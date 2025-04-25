import numpy as np
from struct import unpack_from
from typing import Union, Tuple
import xarray as xr

from acspype.core import LPR, PACKET_TAIL, PACKET_HEAD, PAD_BYTE, WVL_BYTE_OFFSET
from acspype.dev import ACSDev
from acspype.tscor import ACSTSCor
from acspype.structures import ACSPacket, ParsedPacket, CalibratedPacket
from acspype.qaqc import FLAG, internal_temperature_test, syntax_test, elapsed_time_test
from acspype.discontinuity import (find_discontinuity_index, _apply_discontinuity_offset,
                                   _compute_discontinuity_offset, discontinuity_correction)



def parse_packet(acs_packet: ACSPacket,
                 run_syntax_test: bool = True,
                 run_elapsed_time_test: bool = True) -> ParsedPacket:
    """
    Parse an ACSPacket object into a ParsedPacket object.

    :param acs_packet: The ACSPacket object, which contains daq_time and full_packet data.
    :param run_syntax_test:  If True, runs the syntax test on the packet.
    :param run_elapsed_time_test: If True, runs the elapsed time test on the packet.
    :return: A ParsedPacket object containing the parsed data.
    """

    daq_time, full_packet = acs_packet
    if run_syntax_test is True:
        flag_syntax = syntax_test(full_packet)
    else:
        flag_syntax = FLAG.NOT_EVALUATED
    [nwvls] = unpack_from('B', full_packet, offset=WVL_BYTE_OFFSET)
    remaining_packet_size = int(nwvls * 2 * 2)
    full_packet_descriptor = PACKET_HEAD + f"{remaining_packet_size}H" + PACKET_TAIL
    pad_byte_idx = full_packet.rfind(PAD_BYTE)
    checksum_bytes = full_packet[pad_byte_idx - 2:pad_byte_idx]
    checksum = unpack_from('!H', checksum_bytes)
    if isinstance(checksum, tuple) and len(checksum) == 1:
        checksum = checksum[0]
    raw = unpack_from(full_packet_descriptor, full_packet)  # Note: The padbyte is not returned when unpacking.

    # Split and process the raw data.
    reg_bytes = raw[:LPR]
    record_length = raw[LPR + 0]
    packet_type = raw[LPR + 1]
    reserved_1 = raw[LPR + 2]
    sn_int = raw[LPR + 3]
    a_reference_dark = raw[LPR + 4]
    raw_pressure = raw[LPR + 5]
    a_signal_dark = raw[LPR + 6]
    raw_external_temp = raw[LPR + 7]
    raw_internal_temp = raw[LPR + 8]
    c_reference_dark = raw[LPR + 9]
    c_signal_dark = raw[LPR + 10]
    elapsed_time = raw[LPR + 11]
    reserved_2 = raw[LPR + 12]
    nwvls = raw[LPR + 13]
    c_reference = raw[LPR + 14:LPR + 14 + (nwvls * 4):4]
    a_reference = raw[LPR + 15:LPR + 15 + (nwvls * 4):4]
    c_signal = raw[LPR + 16:LPR + 16 + (nwvls * 4):4]
    a_signal = raw[LPR + 17:LPR + 17 + (nwvls * 4):4]

    sn_hexdec = hex(sn_int)
    sn = 'ACS-' + str(int(hex(sn_int)[-6:], 16)).zfill(5)  # sn string becomes ACS-XXXXX. Zero padded to 5 digits.

    if run_elapsed_time_test is True:
        flag_elapsed_time = elapsed_time_test(elapsed_time)
    else:
        flag_elapsed_time = FLAG.NOT_EVALUATED


    parsed_packet = ParsedPacket(daq_time=daq_time, registration_bytes=reg_bytes, record_length=record_length,
                                 packet_type=packet_type, reserved_1=reserved_1, serial_number_int=sn_int,
                                 number_of_output_wavelengths=nwvls, a_reference_dark=a_reference_dark,
                                 c_reference_dark=c_reference_dark, raw_external_temperature=raw_external_temp,
                                 raw_internal_temperature=raw_internal_temp, a_signal_dark=a_signal_dark,
                                 c_signal_dark=c_signal_dark, elapsed_time=elapsed_time, reserved_2=reserved_2,
                                 c_reference=c_reference, a_reference=a_reference, c_signal=c_signal, a_signal=a_signal,
                                 raw_pressure=raw_pressure, checksum=checksum, serial_number_hexdec=sn_hexdec,
                                 serial_number=sn, flag_syntax=flag_syntax, flag_elapsed_time=flag_elapsed_time)
    return parsed_packet


def calibrate_packet(parsed_packet: ParsedPacket, dev: ACSDev,
                     flag_internal_temperature: bool = True) -> CalibratedPacket:
    """
    Calibrate a ParsedPacket object using the supplied ACS device file.
    :param parsed_packet: The ParsedPacket object that contains raw data.
    :param dev: The device file object to apply to the parsed packet.
    :return: Data calibrated using the device file.
    """

    a_uncorr = tuple(map(float,compute_uncorrected(np.array(parsed_packet.a_signal),
                                                   np.array(parsed_packet.a_reference), dev)))
    c_uncorr = tuple(map(float,compute_uncorrected(np.array(parsed_packet.c_signal),
                                                   np.array(parsed_packet.c_reference), dev)))

    internal_temp = float(compute_internal_temperature(parsed_packet.raw_internal_temperature))

    if flag_internal_temperature is True:
        flag_internal_temp = internal_temperature_test(internal_temp, dev)
    else:
        flag_internal_temp = FLAG.NOT_EVALUATED

    external_temp = float(compute_external_temperature(parsed_packet.raw_external_temperature))

    disc_idx = find_discontinuity_index(dev.a_wavelength, dev.c_wavelength)

    a_disc = tuple(map(float,compute_measured(a_uncorr, 'a', internal_temp, dev)))
    c_disc = tuple(map(float,compute_measured(c_uncorr, 'c', internal_temp, dev)))

    a_discontinuity_offset = _compute_discontinuity_offset(a_disc, dev.a_wavelength, disc_idx)
    c_discontinuity_offset = _compute_discontinuity_offset(c_disc, dev.c_wavelength, disc_idx)

    a_m = tuple(map(float,_apply_discontinuity_offset(a_disc, a_discontinuity_offset, disc_idx)))
    c_m = tuple(map(float,_apply_discontinuity_offset(c_disc, c_discontinuity_offset, disc_idx)))

    calibrated_packet = CalibratedPacket(daq_time=parsed_packet.daq_time,
                                         serial_number=parsed_packet.serial_number,
                                         flag_syntax=parsed_packet.flag_syntax,
                                         elapsed_time=parsed_packet.elapsed_time,
                                         flag_elapsed_time=parsed_packet.flag_elapsed_time,
                                         internal_temperature=internal_temp,
                                         flag_internal_temperature= flag_internal_temp,
                                         external_temperature=external_temp,
                                         a_wavelength=tuple(map(float,dev.a_wavelength)),
                                         c_wavelength=tuple(map(float,dev.c_wavelength)),
                                         a_uncorrected=a_uncorr,
                                         c_uncorrected=c_uncorr,
                                         discontinuity_wavelength_index= disc_idx,
                                         a_discontinuity_offset=a_discontinuity_offset,
                                         c_discontinuity_offset= c_discontinuity_offset,
                                         a_m_discontinuity=a_disc,
                                         c_m_discontinuity=c_disc,
                                         a_m=a_m,
                                         c_m=c_m)
    return calibrated_packet


def compute_internal_temperature(counts: Union[int,xr.DataArray, np.array]) -> Union[np.array,xr.DataArray]:
    """
    Compute internal sensor housing temperature in degrees Celsius. The formula for conversion can be found in the
    ACS Manual.
    :param counts: The count output of the internal thermistor.
    :return: The converted data, with a shape/size that matches the input data.
    """

    volts = 5 * counts / 65535
    resistance = 10000 * volts / (4.516 - volts)
    internal_temperature = 1 / (
                0.00093135 + 0.000221631 * np.log(resistance) + 0.000000125741 * np.log(resistance) ** 3) - 273.15
    return internal_temperature


def compute_external_temperature(counts: Union[int,xr.DataArray,np.array]) -> Union[float,np.array,xr.DataArray]:
    """
    Compute external temperature in degrees Celsius. This thermistor has contact with external air or fluid.
    The formula for conversion can be found in the ACS Manual.

    :param counts: The count output of the external thermistor.
    :return: The converted data, with a shape/size that matches the input data.
    """
    a = -7.1023317e-13
    b = 7.09341920e-08
    c = -3.87065673e-03
    d = 95.8241397
    external_temperature = a * counts ** 3 + b * counts ** 2 + c * counts + d
    return external_temperature


def compute_uncorrected(signal_counts: Union[tuple, np.array, xr.DataArray],
                        reference_counts: Union[tuple, np.array, xr.DataArray],
                        dev: ACSDev) -> Union[np.array, xr.DataArray]:
    """
    Compute uncorrected absorption and attenuation from the sensor counts.
    This value should not be used for scientific analysis, but can be used to diagnose issues, such as bubble intrusion
    and blockages.

    :param signal_counts: Raw signal counts from absorption or attenuation channels. Note that it is not necessary for
    users to correct for dark values, as this is already performed internally by the sensor.
    :param reference_counts: The reference absorption and attenuation counts output by the ACS sensor.
    :param dev: The ACSDev object derived from the sensor device file.
    :return: Uncorrected absorption and attenuation values with the units of inverse meters.
    """

    uncorr = (1 / dev.path_length) * np.log(signal_counts / reference_counts)
    return uncorr




def compute_measured(uncorrected: list,
                     channel: str,
                     internal_temperature: float,
                     dev: ACSDev):
    """
    Compute the measured absorption or attenuation from the uncorrected data. The data are corrected using the offsets
    found in the supplied device file. If using the factory device file, data are corrected for the effect of pure water
    on absorption and attenuation and the internal temperature of the sensor.
    :param uncorrected: The uncorrected values in inverse meters.
    :param channel: The channel or tube the data originates from. 'a' for absorption, 'c' for attenuation.
    :param internal_temperature: The internal temperature of the sensor.
    :param dev: A device file represented as a ACSDev object.
    :return: The data corrected for the offset found in the device file.
    """
    if channel.lower() == 'a':
        delta_t = dev.func_a_delta_t(internal_temperature).T
        offsets = dev.a_offset
    elif channel.lower() == 'c':
        delta_t = dev.func_c_delta_t(internal_temperature).T
        offsets = dev.c_offset
    measured = (offsets - uncorrected) - delta_t
    return measured





def ts_correction(m: xr.DataArray,
                  channel: str,
                  temperature: xr.DataArray,
                  salinity: xr.DataArray,
                  dev: ACSDev,
                  tscor: ACSTSCor) -> xr.DataArray:
    """
    Compute temperature-salinity corrected data.

    :param m: The measured signal, corrected for pure water offsets. Does not yet imply the filtration state of the water.
    :param channel: 'a' or 'c', used for obtaining the correct correction coeffs.
    :param temperature: sea water temperature. conservative temperature is acceptable per the ACS manual,
        but does not result in much change
    :param salinity: sea water salinity. absolute salinity is acceptable per the ACS manual,
        but does not result in much change.
    :param dev: The corresponding ACSDev object for the dataset.
    :param tscor: The corresponding ACSTSCor object for the dataset.

    :return: An xr.DataArray containing temperature-salinity corrected data.
    """

    _channel = channel.lower()
    _tcal = dev.tcal

    if _channel == 'a':
        _tscor = tscor.to_xarray().sel(wavelength=m.a_wavelength, method='nearest')
        delta_t, psi_temp = np.meshgrid(temperature - _tcal, _tscor.psi_t)
        s, psi_sal = np.meshgrid(salinity, _tscor.psi_s_a)
    elif _channel == 'c':
        _tscor = tscor.to_xarray().sel(wavelength=m.c_wavelength, method='nearest')
        delta_t, psi_temp = np.meshgrid(temperature - _tcal, _tscor.psi_t)
        s, psi_sal = np.meshgrid(salinity, _tscor.psi_s_c)
    mts = m - ((psi_temp.T * delta_t.T) + (psi_sal.T * s.T))
    return mts


def zero_shift(mts: xr.DataArray) -> xr.DataArray:
    """
    Zero out any values that are between -0.005 and 0. According to the manual, these values are equivalent to 0.

    :param mts: The absorption or attenuation, corrected for temperature and salinity.
    :return: The zeroed out data.
    """

    mts = mts.where((mts > 0) | (mts <= -0.005), 0)
    return mts


def interpolate_common_wavelengths(ds: xr.Dataset, step: int = 1,
                                   wavelength_range: list or str = 'infer') -> xr.Dataset:
    if wavelength_range == 'infer':
        min_wvl = np.ceil(max(ds.a_wavelength.min(), ds.c_wavelength.min())) # Maximum of a and c wavelength minimums.
        max_wvl = np.floor(min(ds.a_wavelength.max(), ds.c_wavelength.max())) # Minimum of a and c wavelength maximums.
    else:
        min_wvl, max_wvl = wavelength_range
    wvls = np.arange(min_wvl, max_wvl, step)
    cds = ds.interp({'a_wavelength': wvls, 'c_wavelength': wvls})  # Interpolate to step size.
    cds = cds.reset_index(['a_wavelength',
                           'c_wavelength'], drop=True).assign_coords(wavelength=wvls).rename(
        {'a_wavelength': 'wavelength',
         'c_wavelength': 'wavelength'})

    cds.attrs['interpolation_step'] = step
    return cds


def scattering_correction_baseline(a_mts: xr.DataArray, reference_wavelength: int = 715):
    """
    METHOD 1
    Perform baseline scattering correction on absorption data.
    Although this function will work with data that is not on a common wavelength bin,
    it is recommended that you apply this function with wavelength interpolated data.

    :param a_mts: TS corrected absorption.
    :param reference_wavelength: The reference wavelength. Usually 715.
    :return: Absorption corrected using the baseline method.
    """
    if 'wavelength' in a_mts.dims:
        ref = a_mts.sel(wavelength=reference_wavelength, method='nearest')
    else:
        ref = a_mts.sel(a_wavelength=reference_wavelength, method='nearest')
    scatcorr = a_mts - ref

    scatcorr.attrs['reference_wavelength'] = reference_wavelength
    scatcorr.attrs['method'] = 'Baseline Scattering Correction - Method 1'
    return scatcorr


def scattering_correction_fixed(a_mts: xr.DataArray, c_mts: xr.DataArray, F: float = 0.18) -> xr.DataArray:
    """
    METHOD 2
    Perform fixed scattering correction on absorption data.
    Use of this method requires that the absorption and attenuation data be on common wavelength bins.

    :param a_mts: TS corrected absorption.
    :param c_mts: TS corrected attenuation.
    :param F: The fixed offset.
    :return: Absorption corrected using the baseline method.
    """
    scatcorr = a_mts - F * (c_mts - a_mts)

    scatcorr.attrs['epsilon'] = F
    scatcorr.attrs['method'] = 'Fixed Scattering Correction - Method 2'
    return scatcorr


def scattering_correction_proportional(a_mts: xr.DataArray, c_mts: xr.DataArray,
                                       reference_wavelength: int = 715) -> xr.DataArray:
    """
    METHOD 3
    Perform proportional scattering correction on absorption data.
    Use of this method requires that the absorption and attenuation data be on common wavelength bins.

    :param a_mts: TS corrected absorption.
    :param c_mts: TS corrected attenuation.
    :param reference_wavelength: The reference wavelength. Usually 715.
    :return: Absorption corrected using the baseline method.
    """

    ref_a = a_mts.sel(wavelength=reference_wavelength, method='nearest')
    ref_c = c_mts.sel(wavelength=reference_wavelength, method='nearest')
    scatcorr = a_mts - ((ref_a / (ref_c - ref_a)) * (c_mts - a_mts))

    scatcorr.attrs['reference_wavelength'] = reference_wavelength
    scatcorr.attrs['method'] = 'Proportional Scattering Correction - Method 3'
    return scatcorr
