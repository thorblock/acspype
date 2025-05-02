import numpy as np
from struct import unpack_from

from acspype.core import WVL_BYTE_OFFSET, PACKET_HEAD, PACKET_TAIL, PAD_BYTE, LPR
from acspype.structures import ACSPacket, ParsedPacket, DeviceCalibratedPacket, TSCorrectedPacket
from acspype.dev import ACSDev
from acspype.tscor import ACSTSCor
from acspype.processing import (convert_sn_hexdec, convert_sn_str, compute_internal_temperature,
                                compute_external_temperature, compute_uncorrected, compute_measured,
                                find_discontinuity_index, compute_discontinuity_offset, apply_discontinuity_offset,
                                ts_correction, zero_shift_correction)


def unpack_packet(full_packet: bytes | bytearray) -> tuple[tuple, int]:
    """
    Unpack a full ACS packet into raw bytes and the checksum. This function is intended for use on a single packet.
        If attempting to implement with Xarray, you will need to create a separate function or attempt to wrap this
        function in an xr.apply_ufunc.
    
    :param full_packet: A bytes or bytesarray object containing the full ACS packet, including the registration
        bytes and the pad byte.
    :return: A tuple containing the raw bytes and the checksum.
    """

    [nwvls] = unpack_from('B', full_packet, offset=WVL_BYTE_OFFSET)  # The number of wavelengths is byte 31.
    remaining_packet_size = int(nwvls * 2 * 2)
    full_packet_descriptor = PACKET_HEAD + f"{remaining_packet_size}H" + PACKET_TAIL
    pad_byte_idx = full_packet.rfind(PAD_BYTE)
    checksum_bytes = full_packet[pad_byte_idx - 2:pad_byte_idx]
    checksum = int(unpack_from('!H', checksum_bytes)[0])
    raw = unpack_from(full_packet_descriptor, full_packet)  # Note: The padbyte is not returned when unpacking.
    return raw, checksum


def parse_packet(acs_packet: ACSPacket) -> ParsedPacket:
    """
    Parse an ACS packet into a ParsedPacket object. This function is intended for use on a single ACSPacket object.
        If attempting to implement with Xarray, you will need to create a separate function or attempt to wrap this
        function in an xr.apply_ufunc.
    
    :param acs_packet: An ACSPacket object containing timestamp and the full packet data. Note: Users looking to 
        only use this to parse the packet without caring about the daq_time can use an ACSPacket with a dummy daq_time.
    :return: A ParsedPacket object containing the parsed data.
    """

    raw, checksum = unpack_packet(acs_packet.full_packet)
    nwvls = raw[LPR + 13]
    parsed_packet = ParsedPacket(
        daq_time=acs_packet.daq_time,
        record_length=raw[LPR + 0],
        packet_type=raw[LPR + 1],
        # reserved_1 = raw[LPR + 2],  # reserved_1 is not used in the ACS.
        sn_int=raw[LPR + 3],
        a_reference_dark=raw[LPR + 4],
        # raw_pressure = raw[LPR + 5], # raw_pressure is deprecated and no longer used.
        a_signal_dark=raw[LPR + 6],
        raw_external_temperature=raw[LPR + 7],
        raw_internal_temperature=raw[LPR + 8],
        c_reference_dark=raw[LPR + 9],
        c_signal_dark=raw[LPR + 10],
        elapsed_time=raw[LPR + 11],
        # reserved_2= raw[LPR + 12], # reserved_2 is not used in the ACS.
        number_of_wavelengths=nwvls,
        c_reference=raw[LPR + 14:LPR + 14 + (nwvls * 4):4],
        a_reference=raw[LPR + 15:LPR + 15 + (nwvls * 4):4],
        c_signal=raw[LPR + 16:LPR + 16 + (nwvls * 4):4],
        a_signal=raw[LPR + 17:LPR + 17 + (nwvls * 4):4],
        checksum=checksum,
    )
    return parsed_packet


def calibrate_packet(parsed_packet: ParsedPacket, dev: ACSDev) -> DeviceCalibratedPacket:
    """
    Calibrate a parsed packet using the device calibration parameters. This function is intended for use on a single
        ParsedPacket object. If attempting to implement with Xarray, you will need to create a separate function or
        attempt to wrap this function in an xr.apply_ufunc.
        
    :param parsed_packet: A ParsedPacket object containing the parsed data. 
    :param dev: An ACSDev-like object containing the device calibration parameters.
    :return: A DeviceCalibratedPacket object containing data corrected/calibrated to the values in the ACSDev object.
    """
    internal_temperature = compute_internal_temperature(parsed_packet.raw_internal_temperature)
    external_temperature = compute_external_temperature(parsed_packet.raw_external_temperature)

    a_uncorrected = compute_uncorrected(parsed_packet.a_signal, parsed_packet.a_reference, dev.path_length)
    c_uncorrected = compute_uncorrected(parsed_packet.c_signal, parsed_packet.c_reference, dev.path_length)
    disc_idx = find_discontinuity_index(dev.a_wavelength, dev.c_wavelength)

    a_m_discontinuity = compute_measured(a_uncorrected, internal_temperature, dev.a_offset, dev.func_a_delta_t)
    c_m_discontinuity = compute_measured(c_uncorrected, internal_temperature, dev.c_offset, dev.func_c_delta_t)

    a_disc_offset = compute_discontinuity_offset(a_m_discontinuity, dev.a_wavelength, disc_idx, 'a_wavelength')
    c_disc_offset = compute_discontinuity_offset(c_m_discontinuity, dev.c_wavelength, disc_idx, 'c_wavelength')

    a_m = apply_discontinuity_offset(a_m_discontinuity, a_disc_offset, disc_idx, 'a_wavelength')
    c_m = apply_discontinuity_offset(c_m_discontinuity, c_disc_offset, disc_idx, 'c_wavelength')

    dev_cal_packet = DeviceCalibratedPacket(
        daq_time=parsed_packet.daq_time,
        a_wavelength=dev.a_wavelength,
        c_wavelength=dev.c_wavelength,
        sn_hexdec=convert_sn_hexdec(parsed_packet.sn_int),
        serial_number=convert_sn_str(parsed_packet.sn_int),
        internal_temperature=internal_temperature,
        external_temperature=external_temperature,
        a_uncorrected=a_uncorrected,
        c_uncorrected=c_uncorrected,
        discontinuity_index=disc_idx,
        a_discontinuity_offset=a_disc_offset,
        c_discontinuity_offset=c_disc_offset,
        a_m_discontinuity=a_m_discontinuity,
        c_m_discontinuity=c_m_discontinuity,
        a_m=a_m,
        c_m=c_m,
    )
    return dev_cal_packet


def ts_correct_packet(device_calibrated_packet: DeviceCalibratedPacket,
                      temperature: float,
                      salinity: float,
                      dev: ACSDev) -> TSCorrectedPacket:
    """
    Apply TS-correction to a device calibrated packet.

    :param device_calibrated_packet: A DeviceCalibratedPacket object containing the device calibrated data.
    :param temperature: The temperature or conservative temperature of the water sample in degrees Celsius.
    :param salinity: The practical or absolute salinity of the water sample.
    :param dev: A ACSDev-like object containing the device calibration parameters.
    :return: A TSCorrectedPacket object containing the TS-corrected data.
    """

    tscor = ACSTSCor().to_xarray()
    psi_s_a = tscor.psi_s_a.sel(wavelength=np.array(device_calibrated_packet.a_wavelength), method='nearest')
    a_psi_t = tscor.psi_t.sel(wavelength=np.array(device_calibrated_packet.a_wavelength), method='nearest')
    a_mts = ts_correction(device_calibrated_packet.a_m, temperature, salinity, a_psi_t.values, psi_s_a.values, dev.tcal)

    psi_s_c = tscor.psi_s_c.sel(wavelength=np.array(device_calibrated_packet.c_wavelength), method='nearest')
    c_psi_t = tscor.psi_t.sel(wavelength=np.array(device_calibrated_packet.c_wavelength), method='nearest')
    c_mts = ts_correction(device_calibrated_packet.c_m, temperature, salinity, c_psi_t.values, psi_s_c.values, dev.tcal)

    a_mts = zero_shift_correction(a_mts)  # Zero shift correction is applied to the a_mts values.
    c_mts = zero_shift_correction(c_mts)

    ts_cor_packet = TSCorrectedPacket(
        daq_time=device_calibrated_packet.daq_time,
        a_wavelength=device_calibrated_packet.a_wavelength,
        c_wavelength=device_calibrated_packet.c_wavelength,
        temperature=temperature,
        salinity=salinity,
        a_mts=a_mts,
        c_mts=c_mts
    )
    return ts_cor_packet
