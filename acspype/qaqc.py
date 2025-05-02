from datetime import timedelta
import numpy as np
from numpy.typing import NDArray
from struct import calcsize
import xarray as xr

from acspype.core import (PACKET_HEAD, PACKET_TAIL, LPR, PACKET_REGISTRATION, PAD_BYTE)
from acspype.dev import ACSDev
from acspype.packet import unpack_packet


class FLAG:
    OK: int = 1
    PASS: int = 1
    NOT_EVALUATED: int = 2
    SUSPECT: int = 3
    HIGH_INTEREST: int = 3
    FAIL: int = 4
    MISSING_DATA: int = 9


def gap_test(now, time_stmp, buffer_length, record_length, time_inc=0.25):
    if now - time_stmp > timedelta(seconds=time_inc):  # Defined in QARTOD Ocean Optics Manual.
        return FLAG.FAIL
    elif buffer_length > record_length:
        """
        This is a custom take on the gap test. 
        If for some reason the buffer length exceeds the previous frame length, that would indicate that
        the buffer is filling up faster than the packets can be unpacked. This would ultimately result in the timestamp
        of the value being off by one or multiples of 250ms periods, depending on the buffer length. This would also 
        indicate an issue with the timing of the data acquisition thread.
        """
        return FLAG.FAIL
    else:
        return FLAG.PASS


def syntax_test(full_packet: bytearray) -> int:
    """
    Test the syntax of an incoming packet. This expands on the QARTOD syntax test.
    Failure of this test indicates that the packet is malformed or incomplete.

    :param full_packet: A full packet from the ACS, including registration bytes and the pad byte.
    :return: Flag indicating pass or fail.
    """
    raw, checksum = unpack_packet(full_packet)
    remaining_packet_size = int(raw[LPR + 13] * 2 * 2)
    full_packet_descriptor = PACKET_HEAD + f"{remaining_packet_size}H" + PACKET_TAIL
    expected_length = raw[LPR + 0]
    if full_packet[:len(PACKET_REGISTRATION)] != PACKET_REGISTRATION:  # If no packet registration...
        return FLAG.FAIL
    elif full_packet[-1] != int.from_bytes(PAD_BYTE):  # If no pad byte...
        return FLAG.FAIL
    elif calcsize(full_packet_descriptor) != len(full_packet):  # If the descriptor doesn't match the packet length...
        return FLAG.FAIL
    elif len(full_packet) - calcsize(PACKET_TAIL) != expected_length:  # If the expected record length doesn't match...
        return FLAG.FAIL
    elif sum(full_packet[:-calcsize(PACKET_TAIL)]) >= 65536:  # If the checksum doesn't add up....
        if checksum != sum(full_packet[:-calcsize(PACKET_TAIL)]) % 65536:
            return FLAG.FAIL
        else:
            return FLAG.PASS
    elif sum(full_packet[:-calcsize(PACKET_TAIL)]) < 65536:
        if checksum != np.uint16(sum(full_packet[:-calcsize(PACKET_TAIL)])):
            return FLAG.FAIL
        else:
            return FLAG.PASS
    else:
        return FLAG.PASS


def elapsed_time_test(elapsed_time: int | xr.DataArray, fail_threshold: int = 1000 * 60,
                      suspect_threshold: int = 1000 * 60 * 3) -> int | xr.DataArray:
    """
    Assess the elapsed time since the instrument was turned on. This is a custom test that is based on information
    in the ACS manual. Generally, the ACS takes 5-10 minutes to warm-up. Waiting this long may not be practical in
    some situations. It is recommended to flag the first minute of data as poor quality and then minutes 2-3 as suspect.

    :param elapsed_time: The elapsed time parsed from an ACS packet.
        Represents the number of milliseconds that have passed since the ACS was turned on.
    :param fail_threshold: The amount of time in milliseconds where data is considered to be of poor quality.
    :param suspect_threshold: The amount of time in milliseconds between the fail and
        suspect threshold where data is considered to be suspect.
    :return: Flag indicating pass, suspect, or fail.
    """
    if not isinstance(elapsed_time, xr.DataArray):
        if elapsed_time <= fail_threshold:
            return FLAG.FAIL
        elif fail_threshold < elapsed_time <= suspect_threshold:
            return FLAG.SUSPECT
        else:
            return FLAG.PASS
    else:
        flags = xr.full_like(elapsed_time, FLAG.NOT_EVALUATED)
        flags = xr.where(elapsed_time <= fail_threshold, FLAG.FAIL, flags)
        flags = xr.where((elapsed_time > fail_threshold) & (elapsed_time <= suspect_threshold), FLAG.SUSPECT, flags)
        flags = xr.where(elapsed_time > suspect_threshold, FLAG.PASS, flags)
        return flags


def internal_temperature_test(internal_temperature: float | xr.DataArray,
                              dev: ACSDev) -> float | xr.DataArray:
    """
    Assess the internal temperature of the ACS. This is a custom test that is based on information within the manual
    and the device file. The ACS is calibrated at specific temperature bins, which are listed in the device file.
    Typically, the internal temperature correction is linearly interpolated from the delta T values in the device file.
    If the internal temperature of the sensor exceed the range of these temperature bins,
    then data are considered to be suspect.

    :param internal_temperature: The internal temperature of the ACS.
    :param dev: The device file for the ACS as an ACSDev object.
    :return: Flag indicating pass or suspect.
    """

    min_t = min(dev.tbin)
    max_t = max(dev.tbin)
    if not isinstance(internal_temperature, xr.DataArray):
        if min_t <= internal_temperature <= max_t:
            return FLAG.PASS
        else:
            return FLAG.SUSPECT
    else:
        flags = xr.full_like(internal_temperature, FLAG.NOT_EVALUATED).astype(int)
        flags = xr.where((min_t > internal_temperature) & (max_t < internal_temperature), FLAG.SUSPECT, flags)
        flags = xr.where((min_t <= internal_temperature) & (max_t >= internal_temperature), FLAG.PASS, flags)
        return flags


def inf_nan_test(uncorrected: NDArray[float] | xr.DataArray) -> int | xr.DataArray:
    """
    Assess if the uncorrected data contains any NaN or Inf values. This is a custom test that can be used to identify
    instances where the reference counts are zero, which produces mathematical errors in the conversion to geophysical
    units. Zero reference counts happen occasionally and can probably be discarded during analysis. Repeated samples
    containing zero reference counts may indicate a problem with the instrument.

    :param uncorrected: Uncorrected a or c data.
    :return: Flag indicating pass or fail.
    """
    if not isinstance(uncorrected, xr.DataArray):
        if np.any(np.isinf(uncorrected)) or np.any(np.isnan(uncorrected)):
            return FLAG.FAIL
        else:
            return FLAG.PASS
    else:
        flags = xr.full_like(uncorrected.time.astype(int), FLAG.PASS).astype(int)
        flags = xr.where((np.any(np.isinf(uncorrected), axis=1)), FLAG.FAIL, flags)
        flags = xr.where((np.any(np.isnan(uncorrected), axis=1)), FLAG.FAIL, flags)
        return flags


def gross_range_test(mts: xr.DataArray,
                     sensor_min: float = -0.005, sensor_max: float = 10.00,
                     op_min: float = 0.001, op_max: float = 8.5) -> xr.DataArray:
    """
    Assess if the TS-corrected data are within the limitations of the instrument. From the ACS manual, the valid range
    of the sensor is 0-10 m^-1.

    The return is a set of flags in the same shape as the input data. Note that fail flags may occur in the
    red wavelengths of absorption spectrum consecutively toward the end of each spectrum, which does not necessarily
    make the entire spectrum of poor quality. Users should consider the pass/fail state of the rest of the spectrum
    before making a decision on the validity of the data.

    :param mts: TS-corrected absorption or attenuation.
    :param sensor_min: The minimum sensor range. Default is -0.005, which is effectively 0 as described in the manual.
    :param sensor_max: The maximum sensor range. Default is 10.00, which is the maximum range of the ACS.
    :param op_min: Operator set minimum for flagging suspect data. Default is 0.001
    :param op_max: Operator set maximum for flagging suspect data. Default is 8.5.
    :return: Flag indicating pass, fail, or suspect.
    """

    if not isinstance(mts, xr.DataArray):
        raise NotImplementedError("Functionality for singletons not yet implemented.")

    else:
        flags = xr.full_like(mts, FLAG.NOT_EVALUATED).astype(int)
        flags = flags.where((mts > sensor_min) & (mts < sensor_max), FLAG.FAIL)
        flags = flags.where((mts > op_min) | (mts < sensor_min), FLAG.SUSPECT)
        flags = flags.where((mts < op_max) | (mts > sensor_max), FLAG.SUSPECT)
        flags = flags.where((mts <= op_min) | (mts >= op_max), FLAG.PASS)
        return flags


# TODO: Test below functions.
#
# def blanket_gross_range_test(mts: xr.DataArray, wavelength_dim: str, percent_unacceptable: float = 10.0,
#                              sensor_min: float = -0.005, sensor_max: float = 10,
#                              op_min:float = 0.001, op_max:float = 9.5):
#     """
#     Assign a blanket gross range flag based on the percentage of values that are flagged as suspect or fail in the spectrum.
#
#     :param mts: Measured and TS corrected absorption or attenuation data
#     :param wavelength_dim: The wavelength dimension of the data.
#     :param percent_unacceptable: The minimum percent of values that can be flagged as suspect or fail in the spectrum.
#     :param sensor_min: The minimum sensor value to flag as fail
#     :param sensor_max: The maximum sensor value to flag as fail.
#     :param op_min: The minimum value to flag as suspect.
#     :param op_max: The maximum value to flag as suspect.
#     :return: The blanket flag of the data, which maintains the same size as the time dimension.
#     """
#
#     ndflag = gross_range_test(mts,sensor_min, sensor_max, op_min, op_max)
#     flag_bool = xr.where(ndflag == 1, 0, 1)
#     num_wvls = len(mts[wavelength_dim])
#     flag_bool_sum = np.sum(flag_bool,axis = 1)
#     flag_percent = (flag_bool_sum / num_wvls) * 100
#     flag = xr.full_like(mts.time.astype(int), FLAG.PASS).astype(int)
#     flag = flag.where(flag_percent <= percent_unacceptable, FLAG.FAIL)
#     return flag
#
#
# def a_gt_c_test(absorption: xr.DataArray, attenuation: xr.DataArray) -> xr.DataArray:
#     """
#     Assess if the absorption is greater than the attenuation. Having absorption greater than attenuation is
#     (mostly) physically impossible.
#     :param absorption: Absorption data.
#     :param attenuation: Attenuation data.
#     :return: Flag indicating if absorption is greater than attenuation.
#     """
#
#     flag = xr.full_like(absorption, FLAG.NOT_EVALUATED).astype(int)
#     flag = flag.where(absorption > attenuation, FLAG.FAIL)
#     flag = flag.where(absorption <= attenuation, FLAG.PASS)
#     return flag
