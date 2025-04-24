from datetime import timedelta
import numpy as np
from struct import unpack_from, calcsize
from typing import Union
import xarray as xr

from acspype.core import WVL_BYTE_OFFSET, PACKET_HEAD, PACKET_TAIL, PACKET_REGISTRATION, PAD_BYTE, LPR, NUM_CHECKSUM_BYTES
from acspype.dev import ACSDev

class FLAG:
    """
    Flag meanings for QAQC tests. This flag syntax follows the QARTOD flag descriptions.

    Users should note that although SUSPECT and HIGH_INTEREST are both flagged as 3, they are not equivalent. Suspect
    flags are used where a returned value is clearly of concern in association with sensor limitations.
    """
    OK: int = 1
    PASS: int = 1
    NOT_EVALUATED: int = 2
    SUSPECT: int = 3
    HIGH_INTEREST: int = 3
    FAIL: int = 4
    MISSING_DATA: int = 9


def gap_test(now, time_stmp, buffer_length, frame_length, time_inc = 0.25):
    if now - time_stmp > timedelta(seconds = time_inc): # Defined in QARTOD Ocean Optics Manual.
        return FLAG.FAIL
    elif buffer_length > frame_length:
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


def syntax_test(packet_data: bytearray) -> int:
    """
    Assess a packet for correct syntax. This function checks if a packet contains the packet registration, the pad byte,
    if the packet length matches the expected length, and if the checksum is correct.

    :param packet_data: The full packet from the ACS. It should include the registration bytes and pad byte.
    :return: The flag indicating the result of the test.
    """

    [nwvls] = unpack_from('B', packet_data, offset=WVL_BYTE_OFFSET)
    remaining_packet_size = int(
        nwvls * 2 * 2)  # Two arrays of signal and two arrays of reference values (four arrays of values total)
    full_packet_descriptor = PACKET_HEAD + f"{remaining_packet_size}H" + PACKET_TAIL
    raw = unpack_from(full_packet_descriptor, packet_data)  # Note: The padbyte is not return when unpacking.
    expected_length = raw[LPR + 0]
    pad_byte_idx = packet_data.rfind(PAD_BYTE)
    checksum_bytes = packet_data[pad_byte_idx - NUM_CHECKSUM_BYTES:pad_byte_idx]
    checksum = unpack_from('!H', checksum_bytes)
    if packet_data[:len(PACKET_REGISTRATION)] != PACKET_REGISTRATION:  # If no packet registration...
        return FLAG.FAIL
    elif packet_data[-1] != int.from_bytes(PAD_BYTE):  # If no pad byte...
        return FLAG.FAIL
    elif calcsize(full_packet_descriptor) != len(packet_data):  # If the descriptor doesn't match the packet length...
        return FLAG.FAIL
    elif len(packet_data) - calcsize(PACKET_TAIL) != expected_length:  # If the expected record length doesn't match.../
        return FLAG.FAIL
    elif checksum != np.uint16(sum(packet_data[:-calcsize(PACKET_TAIL)])):  # If the checksum doesn't add up....
        return FLAG.FAIL
    else:
        return FLAG.PASS


def elapsed_time_test(elapsed_time: Union[int,np.array, xr.DataArray], fail_threshold: int = 1000 * 60,
                      suspect_threshold: int = 1000 * 60 * 3) -> Union[int,np.array, xr.DataArray]:
    """
    Assess if the elapsed time is within acceptable limits. The ACS has a warm-up period after it is turned on
    where data are considered questionable.
    By default, this test assumes any data collected within the first minute of warm up is of poor quality.
    Data are considered suspect if they occur within the first 1-3 minutes of warm up. Users should still review data
    flag as poor and suspect before excluding it from analysis. The ACS manual states the warmup period is 5-10 minutes.

    :param elapsed_time: The time elapsed since the ACS was turned on.
        This data must be an integer representative of the number of milliseconds since the ACS was turned on.
    :param fail_threshold: The number of milliseconds after turning on the ACS to consider data as poor quality.
    :param suspect_threshold: The number of milliseconds after turning on the ACs to consider data as suspect quality.
    :return: The flag indicating the result of the test.
    """
    if isinstance(elapsed_time,int):
        if elapsed_time <= fail_threshold:
            return FLAG.FAIL
        elif fail_threshold < elapsed_time <= suspect_threshold:
            return FLAG.SUSPECT
        else:
            return FLAG.PASS
    elif isinstance(elapsed_time, xr.DataArray):
        flag = xr.full_like(elapsed_time, FLAG.NOT_EVALUATED)
        flag = flag.where(elapsed_time >= fail_threshold, FLAG.FAIL)
        flag = flag.where((elapsed_time < fail_threshold) | (elapsed_time >= suspect_threshold), FLAG.SUSPECT)
        flag = flag.where(elapsed_time < suspect_threshold, FLAG.PASS)
        return flag


def internal_temperature_test(internal_temperature: Union[float, np.array, xr.DataArray],
                              dev: ACSDev) -> Union[float, np.array, xr.DataArray]:
    """
    Assess if the internal temperature is within the range of the ACS device.
    Exceedance of the internal temperature indicates that the environmental conditions outside of
        the sensor are outside what it has been calibrated for.

    :param internal_temperature: The internal temperature of the ACS device in degree Celsius.
    :param dev: The ACSDev object from a corresponding ACS device file.
    :return:The flag indicating the result of the test.
    """

    try:
        min_t = min(dev.tbin)
        max_t = max(dev.tbin)
    except:
        min_t = min(dev.temperature_bin)
        max_t = max(dev.temperature_bin)

    if isinstance(internal_temperature, float):
        if min_t <= internal_temperature <= max_t:
            return FLAG.PASS
        else:
            return FLAG.SUSPECT
    else:
        flag = xr.full_like(internal_temperature, FLAG.NOT_EVALUATED).astype(int)
        flag = flag.where((min_t < internal_temperature) & (max_t > internal_temperature), FLAG.FAIL)
        flag = flag.where((min_t > internal_temperature) & (max_t < internal_temperature), FLAG.PASS)
        return flag


def inf_nan_test(uncorr: Union[tuple, xr.DataArray]):
    """
    Assess if a spectrum contains any infinite or nan values.

    :param uncorr: Uncorrected absorption or attenuation. Really any spectrum in geophysical units will work.
    :return: A flag indicating the result of the test.
    """
    if isinstance(uncorr, tuple):
        if np.any(np.isinf(uncorr)) or np.any(np.isnan(uncorr)):
            return FLAG.FAIL
        else:
            return FLAG.PASS
    elif isinstance(uncorr, xr.DataArray):
        flag = xr.full_like(uncorr.time.astype(int), FLAG.PASS).astype(int)
        flag = flag.where((~np.any(np.isinf(uncorr),axis = 1)), FLAG.FAIL)
        flag = flag.where((~np.any(np.isnan(uncorr), axis = 1)), FLAG.FAIL)
        return flag



def gross_range_test(mts: xr.DataArray,
                     sensor_min: float = -0.005, sensor_max: float = 10.00,
                     op_min: float =0.001,  op_max: float = 9.0) -> xr.DataArray:
    """
    According to the ACS manual, the dynamic range of the sensor is 0.001 to 10. However, in the processing
    protocols, there is a call for values between -0.005 and 0 to be considered equivalent to 0.

    Thus, the decision was made to flag values outside -0.005 and 10.000 as FAIL
        and values outside 0.001 and 10 as SUSPECT.
    It is important to note that a_m and c_m represent data that has not had temperature-salinity-scattering correction.
    The flags associated with this gross range test should be taken lightly and used as an indicator to assess
    data validity after correction has been applied.
    After correction, if the value is within the anticipated sensor range, then it is probably acceptable. Instances of
    values less than 0.001 is more probable than values greater than 10.

    :param m: Absorption or attenuation data.
    :param sensor_min: Anything less than the sensor minimum will be flagged as fail.
    :param sensor_max: Anything greater than the sensor maximum will be flagged as fail.
    :param op_min: Anything between the sensor_min and op_min will be flagged as suspect.
    :param op_max: Anything between the op_max and sensor_max will be flagged as suspect.
    :return: The flag of the data, which maintains the same shape as the input.
    """

    mts = np.array(mts)
    flag = np.where((mts < sensor_min) | (mts > sensor_max), FLAG.FAIL, FLAG.PASS)
    flag = np.where((mts > op_min) | (mts < op_max), flag, FLAG.SUSPECT)
    return flag

def blanket_gross_range_test(mts: xr.DataArray, wavelength_dim: str, percent_unacceptable: float = 10.0,
                             sensor_min: float = -0.005, sensor_max: float = 10,
                             op_min:float = 0.001, op_max:float = 9.5):
    """
    Assign a blanket gross range flag based on the percentage of values that are flagged as suspect or fail in the spectrum.

    :param mts: Measured and TS corrected absorption or attenuation data
    :param wavelength_dim: The wavelength dimension of the data.
    :param percent_unacceptable: The minimum percent of values that can be flagged as suspect or fail in the spectrum.
    :param sensor_min: The minimum sensor value to flag as fail
    :param sensor_max: The maximum sensor value to flag as fail.
    :param op_min: The minimum value to flag as suspect.
    :param op_max: The maximum value to flag as suspect.
    :return: The blanket flag of the data, which maintains the same size as the time dimension.
    """

    ndflag = gross_range_test(mts,sensor_min, sensor_max, op_min, op_max)
    flag_bool = xr.where(ndflag == 1, 0, 1)
    num_wvls = len(mts[wavelength_dim])
    flag_bool_sum = np.sum(flag_bool,axis = 1)
    flag_percent = (flag_bool_sum / num_wvls) * 100
    flag = xr.full_like(mts.time.astype(int), FLAG.PASS).astype(int)
    flag = flag.where(flag_percent <= percent_unacceptable, FLAG.FAIL)
    return flag


def a_gt_c_test(absorption: xr.DataArray, attenuation: xr.DataArray) -> xr.DataArray:
    """
    Assess if the absorption is greater than the attenuation. Having absorption greater than attenuation is
    (mostly) physically impossible.
    :param absorption: Absorption data.
    :param attenuation: Attenuation data.
    :return: Flag indicating if absorption is greater than attenuation.
    """

    flag = xr.full_like(absorption, FLAG.NOT_EVALUATED).astype(int)
    flag = flag.where(absorption > attenuation, FLAG.FAIL)
    flag = flag.where(absorption <= attenuation, FLAG.PASS)
    return flag