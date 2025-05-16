from datetime import datetime, timedelta
import numpy as np
from numpy.typing import NDArray
from struct import calcsize
import xarray as xr

from acspype.core import (PACKET_HEAD, PACKET_TAIL, LPR, PACKET_REGISTRATION, PAD_BYTE)
from acspype.dev import ACSDev
from acspype.packet import unpack_packet


class FLAG:
    """
    Flag values for the ACS QA/QC tests. These flags follow the QARTOD flag meanings.

    :param OK: Indicates that the data passed the test.
    :param PASS: Indicates that the data passed the test.
    :param NOT_EVALUATED: Indicates that the data was not evaluated. Appearance of this flag indicates a problem with
        the test function and not the data.
    :param SUSPECT: Indicates that the data is suspect.
        This flag indicates that the data should be reviewed more closely.
    :param HIGH_INTEREST: Indicates that the data is of high interest. Although it uses the same flag as SUSPECT, the
        context of the test and environmental conditions should be considered when reviewing the data.
    :param FAIL: Indicates that the data failed the test.
    :param BAD: Indicates that the data is bad and failed the test.
    :param MISSING_DATA: Indicates that the input data or ancillary data is missing (NaN) and the test was not run.
    """

    OK: int = 1
    PASS: int = 1
    NOT_EVALUATED: int = 2
    SUSPECT: int = 3
    HIGH_INTEREST: int = 3
    FAIL: int = 4
    BAD: int = 4
    MISSING_DATA: int = 9


def gap_test(now: datetime,
             time_stmp: datetime,
             record_length: int,
             buffer_length: int | None,
             time_inc: float = 0.25) -> int:
    """
    Assess the gap between the current time and the timestamp of the packet. This is a modified form of the generic
    QARTOD gap test.

    :param now: The time at which the test is being run.
    :param time_stmp: The timestamp of the packet.
    :param buffer_length: The number of bytes in the serial buffer. If None, the corresponding section of the test is
        skipped.
    :param record_length: The record length of the packet. This is the expected number of bytes in the packet.
    :param time_inc: The time increment to use for assessing timestamp gaps. Default is 0.25 seconds, which is the
        approximate time that it takes the ACS to send a complete packet.
    :return: A flag indicating pass or fail.
    """
    if now - time_stmp > timedelta(seconds=time_inc):  # Defined in QARTOD Ocean Optics Manual.
        return FLAG.FAIL
    elif buffer_length is not None:
        if buffer_length > record_length:
            """
            This is a custom take on the gap test. 
            If for some reason the buffer length exceeds the previous frame length, that would indicate that
            the buffer is filling up faster than the packets can be unpacked. This would ultimately result in the 
            timestamp of the value being off by one or multiples of 250ms periods, depending on the buffer length. 
            This could also indicate an issue with the timing of the data acquisition thread.
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

        # Assign attributes to the flags if an xarray.DataArray.
        flags.attrs['ancillary_variables'] = elapsed_time.name
        flags.attrs['fail_threshold'] = fail_threshold
        flags.attrs['suspect_threshold'] = suspect_threshold
        flags.attrs['threshold_units'] = 'milliseconds'
        flags.attrs['test_name'] = 'elapsed_time_test'
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

        # Assign attributes to the flags if an xarray.DataArray.
        flags.attrs['ancillary_variables'] = internal_temperature.name
        flags.attrs['minimum_temperature_bin'] = min_t
        flags.attrs['maximum_temperature_bin'] = max_t
        flags.attrs['temperature_bin_units'] = 'degrees Celsius'
        flags.attrs['test_name'] = 'internal_temperature_test'
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

        # Assign attributes to the flags if an xarray.DataArray.
        flags.attrs['ancillary_variables'] = uncorrected.name
        flags.attrs['test_name'] = 'inf_nan_test'
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

        # Assign attributes to the flags if an xarray.DataArray.
        flags.attrs['ancillary_variables'] = mts.name
        flags.attrs['sensor_min'] = sensor_min
        flags.attrs['sensor_max'] = sensor_max
        flags.attrs['operator_min'] = op_min
        flags.attrs['operator_max'] = op_max
        flags.attrs['threshold_units'] = 'm^-1'
        flags.attrs['test_name'] = 'gross_range_test'
        return flags


def discontinuity_offset_test(discontinuity_offset: xr.DataArray,
                              median_multiplier: float = 3,
                              fail_threshold: float = 10) -> xr.DataArray:
    """
    Flag a discontinuity offset value as pass or suspect. This is a custom test that uses a multiple of the median
    of the offset to assess if the offset is acceptable.

    Spectra with significantly large discontinuity offsets are likely to be of poor quality.

    :param discontinuity_offset: The discontinuity offset value for absorption or attenuation.
    :param median_multiplier: The multiplier of the median. Values within the multiplier range are deemed acceptable.
    :param fail_threshold: The threshold for flagging the offset as fail.
        Default is 10, or the maximum value of the ACS range.
    :return: The flag of the discontinuity offset, which maintains the same size as the time dimension.
    """
    flags = xr.full_like(discontinuity_offset, 1).astype(int)
    _median = discontinuity_offset.median(skipna=True)
    flags = flags.where((np.abs(discontinuity_offset) < np.abs(_median) * median_multiplier), FLAG.SUSPECT)
    flags = flags.where((np.abs(discontinuity_offset) < fail_threshold, FLAG.FAIL))

    # Assign attributes to the flags if an xarray.DataArray.
    flags.attrs['ancillary_variables'] = discontinuity_offset.name
    flags.attrs['median_multiplier'] = median_multiplier
    flags.attrs['fail_threshold'] = fail_threshold
    flags.attrs['suspect_threshold'] = np.abs(_median) * median_multiplier
    flags.attrs['threshold_units'] = 'm^-1'
    flags.attrs['test_name'] = 'discontinuity_offset_test'
    return flags


def blanket_gross_range_test(nd_gross_range_results: xr.DataArray,
                             wavelength_dim: str,
                             suspect_threshold: float,
                             fail_threshold: float,
                             include_suspect_flags: bool = False) -> xr.DataArray:
    """
    Assess the validity of the gross range test results as a blanket flag for the spectra.
    This is an experimental test.
    The test considers the number of wavelength bins that were previously flagged as suspect or fail for exceeding the
    gross range limits. If X percent of wavelength bins were flagged as fail (or fail/suspect) and exceed the
    user-defined fail threshold, then the entire spectrum is flagged as fail. If the percentage is between the suspect
    and fail thresholds then the spectra is flagged as suspect.
    The threshold represents the percentage threshold between 0-1 (0-100%).

    :param nd_gross_range_results: The N-dimensional results from the gross range test.
    :param wavelength_dim: The wavelength dimension of the results.
    :param suspect_threshold: The percentage threshold for suspect data.
    :param fail_threshold: The percentage threshold for fail data.
    :param include_suspect_flags: If True, include suspect flags along the spectra in the total fail.
    :return: A blanket flag for the spectra, which maintains the same shape as the time dimension.
    """
    if fail_threshold < suspect_threshold:
        raise ValueError('fail_threshold must be greater than suspect_threshold.')

    num_wvls = len(nd_gross_range_results[wavelength_dim])

    if include_suspect_flags is True:
        _flags_bool = xr.where((nd_gross_range_results == 4) | (nd_gross_range_results == 3), 1, 0)
    else:
        _flags_bool = xr.where(nd_gross_range_results == 4, 1, 0)

    _flags_bool_sum = np.sum(_flags_bool, axis=1)
    ratio = _flags_bool_sum / num_wvls
    flags = xr.full_like(nd_gross_range_results.time.astype(int), FLAG.NOT_EVALUATED).astype(int)

    flags = flags.where((ratio < suspect_threshold) & (ratio > fail_threshold), FLAG.SUSPECT)
    flags = flags.where(ratio < fail_threshold, FLAG.FAIL)

    flags = flags.where(ratio >= suspect_threshold, FLAG.PASS)

    # Assign attributes to the flags if an xarray.DataArray.
    flags.attrs['ancillary_variables'] = nd_gross_range_results.name
    flags.attrs['suspect_threshold'] = suspect_threshold
    flags.attrs['fail_threshold'] = fail_threshold
    flags.attrs['threshold_units'] = 'm^-1'
    flags.attrs['includes_suspect_flags'] = include_suspect_flags
    flags.attrs['test_name'] = 'blanket_gross_range_test'
    return flags


def a_gt_c_test(absorption: xr.DataArray, attenuation: xr.DataArray) -> xr.DataArray:
    """
    Assess if the absorption is greater than the attenuation. Having absorption greater than attenuation is (mostly)
    lawfully impossible. It is recommended that this test be run on scattering corrected absorption with ts-corrected
    attenuation. This test is not included in the QARTOD manual, but is a custom test based on reality checks in the
    ACS protocol document.

    :param absorption: Absorption data.
    :param attenuation: Attenuation data.
    :return: Flag indicating if absorption is greater than attenuation.
    """

    flags = xr.full_like(absorption, FLAG.NOT_EVALUATED).astype(int)
    flags = flags.where(absorption <= attenuation, FLAG.SUSPECT)
    flags = flags.where(absorption > attenuation, FLAG.PASS)

    # Assign attributes to the flags if an xarray.DataArray.
    flags.attrs['ancillary_variables'] = [absorption.name, attenuation.name]
    flags.attrs['test_name'] = 'a_gt_c_test'
    return flags


def rolling_variance_test(mts: xr.DataArray,
                          use_mean: str = 'rolling',
                          window_size: int = 4 * 60,
                          exceedance_percentage: float = 0.25,
                          min_periods: int | None = 1) -> xr.DataArray:
    """
    Apply a rolling variance test to a measured absorption or attenuationcoefficient.

    :param mts: The measured absorption or attenuation coefficient. For absorption the recommendation is to use
        a scattering corrected measurement. For attenuation the recommendation is to us a TS-corrected measurement.
    :param use_mean: Indicates whether to use the timeseries mean or a rolling window mean along the time dimension.
    :param window_size: The centered window size for the variance window
        and for the mean window if use_mean = 'rolling'.
    :param exceedance_percentage: The percentage of the mean that the variance must exceed to be flagged as suspect.
    :param min_periods: The minimum number of periods to perform this test.
        Same functionality as the min_periods argument for xarray.DataArray.rolling.
    :return: A flag indicating SUSPECT or PASS.
    """
    use_mean = use_mean.lower()
    if use_mean == 'timeseries':
        m_mean = mts.mean(dim='time', skipna=True)
        mean_type = 'timeseries_mean'
    elif use_mean == 'rolling':
        m_mean = mts.rolling(time=window_size, center=True, min_periods=min_periods).mean(skipna=True)
        mean_type = 'rolling_mean'
    else:
        raise ValueError("use_mean must be 'timeseries' or 'rolling'.")

    m_var = mts.rolling(time=window_size, center=True, min_periods=min_periods).var(skipna=True)

    flags = xr.full_like(mts, FLAG.NOT_EVALUATED).astype(int)
    flags = flags.where(m_var < exceedance_percentage * m_mean, FLAG.SUSPECT)
    flags = flags.where(m_var >= exceedance_percentage * m_mean, FLAG.PASS)

    # Assign attributes to the flags if an xarray.DataArray.
    flags.attrs['ancillary_variables'] = mts.name
    flags.attrs['rolling_variance_test_window_size'] = window_size
    flags.attrs['exceedance_percentage'] = exceedance_percentage
    flags.attrs['mean_type'] = mean_type
    flags.attrs['min_periods'] = min_periods
    flags.attrs['test_name'] = 'rolling_variance_test'
    return flags
