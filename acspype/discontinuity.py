import numpy as np
from numpy.typing import ArrayLike
from scipy.interpolate import CubicSpline
from typing import Union, Tuple
import xarray as xr


def find_discontinuity_index(a_wavelengths: ArrayLike,
                             c_wavelengths: ArrayLike,
                             min_band: float = 535.0,
                             max_band: float = 600.0) -> int:
    """

    This code is modified from the OPTAA processing utilities in the ooi-data-explorations repo.
    https://github.com/IanTBlack/ooi-data-explorations/blob/master/python/ooi_data_explorations/uncabled/utilities/utilities_optaa.py#L104

    Find the last wavelength index of the first filter based on wavelength differences.
    This function assumes that the discontinuity occurs between 535 nm and 600 nm, which is buffered from the values in the ACS Protocol Document, Rev Q.

    :param a_wavelengths: Absorption wavelengths
    :param c_wavelengths: Attenuation wavelengths
    :return: The last wavelength index at the discontinuity.
    """

    # Copy and convert to numpy arrays because we are paranoid about global variables.
    a_wavelengths = np.array(a_wavelengths).copy()
    c_wavelengths = np.array(c_wavelengths).copy()

    # Set values outside the range to NaN
    a_wavelengths[(a_wavelengths < min_band) | (a_wavelengths > max_band)] = np.nan
    c_wavelengths[(c_wavelengths < min_band) | (c_wavelengths > max_band)] = np.nan

    # Find the index of the discontinuity
    didx = int(np.nanargmin(np.diff(a_wavelengths) + np.diff(c_wavelengths)))
    return didx



def _compute_discontinuity_offset(values: ArrayLike,
                                  wavelength: ArrayLike,
                                  didx: int) -> float:
    """
    This code is modified from the OPTAA processing utilities in the ooi-data-explorations repo.
    https://github.com/IanTBlack/ooi-data-explorations/blob/master/python/ooi_data_explorations/uncabled/utilities/utilities_optaa.py#L212

    Compute the scalar discontinuity offset to be applied to the second half of an ACS spectra.
    NOTE: If the input values contain an inf or NaN value, this function will return a NaN.
    This is to prevent math errors associated with creating a cubic spline on infinite values.
    Spectra with infinite values should be removed at some point in the processing pipeline.

    :param values: The incoming absorption or attenuation values. It is highly recommended that these values be
        representative of ACS data that have been converted to 'geophysical' units (1/m) and corrected for the effects
        of internal temperature on output. That is to say, the recommended input is the measured (a_m and c_m) in the
        ACS protocol documents and manual. acspype inputs would be a_m_discontinuity and c_m_discontinuity.
    :param wavelength: The wavelength bins of the values.
    :param didx: The index of discontinuity.
    :return: The discontinuity offset for the second half of the spectrum.
    """
    _wavelength = np.copy(wavelength)
    _values = np.copy(values)
    _didx = int(np.copy(didx))

    x = _wavelength[_didx - 2:_didx + 1]
    y = _values[_didx - 2:_didx + 1]

    if np.any(np.isinf(y)) or np.any(np.isnan(y)):
        return np.nan
    else:
        cubic_spline = CubicSpline(x, y)
        interp = cubic_spline(_wavelength[_didx + 1], extrapolate=True)
        offset = interp - _values[_didx + 1]
        return offset


def _apply_discontinuity_offset(values: ArrayLike,
                                offset: float,
                                didx: int) -> np.array:

    """
    This code is modified from the OPTAA processing utilities in the ooi-data-explorations repo.
    https://github.com/IanTBlack/ooi-data-explorations/blob/master/python/ooi_data_explorations/uncabled/utilities/utilities_optaa.py#L212

    Apply a pre-determined discontinuity offset to values after the discontinuity index.

    :param values: The measured values to apply the discontinuity offset to.
    :param offset: The scalar offset to apply to the values after the discontinuity index.
    :param didx: The discontinuity index computed from find_discontinuity_index.
    :return: A discontinuity-corrected spectra.
    """

    _values = np.copy(values)
    _offset = np.copy(offset)
    _didx = int(np.copy(didx))
    _values[_didx + 1:] = _values[_didx + 1:] + _offset
    return _values


def compute_discontinuity_offset(measured: xr.DataArray,
                                 wavelength_dim: str,
                                 discontinuity_index: int) -> xr.DataArray:
    """
    This is a wrapper function for _compute_discontinuity_offset that is vectorized for Xarray.

    :param measured: The measured values
    :param wavelength_dim: The dimension to calculate the offset on.
    :param discontinuity_index
    :return: The scalar offset for a given spectrum.
    """

    offset = xr.apply_ufunc(_compute_discontinuity_offset, measured,
                            kwargs = {'wavelength': measured[wavelength_dim].values, 'didx': discontinuity_index},
                            input_core_dims = [[wavelength_dim]],
                            output_core_dims = [[]],
                            vectorize = True)
    return offset


def apply_discontinuity_offset(measured: xr.DataArray,
                               offset: xr.DataArray,
                               wavelength_dim: str,
                               discontinuity_index: int) -> xr.DataArray:
    """
    This is a wrapper function for _apply_discontinuity_offset that is vectorized for Xarray.
    :param measured: The measured values.
    :param offset: The pre-determined discontinuity offset.
    :param wavelength_dim: The wavelength dimension to apply the offset to.
    :param discontinuity_index: The index of discontinuity.
    :return: Spectra with the discontinuity offset applied.
    """

    dc = xr.apply_ufunc(_apply_discontinuity_offset, measured, offset,
                        kwargs = {'didx': discontinuity_index},
                        input_core_dims=[[wavelength_dim],[]],
                        output_core_dims = [[wavelength_dim]],
                        vectorize = True)
    return dc


def discontinuity_correction(measured: xr.DataArray,
                             wavelength_dim: str,
                             discontinuity_index: int) -> Tuple[xr.DataArray, xr.DataArray]:
    """
    This is a convenience function for computing the discontinuity offset and applying it to the measured values in Xarray.

    :param measured: The measured values.
    :param wavelength_dim: The wavelength dimension to apply the correction to.
    :param discontinuity_index: The index of discontinuity.
    :return:
    """

    offset = compute_discontinuity_offset(measured, wavelength_dim, discontinuity_index)
    dc = apply_discontinuity_offset(measured, offset, wavelength_dim, discontinuity_index)
    return (dc, offset)