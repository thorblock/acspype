import numpy as np
from numpy.typing import ArrayLike
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


def ts_correction(m: ArrayLike,
                  a_or_c: str,
                  wavelength_dim: str,
                  temperature: float | ArrayLike,
                  salinity: float | ArrayLike,
                  dev: ACSDev,
                  tscor: ACSTSCor) -> ArrayLike:
    """
    Compute temperature-salinity corrected data.

    :param m: The measured signal, corrected for pure water offsets. Does not yet imply the filtration state of the water.
    :param a_or_c: 'a' or 'c', used for obtaining the correct correction coeffs.
    :param wavelength_dim: The corresponding wavelength dimension for the data.
    :param temperature: sea water temperature. conservative temperature is acceptable per the ACS manual,
        but does not result in much change
    :param salinity: sea water salinity. absolute salinity is acceptable per the ACS manual,
        but does not result in much change.
    :param dev: The corresponding ACSDev object for the dataset.
    :param tscor: The corresponding ACSTSCor object for the dataset.

    :return: An ArrayLike object of the same type as the m input type.
    """

    if a_or_c == 'a':
        _tscor = tscor.to_xarray().sel(wavelength=m[wavelength_dim], method='nearest')
        delta_t, psi_temp = np.meshgrid(temperature - dev.tcal, _tscor.psi_t.values)
        s, psi_sal = np.meshgrid(salinity, _tscor.psi_s_a.values)
    elif a_or_c == 'c':
        _tscor = tscor.to_xarray().sel(wavelength=m[wavelength_dim], method='nearest')
        delta_t, psi_temp = np.meshgrid(temperature - dev.tcal, _tscor.psi_t.values)
        s, psi_sal = np.meshgrid(salinity, _tscor.psi_s_c.values)
    mts = m - ((psi_temp.T * delta_t.T) + (psi_sal.T * s.T))
    if isinstance(m, list):
        mts = mts.tolist()
    elif isinstance(m, tuple):
        mts = tuple(map(tuple, mts))
    return mts


def zero_shift_correction(mts: ArrayLike) -> ArrayLike:
    """
    Zero out any values that are between -0.005 and 0. According to the manual, these values are equivalent to 0.

    :param mts: The absorption or attenuation, corrected for temperature and salinity.
    :return: The zeroed out data.
    """

    if isinstance(mts, xr.DataArray):
        mts = mts.where((mts > 0) | (mts <= -0.005), 0)
    elif isinstance(mts, np.array):
        mts = np.where((mts < 0) & (mts >= -0.005), 0, mts)
    elif isinstance(mts, list):
        mts = np.array(mts)
        mts = np.where((mts < 0) & (mts >= -0.005), 0, mts)
        mts = mts.tolist()
    elif isinstance(mts, tuple):
        mts = np.array(mts)
        mts = np.where((mts < 0) & (mts >= -0.005), 0, mts)
        mts = tuple(map(tuple, mts))
    return mts


def baseline_scattering_correction(a_mts: ArrayLike, reference_wavelength: float = 715.0) -> ArrayLike:
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
