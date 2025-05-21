"""This module contains functions for performing checks on xarray.Datasets that contain ACS data."""

import numpy as np
from numpy.typing import ArrayLike


def _check_dimensions(dimensions: list) -> UserWarning | None:
    """
    Check the dimensions of the xr.Dataset to ensure that 'time' is included.
    acspype relies on the 'time' dimensions to be present in the xr.Dataset.

    :param dimensions: A list of dimension for the xr.Dataset.
    :return: No errors if the dimension checks pass, otherwise raises a ValueError.
    """

    if "time" not in dimensions:
        UserWarning("It is strongly recommended to use 'time' as the time dimension for an ACS dataset.")
    if "a_wavelength" not in dimensions and "wavelength" not in dimensions:
        UserWarning("It is strongly recommended to use 'a_wavelength' as the absorption wavelength dimension for an "
                    "ACS dataset.")
    elif "c_wavelength" not in dimensions and "wavelength" not in dimensions:
        UserWarning("It is strongly recommended to use 'c_wavelength' as the attenuation wavelength dimension for an "
                    "ACS dataset.")
    elif "wavelength" not in dimensions and "a_wavelength" not in dimensions and "c_wavelength" not in dimensions:
        UserWarning("It is strongly recommended to use 'wavelength' as the the common wavelength dimension for "
                    "the absorption and attenuation coefficients for an ACS dataset.")
    else:
        return None


def _compare_wavelengths(wavelength_bins_1: ArrayLike[float],
                         wavelength_bins_2: ArrayLike[float]) -> ValueError | None:
    """
    Compare two sets of wavelengths to ensure they are the same.

    :param wavelength_bins_1: An array of wavelengths that represent wavelengths from an ACS.
    :param wavelength_bins_2: An array of wavelengths that represent wavelengths from an ACS.
    :return: No errors if the wavelengths are the same, otherwise raises a ValueError.
    """

    if np.all(wavelength_bins_1 != wavelength_bins_2):
        raise ValueError("Mismatch between wavelength sets.")
    else:
        return None
