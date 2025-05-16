"""This module contains functions for performing checks on data manipulated through acspype."""

import numpy as np
from numpy.typing import ArrayLike


def _check_dimensions(dimensions: list) -> ValueError | None:
    """
    Check the dimensions of the xr.Dataset to ensure that 'time' is included.
    acspype relies on the 'time' dimensions to be present in the xr.Dataset.

    :param dimensions: A list of dimension for the xr.Dataset.
    :return: No errors if the dimension checks pass, otherwise raises a ValueError.
    """

    if "time" not in dimensions:
        raise ValueError("The xr.Dataset dimensions must include 'time'.")
    else:
        return None


def _compare_wavelengths(wavelength_set_1: ArrayLike[float],
                         wavelength_set_2: ArrayLike[float]) -> ValueError | None:
    """
    Compare two sets of wavelengths to ensure they are the same.

    :param wavelength_set_1: An array of wavelengths that represent wavelengths from an ACS.
    :param wavelength_set_2: An array of wavelengths that represent wavelengths from an ACS.
    :return: No errors if the wavelengths are the same, otherwise raises a ValueError.
    """

    if np.all(wavelength_set_1 != wavelength_set_2):
        raise ValueError("Mismatch between wavelength sets.")
    else:
        return None
