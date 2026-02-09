"""This module contains core constants and functions used across acspype."""

from math import nan, sqrt
from numpy.typing import NDArray
import uncertainties
import xarray as xr

NUM_PAT = r"[+-]?[0-9]*[.]?[0-9]+"  # REGEX for any number, float or int, positive or negative.

# ---------- PACKET HANDLING ---------- #
PACKET_REGISTRATION = b"\xff\x00\xff\x00"  # Start of every ACS packet.
PAD_BYTE = b"\x00"  # End of every ACS packet.
WVL_BYTE_OFFSET = 4 + 2 + 1 + 1 + 1 + 3 + 2 + 2 + 2 + 2 + 2 + 2 + 2 + 4 + 1
NUM_CHECKSUM_BYTES = 2  # Number of bytes in the checksum.
PACKET_HEAD = "!4cHBBl7HIBB"
PACKET_TAIL = "Hx"  # struct descriptor for the last 3 bytes of an ACS packet (checksum + pad byte).
LPR = len(PACKET_REGISTRATION)

# ---------- FILE CREATION ---------- #
NC_ENCODING = {
    "time": {"units": "milliseconds since 1970-01-01"}
}  # Recommended xr.Dataset to netcdf encoding for time.

# ---------- PHYSICAL QUANTITIES ---------- #
EST_FLOW_CELL_VOLUME = 30  # in mL, from the ACS Protocol Document, Rev Q.


# ---------- INSTRUMENT SPECIFICATIONS ---------- #
class SPECS:
    class ACS:
        ACCURACY: float = 0.01
        PRECISION_SHORT: float = 0.012  # Max between 400-449nm
        PRECISION_LONG: float = 0.003  # Max between 450-730nm
        UNC_SHORT: float = sqrt(
            (0.012**2) + (0.01**2)
        )  # Combined uncertainty from base precision/accuracy following Eq. 3b from Csavina et al. 2017.
        UNC_LONG: float = sqrt(
            (0.003**2) + (0.01**2)
        )  # Combined uncertainty from base precision/accuracy following Eq. 3b from Csavina et al. 2017.
        GROSS_MIN: float = 0.001
        GROSS_MAX: float = 10.0

    class INTERNAL_TEMPERATURE:
        ACCURACY: float = 0.1  # From email with Sea-Bird Tech Support.
        UNC: float = 0.1

    class EXTERNAL_TEMPERATURE:
        ACCURACY: float = nan
        UNC: float = nan

    class PRESSURE:
        ACCURACY: float = nan  # Unknown
        UNC: float = nan


def is_uncertainties(values: NDArray, max_checks: int = 5) -> bool:
    """
    Check if an object might be from the uncertainties package.
    :param values: The input to check.
    :param max_checks: The number of times to check in case the object is a nested NDArray.
    :return: A boolean indicating whether or not the object is from the uncertainties package.
    """
    if isinstance(values, xr.DataArray):
        values = values.values
    if isinstance(values, uncertainties.core.Variable) or isinstance(values, uncertainties.unumpy.core.matrix):
        return True
    else:
        for i in range(max_checks):
            try:
                values = values[0]
            except:
                break
        if isinstance(values, uncertainties.core.Variable):
            return True
        else:
            return False
