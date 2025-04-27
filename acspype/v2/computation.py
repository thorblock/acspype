import numpy as np
from numpy.typing import ArrayLike, NDArray
from struct import unpack_from
from typing import Union, Tuple
import xarray as xr

from acspype.dev import ACSDev


def convert_sn_hexdec(sn_int: int) -> str:
    sn_hexdec = hex(sn_int)
    return sn_hexdec

def convert_sn_str(sn_int):
    sn_hexdec = convert_sn_hexdec(sn_int)
    sn_str = 'ACS-' + str(int(sn_hexdec[-6:], 16)).zfill(5)
    return sn_str


def compute_internal_temperature(counts: xr.DataArray) -> xr.DataArray:

    a = 0.00093135
    b = 0.000221631
    c = 0.000000125741
    d = 273.15
    volts = 5 * counts / 65535
    resistance = 10000 * volts / (4.516 - volts)
    internal_temperature = 1 / (a + b * np.log(resistance) + c * np.log(resistance) ** 3) - d

    internal_temperature.attrs['ancillary_variables'] = counts.name
    return internal_temperature


def compute_external_temperature(counts: xr.DataArray) -> xr.DataArray:

    a = -7.1023317e-13
    b = 7.09341920e-08
    c = -3.87065673e-03
    d = 95.8241397
    external_temperature = a * counts ** 3 + b * counts ** 2 + c * counts + d

    external_temperature.attrs['ancillary_variables'] = counts.name
    return external_temperature


def compute_uncorrected(signal_counts: xr.DataArray,
                        reference_counts: xr.DataArray,
                        dev: ACSDev) -> xr.DataArray:

    uncorr = (1 / dev.path_length) * np.log(signal_counts / reference_counts)

    uncorr.attrs['ancillary_variables'] = ', '.join([signal_counts.name, reference_counts.name])
    uncorr.attrs['path_length'] = dev.path_length
    return uncorr


def compute_measured(uncorrected: xr.DataArray,
                     a_or_c: str,
                     internal_temperature: float | ArrayLike,
                     dev: ACSDev) -> xr.DataArray:

    if a_or_c == 'a':
        delta_t = dev.func_a_delta_t(internal_temperature).T
        offsets = dev.a_offset
    elif a_or_c == 'c':
        delta_t = dev.func_c_delta_t(internal_temperature).T
        offsets = dev.c_offset
    measured = (offsets - uncorrected) - delta_t

    measured.attrs['ancillary_variables'] = ', '.join([uncorrected.name, internal_temperature.name, offsets.name])
    measured.attrs['delta_t_interpolation'] = dev.delta_t_interp_method
    return measured




