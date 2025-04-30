from datetime import datetime
import numpy as np
from numpy.typing import ArrayLike
from typing import NamedTuple
import xarray as xr


class ACSPacket(NamedTuple):
    daq_time: datetime
    full_packet: bytes

    def to_dict(self) -> dict:
        """Export the ACSPacket as a dictionary."""
        return self._asdict()

    def to_xarray(self) -> xr.Dataset:
        """Export the ACSPacket as an xr.Dataset."""
        ds = xr.Dataset()
        ds = ds.assign_coords({'daq_time': self.daq_time})
        ds['full_packet'] = bytes(self.full_packet)
        return ds


class ParsedPacket(NamedTuple):
    daq_time: datetime
    record_length: int
    packet_type: int
    # reserved_1: int  # Commented out for unlikely reimplementation in the future.
    sn_int: int
    a_reference_dark: int
    # raw_pressure: int
    a_signal_dark: int
    raw_external_temperature: int
    raw_internal_temperature: int
    c_reference_dark: int
    c_signal_dark: int
    elapsed_time: int
    # reserved_2: int
    number_of_wavelengths: int
    c_reference: tuple[int, ...]
    a_reference: tuple[int, ...]
    c_signal: tuple[int, ...]
    a_signal: tuple[int, ...]
    checksum: int

    def to_dict(self) -> dict:
        """Export the ParsedPacket as a dictionary."""
        return self._asdict()

    def to_xarray(self) -> xr.Dataset:
        """Export the ParsedPacket as an xr.Dataset."""
        ds = xr.Dataset()
        ds = ds.assign_coords({'daq_time': [self.daq_time],
                               'wavelength_index': list(range(self.number_of_wavelengths))})
        ds['record_length'] = self.record_length
        ds['packet_type'] = self.packet_type
        # ds['reserved_1'] = self.reserved_1
        ds['sn_int'] = self.sn_int
        ds['a_reference_dark'] = self.a_reference_dark
        # ds['raw_pressure'] = self.raw_pressure
        ds['a_signal_dark'] = self.a_signal_dark
        ds['raw_external_temperature'] = self.raw_external_temperature
        ds['raw_internal_temperature'] = self.raw_internal_temperature
        ds['c_reference_dark'] = self.c_reference_dark
        ds['c_signal_dark'] = self.c_signal_dark
        ds['elapsed_time'] = self.elapsed_time
        # ds['reserved_2'] = self.reserved_2
        ds['number_of_wavelengths'] = self.number_of_wavelengths
        ds['c_reference'] = (['daq_time', 'wavelength_index'], [self.c_reference])
        ds['a_reference'] = (['daq_time', 'wavelength_index'], [self.a_reference])
        ds['c_signal'] = (['daq_time', 'wavelength_index'], [self.c_signal])
        ds['a_signal'] = (['daq_time', 'wavelength_index'], [self.a_signal])
        ds['checksum'] = self.checksum
        return ds


class DeviceCalibratedPacket(NamedTuple):
    daq_time: datetime
    a_wavelength: tuple[float, ...]
    c_wavelength: tuple[float, ...]
    sn_hexdec: str
    sn_str: str
    internal_temperature: float
    external_temperature: float
    a_uncorrected: ArrayLike
    c_uncorrected: ArrayLike
    discontinuity_index: int
    a_discontinuity_offset: float
    c_discontinuity_offset: float
    a_m_discontinuity: ArrayLike
    c_m_discontinuity: ArrayLike
    a_m: ArrayLike
    c_m: ArrayLike

    def to_dict(self) -> dict:
        return self._asdict()

    def to_xarray(self) -> xr.Dataset:
        ds = xr.Dataset()
        ds = ds.assign_coords({'daq_time': [self.daq_time],
                               'a_wavelength': np.array(self.a_wavelength),
                               'c_wavelength': np.array(self.c_wavelength)})
        ds['sn_hexdec'] = self.sn_hexdec
        ds['sn_str'] = self.sn_str
        ds['internal_temperature'] = self.internal_temperature
        ds['external_temperature'] = self.external_temperature
        ds['a_uncorrected'] = (['daq_time', 'a_wavelength'], [self.a_uncorrected])
        ds['c_uncorrected'] = (['daq_time', 'c_wavelength'], [self.c_uncorrected])
        ds['discontinuity_index'] = self.discontinuity_index
        ds['a_discontinuity_offset'] = self.a_discontinuity_offset
        ds['c_discontinuity_offset'] = self.c_discontinuity_offset
        ds['a_m_discontinuity'] = (['daq_time', 'a_wavelength'], [self.a_m_discontinuity])
        ds['c_m_discontinuity'] = (['daq_time', 'c_wavelength'], [self.c_m_discontinuity])
        ds['a_m'] = (['daq_time', 'a_wavelength'], [self.a_m])
        ds['c_m'] = (['daq_time', 'c_wavelength'], [self.c_m])
        return ds


class TSCorrectedPacket(NamedTuple):
    daq_time: datetime
    a_wavelength: tuple[float, ...]
    c_wavelength: tuple[float, ...]
    temperature: float
    salinity: float
    a_mts: ArrayLike
    c_mts: ArrayLike


    def to_dict(self) -> dict:
        return self._asdict()


    def to_xarray(self) -> xr.Dataset:
        ds = xr.Dataset()
        ds = ds.assign_coords({'daq_time': self.daq_time,
                               'a_wavelength': np.array(self.a_wavelength),
                               'c_wavelength': np.array(self.c_wavelength)})
        ds['temperature'] = self.temperature
        ds['salinity'] = self.salinity
        ds['a_mts'] = (['daq_time', 'a_wavelength'], [self.a_mts])
        ds['c_mts'] = (['daq_time', 'c_wavelength'], [self.c_mts])
        return ds
