"""
This module contains NamedTuple structures that ACS data can be dumped into for easy tracking and exporting.
Each NamedTuple is assigned methods that allows for the export of data to a dictionary or to a formatted
xarray.Dataset.
"""

from datetime import datetime
import numpy as np
from numpy.typing import ArrayLike
from typing import NamedTuple
import xarray as xr


class ACSPacket(NamedTuple):
    """
    A class for representing a time-stamped ACS packet.

    :param daq_time: The host computer time the packet was acquired.
    :param full_packet: The full packet, including the registration bytes and the checksum + pad bytes.
    """
    daq_time: datetime
    full_packet: bytes | bytearray | None

    def to_dict(self) -> dict:
        """
        Export the ACSPacket as a dictionary.

        :return: A dict representation of the ACSPacket, where keys are attributes and values are attribute values.
        """
        return self._asdict()

    def to_xarray(self) -> xr.Dataset:
        """
        Export the ACSPacket as an xarray.Dataset.

        :return: An xarray.Dataset representation of the ACSPacket, where coordinates are the daq_time
            and the data variable is the full_packet.
        """
        ds = xr.Dataset()
        ds = ds.assign_coords({'daq_time': self.daq_time})
        ds['full_packet'] = bytes(self.full_packet)
        return ds


class ParsedPacket(NamedTuple):
    """
    A class for representing a parsed ACS packet.

    :param daq_time: The host computer time the packet was acquired, usually passed from an ACSPacket.
    :param record_length: The length of the packet in bytes.
    :param packet_type: The type of packet.
    :param sn_int: The serial number of the device represented as an integer.
    :param a_reference_dark: The dark reference counts for the absorption channel.
    :param a_signal_dark: The dark signal counts for the absorption channel.
    :param raw_external_temperature: The raw external temperature reading in counts.
    :param raw_internal_temperature: The raw internal temperature reading in counts.
    :param c_reference_dark: The dark reference counts for the attenuation channel.
    :param c_signal_dark: The dark signal counts for the attenuation channel.
    :param elapsed_time: The elapsed time since power was supplied to the ACS, in milliseconds.
    :param number_of_wavelengths: The number of wavelengths for each channel.
    :param c_reference: The attenuation reference in counts.
    :param a_reference: The absorption reference in counts.
    :param c_signal: The attenuation signal in counts.
    :param a_signal: The absorption signal in counts.
    :param checksum: The checksum of the packet.
    """
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
        """
        Export the ParsedPacket as a dictionary.

        :return: A dict representation of the ParsedPacket, where keys are attributes and values are attribute values.
        """
        return self._asdict()

    def to_xarray(self) -> xr.Dataset:
        """
        Export the ParsedPacket as an xr.Dataset.

        :return: An xarray.Dataset representation of the ParsedPacket.
        """
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
    """
    A class for representing a device file-calibrated ACS packet.

    :param daq_time: The host computer time the packet was acquired, usually passed from an ParsedPacket.
    :param a_wavelength: The wavelengths for the absorption channel.
    :param c_wavelength: The wavelengths for the attenuation channel.
    :param sn_hexdec: The serial number of the device represented as a hexadecimal string.
    :param serial_number: The serial number of the device represented as a string.
    :param internal_temperature: The internal temperature of the device in degrees Celsius.
    :param external_temperature: The external temperature of the device in degrees Celsius.
    :param a_uncorrected: The uncorrected absorption signal in inverse meters.
    :param c_uncorrected: The uncorrected attenuation signal in inverse meters.
    :param discontinuity_index: The index of the discontinuity in both channels.
    :param a_discontinuity_offset: The offset for the absorption channel discontinuity.
    :param c_discontinuity_offset: The offset for the attenuation channel discontinuity.
    :param a_m_discontinuity: The absorption coefficient without discontinuity correction applied.
    :param c_m_discontinuity: The attenuation coefficient without discontinuity correction applied.
    :param a_m: The measured absorption coefficient in inverse meters.
    :param c_m: The measured attenuation coefficient in inverse meters.
    """

    daq_time: datetime
    a_wavelength: tuple[float, ...]
    c_wavelength: tuple[float, ...]
    sn_hexdec: str
    serial_number: str
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
        """
        Export the DeviceCalibratedPacket as a dictionary.

        :return: A dict representation of the DeviceCalibratedPacket, where keys are attrs and values are attr values.
        """
        return self._asdict()

    def to_xarray(self) -> xr.Dataset:
        """
        Export the DeviceCalibratedPacket as an xarray.Dataset.

        :return: An xarray.Dataset representation of the DeviceCalibratedPacket.
        """
        ds = xr.Dataset()
        ds = ds.assign_coords({'daq_time': [self.daq_time],
                               'a_wavelength': np.array(self.a_wavelength),
                               'c_wavelength': np.array(self.c_wavelength)})
        ds['sn_hexdec'] = self.sn_hexdec
        ds['serial_number'] = self.serial_number
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
    """
    A class for representing a temperature and salinity-corrected ACS packet.

    :param daq_time: The host computer time the packet was acquired, usually passed from a DeviceCalibratedPacket.
    :param a_wavelength: The wavelengths for the absorption channel.
    :param c_wavelength: The wavelengths for the attenuation channel.
    :param temperature: The temperature of the water in degrees Celsius.
    :param salinity: The salinity of the water in practical salinity units (PSU).
    :param a_mts: The temperature and salinity-corrected absorption coefficient in inverse meters.
    :param c_mts: The temperature and salinity-corrected attenuation coefficient in inverse meters.
    """
    daq_time: datetime
    a_wavelength: tuple[float, ...]
    c_wavelength: tuple[float, ...]
    temperature: float
    salinity: float
    a_mts: ArrayLike
    c_mts: ArrayLike

    def to_dict(self) -> dict:
        """
        Export the TSCorrectedPacket as a dictionary.

        :return: A dict representation of the TSCorrectedPacket, where keys are attributes and values are attr values.
        """
        return self._asdict()

    def to_xarray(self) -> xr.Dataset:
        """
        Export the TSCorrectedPacket as an xarray.Dataset.

        :return: An xarray.Dataset representation of the TSCorrectedPacket.
        """
        ds = xr.Dataset()
        ds = ds.assign_coords({'daq_time': self.daq_time,
                               'a_wavelength': np.array(self.a_wavelength),
                               'c_wavelength': np.array(self.c_wavelength)})
        ds['temperature'] = self.temperature
        ds['salinity'] = self.salinity
        ds['a_mts'] = (['daq_time', 'a_wavelength'], [self.a_mts])
        ds['c_mts'] = (['daq_time', 'c_wavelength'], [self.c_mts])
        return ds
