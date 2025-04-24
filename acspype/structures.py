from datetime import datetime
from typing import NamedTuple


class ACSPacket(NamedTuple):
    daq_time: datetime  # The time that the packet was parsed from the serial stream. The ACS does not provide a timestamp.
    packet_data: bytearray  # The raw byte files of the packet.


class ParsedPacket(NamedTuple):
    daq_time: datetime
    registration_bytes: bytes
    record_length: int
    packet_type: int
    reserved_1: int
    serial_number_int: int
    serial_number_hexdec: str
    serial_number: str
    a_reference_dark: int
    raw_pressure: int
    a_signal_dark: int
    raw_external_temperature: int
    raw_internal_temperature: int
    c_reference_dark: int
    c_signal_dark: int
    elapsed_time: int
    reserved_2: int
    number_of_output_wavelengths: int
    c_reference: tuple
    a_reference: tuple
    c_signal: tuple
    a_signal: tuple
    checksum: int
    flag_syntax: int
    flag_elapsed_time: int


class CalibratedPacket(NamedTuple):
    daq_time: datetime
    serial_number: str
    flag_syntax: int
    elapsed_time: int
    flag_elapsed_time: int
    internal_temperature: float
    flag_internal_temperature: int
    external_temperature: float
    a_wavelength: tuple
    c_wavelength: tuple
    a_uncorrected: tuple
    c_uncorrected: tuple
    discontinuity_wavelength_index: int
    a_discontinuity_offset: int
    c_discontinuity_offset: int
    a_m_discontinuity: tuple
    c_m_discontinuity: tuple
    a_m: tuple
    c_m: tuple