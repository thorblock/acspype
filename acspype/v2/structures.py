from datetime import datetime
from typing import NamedTuple

class ACSPacket:
    """
    A NamedTuple class for holding raw ACS packet data.

    Note A: The daq_time field is a datetime object, in UTC, that represents the host computer time that the packet was
    parsed from the serial buffer. The ACS does not produce a clock time.

    Note B: The full_packet must be a full packet from the ACS, including the registration bytes in the beginning
    and the pad byte at the end.
    """
    daq_time: datetime
    full_packet: bytearray

class SplitPacket(NamedTuple):
    """
    A NamedTuple class for holding parsed ACS packet data.
    Note A: In the ACS Manual (Table 2), reserved_1, raw_pressure, and reserved_2 are listed as data record fields,
    but these outputs do not provide any useful information and are not used anywhere in acspype. They are provided
    in this class and the split_packet() function for completeness and the slight possibility that they are implemented
    in the future.

    Note B: a_reference, a_signal, c_reference, and c_signal are tuples of integers of variable length. The length
    is defined by the number of wavelengths.
    """

    record_length: int
    packet_type: int
    # reserved_1: int
    sn_int: int  # This value is a combination of the meter type and serial number.
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
