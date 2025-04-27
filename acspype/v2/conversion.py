import numpy as np
from numpy.typing import ArrayLike
from struct import unpack_from
from typing import Union, Tuple, NamedTuple
import xarray as xr

from acspype.v2.structures import SplitPacket

from acspype.core import LPR, PACKET_TAIL, PACKET_HEAD, PAD_BYTE, WVL_BYTE_OFFSET
from acspype.dev import ACSDev
from acspype.tscor import ACSTSCor
from acspype.structures import ACSPacket, ParsedPacket, CalibratedPacket
from acspype.qaqc import FLAG, internal_temperature_test, syntax_test, elapsed_time_test
from acspype.discontinuity import (find_discontinuity_index, _apply_discontinuity_offset,
                                   _compute_discontinuity_offset, discontinuity_correction)


def split_packet(full_packet: bytearray):
    """
    Split the full ACS packet into raw data.

    Note A: In the ACS Manual (Table 2), reserved_1, raw_pressure, and reserved_2 are listed as data record fields,
    but these outputs do not provide any useful information and are not used anywhere in acspype. They are provided in
    this function and the SplitPacket structure for completeness and the slight possibility
    that they are implemented in the future by SBS.

    Note B: a_reference, a_signal, c_reference, and c_signal are tuples of integers of variable length.
    The length is defined by the number of wavelengths and will vary between instruments, not within a single instrument.

    :param full_packet: A full packet from the ACS. Must include the registration bytes and the pad byte.
    :return:
    """
    [nwvls] = unpack_from('B', full_packet, offset=WVL_BYTE_OFFSET)  # Get the # of wavelengths from the packet header.
    remaining_packet_size = int(nwvls * 2 * 2) # a_signal, a_reference, c_signal, c_reference
    full_packet_descriptor = PACKET_HEAD + f"{remaining_packet_size}H" + PACKET_TAIL  # Build a full packet descriptor.
    pad_byte_idx = full_packet.rfind(PAD_BYTE) # Get the index of the pad byte so the checksum can be pulled out.
    checksum_bytes = full_packet[pad_byte_idx - 2:pad_byte_idx]
    checksum = int(unpack_from('!H', checksum_bytes))
    raw = unpack_from(full_packet_descriptor, full_packet)  # Note: The padbyte is not returned when unpacking.

    sp = SplitPacket(
        record_length = raw[LPR + 0],
        packet_type = raw[LPR + 1],
        # reserved_1 = raw[LPR + 2],
        sn_int = raw[LPR + 3],
        a_reference_dark = raw[LPR + 4],
        # raw_pressure = raw[LPR + 5],
        a_signal_dark = raw[LPR + 6],
        raw_external_temperature = raw[LPR + 7],
        raw_internal_temperature = raw[LPR + 8],
        c_reference_dark = raw[LPR + 9],
        c_signal_dark = raw[LPR + 10],
        elapsed_time = raw[LPR + 11],
        # reserved_2 = raw[LPR + 12],
        number_of_wavelengths = raw[LPR + 13],
        c_reference = raw[LPR + 14:LPR + 14 + (nwvls * 4):4],
        a_reference = raw[LPR + 15:LPR + 15 + (nwvls * 4):4],
        c_signal = raw[LPR + 16:LPR + 16 + (nwvls * 4):4],
        a_signal = raw[LPR + 17:LPR + 17 + (nwvls * 4):4],
        checksum = checksum
    )
    return sp

