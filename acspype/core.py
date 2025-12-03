"""This module contains core constants and functions used across acspype."""

NUM_PAT = r"[+-]?[0-9]*[.]?[0-9]+"  # REGEX for any number, float or int, positive or negative.

# ---------- PACKET HANDLING ---------- #
PACKET_REGISTRATION = b'\xff\x00\xff\x00'  # Start of every ACS packet.
PAD_BYTE = b'\x00'  # End of every ACS packet.
WVL_BYTE_OFFSET = 4 + 2 + 1 + 1 + 1 + 3 + 2 + 2 + 2 + 2 + 2 + 2 + 2 + 4 + 1  # See Process Data section in ACS manual.
NUM_CHECKSUM_BYTES = 2 # Number of bytes in the checksum.
PACKET_HEAD = '!4cHBBl7HIBB'  # struct descriptor for the first 31 bytes of an ACS packet.
PACKET_TAIL = 'Hx'  # struct descriptor for the last 3 bytes of an ACS packet (checksum + pad byte).
LPR = len(PACKET_REGISTRATION)

# ---------- FILE CREATION ---------- #
NC_ENCODING = {'time': {'units': 'milliseconds since 1970-01-01'}}  # Recommended xr.Dataset to netcdf encoding for time.

# ---------- PHYSICAL QUANTITIES ---------- #
EST_FLOW_CELL_VOLUME = 30  # in mL, from the ACS Protocol Document, Rev Q.
