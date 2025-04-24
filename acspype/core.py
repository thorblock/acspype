NUM_PAT = "[+-]?[0-9]*[.]?[0-9]+"  # REGEX for any number, float or int, positive or negative.

PACKET_REGISTRATION = b'\xff\x00\xff\x00'  # Start of every ACS packet.
PAD_BYTE = b'\x00'  # End of every ACS packet.
WVL_BYTE_OFFSET = 4 + 2 + 1 + 1 + 1 + 3 + 2 + 2 + 2 + 2 + 2 + 2 + 2 + 4 + 1  # See Process Data section in ACS manual.
NUM_CHECKSUM_BYTES = 2
PACKET_HEAD = '!4cHBBl7HIBB'  # struct descriptor for the static header of a packet.
PACKET_TAIL = 'Hx' # struct descriptor for the static tail of a packet.
LPR = len(PACKET_REGISTRATION)

class DefaultSerial:
    BAUDRATE: int = 115200
    BYTESIZE: int = 8
    PARITY: str = 'N'
    STOPBITS: int = 1
    FLOWCONTROL: int = 0
    TIMEOUT: int = 3

# Raw pressure counts are no longer output by an ACS and can be safely ignored. The reserved_1 and reserved_2 variables are single byte variables that are not used by the ACS and can be ignored.
ACS_VARS_TO_IGNORE = ['raw_pressure', 'reserved_1', 'reserved_2']

#---------- File Creation ----------#
ENCODING = {'time': {'units': 'nanoseconds since 1900-01-01'}}  # xr.Dataset to netcdf encoding for time

#---------- PHYSICAL QUANTITIES ----------#
EST_FLOW_CELL_VOLUME = 30 # in mL, from the ACS Protocol Document, Rev Q.