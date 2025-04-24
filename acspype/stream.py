from datetime import datetime, timezone
import serial

from acspype.core import PACKET_REGISTRATION, DefaultSerial
from acspype.structures import ACSPacket

class ACSStream():
    def __init__(self, port: str,
                 baudrate: int = DefaultSerial.BAUDRATE,
                 timeout: float = DefaultSerial.TIMEOUT) -> None:
        """
        Instantiate the ACSStream class, reset the serial port, and clear any data in memory.

        :param port: The COM port the ACS is connected to.
        :param baudrate: The baudrate of the ACS. All ACS sensors come from the factory set to 115200 bps.
        :param timeout: The number of seconds to wait before actions timeout.
        :return: None
        """

        self.clear_packet_buffer()
        self.reset_serial(port = port, baudrate = baudrate, timeout=timeout)


    def clear_serial_buffers(self) -> None:
        """
        Clear the input and output serial buffers on the class assigned port.
        :return: None
        """
        self._serial.reset_input_buffer()
        self._serial.reset_output_buffer()


    def clear_packet_buffer(self):
        """
        Clear the packet buffer that is currently in memory.
        :return: None
        """
        self._buffer = bytearray()


    def reset_serial(self, port: str,
                     baudrate: int,
                     bytesize: int = DefaultSerial.BYTESIZE,
                     parity: str = DefaultSerial.PARITY,
                     stopbits: int = DefaultSerial.STOPBITS,
                     flowcontrol: int = DefaultSerial.FLOWCONTROL,
                     timeout: int = DefaultSerial.TIMEOUT):
        """
        Reset or instantiate serial object using the assigned settings.

        :param port: A string representing the COM port the ACS is connected to. e.g. 'COM1' or '/dev/ttyUSB0'
        :param baudrate: The baudrate the ACS communicates at. All ACS sensors come from the factory set to 115200 bps.
        :param bytesize: The number of bits in a byte. Default is 8.
        :param parity: The parity setting for the ACS. Default is 'N' for no parity.
        :param stopbits: The number of stop bits. Default is 1.
        :param flowcontrol: The flow control setting. Default is 0 for no flow control.
        :param timeout: The number of seconds to wait before timeout. Default is 3 seconds.
        :return: None
        """

        self._serial = serial.Serial()
        self._serial.port = port
        self._serial.baudrate = int(baudrate)
        self._serial.bytesize = int(bytesize)
        self._serial.parity = parity
        self._serial.stopbits = int(stopbits)
        self._serial.xonxoff = int(flowcontrol)
        self._serial.timeout = int(timeout)


    def connect(self) -> None:
        """
        Attempt to connect to the port assigned at class instantiation.

        :return: None
        """

        try:
            self._serial.open()
        except ConnectionError:
            raise serial.PortNotOpenError()


    def disconnect(self) -> None:
        """
        Disconnect from the port assigned at class instantiation.

        :return: None
        """

        self.clear_serial_buffers()
        self._serial.close()
        self.clear_packet_buffer()


    def read_stream(self) -> None:
        """
        Read the current number of bytes waiting in the serial port buffer and append them to a byte array in memory.

        :return: None
        """

        incoming = self._serial.read(self._serial.in_waiting)
        self._buffer.extend(incoming)


    def find_packet(self) -> ACSPacket:
        """
        Find the first full ACS packet within buffered serial data.

        :return: Raw ACS Packet data stored in a custom NamedTuple.
        """

        old_bytes, reg_bytes, remaining_bytes = self._buffer.partition(PACKET_REGISTRATION)
        daq_time = datetime.now(timezone.utc)
        if PACKET_REGISTRATION in remaining_bytes:
            packet_data, next_reg_bytes, next_loop_bytes = remaining_bytes.partition(PACKET_REGISTRATION)
            packet = reg_bytes + packet_data  # Rebuild the complete packet.
            self._buffer = next_reg_bytes + next_loop_bytes  # Redefine buffer.
            acs_packet = ACSPacket(daq_time=daq_time, packet_data=packet)
        else:
            acs_packet = ACSPacket(daq_time=daq_time, packet_data=None)
        return acs_packet


    def __enter__(self) -> object:
        """
        Context Manager: Automatically connect to the assigned port number when used in a "with" statement.

        :return: Class conditions.
        """

        self.connect()
        return self


    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """
        Context Manager: Exit and perform garbage cleanup.

        :param exc_type: The type of exception that was raised.
        :param exc_value: The value of the exception that was raised.
        :param traceback: The traceback of the exception that was raised.
        :return: None
        """

        self.disconnect()