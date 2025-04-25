import serial
import serial.tools.list_ports
import time

from acspype.core import PACKET_REGISTRATION, DefaultSerial


def list_available_ports() -> list:
    """
    List available serial ports for the native operating system.

    :return: A list of available serial port names.
    """

    available_ports = [v.name for v in serial.tools.list_ports.comports()]
    return available_ports


def find_acs_port(baudrate: int = DefaultSerial.BAUDRATE,
                  timeout: float = DefaultSerial.TIMEOUT,
                  check_length: int = 1) -> str:
    """
    Find the serial port where an ACS is connected. Conceptually, this should always connect to the first serial port
    that is avaialble and hosts an ACS.

    :param baudrate: The baudrate of the ACS. The default is 115200 bps.
    :param timeout: The number of seconds to wait for a response from the ACS. The default is 1 second.
    :param check_length: The number of seconds to wait before reading the serial port buffer. The default is 1 second.
    :return: The port of the ACS, if found. Note that this does not keep the port open.
    """

    available_ports = list_available_ports()
    for port in available_ports:
        try:
            with serial.Serial(port = port, baudrate=baudrate, timeout=timeout) as ser:
                time.sleep(check_length)
                incoming = ser.read(ser.in_waiting)
            if PACKET_REGISTRATION in incoming:
                return port
        except:
            continue
    return ConnectionAbortedError('No ACS detected. Is the sensor connected and on?')

