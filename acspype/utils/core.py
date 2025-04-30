import serial
import serial.tools.list_ports
import time

from acspype.core import PACKET_REGISTRATION


def list_available_ports() -> list:
    """
    List available serial ports for the native operating system.

    :return: A list of available serial port names.
    """

    available_ports = [v.name for v in serial.tools.list_ports.comports()]
    return available_ports


def find_acs_port(baudrate: int = 115200,
                  timeout: int = 1,
                  check_length: int = 1) -> str:
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
    raise ConnectionAbortedError('No ACS detected. Is the sensor connected and on?')

