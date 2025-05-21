"""
This module contains functions that are core to utilities that are beneficial to use alongside acspype applications,
but do not support acspype directly.
"""


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
    """
    Iterate through available serial ports and check the incoming data for the ACS registration bytes.
    The first port with ACS registration bytes is returned. This function has not been tested with multiple ACS sensors
    connected to the same computer.

    :param baudrate: The baudrate for the ACS connection. Default is 115200 bps and does not need to be changed.
    :param timeout: The timeout for the serial connection. Default is 1 second.
    :param check_length: The amount of time in seconds to collect data from each serial port in
        hopes of receiving the ACS registration bytes.
    :return: The operating system serial port as a string. Can be used by pyserial for future connections.
    """

    available_ports = list_available_ports()
    for port in available_ports:
        try:
            with serial.Serial(port=port, baudrate=baudrate, timeout=timeout) as ser:
                time.sleep(check_length)
                incoming = ser.read(ser.in_waiting)
            if PACKET_REGISTRATION in incoming:
                return port
        except:  # Bare exceptions don't follow PEP 8: E772, but do we really care for this use case?
            continue
    raise ConnectionAbortedError('No ACS detected. Is the sensor connected and on?')
