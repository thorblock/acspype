"""
This module contains functions that are beneficial to use alongside acspype applications,
but do not support acspype directly.
"""

from numpy.typing import NDArray
import serial
import serial.tools.list_ports
import time
import uncertainties.core

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


def is_uncertainties_object(values: NDArray, max_checks: int = 5) -> bool:
    """
    Check if an input value is an object from the uncertainties package.

    :param values: Any object.
    :param max_checks: The maximum number of times to check if an object is from the uncertainties package.
    :return: If the input is an uncertainties object, return True
        uncertainties supports arrays and matrices. This function will try to whittle down
        an object, if it is an array, until a single element can be tested. If
    """
    if isinstance(values, uncertainties.core.Variable):
        return True
    elif isinstance(values, uncertainties.unumpy.core.matrix):
        return True
    else:
        for i in range(max_checks):
            try:
                values = values[0]
            except:
                break
        if isinstance(values, uncertainties.core.Variable):
            return True
        else:
            return False
