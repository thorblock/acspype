# acpype

acpype provides functions for reading [Sea-Bird Scientific ACS](https://www.seabird.com/ac-s-spectral-absorption-and-attenuation-sensor/product?id=60762467715) data over serial and for performing community-accepted processing on ACS data.

ACS data are inherently complex and difficult to work with, particularly for new users without strong optics backgrounds.
This package attempts to simplify the process of reading and processing ACS data so that users can more quickly get to the more advanced data products for their research.

# CAUTION
If you are using this package to acquire data from an ACS over serial, it is your responsibility to ensure that you are 
using the appropriate equipment. Please contact Sea-Bird Scientific for more information about the relevant equipment and
conditions needed to acquire data from the ACS.


## What does this package need from you (and Sea-Bird Scientific)?
In addition to the physical sensor, one who deploys an ACS will also receive a device file (.dev) and a temperature-salinity correction file (.cor) from Sea-Bird Scientific.

As of publishing of this software package, the current version of the device file received from the factory is version 3. The device file is a text file that contains factory derived pure water offsets.
These offsets are created at select temperature bins, which acpype uses to correct data output for temperature changes within the sensor, which impact the internal optics.
The file provides pure water offsets. Removal of these offsets from a datasets then removes the effect of water on the absorption and attenuation measurements.

The TS4.cor file that is delivered with each ACS provides empirically derived temperature and salinity correction coefficients derived by Sullivan et al. 2006.
These coefficients are used to correct the absorption and attenuation measurements for temperature and salinity effects.

## What does this package not do?
This package does not provide any functionality for logging data from the ACS. The `serial_no_device_file.ipynb` and `serial_with_device_file.ipynb` offer acypype-specific data structures which can then be converted to a format that enables the logging of the ACS data.
Inherently, the custom data structures are NamedTuple objects. To access that information as a Python dictionary, you can use the `._asdict()` method.

In applications where concurrency is required, PostgreSQL is a good option for storing data, since it provides functionality for handling arrays.
If file-based logging is required and single files are desired, SQLite is a good option. The ACS produces a significant amount of data, so logging to text or .netCDF files may best be done as hourly or daily files to prevent excessive memory usage.



## Restrictions
The data streaming and processing functions are intended to be used separately from one another and not as a single pipeline.
The ACS outputs data at 4Hz and processing a serial packet to scattering corrected values is something that can't be done within 250ms for most computers.
Conceptually, this could be resolved by passing data between threads, but for now we will leave that up to the end user to implement.

acpype does not typically enforce a strict naming convention of variables created  using the software, but the examples provided use terminology commonly found in the manual and in literature that uses ACS data.
For obtaining data over serial, the provided functions create predefined names, which users can change out at will.

The only enforced naming conventions are with the dimensions of any ACS data or metadata file imported with Xarray. 
The following dimensions are required for ACS datasets:

- time: The time of the ACS sample, in UTC. Derived from the computer that is reading data over serial.
- a_wavelength: The absorption wavelength bins derived from the device file.
- c_wavelength: The attenuation wavelength bins derived from the device file.


For the device file, the required dimensions are:
- a_wavelength: The absorption wavelength bins derived from the device file.
- c_wavelength: The attenuation wavelength bins derived from the device file.
- temperature_bin: The temperature bin associated with the pure water offset.

For the TS4.cor file, the required dimensions are:
- wavelength



<!--
## Installation

This package is available on [PyPI](https://pypi.org/project/acpype/) and can be installed using pip:
`pip install acspype`

-->

## Suggested Naming Conventions
acpype doesn't enforce a strict naming convention of variables created  using the software, but the examples provided use terminology commonly found in the ACS manual and in literature that uses ACS data.







#### SCRATCH



Serial Pipeline
The ACS communicates over RS232 and asynchronously sends data to produce a new binary packet at 4Hz. A serial pipeline should effectively seek out a full ACS binary packet within a serial port’s buffer, extract it, and then assign it a timestamp based on the host computer clock. Those looking to integrate the ACS into a data acquisition system should know that the instrument does not have a real-time clock and only reports the number of milliseconds that have passed since the system received power. After a full packet is identified, it can be parsed to obtain a set of ACS data with values in engineering units (e.g. counts). Using a corresponding device file (ACS-XXX.dev), the engineering units can then be converted to geophysical units which are indexed by wavelength for the absorption (a) and attenuation (c) channels. Real-time correction of an ACS spectrum in the same thread is limited to what can be done in between receiving packets. A full packet will be received roughly every 250 milliseconds. Anecdotally, computation up through measured absorption (am) and attenuation (cm) can safely be done in the same thread, but temperature-salinity correction may lead to serial buffer pile-up and inaccurate timestamps. Steps 0.0 to 1.5 in Table 1 describe what can be done in single thread data acquisition for the ACS.
acpype also separates itself from software like pyACS by providing functionality to parse an incoming ACS binary packet into engineering values without explicit knowledge of the information in the ACS’ accompanying device file. Conveniently, the first 31 and last 3 bytes of a packet from any ACS (structure version 3) maintain the same format description, and byte 31 provides the number of output wavelengths, which can then be used to calculate the remaining number of bytes within the packet and split it into raw engineering values. This makes it possible, as a recommended last resort, to acquire, identify, and log ACS data from any ACS on a serial port without knowledge of the information within the device file. It is extremely important to keep track of the separate device file that provides the factory (or user) created pure water offsets. A key challenge in processing and merging ACS data is that the number of wavelength outputs and wavelength bins are very rarely the same between instruments, as well as between the same instrument factory refurbishments. This leads to the possibility of incorrectly merging ACS data if not indexed by wavelength. 
After processing the measured value, it is up to the user to determine the storage format of data. acpype does not enforce a storage format or filetype for serial data acquisition. For applications that require concurrency, the use of PostgreSQL (PGDG, 2025) is suggested, as it supports the storage and retrieval of arrays. For file-based applications, using SQLite (SQLite Consortium, 2025) is suggested, but arrays will first need to be converted to character delimited strings. Logging directly to text or netCDF files is discouraged because memory use would increase as the dataset becomes larger, unless functionality is implemented to maintain a consistent file size.

Post-Processing Pipeline
It is presumed that the data processed with acpype will already have been timestamped and split into its engineering parts. acpype does require that the input data be placed in an Xarray.Dataset with the minimum coordinates of time, a_wavelength, and c_wavelength. Steps 0.1 to 5.1 describe the steps that should be taken as part of a standard post-processing pipeline for ACS data. Along this pipeline, several key pieces of information are needed 

A post-processing pipeline should first follow the conversions and corrections implemented in the suggestion serial pipeline. After calculation of the measured value, temperature-salinity correction should then be done for both attenuation and absorption. 




