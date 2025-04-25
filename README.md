# acspype
acspype provides functions for reading [Sea-Bird Scientific ACS](https://www.seabird.com/ac-s-spectral-absorption-and-attenuation-sensor/product?id=60762467715) data over serial and for performing processing on ACS data.
Some of the target audiences for this package include technicians looking to integrate an ACS into an existing data acquisition system, researchers looking to process prior collected ACS data, and data scientists looking to redistribute ACS data in HDF formats.

ACS data are inherently complex and difficult to work with, particularly for new users without strong optics backgrounds that desire to use the data for empirically derived data products like chlorophyll-a and particulate organic carbon.
This package attempts to simplify the process of ingesting and processing ACS data so that users can more quickly get to the science application in their research.

## CAUTION
### Equipment
If you are using this package to acquire data from an ACS over serial, it is your responsibility to ensure that you are 
using the appropriate equipment (power supply, cable, etc.). Please consult Sea-Bird Scientific for more information about the relevant equipment and
power requirements needed to acquire data from the ACS.

If looking for a USB-to-RS232 adapter, the [FTDI US232R-100-BULK](https://ftdichip.com/products/us232r-100-bulk/) is recommended.

### Serial Processing Limitations
Not all processing steps offered in acspype can be done in the same thread as the serial data acquisition. 
Some steps require linear and cubic interpolation, which may take too much time to perform in between packets. 
Processing up through the measured absorption (a_m) and attenuation (c_m) can be done in the same thread in real-time, but temperature-salinity correction may lead to serial buffer pile-up and inaccurate timestamps.
Multi-threading may be a solution, but we will leave that up to you, the user, to implement.

## Installation
This package is available on the [Python Package Index (PyPI)](https://pypi.org/project/acspype/) and can be installed via pip:

`pip install acspype`

## Issues and Discussions
If you experience an issue related to acspype, please use the [GitHub Issues](https://github.com/IanTBlack/acspype/issues) page to report it.

If you would like to propose functionality or to discuss best practices, please start a discussion on the [GitHub Discussions](https://github.com/IanTBlack/acspype/discussions) page.

## Contributions
Contributions to acspype are encouraged via GitHub [pull request](https://github.com/IanTBlack/acspype/pulls). 
acspype will be updated on the PyPI as patches and major changes are implemented.


# Overview
## What does acspype need from you (and Sea-Bird Scientific)?
In addition to the physical sensor and a test cable, one who deploys an ACS will also receive a device file (ACS-XXX.dev) and a temperature-salinity correction file (TS4.cor) from Sea-Bird Scientific.

The device file is a text file that contains pure water offsets, delta T values, and the wavelength bins for your specific ACS.
No two device files are alike, making it critical to keep track of this information in order to compute advanced products.

The TS4.cor file is a text file that contains empirically derived temperature and salinity correction coefficients needed to correct absorption and attenuation data for changes in temperature and salinity.
These coefficients have been derived and interpolated to 0.1 nm bins from the original data presented in [Sullivan et al. 2006](https://doi.org/10.1364/AO.45.005294).
The TS4.cor file received from Sea-Bird Scientific is not unique. This file has been converted and stored within acspype for convenience.

## What does this package **not** do?
This package does not provide any functionality for logging data from the ACS. The file or database type is entirely left up to the user.
For file-based logging, we recommend using SQLite, which can be used to store data within multiple tables in a single file. 
For application where concurrency is required, PostgreSQL is a good option. It provides functionality for handling arrays, which is how ACS data are best represented.
The ACS produces a significant amount of data, so logging to text (.txt, .csv) or netCDF (.nc) files may best be done as hourly or daily files to prevent excessive memory usage.

## Naming and Processing Conventions
Although it was designed to not enforce strict naming conventions, there are some functions that require the specific coordinates of 'a_wavelength', 'c_wavelength', or 'wavelength'.
It is strongly encouraged that users follow the naming conventions used in the examples and within the metadata file. acspype contains pre-defined and dynamic metadata that clearly identifies products and the ancillary data used to create them.




### Serial Pipeline Example
The ACS communicates over RS232 and asynchronously sends data to produce a new binary packet at 4Hz. A serial pipeline should effectively seek out a full ACS binary packet within a serial port’s buffer, extract it, and then assign it a timestamp based on the host computer clock. Those looking to integrate the ACS into a data acquisition system should know that the instrument does not have a real-time clock and only reports the number of milliseconds that have passed since the system received power. After a full packet is identified, it can be parsed to obtain a set of ACS data with values in engineering units (e.g. counts). Using a corresponding device file (ACS-XXX.dev), the engineering units can then be converted to geophysical units which are indexed by wavelength for the absorption (a) and attenuation (c) channels. Real-time correction of an ACS spectrum in the same thread is limited to what can be done in between receiving packets. A full packet will be received roughly every 250 milliseconds. Anecdotally, computation up through measured absorption (am) and attenuation (cm) can safely be done in the same thread, but temperature-salinity correction may lead to serial buffer pile-up and inaccurate timestamps. Steps 0.0 to 1.5 in Table 1 describe what can be done in single thread data acquisition for the ACS.
acspype also separates itself from software like pyACS by providing functionality to parse an incoming ACS binary packet into engineering values without explicit knowledge of the information in the ACS’ accompanying device file. Conveniently, the first 31 and last 3 bytes of a packet from any ACS (structure version 3) maintain the same format description, and byte 31 provides the number of output wavelengths, which can then be used to calculate the remaining number of bytes within the packet and split it into raw engineering values. This makes it possible, as a recommended last resort, to acquire, identify, and log ACS data from any ACS on a serial port without knowledge of the information within the device file. It is extremely important to keep track of the separate device file that provides the factory (or user) created pure water offsets. A key challenge in processing and merging ACS data is that the number of wavelength outputs and wavelength bins are very rarely the same between instruments, as well as between the same instrument factory refurbishments. This leads to the possibility of incorrectly merging ACS data if not indexed by wavelength. 
After processing the measured value, it is up to the user to determine the storage format of data. acspype does not enforce a storage format or filetype for serial data acquisition. For applications that require concurrency, the use of PostgreSQL (PGDG, 2025) is suggested, as it supports the storage and retrieval of arrays. For file-based applications, using SQLite (SQLite Consortium, 2025) is suggested, but arrays will first need to be converted to character delimited strings. Logging directly to text or netCDF files is discouraged because memory use would increase as the dataset becomes larger, unless functionality is implemented to maintain a consistent file size.


### Post-Processsing Pipeline Example

A post-processing pipeline should first follow the conversions and corrections implemented in the suggestion serial pipeline. After calculation of the measured value, temperature-salinity correction should then be done for both attenuation and absorption. 
It is presumed that the data processed with acspype will already have been timestamped and split into its engineering parts. acspype does require that the input data be placed in an Xarray.Dataset with the minimum coordinates of time, a_wavelength, and c_wavelength. Steps 0.1 to 5.1 describe the steps that should be taken as part of a standard post-processing pipeline for ACS data. Along this pipeline, several key pieces of information are needed 






