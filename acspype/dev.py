from datetime import datetime
import numpy as np
import re
from scipy import interpolate
import xarray as xr

from acspype.core import NUM_PAT


class ACSDev:
    """
    A class for parsing ACS .dev files and putting them into a format that is easier to work with for larger or
    multiple file datasets.

    Generally, users will not call individual functions, but rather use the class to obtain attributes, which are
    created at class instantiation or convert the data to an xr.Dataset using the to_xarray function.
    """

    def __init__(self, filepath: str) -> None:
        """
        Run the following functions at instantiation to parse the .dev file and store the info as class attributes.

        :param filepath: The filepath to the .dev file.
        :return: None
        """

        self._filepath = filepath
        self.__read_dev()
        self.__parse_metadata()
        self.__parse_tbins()
        self.__parse_offsets()
        self.__build_interp_funcs()
        self.__check_parse()


    def __read_dev(self) -> None:
        """
        Import the .dev file as a text file.
        The file contents are stored as a list of strings in the class attribute self._content.

        :return: None
        """

        with open(self._filepath, 'r') as _file:
            self._content = _file.readlines()


    def __parse_metadata(self) -> None:
        """
        Parse the .dev file for individual sensor metadata.
        Sensor specific metadata are stored as class attributes.

        :return: None
        """

        metadata_lines = [line for line in self._content if 'C and A offset' not in line]
        for line in metadata_lines:
            if 'ACS Meter' in line:
                self.sensor_type = re.findall('(.*?)\n', line)[0]
            elif 'Serial' in line:
                self.sn_hexdec = re.findall('(.*?)\t', line)[0]
                self.sn = 'ACS-' + str(int(self.sn_hexdec[-6:], 16)).zfill(5)  # Convert to sn shown on product sticker.
            elif 'structure version' in line:
                self.structure_version = int(re.findall(f'({NUM_PAT})\t', line)[0])
            elif 'tcal' in line or 'Tcal' in line:
                self.tcal, self.ical = [float(v) for v in re.findall(f': ({NUM_PAT}) C', line)]
                cal_date_str = re.findall('file on (.*?)[.]', line)[0].replace(' ', '')
                try:  # Sometimes the file date is entered as yyyy or yy. This should handle both cases.
                    self.cal_date = datetime.strptime(cal_date_str, '%m/%d/%Y').strftime('%Y-%m-%d')
                except:
                    self.cal_date = datetime.strptime(cal_date_str, '%m/%d/%y').strftime('%Y-%m-%d')
            elif 'Depth calibration' in line:
                (self.depth_cal_1,
                 self.depth_cal_2) = [float(v) for v in re.findall(f'({NUM_PAT})', line)]
            elif 'Baud' in line:
                self.baudrate = int(re.findall(f'({NUM_PAT})\t', line)[0])
            elif 'Path' in line:
                self.path_length = float(re.findall(f'({NUM_PAT})\t', line)[0])
            elif 'wavelengths' in line:
                self.num_wavelength = int(re.findall(f'({NUM_PAT})\t', line)[0])
            elif 'number of temperature bins' in line:
                self.num_tbin = int(re.findall(f'({NUM_PAT})\t', line)[0])
            elif 'maxANoise' in line:
                (self.max_a_noise, self.max_c_noise, self.max_a_nonconform, self.max_c_nonconform,
                 self.max_a_difference, self.max_c_difference, self.min_a_counts,
                 self.min_c_counts, self.min_r_counts, self.max_temp_sd,
                 self.max_depth_sd) = [float(v) for v in re.findall(f'({NUM_PAT})\t', line)]


    def __parse_tbins(self) -> None:
        """
        Parse the .dev file for temperature bin information.

        :return: None
        """
        tbin_line = [line for line in self._content if '; temperature bins' in line][0]
        tbins = tbin_line.split('\t')
        tbins = [v for v in tbins if v]  # Toss empty strings.
        tbins = [v for v in tbins if v != '\n']  # Toss newline characters.
        self.tbin = np.array([float(v) for v in tbins if 'temperature bins' not in v])  # Convert to float and toss comment.


    def __parse_offsets(self) -> None:
        """
        Parse the .dev file for a and c offsets. Data are saved as class attributes for access at a later time.

        :return: None
        """

        offset_lines = [line for line in self._content if 'C and A offset' in line]

        # Create holder arrays to loop over and append data to.
        c_wvls = []
        a_wvls = []
        c_offs = []
        a_offs = []
        c_deltas = []
        a_deltas = []
        wavelength_color_schemes = []

        for line in offset_lines:
            offsets, c_delta, a_delta = line.split('\t\t')[:-1]
            c_wvl, a_wvl, wvl_color, c_off, a_off = offsets.split('\t')

            # Convert strings to proper pythonic datatypes.
            c_wvl = float(c_wvl.replace('C', ''))
            a_wvl = float(a_wvl.replace('A', ''))
            c_off = float(c_off)
            a_off = float(a_off)
            c_delta = [float(v) for v in c_delta.split('\t')]
            a_delta = [float(v) for v in a_delta.split('\t')]

            # Append files to holder arrays.
            c_wvls.append(c_wvl)
            a_wvls.append(a_wvl)
            c_offs.append(c_off)
            a_offs.append(a_off)
            c_deltas.append(c_delta)
            a_deltas.append(a_delta)
            wavelength_color_schemes.append(wvl_color)

        # Convert holder arrays to numpy arrays.
        self.c_wavelength = np.array(c_wvls)
        self.a_wavelength = np.array(a_wvls)
        self.c_offset = np.array(c_offs)
        self.a_offset = np.array(a_offs)
        self.c_delta_t = np.array(c_deltas)
        self.a_delta_t = np.array(a_deltas)
        self.wavelength_color_schemes = wavelength_color_schemes


    def __build_interp_funcs(self) -> None:
        """
        Build interpolation functions for the a and c delta_t values and store as class attributes.

        :return: None
        """
        self.func_a_delta_t = interpolate.interp1d(self.tbin, self.a_delta_t, axis=1)
        self.func_c_delta_t = interpolate.interp1d(self.tbin, self.c_delta_t, axis=1)
        self.delta_t_interp_method = 'scipy.interpolate.interp1d'


    def __check_parse(self) -> None:
        """
        Verify that the shape of the data is as expected.

        :return: None
        """

        if len(self.a_wavelength) != self.num_wavelength:
            raise ValueError('Mismatch between number of wavelengths extracted for A and expected from file.'
                             'Please verify the .dev file integrity.')
        if len(self.c_wavelength) != self.num_wavelength:
            raise ValueError('Mismatch between number of wavelengths extracted for C and expected from file.'
                             'Please verify the .dev file integrity.')
        if len(self.c_wavelength) != len(self.a_wavelength):
            raise ValueError('Mismatch between number of wavelengths extracted for A and C.'
                             'Please verify the .dev file integrity.')
        if np.array(self.a_delta_t).shape != (len(self.a_wavelength), self.num_tbin):
            raise ValueError('Mismatch between length of A wavelengths and number of temperature bins.'
                             'Please verify the .dev file integrity.')
        if np.array(self.c_delta_t).shape != (len(self.a_wavelength), self.num_tbin):
            raise ValueError('Mismatch between length of C wavelengths and number of temperature bins.'
                             'Please verify the .dev file integrity.')


    def to_xarray(self) -> xr.Dataset:
        """
        Convert the parsed .dev file files to an xarray dataset

        Returns: An appropriately dimensioned xarray dataset containing device file files.
        """
        ds = xr.Dataset()
        ds = ds.assign_coords({'a_wavelength': self.a_wavelength,
                               'c_wavelength': self.c_wavelength,
                               'temperature_bin': self.tbin})

        ds['a_offset'] = (['a_wavelength'], self.a_offset)
        ds['a_delta_t'] = (['a_wavelength', 'temperature_bin'], self.a_delta_t)

        ds['c_offset'] = (['c_wavelength'], np.array(self.c_offset))
        ds['c_delta_t'] = (['c_wavelength', 'temperature_bin'], self.c_delta_t)

        ds.attrs['device_filepath'] = self._filepath
        ds.attrs['sensor_type'] = self.sensor_type
        ds.attrs['serial_number_hexdec'] = self.sn_hexdec
        ds.attrs['serial_number'] = self.sn
        ds.attrs['device_file_structure_version'] = self.structure_version
        ds.attrs['tcal'] = self.tcal
        ds.attrs['ical'] = self.ical
        ds.attrs['calibration_date'] = self.cal_date
        ds.attrs['depth_cal_1'] = self.depth_cal_1
        ds.attrs['depth_cal_2'] = self.depth_cal_2
        ds.attrs['baudrate'] = self.baudrate
        ds.attrs['path_length'] = self.path_length
        ds.attrs['number_of_wavelength_bins'] = self.num_wavelength
        ds.attrs['number_of_temperature_bins'] = self.num_tbin
        ds.attrs['max_a_noise'] = self.max_a_noise
        ds.attrs['max_c_noise'] = self.max_c_noise
        ds.attrs['max_a_nonconform'] = self.max_a_nonconform
        ds.attrs['max_c_nonconform'] = self.max_c_nonconform
        ds.attrs['max_a_difference'] = self.max_a_difference
        ds.attrs['max_c_difference'] = self.max_c_difference
        ds.attrs['min_a_counts'] = self.min_a_counts
        ds.attrs['min_c_counts'] = self.min_c_counts
        ds.attrs['min_r_counts'] = self.min_r_counts
        ds.attrs['max_temp_sd'] = self.max_temp_sd
        ds.attrs['max_depth_sd'] = self.max_depth_sd
        return ds
