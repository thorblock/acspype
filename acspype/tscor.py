import numpy as np
import os
import xarray as xr

from acspype.ts4cor import TS4COR


class ACSTSCor:
    """A wrapper class for the TS4COR dictionary that is included with the package."""

    def __init__(self):
        self.wavelengths = tuple(TS4COR.keys())
        self.psi_s_a = list(map(lambda x: x[1]['psi_s_a'], TS4COR.items()))
        self.psi_s_c = list(map(lambda x: x[1]['psi_s_c'], TS4COR.items()))
        self.psi_t = list(map(lambda x: x[1]['psi_t'], TS4COR.items()))
        self.method = 'Sullivan et al., 2006'

    def to_xarray(self):
        """
        Convert the TS4COR dictionary to an xarray dataset.

        :return: A dimensioned xarray.Dataset.
        """
        ds = xr.Dataset()
        ds = ds.assign_coords({'wavelength': np.array(self.wavelengths)})
        ds['psi_t'] = (['wavelength'], np.array(self.psi_t))
        ds['psi_s_c'] = (['wavelength'], np.array(self.psi_s_c))
        ds['psi_s_a'] = (['wavelength'], np.array(self.psi_s_a))

        ds.attrs['tscor_data'] = self.method
        return ds


class ACSTS4CorReader:
    """
    A class for parsing ACS TS4.cor files and putting them into a format that is easier to work with for larger or
    multiple file datasets.

    This class can be used if users want to import the TS4.cor file themselves.
    """

    def __init__(self, filepath: str) -> None:
        """
        Parse the .cor file and assign data as attributes.

        :param filepath: The filepath of the TS4.cor file.
        """

        self.filepath = os.path.normpath(filepath)
        self.__read_cor()
        self.__parse_lines()
        self.method = 'Sullivan et al., 2006'

    def __read_cor(self) -> None:
        """
        Read .cor file and store lines as a class attribute.

        :return: None
        """

        with open(self.filepath, 'r') as _file:
            self._lines = _file.readlines()

    def __parse_lines(self) -> None:
        """
        Parse the lines of the .cor file to get correction information.

        :return: None
        """

        wavelengths = []
        psi_t = []
        psi_s_c = []
        psi_s_a = []
        for line in self._lines:
            line_data = line.split('\t')
            line_data = [v.replace('\n', '') for v in line_data]
            line_data = [v.replace(' ', '') for v in line_data]
            if line_data == ['']:
                break
            line_data = [float(v) for v in line_data]
            wavelengths.append(line_data[0])
            psi_t.append(line_data[1])
            psi_s_c.append(line_data[2])
            psi_s_a.append(line_data[3])
        if len(wavelengths) != len(psi_t) != len(psi_s_c) != len(psi_s_a):
            raise ValueError('Mismatch in length of TS4cor file.')
        else:
            self.wavelengths = np.array(wavelengths)
            self.psi_t = np.array(psi_t)
            self.psi_s_c = np.array(psi_s_c)
            self.psi_s_a = np.array(psi_s_a)

    def to_xarray(self):
        """
        Convert the imported TS4 data to an xarray dataset.

        :return: A dimensioned xarray.Dataset.
        """
        ds = xr.Dataset()
        ds = ds.assign_coords({'wavelength': np.array(self.wavelengths)})
        ds['psi_t'] = (['wavelength'], np.array(self.psi_t))
        ds['psi_s_c'] = (['wavelength'], np.array(self.psi_s_c))
        ds['psi_s_a'] = (['wavelength'], np.array(self.psi_s_a))

        ds.attrs['tscor_data'] = self.method
        return ds
