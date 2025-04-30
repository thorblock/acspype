from datetime import datetime
import numpy as np
import xarray as xr

from acspype.processing import compute_uncorrected

TEST_DATA = xr.Dataset()
TEST_DATA = TEST_DATA.assign_coords({'time': [datetime(2025,4,28,0,0,0,100000),
                                              datetime(2025,4,28,0,0,0,350000),
                                              datetime(2025,4,28,0,0,0,600000),
                                              datetime(2025,4,28,0,0,0,850000),
                                              datetime(2025,4,28,0,0,1,100000)],
                                     'a_wavelength': [400,404.3, 408.6, 412.9, 417.2]})
TEST_DATA['a_signal'] = (['time','a_wavelength'],
                         (np.random.rand(len(TEST_DATA.time.values),
                                         len(TEST_DATA.a_wavelength.values))*10000).astype(int))
TEST_DATA['a_reference'] = (0.5 * TEST_DATA.a_signal).astype(int)


def test():
    assert isinstance(TEST_DATA['a_signal'], xr.DataArray)
    assert isinstance(compute_uncorrected(TEST_DATA.a_signal, TEST_DATA.a_reference),xr.DataArray)


