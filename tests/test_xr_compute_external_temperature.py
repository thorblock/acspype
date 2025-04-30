from datetime import datetime
import xarray as xr
from acspype.processing import compute_external_temperature

TEST_DATA = xr.Dataset()
TEST_DATA = TEST_DATA.assign_coords({'time': [datetime(2025,4,28,0,0,0,100000),
                                              datetime(2025,4,28,0,0,0,350000),
                                              datetime(2025,4,28,0,0,0,600000),
                                              datetime(2025,4,28,0,0,0,850000),
                                              datetime(2025,4,28,0,0,1,100000)]})
TEST_DATA['raw_external_temperature'] = (['time'], [38000, 39000,40000,41000,42000])

def test():
    assert isinstance(TEST_DATA.raw_external_temperature, xr.DataArray)
    assert isinstance(compute_external_temperature(TEST_DATA.raw_external_temperature), xr.DataArray)