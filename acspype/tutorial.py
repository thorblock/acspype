from importlib.resources import files
import xarray as xr

import acspype
from acspype import ACSDev

def load_shimada_tutorial_data():
    source = files(acspype).joinpath("tutorial_data/EXAMPLE_SHIMADA_202405.nc")
    return xr.open_dataset(source)

def load_shimada_tutorial_device_file():
    source = files(acspype).joinpath("tutorial_data/ACS-00412_2023-05-10.dev")
    return ACSDev(source)