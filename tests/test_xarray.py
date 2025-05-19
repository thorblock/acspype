import numpy as np
import xarray as xr

from src.acspype import ACSDev, ACSTSCor
import src.acspype.processing as acsproc

TEST_FILEPATH = '../dev_tools/test_files/TEST_SHIMADA_202405.nc'
TEST_DEV = '../dev_tools/test_files/ACS-00412_2023-05-10.dev'


def test():
    dev = ACSDev(TEST_DEV)
    acs = xr.open_dataset(TEST_FILEPATH)
    tscor = ACSTSCor().to_xarray()
    assert np.all(dev.a_wavelength == acs.a_wavelength)

    # Compute Temperatures from ACS Thermistors
    acs['internal_temperature'] = acsproc.compute_internal_temperature(acs.raw_internal_temperature)
    acs['external_temperature'] = acsproc.compute_external_temperature(acs.raw_external_temperature)

    # Compute Uncorrected Absorption and Attenuation
    acs['a_uncorrected'] = acsproc.compute_uncorrected(acs.a_signal, acs.a_reference, dev.path_length)
    acs['c_uncorrected'] = acsproc.compute_uncorrected(acs.c_signal, acs.c_reference, dev.path_length)

    # Compute Measured Absorption and Attenuation
    acs['a_m'] = acsproc.compute_measured(acs.a_uncorrected, acs.internal_temperature, dev.a_offset, dev.func_a_delta_t)
    acs['c_m'] = acsproc.compute_measured(acs.c_uncorrected, acs.external_temperature, dev.c_offset, dev.func_c_delta_t)

    # Discontinuity Correction
    discontinuity_index = acsproc.find_discontinuity_index(acs.a_wavelength, acs.c_wavelength)
    acs['a_m'], acs['a_discontinuity_offset'] = acsproc.discontinuity_correction(acs.a_m,
                                                                                 discontinuity_index,
                                                                                 'a_wavelength')
    acs['c_m'], acs['c_discontinuity_offset'] = acsproc.discontinuity_correction(acs.c_m,
                                                                                 discontinuity_index,
                                                                                 'c_wavelength')

    # Temperature and Salinity Correction
    psi_s_a = tscor.psi_s_a.sel(wavelength=np.array(acs.a_wavelength), method='nearest')
    psi_t = tscor.psi_t.sel(wavelength=np.array(acs.a_wavelength), method='nearest')
    acs['a_mts'] = acsproc.ts_correction(acs.a_m, acs.sea_water_temperature, acs.sea_water_practical_salinity, psi_t,
                                         psi_s_a, dev.tcal)

    psi_s_c = tscor.psi_s_c.sel(wavelength=np.array(acs.c_wavelength), method='nearest')
    psi_t = tscor.psi_t.sel(wavelength=np.array(acs.c_wavelength), method='nearest')
    acs['c_mts'] = acsproc.ts_correction(acs.c_m, acs.sea_water_temperature, acs.sea_water_practical_salinity, psi_t,
                                         psi_s_c, dev.tcal)

    # Zero Shift Correction
    acs['a_mts'] = acsproc.zero_shift_correction(acs.a_mts)
    acs['c_mts'] = acsproc.zero_shift_correction(acs.c_mts)

    acs = acsproc.interpolate_common_wavelengths(acs, 'a_wavelength', 'c_wavelength', 'wavelength')

    ref_a = acs.a_mts.sel(wavelength=715)
    acs['a_mts_baseline'] = acsproc.baseline_scattering_correction(acs.a_mts, ref_a)

    acs['a_mts_fixed'] = acsproc.fixed_scattering_correction(acs.a_mts, acs.c_mts)

    ref_a = acs.a_mts.sel(wavelength=715)
    ref_c = acs.c_mts.sel(wavelength=715)
    acs['a_mts_proportional'] = acsproc.proportional_scattering_correction(acs.a_mts, acs.c_mts, ref_a, ref_c)

    gelbstoff = acs.where(acs.seawater_state == 1, drop=True)

    assert gelbstoff.seawater_state.max() == 1
