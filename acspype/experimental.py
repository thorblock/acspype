import xarray as xr


def compute_chl_alh(absorption: xr.DataArray, alh_star: float = 0.0104) -> xr.DataArray:
    """
    Compute chlorophyll-a from absorption line height via Roesler and Barnard 2013.
    https://www.sciencedirect.com/science/article/pii/S2211122013000509

    :param absorption: Absorption data with wavelength as a coordinate.
    :param alh_star: Absorption line height coefficient value, default is 0.0104,
        which is the average from Table 1 in Roesler and Barnard 2013.
    :return: Chlorophyll-a concentration in mg/m^3.
    """

    a650 = absorption.sel(wavelength=650, method='nearest')
    a676 = absorption.sel(wavelength=676, method='nearest')
    a715 = absorption.sel(wavelength=715, method='nearest')

    abl = ((a715 - a650) / (715 - 650)) * (676 - 650) + a650  #EQ 1 in Roesler and Barnard 2013
    alh = a676 - abl  #EQ 2 in Roesler and Barnard 2013
    chl_alh = alh / alh_star #EQ 3 in Roesler and Barnard 2013

    chl_alh.attrs['alh_star'] = alh_star
    chl_alh.attrs['method'] = 'Roesler and Barnard, 2013'
    chl_alh.attrs['ancillary_variables'] = str(absorption.name)
    return chl_alh

#
# def compute_poc(attenuation: xr.DataArray, method = 'gardner_et_al_2006'):
#     """
#     Compute particulate organic carbon (POC) from attenuation based on linear models defined by the issued method.
#
#     :param attenuation: The particulate attenuation coefficient (c_p).
#     :param method: The method to use for POC calculation. Options are:
#         gardner_et_al_2006
#         behrenfeld_and_boss_2006
#         stramski_et_al_2008
#         cetenic_et_al_2012
#         goni_et_al_2021 -> Oregon Coast, August, 2011, Sigma 22-23
#     :return: POC in mg/m3
#     """
#
#     if method == 'gardner_et_al_2006':
#         m = 381  #POC to c_p slope in units of mgC/m2
#         c_p = attenuation.sel(wavelength=660, method='nearest')
#         b = 9.4
#     elif method == 'behrenfeld_and_boss_2006':
#         m = 585
#         c_p = attenuation.sel(wavelength=660, method='nearest')
#         b = 7.6
#     elif method == 'stramski_et_al_2008':
#         m = 458
#         c_p = attenuation.sel(wavelength=660, method='nearest')
#         b = 0
#     elif method == 'cetenic_et_al_2012':
#         m = 391
#         c_p = attenuation.sel(wavelength=660, method='nearest')
#         b = -5.8
#     elif method == 'goni_et_al_2021':
#         m = 38.9* 12.01
#         c_p = attenuation.sel(wavelength=650, method='nearest')
#         b = 0
#     poc = m * c_p + b
#     return poc