import xarray as xr


def estimate_chl_alh(a_p: xr.DataArray, alh_star: float = 0.0104) -> xr.DataArray:
    """
    Compute chlorophyll-a from absorption line height via Roesler and Barnard 2013.
    https://www.sciencedirect.com/science/article/pii/S2211122013000509

    :param a_p: Particulate absorption data with wavelength as a coordinate.
    :param alh_star: Absorption line height coefficient value, default is 0.0104,
        which is the average from Table 1 in Roesler and Barnard 2013.
    :return: Chlorophyll-a concentration in mg/m^3.
    """

    a650 = a_p.sel(wavelength=650, method='nearest')
    a676 = a_p.sel(wavelength=676, method='nearest')
    a715 = a_p.sel(wavelength=715, method='nearest')

    abl = ((a715 - a650) / (715 - 650)) * (676 - 650) + a650  #EQ 1 in Roesler and Barnard 2013
    alh = a676 - abl  #EQ 2 in Roesler and Barnard 2013
    chl_alh = alh / alh_star #EQ 3 in Roesler and Barnard 2013

    chl_alh.attrs['alh_star'] = alh_star
    chl_alh.attrs['method'] = 'Roesler and Barnard, 2013'
    chl_alh.attrs['ancillary_variables'] = str(a_p.name)
    return chl_alh


def estimate_poc(c_p: xr.DataArray, method = 'gardner_etal_2006') -> xr.DataArray:
    """
    Estimate particulate organic carbon (POC) from attenuation based on linear models defined by the issued method.

    Possible methods include:
    gardner_etal_2006 -> From Gardner et al. 2006
    cetenic_etal_2012 -> From Cetinic et al. 2012
    stramski_etal_2008 -> From Stramski et al. 2008
    behrenfeld_and_boss_2006 -> From Behrenfeld and Boss 2006

    :param c_p: Particulate beam attenuation coefficient data with wavelength as a coordinate.
    :param method:
    :return:
    """
    method = method.lower()
    if method == 'gardner_etal_2006':
        wvl = 660
        _c_p = c_p.sel(wavelength=wvl, method='nearest')
        m = 381
        b = 9.4
        poc = m * _c_p + b
        return poc
    elif method == 'cetenic_etal_2012':
        wvl = 660
        _c_p = c_p.sel(wavelength=wvl, method='nearest')
        m = 391
        b = -5.8
        poc = m * _c_p + b
        return poc
    elif method == 'stramski_etal_2008':
        wvl = 660
        _c_p = c_p.sel(wavelength=wvl, method='nearest')
        m = 458
        b = 0
        poc = m * _c_p + b
        return poc
    elif method == 'behrenfeld_and_boss_2006':
        wvl = 660
        _c_p = c_p.sel(wavelength=wvl, method='nearest')
        m = 585
        b = 7.6
        poc = m * _c_p + b
        return poc
    