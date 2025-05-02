import xarray as xr


def compute_alh_676(a_p_650: xr.DataArray, a_p_676: xr.DataArray, a_p_715: xr.DataArray) -> xr.DataArray:
    """
    Compute absorption line height at 676 nm via Boss et al, 2007.
    https://link.springer.com/chapter/10.1007/978-1-4020-5824-0_9

    :param a_p_650: a_p data at 650 nm.
    :param a_p_676: a_p data at 676 nm.
    :param a_p_715: a_p data at 715 nm.
    :return: Absorption line height at 676 nm in m^-1.
    """

    alh = (a_p_676 - (39/65 * a_p_650) + (26/65 * a_p_715))
    return alh


def compute_chl_alh(a_p_650: xr.DataArray,
                    a_p_676: xr.DataArray,
                    a_p_715: xr.DataArray,
                    alh_star: float = 0.0104) -> xr.DataArray:
    """
    Compute chlorophyll-a from absorption line height via Roesler and Barnard 2013.
    https://www.sciencedirect.com/science/article/pii/S2211122013000509

    :param a_p_650: a_p data at 650 nm.
    :param a_p_676: a_p data at 676 nm.
    :param a_p_715: a_p data at 715 nm.
    :param alh_star: a_p line height coefficient value, default is 0.0104,
        which is the average from Table 1 in Roesler and Barnard 2013.
    :return: Chlorophyll-a concentration in mg/m^3.
    """

    abl = ((a_p_715 - a_p_650) / (715 - 650)) * (676 - 650) + a_p_650  #EQ 1 in Roesler and Barnard 2013
    alh = a_p_676 - abl  #EQ 2 in Roesler and Barnard 2013
    chl_alh = alh / alh_star #EQ 3 in Roesler and Barnard 2013
    return chl_alh
