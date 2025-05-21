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

    alh = (a_p_676 - (39 / 65 * a_p_650) + (26 / 65 * a_p_715))

    if isinstance(alh, xr.DataArray):
        alh.attrs['ancillary_variables'] = [a_p_650.name, a_p_676.name, a_p_715.name]
        alh.attrs['long_name'] = 'Absorption Line Height at 676 nm'
        alh.attrs['units'] = 'm^-1'
        alh.attrs['reference'] = 'Boss et al. (2007)'
        alh.attrs['reference_doi'] = 'https://doi.org/10.1007/978-1-4020-5824-0_9'
    return alh


def estimate_chl(a_p_650: xr.DataArray,
                 a_p_676: xr.DataArray,
                 a_p_715: xr.DataArray,
                 alh_coeff: float = 0.0104) -> xr.DataArray:
    """
    Estimate chlorophyll-a from absorption line height via Roesler and Barnard 2013.
    https://doi.org/10.1016/j.mio.2013.12.003

    This function has only been tested on xarray.Datasets.

    :param a_p_650: a_p data at 650 nm.
    :param a_p_676: a_p data at 676 nm.
    :param a_p_715: a_p data at 715 nm.
    :param alh_coeff: (a_p line height coefficient value, default is 0.0104,
        which is the average from Table 1 in Roesler and Barnard 2013.
    :return: Chlorophyll-a concentration in mg/m^3.
    """

    abl = ((a_p_715 - a_p_650) / (715 - 650)) * (676 - 650) + a_p_650  # EQ 1 in Roesler and Barnard 2013
    alh = a_p_676 - abl  # EQ 2 in Roesler and Barnard 2013
    chl_alh = alh / alh_coeff  # EQ 3 in Roesler and Barnard 2013

    if isinstance(chl_alh, xr.DataArray):
        chl_alh.attrs['ancillary_variables'] = a_p_650.name
        chl_alh.attrs['alh_coeff'] = alh_coeff
        chl_alh.attrs['long_name'] = 'Chlorophyll-a Concentration from Absorption Line Height'
        chl_alh.attrs['units'] = 'mg/m^3'
        chl_alh.attrs['estimation_method'] = 'Roesler and Barnard (2013)'
        chl_alh.attrs['method_data_source'] = 'Cultures and Gulf of Maine'
        chl_alh.attrs['reference'] = 'Roesler and Barnard (2013)'
        chl_alh.attrs['reference_doi'] = 'https://doi.org/10.1016/j.mio.2013.12.003'
    return chl_alh


def estimate_poc(c_p_660: xr.DataArray, slope_offset: str | tuple | list) -> xr.DataArray:
    """
    Compute particulate organic carbon from particle beam attenuation at 660nm using a provided slope and offset.
    If the slope_offset value is a tuple or list, the first value is the slope and the second value is the offset.

    If the slope_offset value is a string, a predefined slope and offset will be used. Strings represent the first
    author and the year of the publication.

    Note: If supplying your own slope and offset, please ensure that the units are correct.

    :param c_p_660: particulate attenuation (c_p) at 660 nm.
    :param slope_offset: A string indicating the literary source for the slope and offset or
        a tuple or list containing the slope and offset values.
    :return: Estimated particulate organic carbon in mg/m3.
    """

    if isinstance(slope_offset, str):
        if slope_offset == 'gardner2006':  # Gardner et al. 2006, Global, All Seasons
            slope = 381
            offset = 9.4
            attrs = {'estimation_method': 'Gardner et al. (2006)',
                     'method_data_source': 'Global, All Seasons',
                     'reference': 'Gardner et al. (2006)',
                     'reference_doi': 'https://doi.org/10.1016/j.dsr2.2006.01.029'}
        elif slope_offset == 'stramski2008':  # Stramski et al. 2008, Fall, Pacific, Atlantic
            slope = 458
            offset = 10.7
            attrs = {'estimation_method': 'Stramski et al. (2008)',
                     'method_data_source': 'Fall, Pacific, Atlantic',
                     'reference': 'Stramski et al. (2008)',
                     'reference_doi': 'https://doi.org/10.5194/bg-5-171-2008'}
        elif slope_offset == 'behrenfeld2006':  # Behrenfeld and Boss 2006, Fall, Equatorial Pacific
            slope = 585
            offset = 7.6
            attrs = {'estimation_method': 'Behrenfeld and Boss (2006)',
                     'method_data_source': 'Fall, Equatorial Pacific',
                     'reference': 'Behrenfeld and Boss (2006)',
                     'reference_doi': 'https://doi.org/10.1357/002224006778189563'}
        elif slope_offset == 'cetenic2008':  # Cetinić et al. 2008, Spring, NE Atlantic
            slope = 391
            offset = -5.8
            attrs = {'estimation_method': 'Cetinić et al. (2008)',
                     'method_data_source': 'Spring, NE Atlantic',
                     'reference': 'Cetinić et al. (2008)',
                     'reference_doi': ' https://doi.org/10.1029/2011JC007771 '}
    elif isinstance(slope_offset, tuple | list):
        slope, offset = slope_offset

    poc = slope * c_p_660 + offset  # y = m * x + b

    if isinstance(poc, xr.DataArray):
        for k, v in attrs.items():
            poc.attrs[k] = v
        poc.attrs['ancillary_variables'] = c_p_660.name
        poc.attrs['slope'] = slope
        poc.attrs['offset'] = offset
        poc.attrs['long_name'] = 'Particulate Organic Carbon from Particle Beam Attenuation (c_p) at 660nm'
        poc.attrs['units'] = 'mg/m^3'
    return poc
