import numpy as np
import os
from scipy.interpolate import interp1d

def main():
    new_wavelengths = np.arange(395, 755.1, 0.1).tolist()

    # Values from Sullivan et al. 2006.
    # Sorry to whomever looks at this code. It was late and I couldn't think of a better way to write out the values from Table 1.
    sigma_psi_t = {0.0002: [400, 700, 702],
                   0.0001: list(np.arange(402, 434, 0.1)) +
                           list(np.arange(444, 448, 0.1)) +
                           list(np.arange(506, 536)) +
                           list(np.arange(536, 564, 1)) + [544] + list(
                       np.arange(570, 698, 0.1)),
                   0.0000: list(np.arange(436, 442, 0.1)) +
                           list(np.arange(450, 504, 0.1)) +
                           list(np.arange(538, 542, 0.1)) + list(
                       np.arange(546, 568, 0.1)),
                   0.0003: [704, 706, 736, 738, ],
                   0.0004: [708, 710, 712, 714, 734, 740],
                   0.0005: [716, 718, 732, 742],
                   0.0006: [720, 722, 730, 744],
                   0.0007: [724, 726, 728, 746],
                   0.0008: [748],
                   0.0009: [750]}

    sigma_psi_s_c = {0.0004: list(np.arange(400, 412, 0.1)),
                     0.0003: list(np.arange(414, 462, 0.1)) + list(
                         np.arange(746, 750, 0.1)),
                     0.0002: list(np.arange(464, 548, 0.1)) + list(
                         np.arange(740, 744, 0.1)),
                     0.0001: list(np.arange(550, 738, 0.1)), }

    sigma_psi_s_a = {0.00003: list(np.arange(400, 438, 0.1)) + [750],
                     0.00002: list(np.arange(440, 506, 0.1)) + list(
                         np.arange(744, 748, 0.1)),
                     0.00001: list(np.arange(508, 742, 0.1))}


    interp_sigma_psi_t = interp_sigma(sigma_psi_t, new_wavelengths)
    interp_sigma_psi_s_c = interp_sigma(sigma_psi_s_c, new_wavelengths)
    interp_sigma_psi_s_a = interp_sigma(sigma_psi_s_a, new_wavelengths)

    sigma_dict = {}
    for wavelength in new_wavelengths:
        idx = new_wavelengths.index(wavelength)
        sigma_dict[round(wavelength,1)] = {'sigma_psi_s_a': float(interp_sigma_psi_s_a[idx]),
                                  'sigma_psi_s_c': float(interp_sigma_psi_s_c[idx]),
                                  'sigma_psi_t': float(interp_sigma_psi_t[idx])}


    filename = 'ts4cor_sigma.py'
    filepath = f'../../acspype/{filename}'
    with open(filepath, 'w') as f:
        f.write('"""\nThis file was automatically generated with the standard deviations found in Table 1 of\n'
                'Sullivan et al. 2006. (10.1364/AO.45.005294)\n')
        f.write('The TS4COR_SIGMA variable is a dictionary with the wavelength as the key and the '
                '\nstandard deviations as the '
                'value in the form of a dictionary.\n')
        f.write('Sullivan et al. 2006 provide these values in 2nm bins between 400nm and 750nm,\n'
                'so values have been nearest neighbor interpolated.')
        f.write('Values less than 400nm have been nearest\nneighbor interpolated using a '
                'SciPy interp1d and are not from\nSBS or Sullivan et al. 2006.\n"""')
        f.write('\n\n')
        f.write('TS4COR_SIGMA = {\n')
        for wavelength in new_wavelengths:
            f.write(f'\t{round(wavelength,1)}: {sigma_dict[round(wavelength,1)]},\n')
        f.write('}\n')

    if os.path.isfile(filepath):
        print(r"Created ts4cor_sigma.py in acspype package.")
    else:
        raise FileNotFoundError(f"File {filepath} not found.")



def interp_sigma(vals: dict, new_wavelengths: list) -> list:
    """
    Interpolate sigma values using the nearest value.

    :param vals: A mish-mashed dictionary compiled haphazardly and manually from Sullivan et al 2006.
    :param new_wavelengths: The wavelengths to do nearest interpolation to.
    :return: A list of sigma values.
    """
    data = {}
    for k, v in vals.items():
        for wvl in v:
            data[float(round(wvl, 1))] = k
    data = dict(sorted(data.items()))
    x, y = zip(*data.items())
    f = interp1d(x, y, kind='nearest', fill_value='extrapolate')
    interp_sigma = f(new_wavelengths)
    return interp_sigma


if __name__ == '__main__':
    main()