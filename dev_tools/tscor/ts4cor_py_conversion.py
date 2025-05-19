import numpy as np
import os
from scipy.interpolate import make_interp_spline

from src.acspype import ACSTS4CorReader


def main():
    tscor = ACSTS4CorReader('TS4.cor')
    wavelengths = tscor.wavelengths
    psi_s_a = tscor.psi_s_a
    psi_s_c = tscor.psi_s_c
    psi_t = tscor.psi_t

    # Build Splines
    psi_s_a_spl = make_interp_spline(wavelengths, psi_s_a, k=1)
    psi_s_c_spl = make_interp_spline(wavelengths, psi_s_c, k=1)
    psi_t_spl = make_interp_spline(wavelengths, psi_t, k=1)

    # Extrapolate
    extrap_wvls = np.arange(395, 399.9, 0.1)
    extrap_psi_s_a = psi_s_a_spl(extrap_wvls)
    extrap_psi_s_c = psi_s_c_spl(extrap_wvls)
    extrap_psi_t = psi_t_spl(extrap_wvls)

    # Build new wavelengths and coefficients
    new_wavelengths = tuple(map(float, np.round(np.concatenate([extrap_wvls, wavelengths]), 6)))
    new_psi_s_a = tuple(map(float, np.round(np.concatenate([extrap_psi_s_a, psi_s_a]), 6)))
    new_psi_s_c = tuple(map(float, np.round(np.concatenate([extrap_psi_s_c, psi_s_c]), 6)))
    new_psi_t = tuple(map(float, np.round(np.concatenate([extrap_psi_t, psi_t]), 6)))

    tscor_dict = {}
    for wavelength in new_wavelengths:
        tscor_dict[wavelength] = {
            'psi_s_a': new_psi_s_a[new_wavelengths.index(wavelength)],
            'psi_s_c': new_psi_s_c[new_wavelengths.index(wavelength)],
            'psi_t': new_psi_t[new_wavelengths.index(wavelength)]
        }


    filename = 'ts4cor.py'
    filepath = f'../../acspype/{filename}'
    with open(filepath, 'w') as f:
        f.write('"""\nThis file was automatically generated from the TS4.cor file that accompanies all ACS sensors \nusing ts4cor_py_conversion.py. ')
        f.write('The TS4COR variable is a dictionary with the wavelength as the key and the \ncoefficients as the value in the form of a dictionary.\n')
        f.write('Values less than 400nm have been linearly extrapolated and rounded to 6 decimal places \nusing a SciPy BSpline and are not from SBS or Sullivan et al. 2006.\n')
        f.write('The TS4COR variable takes up approximately 150 kB of memory.\n"""')
        f.write('\n\n')
        f.write('TS4COR = {\n')
        for wavelength in new_wavelengths:
            f.write(f'\t{wavelength}: {tscor_dict[wavelength]},\n')
        f.write('}\n')

    if os.path.isfile(filepath):
        print(r"Created ts4cor.py in acspype package.")
    else:
        raise FileNotFoundError(f"File {filepath} not found.")

if __name__ == "__main__":
    main()