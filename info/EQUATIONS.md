
## Absorption, Scattering, and Attenuation
#### Equation 1
$$a+b=c$$
where $a$ is the absorption coefficient, $b$ is the scattering coefficient, and $c$ is the attenuation coefficient.


## Absorption Processing
#### Equation 2a - Total Absorption Constituents
$$a_t = a_w + \Sigma a_p + \Sigma a_g$$
where $a_t$ is the total absorption of seawater and its constituents, $a_w$ is the absorption by water, $\Sigma a_p$ is the absorption by all particulates, and $\Sigma a_g$ is the absorption by all gelbstoff (dissolved) constituents.

#### Equation 3a - Uncorrected Absorption from Raw Counts
$$
a(\lambda)_{\mathrm{uncorr}} = \frac{1}{x} \cdot \ln\!\left( 
\frac{a(\lambda)_{\mathrm{signal}}}{a(\lambda)_{\mathrm{reference}}}
\right)
$$
where $a(\lambda)_{\mathrm{uncorr}}$ is the uncorrected absorption, x is the path length, $a(\lambda)_{\mathrm{signal}}$ is raw signal counts for absorption, and $a(\lambda)_{\mathrm{reference}}$ is raw reference counts for absorption.

#### Equation 4a - Measured Absorption, Corrected for Internal Temperature, Uncorrected for Spectra Discontinuity
$$a(\lambda)_{m\_discontinuity} =  a(\lambda)_{offset}-a(\lambda)_{uncorr}-\Delta T_a$$
where $a(\lambda)_{m\_discontinuity}$ is the measured absorption not corrected for the inherent discontinuity,  $a(\lambda)_{offset}$ is the absorption offset from the device file, $a(\lambda)_{uncorr}$ is the uncorrected absorption from Equation 3a, and $\Delta T_a$ is the linearly interpolated internal temperature correction. 


#### Equation 5a - Measured Absorption, Corrected for Spectra Discontinuity
$$
a(\lambda)_{\text{m}} = a(\lambda)_{\text{m_discontinuity}} + \mathrm{offset}_d
\left\{
\begin{array}{lll}
\lambda_0 \rightarrow \lambda_{\text{discontinuity_idx}} &,& 0 \\
\lambda_{\text{discontinuity_idx}+1} \rightarrow \lambda_n &,& \text{offset}_d
\end{array}
\right.
$$
where $a(\lambda)_{m\_discontinuity}$ is from Equation 4a,  and $offset_d$ is a scalar offset applied to the latter half of a spectra, derived using a cubic spline. 


#### Equation 6a - Temperature Salinity Corrected Absorption
$$
a(\lambda)_{mts} = a(\lambda)_m - \bigl[\psi_t \, \cdot (T - T_{\mathrm{ref}}) + \psi_{sa}\,\cdot S\bigr]
$$

where $a(\lambda)_{mts}$ is the measured absorption, corrected for temperature and salinity dependence, $a(\lambda)_m$ is derived from Equation 5a, $\psi_t$ and $\psi_{sa}$ are derived from the TS4.cor file, and $T$ and $S$ are derived from an ancillary instrument, $T_{ref}$ is from the device file.

### Absorption Scattering Corrections
#### Equation 7.1a - Baseline Scattering Correction
$$
a(\lambda)_{mts\_{baseline}} = a(\lambda)_{mts} - a(\lambda_{ref})_{mts}
$$
where $a(\lambda)_{mts\_{baseline}}$ is the absorption corrected for scattering using the baseline method, $a(\lambda)_{mts}$ is from Equation 6a and $\lambda_{ref}$ is the absorption at the reference wavelength, typically 715nm.

#### Equation 7.2a - Fixed Scattering Correction
$$
a(\lambda)_{mts\_{fixed}} = a(\lambda)_{mts} - \varepsilon \cdot (c(\lambda)_{mts} - a(\lambda)_{mts})
$$
where $a(\lambda)_{mts\_{fixed}}$ is the absorption corrected for scattering using the fixed method, $a(\lambda)_{mts}$ and $c(\lambda)_{mts}$ are from Equation 6a/6c and $\varepsilon$ is an empirically derived coefficient.

#### Equation 7.3a - Proportional Scattering Correction
$$
a(\lambda)_{\text{mts\_proportional}} = a(\lambda)_{mts} - \left(\left(
\frac{a(\mathrm{\lambda_ref})_{mts}}
     {c(\mathrm{\lambda_ref})_{mts} - a(\mathrm{\lambda_ref})_{mts}}
\right) \cdot \bigl( c(\lambda)_{mts} - a(\lambda)_{mts} \bigr)\right)
$$
where $a(\lambda)_{\text{mts\_proportional}}$ is the absorption corrected for scattering using the proportional method, $\lambda_{ref}$ is the absorption/attenuation at the reference wavelength, typically 715nm.



## Attenuation Processing
#### Equation 2c - Total Attenuation Constituents
$$c_t = c_w + \Sigma c_p + \Sigma c_g$$
where $c_t$ is the total attenuation of seawater and its constituents, $c_w$ is the attenuation by water, $\Sigma c_p$ is the attenuation by all particulates, and $\Sigma c_g$ is the attenuation by all gelbstoff (dissolved) constituents.

#### Equation 3c - Uncorrected Attenuation From Raw Counts
$$
c(\lambda)_{\mathrm{uncorr}} = \frac{1}{x} \cdot \ln\!\left( 
\frac{c(\lambda)_{\mathrm{signal}}}{c(\lambda)_{\mathrm{reference}}}
\right)
$$
where $c(\lambda)_{\mathrm{uncorr}}$ is the uncorrected attenuation, x is the path length, $c(\lambda)_{\mathrm{signal}}$ is raw signal counts for attenuation, and $c(\lambda)_{\mathrm{reference}}$ is raw reference counts for attenuation.


#### Equation 4c - Measured Attenuation, Corrected for Internal Temperature, Uncorrected for Spectra Discontinuity
$$c(\lambda)_{m\_discontinuity} =  c(\lambda)_{offset}-c(\lambda)_{uncorr}-\Delta T_a$$
where $c(\lambda)_{m\_discontinuity}$ is the measured attenuation not corrected for the inherent discontinuity,  $c(\lambda)_{offset}$ is the attenuation offset from the device file, $c(\lambda)_{uncorr}$ is the uncorrected attenuation from Equation 3c, and $\Delta T_c$ is the linearly interpolated internal temperature correction. 


#### Equation 5c - Measured Attenuation, Corrected for Spectra Discontinuity
$$
c(\lambda)_{\text{m}} = c(\lambda)_{\text{m_discontinuity}} + \mathrm{offset}_d
\left\{
\begin{array}{lll}
\lambda_0 \rightarrow \lambda_{\text{discontinuity_idx}} &,& 0 \\
\lambda_{\text{discontinuity_idx}+1} \rightarrow \lambda_n &,& \text{offset}_d
\end{array}
\right.
$$
where $c(\lambda)_{m\_discontinuity}$ is from Equation 4c,  and $offset_d$ is a scalar offset applied to the latter half of a spectra, derived using a cubic spline. 

#### Equation 6c - Temperature Salinity Corrected Attenuation
$$
c(\lambda)_{mts} = c(\lambda)_m - \bigl[\psi_t \, \cdot (T - T_{\mathrm{ref}}) + \psi_{sc}\,\cdot S\bigr]
$$
where $c(\lambda)_{mts}$ is the measured attenuation, corrected for temperature and salinity dependence, $c(\lambda)_m$ is derived from Equation 5c, $\psi_t$ and $\psi_{sc}$ are derived from the TS4.cor file, and $T$ and $S$ are derived from an ancillary instrument, $T_{ref}$ is from the device file.
