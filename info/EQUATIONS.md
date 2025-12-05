
$a+b=c$

## Absorption
$a_t = a_w + \Sigma a_p + \Sigma a_g$

$
a(\lambda)_{\mathrm{uncorr}} = \frac{1}{x} \cdot \ln\!\left( 
\frac{a(\lambda)_{\mathrm{signal}}}{a(\lambda)_{\mathrm{reference}}}
\right)
$

$a(\lambda)_{m\_discontinuity} =  a(\lambda)_{offset}-a(\lambda)_{uncorr}-\Delta T_a$

## Attenuation
$c_t = c_w + \Sigma c_p + \Sigma c_g$

$
a(\lambda)_{\mathrm{uncorr}} = \frac{1}{x} \cdot \ln\!\left( 
\frac{a(\lambda)_{\mathrm{signal}}}{a(\lambda)_{\mathrm{reference}}}
\right)
$

## Spectrum Discontinuity Correction
$
a(\lambda)_{\text{m}} = a(\lambda)_{\text{m\_discontinuity}} + \mathrm{offset}_d
\left\{
\begin{array}{lll}
\lambda_0 \rightarrow \lambda_{\text{discontinuity\_idx}} &,& 0 \\
\lambda_{\text{discontinuity\_idx}+1} \rightarrow \lambda_n &,& \text{offset}_d
\end{array}
\right.
$

## Temperature Salinity Dependence
$
a(\lambda)_{mts} = a(\lambda)_m - \bigl[\psi_t \, \cdot (T - T_{\mathrm{ref}}) + \psi_{sa}\,\cdot S\bigr]
$

## Absorption Scattering Correction
$
a(\lambda)_{\text{mts\_baseline}} = a(\lambda)_{mts} - a(\mathrm{ref})_{mts}
$

$
a(\lambda)_{\text{mts\_fixed}} = a(\lambda)_{mts} - \varepsilon \cdot \bigl(c(\lambda)_{mts} - a(\lambda)_{mts}\bigr)
$

$
a(\lambda)_{\text{mts\_proportional}} = a(\lambda)_{mts} - \left(\left(
\frac{a(\mathrm{ref})_{mts}}
     {c(\mathrm{ref})_{mts} - a(\mathrm{ref})_{mts}}
\right) \cdot \bigl( c(\lambda)_{mts} - a(\lambda)_{mts} \bigr)\right)
$