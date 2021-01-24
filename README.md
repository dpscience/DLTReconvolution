Support this project and keep always updated about recent software releases, bug fixes and major improvements by [following on github](https://github.com/dpscience?tab=followers).

![badge-followers](https://img.shields.io/github/followers/dpscience?style=social)
![badge-stars](https://img.shields.io/github/stars/dpscience/DLTReconvolution?style=social)
![badge-forks](https://img.shields.io/github/forks/dpscience/DLTReconvolution?style=social)

# DLTReconvolution

![badge-language](https://img.shields.io/badge/language-Python-blue)
![badge-license](https://img.shields.io/badge/license-BSD-green)

Copyright (c) 2017-2019 Danny Petschke (danny.petschke@uni-wuerzburg.de). All rights reserved.<br><br>
<b>DLTReconvolution</b> - A Python based software for the analysis of lifetime spectra using the iterative least-square reconvolution method.<br>

This program allows the analysis of lifetime spectra using ... 

* either the experimentally obtained instrumental response (IRF)  
* or a recorded mono-exponential decay spectrum such as from 207-Bi.

![DLTReconvolution output](/testData/demo.png)

# Quickstart Guide

`DLTReconvolution` consists of 3 files ...<br>

`DReconvolutionInput.py`<br>
`DReconvolutionModel.py`<br>
`DReconvolutionProc.py`<br>

* <b>edit</b> the input file `DReconvolutionInput.py` or just start the demonstration using the provided test data set:

```python
# note: spectrum and IRF (or mono-exponential decay spectrum) data vectors require equal length!

# file path (and name) to the SPECTRUM data:
__filePathSpec          = 'testData/spectrum_5ps.dat'
__specDataDelimiter     = '\t'

# file path (and name) to the IRF data:
__filePathIRF           = 'testData/irf_5ps.dat'
__irfDataDelimiter      = '\t'

# if TRUE, the fitted model function according to '__modelType' will be used as IRF data
__bUsingModel           = False

# when using model function ('__bUsingModel' = True) choose type of model (defined in DReconvolutionModel.py)
# ------------------
# Gaussian       = 1
# Lorentz_Cauchy = 2
# Pseudovoigt1   = 3
# Pearson7       = 4
# ------------------
__modelType             = reconvModel.Gaussian

# define the number of rows to be skipped during the import (e.g. for ignoring the header entries)
__skipRows              = 5;

# channel/bin resolution [ps]
__channelResolutionInPs = 5.0

# binning factor:
__binningFactor         = 1;

# expected number of components (number of exponential decay functions - LIMITED to 4)
__numberOfExpDec        = 3

# expected discrete characteristic lifetimes (decay = 1/tau) -> start values in units of picoseconds [ps]
# note: the values are considered in top-down order (e.g.: for '__numberOfExpDec' = 2 --> '__expectedTau_1_in_ps' AND '__expectedTau_2_in_ps' will be considered)
__expectedTau_1_in_ps   = 110.0;
__expectedTau_2_in_ps   = 375.0;
__expectedTau_3_in_ps   = 2200.0;
__expectedTau_4_in_ps   = 160.0;

# fit weighting: y variance? w = 1/sqrt(y) for poisson noise assumption otherwise the weighting is equally distributed w = 1.0
__bUsingYVarAsWeighting = True

# background estimation
__bkgrd_startIndex      = 8900;
__bkgrd_count           = 1000; # number of channels with respect to the 'startIndex'

# fixed background? if set True, the value of the estimated background based on the calculated mean ['__bkgrd_startIndex':'__bkgrd_startIndex' + '__bkgrd_count'] will be used
__bkgrdFixed            = False;

# set TRUE if the IRF should be retrieved from a mono-exponential decay spectrum such as well annealed metals (Al, Fe, ..) or the 207-Bi isotope using the 'graphical deconvolution' technique presented by Koechlin & Raviart (1964) (in this case the IRF data will be ignored)
__bUsingMonoDecaySpecForIRF               = False

# fixed mono-decay component in units of picoseconds [ps] (1/lambda = tau):
__tau_monoDecaySpec_in_ps                 = 182.0 #[ps]

__filePathMonoDecaySpec                   = 'C:/Users/.../207_Bi.dat'
__monoDecaySpecDataDelimiter              = '\t'

# data pre-processing for indirect IRF extraction from a mono-exponential decay spectrum using the 'graphical deconvolution' technique presented by Koechlin & Raviart (1964):
# 1. stage: re-binning >> 2. stage: smoothing
# Note: re-binning is only applied in case of '__bSmoothMonoDecaySpecForIRF' = True

# 1. stage: re-binning
__bReBinMonoDecaySpecForIRF               = False
__bReBinFacMonoDecaySpecForIRF            = 4

# 2. stage: smoothing by Savitzky-Golay filtering
__bSmoothMonoDecaySpecForIRF              = False
__SmoothingWindowDecaySpecForIRF          = 11
__SmoothingPolynomialOrderDecaySpecForIRF = 3

# set TRUE if the IRF data should be artificially broadened (~FWHM) applying an additional convolution using a Gaussian kernel (e.g. for compensation of energy differences)
__bUsingAdditionalGaussianKernel          = False

__gaussianKernelFWHM                      = 90.2  #[ps]
__bVaryGaussianKernelFWHM                 = False #if TRUE, this values will be used a an additional fitting parameter

# save output as *.txt file after success?
__saveReconvolutionSpectrum             = False
__saveReconvolutionSpectrumPath         = 'C:/Users/.../Bi_207_analytical_additionalConvKernel_fitdata.txt'
__saveReconvolutionSpectrumResidualPath = 'C:/Users/.../Bi_207_analytical_additionalConvKernel_residuals.txt'

__saveIRFSpectrum                       = False
__saveIRFSpectrumPath                   = 'C:/Users/.../Bi_207_analytical_additionalConvKernel_irfdata.txt'

__saveReconvolutionResults              = False
__saveReconvolutionResultsPath          = 'C:/Users/.../Bi_207_analytical_additionalConvKernel_results.txt'

# note: IRF output is only saved if the model function is used, meaning '__bUsingModel' = True
__saveReconvolutionIRF                  = False
__saveReconvolutionIRFPath              = 'output/...*txt'
__saveReconvolutionIRFResidualPath      = 'output/...*txt'
```

* <b>execute</b> `DReconvolutionProc.py`.

* <b>finished</b>. You should see the results as shown above in the figures.

# Related Publications using this Program

[![DOI](https://img.shields.io/badge/DOI-10.12693%2FAPhysPolA.137.171-yellowgreen)](http://doi.org/10.12693/APhysPolA.137.171)

[Limitations on the positron lifetime spectra decomposability applying the iterative least-square re-convolution method using the instrumental responses obtained from 207-Bi and 60-Co, Acta Physica Polonica A](http://doi.org/10.12693/APhysPolA.137.171).<br>

# How to cite this Software?

* <b>You should at least cite the applied version of this program in your study.</b><br>

You can cite all versions by using the <b>DOI 10.5281/zenodo.1255105</b>. This DOI represents all versions, and will always resolve to the latest one.<br>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1255105.svg)](https://doi.org/10.5281/zenodo.1255105)

## ``v1.x``
<b>DLTReconvolution v1.0</b><br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1255106.svg)](https://doi.org/10.5281/zenodo.1255106)<br>
<b>DLTReconvolution v1.1</b><br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1414107.svg)](https://doi.org/10.5281/zenodo.1414107)<br>
<b>DLTReconvolution v1.2</b><br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3464523.svg)](https://doi.org/10.5281/zenodo.3464523)<br>

# License of DLTReconvolution (BSD-3-Clause)

Copyright (c) 2017-2019 Danny Petschke (danny.petschke@uni-wuerzburg.de). All rights reserved.<br><br>

Redistribution and use in source and binary forms, with or without modification,<br> 
are permitted provided that the following conditions are met:<br><br>

 1. Redistributions of source code must retain the above copyright notice<br>
    this list of conditions and the following disclaimer.<br><br>

 2. Redistributions in binary form must reproduce the above copyright notice,<br> 
    this list of conditions and the following disclaimer in the documentation<br> 
    and/or other materials provided with the distribution.<br><br>

 3. Neither the name of the copyright holder "Danny Petschke" nor the names of<br> 
    its contributors may be used to endorse or promote products derived from <br>
    this software without specific prior written permission.<br><br>


 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS<br> 
 OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF<br> 
 MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE<br> 
 COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,<br> 
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF<br> 
 SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)<br> 
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR<br> 
 TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,<br> 
 EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.<br>
 
 see also [BSD-3-Clause License](https://opensource.org/licenses/BSD-3-Clause)
