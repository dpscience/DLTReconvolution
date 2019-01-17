# DLTReconvolution
Copyright (c) 2017-2019 Danny Petschke (danny.petschke@uni-wuerzburg.de). All rights reserved.<br><br>
DLTReconvolution - A Python based Software for the Analysis of Lifetime Spectra using the iterative least-square Reconvolution Method

# Short Introduction

![DLTReconvolution output](/testData/demo.png)

### DLTReconvolution consists of 3 python files:

- DReconvolutionInput.py
- DReconvolutionModel.py
- DReconvolutionProc.py

### How to start?

1. edit <b>DReconvolutionInput.py (or simply start with the provided test data...)</b>:

```python
#save output as *.txt file after success?
__saveReconvolutionSpectrum             = False
__saveReconvolutionSpectrumPath         = 'output/...*txt'
__saveReconvolutionSpectrumResidualPath = 'output/...*txt'

#!note: IRF output is only saved if the model function is used, meaning--> (__bUsingModel = True)
__saveReconvolutionIRF                  = False
__saveReconvolutionIRFPath              = 'output/...*txt'
__saveReconvolutionIRFResidualPath      = 'output/...*txt'

#channel/bin resolution [ps]
__channelResolutionInPs = 5.0

#binning factor:
__binningFactor = 1;

#expected number of components (number of exponential decay functions - LIMITED to MAX: 4):
__numberOfExpDec = 3

#expected lifetimes (tau) -> start values in [ps] (required for the levenberg marquardt fit using lmfit library)
#note: only the first '__numberOfExpDec' related values are considered (e.g.: for __numberOfExpDec = 2 --> __expectedTau_1_in_ps AND __expectedTau_2_in_ps)
__expectedTau_1_in_ps = 190.0;
__expectedTau_2_in_ps = 385.0;
__expectedTau_3_in_ps = 1200.0;
__expectedTau_4_in_ps = 160.0;

#background estimation (right side of spectrum data):
__bkgrd_startIndex = 8000;
__bkgrd_count = 999;

#fixed background? (value of estimated background is used)
__bkgrdFixed = False;

#NOTE: Spectrum and IRF (or mono-exponential decay spectrum) data vectors require equal length!!!

#define the number of rows which should be skipped during the import:
__skipRows = 5;

#file path which contains the SPECTRUM data:
__filePathSpec = 'testData/bi207_ref/sn_5ps.dat'
__specDataDelimiter = '\t'

#file path which contains the IRF data:
__filePathIRF = 'path to irf data'
__irfDataDelimiter = '\t'


#set TRUE if the irf is retrieved from a mono-exponential decay spectrum such as well annealed metals or 207-Bi (in this case, the IRF data will be ignored):
__bUsingMonoDecaySpecForIRF = True

__tau_monoDecaySpec_in_ps = 182.0 #[ps] 207-Bi

__filePathMonoDecaySpec = 'testData/bi207_ref/bi207_5ps.dat'
__monoDecaySpecDataDelimiter = '\t'

#using model function for IRF?
__bUsingModel = False

#fit weighting: y variance? w = 1/sqrt(y) <--- <poisson noise> otherwise the weighting is equally distributed: w = 1
__bUsingYVarAsWeighting = True

#if using model function? choose type of model (defined in DReconvolutionModel.py):
#------------------
#Gaussian       = 1
#Lorentz_Cauchy = 2
#Pseudovoigt1   = 3
#Pearson7       = 4
#------------------
__modelType = reconvModel.Pearson7
```
2. run <b>DReconvolutionProc.py</b>
3. <b>finished!</b>

Note: all parameter estimations required for the best model fit using the levenberg marquardt algorithm are done automatically. 

# How to cite this Software?

You can cite all versions by using the <b>DOI 10.5281/zenodo.1255105</b>. This DOI represents all versions, and will always resolve to the latest one.<br>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1255105.svg)](https://doi.org/10.5281/zenodo.1255105)

## v1.x
DLTReconvolution v1.0:<br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1255106.svg)](https://doi.org/10.5281/zenodo.1255106)<br>
DLTReconvolution v1.1:<br>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1414107.svg)](https://doi.org/10.5281/zenodo.1414107)<br>
DLTReconvolution v1.2:<br><b>not officialy released yet...<br>

# Requirements
- [LMFIT](https://lmfit.github.io/lmfit-py/model.html) 
- [NumPy](http://www.numpy.org/) 
- [matplotlib](https://matplotlib.org/)<br>

#### [WinPython](https://sourceforge.net/projects/winpython/) meets all requirements. 

# License (BSD-3-Clause)

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