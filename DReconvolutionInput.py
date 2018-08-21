#*************************************************************************************************
#**
#** Copyright (c) 2017, 2018 Danny Petschke. All rights reserved.
#** 
#** Redistribution and use in source and binary forms, with or without modification, 
#** are permitted provided that the following conditions are met:
#**
#** 1. Redistributions of source code must retain the above copyright notice, 
#**    this list of conditions and the following disclaimer.
#**
#** 2. Redistributions in binary form must reproduce the above copyright notice, 
#**    this list of conditions and the following disclaimer in the documentation 
#**    and/or other materials provided with the distribution.
#**
#** 3. Neither the name of the copyright holder "Danny Petschke" nor the names of its  
#**    contributors may be used to endorse or promote products derived from this software  
#**    without specific prior written permission.
#**
#**
#** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS 
#** OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
#** MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
#** COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
#** EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
#** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
#** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
#** TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
#** EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#**
#** Contact: danny.petschke@uni-wuerzburg.de
#**
#*************************************************************************************************

from DReconvolutionModel import ReconvolutionModel as reconvModel

#save output as *.txt file after success?
__saveReconvolutionSpectrum             = False
__saveReconvolutionSpectrumPath         = 'testData/recovolutionSpectrumOutput.txt'
__saveReconvolutionSpectrumResidualPath = 'testData/recovolutionSpectrumResidualsOutput.txt'

#!note: IRF output is only saved if the model function is used, meaning--> (__bUsingModel = True)
__saveReconvolutionIRF                  = False
__saveReconvolutionIRFPath              = 'testData/recovolutionIRFOutput.txt'
__saveReconvolutionIRFResidualPath      = 'testData/recovolutionIRFResidualsOutput.txt'


#channel/bin resolution [ps]
__channelResolutionInPs = 5.0

#binning factor:
__binningFactor = 1;

#expected number of components (number of exponential decay functions - LIMITED to MAX: 4):
__numberOfExpDec = 2

#expected lifetimes (tau) -> start values in [ps] (required for the levenberg marquardt fit using lmfit library)
#note: only the first '__numberOfExpDec' related values are considered (e.g.: for __numberOfExpDec = 2 --> __expectedTau_1_in_ps AND __expectedTau_2_in_ps)
__expectedTau_1_in_ps = 240.0;
__expectedTau_2_in_ps = 1200.0;
__expectedTau_3_in_ps = 2800.0;
__expectedTau_4_in_ps = 160.0;

#background estimation (right side of spectrum data):
__bkgrd_startIndex = 8000;
__bkgrd_count = 999;

#fixed background? (value of estimated background is used)
__bkgrdFixed = False;


#NOTE: Spectrum and IRF data vectors require equal length!!!

#file path which contains the SPECTRUM data:
__filePathSpec = 'testData/spectrum_5ps.dat'
__specDataDelimiter = '\t'

#file path which contains the IRF data:
__filePathIRF = 'testData/irf_5ps.dat'
__irfDataDelimiter = '\t'

#define the number of rows which should be skipped during the import:
__skipRows = 0;

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
