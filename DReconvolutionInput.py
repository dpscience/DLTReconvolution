#*************************************************************************************************
#**
#** DLTReconvolution v1.3 (06.12.2019)
#**
#**
#** Copyright (c) 2017 - 2019 Danny Petschke. All rights reserved.
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

#NOTE: spectrum and IRF (or mono-exponential decay spectrum) data vectors require equal length!

#file path (and name) to the SPECTRUM data:
__filePathSpec          = 'testData/spectrum_5ps.dat'
__specDataDelimiter     = '\t'

#file path (and name) to the IRF data:
__filePathIRF           = 'testData/irf_5ps.dat'
__irfDataDelimiter      = '\t'

#if TRUE, the fitted model function according to '__modelType' will be used as IRF data
__bUsingModel           = False

#if using model function? choose type of model (defined in DReconvolutionModel.py):
#------------------
#Gaussian       = 1
#Lorentz_Cauchy = 2
#Pseudovoigt1   = 3
#Pearson7       = 4
#------------------
__modelType             = reconvModel.Gaussian

#define the number of rows, which should be skipped during the import (e.g. for ignoring the header entries):
__skipRows              = 5;

#channel/bin resolution [ps]
__channelResolutionInPs = 5.0

#binning factor:
__binningFactor         = 1;

#expected number of components (number of exponential decay functions - LIMITED to MAX: 4):
__numberOfExpDec        = 3

#expected discrete characteristic lifetimes (tau) -> start values in units of picoseconds [ps]
#note: the values are considered in top-down order (e.g.: for __numberOfExpDec = 2 --> __expectedTau_1_in_ps AND __expectedTau_2_in_ps are considered)
__expectedTau_1_in_ps   = 108.0;
__expectedTau_2_in_ps   = 385.0;
__expectedTau_3_in_ps   = 2200.0;
__expectedTau_4_in_ps   = 160.0;

#fit weighting: y variance? w = 1/sqrt(y) <--- <assumption: poisson noise> otherwise the weighting is equally distributed: w = 1.0
__bUsingYVarAsWeighting = True

#background estimation:
__bkgrd_startIndex      = 9000;
__bkgrd_count           = 1000; # number of channels with respect to the 'startIndex'

#fixed background? >> if True, the value of the estimated background based on the calculated mean [__bkgrd_startIndex:__bkgrd_startIndex + __bkgrd_count] will be used
__bkgrdFixed            = False;


#set TRUE if the irf should be retrieved from a mono-exponential decay spectrum such as well annealed metals (Al, Fe, ..) or the 207-Bi isotope using the 'graphical deconvolution' technique presented by Koechlin & Raviart (1964) (in this case, the IRF data will be ignored):
__bUsingMonoDecaySpecForIRF               = False

#fixed mono-decay component in units of picoseconds [ps] (1/lambda = tau):
__tau_monoDecaySpec_in_ps                 = 182.0 #[ps]

__filePathMonoDecaySpec                   = 'C:/Users/.../207_Bi.dat'
__monoDecaySpecDataDelimiter              = '\t'

#data pre-processing for indirect IRF extraction from a mono-exponential decay spectrum using the 'graphical deconvolution' technique presented by Koechlin & Raviart (1964):
#1. stage: re-binning >> 2. stage: smoothing
#Note: re-binning is only applied in case of '__bSmoothMonoDecaySpecForIRF = True'

#1. stage: re-binning:
__bReBinMonoDecaySpecForIRF               = False
__bReBinFacMonoDecaySpecForIRF            = 4

#2. stage: smoothing by Savitzky-Golay filtering:
__bSmoothMonoDecaySpecForIRF              = False
__SmoothingWindowDecaySpecForIRF          = 11
__SmoothingPolynomialOrderDecaySpecForIRF = 3


#set TRUE if the irf data should be artificially broadened (~FWHM) applying an additional convolution using a Gaussian kernel (e.g. for compensation of energy differences)
__bUsingAdditionalGaussianKernel          = False

__gaussianKernelFWHM                      = 90.2  #[ps]
__bVaryGaussianKernelFWHM                 = False #if TRUE, this values will be used a an additional fitting parameter

#set TRUE if synthetically generated data should be used:
__bUsingSimSpectra                        = True

__gaussianIRFFWHM                         = 235.0 #[ps]

__numberOfChannels                        = 10000
__integralCountsSpectrum                  = 5000000
__integralCountsGaussianIRF               = 5000000
__constBackground                         = 10

__simI_1                                  = 0.85
__simI_2                                  = 0.147
__simI_3                                  = 0.003
__simI_4                                  = 0.0
__simI_5                                  = 0.0


#save output as *.txt file after success?
__saveReconvolutionSpectrum             = False
__saveReconvolutionSpectrumPath         = 'C:/Users/.../Bi_207_analytical_additionalConvKernel_fitdata.txt'
__saveReconvolutionSpectrumResidualPath = 'C:/Users/.../Bi_207_analytical_additionalConvKernel_residuals.txt'

__saveIRFSpectrum                       = False
__saveIRFSpectrumPath                   = 'C:/Users/.../Bi_207_analytical_additionalConvKernel_irfdata.txt'

__saveReconvolutionResults              = False
__saveReconvolutionResultsPath          = 'C:/Users/.../Bi_207_analytical_additionalConvKernel_results.txt'

#Note: IRF output is only saved if the model function is used, meaning--> (__bUsingModel = True)
__saveReconvolutionIRF                  = False
__saveReconvolutionIRFPath              = 'output/...*txt'
__saveReconvolutionIRFResidualPath      = 'output/...*txt'
