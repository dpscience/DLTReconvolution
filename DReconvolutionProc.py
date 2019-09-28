#*************************************************************************************************
#**
#** DLTReconvolution v1.2 (21.09.2019)
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

import DReconvolutionModel as functionModelList
import DReconvolutionInput as userInput

from lmfit import Model
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.signal import savgol_filter
from scipy.stats import norm
import sys,os

def wiener_deconvolution(signal, kernel, SNR):
    kernel = np.hstack((kernel, np.zeros(len(signal) - len(kernel)))) # zero pad the kernel to same length
    H = np.fft.fft(kernel)
    deconvolved = np.real(np.fft.ifft(np.fft.fft(signal)*np.conj(H)/(H*np.conj(H) + SNR**2)))
    return deconvolved

# only applicable on noiseless data
def deconvolveData(ab, b):
    AB = np.fft.fft(ab);
    B = np.fft.fft(b);
    deconvA = np.real(np.fft.ifft(AB/B));
    return deconvA;

print("importing spectrum & irf data...");

xSpec,ySpec = np.loadtxt(userInput.__filePathSpec, delimiter=userInput.__specDataDelimiter, skiprows=userInput.__skipRows, unpack=True, dtype='float');

if not userInput.__bUsingMonoDecaySpecForIRF:
    xIRF,yIRF   = np.loadtxt(userInput.__filePathIRF, delimiter=userInput.__irfDataDelimiter, skiprows=userInput.__skipRows, unpack=True, dtype='float');
else:
    xIRF,yMonoDecaySpec = np.loadtxt(userInput.__filePathMonoDecaySpec, delimiter=userInput.__monoDecaySpecDataDelimiter, skiprows=userInput.__skipRows, unpack=True, dtype='float');

    #plt.semilogy(ySpec, 'ro', yMonoDecaySpec, 'bo');
    #plt.show()

    if userInput.__tau_monoDecaySpec_in_ps <= 0:
        print("Input error: negative lifetime for mono-decay spectrum detected.");
        quit(); #kill the process on error.

    yIRF = np.zeros(xSpec.size);
    yIRF += 1 #prevent zero values

    if userInput.__bSmoothMonoDecaySpecForIRF:
        if userInput.__bReBinMonoDecaySpecForIRF:
            newLength = 1+((int)(len(xIRF)/userInput.__bReBinFacMonoDecaySpecForIRF))
            
            rebin = 0;
            yMonoSpecIRFRebinned = np.zeros(newLength)
            xValRebinned         = np.zeros(newLength)

            for i in range(0, len(xIRF)):
                xValRebinned[rebin]          = rebin*userInput.__bReBinFacMonoDecaySpecForIRF;
                yMonoSpecIRFRebinned[rebin] += yMonoDecaySpec[i]

                if (((i+1)%userInput.__bReBinFacMonoDecaySpecForIRF) == 1):
                    rebin += 1

            yMonoSpecIRFRebinned[0] = yMonoSpecIRFRebinned[1]
            
            sv = savgol_filter(yMonoSpecIRFRebinned, userInput.__SmoothingWindowDecaySpecForIRF, userInput.__SmoothingPolynomialOrderDecaySpecForIRF) 
            cs = CubicSpline(xValRebinned, sv)

            #plt.semilogy(yMonoSpecIRFRebinned, 'r');
            #plt.semilogy(cs(xValRebinned), 'bo');
            #plt.show()

            for i in range(0, len(yMonoDecaySpec)):
                if cs(xIRF[i]) < 0.0:
                    yMonoDecaySpec[i] = 0.0
                else:
                    yMonoDecaySpec[i] = cs(xIRF[i])

            # formula according to Koechlin & Raviart (1964)
            for i in range(1, len(xIRF)-1):
                yIRF[i] = yMonoDecaySpec[i] + (userInput.__tau_monoDecaySpec_in_ps/(2*userInput.__channelResolutionInPs))*(yMonoDecaySpec[i+1] - yMonoDecaySpec[i-1])
                
                if yIRF[i] < 0:
                    yIRF[i] = 0

            yIRF[0] = yIRF[1]
            yIRF[len(yIRF)-1] = yIRF[len(yIRF)-2]
            
            #plt.semilogy(yIRF, 'r');
            #plt.semilogy(yMonoDecaySpec, 'bo');
            #plt.show()
        else:
            sv = savgol_filter(yMonoDecaySpec, userInput.__SmoothingWindowDecaySpecForIRF, userInput.__SmoothingPolynomialOrderDecaySpecForIRF) 
            cs = CubicSpline(xIRF, sv)
            
            # formula according to Koechlin & Raviart (1964)
            for i in range(1, len(xIRF)-1):
                yIRF[i] = cs(xIRF[i]) + (userInput.__tau_monoDecaySpec_in_ps/(2*userInput.__channelResolutionInPs))*(cs(xIRF[i+1]) - cs(xIRF[i-1]))

                if yIRF[i] < 0:
                    yIRF[i] = 0

            yIRF[0] = yIRF[1]
            yIRF[len(yIRF)-1] = yIRF[len(yIRF)-2]
            
            #plt.semilogy(xIRF, cs(xIRF), 'r');
            #plt.semilogy(yIRF, 'ro', yMonoDecaySpec, 'bo');
            #plt.show()
    else:
        # formula according to Koechlin & Raviart (1964)
        for i in range(1, len(xIRF)-1):
            yIRF[i] = yMonoDecaySpec[i] + (userInput.__tau_monoDecaySpec_in_ps/(2*userInput.__channelResolutionInPs))*(yMonoDecaySpec[i+1] - yMonoDecaySpec[i-1])

            if yIRF[i] < 0:
                    yIRF[i] = 0

        yIRF[0] = yIRF[1]
        yIRF[len(yIRF)-1] = yIRF[len(yIRF)-2]
            
        #plt.semilogy(yIRF, 'bo', yMonoDecaySpec, 'r');
        #plt.show()


yIRFOrigin = yIRF;

#shift data in the way that: vec[0] = 0:
xSpec -= xSpec[0];
xIRF  -= xIRF[0];

if (len(xSpec) != len(xIRF) or len(ySpec) != len(yIRF)):
    print("Error: same length is required for the data vectors (x,y) of IRF and SPECTRUM.");
    quit(); #kill the process on error.

#merging into one vector (xVal) for x values (Fourier transforms are 1D):
xVal = np.zeros(xSpec.size);
xVal = xSpec;

#re-bin data if required: 
if userInput.__binningFactor > 1:
    print("data re-binning by factor: {0} ...".format(userInput.__binningFactor))
    
    newLength = 1+((int)(len(xVal)/userInput.__binningFactor));
    
    xValRebinned = np.zeros(newLength);

    rebin = 0;
    ySpecRebinned = np.zeros(newLength);
    yIRFRebinned = np.zeros(newLength);

    for i in range(0, len(xVal)):
        xValRebinned[rebin] = rebin;
        ySpecRebinned[rebin] += ySpec[i]
        yIRFRebinned[rebin] += yIRF[i]

        if (((i+1)%userInput.__binningFactor) == 1):
            rebin += 1;

    xVal = np.zeros(newLength);
    yIRF = np.zeros(newLength);
    ySpec = np.zeros(newLength);
    
    xVal = xValRebinned;
    yIRF = yIRFRebinned;
    ySpec = ySpecRebinned;

    userInput.__channelResolutionInPs *= userInput.__binningFactor;

    userInput.__bkgrd_startIndex /= userInput.__binningFactor;
    userInput.__bkgrd_count /= userInput.__binningFactor;

#calculate fit weightings:
fitWeightingIRF  = 1.0;
fitWeightingSpec = 1.0;

print("calculating individual fit weights...")

if userInput.__bUsingYVarAsWeighting:
    fitWeightingIRF  = 1.0/np.sqrt(yIRF +1); #prevent zero devision
    fitWeightingSpec = 1.0/np.sqrt(ySpec+1); #prevent zero devision

#estimate IRF start values (amplitude, xOffset/mean, stddev) <<---- based on Gaussian and should work for the most symmetric distribution functions:
yMaxIRF        = np.amax(yIRF)
yMaxSpec       = np.amax(ySpec)
 
xWhereYMaxIRF  = np.argmax(yMaxIRF)  
xWhereYMaxSpec = np.argmax(yMaxSpec) 


stddevIRF          = 1.0
invStddevToFWHMFac = 1.0/(2*np.sqrt(2*np.log(2)));

for i in range(0, len(xVal)):
    if yIRF[i] > yMaxIRF*0.5:
        stddevIRF = np.abs((xWhereYMaxIRF-xIRF[i]))*invStddevToFWHMFac;
        break;

estimatedBkgrd    = 0.0;
estimatedBkgrdIRF = 0.0;
    
for i in range((int)(userInput.__bkgrd_startIndex), (int)(userInput.__bkgrd_startIndex + userInput.__bkgrd_count)):
    estimatedBkgrd    += ySpec[i];
    estimatedBkgrdIRF += yIRF[i];
    
estimatedBkgrd    /= userInput.__bkgrd_count;
estimatedBkgrdIRF /= userInput.__bkgrd_count;

#substract background from SPECTRUM and IRF - not used yet -:
ySpecMinusBkgrd = ((ySpec-estimatedBkgrd)+abs(ySpec-estimatedBkgrd))/2;
yIRFMinusBkgrd  = ((yIRF-estimatedBkgrdIRF)+abs(yIRF-estimatedBkgrdIRF))/2;

def convolveData(a, b):
    A = np.fft.fft(a);
    B = np.fft.fft(b);
    convAB = np.real(np.fft.ifft(A*B));
    return convAB;

# gaussian distribution function: G(mu = t_zero_in_ps, fwhm)
def generateGaussian(binWidth_in_ps=5.0, 
                     numberOfIntegralCounts=5000000, 
                     numberOfBins=10000,
                     tZero_in_ps = 0.0,
                     fwhm_in_ps = 230.0):
    
    # may change later for multiple gaussian functions
    numberOfComponents = 1
    intensitiesOfGaussian = [1.0]
    
    timeBin_in_ps = np.zeros(numberOfBins)
    counts_y      = np.zeros(numberOfBins)
    
    countsInitial = np.zeros(numberOfComponents)
    areaInitial   = np.zeros(numberOfComponents) 
    
    sumOfCounts   = 0
    
    sigma         = fwhm_in_ps/(2*np.sqrt(2*np.log(2)))

    for i in range(0, numberOfComponents):
        countsInitial[i] = numberOfIntegralCounts*intensitiesOfGaussian[i]
    
    for bin in range(0, numberOfBins - 1):
        timeBin_in_ps[bin] = (2*bin + 1)*binWidth_in_ps*0.5 
        
        for i in range(0, numberOfComponents):
            areaInitial[i] += (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-0.5*((timeBin_in_ps[bin]-tZero_in_ps)/sigma)**2)
        
    for i in range(0, numberOfComponents):
        areaInitial[i] *= intensitiesOfGaussian[i]*numberOfComponents
            
    for bin in range(0, numberOfBins):
        for i in range(0, numberOfComponents):
            counts_y[bin] += (countsInitial[i]/areaInitial[i])*(1/(sigma*np.sqrt(2*np.pi)))*np.exp(-0.5*((timeBin_in_ps[bin]-tZero_in_ps)/sigma)**2)      
            
        sumOfCounts += (int)(counts_y[bin])
     
    return counts_y

#fit the IRF model function on data (xIRF, yIRF):
if userInput.__bUsingModel:
    print("fitting irf model...")
    
    #Gaussian:
    if userInput.__modelType == functionModelList.ReconvolutionModel.Gaussian:
        fitModelIRF = Model(functionModelList.Gaussian);
        fitModelIRF.set_param_hint('sigma', min=0.0);
        fitModelIRF.set_param_hint('ampl', min=0.0);
        fitModelIRF.set_param_hint('y0', min=0.0);

        parameterListIRFFit = fitModelIRF.make_params(x=xVal, ampl=yMaxIRF, sigma=stddevIRF, y0=estimatedBkgrdIRF, x0=xWhereYMaxIRF, args=yIRF);
        #change here if you want to fix x0 and/or y0:
        parameterListIRFFit['x0'].vary = True; 
        parameterListIRFFit['y0'].vary = True;

        #run the fit:
        resultsOfModelIRF = fitModelIRF.fit(yIRF, params=parameterListIRFFit, weights=fitWeightingIRF, method='leastsq', x=xVal);

        chiSquare = resultsOfModelIRF.chisqr;
        redChiSquare = resultsOfModelIRF.redchi;

        sigma = (float)(resultsOfModelIRF.params['sigma'].value*userInput.__channelResolutionInPs);
        sigma_err = (float)(resultsOfModelIRF.params['sigma'].stderr*userInput.__channelResolutionInPs);

        amplitude = (float)(resultsOfModelIRF.params['ampl'].value);
        amplitude_err = (float)(resultsOfModelIRF.params['ampl'].stderr);

        yRes = (float)(resultsOfModelIRF.params['y0'].value);
        yRes_err = (float)(resultsOfModelIRF.params['y0'].stderr);


        print("\nFit results: IRF - Model Function (Type: Gaussian):");
        print("-----------------------------------------------------");
        print("X²          = {0}".format(redChiSquare));
        print("");
        print("stddev [ps] = {0} ({1})".format(sigma, sigma_err));
        print("FWHM   [ps] = {0} ({1})".format(sigma*(2*np.sqrt(2*np.log(2))), sigma_err*(2*np.sqrt(2*np.log(2)))));
        print("");
        print("amplitude   = {0} ({1})".format(amplitude, amplitude_err));
        print("");
        print("background  = {0} ({1})".format(yRes, yRes_err));
        print("-----------------------------------------------------");


        plt.figure(1);
        ax = plt.subplot(2,1,1);
        ax.set_title("Best fit: IRF Gaussian model fit");
        plt.semilogy(xVal, yIRF,'o', xVal , resultsOfModelIRF.best_fit, 'b');
        ax2 = plt.subplot(2,1,2);
        ax2.set_title("Best fit: residuals");
        plt.plot(xVal, resultsOfModelIRF.residual);
        
        plt.figure(2)
        ax3 = plt.subplot(1,1,1);
        ax3.set_title("lifetime spectrum + IRF + IRF best fit");
        plt.semilogy(xVal, ySpec, xVal, yIRF, xVal, resultsOfModelIRF.best_fit)
        
        #replace input by model fit data:
        yIRF = resultsOfModelIRF.best_fit; 


    #Lorentz/Cauchy:
    if userInput.__modelType == functionModelList.ReconvolutionModel.Lorentz_Cauchy:
        fitModelIRF = Model(functionModelList.Lorentz_Cauchy);
        fitModelIRF.set_param_hint('s', min=0.0);
        fitModelIRF.set_param_hint('ampl', min=0.0);
        fitModelIRF.set_param_hint('y0', min=0.0);
        
        parameterListIRFFit = fitModelIRF.make_params(x=xVal, ampl=yMaxIRF, s=stddevIRF, y0=estimatedBkgrdIRF, x0=xWhereYMaxIRF, args=yIRF);
        #change here if you want to fix x0 and/or y0:
        parameterListIRFFit['x0'].vary = True; 
        parameterListIRFFit['y0'].vary = True;

        #run the fit:
        resultsOfModelIRF = fitModelIRF.fit(yIRF, params=parameterListIRFFit, weights=fitWeightingIRF, method='leastsq', x=xVal);

        chiSquare = resultsOfModelIRF.chisqr;
        redChiSquare = resultsOfModelIRF.redchi;

        s = (float)(resultsOfModelIRF.params['s'].value*userInput.__channelResolutionInPs);
        s_err = (float)(resultsOfModelIRF.params['s'].stderr*userInput.__channelResolutionInPs);

        amplitude = (float)(resultsOfModelIRF.params['ampl'].value);
        amplitude_err = (float)(resultsOfModelIRF.params['ampl'].stderr);

        yRes = (float)(resultsOfModelIRF.params['y0'].value);
        yRes_err = (float)(resultsOfModelIRF.params['y0'].stderr);


        print("\nFit results: IRF - Model Function (Type: Lorentzian/Cauchy):");
        print("--------------------------------------------------------------");
        print("X²          = {0}".format(redChiSquare));
        print("");
        print("s [ps]      = {0} ({1})".format(s, s_err));
        print("");
        print("amplitude   = {0} ({1})".format(amplitude, amplitude_err));
        print("");
        print("background  = {0} ({1})".format(yRes, yRes_err));
        print("--------------------------------------------------------------");

        
        plt.figure(1);
        ax = plt.subplot(2,1,1);
        ax.set_title("Best fit: IRF Lorentzian/Cauchy model fit");
        plt.semilogy(xVal, yIRF,'o', xVal , resultsOfModelIRF.best_fit, 'b');
        ax2 = plt.subplot(2,1,2);
        ax2.set_title("Best fit: residuals");
        plt.plot(xVal, resultsOfModelIRF.residual);
        
        plt.figure(2)
        ax3 = plt.subplot(1,1,1);
        ax3.set_title("lifetime spectrum + IRF + IRF best fit");
        plt.semilogy(xVal, ySpec, xVal, yIRF, xVal, resultsOfModelIRF.best_fit)
        
        #replace input by model fit data:
        yIRF = resultsOfModelIRF.best_fit; 


    #PseudoVoigt-1:
    if userInput.__modelType == functionModelList.ReconvolutionModel.Pseudovoigt1:
        fitModelIRF = Model(functionModelList.Pseudovoigt1);
        fitModelIRF.set_param_hint('a', min=0.0, max=1.0);
        fitModelIRF.set_param_hint('sigma', min=0.0);
        fitModelIRF.set_param_hint('s', min=0.0);
        fitModelIRF.set_param_hint('ampl', min=0.0);
        fitModelIRF.set_param_hint('y0', min=0.0);

        parameterListIRFFit = fitModelIRF.make_params(x=xVal, ampl=yMaxIRF, a=0.8, sigma=stddevIRF, s=stddevIRF*0.5, y0=estimatedBkgrdIRF, x0=xWhereYMaxIRF, args=yIRF); 
        #change here if you want to fix x0 and/or y0:
        parameterListIRFFit['x0'].vary = True; 
        parameterListIRFFit['y0'].vary = True;

        #run the fit:
        resultsOfModelIRF = fitModelIRF.fit(yIRF, params=parameterListIRFFit, weights=fitWeightingIRF, method='leastsq', x=xVal);

        chiSquare = resultsOfModelIRF.chisqr;
        redChiSquare = resultsOfModelIRF.redchi;

        fract = (float)(resultsOfModelIRF.params['a'].value);
        fract_err = (float)(resultsOfModelIRF.params['a'].stderr);
        
        sigma = (float)(resultsOfModelIRF.params['sigma'].value*userInput.__channelResolutionInPs);
        sigma_err = (float)(resultsOfModelIRF.params['sigma'].stderr*userInput.__channelResolutionInPs);
        
        s = (float)(resultsOfModelIRF.params['s'].value*userInput.__channelResolutionInPs);
        s_err = (float)(resultsOfModelIRF.params['s'].stderr*userInput.__channelResolutionInPs);

        amplitude = (float)(resultsOfModelIRF.params['ampl'].value);
        amplitude_err = (float)(resultsOfModelIRF.params['ampl'].stderr);

        yRes = (float)(resultsOfModelIRF.params['y0'].value);
        yRes_err = (float)(resultsOfModelIRF.params['y0'].stderr);


        print("\nFit results: IRF - Model Function (Type: PseudoVoigt type 1):");
        print("---------------------------------------------------------------");
        print("X²          = {0}".format(redChiSquare));
        print("");
        print("amplitude   = {0} ({1})".format(amplitude, amplitude_err));
        print("background  = {0} ({1})".format(yRes, yRes_err));
        print("");
        print("---------------------------------------------------------------");
        print("G - Gaussian:  a        = {0} ({1})".format(fract, fract_err));
        print("---------------------------------------------------------------");
        print("stddev [ps] = {0} ({1})".format(sigma, sigma_err));
        print("FWHM   [ps] = {0} ({1})".format(sigma*(2*np.sqrt(2*np.log(2))), sigma_err*(2*np.sqrt(2*np.log(2)))));
        print("");
        print("---------------------------------------------------------------");
        print("L - Lorentzian: (1 - a) = {0} ({1})".format(1-fract, fract_err));
        print("---------------------------------------------------------------");
        print("s [ps]      = {0} ({1})".format(s, s_err));
        print("---------------------------------------------------------------");


        plt.figure(1);
        ax = plt.subplot(2,1,1);
        ax.set_title("Best fit: IRF Pseudo-Voigt 1 model fit");
        plt.semilogy(xVal, yIRF,'o', xVal , resultsOfModelIRF.best_fit, 'b');
        ax2 = plt.subplot(2,1,2);
        ax2.set_title("Best fit: residuals");
        plt.plot(xVal, resultsOfModelIRF.residual);
        
        plt.figure(2)
        ax3 = plt.subplot(1,1,1);
        ax3.set_title("lifetime spectrum + IRF + IRF best fit");
        plt.semilogy(xVal, ySpec, xVal, yIRF, xVal, resultsOfModelIRF.best_fit)
        
        #replace input by model fit data:
        yIRF = resultsOfModelIRF.best_fit; 

        
    #Pearson Type VII:
    if userInput.__modelType == functionModelList.ReconvolutionModel.Pearson7:
        fitModelIRF = Model(functionModelList.Pearson7);
        fitModelIRF.set_param_hint('alpha', min=0.0);
        fitModelIRF.set_param_hint('m', min=0.1);
        fitModelIRF.set_param_hint('ampl', min=0.0);
        fitModelIRF.set_param_hint('y0', min=0.0);

        
        parameterListIRFFit = fitModelIRF.make_params(x=xVal, ampl=yMaxIRF, alpha=stddevIRF, m=2, y0=estimatedBkgrdIRF, x0=xWhereYMaxIRF, args=yIRF);
        #change here if you want to fix x0 and/or y0:
        parameterListIRFFit['x0'].vary = True; 
        parameterListIRFFit['y0'].vary = True;

        #run the fit:
        resultsOfModelIRF = fitModelIRF.fit(yIRF, params=parameterListIRFFit, weights=fitWeightingIRF, method='leastsq', x=xVal);

        chiSquare = resultsOfModelIRF.chisqr;
        redChiSquare = resultsOfModelIRF.redchi;

        alpha = (float)(resultsOfModelIRF.params['alpha'].value*userInput.__channelResolutionInPs);
        alpha_err = (float)(resultsOfModelIRF.params['alpha'].stderr*userInput.__channelResolutionInPs);
        
        m = (float)(resultsOfModelIRF.params['m'].value);
        m_err = (float)(resultsOfModelIRF.params['m'].stderr);
        
        amplitude = (float)(resultsOfModelIRF.params['ampl'].value);
        amplitude_err = (float)(resultsOfModelIRF.params['ampl'].stderr);

        yRes = (float)(resultsOfModelIRF.params['y0'].value);
        yRes_err = (float)(resultsOfModelIRF.params['y0'].stderr);


        print("\nFit results: IRF - Model Function (Type: Pearson type 7):");
        print("-----------------------------------------------------------");
        print("X²          = {0}".format(redChiSquare));
        print("");
        print("alpha [ps]  = {0} ({1})".format(alpha, alpha_err));
        print("m           = {0} ({1})".format(m, m_err));
        print("");
        print("amplitude   = {0} ({1})".format(amplitude, amplitude_err));
        print("");
        print("background  = {0} ({1})".format(yRes, yRes_err));
        print("-----------------------------------------------------------");
        

        plt.figure(1);
        ax = plt.subplot(2,1,1);
        ax.set_title("Best fit: IRF Pearson Type VII model fit");
        plt.semilogy(xVal, yIRF,'o', xVal , resultsOfModelIRF.best_fit, 'b');
        ax2 = plt.subplot(2,1,2);
        ax2.set_title("Best fit: residuals");
        plt.plot(xVal, resultsOfModelIRF.residual);
        
        plt.figure(2)
        ax3 = plt.subplot(1,1,1);
        ax3.set_title("lifetime spectrum + IRF + IRF best fit");
        plt.semilogy(xVal, ySpec, xVal, yIRF, xVal, resultsOfModelIRF.best_fit)
        
        #replace input by model fit data:
        yIRF = resultsOfModelIRF.best_fit; 

#1 component expontential distribution (with convoluted Gaussian kernel):
def ExpDecay_1_conv(x, ampl1, tau1, y0, x0, addKernelFWHM, args=(yIRF)):  
    h = np.zeros(x.size) 
    lengthVec = len(ySpec)
    
    irf_norm        = yIRF/sum(yIRF)

    x               = np.arange(0, lengthVec, 1)
    centerOfMass    = np.argmax(x*yIRF/sum(yIRF))
    
    sigma           = addKernelFWHM/(2*np.sqrt(2*np.log(2)))
    kernelData      = norm.pdf(x, 0.5*centerOfMass, sigma/userInput.__channelResolutionInPs) 
    kernelData_norm = kernelData/sum(kernelData)
    
    convK           = convolveData(irf_norm, kernelData_norm)
    centerOfMass_c  = np.argmax(x*convK/sum(convK))

    convK_norm      = np.roll(convK, -centerOfMass_c + centerOfMass)

    shift_1 = np.remainder(np.remainder(x-np.floor(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr1 = (1 - x0 + np.floor(x0))*convK_norm[shift_1.astype(int)]
    
    shift_2 = np.remainder(np.remainder(x-np.ceil(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr2 = (x0 - np.floor(x0))*convK_norm[shift_2.astype(int)]

    irf_shifted = (shift_Incr1 + shift_Incr2)
    convK_norm = irf_shifted/sum(irf_shifted)

    #plt.semilogy(irf_norm, 'go', kernelData_norm, 'r', convK_norm, 'b');
    #plt.show()

    h = ampl1*np.exp(-(x)/tau1) 
    hConvIrf_norm = convolveData(h, convK_norm)

    return hConvIrf_norm + y0

#1 component expontential distribution:
def ExpDecay_1(x, ampl1, tau1, y0, x0, args=(yIRF)):
    h = np.zeros(x.size) 
    lengthVec = len(ySpec)

    shift_1 = np.remainder(np.remainder(x-np.floor(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr1 = (1 - x0 + np.floor(x0))*yIRF[shift_1.astype(int)]
    
    shift_2 = np.remainder(np.remainder(x-np.ceil(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr2 = (x0 - np.floor(x0))*yIRF[shift_2.astype(int)]
    
    irf_shifted = (shift_Incr1 + shift_Incr2)
    irf_norm = irf_shifted/sum(irf_shifted)
    
    h = ampl1*np.exp(-(x)/tau1)
    hConvIrf_norm = convolveData(h, irf_norm)
    return hConvIrf_norm + y0

#2 component expontential distribution (with convoluted Gaussian kernel):
def ExpDecay_2_conv(x, ampl1, tau1, ampl2, tau2, y0, x0, addKernelFWHM, args=(yIRF)):  
    h = np.zeros(x.size) 
    lengthVec = len(ySpec)
    
    irf_norm        = yIRF/sum(yIRF)

    x               = np.arange(0, lengthVec, 1)
    centerOfMass    = np.argmax(x*yIRF/sum(yIRF))
    
    sigma           = addKernelFWHM/(2*np.sqrt(2*np.log(2)))
    kernelData      = norm.pdf(x, 0.5*centerOfMass, sigma/userInput.__channelResolutionInPs) 
    kernelData_norm = kernelData/sum(kernelData)
    
    convK           = convolveData(irf_norm, kernelData_norm)
    centerOfMass_c  = np.argmax(x*convK/sum(convK))

    convK_norm      = np.roll(convK, -centerOfMass_c + centerOfMass)

    shift_1 = np.remainder(np.remainder(x-np.floor(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr1 = (1 - x0 + np.floor(x0))*convK_norm[shift_1.astype(int)]
    
    shift_2 = np.remainder(np.remainder(x-np.ceil(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr2 = (x0 - np.floor(x0))*convK_norm[shift_2.astype(int)]

    irf_shifted = (shift_Incr1 + shift_Incr2)
    convK_norm = irf_shifted/sum(irf_shifted)

    #plt.semilogy(irf_norm, 'go', kernelData_norm, 'r', convK_norm, 'b');
    #plt.show()

    h = ampl1*np.exp(-(x)/tau1) + ampl2*np.exp(-(x)/tau2) 
    hConvIrf_norm = convolveData(h, convK_norm)

    return hConvIrf_norm + y0

#2 component expontential distribution:
def ExpDecay_2(x, ampl1, tau1, ampl2, tau2, y0, x0, args=(yIRF)):  
    h = np.zeros(x.size) 
    lengthVec = len(ySpec)

    shift_1 = np.remainder(np.remainder(x-np.floor(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr1 = (1 - x0 + np.floor(x0))*yIRF[shift_1.astype(int)]
    
    shift_2 = np.remainder(np.remainder(x-np.ceil(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr2 = (x0 - np.floor(x0))*yIRF[shift_2.astype(int)]
    
    irf_shifted = (shift_Incr1 + shift_Incr2)
    irf_norm = irf_shifted/sum(irf_shifted)

    h = ampl1*np.exp(-(x)/tau1) + ampl2*np.exp(-(x)/tau2);
    hConvIrf_norm = convolveData(h, irf_norm)
    return hConvIrf_norm + y0

#3 component expontential distribution (with convoluted Gaussian kernel):
def ExpDecay_3_conv(x, ampl1, tau1, ampl2, tau2, ampl3, tau3, y0, x0, addKernelFWHM, args=(yIRF)):  
    h = np.zeros(x.size) 
    lengthVec = len(ySpec)
    
    irf_norm        = yIRF/sum(yIRF)

    x               = np.arange(0, lengthVec, 1)
    centerOfMass    = np.argmax(x*yIRF/sum(yIRF))
    
    sigma           = addKernelFWHM/(2*np.sqrt(2*np.log(2)))
    kernelData      = norm.pdf(x, 0.5*centerOfMass, sigma/userInput.__channelResolutionInPs) 
    kernelData_norm = kernelData/sum(kernelData)
    
    convK           = convolveData(irf_norm, kernelData_norm)
    centerOfMass_c  = np.argmax(x*convK/sum(convK))

    convK_norm      = np.roll(convK, -centerOfMass_c + centerOfMass)

    shift_1 = np.remainder(np.remainder(x-np.floor(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr1 = (1 - x0 + np.floor(x0))*convK_norm[shift_1.astype(int)]
    
    shift_2 = np.remainder(np.remainder(x-np.ceil(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr2 = (x0 - np.floor(x0))*convK_norm[shift_2.astype(int)]

    irf_shifted = (shift_Incr1 + shift_Incr2)
    convK_norm = irf_shifted/sum(irf_shifted)

    #plt.semilogy(irf_norm, 'go', kernelData_norm, 'r', convK_norm, 'b');
    #plt.show()

    h = ampl1*np.exp(-(x)/tau1) + ampl2*np.exp(-(x)/tau2) + ampl3*np.exp(-(x)/tau3) 
    hConvIrf_norm = convolveData(h, convK_norm)

    return hConvIrf_norm + y0

#3 component expontential distribution
def ExpDecay_3(x, ampl1, tau1, ampl2, tau2, ampl3, tau3, y0, x0, args=(yIRF)):  
    h = np.zeros(x.size) 
    lengthVec = len(ySpec)
    
    shift_1 = np.remainder(np.remainder(x-np.floor(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr1 = (1 - x0 + np.floor(x0))*yIRF[shift_1.astype(int)]
    
    shift_2 = np.remainder(np.remainder(x-np.ceil(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr2 = (x0 - np.floor(x0))*yIRF[shift_2.astype(int)]
    
    irf_shifted = (shift_Incr1 + shift_Incr2)
    irf_norm = irf_shifted/sum(irf_shifted)

    h = ampl1*np.exp(-(x)/tau1) + ampl2*np.exp(-(x)/tau2) + ampl3*np.exp(-(x)/tau3) 
    hConvIrf_norm = convolveData(h, irf_norm)
    
    return hConvIrf_norm + y0

#4 component expontential distribution (with convoluted Gaussian kernel):
def ExpDecay_4_conv(x, ampl1, tau1, ampl2, tau2, ampl3, tau3, ampl4, tau4, y0, x0, addKernelFWHM, args=(yIRF)):  
    h = np.zeros(x.size) 
    lengthVec = len(ySpec)
    
    irf_norm        = yIRF/sum(yIRF)

    x               = np.arange(0, lengthVec, 1)
    centerOfMass    = np.argmax(x*yIRF/sum(yIRF))
    
    sigma           = addKernelFWHM/(2*np.sqrt(2*np.log(2)))
    kernelData      = norm.pdf(x, 0.5*centerOfMass, sigma/userInput.__channelResolutionInPs) 
    kernelData_norm = kernelData/sum(kernelData)
    
    convK           = convolveData(irf_norm, kernelData_norm)
    centerOfMass_c  = np.argmax(x*convK/sum(convK))

    convK_norm      = np.roll(convK, -centerOfMass_c + centerOfMass)

    shift_1 = np.remainder(np.remainder(x-np.floor(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr1 = (1 - x0 + np.floor(x0))*convK_norm[shift_1.astype(int)]
    
    shift_2 = np.remainder(np.remainder(x-np.ceil(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr2 = (x0 - np.floor(x0))*convK_norm[shift_2.astype(int)]

    irf_shifted = (shift_Incr1 + shift_Incr2)
    convK_norm = irf_shifted/sum(irf_shifted)

    #plt.semilogy(irf_norm, 'go', kernelData_norm, 'r', convK_norm, 'b');
    #plt.show()

    h = ampl1*np.exp(-(x)/tau1) + ampl2*np.exp(-(x)/tau2) + ampl3*np.exp(-(x)/tau3) + ampl4*np.exp(-(x)/tau4) 
    hConvIrf_norm = convolveData(h, convK_norm)

    return hConvIrf_norm + y0

#4 component expontential distribution:
def ExpDecay_4(x, ampl1, tau1, ampl2, tau2, ampl3, tau3, ampl4, tau4, y0, x0, args=(yIRF)):  
    h = np.zeros(x.size) 
    lengthVec = len(ySpec)
    
    shift_1 = np.remainder(np.remainder(x-np.floor(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr1 = (1 - x0 + np.floor(x0))*yIRF[shift_1.astype(int)]
    
    shift_2 = np.remainder(np.remainder(x-np.ceil(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr2 = (x0 - np.floor(x0))*yIRF[shift_2.astype(int)]
    
    irf_shifted = (shift_Incr1 + shift_Incr2)
    irf_norm = irf_shifted/sum(irf_shifted)
    
    h = ampl1*np.exp(-(x)/tau1) + ampl2*np.exp(-(x)/tau2) + ampl3*np.exp(-(x)/tau3) + ampl4*np.exp(-(x)/tau4)
    hConvIrf_norm = convolveData(h, irf_norm)
    return hConvIrf_norm + y0

#applying reconvolution:
if userInput.__numberOfExpDec == 1:
    print("\n\nreconvolution fitting with 1 component...\n")

    if userInput.__bUsingAdditionalGaussianKernel:
        fitModelDecay = Model(ExpDecay_1_conv);
    else:
        fitModelDecay = Model(ExpDecay_1)
        
    fitModelDecay.set_param_hint('ampl1', min=0.0);
    fitModelDecay.set_param_hint('tau1', min=0.00001);
    fitModelDecay.set_param_hint('y0', min=0.0);
    
    # convolve IRF data with an additional Gaussian kernel for producing artificial broadening:
    if userInput.__bUsingAdditionalGaussianKernel:
        fitModelDecay.set_param_hint('addKernelFWHM', min=0.01);
            
        parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=(userInput.__expectedTau_1_in_ps/userInput.__channelResolutionInPs), y0=estimatedBkgrd, x0=xWhereYMaxSpec, addKernelFWHM=userInput.__gaussianKernelFWHM, args=ySpec)
        parameterListDecayFit['addKernelFWHM'].vary = userInput.__bVaryGaussianKernelFWHM;
    else:
        parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=(userInput.__expectedTau_1_in_ps/userInput.__channelResolutionInPs), y0=estimatedBkgrd, x0=xWhereYMaxSpec, args=ySpec)
    
    #change here if you want to fix x0 and/or y0:
    parameterListDecayFit['x0'].vary = True; 
    parameterListDecayFit['y0'].vary = not userInput.__bkgrdFixed;

    #run the fit:
    resultsOfModelDecay = fitModelDecay.fit(ySpec, params=parameterListDecayFit, weights=fitWeightingSpec, method='leastsq', x=xVal);

    #calculate results:
    chiSquare = resultsOfModelDecay.chisqr;
    redChiSquare = resultsOfModelDecay.redchi;
    
    t1 = (float)(resultsOfModelDecay.params['tau1'].value*userInput.__channelResolutionInPs);
    t1_err = (float)(resultsOfModelDecay.params['tau1'].stderr*userInput.__channelResolutionInPs);

    yRes = (float)(resultsOfModelDecay.params['y0'].value);
    yRes_err = (float)(resultsOfModelDecay.params['y0'].stderr);

    amplitude1 = (float)(resultsOfModelDecay.params['ampl1'].value);
    amplitude1_err = (float)(resultsOfModelDecay.params['ampl1'].stderr);
    
    counts_1 = 0;
    
    counts_1_stddev = 0;
    
    for i in range(0, len(xVal)):
        counts_1 += amplitude1*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t1) 
        
        counts_1_stddev += np.exp(-2*xVal[i]*userInput.__channelResolutionInPs/t1)*(amplitude1_err**2 + (t1_err**2/t1**4)) + yRes_err**2
        
    counts_1_err = np.sqrt(counts_1_stddev)/np.sqrt(len(xVal));
    
    counts_sum = (counts_1);
    counts_sum_err = np.sqrt(counts_1_err*counts_1_err);

    I1 = (counts_1/counts_sum)
    I1_err = np.sqrt((1/counts_sum)**2*counts_1_err**2 + counts_1**2*counts_sum_err**2/counts_sum**4); 

    if userInput.__bUsingAdditionalGaussianKernel:
            kernelFWHM = (float)(resultsOfModelDecay.params['addKernelFWHM'].value);
            kernelFWHM_err = (float)(resultsOfModelDecay.params['addKernelFWHM'].stderr);
    
            str1  = "Fit results: Reconvolution 1 component (with convoluted Gaussian kernel):\n"
            str2  = "--------------------------------------------------------------------------\n"
            str3  = "X²                   = {0}\n\n".format(redChiSquare)
            str5  = "tau  (1)             = {0} ({1})\n".format(t1, t1_err)
            str6  = "I    (1)             = {0} ({1})\n".format(I1, I1_err)
            str14 = "background           = {0} ({1})\n".format(yRes, yRes_err)
            str15 = "\n"
            str16 = "conv. kernel (FWHM)  = {0} ({1})\n".format(kernelFWHM, kernelFWHM_err)
            str17 = "--------------------------------------------------------------------------"

            resStr = str1+str2+str3+str5+str6+str14+str15+str16+str17
            print(resStr)

            if userInput.__saveReconvolutionResults:
                resultsFile = open(userInput.__saveReconvolutionResultsPath,"a")
                resultsFile.write(resStr) 
                resultsFile.close()
    else:
        str1 = "Fit results: Reconvolution 1 component:\n"
        str2 = "----------------------------------------\n"
        str3 = "X²                   = {0}\n\n".format(redChiSquare)
        str5 = "tau  (1)             = {0} ({1})\n".format(t1, t1_err)
        str6 = "I    (1)             = {0} ({1})\n".format(I1, I1_err)
        str14 = "background          = {0} ({1})\n".format(yRes, yRes_err)
        str17 = "----------------------------------------"

        resStr = str1+str2+str3+str5+str6+str14+str17
        print(resStr)

        if userInput.__saveReconvolutionResults:
            resultsFile = open(userInput.__saveReconvolutionResultsPath,"a")
            resultsFile.write(resStr) 
            resultsFile.close()

    plt.figure(3)
    ax = plt.subplot(2,1,1);
    ax.set_title("Best fit: Reconvolution with 1 lifetime component");
    plt.semilogy(xVal, ySpec,'o', xVal ,resultsOfModelDecay.best_fit, 'b');
    ax2 = plt.subplot(2,1,2);
    ax2.set_title("Best fit: Residuals");
    plt.plot(xVal, resultsOfModelDecay.residual);

if userInput.__numberOfExpDec == 2:
    print("\n\nreconvolution fitting with 2 components...\n")

    if userInput.__bUsingAdditionalGaussianKernel:
        fitModelDecay = Model(ExpDecay_2_conv);
    else:
        fitModelDecay = Model(ExpDecay_2)
        
    fitModelDecay.set_param_hint('ampl1', min=0.0);
    fitModelDecay.set_param_hint('tau1', min=0.00001);
    fitModelDecay.set_param_hint('ampl2', min=0.0);
    fitModelDecay.set_param_hint('tau2', min=0.00001);
    fitModelDecay.set_param_hint('y0', min=0.0);
    
    # convolve IRF data with an additional Gaussian kernel for producing artificial broadening:
    if userInput.__bUsingAdditionalGaussianKernel:
        fitModelDecay.set_param_hint('addKernelFWHM', min=0.01);
            
        parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=(userInput.__expectedTau_1_in_ps/userInput.__channelResolutionInPs), ampl2=yMaxSpec, tau2=(userInput.__expectedTau_2_in_ps/userInput.__channelResolutionInPs), y0=estimatedBkgrd, x0=xWhereYMaxSpec, addKernelFWHM=userInput.__gaussianKernelFWHM, args=ySpec)
        parameterListDecayFit['addKernelFWHM'].vary = userInput.__bVaryGaussianKernelFWHM;
    else:
        parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=(userInput.__expectedTau_1_in_ps/userInput.__channelResolutionInPs), ampl2=yMaxSpec, tau2=(userInput.__expectedTau_2_in_ps/userInput.__channelResolutionInPs), y0=estimatedBkgrd, x0=xWhereYMaxSpec, args=ySpec)
    
    #change here if you want to fix x0 and/or y0:
    parameterListDecayFit['x0'].vary = True; 
    parameterListDecayFit['y0'].vary = not userInput.__bkgrdFixed;

    #run the fit:
    resultsOfModelDecay = fitModelDecay.fit(ySpec, params=parameterListDecayFit, weights=fitWeightingSpec, method='leastsq', x=xVal);

    #calculate results:
    chiSquare = resultsOfModelDecay.chisqr;
    redChiSquare = resultsOfModelDecay.redchi;
    
    t1 = (float)(resultsOfModelDecay.params['tau1'].value*userInput.__channelResolutionInPs);
    t1_err = (float)(resultsOfModelDecay.params['tau1'].stderr*userInput.__channelResolutionInPs);

    t2 = (float)(resultsOfModelDecay.params['tau2'].value*userInput.__channelResolutionInPs);
    t2_err = (float)(resultsOfModelDecay.params['tau2'].stderr*userInput.__channelResolutionInPs);

    yRes = (float)(resultsOfModelDecay.params['y0'].value);
    yRes_err = (float)(resultsOfModelDecay.params['y0'].stderr);

    amplitude1 = (float)(resultsOfModelDecay.params['ampl1'].value);
    amplitude1_err = (float)(resultsOfModelDecay.params['ampl1'].stderr);
    amplitude2 = (float)(resultsOfModelDecay.params['ampl2'].value);
    amplitude2_err = (float)(resultsOfModelDecay.params['ampl2'].stderr);
    
    counts_1 = 0;
    counts_2 = 0;
    
    counts_1_stddev = 0;
    counts_2_stddev = 0;
    
    for i in range(0, len(xVal)):
        counts_1 += amplitude1*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t1) 
        counts_2 += amplitude2*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t2) 
       
        counts_1_stddev += np.exp(-2*xVal[i]*userInput.__channelResolutionInPs/t1)*(amplitude1_err**2 + (t1_err**2/t1**4)) + yRes_err**2
        counts_2_stddev += np.exp(-2*xVal[i]*userInput.__channelResolutionInPs/t2)*(amplitude2_err**2 + (t2_err**2/t2**4)) + yRes_err**2
        
    counts_1_err = np.sqrt(counts_1_stddev)/np.sqrt(len(xVal));
    counts_2_err = np.sqrt(counts_2_stddev)/np.sqrt(len(xVal));
    
    counts_sum = (counts_1 + counts_2);
    counts_sum_err = np.sqrt(counts_1_err*counts_1_err + counts_2_err*counts_2_err);

    I1 = (counts_1/counts_sum)
    I1_err = np.sqrt((1/counts_sum)**2*counts_1_err**2 + counts_1**2*counts_sum_err**2/counts_sum**4); 

    I2 = (counts_2/counts_sum)
    I2_err = np.sqrt((1/counts_sum)**2*counts_2_err**2 + counts_2**2*counts_sum_err**2/counts_sum**4);

    if userInput.__bUsingAdditionalGaussianKernel:
            kernelFWHM = (float)(resultsOfModelDecay.params['addKernelFWHM'].value);
            kernelFWHM_err = (float)(resultsOfModelDecay.params['addKernelFWHM'].stderr);
    
            str1  = "Fit results: Reconvolution 2 components (with convoluted Gaussian kernel):\n"
            str2  = "--------------------------------------------------------------------------\n"
            str3  = "X²                   = {0}\n\n".format(redChiSquare)
            str5  = "tau  (1)             = {0} ({1})\n".format(t1, t1_err)
            str6  = "I    (1)             = {0} ({1})\n".format(I1, I1_err)
            str8  = "tau  (2)             = {0} ({1})\n".format(t2, t2_err)
            str9  = "I    (2)             = {0} ({1})\n".format(I2, I2_err)
            str14 = "background          = {0} ({1})\n".format(yRes, yRes_err)
            str15 = "\n"
            str16 = "conv. kernel (FWHM) = {0} ({1})\n".format(kernelFWHM, kernelFWHM_err)
            str17 = "--------------------------------------------------------------------------"

            resStr = str1+str2+str3+str5+str6+str8+str9+str14+str15+str16+str17
            print(resStr)

            if userInput.__saveReconvolutionResults:
                resultsFile = open(userInput.__saveReconvolutionResultsPath,"a")
                resultsFile.write(resStr) 
                resultsFile.close()
    else:
        str1 = "Fit results: Reconvolution 2 components:\n"
        str2 = "----------------------------------------\n"
        str3 = "X²                   = {0}\n\n".format(redChiSquare)
        str5 = "tau  (1)             = {0} ({1})\n".format(t1, t1_err)
        str6 = "I    (1)             = {0} ({1})\n".format(I1, I1_err)
        str8 = "tau  (2)             = {0} ({1})\n".format(t2, t2_err)
        str9 = "I    (2)             = {0} ({1})\n".format(I2, I2_err)
        str14 = "background          = {0} ({1})\n".format(yRes, yRes_err)
        str17 = "----------------------------------------"

        resStr = str1+str2+str3+str5+str6+str8+str9+str14+str17
        print(resStr)

        if userInput.__saveReconvolutionResults:
            resultsFile = open(userInput.__saveReconvolutionResultsPath,"a")
            resultsFile.write(resStr) 
            resultsFile.close()

    plt.figure(3)
    ax = plt.subplot(2,1,1);
    ax.set_title("Best fit: Reconvolution with 2 lifetime components");
    plt.semilogy(xVal, ySpec,'o', xVal ,resultsOfModelDecay.best_fit, 'b');
    ax2 = plt.subplot(2,1,2);
    ax2.set_title("Best fit: Residuals");
    plt.plot(xVal, resultsOfModelDecay.residual);

if userInput.__numberOfExpDec == 3:
    print("\n\nreconvolution fitting with 3 components...\n")

    if userInput.__bUsingAdditionalGaussianKernel:
        fitModelDecay = Model(ExpDecay_3_conv);
    else:
        fitModelDecay = Model(ExpDecay_3)
        
    fitModelDecay.set_param_hint('ampl1', min=0.0);
    fitModelDecay.set_param_hint('tau1', min=0.00001);
    fitModelDecay.set_param_hint('ampl2', min=0.0);
    fitModelDecay.set_param_hint('tau2', min=0.00001);
    fitModelDecay.set_param_hint('ampl3', min=0.0);
    fitModelDecay.set_param_hint('tau3', min=0.00001);
    fitModelDecay.set_param_hint('y0', min=0.0);
    
    # (de)convolve IRF data with an additional Gaussian kernel for producing artificial broadening:
    if userInput.__bUsingAdditionalGaussianKernel:
        fitModelDecay.set_param_hint('addKernelFWHM', min=0.01);
            
        parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=(userInput.__expectedTau_1_in_ps/userInput.__channelResolutionInPs), ampl2=yMaxSpec, tau2=(userInput.__expectedTau_2_in_ps/userInput.__channelResolutionInPs), ampl3=yMaxSpec, tau3=(userInput.__expectedTau_3_in_ps/userInput.__channelResolutionInPs), y0=estimatedBkgrd, x0=xWhereYMaxSpec, addKernelFWHM=userInput.__gaussianKernelFWHM, args=ySpec)
        parameterListDecayFit['addKernelFWHM'].vary = userInput.__bVaryGaussianKernelFWHM;
    else:
        parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=(userInput.__expectedTau_1_in_ps/userInput.__channelResolutionInPs), ampl2=yMaxSpec, tau2=(userInput.__expectedTau_2_in_ps/userInput.__channelResolutionInPs), ampl3=yMaxSpec, tau3=(userInput.__expectedTau_3_in_ps/userInput.__channelResolutionInPs), y0=estimatedBkgrd, x0=xWhereYMaxSpec, args=ySpec)
    
    #change here if you want to fix x0 and/or y0:
    parameterListDecayFit['x0'].vary = True; 
    parameterListDecayFit['y0'].vary = not userInput.__bkgrdFixed;

    #run the fit:
    resultsOfModelDecay = fitModelDecay.fit(ySpec, params=parameterListDecayFit, weights=fitWeightingSpec, method='leastsq', x=xVal);

    #calculate results:
    chiSquare = resultsOfModelDecay.chisqr;
    redChiSquare = resultsOfModelDecay.redchi;
    
    t1 = (float)(resultsOfModelDecay.params['tau1'].value*userInput.__channelResolutionInPs);
    t1_err = (float)(resultsOfModelDecay.params['tau1'].stderr*userInput.__channelResolutionInPs);

    t2 = (float)(resultsOfModelDecay.params['tau2'].value*userInput.__channelResolutionInPs);
    t2_err = (float)(resultsOfModelDecay.params['tau2'].stderr*userInput.__channelResolutionInPs);

    t3 = (float)(resultsOfModelDecay.params['tau3'].value*userInput.__channelResolutionInPs);
    t3_err = (float)(resultsOfModelDecay.params['tau3'].stderr*userInput.__channelResolutionInPs);

    yRes = (float)(resultsOfModelDecay.params['y0'].value);
    yRes_err = (float)(resultsOfModelDecay.params['y0'].stderr);

    amplitude1 = (float)(resultsOfModelDecay.params['ampl1'].value);
    amplitude1_err = (float)(resultsOfModelDecay.params['ampl1'].stderr);
    amplitude2 = (float)(resultsOfModelDecay.params['ampl2'].value);
    amplitude2_err = (float)(resultsOfModelDecay.params['ampl2'].stderr);
    amplitude3 = (float)(resultsOfModelDecay.params['ampl3'].value);
    amplitude3_err = (float)(resultsOfModelDecay.params['ampl3'].stderr);
    
    counts_1 = 0;
    counts_2 = 0;
    counts_3 = 0;
    
    counts_1_stddev = 0;
    counts_2_stddev = 0;
    counts_3_stddev = 0;
    
    for i in range(0, len(xVal)):
        counts_1 += amplitude1*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t1) 
        counts_2 += amplitude2*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t2) 
        counts_3 += amplitude3*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t3) 

        counts_1_stddev += np.exp(-2*xVal[i]*userInput.__channelResolutionInPs/t1)*(amplitude1_err**2 + (t1_err**2/t1**4)) + yRes_err**2
        counts_2_stddev += np.exp(-2*xVal[i]*userInput.__channelResolutionInPs/t2)*(amplitude2_err**2 + (t2_err**2/t2**4)) + yRes_err**2
        counts_3_stddev += np.exp(-2*xVal[i]*userInput.__channelResolutionInPs/t3)*(amplitude3_err**2 + (t3_err**2/t3**4)) + yRes_err**2

    counts_1_err = np.sqrt(counts_1_stddev)/np.sqrt(len(xVal));
    counts_2_err = np.sqrt(counts_2_stddev)/np.sqrt(len(xVal));
    counts_3_err = np.sqrt(counts_3_stddev)/np.sqrt(len(xVal));
    
    counts_sum = (counts_1 + counts_2 + counts_3);
    counts_sum_err = np.sqrt(counts_1_err*counts_1_err + counts_2_err*counts_2_err + counts_3_err*counts_3_err);

    I1 = (counts_1/counts_sum)
    I1_err = np.sqrt((1/counts_sum)**2*counts_1_err**2 + counts_1**2*counts_sum_err**2/counts_sum**4); 

    I2 = (counts_2/counts_sum)
    I2_err = np.sqrt((1/counts_sum)**2*counts_2_err**2 + counts_2**2*counts_sum_err**2/counts_sum**4);

    I3 = (counts_3/counts_sum)
    I3_err = np.sqrt((1/counts_sum)**2*counts_3_err**2 + counts_3**2*counts_sum_err**2/counts_sum**4);

    if userInput.__bUsingAdditionalGaussianKernel:
            kernelFWHM = (float)(resultsOfModelDecay.params['addKernelFWHM'].value);
            kernelFWHM_err = (float)(resultsOfModelDecay.params['addKernelFWHM'].stderr);
    
            str1  = "Fit results: Reconvolution 3 components (with convoluted Gaussian kernel):\n"
            str2  = "--------------------------------------------------------------------------\n"
            str3  = "X²                   = {0}\n\n".format(redChiSquare)
            str5  = "tau  (1)             = {0} ({1})\n".format(t1, t1_err)
            str6  = "I    (1)             = {0} ({1})\n".format(I1, I1_err)
            str8  = "tau  (2)             = {0} ({1})\n".format(t2, t2_err)
            str9  = "I    (2)             = {0} ({1})\n".format(I2, I2_err)
            str11 = "tau (3)             = {0} ({1})\n".format(t3, t3_err)
            str12 = "I   (3)             = {0} ({1})\n".format(I3, I3_err)
            str14 = "background          = {0} ({1})\n".format(yRes, yRes_err)
            str15 = "\n"
            str16 = "conv. kernel (FWHM) = {0} ({1})\n".format(kernelFWHM, kernelFWHM_err)
            str17 = "--------------------------------------------------------------------------"

            resStr = str1+str2+str3+str5+str6+str8+str9+str11+str12+str14+str15+str16+str17
            print(resStr)

            if userInput.__saveReconvolutionResults:
                resultsFile = open(userInput.__saveReconvolutionResultsPath,"a")
                resultsFile.write(resStr) 
                resultsFile.close()
    else:
        str1 = "Fit results: Reconvolution 3 components:\n"
        str2 = "----------------------------------------\n"
        str3 = "X²                   = {0}\n\n".format(redChiSquare)
        str5 = "tau  (1)             = {0} ({1})\n".format(t1, t1_err)
        str6 = "I    (1)             = {0} ({1})\n".format(I1, I1_err)
        str8 = "tau  (2)             = {0} ({1})\n".format(t2, t2_err)
        str9 = "I    (2)             = {0} ({1})\n".format(I2, I2_err)
        str11 = "tau (3)             = {0} ({1})\n".format(t3, t3_err)
        str12 = "I   (3)             = {0} ({1})\n".format(I3, I3_err)
        str14 = "background          = {0} ({1})\n".format(yRes, yRes_err)
        str17 = "----------------------------------------"

        resStr = str1+str2+str3+str5+str6+str8+str9+str11+str12+str14+str17
        print(resStr)

        if userInput.__saveReconvolutionResults:
            resultsFile = open(userInput.__saveReconvolutionResultsPath,"a")
            resultsFile.write(resStr) 
            resultsFile.close()

    plt.figure(3)
    ax = plt.subplot(2,1,1);
    ax.set_title("Best fit: Reconvolution with 3 lifetime components");
    plt.semilogy(xVal, ySpec,'o', xVal ,resultsOfModelDecay.best_fit, 'b');
    ax2 = plt.subplot(2,1,2);
    ax2.set_title("Best fit: Residuals");
    plt.plot(xVal, resultsOfModelDecay.residual);

if userInput.__numberOfExpDec == 4:
    print("\n\nreconvolution fitting with 4 components...\n")

    if userInput.__bUsingAdditionalGaussianKernel:
        fitModelDecay = Model(ExpDecay_4_conv);
    else:
        fitModelDecay = Model(ExpDecay_4)
        
    fitModelDecay.set_param_hint('ampl1', min=0.0);
    fitModelDecay.set_param_hint('tau1', min=0.00001);
    fitModelDecay.set_param_hint('ampl2', min=0.0);
    fitModelDecay.set_param_hint('tau2', min=0.00001);
    fitModelDecay.set_param_hint('ampl3', min=0.0);
    fitModelDecay.set_param_hint('tau3', min=0.00001);
    fitModelDecay.set_param_hint('ampl4', min=0.0);
    fitModelDecay.set_param_hint('tau4', min=0.00001);
    fitModelDecay.set_param_hint('y0', min=0.0);
    
    # (de)convolve IRF data with an additional Gaussian kernel for producing artificial broadening:
    if userInput.__bUsingAdditionalGaussianKernel:
        fitModelDecay.set_param_hint('addKernelFWHM', min=0.01);
            
        parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=(userInput.__expectedTau_1_in_ps/userInput.__channelResolutionInPs), ampl2=yMaxSpec, tau2=(userInput.__expectedTau_2_in_ps/userInput.__channelResolutionInPs), ampl3=yMaxSpec, tau3=(userInput.__expectedTau_3_in_ps/userInput.__channelResolutionInPs), ampl4=yMaxSpec, tau4=(userInput.__expectedTau_4_in_ps/userInput.__channelResolutionInPs), y0=estimatedBkgrd, x0=xWhereYMaxSpec, addKernelFWHM=userInput.__gaussianKernelFWHM, args=ySpec)
        parameterListDecayFit['addKernelFWHM'].vary = userInput.__bVaryGaussianKernelFWHM;
    else:
        parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=(userInput.__expectedTau_1_in_ps/userInput.__channelResolutionInPs), ampl2=yMaxSpec, tau2=(userInput.__expectedTau_2_in_ps/userInput.__channelResolutionInPs), ampl3=yMaxSpec, tau3=(userInput.__expectedTau_3_in_ps/userInput.__channelResolutionInPs), ampl4=yMaxSpec, tau4=(userInput.__expectedTau_4_in_ps/userInput.__channelResolutionInPs), y0=estimatedBkgrd, x0=xWhereYMaxSpec, args=ySpec)
    
    #change here if you want to fix x0 and/or y0:
    parameterListDecayFit['x0'].vary = True; 
    parameterListDecayFit['y0'].vary = not userInput.__bkgrdFixed;

    #run the fit:
    resultsOfModelDecay = fitModelDecay.fit(ySpec, params=parameterListDecayFit, weights=fitWeightingSpec, method='leastsq', x=xVal);

    #calculate results:
    chiSquare = resultsOfModelDecay.chisqr;
    redChiSquare = resultsOfModelDecay.redchi;
    
    t1 = (float)(resultsOfModelDecay.params['tau1'].value*userInput.__channelResolutionInPs);
    t1_err = (float)(resultsOfModelDecay.params['tau1'].stderr*userInput.__channelResolutionInPs);

    t2 = (float)(resultsOfModelDecay.params['tau2'].value*userInput.__channelResolutionInPs);
    t2_err = (float)(resultsOfModelDecay.params['tau2'].stderr*userInput.__channelResolutionInPs);

    t3 = (float)(resultsOfModelDecay.params['tau3'].value*userInput.__channelResolutionInPs);
    t3_err = (float)(resultsOfModelDecay.params['tau3'].stderr*userInput.__channelResolutionInPs);

    t4 = (float)(resultsOfModelDecay.params['tau4'].value*userInput.__channelResolutionInPs);
    t4_err = (float)(resultsOfModelDecay.params['tau4'].stderr*userInput.__channelResolutionInPs);

    yRes = (float)(resultsOfModelDecay.params['y0'].value);
    yRes_err = (float)(resultsOfModelDecay.params['y0'].stderr);

    amplitude1 = (float)(resultsOfModelDecay.params['ampl1'].value);
    amplitude1_err = (float)(resultsOfModelDecay.params['ampl1'].stderr);
    amplitude2 = (float)(resultsOfModelDecay.params['ampl2'].value);
    amplitude2_err = (float)(resultsOfModelDecay.params['ampl2'].stderr);
    amplitude3 = (float)(resultsOfModelDecay.params['ampl3'].value);
    amplitude3_err = (float)(resultsOfModelDecay.params['ampl3'].stderr);
    amplitude4 = (float)(resultsOfModelDecay.params['ampl4'].value);
    amplitude4_err = (float)(resultsOfModelDecay.params['ampl4'].stderr);
    
    counts_1 = 0;
    counts_2 = 0;
    counts_3 = 0;
    counts_4 = 0;
    
    counts_1_stddev = 0;
    counts_2_stddev = 0;
    counts_3_stddev = 0;
    counts_4_stddev = 0;
    
    for i in range(0, len(xVal)):
        counts_1 += amplitude1*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t1) 
        counts_2 += amplitude2*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t2) 
        counts_3 += amplitude3*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t3)
        counts_4 += amplitude4*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t4)

        counts_1_stddev += np.exp(-2*xVal[i]*userInput.__channelResolutionInPs/t1)*(amplitude1_err**2 + (t1_err**2/t1**4)) + yRes_err**2
        counts_2_stddev += np.exp(-2*xVal[i]*userInput.__channelResolutionInPs/t2)*(amplitude2_err**2 + (t2_err**2/t2**4)) + yRes_err**2
        counts_3_stddev += np.exp(-2*xVal[i]*userInput.__channelResolutionInPs/t3)*(amplitude3_err**2 + (t3_err**2/t3**4)) + yRes_err**2
        counts_4_stddev += np.exp(-2*xVal[i]*userInput.__channelResolutionInPs/t4)*(amplitude4_err**2 + (t4_err**2/t4**4)) + yRes_err**2

    counts_1_err = np.sqrt(counts_1_stddev)/np.sqrt(len(xVal));
    counts_2_err = np.sqrt(counts_2_stddev)/np.sqrt(len(xVal));
    counts_3_err = np.sqrt(counts_3_stddev)/np.sqrt(len(xVal));
    counts_4_err = np.sqrt(counts_4_stddev)/np.sqrt(len(xVal));
    
    counts_sum = (counts_1 + counts_2 + counts_3 + counts_4);
    counts_sum_err = np.sqrt(counts_1_err*counts_1_err + counts_2_err*counts_2_err + counts_3_err*counts_3_err + counts_4_err*counts_4_err);

    I1 = (counts_1/counts_sum)
    I1_err = np.sqrt((1/counts_sum)**2*counts_1_err**2 + counts_1**2*counts_sum_err**2/counts_sum**4); 

    I2 = (counts_2/counts_sum)
    I2_err = np.sqrt((1/counts_sum)**2*counts_2_err**2 + counts_2**2*counts_sum_err**2/counts_sum**4);

    I3 = (counts_3/counts_sum)
    I3_err = np.sqrt((1/counts_sum)**2*counts_3_err**2 + counts_3**2*counts_sum_err**2/counts_sum**4);

    I4 = (counts_4/counts_sum)
    I4_err = np.sqrt((1/counts_sum)**2*counts_4_err**2 + counts_4**2*counts_sum_err**2/counts_sum**4);

    if userInput.__bUsingAdditionalGaussianKernel:
            kernelFWHM = (float)(resultsOfModelDecay.params['addKernelFWHM'].value);
            kernelFWHM_err = (float)(resultsOfModelDecay.params['addKernelFWHM'].stderr);
    
            str1  = "Fit results: Reconvolution 4 components (with convoluted Gaussian kernel):\n"
            str2  = "--------------------------------------------------------------------------\n"
            str3  = "X²                   = {0}\n\n".format(redChiSquare)
            str5  = "tau  (1)             = {0} ({1})\n".format(t1, t1_err)
            str6  = "I    (1)             = {0} ({1})\n".format(I1, I1_err)
            str8  = "tau  (2)             = {0} ({1})\n".format(t2, t2_err)
            str9  = "I    (2)             = {0} ({1})\n".format(I2, I2_err)
            str11 = "tau  (3)             = {0} ({1})\n".format(t3, t3_err)
            str12 = "I    (3)             = {0} ({1})\n".format(I3, I3_err)
            str13 = "tau  (4)             = {0} ({1})\n".format(t4, t4_err)
            str18 = "I    (4)             = {0} ({1})\n".format(I4, I4_err)
            str14 = "background           = {0} ({1})\n".format(yRes, yRes_err)
            str15 = "\n"
            str16 = "conv. kernel (FWHM)  = {0} ({1})\n".format(kernelFWHM, kernelFWHM_err)
            str17 = "--------------------------------------------------------------------------"

            resStr = str1+str2+str3+str5+str6+str8+str9+str11+str12+str13+str18+str14+str15+str16+str17
            print(resStr)

            if userInput.__saveReconvolutionResults:
                resultsFile = open(userInput.__saveReconvolutionResultsPath,"a")
                resultsFile.write(resStr) 
                resultsFile.close()
    else:
        str1 = "Fit results: Reconvolution 4 components:\n"
        str2 = "----------------------------------------\n"
        str3 = "X²                   = {0}\n\n".format(redChiSquare)
        str5 = "tau  (1)             = {0} ({1})\n".format(t1, t1_err)
        str6 = "I    (1)             = {0} ({1})\n".format(I1, I1_err)
        str8 = "tau  (2)             = {0} ({1})\n".format(t2, t2_err)
        str9 = "I    (2)             = {0} ({1})\n".format(I2, I2_err)
        str11 = "tau (3)             = {0} ({1})\n".format(t3, t3_err)
        str12 = "I   (3)             = {0} ({1})\n".format(I3, I3_err)
        str13 = "tau (4)             = {0} ({1})\n".format(t4, t4_err)
        str18 = "I   (4)             = {0} ({1})\n".format(I4, I4_err)
        str14 = "background          = {0} ({1})\n".format(yRes, yRes_err)
        str17 = "----------------------------------------"

        resStr = str1+str2+str3+str5+str6+str8+str9+str11+str12+str13+str18+str14+str17
        print(resStr)

        if userInput.__saveReconvolutionResults:
            resultsFile = open(userInput.__saveReconvolutionResultsPath,"a")
            resultsFile.write(resStr) 
            resultsFile.close()

    plt.figure(3)
    ax = plt.subplot(2,1,1);
    ax.set_title("Best fit: Reconvolution with 4 lifetime components");
    plt.semilogy(xVal, ySpec,'o', xVal ,resultsOfModelDecay.best_fit, 'b');
    ax2 = plt.subplot(2,1,2);
    ax2.set_title("Best fit: Residuals");
    plt.plot(xVal, resultsOfModelDecay.residual);

#save data if required:
if userInput.__saveReconvolutionSpectrum:
    ab = np.zeros(len(xVal), dtype=[('time_[ps]', float), ('raw_counts_[#]', float), ('best_fit_[#]', float), ('raw_counts_area_normalized_[#]', float), ('best_fit_area_normalized_[#]', float)]);
    ab['time_[ps]']      = xVal*userInput.__channelResolutionInPs;
    ab['raw_counts_[#]'] = ySpec;
    ab['best_fit_[#]']   = resultsOfModelDecay.best_fit;
    ab['raw_counts_area_normalized_[#]'] = ySpec/sum(ySpec);
    ab['best_fit_area_normalized_[#]']   = resultsOfModelDecay.best_fit/sum(resultsOfModelDecay.best_fit);
    np.savetxt(userInput.__saveReconvolutionSpectrumPath, ab, fmt='%10.9f\t%10.9f\t%10.9f\t%10.9f\t%10.9f', delimiter='\t', newline='\n', header='time [ps]\traw_counts_[#]\tbest_fit_[#]\nraw_counts_area_normalized_[#]\nbest_fit_area_normalized_[#]\n');

    abRes = np.zeros(len(xVal), dtype=[('time_[ps]', float), ('residuals_[stddev]', float)]);
    abRes['time_[ps]']               = xVal*userInput.__channelResolutionInPs;
    abRes['residuals_[stddev]'] = resultsOfModelDecay.residual;
    np.savetxt(userInput.__saveReconvolutionSpectrumResidualPath, abRes, fmt='%10.9f\t%10.9f', delimiter='\t', newline='\n', header='time [ps]\traw_counts_[#]\tresiduals_[stddev]\n');

if userInput.__saveIRFSpectrum:
    ab = np.zeros(len(xVal), dtype=[('time_[ps]', float), ('raw_counts_[#]', float), ('raw_counts_area_normalized[#]', float)]);
    ab['time_[ps]']      = xVal*userInput.__channelResolutionInPs;
    ab['raw_counts_[#]'] = yIRF;
    ab['raw_counts_area_normalized[#]'] = yIRF/sum(yIRF);
    np.savetxt(userInput.__saveIRFSpectrumPath, ab, fmt='%10.9f\t%10.9f\t%10.9f', delimiter='\t', newline='\n', header='time_[ps]\traw_counts_[#]\nraw_counts_area_normalized[#]\n');

if userInput.__saveReconvolutionIRF and userInput.__bUsingModel:
    ab = np.zeros(len(xVal), dtype=[('time [ps]', float), ('raw counts [#]', float), ('best fit [#]', float)]);
    ab['time [ps]']      = xVal*userInput.__channelResolutionInPs;
    ab['raw counts [#]'] = yIRFOrigin;
    ab['best fit [#]']   = yIRF;
    np.savetxt(userInput.__saveReconvolutionIRFPath, ab, fmt='%10.9f\t%10.9f\t%10.9f', delimiter='\t', newline='\n', header='time [ps]\traw counts [#]\tbest fit [#]\n');

    abRes = np.zeros(len(xVal), dtype=[('time [ps]', float), ('residuals [conv. level]', float)]);
    abRes['time [ps]']               = xVal*userInput.__channelResolutionInPs;
    abRes['residuals [conv. level]'] = resultsOfModelIRF.residual;
    np.savetxt(userInput.__saveReconvolutionIRFResidualPath, abRes, fmt='%10.9f\t%10.9f', delimiter='\t', newline='\n', header='time [ps]\traw counts [#]\tresiduals [conv. level]\n');

plt.show();
    
