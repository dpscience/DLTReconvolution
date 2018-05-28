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

import DReconvolutionModel as functionModelList
import DReconvolutionInput as userInput

from lmfit import Model
import numpy as np
import matplotlib.pyplot as plt

print("....started....");

xSpec,ySpec = np.loadtxt(userInput.__filePathSpec, delimiter=userInput.__specDataDelimiter, unpack=True, dtype='float');
xIRF,yIRF   = np.loadtxt(userInput.__filePathIRF, delimiter=userInput.__irfDataDelimiter, unpack=True, dtype='float');

print("shifting x to x = 0...")

#shift that: [0] = 0:
xSpec -= xSpec[0];
xIRF  -= xIRF[0];

if (len(xSpec) != len(xIRF) or len(ySpec) != len(yIRF)):
    print("Error: the data vectors (x,y) of IRF and SPECTRUM require the same length.");
    quit(); #kill the process on error.

#merge into one vector (xVal) for x values:
xVal = np.zeros(xSpec.size);
xVal = xSpec;

#calculate the vector which contains the fit weightings:
fitWeightingIRF  = 1.0;
fitWeightingSpec = 1.0;

print("calculating fit weights...")

if userInput.__bUsingYVarAsWeighting:
    fitWeightingIRF  = 1.0/np.sqrt(yIRF+1);  #prevent zero devision
    fitWeightingSpec = 1.0/np.sqrt(ySpec+1); #prevent zero devision

print("estimate start values...")

#estimate IRF start values (amplitude, xOffset/mean, stddev):
yMaxIRF  = np.amax(yIRF);
yMaxSpec = np.amax(ySpec);
 
xWhereYMaxIRF  = np.argmax(yMaxIRF);
xWhereYMaxSpec = np.argmax(yMaxSpec);

stddevIRF = 1.0

for i in range(0, len(xVal)-1):
    if yIRF[i] > yMaxIRF*0.5:
        stddevIRF = np.abs((xWhereYMaxIRF-xIRF[i]))/(2*np.sqrt(2*np.log(2)));
        break;

estimatedBkgrd = 0;
for i in range(userInput.__bkgrd_startIndex, userInput.__bkgrd_startIndex + userInput.__bkgrd_count):
    estimatedBkgrd += ySpec[i];

estimatedBkgrd /= userInput.__bkgrd_count;
    
#fit the IRF model function on data (xIRF, yIRF):
if userInput.__bUsingModel:
    print("IRF: running model fit...")
    #Gaussian:
    if userInput.__modelType == functionModelList.ReconvolutionModel.Gaussian:
        fitModelIRF = Model(functionModelList.Gaussian);
        fitModelIRF.set_param_hint('sigma', min=0.0);
        fitModelIRF.set_param_hint('ampl', min=0.0);
        fitModelIRF.set_param_hint('y0', min=0.0);

        parameterListIRFFit = fitModelIRF.make_params(x=xVal, ampl=yMaxIRF, sigma=stddevIRF, y0=0, x0=xWhereYMaxIRF, args=yIRF);
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
        
        parameterListIRFFit = fitModelIRF.make_params(x=xVal, ampl=yMaxIRF, s=stddevIRF, y0=0, x0=xWhereYMaxIRF, args=yIRF);
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


        print("\nFit results: IRF - Model Function (Type: Lorentz/Cauchy):");
        print("-----------------------------------------------------------");
        print("X²          = {0}".format(redChiSquare));
        print("");
        print("s [ps]      = {0} ({1})".format(s, s_err));
        print("");
        print("amplitude   = {0} ({1})".format(amplitude, amplitude_err));
        print("");
        print("background  = {0} ({1})".format(yRes, yRes_err));
        print("-----------------------------------------------------------");

        
        plt.figure(1);
        ax = plt.subplot(2,1,1);
        ax.set_title("Best fit: IRF Lorentz/Cauchy model fit");
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

        parameterListIRFFit = fitModelIRF.make_params(x=xVal, ampl=yMaxIRF, a=1, sigma=stddevIRF, s=stddevIRF, y0=0, x0=xWhereYMaxIRF, args=yIRF);
        #change here if you want to fix x0 and/or y0:
        parameterListIRFFit['x0'].vary = True; 
        parameterListIRFFit['y0'].vary = True;

        #run the fit:
        resultsOfModelIRF = fitModelIRF.fit(yIRF, params=parameterListIRFFit, weights=fitWeightingIRF, method='leastsq', x=xVal);

        chiSquare = resultsOfModelIRF.chisqr;
        redChiSquare = resultsOfModelIRF.redchi;

        fract = (float)(resultsOfModelIRF.params['a'].value*userInput);
        fract_err = (float)(resultsOfModelIRF.params['a'].stderr*userInput);
        
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
        print("G - Gaussian:  a        = {0} ({})".format(fract, fract_err));
        print("---------------------------------------------------------------");
        print("stddev [ps] = {0} ({1})".format(sigma, sigma_err));
        print("FWHM   [ps] = {0} ({1})".format(sigma*(2*np.sqrt(2*np.log(2))), sigma_err*(2*np.sqrt(2*np.log(2)))));
        print("");
        print("---------------------------------------------------------------");
        print("L - Lorentzian: (1 - a) = {0} ({})".format(1-fract, fract_err));
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

        parameterListIRFFit = fitModelIRF.make_params(x=xVal, ampl=yMaxIRF, alpha=stddevIRF, m=2, y0=0, x0=xWhereYMaxIRF, args=yIRF);
        #change here if you want to fix x0 and/or y0:
        parameterListIRFFit['x0'].vary = True; 
        parameterListIRFFit['y0'].vary = True;

        #run the fit:
        resultsOfModelIRF = fitModelIRF.fit(yIRF, params=parameterListIRFFit, weights=fitWeightingIRF, method='leastsq', x=xVal);

        chiSquare = resultsOfModelIRF.chisqr;
        redChiSquare = resultsOfModelIRF.redchi;

        alpha = (float)(resultsOfModelIRF.params['alpha'].value*userInput.__channelResolutionInPs);
        alpha_err = (float)(resultsOfModelIRF.params['alpha'].stderr*userInput.__channelResolutionInPs);
        
        m = (float)(resultsOfModelIRF.params['m'].value*userInput);
        m_err = (float)(resultsOfModelIRF.params['m'].stderr*userInput);
        
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


def convolveData(a, b):
    A = np.fft.fft(a);
    B = np.fft.fft(b);
    convAB = np.real(np.fft.ifft(A*B));
    return convAB;

#1 component expontential distribution:
def ExpDecay_1(x, ampl1, tau1, y0, x0, args=(yIRF)):
    h = np.zeros(x.size) 
    lengthVec = len(ySpec)

    shift_1 = np.remainder(np.remainder(x-np.floor(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr1 = (1 - x0 + np.floor(x0))*yIRF[shift_1.astype(int)]
    
    shift_2 = np.remainder(np.remainder(x-np.ceil(c)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr2 = (x0 - np.floor(x0))*yIRF[shift_2.astype(int)]
    
    irf_shifted = (shift_Incr1 + shift_Incr2)
    irf_norm = irf_shifted/sum(irf_shifted)
    
    h = ampl1*np.exp(-(x)/tau1)
    hConvIrf_norm = convolveData(h, irf_norm)
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

#3 component expontential distribution:
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
    print("\nrunning reconvolution with 1 component...\n")
    fitModelDecay = Model(ExpDecay_1);
    fitModelDecay.set_param_hint('ampl1', min=0.0);
    fitModelDecay.set_param_hint('tau1', min=0.00001);
    fitModelDecay.set_param_hint('y0', min=0.0);

    parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=(userInput.__expectedTau_1_in_ps/userInput.__channelResolutionInPs), y0=estimatedBkgrd, x0=xWhereYMaxSpec, args=ySpec);
    #change here if you want to fix x0 and/or y0:
    parameterListDecayFit['x0'].vary = True; 
    parameterListDecayFit['y0'].vary = True;

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
        counts_1 += amplitude1*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t1) + yRes;
        
        counts_1_stddev += np.exp(-2*xVal[i]*userInput.__channelResolutionInPs/t1)*(amplitude1_err**2 + (t1_err**2/t1**4)) + yRes_err**2
        
    counts_1_err = np.sqrt(counts_1_stddev)/np.sqrt(len(xVal));
    
    counts_sum = (counts_1);
    counts_sum_err = np.sqrt(counts_1_err*counts_1_err);

    I1 = (counts_1/counts_sum)
    I1_err = np.sqrt((1/counts_sum)**2*counts_1_err**2 + counts_1**2*counts_sum_err**2/counts_sum**4); 


    print("Fit results: Reconvolution 1 component:");
    print("---------------------------------------");
    print("X²         = {0}".format(redChiSquare));
    print("");
    print("tau (1)    = {0} ({1})".format(t1, t1_err));
    print("I   (1)    = {0} ({1})".format(I1, I1_err));
    print("");
    print("background = {0} ({1})".format(yRes, yRes_err));
    print("---------------------------------------");

    print("\nplotting results...")
    plt.figure(3)
    ax = plt.subplot(2,1,1);
    ax.set_title("Best fit: Reconvolution with 1 lifetime component");
    plt.semilogy(xVal, ySpec,'o', xVal ,resultsOfModelDecay.best_fit, 'b');
    ax2 = plt.subplot(2,1,2);
    ax2.set_title("Best fit: Residuals");
    plt.plot(xVal, resultsOfModelDecay.residual);
    plt.show();

if userInput.__numberOfExpDec == 2:
    print("\nrunning reconvolution with 2 component...\n")
    fitModelDecay = Model(ExpDecay_2);
    fitModelDecay.set_param_hint('ampl1', min=0.0);
    fitModelDecay.set_param_hint('tau1', min=0.00001);
    fitModelDecay.set_param_hint('ampl2', min=0.0);
    fitModelDecay.set_param_hint('tau2', min=0.00001);
    fitModelDecay.set_param_hint('y0', min=0.0);

    parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=(userInput.__expectedTau_1_in_ps/userInput.__channelResolutionInPs), ampl2=yMaxSpec, tau2=(userInput.__expectedTau_2_in_ps/userInput.__channelResolutionInPs), y0=estimatedBkgrd, x0=xWhereYMaxSpec, args=ySpec);
    #change here if you want to fix x0 and/or y0:
    parameterListDecayFit['x0'].vary = True; 
    parameterListDecayFit['y0'].vary = True;
    
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
        counts_1 += amplitude1*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t1) + yRes;
        counts_2 += amplitude2*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t2) + yRes;

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

    
    print("Fit results: Reconvolution 2 components:");
    print("----------------------------------------");
    print("X²         = {0}".format(redChiSquare));
    print("");
    print("tau (1)    = {0} ({1})".format(t1, t1_err));
    print("I   (1)    = {0} ({1})".format(I1, I1_err));
    print("");
    print("tau (2)    = {0} ({1})".format(t2, t2_err));
    print("I   (2)    = {0} ({1})".format(I2, I2_err));
    print("");
    print("background = {0} ({1})".format(yRes, yRes_err));
    print("----------------------------------------");
    

    print("\nplotting results...")
    plt.figure(3)
    ax = plt.subplot(2,1,1);
    ax.set_title("Best fit: Reconvolution with 2 lifetime components");
    plt.semilogy(xVal, ySpec,'o', xVal ,resultsOfModelDecay.best_fit, 'b');
    ax2 = plt.subplot(2,1,2);
    ax2.set_title("Best fit: Residuals");
    plt.plot(xVal, resultsOfModelDecay.residual);
    plt.show();

if userInput.__numberOfExpDec == 3:
    print("\nrunning reconvolution with 3 component...\n")
    fitModelDecay = Model(ExpDecay_3);
    fitModelDecay.set_param_hint('ampl1', min=0.0);
    fitModelDecay.set_param_hint('tau1', min=0.00001);
    fitModelDecay.set_param_hint('ampl2', min=0.0);
    fitModelDecay.set_param_hint('tau2', min=0.00001);
    fitModelDecay.set_param_hint('ampl3', min=0.0);
    fitModelDecay.set_param_hint('tau3', min=0.00001);
    fitModelDecay.set_param_hint('y0', min=0.0);

    parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=(userInput.__expectedTau_1_in_ps/userInput.__channelResolutionInPs), ampl2=yMaxSpec, tau2=(userInput.__expectedTau_2_in_ps/userInput.__channelResolutionInPs), ampl3=yMaxSpec, tau3=(userInput.__expectedTau_3_in_ps/userInput.__channelResolutionInPs), y0=estimatedBkgrd, x0=xWhereYMaxSpec, args=ySpec);
    #change here if you want to fix x0 and/or y0:
    parameterListDecayFit['x0'].vary = True; 
    parameterListDecayFit['y0'].vary = True;

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
        counts_1 += amplitude1*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t1) + yRes;
        counts_2 += amplitude2*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t2) + yRes;
        counts_3 += amplitude3*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t3) + yRes;

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
    
    
    print("Fit results: Reconvolution 2 components:");
    print("----------------------------------------");
    print("X²         = {0}".format(redChiSquare));
    print("");
    print("tau (1)    = {0} ({1})".format(t1, t1_err));
    print("I   (1)    = {0} ({1})".format(I1, I1_err));
    print("");
    print("tau (2)    = {0} ({1})".format(t2, t2_err));
    print("I   (2)    = {0} ({1})".format(I2, I2_err));
    print("");
    print("tau (3)    = {0} ({1})".format(t3, t3_err));
    print("I   (3)    = {0} ({1})".format(I3, I3_err));
    print("");
    print("background = {0} ({1})".format(yRes, yRes_err));
    print("----------------------------------------");

    print("\nplotting results...")
    plt.figure(3)
    ax = plt.subplot(2,1,1);
    ax.set_title("Best fit: Reconvolution with 3 lifetime components");
    plt.semilogy(xVal, ySpec,'o', xVal ,resultsOfModelDecay.best_fit, 'b');
    ax2 = plt.subplot(2,1,2);
    ax2.set_title("Best fit: Residuals");
    plt.plot(xVal, resultsOfModelDecay.residual);
    plt.show();

if userInput.__numberOfExpDec == 4:
    print("\nrunning reconvolution with 4 component...\n")
    fitModelDecay = Model(ExpDecay_4);
    fitModelDecay.set_param_hint('ampl1', min=0.0);
    fitModelDecay.set_param_hint('tau1', min=0.00001);
    fitModelDecay.set_param_hint('ampl2', min=0.0);
    fitModelDecay.set_param_hint('tau2', min=0.00001);
    fitModelDecay.set_param_hint('ampl3', min=0.0);
    fitModelDecay.set_param_hint('tau3', min=0.00001);
    fitModelDecay.set_param_hint('ampl4', min=0.0);
    fitModelDecay.set_param_hint('tau4', min=0.00001);
    fitModelDecay.set_param_hint('y0', min=0.0);

    parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=(userInput.__expectedTau_1_in_ps/userInput.__channelResolutionInPs), ampl2=yMaxSpec, tau2=(userInput.__expectedTau_2_in_ps/userInput.__channelResolutionInPs), ampl3=yMaxSpec, tau3=(userInput.__expectedTau_3_in_ps/userInput.__channelResolutionInPs), ampl4=yMaxSpec, tau4=(userInput.__expectedTau_4_in_ps/userInput.__channelResolutionInPs), y0=estimatedBkgrd, x0=xWhereYMaxSpec, args=ySpec);
    #change here if you want to fix x0 and/or y0:
    parameterListDecayFit['x0'].vary = True; 
    parameterListDecayFit['y0'].vary = True;

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
        counts_1 += amplitude1*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t1) + yRes;
        counts_2 += amplitude2*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t2) + yRes;
        counts_3 += amplitude3*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t3) + yRes;
        counts_4 += amplitude4*np.exp(-(xVal[i]*userInput.__channelResolutionInPs)/t4) + yRes;

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
    
    
    print("Fit results: Reconvolution 2 components:");
    print("----------------------------------------");
    print("X²         = {0}".format(redChiSquare));
    print("");
    print("tau (1)    = {0} ({1})".format(t1, t1_err));
    print("I   (1)    = {0} ({1})".format(I1, I1_err));
    print("");
    print("tau (2)    = {0} ({1})".format(t2, t2_err));
    print("I   (2)    = {0} ({1})".format(I2, I2_err));
    print("");
    print("tau (3)    = {0} ({1})".format(t3, t3_err));
    print("I   (3)    = {0} ({1})".format(I3, I3_err));
    print("");
    print("tau (4)    = {0} ({1})".format(t4, t4_err));
    print("I   (4)    = {0} ({1})".format(I4, I4_err));
    print("");
    print("background = {0} ({1})".format(yRes, yRes_err));
    print("----------------------------------------");

    print("\nplotting results...")
    plt.figure(3)
    ax = plt.subplot(2,1,1);
    ax.set_title("Best fit: Reconvolution with 2 lifetime components");
    plt.semilogy(xVal, ySpec,'o', xVal ,resultsOfModelDecay.best_fit, 'b');
    ax2 = plt.subplot(2,1,2);
    ax2.set_title("Best fit: Residuals");
    plt.plot(xVal, resultsOfModelDecay.residual);
    plt.show();
    
