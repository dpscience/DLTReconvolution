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

#shift that: [0] = 0:
xSpec -= xSpec[0];
xIRF  -= xIRF[0];

print("shifting x to x = 0...")

#convert x bins to time scale [ps]:
xIRF  *= userInput.__channelResInPicoseconds;
xSpec *= userInput.__channelResInPicoseconds;

print("converting bins into time scale...")

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
 
xWhereYMaxIRF  = np.argmax(yMaxIRF)*userInput.__channelResInPicoseconds;
xWhereYMaxSpec = np.argmax(yMaxSpec)*userInput.__channelResInPicoseconds;

stddevIRF = 1.0

for i in range(0, len(xVal)-1):
    if yIRF[i] > yMaxIRF*0.5:
        stddevIRF = np.abs((xWhereYMaxIRF-xIRF[i]))/(2*np.sqrt(2*np.log(2)));
        break;
    
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

        plt.figure(1);
        plt.subplot(2,1,1);
        plt.semilogy(xVal, yIRF,'o', xVal , resultsOfModelIRF.best_fit, 'b');
        plt.subplot(2,1,2);
        plt.plot(xVal, resultsOfModelIRF.residual);
        
        plt.figure(2)
        plt.semilogy(xVal, ySpec, xVal, yIRF, xVal, resultsOfModelIRF.best_fit)
        
        #replace input by model fit data:
        yIRF = resultsOfModelIRF.best_fit; 

        print("\nFit results: IRF - Model Function (Type: Gaussian):");
        print("---------------------------------------------------");
        print(resultsOfModelIRF.fit_report());
        print("---------------------------------------------------");

    #Loretz/Cauchy:
    if userInput.__modelType == functionModelList.ReconvolutionModel.Lorentz_Cauchy:
        fitModelIRF = Model(functionModelList.Lorentz_Cauchy);
        fitModelIRF.set_param_hint('a', min=0.0);
        fitModelIRF.set_param_hint('ampl', min=0.0);
        fitModelIRF.set_param_hint('y0', min=0.0);
        fitModelIRF.set_param_hint('wing', min=0.0);

        parameterListIRFFit = fitModelIRF.make_params(x=xVal, ampl=yMaxIRF, a=1, wing=stddevIRF, y0=0, x0=xWhereYMaxIRF, args=yIRF);
        #change here if you want to fix x0 and/or y0:
        parameterListIRFFit['x0'].vary = True; 
        parameterListIRFFit['y0'].vary = True;

        #run the fit:
        resultsOfModelIRF = fitModelIRF.fit(yIRF, params=parameterListIRFFit, weights=fitWeightingIRF, method='leastsq', x=xVal);

        plt.figure(1);
        plt.subplot(2,1,1);
        plt.semilogy(xVal, yIRF,'o', xVal , resultsOfModelIRF.best_fit, 'b');
        plt.subplot(2,1,2);
        plt.plot(xVal, resultsOfModelIRF.residual);
        
        plt.figure(2)
        plt.semilogy(xVal, ySpec, xVal, yIRF, xVal, resultsOfModelIRF.best_fit)
        
        #replace input by model fit data:
        yIRF = resultsOfModelIRF.best_fit; 

        print("\nFit results: IRF - Model Function (Type: Lorentz/Cauchy):");
        print("---------------------------------------------------------");
        print(resultsOfModelIRF.fit_report());
        print("---------------------------------------------------------");

    #PseudoVoigt-1:
    if userInput.__modelType == functionModelList.ReconvolutionModel.Pseudovoigt1:
        fitModelIRF = Model(functionModelList.Pseudovoigt1);
        fitModelIRF.set_param_hint('a', min=0.0);
        fitModelIRF.set_param_hint('sigma', min=0.0);
        fitModelIRF.set_param_hint('wing', min=0.0);
        fitModelIRF.set_param_hint('ampl', min=0.0);
        fitModelIRF.set_param_hint('y0', min=0.0);

        parameterListIRFFit = fitModelIRF.make_params(x=xVal, ampl=yMaxIRF, a=1, sigma=stddevIRF, wing=stddevIRF, y0=0, x0=xWhereYMaxIRF, args=yIRF);
        #change here if you want to fix x0 and/or y0:
        parameterListIRFFit['x0'].vary = True; 
        parameterListIRFFit['y0'].vary = True;

        #run the fit:
        resultsOfModelIRF = fitModelIRF.fit(yIRF, params=parameterListIRFFit, weights=fitWeightingIRF, method='leastsq', x=xVal);

        plt.figure(1);
        plt.subplot(2,1,1);
        plt.semilogy(xVal, yIRF,'o', xVal , resultsOfModelIRF.best_fit, 'b');
        plt.subplot(2,1,2);
        plt.plot(xVal, resultsOfModelIRF.residual);
        
        plt.figure(2)
        plt.semilogy(xVal, ySpec, xVal, yIRF, xVal, resultsOfModelIRF.best_fit)
        
        #replace input by model fit data:
        yIRF = resultsOfModelIRF.best_fit; 

        print("\nFit results: IRF - Model Function (Type: PseudoVoigt type 1):");
        print("-------------------------------------------------------------");
        print(resultsOfModelIRF.fit_report());
        print("-------------------------------------------------------------");

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

        plt.figure(1);
        plt.subplot(2,1,1);
        plt.semilogy(xVal, yIRF,'o', xVal , resultsOfModelIRF.best_fit, 'b');
        plt.subplot(2,1,2);
        plt.plot(xVal, resultsOfModelIRF.residual);
        
        plt.figure(2)
        plt.semilogy(xVal, ySpec, xVal, yIRF, xVal, resultsOfModelIRF.best_fit)
        
        #replace input by model fit data:
        yIRF = resultsOfModelIRF.best_fit; 

        print("\nFit results: IRF - Model Function (Type: Pearson type 7):");
        print("---------------------------------------------------------");
        print(resultsOfModelIRF.fit_report());
        print("---------------------------------------------------------");


#applying reconvolution:
if userInput.__numberOfExpDec == 1:
    print("\nrunning reconvolution with 1 component...")
    fitModelDecay = Model(functionModelList.ExpDecay_1);
    fitModelDecay.set_param_hint('ampl1', min=0.0);
    fitModelDecay.set_param_hint('tau1', min=0.00001);

    parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=200.0, y0=0, x0=xWhereYMaxSpec, args=ySpec);
    #change here if you want to fix x0 and/or y0:
    parameterListDecayFit['x0'].vary = True; 
    parameterListDecayFit['y0'].vary = True;

    #run the fit:
    resultsOfModelDecay = fitModelDecay.fit(ySpec, params=parameterListDecayFit, weights=fitWeightingSpec, method='leastsq', x=xVal);

    print("Fit results: Reconvolution 1 component:");
    print("---------------------------------------");
    print(resultsOfModelDecay.fit_report());
    print("---------------------------------------");

    print("\plotting results...")
    plt.figure(3)
    plt.subplot(2,1,1);
    plt.semilogy(xVal, ySpec,'o', xVal ,resultsOfModelDecay.best_fit, 'b');
    plt.subplot(2,1,2);
    plt.plot(xVal, resultsOfModelDecay.residual);
    plt.show();

if userInput.__numberOfExpDec == 2:
    print("\nrunning reconvolution with 2 component...")
    fitModelDecay = Model(functionModelList.ExpDecay_2);
    fitModelDecay.set_param_hint('ampl1', min=0.0);
    fitModelDecay.set_param_hint('tau1', min=0.00001);
    fitModelDecay.set_param_hint('ampl2', min=0.0);
    fitModelDecay.set_param_hint('tau2', min=0.00001);

    parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=200.0, ampl2=yMaxSpec, tau2=200.0, y0=0, x0=xWhereYMaxSpec, args=ySpec);
    #change here if you want to fix x0 and/or y0:
    parameterListDecayFit['x0'].vary = True; 
    parameterListDecayFit['y0'].vary = True;

    #run the fit:
    resultsOfModelDecay = fitModelDecay.fit(ySpec, params=parameterListDecayFit, weights=fitWeightingSpec, method='leastsq', x=xVal);

    print("Fit results: Reconvolution 2 component:");
    print("---------------------------------------");
    print(resultsOfModelDecay.fit_report());
    print("---------------------------------------");

    print("\nplotting results...")
    plt.figure(3)
    plt.subplot(2,1,1);
    plt.semilogy(xVal, ySpec,'o', xVal ,resultsOfModelDecay.best_fit, 'b');
    plt.subplot(2,1,2);
    plt.plot(xVal, resultsOfModelDecay.residual);
    plt.show();

if userInput.__numberOfExpDec == 3:
    print("\nrunning reconvolution with 3 component...")
    fitModelDecay = Model(functionModelList.ExpDecay_3);
    fitModelDecay.set_param_hint('ampl1', min=0.0);
    fitModelDecay.set_param_hint('tau1', min=0.00001);
    fitModelDecay.set_param_hint('ampl2', min=0.0);
    fitModelDecay.set_param_hint('tau2', min=0.00001);
    fitModelDecay.set_param_hint('ampl3', min=0.0);
    fitModelDecay.set_param_hint('tau3', min=0.00001);

    parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=200.0, ampl2=yMaxSpec, tau2=200.0, ampl3=yMaxSpec, tau3=200.0, y0=0, x0=xWhereYMaxSpec, args=ySpec);
    #change here if you want to fix x0 and/or y0:
    parameterListDecayFit['x0'].vary = True; 
    parameterListDecayFit['y0'].vary = True;

    #run the fit:
    resultsOfModelDecay = fitModelDecay.fit(ySpec, params=parameterListDecayFit, weights=fitWeightingSpec, method='leastsq', x=xVal);

    print("Fit results: Reconvolution 3 component:");
    print("---------------------------------------");
    print(resultsOfModelDecay.fit_report());
    print("---------------------------------------");

    print("\nplotting results...")
    plt.figure(3)
    plt.subplot(2,1,1);
    plt.semilogy(xVal, ySpec,'o', xVal ,resultsOfModelDecay.best_fit, 'b');
    plt.subplot(2,1,2);
    plt.plot(xVal, resultsOfModelDecay.residual);
    plt.show();

if userInput.__numberOfExpDec == 4:
    print("\nrunning reconvolution with 4 component...")
    fitModelDecay = Model(functionModelList.ExpDecay_4);
    fitModelDecay.set_param_hint('ampl1', min=0.0);
    fitModelDecay.set_param_hint('tau1', min=0.00001);
    fitModelDecay.set_param_hint('ampl2', min=0.0);
    fitModelDecay.set_param_hint('tau2', min=0.00001);
    fitModelDecay.set_param_hint('ampl3', min=0.0);
    fitModelDecay.set_param_hint('tau3', min=0.00001);
    fitModelDecay.set_param_hint('ampl4', min=0.0);
    fitModelDecay.set_param_hint('tau4', min=0.00001);

    parameterListDecayFit = fitModelDecay.make_params(x=xVal, ampl1=yMaxSpec, tau1=200.0, ampl2=yMaxSpec, tau2=200.0, ampl3=yMaxSpec, tau3=200.0, ampl4=yMaxSpec, tau4=200.0, y0=0, x0=xWhereYMaxSpec, args=ySpec);
    #change here if you want to fix x0 and/or y0:
    parameterListDecayFit['x0'].vary = True; 
    parameterListDecayFit['y0'].vary = True;

    #run the fit:
    resultsOfModelDecay = fitModelDecay.fit(ySpec, params=parameterListDecayFit, weights=fitWeightingSpec, method='leastsq', x=xVal);

    print("Fit results: Reconvolution 4 component:");
    print("---------------------------------------");
    print(resultsOfModelDecay.fit_report());
    print("---------------------------------------");

    print("\nplotting results...")
    plt.figure(3)
    plt.subplot(2,1,1);
    plt.semilogy(xVal, ySpec,'o', xVal ,resultsOfModelDecay.best_fit, 'b');
    plt.subplot(2,1,2);
    plt.plot(xVal, resultsOfModelDecay.residual);
    plt.show();
    
