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

import numpy as np
from scipy.special import beta
from enum import Enum

#list of available models:
class ReconvolutionModel(Enum):
     Gaussian       = 1
     Lorentz_Cauchy = 2
     Pseudovoigt1   = 3
     Pearson7       = 4

#definition of available models:
def Gaussian(x, ampl, sigma, y0, x0, args=()):
    h=np.zeros(x.size)
    N=np.zeros(x.size)
    N=1.0/(sigma*np.sqrt(2*np.pi))
    h=N*np.exp(-0.5*((x-x0)/sigma)**2);
    return ampl*h+y0

def Lorentz_Cauchy(x, ampl, a, wing, y0, x0, args=()):
    h=np.zeros(x.size)
    h=wing/(np.pi*((x-x0)*(x-x0) + wing*wing))
    return ampl*h+y0

def Pseudovoigt1(x, ampl, a, sigma, wing, y0, x0, args=()):
    G=np.zeros(x.size)
    L=np.zeros(x.size)
    G=(1.0/(sigma*np.sqrt(2*np.pi)))*np.exp(-0.5*((x-x0)/sigma)*((x-x0)/sigma))
    L=wing/(np.pi*((x-x0)*(x-x0) + wing*wing))
    return ampl*(a*G+(1-a)*L)+y0

def Pearson7(x, ampl, alpha, m, y0, x0, args=()):
    h=np.zeros(x.size)
    h=(1/(alpha*beta(m-0.5,0.5)))*(1+((x-x0)/alpha)**2)**(-m)
    return ampl*h+y0


def convolveData(a, b): 
    return np.real(np.fft.ifft(np.fft.fft(a)*np.fft.fft(b)))

#1 component expontential distribution:
def ExpDecay_1(x, ampl1, tau1, y0, x0, args=()):  
    h = np.zeros(x.size) 
    lengthVec = len(args)
    
    shift_1 = np.remainder(np.remainder(x-np.floor(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr1 = (1 - x0 + np.floor(x0))*args[shift_1.astype(int)]
    
    shift_2 = np.remainder(np.remainder(x-np.ceil(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr2 = (x0 - np.floor(x0))*args[shift_2.astype(int)]
    
    irf_shifted = (shift_Incr1 + shift_Incr2)
    irf_norm = irf_shifted/sum(irf_shifted)
    
    h = ampl1*np.exp(-(x)/tau1)
    hConvIrf_norm = convolveData(h, irf_norm) + y0
    return hConvIrf_norm

#2 component expontential distribution:
def ExpDecay_2(x, ampl1, tau1, ampl2, tau2, y0, x0, args=()):  
    h = np.zeros(x.size) 
    lengthVec = len(x)

    shift_1 = np.remainder(np.remainder(x-np.floor(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr1 = (1 - x0 + np.floor(x0))*args[shift_1.astype(int)]
    
    shift_2 = np.remainder(np.remainder(x-np.ceil(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr2 = (x0 - np.floor(x0))*args[shift_2.astype(int)]
    
    irf_shifted = (shift_Incr1 + shift_Incr2)
    irf_norm = irf_shifted/sum(irf_shifted)
    
    h = ampl1*np.exp(-(x)/tau1) + ampl2*np.exp(-(x)/tau2)
    hConvIrf_norm = convolveData(h, irf_norm) + y0
    return hConvIrf_norm

#3 component expontential distribution:
def ExpDecay_3(x, ampl1, tau1, ampl2, tau2, ampl3, tau3, y0, x0, args=()):  
    h = np.zeros(x.size) 
    lengthVec = len(args)
    
    shift_1 = np.remainder(np.remainder(x-np.floor(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr1 = (1 - x0 + np.floor(x0))*args[shift_1.astype(int)]
    
    shift_2 = np.remainder(np.remainder(x-np.ceil(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr2 = (x0 - np.floor(x0))*args[shift_2.astype(int)]
    
    irf_shifted = (shift_Incr1 + shift_Incr2)
    irf_norm = irf_shifted/sum(irf_shifted)
    
    h = ampl1*np.exp(-(x)/tau1) + ampl2*np.exp(-(x)/tau2) + ampl3*np.exp(-(x)/tau3)
    hConvIrf_norm = convolveData(h, irf_norm) + y0
    return hConvIrf_norm

#4 component expontential distribution:
def ExpDecay_4(x, ampl1, tau1, ampl2, tau2, ampl3, tau3, ampl4, tau4, y0, x0, args=()):  
    h = np.zeros(x.size) 
    lengthVec = len(args)
    
    shift_1 = np.remainder(np.remainder(x-np.floor(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr1 = (1 - x0 + np.floor(x0))*args[shift_1.astype(int)]
    
    shift_2 = np.remainder(np.remainder(x-np.ceil(x0)-1, lengthVec) + lengthVec, lengthVec)
    shift_Incr2 = (x0 - np.floor(x0))*args[shift_2.astype(int)]
    
    irf_shifted = (shift_Incr1 + shift_Incr2)
    irf_norm = irf_shifted/sum(irf_shifted)
    
    h = ampl1*np.exp(-(x)/tau1) + ampl2*np.exp(-(x)/tau2) + ampl3*np.exp(-(x)/tau3) + ampl4*np.exp(-(x)/tau4)
    hConvIrf_norm = convolveData(h, irf_norm) + y0
    return hConvIrf_norm





