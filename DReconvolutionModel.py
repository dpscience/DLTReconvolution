#*************************************************************************************************
#**
#** DLTReconvolution v1.2 (17.01.2019)
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

def Lorentz_Cauchy(x, ampl, s, y0, x0, args=()):
    h=np.zeros(x.size)
    h=s/(np.pi*((x-x0)*(x-x0) + s*s))
    return ampl*h+y0

def Pseudovoigt1(x, ampl, a, sigma, s, y0, x0, args=()):
    G=np.zeros(x.size)
    L=np.zeros(x.size)

    G=(1.0/(sigma*np.sqrt(2*np.pi)))*np.exp(-0.5*((x-x0)/sigma)**2)
    L=s/(np.pi*((x-x0)**2 + s*s))
    return ampl*(a*G+(1-a)*L)+y0

def Pearson7(x, ampl, alpha, m, y0, x0, args=()):
    h=np.zeros(x.size)
    h=(1/(alpha*beta(m-0.5,0.5)))*(1+((x-x0)/alpha)**2)**(-m)
    return ampl*h+y0


