from primordial_functions import P, pps
#Dependencies
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from IPython.display import set_matplotlib_formats

import scipy.special as sp
import numpy as np
import sys
import math
import os
import warnings

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import rc
import classy


#
# Bulk of code based on https://zenodo.org/record/4024321#.YyHBViHMIbk
#

class universe():
    def __init__(self,K):
        self.K=K
        self.kp = 0.05
        klow = -7.0
        khigh = 2.0
        knum = 100000
        self.ks = np.logspace(klow,khigh,knum)
        ###Parameters estimated using anesthetic and/or taken from original code
        if(K==0):
            self.H=67.248066
            self.omega_k=0
            self.omega_cdm=0.120218
            self.omega_b=0.022342
            self.tau_reio=0.055270
            self.ns=0.964235
            self.As=2.104694*(10**(-9))

            self.h=self.H/100
            self.a= 43000 #Can be any arbitrary number really, as long as you rescale nt correctly as well.
            self.kcom=(self.ks*self.a)
        elif(K==1):
            self.H = 64.03
            self.omega_k = -0.0092
            self.omega_cdm = 0.11839
            self.omega_b = 0.022509
            self.tau_reio = 0.0515
            self.ns = 0.9699
            self.As = 2.077187691*(10**(-9))

            self.h=self.H/100
            self.a= 43000*((self.h/0.7)**(-1))*(np.abs(self.omega_k/0.01)**(-1/2))
            self.kcom=(self.ks*self.a)+3 ### Why is there a +3 here???


        elif(K==-1):
            self.H=68.165888
            self.omega_k=0.002252
            self.omega_cdm=0.120733
            self.omega_b=0.022326
            self.tau_reio=0.057484
            self.ns=0.963318
            self.As=2.116484*(10**(-9))

            self.h=self.H/100
            self.a= 43000*((self.h/0.7)**(-1))*(np.abs(self.omega_k/0.01)**(-1/2))
            self.kcom=(self.ks*self.a)

    def gen_CMB(self,nt,cs2m,cs2p,label,baseline=False):
        #Print out in required format
        Ps = np.nan_to_num(pps(self,self.ks,self.kcom,nt,cs2m,cs2p,baseline))
        np.savetxt(label,np.vstack([self.ks,Ps]).T)
        print(Ps.min(),Ps.max())

        cl = {}
        l = np.zeros(1999)
        TT = np.zeros(1999)


        curved = classy.Class()
        params = {
                 'output': 'tCl lCl',
                 'l_max_scalars': 2000,
                 'lensing': 'yes',
                 'P_k_ini type': 'external_Pk',
                 'command': 'cat '+label,
                 'tau_reio': self.tau_reio,
                 'h': self.h,
                 'omega_b': self.omega_b,
                 'Omega_k': self.omega_k,
                 'omega_cdm':self.omega_cdm}
        curved.set(params)
        curved.compute()
        cl.update(curved.lensed_cl(2000))
        l[:] = cl['ell'][2:]
        TT[:] = l*(l+1)*cl['tt'][2:] * (1e6 * 2.7255)**2 / (2*np.pi)


        return (l,TT)

    def comk(self,k):
        return k/self.a

    def physk(self,k):
        return k*self.a
