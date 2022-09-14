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
from primordial_functions import P, pps


from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import rc
import classy

filename = "COM_PowerSpect_CMB-TT-full_R3.01.txt"

file = open(filename, "r")

#Remove the Header
dump = file.readline()

#Read in the actual data
exp = np.empty(2507)
errors = np.empty((2,2507))
count = 0
no_l=1000#29
for line in file:
    exp[count] = float(line.split()[1])
    errors[0,count] = float(line.split()[2])
    errors[1,count] = float(line.split()[3])
    count += 1

def plot_1(universes,nts,cs2m,cs2p,experim):
    #REVTEX Params
    set_matplotlib_formats('pdf', 'svg')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    matplotlib.use('PDF')

    plt.rcParams['axes.labelsize']=10
    plt.rcParams['legend.fontsize']=8
    plt.rcParams['xtick.labelsize']=8
    plt.rcParams['ytick.labelsize']=8
    plt.rcParams['text.usetex']=True
    plt.rcParams['font.family']="serif"
    plt.rcParams['font.serif']="cm"

    errorlinewidth = 0.5
    errorwidth = 0.2

    #Planck Results

    ## Plotting
    m=(3*cs2m+1)/2;
    eta_max=np.pi/(2*m)
    nts=eta_max*np.array(nts)

    Cherry_Red = '#900001'
    Bright_Red = '#fe3533'
    Warm_Green = '#8db046'
    Ultramarine = '#4166f5'
    Deep_Blue = '#1d3e6e'
    Black = 'k'
    Light_Grey = '#A0A0A0BF'

    # lclosed = (closed00_l,closed03_l,closed0_l,closed_l,closedAs_l)
    # TTclosed = (closed00_TT,closed03_TT,closed0_TT,closed_TT,closedAs_TT)
    # lopen = (open00_l,open03_l,open0_l,open_l,openAs_l)
    # TTopen = (open00_TT,open03_TT,open0_TT,open_TT,openAs_TT)
    # lflat = (flat00_l,flat03_l,flat0_l,flat_l,flatAs_l)
    # TTflat = (flat00_TT,flat03_TT,flat0_TT,flat_TT,flatAs_TT)

    zorders = (1.0,1.0,1.0,1.0,1.0)

    cols = ('C0','C1','C2','C3','grey')

    errorlinewidth = 0.5
    errorwidth = 0.2
    fig=plt.figure(figsize = (7.0, 7.0))
    plt.plot([])
    #fig.axis('off')
    ax1=fig.axes[0]
    ax2=ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.plot([])

    # ax2.tick_params(axis='y')
    ax1.set_ylabel(r'$\log(10^{10} \mathcal{P}_\mathcal{R})$', fontsize = 10)  # we already handled the x-label with ax1
    ax2.set_ylabel(r'$\mathcal{D}_{\ell}^{TT}$ [$\mu K^2$]', fontsize = 10)  # we already handled the x-label with ax1

    ##Remove ugly things
    ax1.spines['top'].set_color('white')
    ax1.spines['right'].set_color('white')
    ax1.spines['bottom'].set_color('white')
    ax1.spines['left'].set_color('white')

    ax2.spines['top'].set_color('white')
    ax2.spines['right'].set_color('white')
    ax2.spines['bottom'].set_color('white')
    ax2.spines['left'].set_color('white')

    ax1.tick_params(axis='y', colors='white')
    ax1.tick_params(axis='x', colors='white')

    ax2.tick_params(axis='y', colors='white')
    ax2.tick_params(axis='x', colors='white')


    axs = fig.subplots(len(universes), 2)
    # color = 'tab:blue'


    #
    if((cs2m==cs2p)):
        fig.suptitle(r'$c_-^2 = c_+^2 = '+r'{}'.format('{0:2}'.format(cs2m)+r'$'))
    else:
        fig.suptitle(r'$c_-^2 = '+r'{}'.format('{0:2}'.format(cs2m)+r'$')+r'; $c_+^2 = '+r'{}'.format('{0:2}'.format(cs2p)+r'$'))
    for i,universe in enumerate(universes):

        secax = axs[i, 0].secondary_xaxis('top', functions = (universe.comk, universe.physk))
        if(i==0):secax.set_xlabel(r'physical $k$ [Mpc$^{-1}$]', fontsize = 10, labelpad = 10)

        secax = axs[i, 1].secondary_xaxis('top', functions = (universe.comk, universe.physk))
        if(i==0):secax.set_xlabel(r'physical $k$ [Mpc$^{-1}$]', fontsize = 10, labelpad = 10)

        #Left plot
        if(universe.K==1): k = np.logspace(np.log10(3), 4, 1000)
        else: k = np.logspace(np.log10(2), 4, 1000)

        for nt in np.array(nts):
            axs[i, 0].semilogx(k, np.log((10**10)*pps(universe,k,k,nt,cs2m,cs2p)), label = r'$\eta_{\mathrm{t}} = '+r'{}'.format('{}'.format(round(nt/eta_max,2)))+r'\eta_{\mathrm{max}}$', alpha = 0.5)

        k = np.logspace(np.log10(2), 4, 1000)

        axs[i, 0].semilogx(k, np.log(pps(universe,k,k,nt,cs2m,cs2p,baseline=True)*(10**10)), linestyle = '--', color = 'grey', label = r'$\mathcal{P}_\mathcal{R}^{\,\,K\Lambda\mathrm{CDM}}$')
        if(i==len(universes)-1):axs[i, 0].set_xlabel(r'comoving $k$', fontsize = 10)
        # fig.supxlabel(r'comoving $k$', fontsize = 10)
        # axs[i, 0].set_ylabel(r'$\log(10^{10} \mathcal{P}_\mathcal{R})$', fontsize = 10)
        # fig.supylabel(r'$\log(10^{10} \mathcal{P}_\mathcal{R})$', fontsize = 10)
        axs[i, 0].set_xlim([2, 10**4])
        # axs[i, 0].set_ylim([2, 3.1])
        # axs[i, 0].set_ylim([2,3.6])
        axs[i, 0].set_ylim([0.5,5.0])
        axs[i, 0].set_yticks([2.2, 2.6, 3, 3.4])
        # axs[i, 0].set_yticks([1.2,1.6,2.0,2.4,2.8,3.2,3.6])
        axs[i, 0].text(x = 1000, y = 2.9, s = r'$K = {}$'.format(universe.K), fontsize = 10)

        if(universe.K==1):
            k = np.arange(3, 500, 1)
            for nt,col in zip(nts,cols[:-1]):
                axs[i, 0].semilogx(k, np.log((10**10)*pps(universe,k,k,nt,cs2m,cs2p)), ".", markersize = 5, alpha = 0.7, color = col)
            axs[i, 0].vlines(x = 3, ymin = 0.8, ymax = 4, linestyle = '--', lw = 0.5, color = 'black')
        ls=[]
        TTs=[]
        for j,nt in enumerate(nts):
            l_1,TT_1=universe.gen_CMB(nt,cs2m,cs2p,'uebak'+str(j),baseline=False)
            ls.append(l_1)
            TTs.append(TT_1)
        base_l,base_TT=universe.gen_CMB(nt,cs2m,cs2p,'uebak'+str(j),baseline=True)
        ls.append(base_l)
        TTs.append(base_TT)
        #Right plot
        axs[i, 1].errorbar(base_l[0:no_l], exp[0:no_l], errors[:,0:no_l], fmt = 'o', mfc = Light_Grey, elinewidth=errorlinewidth,markeredgewidth = errorwidth, color = Light_Grey, ecolor = Light_Grey, label = 'Planck 2018', zorder = 0.5)
        # fig.supxlabel('commonx')
        for l,TT,col,zo in zip(ls,TTs,cols,zorders):
            axs[i, 1].semilogx(l[0:no_l],TT[0:no_l],col,zorder = zo,alpha=0.7)
            print(l[0:no_l].min(),l[0:no_l].max())
            print(TT[0:no_l].min(),TT[0:no_l].max())
        axs[i, 1].plot(100, 100, color = "grey", label = r'$\mathcal{P}_{\mathcal{R}}^{K\Lambda CDM}$')
        # axs[i, 1].legend(loc = 'best', prop = {'size': 6})#1
        axs[i, 1].set_xticks([2,10,30])
        axs[i, 1].set_xticklabels([2,10,30])
        # axs[i, 1].set_ylim([0,2000])
        axs[i, 1].set_ylim([0,8000])
        # axs[i, 1].set_yticks([500,1000,1500,2000])
        axs[i, 1].set_yticks([500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000])
        axs[i, 1].set_xlim([1.5,no_l+1])
        if(i==len(universes)-1):axs[i, 1].set_xlabel(r'$\ell$', fontsize = 10)
        # axs[i, 1].set_ylabel(r'$\mathcal{D}_{\ell}^{TT}$ [$\mu K^2$]', fontsize = 10)
        axs[i,1].yaxis.set_label_position("right")
        axs[i,1].yaxis.tick_right()

        # fig.set_tight_layout({'rect': [0, 0.03, 1, 0.95]})
        # fig.text(0.5, 0.03, 'Position', ha='center')

        # fig.supylabel(r'$\mathcal{D}_{\ell}^{TT}$ [$\mu K^2$]', fontsize = 10)
        axs[i, 1].text(2, 1650, r'$K = {}$'.format(universe.K), fontsize = 10)
        if(i==1):
            axs[i, 0].legend(prop = {'size': 6},loc='upper left')#4
            axs[i, 1].legend(prop = {'size': 6},loc='upper left')#1

    plt.subplots_adjust(hspace = 0.6, wspace = 0.4)
    plt.tight_layout()
    os.makedirs('./figures', exist_ok=True)
    plt.savefig('./figures/cl_{}.pdf'.format(experim))
