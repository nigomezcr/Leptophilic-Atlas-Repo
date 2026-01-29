import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import PlotSettings as ps

from matplotlib.lines import Line2D

"""
////////// Constants and Patameters /////////////////////////////////
"""

# Constants
mm = .105 # muon mass in GeV
me = 0.000511 # electron mass in GeV
Delta_amu = 38 # discrepancy in a_mu 10^{-11} units
sigma_amu = 63
Delta_ae = 0.087 # discrepancy in a_e in 10^{-11} units
sigma_ae = 0.036

#Parameters for the plots (charges)
rho = np.array([0, -1/6, -1/2, -1/20])
#nu = np.array([1, 1/3, 1/2, 1/4])

Qe = rho
Qmu = np.array([1, 1/2, 1/4, 1])#nu - rho



"""
////////// Import data on constrains //////////////////////////
"""

# BaBar
BaBar_data = np.loadtxt("Data-Sets/Babar_exclusion.csv", delimiter=',')
Babar_x = BaBar_data[:,0]
Babar_y = BaBar_data[:,1]

# CCFR
CCFR_data = np.loadtxt("Data-Sets/Nu-trident_exclusion.csv", delimiter=',')
CCFR_x = CCFR_data[:,0]
CCFR_y = CCFR_data[:,1]

# Borexino
Borexino_data = np.loadtxt("Data-Sets/Borexino_exclusion.csv", delimiter=',')
Borexino_x = Borexino_data[:,0]
Borexino_y = Borexino_data[:,1]

"""
////////// Coupling constant function /////////////////////////////////
"""

def gp(charge, delta_a, M, m):
    if charge == 0:
        return 0

    integ = quad(lambda x: m**2*x**2*(1-x)/(x**2*m**2 + (1-x)*M**2), 0, 1)[0]
    gp2 = 4*np.pi**2*(delta_a)/( np.abs(charge)*integ )
    return np.sqrt(gp2)

"""
///////////  Arrays For the Plots  //////////////////////////
"""

#Z' mass array
N=1000
mZp = np.linspace(1, 7000, N)*10**(-3) # Mass in GeV

#Create arrays
g_muon = np.full([4, len(mZp)], 0.)
g_electron = np.full([4, len(mZp)], 0.)

for i in range(1, 4):
    g_muon[i] =  [gp(Qmu[i], Delta_amu, M, mm)*10**(-5.5) for M in mZp]
    g_electron[i] =  [gp(Qe[i], Delta_ae, M, me)*10**(-5.5) for M in mZp]

"""
/////////// Plotting Functions ///////////////////////////////////
"""
def plot_Settings():
    plt.ylim(1e-4, 1e-1)
    plt.xlim(1e-3, mZp[-1])
    plt.xlabel('$m_{Z\'} ~[\mathrm{GeV}]$')
    plt.ylabel(' ${g\'}$ ')
    plt.yscale('log')
    plt.xscale('log')
    plt.tight_layout()


def plot_Constraints():
    # Borexino
    #BorexinoColor = ps.Gray2
    #plt.plot(Borexino_x,Borexino_y,'-',color=BorexinoColor, alpha=1)
    #plt.fill_between(Borexino_x, Borexino_y, 1e-1, color=BorexinoColor,alpha=ps.RegionAlpha)
    #plt.text(1e-2,5e-3, 'Borexino', fontsize=12, color=BorexinoColor)

    # CCFR
    CCFRColor = ps.MainColor2
    plt.plot(CCFR_x,CCFR_y,'-',color=CCFRColor, alpha=ps.RegionAlpha)
    plt.fill_between(CCFR_x, CCFR_y, 1e-1, color=CCFRColor,alpha=ps.RegionAlpha)
    plt.text(1.5e-2,3e-2, 'CCFR', fontsize=12, color=CCFRColor)

    #  BaBar
    BaBarColor = ps.Gray1
    plt.fill_between(Babar_x, Babar_y, 1e-1, color=BaBarColor, alpha=ps.RegionAlpha)
    plt.loglog(Babar_x,Babar_y, color=BaBarColor, alpha=0.7)
    plt.text(1.3, 3e-2, 'BaBar', fontsize=12, color=ps.Gray2)

    #  Neff
    NeffColor = ps.Gray2
    plt.fill_between(([1e-3,5.3e-3]),1e-4,1e-1, color=NeffColor, alpha=ps.RegionAlpha)
    plt.text(1.4e-3 , 3e-2, r'$\Delta N_{eff}$', fontsize= 12,color=ps.Gray1)

    


def plot_Atlas():

    gColor2 = ps.MainColor1
    gColor1 = ps.BackgroundColor1

    #Fill g-2 region for Lmu - Ltau, left plot
    g_sigma1p = [gp(Qmu[0], Delta_amu , M, mm)*10**(-5.5) for M in mZp ]
    #g_sigma1m = [gp(Qmu[0], Delta_amu - 1*sigma_amu, M, mm)*10**(-5.5) for M in mZp ]
    


    #plt.fill_between(mZp, g_sigma1m, g_sigma1p, color=gColor1)

    # Plot the other models
    lines = ['solid', 'dashed', 'dotted']

    


    custom_lines = [Line2D([0], [0]),
                    Line2D([0], [0]),
                    Line2D([0], [0])]


    for i in range(1, 4):
        custom_lines[i-1] = Line2D([0], [0], color='k', lw=3, linestyle=lines[i-1], label='$(g-2)_e \; BM{:n}$'.format(i+1))
        plt.loglog(mZp, g_muon[i],  color=gColor1, linestyle=lines[i-1], label=r'$(g-2)_e \; BM{:n}$'.format(i+1))
        plt.loglog(mZp, g_electron[i],  color='k', linestyle=lines[i-1], label=r'$(g-2)_e \; BM{:n}$'.format(i+1))


    g_patch = ps.mpatches.Patch(color=gColor1, label=r'$1\sigma ~(L_{\mu} - L_{{\tau}})$')
    
 


    plt.legend(handles=[custom_lines[0], custom_lines[1], custom_lines[2], g_patch, g_patch2], loc='lower right') 

    plt.savefig("Plots/g_minus_2-final.svg")
    plt.savefig("Plots/g_minus_2-final.pdf")


"""
/////////// Generate Plots ///////////////////////////
"""


plot_Settings()
plot_Constraints()
plot_Atlas()

#PrintContributions()
