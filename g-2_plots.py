import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import PlotSettings as ps

"""
////////// Constants and Patameters /////////////////////////////////
"""

# Constants
mm = .105 # muon mass in GeV
me = 0.000511 # electron mass in GeV
Delta_amu = 248 # discrepancy in a_mu 10^{-11} units
sigma_amu = 49
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
    CCFRColor = ps.MainColor1
    plt.plot(CCFR_x,CCFR_y,'-',color=CCFRColor, alpha=ps.RegionAlpha)
    plt.fill_between(CCFR_x, CCFR_y, 1e-1, color=CCFRColor,alpha=ps.RegionAlpha)
    plt.text(1e-2,3e-2, 'CCFR', fontsize=12, color=CCFRColor)

    #  BaBar
    BaBarColor = ps.Gray1
    plt.fill_between(Babar_x, Babar_y, 1e-1, color=BaBarColor, alpha=ps.RegionAlpha)
    plt.loglog(Babar_x,Babar_y, color=BaBarColor, alpha=0.7)
    plt.text(1, 3e-2, 'BaBar', fontsize=12, color=BaBarColor)

    #  Neff
    NeffColor = ps.BackgroundColor2
    plt.fill_between(([1e-3,5.3e-3]),1e-4,1e-1, color=NeffColor, alpha=ps.RegionAlpha)
    plt.text(1.4e-3 , 3e-2, r'$\Delta N_{eff}$', fontsize= 12,color=(155/255, 8/255, 0/255))

    


def plot_Atlas(mass_arr, g_arr, color='k', lsty='solid', label=None):

    gColor1 = ps.MainColor3
    gColor2 = ps.MainColor2

    #Fill g-2 region for Lmu - Ltau, left plot
    g_sigma1p = [gp(Qmu[0], Delta_amu + 1*sigma_amu, M, mm)*10**(-5.5) for M in mZp ]
    g_sigma1m = [gp(Qmu[0], Delta_amu - 1*sigma_amu, M, mm)*10**(-5.5) for M in mZp ]
    g_sigma2p = [gp(Qmu[0], Delta_amu + 2*sigma_amu, M, mm)*10**(-5.5) for M in mZp ]
    g_sigma2m = [gp(Qmu[0], Delta_amu - 2*sigma_amu, M, mm)*10**(-5.5) for M in mZp ]

    plt.title(r'$L_{\mu} - L_{\tau}$ model')
    plt.fill_between(mZp, g_sigma2m, g_sigma2p, color=gColor2)
    plt.fill_between(mZp, g_sigma1m, g_sigma1p, color=gColor1)

    line = plt.plot(mass_arr, g_arr, color=color, linestyle=lsty, label=label)

    g_patch = ps.mpatches.Patch(color=gColor1, label=r'$1\sigma ~(L_{\mu} - L_{\tau})$')
    g_patch2 = ps.mpatches.Patch(color=gColor2, label=r'$2\sigma  ~(L_{\mu} - L_{\tau})$')
    label_patch = ps.mpatches.Patch(color='k', label=label)


    plt.legend(handles=[g_patch, g_patch2, line], loc='lower right')
    plt.savefig("Plots/g_minus_2-final.svg")
    plt.savefig("Plots/g_minus_2-final.pdf")


def PrintContributions():
    print("The Contribution for the charged fermions anomalous magnetic moment:\n")
    print("$(g-2)_{e}$=\t", gp(Qe[2], 0, 10, me))
    print("$(g-2)_{\mu}$=\t",  gp(Qmu[2], 0, 10, mm))
    #print("$(g-2)_{\tau}$=\t",  gp(Qtau, 0, 10, mtau))

"""
/////////// Generate Plots ///////////////////////////////////
"""

lines = ['solid', 'dashed', 'dotted']
colors = [ps.MainColor1, ps.MainColor2, ps.MainColor3]

plot_Settings()
#plot_Constraints()
#plot_Lmu_minus_Ltau()
for i in range(1, 4):
    plot_Atlas(mZp, g_muon[i],  color=ps.MainColor2, lsty=lines[i-1], label='BM{:n}'.format(i+1))
#    plot_Atlas(mZp, g_electron[i],  color='k', lsty=lines[i-1])

#PrintContributions()
print('----> Done!')
