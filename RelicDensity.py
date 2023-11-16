import numpy as np
import PlotSettings as ps
from PlotSettings import plt as plt
import LatexSettings


#Arrays
RelicData = np.loadtxt('Data-Sets/RelicDensityData.csv', delimiter=',')
DMmass = RelicData[:,0]
Density = RelicData[:,1:5].T
ColorArray = [ps.MainColor1, ps.MainColor2, ps.MainColor3, ps.Gray1]
LabelArray = [r'$g^{\prime}_{\chi} = 8\times 10^{-3}$', r'$g^{\prime}_{\chi} = 1\times 10^{-2}$', r'$g^{\prime}_{\chi} = 2\times 10^{-2}$', r'$g^{\prime}_{\chi} = 4 \times 10^{-2}$']


for i in range(0,4):
    plt.plot(DMmass, Density[i], color=ColorArray[i], label=LabelArray[i])

plt.hlines(1, xmin=DMmass[0], xmax=DMmass[-1], colors='k')

plt.yscale('log')
plt.xscale('log')
plt.xlabel('$m_{\chi} [\mathrm{GeV}]$')
plt.ylabel('$ \Omega_{\chi} h^2/ \Omega_{DM} h^2$')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig('Plots/RelicDensityPlotpy.pdf')
