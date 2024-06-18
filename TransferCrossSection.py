import numpy as np
import matplotlib.pyplot as plt
from PlotSettings import MainColor1, MainColor2, MainColor3, BackgroundColor1, BackgroundColor2, BackgroundColor3, Gray1

# Arrays
mDM_array = np.logspace(-1, 2, 100)
v0 = (10, 1000)

alphaDM, alphaDM2, alphaDM3 = 0.005, 0.01, 0.03
mZp, mZp2, mZp3 = 9.8, 20, 15
mDM, mDM2, mDM3 = 6.9, 165, 8


def sigmatransfer(V, mphi, mchi, alphax):
    w = 300 * (mphi / 10) * (10 / mchi)
    st = (275.73) * (alphax / 1e-2) ** 2 * (mchi / 10) * (10 / mphi) ** (4)
    sv = 2 * st * (w ** 4 / V ** 4) * (2 * np.log(1.0 + V ** 2 / (2 * w ** 2)) - np.log(1.0 + V ** 2 / (w ** 2)))
    return sv

###### For the Plots ######

line_colors = (MainColor1, MainColor3)
line_styles = ('solid', 'dashed', 'dotted')

def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(c1)
    c2=np.array(c2)
    return (1-mix)*c1 + mix*c2

n = 500
c1 = BackgroundColor2  # yellow
c2 = (1, 1, 1)  # white
y_values = np.logspace(-1.6, 0, n+1)  # Generate logarithmically spaced values


#Plots
fig, ax = plt.subplots(1, 1)

ax.text(mDM_array[5], 0.2, 'Clusters', fontsize=14)
for y in y_values:
    ax.axhline(y, 0.0, 1, color=colorFader(c2, c1, (np.log10(y) + 3) / 4), linewidth=4)

ax.text(mDM_array[5], 20, 'Dwarfs')
ax.fill_between([mDM_array[0], mDM_array[-1]], 10, 100, color=BackgroundColor1, alpha=0.7)

ax.set_ylabel(r'$\sigma_T /m_{\chi} ~[ \mathrm{cm^2/g}]$')
ax.set_xlabel(r'$m_{{\chi}} ~[\mathrm{GeV} ] $')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([mDM_array[0], mDM_array[-1]])
ax.set_ylim(1e-2, 200)

Dwarf1 = [sigmatransfer(v0[0], mZp, m, alphaDM) for m in mDM_array]
Dwarf2 = [sigmatransfer(v0[0], mZp2, m, alphaDM2) for m in mDM_array]
Dwarf3 = [sigmatransfer(v0[0], mZp3, m, alphaDM3) for m in mDM_array]

Cluster1 = [sigmatransfer(v0[1], mZp, m, alphaDM) for m in mDM_array]
Cluster2 = [sigmatransfer(v0[1], mZp2, m, alphaDM2) for m in mDM_array]
Cluster3 = [sigmatransfer(v0[1], mZp3, m, alphaDM3) for m in mDM_array]


DwPlot1 = ax.plot(mDM_array, Dwarf1, linestyle='solid', color=MainColor1, label='DM1')
DwPlot2 = ax.plot(mDM_array, Dwarf2, linestyle='dashed', color=MainColor1, label='DM2')
DwPlot3 = ax.plot(mDM_array, Dwarf3, linestyle='dotted', color=MainColor1, label='DM3')

ClPlot1 = ax.plot(mDM_array, Cluster1, linestyle='solid',  color=MainColor2, label='DM1')
ClPlot2 = ax.plot(mDM_array, Cluster2, linestyle='dashed', color=MainColor2, label='DM2')
ClPlot3 = ax.plot(mDM_array, Cluster3, linestyle='dotted', color=MainColor2, label='DM3')

first_legend = ax.legend(handles=[DwPlot1[0], DwPlot2[0], DwPlot3[0]], loc=4, title='Models')
ax.add_artist(first_legend)

ax.legend(handles=[DwPlot1[0], ClPlot1[0]], loc=2, title='Velocity')




fig.tight_layout()
fig.savefig('Plots/TransferCrossSections.pdf')
