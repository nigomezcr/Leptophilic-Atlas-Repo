import numpy as np
import matplotlib.pyplot as plt
from PlotSettings import MainColor1, MainColor2, MainColor3, BackgroundColor1, BackgroundColor2, BackgroundColor3, Gray1

# Arrays
mDM_array = np.logspace(0, 3, 100)
v0 = (30, 1000)

alphaDM = 0.01
mZp, mZp2, mZp3 = 19.9, 25.1, 22
mDM, mDM2, mDM3 = 161.5, 209.6, 25

######## Functions ############
def sigmatransfer(V, mphi, mchi, alphax):
    w = 300 * (mphi / 10) * (10 / mchi)
    st = (275.73) * (alphax / 1e-2) ** 2 * (mchi / 10) * (10 / mphi) ** (4)
    sv = 2 * st * (w ** 4 / V ** 4) * (2 * np.log(1.0 + V ** 2 / (2 * w ** 2)) - np.log(1.0 + V ** 2 / (w ** 2)))
    return sv

def sigmaviscosity(V, mphi, mchi, alphax):
    w  = 300*(mphi/(10))*(10/(mchi))
    beta = V/w
    st = (275.73)*(alphax/0.01)**2*(mchi/10.0)*(10.0/(mphi))**(4)

    if beta < 0.01:
        sv = 3/6*st
    else:
        sv = 3*st*(1/beta**6)*( (2+beta**2)*np.log(1+ beta**2)- 2*beta**2 )
    return sv

############# Filling arrays for the plots #############
Dwarfv = [sigmaviscosity(v0[0], mZp, m, alphaDM) for m in mDM_array]
Dwarft = [sigmatransfer(v0[0], mZp2, m, alphaDM) for m in mDM_array]

Clusterv = [sigmaviscosity(v0[1], mZp, m, alphaDM) for m in mDM_array]
Clustert = [sigmatransfer(v0[1], mZp2, m, alphaDM) for m in mDM_array]

###### For the Plots ######
line_colors = (MainColor1, MainColor3)
line_styles = ('solid', 'dashed', 'dotted')

def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(c1)
    c2=np.array(c2)
    return (1-mix)*c1 + mix*c2

n = 500
c1 = MainColor2  # blue
c2 = (1, 1, 1)  # white
y_values = np.logspace(-2, 0, n+1)  # Generate logarithmically spaced values

#Plots
fig, ax = plt.subplots(1, 1)

ax.text(mDM_array[5], 0.2, 'Clusters', fontsize=16, color=BackgroundColor2)
for y in y_values:
    ax.axhline(y, 0.0, 1, color=colorFader(c2, c1, (np.log10(y) + 3) / 4), linewidth=1)

ax.text(mDM_array[5], 22, 'Dwarfs', fontsize=16, color=BackgroundColor1)
ax.fill_between([mDM_array[0], mDM_array[-1]], 10, 100, color=BackgroundColor1, alpha=0.7)

ax.set_ylabel(r'$\sigma /m_{\chi} ~[ \mathrm{cm^2/g}]$')
ax.set_xlabel(r'$m_{{\chi}} ~[\mathrm{GeV} ] $')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([mDM_array[0], mDM_array[-1]])
ax.set_ylim(1e-2, 200)


DwPlotv = ax.plot(mDM_array, Dwarfv, linestyle='solid', color=BackgroundColor3)
DwPlott = ax.plot(mDM_array, Dwarft, linestyle='solid', color=BackgroundColor1)

ClPlotv = ax.plot(mDM_array, Clusterv, linestyle='dashed',  color=BackgroundColor3)
ClPlott = ax.plot(mDM_array, Clustert, linestyle='dashed', color=BackgroundColor2)


ax.vlines(([20, 200]), ([1e-2, 1e-2]), ([200, 200]), colors=('k', 'k'), linestyles=['dotted', 'dotted'])
first_legend = ax.legend(handles=[DwPlotv[0], DwPlott[0]], labels=['Viscosity', 'Transfer'], loc=4, title='Models')
ax.add_artist(first_legend)

ax.legend(handles=[DwPlott[0], ClPlott[0]], labels=[r'$v = 30~km/s $', r'$v=1000~km/s$'], loc=3, title='Velocity')




fig.tight_layout()
fig.savefig('Plots/TransferCrossSections.pdf')
