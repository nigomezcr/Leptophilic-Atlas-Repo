import  numpy as np
from Constants import c, fc
import matplotlib.pyplot as plt

def sigmaviscosity(V, mphi=19.9, mchi=16, alphax=0.01):
    w  = 300*(mphi/(10))*(10/(mchi))
    st = (275.73)*(alphax/0.01)**2*(mchi/10.0)*(10.0/(mphi))**(4)
    sv = 3*st*(w**6/(V**6))*( (2+V**2/(w**2))*np.log(1+V**2/(w**2))-2*V**2/w**2)
    return sv

def sigmatransfer(V, mphi=19.9, mchi=160, alphax=0.01):
    w = 300 * (mphi / 10) * (10 / mchi)
    st = (275.73) * (alphax / 1e-2) ** 2 * (mchi / 10) * (10 / mphi) ** (4)
    sv = 2 * st * (w ** 4 / V ** 4) * (2 * np.log(1.0 + V ** 2 / (2 * w ** 2)) - np.log(1.0 + V ** 2 / (w ** 2)))
    return sv

def sigmaviscosity_alter(V, mphi=19.9, mchi=16, alphax=0.01):
    w  = 300*(mphi/(10))*(10/(mchi))
    beta = V/w
    st = (275.73)*(alphax/0.01)**2*(mchi/10.0)*(10.0/(mphi))**(4)

    if beta < 0.01:
        sv = 3/6*st
    else:
        sv = 3*st*(1/beta**6)*( (2+beta**2)*np.log(1+ beta**2)- 2*beta**2 )
    return sv



v_array = np.logspace(0.2, 1, 100)
crossSection_array_transfer = [sigmaviscosity_alter(v) for v in v_array]
crossSection_array_viscosity = [sigmaviscosity(v) for v in v_array]

plt.plot(v_array, crossSection_array_transfer, 'k')
plt.plot(v_array, crossSection_array_viscosity, 'g')



plt.xscale('log')
plt.yscale('log')
plt.ylabel('Cross Section')
plt.xlabel('velocity')
plt.savefig('Plots/CrossSectionsTest.pdf')
#plt.show()
