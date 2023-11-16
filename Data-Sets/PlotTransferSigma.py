import  numpy as np
from Constants import c, fc
import matplotlib.pyplot as plt
from CrossSections import Transfer_Sigma_Low_Energy, Transfer_sigma


v_array = np.logspace(1, 3, 100)
crossSection_array_low = [Transfer_Sigma_Low_Energy(v/c) for v in v_array]
crossSection_array = [Transfer_sigma(v/c) for v in v_array]


plt.plot(v_array, crossSection_array_low, 'k')
plt.plot(v_array, crossSection_array)


plt.xscale('log')
plt.yscale('log')
plt.ylabel('Transfer Cross Section')
plt.xlabel('velocity')
plt.show()
