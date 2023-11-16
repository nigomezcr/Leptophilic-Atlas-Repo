
"""
///////////// constants and Parameters ///////////////
"""
#Universal constants
c = 3e+5 #Light speed in km/s
G = 4.3e-3 # In units of pc/Ms * (km/s)^2
h = 0.674 # Hubble constant in units of km/s/Mpc
Mp = 1.22E19 # Planck mass in GeV
me = 0.5e-3 # Electron Mass in GeV
mm = 0.105 # Muon Mass in GeV
DensityFactor = 1.5E8 # Product of s0/rho_c h^2

#Conversion Factors
rho_c = 2.7754e-7*h**2 # In units of Ms/Mpc-3
GeVtocm2 = (1/5.06e13)**2 # cm^2
GeVtog = (1/1.78e-24)
fc =  GeVtocm2*GeVtog
