import  numpy as np
from scipy.integrate import quad
from Constants import c, fc, GeVtocm2

# Defult values of nice-working example
gp = 0.6
mZp = 20 # in MeV
mDM = 1 # in GeV


"""
///////////// Annihilation Cross Sections  ///////////////
"""


#Annihilation cross-section to leptons
def sigma_to_LL(s, g, M, m):
    Qc = 1
    Ql = 1
    Prefac = 4*np.pi*Qc*Qc*Ql
    alpha = g**2/(4*np.pi)
    return Prefac*alpha**2*s*(1 + 2*m**2/s)/((s - M**2)**2*np.sqrt(1 - 4*m**2/s))

#Thermal Cross section to leptons
def sigmaV_to_LL(gl, gc, ml, mc, M):
    term1 = gl**2*gc**2/(2*np.pi)
    term2 = np.sqrt(1 - ml**2/mc**2)
    term3 = ( 2*mc**2 + ml**2)/( (4*mc**2 - M**2)**2)
    return term1 * term2 * term3

#Thermal Cross section to Z's
def sigmaV_to_ZZ(gDM, mDM, M):
    term1 = gDM**4/(16*np.pi*mDM**2)
    term2 = (1 - M**2/mDM**2)**(3/2)
    term3 = (1 - M**2/(2*mDM**2))**(-2)
    return term1 * term2 * term3

# Approximated Thermal Cross section to Z's
def sigmaV_to_ZZ_Approx(g, m):
    term1 = g**4/(16*np.pi*m**2)
    return term1

#Non relativistic thermal cross-section to neutrinos
def sigvNR(g, M, m):
    Qc = 1
    return (Qc**4*g**4/np.pi)*(m**2/(4*m**2 - M**2)**2)


"""
///////////// Scattering Cross Sections Computed Analytically  ///////////////
"""


#### Total cross section, analitical #####
#  Particle-Particle: chi chi -> chi chi
def total_sigma_repulsive(v, g, M, m):
    d = 4 #Dirac
    s = 4*m**2/(1 - v**2)
    Prefactor = (2*g**4/(d*np.pi*s*m))
    Term1 = ((2*s + 3*M**2)*s*v**2 + 2*(M**2 + 2*m**2)**2)/(2*M**2*(M**2 + s*v**2))
    LogTerm1 = -(s*v**2*(3*M**2 + 4*m**2) + 2*(M**2 + 2*m**2)**2 - 4*m**4)/(s*v**2*(2*M**2 + s*v**2))
    return  fc*Prefactor*(Term1 + LogTerm1*np.log(1 + s*v**2/M**2) )


# Particle-Antiparticle: chi chibar -> chi chibar
def total_sigma_attractive(v, g, M, m):
    d = 4 #Dirac
    s = 4*m**2/(1 - v**2)
    Prefactor = g**4/(d*np.pi*s*m)
    sigma_t = Prefactor*( (s*v**2*(2*s + 3*M**2)+2*(M**2 + 2*m**2)**2)/(2*M**2*(M**2+s*v**2)) - (M**2 + s)/(s*v**2)*np.log(1 + s*v**2/M**2) )
    sigma_s = Prefactor/3*( (12*m**4 + 6*(v**2 + 1)*s*m**2 + s**2*v**4 ) / (s - M**2)**2 )
    sigma_st = -Prefactor/2*( (16*m**2+2*M**2+ 3*s*v**2)/(s-M**2) - (4*m**2 + 2*m**2*(4*M**2+3*s*v**2+s ) + (M**2+s*v**2)**2 )/(s*v**2*(s-M**2))*np.log(1 + s*v**2/M**2) )

    return fc*( sigma_s + sigma_t + sigma_st)

def total_sigma(v, g=gp, M=mZp ,m=mDM ):
    return (total_sigma_attractive(v, g, M, m) + total_sigma_repulsive(v, g, M, m))/2



def Normalized_Transer(v,  g=gp, M=mZp, m=mDM):
    alpha = g**2/(4*np.pi)
    w = 300*(M/10)*(10/m)
    sigma0T = 137.86*(alpha/0.01)**2*(m/10)*(10/M)**4
    return sigma0T*4*w**4/v**4 * (np.log(1 + v**2/(w**2)) - (v/w)**2/(1 + (v/w)**2) )


def Transfer_sigma_repulsive(v, g=gp, M=mZp ,m=mDM ):
    Beta =  v /(2)
    s = 4 * m**2 / (1 - Beta**2)
    sigma0 = -g**4/(4*np.pi*s**3*Beta**4)

    t_chn = -((s * Beta**2 * (-16 * m**4 - 6 * M**4 + 16 * m**2 * s - M**2 * s * (8 + 3 * Beta**2) + s**2 * (-4 - 4 * Beta**2 + Beta**4))) / (2 * (M**2 + s * Beta**2)))
    t_chn_log = - (8 * m**4 + 3 * M**4 - 8 * m**2 * s + 4 * M**2 * s + 2 * s**2)
    TR_t = t_chn + t_chn_log*np.log(1 + s*Beta**2/M**2)
    
    u_chn = -((s * Beta**2 * (48 * m**4 + 6 * M**4 + 9 * M**2 * s * Beta**2 + 2 * s**2 * (1 + Beta**4) + 16 * m**2 * (2 * M**2 + s * (-1 + Beta**2)))) / (2 * M**2))
    u_chn_log = (24 * m**4 + s**2 + 3 * (M**2 + s * Beta**2)**2 + 8 * m**2 * (2 * M**2 + s * (-1 + 2 * Beta**2)) )
    TR_u = u_chn + u_chn_log*np.log(1 + s*Beta**2/M**2)
    
    tu_chn =  (s * Beta**2) / (2 * M**2 + s * Beta**2) * (12 * m**4 - 8 * m**2 * s + s**2)
    TR_tu = tu_chn*np.log(1 + s*Beta**2/M**2)


    return fc* sigma0 /m  * ( TR_t )# + TR_u - 2*TR_tu  ) 

def Transfer_sigma_attractive(v, g=gp, M=mZp ,m=mDM ):
    Beta =  v /(2)
    s = 4 * m**2 / (1 - Beta**2)
    sigma0 = -g**4/(4*np.pi*s**3*Beta**4)

    t_chn = -((s * Beta**2 * (-16 * m**4 - 6 * M**4 + 16 * m**2 * s - M**2 * s * (8 + 3 * Beta**2) + s**2 * (-4 - 4 * Beta**2 + Beta**4))) / (2 * (M**2 + s * Beta**2)))
    t_chn_log = - (8 * m**4 + 3 * M**4 - 8 * m**2 * s + 4 * M**2 * s + 2 * s**2)
    TA_t = t_chn + t_chn_log*np.log(1 + s*Beta**2/M**2)
    
    TA_s = -s**2 * Beta**4 * (24 * m**4 + 4 * m**2 * s * Beta**2 + s**2 * (3 - 2 * Beta**2)) / (6 * (M**2 - s)**2)

    ts_chn = - s*Beta**2 * (-24*m**4 + 6*(M**2 + s)**2 - 3*s*Beta**2*(2*s + M**2) + 2*s**2*Beta**4)/(6*(s-M**2))
    ts_chn_log = M**2 * ((s + M**2)**2 - 4*m**4)/(s-M**2)
    TA_st = ts_chn + ts_chn_log*np.log(1 + s*Beta**2/M**2)

    return fc* sigma0 / m  * (TA_t )# + TA_s+ TA_st)


def Transfer_sigma(v, g=gp, M=mZp ,m=mDM ):
    M = M/1000
    return (Transfer_sigma_repulsive(v/c, g, M, m) + Transfer_sigma_attractive(v/c, g, M, m))/2


"""
///////////// Velocity Weighted Cross sections  ///////////////
"""

# Standard Transfer Cross Section Wighter

def sigv_integrand_low(v, v0, g, M, m):
    return Transfer_Sigma_Low_Energy(v/c, g, M, m)*(v/v0)*np.exp(-0.5*v**2/v0**2)*(v/v0)**2


def Transfer_Sigmav_Low_Energy(v0, g=gp, M=mZp, m=mDM):
    sigma2_MB = v0**2*np.pi*(3*np.pi - 8)/np.pi
    vmax = 2*np.sqrt(sigma2_MB)

    Prefactor = 4*np.pi/((2*np.pi)**1.5)
    Integral = quad(sigv_integrand_low, 0.0, vmax, args=(v0, g, M, m))[0]
    
    return Prefactor*Integral




# Normalized cross section

def Normalized_SigmaV(v0, g=gp, M=mZp, m=mDM):
    alpha = g**2/(4*np.pi)
    w = 300*(M/10)*(10/m)


    sigma0T = 137.86*(alpha/0.01)**2*(m/10)*(10/M)**4
    sigmaT = lambda v: 4*w**4/v**4 * (np.log(1 + v**2/(w**2)) - (v/w)**2/(1 + (v/w)**2) )

    sigma2_MB = v0**2*np.pi*(3*np.pi - 8)/np.pi
    vmax = 2*np.sqrt(sigma2_MB)
    integrand = lambda v: sigmaT(v)*(v/v0)*np.exp(-0.5*v**2/v0**2)*(v/v0)**2


    Prefactor = 4*np.pi*sigma0T/((2*np.pi)**1.5 )
    Integral = Prefactor*quad(integrand, 0.0, vmax)[0]

    
    return Integral


# PP Transfer Cross Section

def Transfer_SigmaV(v0, g=gp, M=mZp, m=mDM):
    M = M/1000
    Beta = v0/(2*c)
    s = 4 * m**2 / (1 - Beta**2)
    sigma0 = -g**4/(4*np.pi*s**3*Beta**4)
    
    # Attractive
    t_chn = -((s * Beta**2 * (-16 * m**4 - 6 * M**4 + 16 * m**2 * s - M**2 * s * (8 + 3 * Beta**2) + s**2 * (-4 - 4 * Beta**2 + Beta**4))) / (2 * (M**2 + s * Beta**2)))
    t_chn_log = - (8 * m**4 + 3 * M**4 - 8 * m**2 * s + 4 * M**2 * s + 2 * s**2)
    TR_t = t_chn + t_chn_log*np.log(1 + s*Beta**2/M**2)
    
    TA_s = -s**2 * Beta**4 * (24 * m**4 + 4 * m**2 * s * Beta**2 + s**2 * (3 - 2 * Beta**2)) / (6 * (M**2 - s)**2)

    ts_chn = - s*Beta**2 * (-24*m**4 + 6*(M**2 + s)**2 - 3*s*Beta**2*(2*s + M**2) + 2*s**2*Beta**4)/(6*(s-M**2))
    ts_chn_log = M**2 * ((s + M**2)**2 - 4*m**4)/(s-M**2)
    TA_st = ts_chn + ts_chn_log*np.log(1 + s*Beta**2/M**2)

    Transfer_sigma_attractive =  sigma0 / m  * (TR_t  + TA_s+ TA_st)

    # Repulsive
    u_chn = -((s * Beta**2 * (48 * m**4 + 6 * M**4 + 9 * M**2 * s * Beta**2 + 2 * s**2 * (1 + Beta**4) + 16 * m**2 * (2 * M**2 + s * (-1 + Beta**2)))) / (2 * M**2))
    u_chn_log = (24 * m**4 + s**2 + 3 * (M**2 + s * Beta**2)**2 + 8 * m**2 * (2 * M**2 + s * (-1 + 2 * Beta**2)) )
    TR_u = u_chn + u_chn_log*np.log(1 + s*Beta**2/M**2)
    
    tu_chn =  (s * Beta**2) / (2 * M**2 + s * Beta**2) * (12 * m**4 - 8 * m**2 * s + s**2)
    TR_tu = tu_chn*np.log(1 + s*Beta**2/M**2)

    Transfer_sigma_repulsive =  sigma0 / m  * ( TR_t  + TR_u - 2*TR_tu  ) 

    
    Transfer_sigma =  (Transfer_sigma_repulsive + Transfer_sigma_attractive)/2

    Integrand = lambda v: Transfer_sigma*(v/v0)*np.exp(-0.5*v**2/v0**2)*(v/v0)**2


    sigma2_MB = v0**2*np.pi*(3*np.pi - 8)/np.pi
    vmax = 2*np.sqrt(sigma2_MB)

    Prefactor = 4*np.pi/((2*np.pi)**1.5 )
    Integral = quad(Integrand, 0.0, vmax)[0]
    return fc* Prefactor*Integral


# Total Cross Section Weighted
def sigv_integrand(v, v0, g, M, m):
    return total_sigma(v, g, M, m)*v*np.exp(-0.5*v**2/v0**2)*v**2


def sigv_T(v0, g, M, m):
    sigma2_MB = v0**2*np.pi*(3*np.pi - 8)/np.pi
    vmax = 2*np.sqrt(sigma2_MB)

    Prefactor = 4*np.pi/((2*np.pi*v0**2)**1.5 * m)
    Integral = quad(sigv_integrand, 0.1, vmax, args=(v0, g, M, m))[0]
    return Prefactor*Integral



