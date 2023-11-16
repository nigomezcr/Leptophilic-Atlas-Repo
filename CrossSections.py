import  numpy as np
from scipy.integrate import quad
from Constants import c, fc, GeVtocm2

# Defult values of nice-working example
gp = 0.6
mZp = 0.02
mDM = 1


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
def sigmaV_to_ZZ_Approx(g, m, M):
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


def Transfer_Sigma_Low_Energy(v, g=gp, M=mZp, m=mDM):
    R = m*v/M
    alph = g**2/(4*np.pi)
    return fc*8*np.pi*alph**2/(m**3 * v**4) * (np.log(1 + R**2) - R**2/(1 + R**2) )


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


    return fc* sigma0 * ( TR_t  + TR_u - 2*TR_tu  ) / m

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

    return fc* sigma0 * (TA_t + TA_s+ TA_st)/m


def Transfer_sigma(v, g=gp, M=mZp ,m=mDM ):
    return (Transfer_sigma_repulsive(v, g, M, m) + Transfer_sigma_attractive(v, g, M, m))/2

"""
///////////// Scattering Cross Sections Computed Numerically  ///////////////
"""

#### Total cross section, numerical #####
# t channel
def sigma_t_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(2*d*np.pi*s**2*v**2*m))
    integrand = lambda t: (t**2 + 2*s*t + 2*s**2*v**2 + 8*m**4)/(t-M**2)**2
    return fc*Prefactor*quad(integrand, -s*v**2, 0)[0]

# u channel
def sigma_u_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(2*d*np.pi*s**2*v**2*m))
    integrand = lambda t: (t**2 - 8*m**2*(s + t) + s**2 + 24*m**4)/(t + s*v**2 + M**2)**2
    return fc*Prefactor*quad(integrand, -s*v**2, 0)[0]

# t-u interference
def sigma_tu_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(2*d*np.pi*s**2*v**2*m))
    integrand = lambda t: (4*m**4 - 2*(v**2 + 3)*m**2*s + s**2)/((t-M**2)*(t + s*v**2 + M**2))
    return - fc *Prefactor * quad(integrand, -s*v**2, 0)[0]


# s-channel
def sigma_s_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(16*d*np.pi*s**2*v**2*m))
    integrand = lambda t:  (8*((2*m**2 + s*v**2 + t)**2 - s*(s-2*m**2) + (2*m**2 - t)**2) - 2*s*(4*(s-2*m**2) - 8*s))/(s - M**2)**2
    return fc *Prefactor * quad(integrand, -s*v**2, 0)[0]

# s-t interference
def sigma_st_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(16*d*np.pi*s**2*v**2*m))
    integrand = lambda t:  - 8*(4*m**4 + 6*m**2*s*v**2 + 2*m**2*s + 8*m**2*t + s**2*v**4 + 2*s*v**2*t + t**2)/((s-M**2)*(t-M**2))
    return fc *Prefactor * quad(integrand, -s*v**2, 0)[0]


# Attractive
def sigma_a_Num(v, g, M, m):
    return sigma_t_num(v, g, M, m) + sigma_s_num(v, g, M, m) + sigma_st_num(v, g, M, m)

# Repulsive
def sigma_r_Num(v, g, M, m):
    return sigma_t_num(v, g, M, m) + sigma_u_num(v, g, M, m) + 2*sigma_tu_num(v, g, M, m)


# Total effective
def sigma_eff_Num(v, g, M, m):
    return (sigma_a_Num(v, g, M, m))# + sigma_r_Num(v, g, M, m) )/2

#### Transfer cross section, numerical #####
# t channel
def sigmaT_t_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(2*d*np.pi*s**2*v**2*m))
    integrand = lambda t: t*(t**2 + 2*s*t + 2*s**2*v**2 + 8*m**4)/(t-M**2)**2
    return -2/(s*v**2) * fc*Prefactor*quad(integrand, -s*v**2, 0)[0]

# u channel
def sigmaT_u_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(2*d*np.pi*s**2*v**2*m))
    integrand = lambda t: t*(t**2 - 8*m**2*(s + t) + s**2 + 24*m**4)/(t + s*v**2 + M**2)**2
    return -2/(s*v**2) * fc*Prefactor*quad(integrand, -s*v**2, 0, )[0]

# t-u interference
def sigmaT_tu_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = -(g**4/(2*d*np.pi*s**2*v**2*m))
    integrand = lambda t: t*(4*m**4 - 2*(v**2 + 3)*m**2*s + s**2)/((t-M**2)*(t + s*v**2 + M**2))
    return  -2/(s*v**2) * fc *Prefactor * quad(integrand, -s*v**2, 0)[0]


# s-channel
def sigmaT_s_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(16*d*np.pi*s**2*v**2*m))
    integrand = lambda t:  t*(8*((2*m**2 + s*v**2 + t)**2 - s*(s-2*m**2) + (2*m**2 - t)**2) - 2*s*(4*(s-2*m**2) - 8*s))/(s - M**2)**2
    return -2/(s*v**2) * fc *Prefactor * quad(integrand, -s*v**2, 0)[0]

# s-t interference
def sigmaT_st_num(v, g, M, m):
    d = 4 #Dirac
    v = v/2
    s = 4*m**2/(1 - v**2)
    Prefactor = (g**4/(16*d*np.pi*s**2*v**2*m))
    integrand = lambda t:  - t*8*(4*m**4 + 6*m**2*s*v**2 + 2*m**2*s + 8*m**2*t + s**2*v**4 + 2*s*v**2*t + t**2)/((s-M**2)*(t-M**2))
    return -2/(s*v**2) * fc *Prefactor * quad(integrand, -s*v**2, 0)[0]


# Attractive Transfer
def sigmaT_a_Num(v, g, M, m):
    return sigmaT_t_num(v, g, M, m) + sigmaT_s_num(v, g, M, m) + sigmaT_st_num(v, g, M, m)

# Repulsive Transfer
def sigmaT_r_Num(v, g, M, m):
    return sigmaT_t_num(v, g, M, m) + sigmaT_u_num(v, g, M, m) + 2*sigmaT_tu_num(v, g, M, m)


# Transfer effective
def sigmaT_Num(v, g, M, m):
    return (sigmaT_a_Num(v, g, M, m) + sigmaT_r_Num(v, g, M, m) )/2



"""
///////////// Velocity Weighted Cross sections  ///////////////
"""



# Ensemble average cross section
def sigv_integrand(v, v0, g, M, m):
    return total_sigma(v, g, M, m)*v*np.exp(-0.5*v**2/v0**2)*v**2


def sigv_T(v0, g, M, m):
    sigma2_MB = v0**2*np.pi*(3*np.pi - 8)/np.pi
    vmax = 2*np.sqrt(sigma2_MB)

    Prefactor = 4*np.pi/((2*np.pi*v0**2)**1.5 * m)
    Integral = quad(sigv_integrand, 0.1, vmax, args=(v0, g, M, m))[0]
    return Prefactor*Integral
