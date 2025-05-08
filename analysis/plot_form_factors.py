import numpy as np
from scipy.integrate import simpson, quad
from scipy.special import spherical_jn
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.patheffects as pe
import csv

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

# Constants
Mproton = 0.938272
fm_to_invGeV = 5.068 # 1 fm = 5.068 GeV^-1
GeV_to_invfm = 5.068 # 1 GeV = 5.068 fm^-1

### Data from CXX file ###
Q_p = []
FF_p = []

Q_n = []
FF_n = []

Q_p_e = []
Q_p_m = []
FF_p_e = []
FF_p_m = []

Q_n_e = []
Q_n_m = []
FF_n_e = []
FF_n_m = []

Q2s_Fe = []
Us_Fe = []
Q_Fe_CXX = []
FF_Fe_CXX = []
Q_Fe = []
FF_Fe = []

with open('../csv/form_factors/iron.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        u = float(row[1])
        Q2s_Fe.append(q2)
        Us_Fe.append(u)

FF2_array_Fe = np.gradient(Us_Fe, Q2s_Fe)
for i in range(len(Q2s_Fe)):
    q = np.sqrt(Q2s_Fe[i])
    ff = np.sqrt(-FF2_array_Fe[i])
    Q_Fe_CXX.append(q)
    FF_Fe_CXX.append(ff)

with open('../csv/form_factors/FF_p_electric.csv') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        Q_p_e.append(float(row[0]))
        FF_p_e.append(float(row[1]))

with open('../csv/form_factors/FF_p_magnetic.csv') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        Q_p_m.append(float(row[0]))
        FF_p_m.append(float(row[1]))

with open('../csv/form_factors/FF_n_electric.csv') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        Q_n_e.append(float(row[0]))
        FF_n_e.append(float(row[1]))

with open('../csv/form_factors/FF_n_magnetic.csv') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        Q_n_m.append(float(row[0]))
        FF_n_m.append(float(row[1]))

with open('../csv/form_factors/FF_Fe_3Gp-orange_Alt.csv') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        Q_Fe.append(float(row[0]))
        FF_Fe.append(float(row[1]))

###################################
##### Woods-Saxon Form Factor #####
###################################

def form_factor_WoodsSaxon(A, Q2, a=0.523):
    """
    Woods-Saxon form factor obtained analytically using a symmetrized Fermi function. See Appendix A in 1807.10973.
    The parameter a is set to 0.523 fm.
    """
    r0 = 1.126 * A**(1/3)
    Q = np.sqrt(Q2) * GeV_to_invfm

    WS = (3 * np.pi * a) / (r0**2 + np.pi**2 * a**2)
    if Q == 0:
        WS *= (r0**2 + np.pi**2 * a**2) / (3 * np.pi * a)
    else:
        WS *= np.pi * a * np.tanh(np.pi * Q * a)**(-1) * np.sin(Q * r0) - r0 * np.cos(Q * r0) # The Taylor expansion of this and the next line give the result within the if statement
        WS /= Q * r0 * np.sinh(np.pi * Q * a)
    return WS


############################
##### 3pF Form Factor ######
############################

# Cutoff radius for the radial integral of the Woods-Saxon form factor
R_cutoff = 9.0 # fm

def rho_3pF(r,w=-0.19,c=3.73,z=0.62):
    """
    3 parameter Fermi nuclear density distribution. The parameters c and z and the input radius are in femtometers. The parameter w is dimensionless. The nuclear density distribution is also dimensionless.
    """
    result = 1+w*r**2/c**2
    result /= 1+np.exp((r-c)/z)
    return result

def form_factor_integrand(r,Q2,rho):
    """
    The Woods-Saxon form factor integrand dr r^2 sin(qr)/(qr) rho(r). r is in fm and Q2 is in GeV^2.
    """
    Q = np.sqrt(Q2) # [GeV]
    if Q == 0:
        return r**2 * rho(r)
    else:
        return r**2 * np.sin(Q*GeV_to_invfm*r) * rho(r) / (Q*GeV_to_invfm*r)

def form_factor(Q2,rho):
    num, _ = quad(form_factor_integrand, 0, R_cutoff, args=(Q2,rho))
    denom, _ = quad(form_factor_integrand, 0, R_cutoff, args=(0,rho))
    return num/denom

def form_factor_3pF(Q2):
    return form_factor(Q2,rho_3pF)

############################
##### Helm Form Factor #####
############################


def form_factor_Helm(Q2,A=40,r0=0.52,s=0.9):
    """
    Helm form factor calculation assuming r0 and s are in fm, and the input q is in GeV
    """
    R0 = np.sqrt((1.23*A**(1/3)-0.6)**2 + 7/3*np.pi**2*r0**2 - 5*s**2)
    Q = np.sqrt(Q2)*GeV_to_invfm
    if Q == 0:
        Helm = 1.0
    else:
        Helm = 3*spherical_jn(1,Q*R0)/(Q*R0)  # The Taylor expansion of j1(x)/x around x=0 is 1/3
    Helm *= np.exp(-Q2*(s*fm_to_invGeV)**2/2)
    return Helm


##############################################
##### Adapted Klein-Nystrand Form Factor #####
##############################################

def form_factor_KNad(Q2,r0=3.427,a=0.7):
    """
    Klein-Nystrand form factor calculation assuming r0 and a are in fm, and the input q is in GeV.
    For argon 40, r0 = 3.427 fm.
    """
    RA = np.sqrt(5/3*r0**2 - 10*a**2)
    RA = RA * fm_to_invGeV
    q = np.sqrt(Q2)
    if q == 0:
        KN = 1.0
    else:
        KN = 3*spherical_jn(1,q*RA)/(q*RA) # The Taylor expansion of j1(x)/x around x=0 is 1/3
    KN /= 1 + Q2*(a*fm_to_invGeV)**2
    return KN


################################
#### Calculate form factors ####
################################

form_factor_3pF = np.vectorize(form_factor_3pF)
form_factor_WoodsSaxon = np.vectorize(form_factor_WoodsSaxon)
form_factor_Helm = np.vectorize(form_factor_Helm)
form_factor_KNad = np.vectorize(form_factor_KNad)

Q_array = np.linspace(0,0.7,1000)
Q2_array = Q_array**2

FF_W_WS = np.abs(form_factor_WoodsSaxon(184,Q2_array))
FF_Ar_3pF = np.abs(form_factor_3pF(Q2_array))
FF_Ar_Helm = np.abs(form_factor_Helm(Q2_array))
FF_Ar_KNad = np.abs(form_factor_KNad(Q2_array))

##########################################
#### Calculate form factor (u1) grids ####
##########################################

def u1(Q2,form_factor_func):
    f2 = lambda x: form_factor_func(x)**2 
    uval, _ = quad(f2,Q2,np.inf)
    return uval

u1 = np.vectorize(u1)

Q2_grid = []
### Read Q2 grid values from CXX file ###
with open('../csv/form_factors/argon.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q2_grid.append(q2)

Q2_grid = np.array(Q2_grid)

U1_Ar_3pF = u1(Q2_grid,form_factor_3pF)
U1_Ar_Helm = u1(Q2_grid,form_factor_Helm)
U1_Ar_KNad = u1(Q2_grid,form_factor_KNad)

with open('../csv/form_factors/u1_grids/argon/3pF.csv','w',newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
    for q2, u1 in zip(Q2_grid,U1_Ar_3pF):
        q2 = '{:.17f}'.format(q2)
        u1 = '{:.30f}'.format(u1)
        writer.writerow([q2, u1])

with open('../csv/form_factors/u1_grids/argon/Helm.csv','w',newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
    for q2, u1 in zip(Q2_grid,U1_Ar_Helm):
        q2 = '{:.17f}'.format(q2)
        u1 = '{:.30f}'.format(u1)
        writer.writerow([q2, u1])

with open('../csv/form_factors/u1_grids/argon/KNad.csv','w',newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
    for q2, u1 in zip(Q2_grid,U1_Ar_KNad):
        q2 = '{:.17f}'.format(q2)
        u1 = '{:.30f}'.format(u1)
        writer.writerow([q2, u1])

##################
#### Plotting ####
##################

fig1, ax1 = plt.subplots(1, 1, tight_layout=True) # Argon form factor
fig2, ax2 = plt.subplots(1, 2, figsize=(26,13), tight_layout=True, sharex=True, sharey=True) # Proton and neutron form factors
fig3, ax3 = plt.subplots(1, 1, tight_layout=True) # Tungsten form factor
fig4, ax4 = plt.subplots(1, 1, tight_layout=True) # Iron form factor

### Plot 3-Fermi parameter form factor from CXX code
ax1.plot(Q_array, FF_Ar_Helm, '-', color = 'crimson', label='Helm', path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
ax1.plot(Q_array, FF_Ar_KNad, '-', color = 'orange', label='(ad.) KN', path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
ax1.plot(Q_array, FF_Ar_3pF, '-', color = 'royalblue', label='3-Fermi Param.', path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])

### Plot proton and neutron magnetic form factors
ax2[0].plot(Q_p_e, FF_p_e, '-', color = '#A71342', label='Electric', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax2[0].plot(Q_p_m, FF_p_m, '-', color = 'orange', label='Magnetic', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])

ax2[1].plot(Q_n_e, FF_n_e, '-', color = '#A71342', label='Electric', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax2[1].plot(Q_n_m, FF_n_m, '-', color = 'orange', label='Magnetic', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])

### Plot Woods-Saxon tungsten form factor 
ax3.plot(Q_array, FF_W_WS, '-', color = 'royalblue', label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])

### Plot CXX iron form factor 
ax4.plot(Q_Fe, FF_Fe, '-', color = 'royalblue', label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])

### Styling and save ###
# Scales #
ax1.set_yscale('log')
ax2[0].set_yscale('log')
ax2[1].set_yscale('log')
ax3.set_yscale('log')
ax4.set_yscale('log')

#Fixing ticks
locmaj1 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin1 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)
locmaj2 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin2 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)
locmaj3 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin3 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)
locmaj4 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin4 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)

ax1.yaxis.set_major_locator(locmaj1)
ax1.yaxis.set_minor_locator(locmin1)
ax2[0].yaxis.set_major_locator(locmaj2)
ax2[0].yaxis.set_minor_locator(locmin2)
ax3.yaxis.set_major_locator(locmaj3)
ax3.yaxis.set_minor_locator(locmin3)
ax4.yaxis.set_major_locator(locmaj4)
ax4.yaxis.set_minor_locator(locmin4)

ax1.set_ylim(1e-4,1)
ax1.set_xlim(0,0.7)

ax2[0].set_xlim(0,3)
ax2[1].set_xlim(0,3)

ax3.set_xlim(0,0.7)
ax3.set_ylim(1e-8,1)

ax4.set_xlim(0,0.7)
ax4.set_ylim(1e-4,1)

# Axis labels #
ax1.set_xlabel(r'\textbf{Momentum Transfer} $Q$ (GeV)')
ax1.set_ylabel(r'\textbf{Form Factor} $|F(Q^2)|$')
ax1.set_title(r'\textbf{Argon Nuclear Form Factor}')

ax2[0].set_xlabel(r'\textbf{Momentum Transfer} $Q$ (GeV)')
ax2[1].set_xlabel(r'\textbf{Momentum Transfer} $Q$ (GeV)')
ax2[0].set_ylabel(r'\textbf{Form Factor} $|F(Q^2)|$')
ax2[0].set_title(r'\textbf{Proton Form Factors}')
ax2[1].set_title(r'\textbf{Neutron Form Factors}')

ax2[0].legend(loc='upper right')
ax2[1].legend(loc='upper right')
ax1.legend()

ax3.set_xlabel(r'\textbf{Momentum Transfer} $Q$ (GeV)')
ax3.set_ylabel(r'\textbf{Form Factor} $|F(Q^2)|$')
ax3.set_title(r'\textbf{Tungsten Nuclear Form Factor}')

ax4.set_xlabel(r'\textbf{Momentum Transfer} $Q$ (GeV)')
ax4.set_ylabel(r'\textbf{Form Factor} $|F(Q^2)|$')
ax4.set_title(r'\textbf{Iron Nuclear Form Factor}')

# Grids #
ax1.grid(which = 'major', axis = 'both')
ax2[0].grid(which= 'major', axis = 'both')
ax2[1].grid(which= 'major', axis = 'both')
ax3.grid(which = 'major', axis = 'both')
ax4.grid(which = 'major', axis = 'both')

# Save #
fig1.savefig("../plots/form_factor_argon.png", dpi=100, bbox_inches='tight')
fig2.savefig("../plots/form_factor_nucleons.png", dpi=100, bbox_inches='tight')
fig3.savefig("../plots/form_factor_tungsten.png", dpi=100, bbox_inches='tight')
fig4.savefig("../plots/form_factor_iron.png", dpi=100, bbox_inches='tight')
