import numpy as np
from scipy.integrate import simpson
import matplotlib.pyplot as plt
import csv

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

##################################
#### Fourier-Bessel Expansion ####
##################################

### Q2 for argon ###
Q_FB_Alt = []
Q_FB_Alt_digitized = []
Q_3Fp_Alt_digitized = []

### Argon form factors ###
Ar_FM_FB_Alt = [] # copied from CXX file -- Fourier Bessel expansion
Ar_FM_FB_Alt_digitized = [] # digitized from paper; should be the same as above.
Ar_FM_3Fp_Alt_digitized = []

### Data from CXX file ###
Q2s = []   # data from CXX file is expressed as [Q^2,
Us = []    # U] where U is int_Q^2^infty FF^2 dQ^2. See Eq. 24 in Lovseth and Radomski.



### Read CXX file ###
with open('../csv/form_factors/argon.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        u = float(row[1])
        Q2s.append(q2)
        Us.append(u)

### Calculate form factor from data and store ###
FF2_array = np.gradient(Us, Q2s)

for i in range(len(Q2s)):
    q = np.sqrt(Q2s[i])
    FF = np.sqrt(-FF2_array[i])
    Q_FB_Alt.append(q)
    Ar_FM_FB_Alt.append(FF)

### Read digitized data and store ###
with open('../csv/form_factors/FF_Ar_FourierBessel_Alt.csv') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        Q_FB_Alt_digitized.append(float(row[0]))
        Ar_FM_FB_Alt_digitized.append(float(row[1]))

with open('../csv/form_factors/FF_Ar_3Fp-red_Alt.csv') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        Q_3Fp_Alt_digitized.append(float(row[0]))
        Ar_FM_3Fp_Alt_digitized.append(float(row[1]))

### Normalize form factor to F(0) = 1 and save abs value ###
Ar_FM_FB_Alt_norm = Ar_FM_FB_Alt[0]
Ar_FM_FB_abs_Alt = [np.abs(fm / Ar_FM_FB_Alt_norm) for fm in Ar_FM_FB_Alt]

###################################
##### Woods-Saxon Form Factor #####
###################################

Ar_FM_WS_Alt_digitized = [] # dashed gray line in FF plot (Wood-Saxon charge distribution), digitized from Alt. paper

Q_WS_Alt_digitized = []

### Read digitized data and store ###
with open('../csv/form_factors/FF_Ar_2Fp-dashed_Alt.csv') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        Q_WS_Alt_digitized.append(float(row[0]))
        Ar_FM_WS_Alt_digitized.append(float(row[1]))

### Form factor function from Woods-Saxon distribution ###
fm_to_invGeV = 5.068 # 1 fm = 5.068 GeV^-1
a = 0.523 * fm_to_invGeV # fm -> GeV^-1

def WS_form_factor(A, Q):
    r0 = 1.126 * A**(1/3) * fm_to_invGeV   # Form factor should be dimensionless. Convert fm to GeV^-1.
    Q2 = Q*Q

    # WS form factor can be written analytically; see (A.3) in Ballett et al.
    WS = (3 * np.pi * a) / (r0**2 + np.pi**2 * a**2)

    WS *= np.pi * a * np.tanh(np.pi * Q * a)**(-1) * np.sin(Q * r0) - r0 * np.cos(Q * r0)

    WS /= Q * r0 * np.sinh(np.pi * Q * a)

    return WS

### Calculate form factors for argon and tungsten ###
Ar_A = 40
W_A = 184

# Normalize to F(0) = 1; use 0.001 instead.
Ar_norm = WS_form_factor(Ar_A, 0.001)
W_norm  = WS_form_factor(W_A, 0.001)

# Q is in GeV
Q_array = np.linspace(0.001, 3, 1000)

Ar_FM = [np.abs( WS_form_factor(Ar_A, Q) / Ar_norm ) for Q in Q_array]
W_FM = [np.abs( WS_form_factor(W_A, Q) / W_norm ) for Q in Q_array]

###################################
## 3-Fermi Parameter Form Factor ##
###################################

### Analytic 3Fp Argon ###
parameter_c = 3.73*fm_to_invGeV # fm
parameter_z = 0.62*fm_to_invGeV # fm
parameter_w = -0.19 # 1

def ThreeFp_form_factor(Q):
    r_array = np.linspace(0.00001*fm_to_invGeV,1000*fm_to_invGeV,num=1000) # in fm
    integrand_array = []
    rho0_integrand_array = []
    for r in r_array:
        rho = (1+parameter_w*r**2/parameter_c**2) / (1+np.exp((r-parameter_c)/parameter_z))
        integrand = r**2*np.sin(Q*r)/(Q*r)*rho
        rho0_integrand = r**2*rho
        integrand_array.append(integrand)
        rho0_integrand_array.append(rho0_integrand)

    integral = simpson(integrand_array, x=r_array)
    rho0_integral = simpson(rho0_integrand_array, x=r_array)
    return integral / rho0_integral

Ar_FM_3Fp = [np.abs(ThreeFp_form_factor(Q)) for Q in Q_array]

##################
#### Plotting ####
##################

fig, ax = plt.subplots(2, 1, figsize=(15, 12), tight_layout=True)

### Plot Woods-Saxon form factors for argon and tungsten ###
ax[1].plot(Q_array, Ar_FM, color = 'c', label = r'$^{40}$Ar')
ax[1].plot(Q_array, W_FM, color = 'firebrick', label = r'$^{184}$W')

### Plot Fourier-Bessel form factor from CXX code
ax[0].scatter(Q_FB_Alt, Ar_FM_FB_abs_Alt, color = 'royalblue', s=45, edgecolor='k', linewidths=1, alpha=0.75, marker='o', label = r'$^{40}$Ar (CXX)')

### Plot digitized Fourier-Bessel and Woods-Saxon form factors from Alt. paper ###
ax[0].scatter(Q_FB_Alt_digitized, Ar_FM_FB_Alt_digitized, color = 'orange', s=45, marker='^', edgecolor='k', linewidths=1, alpha=0.75, label= r'$^{40}$Ar (digit-FB)')
ax[0].scatter(Q_3Fp_Alt_digitized, Ar_FM_3Fp_Alt_digitized, color = 'green', s=45, marker='*', edgecolor='k', linewidth=1, alpha=0.75, label= r'$^{40}$Ar (digit-3Fp)')
ax[1].scatter(Q_WS_Alt_digitized, Ar_FM_WS_Alt_digitized, color = 'grey', s=20, marker='o', edgecolor='k', linewidths=1, alpha=0.75, label = r'$^{40}$Ar (digit)')

### Plot analytic 3-Fermi parameter form factor for argon 
ax[0].plot(Q_array, Ar_FM_3Fp, color = 'darkmagenta', label = r'$^{40}$Ar (3Fp)')

### Styling and save ###
# Axis labels #
ax[0].set_xlabel('$Q$ (GeV)')
ax[0].set_ylabel(r'$|F(Q^2)|$')
ax[0].set_title('Fourier-Bessel')

ax[1].set_xlabel('$Q$ (GeV)')
ax[1].set_ylabel(r'$|F(Q^2)|$')
ax[1].set_title('Woods-Saxon')

# Axis limits #
ax[0].set_ylim(1e-4, 1)
ax[0].set_xlim(0,0.75)

ax[1].set_ylim(1e-7, 1)
ax[1].set_xlim(0,1)

# Legends #
ax[0].legend(fontsize='25')
ax[1].legend(fontsize='25')

# Scales #
ax[0].set_yscale('log')
ax[1].set_yscale('log')

# Grids #
ax[0].grid(which = 'major', axis = 'both')
ax[1].grid(which = 'major', axis = 'both')

# Save #
fig.savefig("../plots/form_factors.png", dpi=400, bbox_inches='tight')
