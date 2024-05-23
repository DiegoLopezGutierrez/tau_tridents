import numpy as np
from scipy.integrate import simpson
import matplotlib.pyplot as plt
import csv

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

##################################
#### Fourier-Bessel Expansion ####
##################################

### Data from CXX file ###
Q2s = []   # data from CXX file is expressed as [Q^2,
Ar_Us = []    # U] where U is int_Q^2^infty FF^2 dQ^2. See Eq. 24 in Lovseth and Radomski.

Qs = []
Us = []

### Read CXX file ###
with open('../csv/form_factors/argon.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        q = np.sqrt(q2)
        Q2s.append(q2)
        Qs.append(q)
        Ar_Us.append(float(row[1]))

### Calculate form factor from data and store ###
#FF2_array = np.gradient(Us, Q2s)

#for i in range(len(Q2s)):
#    q = np.sqrt(Q2s[i])
#    FF = np.sqrt(-FF2_array[i])
#    Q_FB_Alt.append(q)
#    Ar_FM_FB_Alt.append(FF)


###################################
##### Woods-Saxon Form Factor #####
###################################

### Form factor function from Woods-Saxon distribution ###
fm_to_invGeV = 5.068 # 1 fm = 5.068 GeV^-1
a = 0.523 * fm_to_invGeV # fm -> GeV^-1

def WS_form_factor(A, Q):
    if Q == 0:   # First value of Q2s is 0.0 which will lead to division by zero when integrating
        Q = 0.00001
    r0 = 1.126 * A**(1/3) * fm_to_invGeV   # Form factor should be dimensionless. Convert fm to GeV^-1.
    Q2 = Q*Q

    # WS form factor can be written analytically; see (A.3) in Ballett et al.
    WS = (3 * np.pi * a) / (r0**2 + np.pi**2 * a**2)

    WS *= np.pi * a * np.tanh(np.pi * Q * a)**(-1) * np.sin(Q * r0) - r0 * np.cos(Q * r0)

    WS /= Q * r0 * np.sinh(np.pi * Q * a)

    return WS

### Calculate form factors for argon and tungsten ###
W_A = 184

# Normalize to F(0) = 1; use 0.001 instead.
W_norm  = WS_form_factor(W_A, 0.001)

# Q is in GeV
W_FM = [np.abs( WS_form_factor(W_A, Q) / W_norm)**2 for Q in Qs]

def U1(q2_array, F2_array):
    return simpson(F2_array, x=q2_array) 

for i, q2 in enumerate(Q2s):
    q2_array = Q2s[i:]
    F2_array = W_FM[i:]
    u1 = U1(q2_array, F2_array)
    Us.append(u1)

print(len(Q2s), len(Us))

with open('tungsten.csv','w',newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for q2, u1 in zip(Q2s, Us):
        q2 = '{:.17f},'.format(q2)
        u1 = '{:.30f},'.format(u1)
        writer.writerow([q2, ",", u1, ","])

##################
#### Plotting ####
##################

fig, ax = plt.subplots(1, 1, tight_layout=True)

### Plot Woods-Saxon form factors for argon and tungsten ###
ax.scatter(Q2s, Ar_Us, color = 'c', s=25, marker='o', edgecolor='k', linewidths=1, alpha=0.75, label = r'$^{40}$Ar (digit)')
ax.scatter(Q2s, Us, color = 'firebrick', s=25, marker='^', edgecolor='k', linewidths=1, alpha=0.75, label = r'$^{184}$W')

### Styling and save ###
# Axis labels #
ax.set_xlabel('$Q^2$ (GeV$^2$)')
ax.set_ylabel(r'$u_1$')
ax.set_title('Woods-Saxon')

# Axis limits #
#ax[1].set_ylim(1e-7, 1)
#ax[1].set_xlim(0,1)

# Legends #
ax.legend(fontsize='25')

# Scales #
ax.set_yscale('log')

# Grids #
ax.grid(which = 'major', axis = 'both')

# Save #
fig.savefig("../plots/u1s.png", dpi=400, bbox_inches='tight')
