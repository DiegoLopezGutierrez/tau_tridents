import numpy as np
from scipy.integrate import simpson
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.patheffects as pe
import csv

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

##################################
#### Fourier-Bessel Expansion ####
##################################

### Q2 for argon ###
#Q_3Fp_Alt_digitized = []

### Argon form factors ###
#Ar_FM_3Fp_Alt_digitized = []

Mproton = 0.938272

### Data from CXX file ###
Q2s_Ar = []   # data from CXX file is expressed as [Q^2,
Us_Ar = []    # U] where U is int_Q^2^infty FF^2 dQ^2. See Eq. 24 in Lovseth and Radomski.

Q2s_p = []
Us_p = []

Q2s_n = []
Us_n = []

Q_Ar = []
FF_Ar = []

Q_p = []
FF_p = []

Q_n = []
FF_n = []

Q_p_e = []
Q_p_m = []
FF_p_e = []
FF_p_m = []
FF_p_calc = []

Q_n_e = []
Q_n_m = []
FF_n_e = []
FF_n_m = []
FF_n_calc = []

### Read CXX file ###
#with open('../csv/form_factors/argon.csv','r') as csvfile:
#    data = csv.reader(csvfile, delimiter = ',')
#    for row in data:
#        q2 = float(row[0])
#        u = float(row[1])
#        Q2s_Ar.append(q2)
#        Us_Ar.append(u)
#
#with open('../csv/form_factors/proton.csv','r') as csvfile:
#    data = csv.reader(csvfile, delimiter = ',')
#    for row in data:
#        q2 = float(row[0])
#        u = float(row[1])
#        Q2s_p.append(q2)
#        Us_p.append(u)
#
#with open('../csv/form_factors/neutron.csv','r') as csvfile:
#    data = csv.reader(csvfile, delimiter = ',')
#    for row in data:
#        q2 = float(row[0])
#        u = float(row[1])
#        Q2s_n.append(q2)
#        Us_n.append(u)
#
#### Calculate form factor from data and store ###
#FF2_array_Ar = np.gradient(Us_Ar, Q2s_Ar)
#FF2_array_p = np.gradient(Us_p, Q2s_p)
#FF2_array_n = np.gradient(Us_n, Q2s_n)
#
#for i in range(len(Q2s_Ar)):
#    q = np.sqrt(Q2s_Ar[i])
#    ff = np.sqrt(-FF2_array_Ar[i])
#    Q_Ar.append(q)
#    FF_Ar.append(ff)

#for i in range(len(Q2s_p)):
#    q = np.sqrt(Q2s_p[i])
#    ff = np.sqrt(-FF2_array_p[i])
#    Q_p.append(q)
#    FF_p.append(ff)
#
#for i in range(len(Q2s_n)):
#    q = np.sqrt(Q2s_n[i])
#    ff = np.sqrt(-FF2_array_n[i])
#    Q_n.append(q)
#    FF_n.append(ff)

### Read digitized data and store ###
#with open('../csv/form_factors/FF_Ar_FourierBessel_Alt.csv') as csvfile:
#    data = csv.reader(csvfile, delimiter = ',')
#    for row in data:
#        Q_FB_Alt_digitized.append(float(row[0]))
#        Ar_FM_FB_Alt_digitized.append(float(row[1]))
#
#with open('../csv/form_factors/FF_Ar_3Fp-red_Alt.csv') as csvfile:
#    data = csv.reader(csvfile, delimiter = ',')
#    for row in data:
#        Q_3Fp_Alt_digitized.append(float(row[0]))
#        Ar_FM_3Fp_Alt_digitized.append(float(row[1]))

### Normalize form factor to F(0) = 1 and save abs value ###
#Ar_FM_FB_Alt_norm = Ar_FM_FB_Alt[0]
#Ar_FM_FB_abs_Alt = [np.abs(fm / Ar_FM_FB_Alt_norm) for fm in Ar_FM_FB_Alt]

with open('../csv/form_factors/FF_Ar_3Fp-red_Alt.csv') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        Q_Ar.append(float(row[0]))
        FF_Ar.append(float(row[1]))

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

#for i, Q in enumerate(Q_p_e):
#    FF_p_calc.append(np.sqrt(FF_p_e[i]**2 + Q**2/(4*Mproton**2)*FF_p_n[i]**2))
#
#for i, Q in enumerate(Q_p_e):
#    FF_p_calc.append(np.sqrt(FF_p_e[i]**2 + Q**2/(4*Mproton**2)*FF_p_n[i]**2))

##################
#### Plotting ####
##################

fig1, ax1 = plt.subplots(1, 1, tight_layout=True)
fig2, ax2 = plt.subplots(1, 2, figsize=(26,13), tight_layout=True, sharex=True, sharey=True)

### Plot 3-Fermi parameter form factor from CXX code
ax1.plot(Q_Ar, FF_Ar, '-', color = 'royalblue', label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])

### Plot proton and neutron magnetic form factors
ax2[0].plot(Q_p_e, FF_p_e, '-', color = 'green', label='Electric', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax2[0].plot(Q_p_m, FF_p_m, '-', color = 'orange', label='Magnetic', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
#ax2[0].plot(Q_p, FF_p, '--', color = 'royalblue', label='CXX', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])

ax2[1].plot(Q_n_e, FF_n_e, '-', color = 'green', label='Electric', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax2[1].plot(Q_n_m, FF_n_m, '-', color = 'orange', label='Magnetic', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
#ax2[1].plot(Q_n, FF_n, '--', color = 'royalblue', label='CXX', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])

### Styling and save ###
# Scales #
ax1.set_yscale('log')
ax2[0].set_yscale('log')
ax2[1].set_yscale('log')

#Fixing ticks
locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)
ax1.yaxis.set_major_locator(locmaj)
ax1.yaxis.set_minor_locator(locmin)
ax2[0].yaxis.set_major_locator(locmaj)
ax2[0].yaxis.set_minor_locator(locmin)

ax1.set_ylim(1e-4,1)
ax1.set_xlim(0,0.7)

# Axis labels #
ax1.set_xlabel(r'\textbf{Momentum Transfer} $Q$ (GeV)')
ax1.set_ylabel(r'$|F(Q^2)|$')
ax1.set_title(r'\textbf{Argon Nuclear Form Factor}')

ax2[0].set_xlabel(r'\textbf{Momentum Transfer} $Q$ (GeV)')
ax2[1].set_xlabel(r'\textbf{Momentum Transfer} $Q$ (GeV)')
ax2[0].set_ylabel(r'$|F(Q^2)|$')
ax2[0].set_title(r'\textbf{Proton Form Factors}')
ax2[1].set_title(r'\textbf{Neutron Form Factors}')

ax2[0].legend(loc='upper right')
ax2[1].legend(loc='upper right')

# Grids #
ax1.grid(which = 'major', axis = 'both')
ax2[0].grid(which= 'major', axis = 'both')
ax2[1].grid(which= 'major', axis = 'both')

# Save #
fig1.savefig("../plots/Ar_form_factor.png", dpi=400, bbox_inches='tight')
fig2.savefig("../plots/nucleon_form_factors.png", dpi=400, bbox_inches='tight')
