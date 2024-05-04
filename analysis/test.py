import numpy as np
import matplotlib.pyplot as plt
import csv
import math

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

Q2DIR_1TAU = '../csv/distributions/vmu_to_vtau_tau+_mu-'

# Changed to proton instead but will leave names for now
Q_1tau_coh_Ar_5GeV = []

with open(Q2DIR_1TAU+'/nucleon/proton/1tau_q2dist_nucleon_p_5GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_1tau_coh_Ar_5GeV.append(np.sqrt(q2))

Q_1tau_coh_Ar_5GeV_max = max(Q_1tau_coh_Ar_5GeV)
Q_1tau_coh_Ar_5GeV_min = min(Q_1tau_coh_Ar_5GeV)

print("Min Q: ", Q_1tau_coh_Ar_5GeV_min)
print("Max Q: ", Q_1tau_coh_Ar_5GeV_max)

bin_max = math.ceil(np.log10(Q_1tau_coh_Ar_5GeV_max))
bin_min = math.floor(np.log10(Q_1tau_coh_Ar_5GeV_min))

print("Min logspace: ", bin_min)
print("Max logspace: ", bin_max)

bins = np.logspace(bin_min, bin_max, num=100)

binsgeom = np.geomspace(Q_1tau_coh_Ar_5GeV_min, Q_1tau_coh_Ar_5GeV_max, num=100)

fig, ax = plt.subplots(1, 2, figsize=(28, 12), tight_layout=True)

### Plot Ar for 5, 20, 50 and 200 GeV ###
ax[0].hist(Q_1tau_coh_Ar_5GeV, density=False, bins=bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
ax[0].hist(Q_1tau_coh_Ar_5GeV, density=False, bins=bins, histtype='step', color='black',lw=2, alpha=1)
ax[0].text(0.905,0.85,r'{\bf p - 5 GeV - log}',ha='right',fontsize=30, transform=ax[0].transAxes)

ax[0].set_xlabel('$Q$ (GeV)')
ax[0].set_xscale('log')
ax[0].set_yscale('log')

ax[1].hist(Q_1tau_coh_Ar_5GeV, density=False, bins=binsgeom, color='green',alpha=0.3, edgecolor='black', lw=0.5)
ax[1].hist(Q_1tau_coh_Ar_5GeV, density=False, bins=binsgeom, histtype='step', color='black',lw=2, alpha=1)
ax[1].text(0.905,0.85,r'{\bf p - 5 GeV - geom}',ha='right',fontsize=30, transform=ax[1].transAxes)

ax[1].set_xlabel('$Q$ (GeV)')
ax[1].set_xscale('log')
ax[1].set_yscale('log')

fig.suptitle(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$', fontsize=50)
# Save #
fig.savefig("test_q2.png", dpi=400, bbox_inches='tight')
