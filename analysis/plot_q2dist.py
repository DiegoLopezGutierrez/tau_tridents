import numpy as np
import matplotlib.pyplot as plt
import csv

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

Q2DIR_1TAU = '../csv/distributions/vmu_to_vtau_tau+_mu-'

##########################
#### Q2 Distributions ####
##########################

###### Coherent Argon ######
### 1tau ; coherent ; argon ; 5 GeV ###
Q_1tau_coh_Ar_5GeV = []

### 1tau ; coherent ; argon ; 20 GeV ###
Q_1tau_coh_Ar_20GeV = []

### 1tau ; coherent ; argon ; 50 GeV ###
Q_1tau_coh_Ar_50GeV = []

### 1tau ; coherent ; argon ; 200 GeV ###
Q_1tau_coh_Ar_200GeV = []

###### Nucleons ######
### 1tau ; nucleon ; proton ; 5 GeV ###
Q_1tau_nucleon_p_5GeV = []

### 1tau ; nucleon ; proton ; 20 GeV ###
Q_1tau_nucleon_p_20GeV = []

### 1tau ; nucleon ; proton ; 50 GeV ###
Q_1tau_nucleon_p_50GeV = []

### 1tau ; nucleon ; proton ; 200 GeV ###
Q_1tau_nucleon_p_200GeV = []

############################
#### Read Distributions ####
############################

###### Coherent Argon ######
### Read 1tau ; coherent ; argon ; 5 Gev ###
with open(Q2DIR_1TAU+'/coherent/argon/1tau_q2dist_coh_Ar_5GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_1tau_coh_Ar_5GeV.append(np.sqrt(q2))

### Read 1tau ; coherent ; argon ; 20 Gev ###
with open(Q2DIR_1TAU+'/coherent/argon/1tau_q2dist_coh_Ar_20GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_1tau_coh_Ar_20GeV.append(np.sqrt(q2))

### Read 1tau ; coherent ; argon ; 50 Gev ###
with open(Q2DIR_1TAU+'/coherent/argon/1tau_q2dist_coh_Ar_50GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_1tau_coh_Ar_50GeV.append(np.sqrt(q2))

### Read 1tau ; coherent ; argon ; 200 Gev ###
with open(Q2DIR_1TAU+'/coherent/argon/1tau_q2dist_coh_Ar_200GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_1tau_coh_Ar_200GeV.append(np.sqrt(q2))


###### Nucleons ######
### Read 1tau ; nucleon ; proton ; 5 Gev ###
with open(Q2DIR_1TAU+'/nucleon/proton/1tau_q2dist_nucleon_p_5GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_1tau_nucleon_p_5GeV.append(np.sqrt(q2))

### Read 1tau ; nucleon ; proton ; 20 Gev ###
with open(Q2DIR_1TAU+'/nucleon/proton/1tau_q2dist_nucleon_p_20GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_1tau_nucleon_p_20GeV.append(np.sqrt(q2))

### Read 1tau ; nucleon ; proton ; 50 Gev ###
with open(Q2DIR_1TAU+'/nucleon/proton/1tau_q2dist_nucleon_p_50GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_1tau_nucleon_p_50GeV.append(np.sqrt(q2))

### Read 1tau ; nucleon ; proton ; 200 Gev ###
with open(Q2DIR_1TAU+'/nucleon/proton/1tau_q2dist_nucleon_p_200GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_1tau_nucleon_p_200GeV.append(np.sqrt(q2))

##################
#### Binning ####
##################
NBINS = 100

Q_1tau_coh_Ar_5GeV_max = max(Q_1tau_coh_Ar_5GeV)
Q_1tau_coh_Ar_5GeV_min = min(Q_1tau_coh_Ar_5GeV)
Q_1tau_coh_Ar_5GeV_bins = np.geomspace(Q_1tau_coh_Ar_5GeV_min, Q_1tau_coh_Ar_5GeV_max, num=NBINS)
print("Q ; Ar - 5 GeV:")
print("\tmin: ", Q_1tau_coh_Ar_5GeV_min)
print("\tmax: ", Q_1tau_coh_Ar_5GeV_max)

Q_1tau_coh_Ar_20GeV_max = max(Q_1tau_coh_Ar_20GeV)
Q_1tau_coh_Ar_20GeV_min = min(Q_1tau_coh_Ar_20GeV)
Q_1tau_coh_Ar_20GeV_bins = np.geomspace(Q_1tau_coh_Ar_20GeV_min, Q_1tau_coh_Ar_20GeV_max, num=NBINS)
print("Q ; Ar - 20 GeV:")
print("\tmin: ", Q_1tau_coh_Ar_20GeV_min)
print("\tmax: ", Q_1tau_coh_Ar_20GeV_max)

Q_1tau_coh_Ar_50GeV_max = max(Q_1tau_coh_Ar_50GeV)
Q_1tau_coh_Ar_50GeV_min = min(Q_1tau_coh_Ar_50GeV)
Q_1tau_coh_Ar_50GeV_bins = np.geomspace(Q_1tau_coh_Ar_50GeV_min, Q_1tau_coh_Ar_50GeV_max, num=NBINS) 
print("Q ; Ar - 50 GeV:")
print("\tmin: ", Q_1tau_coh_Ar_50GeV_min)
print("\tmax: ", Q_1tau_coh_Ar_50GeV_max)

Q_1tau_coh_Ar_200GeV_max = max(Q_1tau_coh_Ar_200GeV)
Q_1tau_coh_Ar_200GeV_min = min(Q_1tau_coh_Ar_200GeV)
Q_1tau_coh_Ar_200GeV_bins = np.geomspace(Q_1tau_coh_Ar_200GeV_min, Q_1tau_coh_Ar_200GeV_max, num=NBINS) 
print("Q ; Ar - 200 GeV:")
print("\tmin: ", Q_1tau_coh_Ar_200GeV_min)
print("\tmax: ", Q_1tau_coh_Ar_200GeV_max)

Q_1tau_nucleon_p_5GeV_max = max(Q_1tau_nucleon_p_5GeV)
Q_1tau_nucleon_p_5GeV_min = min(Q_1tau_nucleon_p_5GeV)
Q_1tau_nucleon_p_5GeV_bins = np.geomspace(Q_1tau_nucleon_p_5GeV_min, Q_1tau_nucleon_p_5GeV_max, num=NBINS) 
print("Q ; p - 5 GeV:")
print("\tmin: ", Q_1tau_nucleon_p_5GeV_min)
print("\tmax: ", Q_1tau_nucleon_p_5GeV_max)

Q_1tau_nucleon_p_20GeV_max = max(Q_1tau_nucleon_p_20GeV)
Q_1tau_nucleon_p_20GeV_min = min(Q_1tau_nucleon_p_20GeV)
Q_1tau_nucleon_p_20GeV_bins = np.geomspace(Q_1tau_nucleon_p_20GeV_min, Q_1tau_nucleon_p_20GeV_max, num=NBINS) 
print("Q ; p - 20 GeV:")
print("\tmin: ", Q_1tau_nucleon_p_20GeV_min)
print("\tmax: ", Q_1tau_nucleon_p_20GeV_max)

Q_1tau_nucleon_p_50GeV_max = max(Q_1tau_nucleon_p_50GeV)
Q_1tau_nucleon_p_50GeV_min = min(Q_1tau_nucleon_p_50GeV)
Q_1tau_nucleon_p_50GeV_bins = np.geomspace(Q_1tau_nucleon_p_50GeV_min, Q_1tau_nucleon_p_50GeV_max, num=NBINS)
print("Q ; p - 50 GeV:")
print("\tmin: ", Q_1tau_nucleon_p_50GeV_min)
print("\tmax: ", Q_1tau_nucleon_p_50GeV_max)

Q_1tau_nucleon_p_200GeV_max = max(Q_1tau_nucleon_p_200GeV)
Q_1tau_nucleon_p_200GeV_min = min(Q_1tau_nucleon_p_200GeV)
Q_1tau_nucleon_p_200GeV_bins = np.geomspace(Q_1tau_nucleon_p_200GeV_min, Q_1tau_nucleon_p_200GeV_max, num=NBINS)
print("Q ; p - 200 GeV:")
print("\tmin: ", Q_1tau_nucleon_p_200GeV_min)
print("\tmax: ", Q_1tau_nucleon_p_200GeV_max)

##################
#### Plotting ####
##################

fig, ax = plt.subplots(4, 2, figsize=(28, 36), tight_layout=True)

### Plot Ar for 5, 20, 50 and 200 GeV ###
ax[0,0].hist(Q_1tau_coh_Ar_5GeV, density=False, bins=Q_1tau_coh_Ar_5GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
ax[0,0].hist(Q_1tau_coh_Ar_5GeV, density=False, bins=Q_1tau_coh_Ar_5GeV_bins, histtype='step', color='black',lw=2, alpha=1)
ax[0,0].text(0.905,0.85,r'{\bf Ar - 5 GeV}',ha='right',fontsize=30, transform=ax[0,0].transAxes)

ax[1,0].hist(Q_1tau_coh_Ar_20GeV, density=False, bins=Q_1tau_coh_Ar_20GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
ax[1,0].hist(Q_1tau_coh_Ar_20GeV, density=False, bins=Q_1tau_coh_Ar_20GeV_bins, histtype='step', color='black',lw=2, alpha=1)
ax[1,0].text(0.905,0.85,r'{\bf Ar - 20 GeV}',ha='right',fontsize=30, transform=ax[1,0].transAxes)

ax[2,0].hist(Q_1tau_coh_Ar_50GeV, density=False, bins=Q_1tau_coh_Ar_50GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
ax[2,0].hist(Q_1tau_coh_Ar_50GeV, density=False, bins=Q_1tau_coh_Ar_50GeV_bins, histtype='step', color='black',lw=2, alpha=1)
ax[2,0].text(0.905,0.85,r'{\bf Ar - 50 GeV}',ha='right',fontsize=30, transform=ax[2,0].transAxes)

ax[3,0].hist(Q_1tau_coh_Ar_200GeV, density=False, bins=Q_1tau_coh_Ar_200GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
ax[3,0].hist(Q_1tau_coh_Ar_200GeV, density=False, bins=Q_1tau_coh_Ar_200GeV_bins, histtype='step', color='black',lw=2, alpha=1)
ax[3,0].text(0.905,0.85,r'{\bf Ar - 200 GeV}',ha='right',fontsize=30, transform=ax[3,0].transAxes)

### Plot proton for 5, 20, 50 and 200 GeV ###
ax[0,1].hist(Q_1tau_nucleon_p_5GeV, density=False, bins=Q_1tau_nucleon_p_5GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
ax[0,1].hist(Q_1tau_nucleon_p_5GeV, density=False, bins=Q_1tau_nucleon_p_5GeV_bins, histtype='step', color='black',lw=2, alpha=1)
ax[0,1].text(0.905,0.85,r'{\bf $p$ - 5 GeV}',ha='right',fontsize=30, transform=ax[0,1].transAxes)

ax[1,1].hist(Q_1tau_nucleon_p_20GeV, density=False, bins=Q_1tau_nucleon_p_20GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
ax[1,1].hist(Q_1tau_nucleon_p_20GeV, density=False, bins=Q_1tau_nucleon_p_20GeV_bins, histtype='step', color='black',lw=2, alpha=1)
ax[1,1].text(0.225,0.85,r'{\bf $p$ - 20 GeV}',ha='right',fontsize=30, transform=ax[1,1].transAxes)

ax[2,1].hist(Q_1tau_nucleon_p_50GeV, density=False, bins=Q_1tau_nucleon_p_50GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
ax[2,1].hist(Q_1tau_nucleon_p_50GeV, density=False, bins=Q_1tau_nucleon_p_50GeV_bins, histtype='step', color='black',lw=2, alpha=1)
ax[2,1].text(0.225,0.85,r'{\bf $p$ - 50 GeV}',ha='right',fontsize=30, transform=ax[2,1].transAxes)

ax[3,1].hist(Q_1tau_nucleon_p_200GeV, density=False, bins=Q_1tau_nucleon_p_200GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
ax[3,1].hist(Q_1tau_nucleon_p_200GeV, density=False, bins=Q_1tau_nucleon_p_200GeV_bins, histtype='step', color='black',lw=2, alpha=1)
ax[3,1].text(0.225,0.85,r'{\bf $p$ - 200 GeV}',ha='right',fontsize=30, transform=ax[3,1].transAxes)

### Styling and save ###
# Axis labels #
ax[0,0].set_xlabel('$Q$ (GeV)')
#ax[0,0].set_ylabel('Counts')

ax[1,0].set_xlabel('$Q$ (GeV)')
#ax[1,0].set_ylabel('Counts')

ax[2,0].set_xlabel('$Q$ (GeV)')
#ax[2,0].set_ylabel('Counts')

ax[3,0].set_xlabel('$Q$ (GeV)')
#ax[3,0].set_ylabel('Counts')

ax[0,1].set_xlabel('$Q$ (GeV)')
#ax[0,1].set_ylabel('Counts')

ax[1,1].set_xlabel('$Q$ (GeV)')
#ax[1,1].set_ylabel('Counts')

ax[2,1].set_xlabel('$Q$ (GeV)')
#ax[2,1].set_ylabel('Counts')

ax[3,1].set_xlabel('$Q$ (GeV)')
#ax[3,1].set_ylabel('Counts')

# Scales #
ax[0,0].set_xscale('log')
ax[0,0].set_yscale('log')

ax[1,0].set_xscale('log')
ax[1,0].set_yscale('log')

ax[2,0].set_xscale('log')
ax[2,0].set_yscale('log')

ax[3,0].set_xscale('log')
ax[3,0].set_yscale('log')

ax[0,1].set_xscale('log')
ax[0,1].set_yscale('log')

ax[1,1].set_xscale('log')
ax[1,1].set_yscale('log')

ax[2,1].set_xscale('log')
ax[2,1].set_yscale('log')

ax[3,1].set_xscale('log')
ax[3,1].set_yscale('log')

fig.suptitle(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$', fontsize=50)
# Save #
fig.savefig("../plots/1tau_q2dist.png", dpi=400, bbox_inches='tight')
