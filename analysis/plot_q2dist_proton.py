import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import csv

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

Q2DIR_1TAU = '../csv/distributions/vmu_to_vtau_tau+_mu-/nucleon/proton'
Q2DIR_2TAU = '../csv/distributions/vmu_to_vmu_tau+_tau-/nucleon/proton'

vertical_flag = False
horizontal_flag = True

##########################
#### Q2 Distributions ####
##########################

###### Nucleon Proton ######
### 1tau ; proton ; 5 GeV ###
Q_1tau_nucleon_p_20GeV = []

### 1tau ; proton ; 20 GeV ###
Q_1tau_nucleon_p_50GeV = []

### 1tau ; proton ; 50 GeV ###
Q_1tau_nucleon_p_200GeV = []

### 1tau ; proton ; 200 GeV ###
Q_1tau_nucleon_p_1000GeV = []

### 2tau ; proton ; 5 GeV ###
Q_2tau_nucleon_p_20GeV = []

### 2tau ; proton ; 20 GeV ###
Q_2tau_nucleon_p_50GeV = []

### 2tau ; proton ; 50 GeV ###
Q_2tau_nucleon_p_200GeV = []

### 2tau ; proton ; 200 GeV ###
Q_2tau_nucleon_p_1000GeV = []


############################
#### Read Distributions ####
############################

###### Coherent Argon ######
#### vmu -> vtau tau+ mu- ####
### Read 1tau ; nucleon ; proton ; 10 Gev ###
with open(Q2DIR_1TAU+'/1tau_q2dist_nucleon_p_20GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_1tau_nucleon_p_20GeV.append(np.sqrt(q2))

### Read 1tau ; nucleon ; proton ; 20 Gev ###
with open(Q2DIR_1TAU+'/1tau_q2dist_nucleon_p_50GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_1tau_nucleon_p_50GeV.append(np.sqrt(q2))

### Read 1tau ; nucleon ; proton ; 50 Gev ###
with open(Q2DIR_1TAU+'/1tau_q2dist_nucleon_p_200GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_1tau_nucleon_p_200GeV.append(np.sqrt(q2))

### Read 1tau ; nucleon ; proton ; 200 Gev ###
with open(Q2DIR_1TAU+'/1tau_q2dist_nucleon_p_1000GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_1tau_nucleon_p_1000GeV.append(np.sqrt(q2))

#### vmu -> vmu tau+ tau- ####
### Read 2tau ; nucleon ; proton ; 10 Gev ###
with open(Q2DIR_2TAU+'/2tau_q2dist_nucleon_p_20GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_2tau_nucleon_p_20GeV.append(np.sqrt(q2))

### Read 2tau ; nucleon ; proton ; 20 Gev ###
with open(Q2DIR_2TAU+'/2tau_q2dist_nucleon_p_50GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_2tau_nucleon_p_50GeV.append(np.sqrt(q2))

### Read 2tau ; nucleon ; proton ; 50 Gev ###
with open(Q2DIR_2TAU+'/2tau_q2dist_nucleon_p_200GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_2tau_nucleon_p_200GeV.append(np.sqrt(q2))

### Read 2tau ; nucleon ; proton ; 200 Gev ###
with open(Q2DIR_2TAU+'/2tau_q2dist_nucleon_p_1000GeV.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        q2 = float(row[0])
        Q_2tau_nucleon_p_1000GeV.append(np.sqrt(q2))


##################
#### Binning ####
##################
NBINS = 100

## Argon ##
# vmu -> vtau tau+ mu- #
print("vmu -> vtau tau+ mu-")
Q_1tau_nucleon_p_20GeV_max = max(Q_1tau_nucleon_p_20GeV)
Q_1tau_nucleon_p_20GeV_min = min(Q_1tau_nucleon_p_20GeV)
Q_1tau_nucleon_p_20GeV_bins = np.geomspace(Q_1tau_nucleon_p_20GeV_min, Q_1tau_nucleon_p_20GeV_max, num=NBINS)
print("Q ; Ar - 20 GeV:")
print("\tmin: ", Q_1tau_nucleon_p_20GeV_min)
print("\tmax: ", Q_1tau_nucleon_p_20GeV_max)

Q_1tau_nucleon_p_50GeV_max = max(Q_1tau_nucleon_p_50GeV)
Q_1tau_nucleon_p_50GeV_min = min(Q_1tau_nucleon_p_50GeV)
Q_1tau_nucleon_p_50GeV_bins = np.geomspace(Q_1tau_nucleon_p_50GeV_min, Q_1tau_nucleon_p_50GeV_max, num=NBINS)
print("Q ; Ar - 50 GeV:")
print("\tmin: ", Q_1tau_nucleon_p_50GeV_min)
print("\tmax: ", Q_1tau_nucleon_p_50GeV_max)

Q_1tau_nucleon_p_200GeV_max = max(Q_1tau_nucleon_p_200GeV)
Q_1tau_nucleon_p_200GeV_min = min(Q_1tau_nucleon_p_200GeV)
Q_1tau_nucleon_p_200GeV_bins = np.geomspace(Q_1tau_nucleon_p_200GeV_min, Q_1tau_nucleon_p_200GeV_max, num=NBINS) 
print("Q ; Ar - 200 GeV:")
print("\tmin: ", Q_1tau_nucleon_p_200GeV_min)
print("\tmax: ", Q_1tau_nucleon_p_200GeV_max)

Q_1tau_nucleon_p_1000GeV_max = max(Q_1tau_nucleon_p_1000GeV)
Q_1tau_nucleon_p_1000GeV_min = min(Q_1tau_nucleon_p_1000GeV)
Q_1tau_nucleon_p_1000GeV_bins = np.geomspace(Q_1tau_nucleon_p_1000GeV_min, Q_1tau_nucleon_p_1000GeV_max, num=NBINS) 
print("Q ; Ar - 1000 GeV:")
print("\tmin: ", Q_1tau_nucleon_p_1000GeV_min)
print("\tmax: ", Q_1tau_nucleon_p_1000GeV_max)

# vmu -> vmu tau+ tau- #
print("vmu -> vmu tau+ tau-")
Q_2tau_nucleon_p_20GeV_max = max(Q_2tau_nucleon_p_20GeV)
Q_2tau_nucleon_p_20GeV_min = min(Q_2tau_nucleon_p_20GeV)
Q_2tau_nucleon_p_20GeV_bins = np.geomspace(Q_2tau_nucleon_p_20GeV_min, Q_2tau_nucleon_p_20GeV_max, num=NBINS)
print("Q ; Ar - 20 GeV:")
print("\tmin: ", Q_2tau_nucleon_p_20GeV_min)
print("\tmax: ", Q_2tau_nucleon_p_20GeV_max)

Q_2tau_nucleon_p_50GeV_max = max(Q_2tau_nucleon_p_50GeV)
Q_2tau_nucleon_p_50GeV_min = min(Q_2tau_nucleon_p_50GeV)
Q_2tau_nucleon_p_50GeV_bins = np.geomspace(Q_2tau_nucleon_p_50GeV_min, Q_2tau_nucleon_p_50GeV_max, num=NBINS)
print("Q ; Ar - 50 GeV:")
print("\tmin: ", Q_2tau_nucleon_p_50GeV_min)
print("\tmax: ", Q_2tau_nucleon_p_50GeV_max)

Q_2tau_nucleon_p_200GeV_max = max(Q_2tau_nucleon_p_200GeV)
Q_2tau_nucleon_p_200GeV_min = min(Q_2tau_nucleon_p_200GeV)
Q_2tau_nucleon_p_200GeV_bins = np.geomspace(Q_2tau_nucleon_p_200GeV_min, Q_2tau_nucleon_p_200GeV_max, num=NBINS) 
print("Q ; Ar - 200 GeV:")
print("\tmin: ", Q_2tau_nucleon_p_200GeV_min)
print("\tmax: ", Q_2tau_nucleon_p_200GeV_max)

Q_2tau_nucleon_p_1000GeV_max = max(Q_2tau_nucleon_p_1000GeV)
Q_2tau_nucleon_p_1000GeV_min = min(Q_2tau_nucleon_p_1000GeV)
Q_2tau_nucleon_p_1000GeV_bins = np.geomspace(Q_2tau_nucleon_p_1000GeV_min, Q_2tau_nucleon_p_1000GeV_max, num=NBINS) 
print("Q ; Ar - 1000 GeV:")
print("\tmin: ", Q_2tau_nucleon_p_1000GeV_min)
print("\tmax: ", Q_2tau_nucleon_p_1000GeV_max)

##################
#### Plotting ####
##################

# vertical #
if vertical_flag:
    fig, ax = plt.subplots(4, 2, figsize=(28, 36), tight_layout=True)

    ### Plot W for 5, 20, 50 and 200 GeV ###
    #ax[0,1].hist(Q_1tau_coh_W_5GeV, density=False, bins=Q_1tau_coh_W_5GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    #ax[0,1].hist(Q_1tau_coh_W_5GeV, density=False, bins=Q_1tau_coh_W_5GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    #ax[0,1].text(0.905,0.85,r'{\bf W - 5 GeV}',ha='right',fontsize=30, transform=ax[0,1].transAxes)
    #
    #ax[1,1].hist(Q_1tau_coh_W_50GeV, density=False, bins=Q_1tau_coh_W_50GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    #ax[1,1].hist(Q_1tau_coh_W_50GeV, density=False, bins=Q_1tau_coh_W_50GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    #ax[1,1].text(0.905,0.85,r'{\bf W - 50 GeV}',ha='right',fontsize=30, transform=ax[1,1].transAxes)
    #
    #ax[2,1].hist(Q_1tau_coh_W_200GeV, density=False, bins=Q_1tau_coh_W_200GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    #ax[2,1].hist(Q_1tau_coh_W_200GeV, density=False, bins=Q_1tau_coh_W_200GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    #ax[2,1].text(0.905,0.85,r'{\bf W - 200 GeV}',ha='right',fontsize=30, transform=ax[2,1].transAxes)
    #
    #ax[3,1].hist(Q_1tau_coh_W_1000GeV, density=False, bins=Q_1tau_coh_W_1000GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    #ax[3,1].hist(Q_1tau_coh_W_1000GeV, density=False, bins=Q_1tau_coh_W_1000GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    #ax[3,1].text(0.905,0.85,r'{\bf W - 1000 GeV}',ha='right',fontsize=30, transform=ax[3,1].transAxes)

# horizontal #
if horizontal_flag:
    fig, ax = plt.subplots(2, 4, figsize=(50, 26), tight_layout=True)

    ### vmu -> vtau tau+ mu- ###
    # Plot Ar for 10, 20, 50 and 200 GeV #
    ax[0,0].hist(Q_1tau_nucleon_p_20GeV, density=False, bins=Q_1tau_nucleon_p_20GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[0,0].hist(Q_1tau_nucleon_p_20GeV, density=False, bins=Q_1tau_nucleon_p_20GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[0,0].text(0.935,0.90,r'$E_\nu=$\textbf{ 20 GeV}',ha='right',fontsize=30, transform=ax[0,0].transAxes)

    ax[0,1].hist(Q_1tau_nucleon_p_50GeV, density=False, bins=Q_1tau_nucleon_p_50GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[0,1].hist(Q_1tau_nucleon_p_50GeV, density=False, bins=Q_1tau_nucleon_p_50GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[0,1].text(0.935,0.90,r'$E_\nu=$\textbf{ 50 GeV}',ha='right',fontsize=30, transform=ax[0,1].transAxes)

    ax[0,2].hist(Q_1tau_nucleon_p_200GeV, density=False, bins=Q_1tau_nucleon_p_200GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[0,2].hist(Q_1tau_nucleon_p_200GeV, density=False, bins=Q_1tau_nucleon_p_200GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[0,2].text(0.935,0.90,r'$E_\nu=$\textbf{ 200 GeV}',ha='right',fontsize=30, transform=ax[0,2].transAxes)

    ax[0,3].hist(Q_1tau_nucleon_p_1000GeV, density=False, bins=Q_1tau_nucleon_p_1000GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[0,3].hist(Q_1tau_nucleon_p_1000GeV, density=False, bins=Q_1tau_nucleon_p_1000GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[0,3].text(0.935,0.90,r'$E_\nu=$\textbf{ 1000 GeV}',ha='right',fontsize=30, transform=ax[0,3].transAxes)

    ### vmu -> vmu tau+ tau- ###
    # Plot Ar for 10, 20, 50 and 200 GeV #
    ax[1,0].hist(Q_2tau_nucleon_p_20GeV, density=False, bins=Q_2tau_nucleon_p_20GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[1,0].hist(Q_2tau_nucleon_p_20GeV, density=False, bins=Q_2tau_nucleon_p_20GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[1,0].text(0.935,0.90,r'$E_\nu=$\textbf{ 20 GeV}',ha='right',fontsize=30, transform=ax[1,0].transAxes)

    ax[1,1].hist(Q_2tau_nucleon_p_50GeV, density=False, bins=Q_2tau_nucleon_p_50GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[1,1].hist(Q_2tau_nucleon_p_50GeV, density=False, bins=Q_2tau_nucleon_p_50GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[1,1].text(0.935,0.90,r'$E_\nu=$\textbf{ 50 GeV}',ha='right',fontsize=30, transform=ax[1,1].transAxes)

    ax[1,2].hist(Q_2tau_nucleon_p_200GeV, density=False, bins=Q_2tau_nucleon_p_200GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[1,2].hist(Q_2tau_nucleon_p_200GeV, density=False, bins=Q_2tau_nucleon_p_200GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[1,2].text(0.935,0.90,r'$E_\nu=$\textbf{ 200 GeV}',ha='right',fontsize=30, transform=ax[1,2].transAxes)

    ax[1,3].hist(Q_2tau_nucleon_p_1000GeV, density=False, bins=Q_2tau_nucleon_p_1000GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[1,3].hist(Q_2tau_nucleon_p_1000GeV, density=False, bins=Q_2tau_nucleon_p_1000GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[1,3].text(0.935,0.90,r'$E_\nu=$\textbf{ 1000 GeV}',ha='right',fontsize=30, transform=ax[1,3].transAxes)

 
    ### Styling and save ###
    # Axis labels #
    ax[0,0].set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$', fontsize=60)

    ax[1,0].set_ylabel(r'$\nu_\mu \to \nu_\mu \tau^+ \tau^-$', fontsize=60)
    
    locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
    locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)

    xmin_values = [Q_1tau_nucleon_p_20GeV_min, Q_1tau_nucleon_p_50GeV_min, Q_1tau_nucleon_p_200GeV_min, Q_1tau_nucleon_p_1000GeV_min,
                   Q_2tau_nucleon_p_20GeV_min, Q_2tau_nucleon_p_50GeV_min, Q_2tau_nucleon_p_200GeV_min, Q_2tau_nucleon_p_1000GeV_min]

    k = 0
    for i in range(2):
        for j in range(4):        
            xmin = xmin_values[k]
            print(i, j, k, xmin)
            ax[i,j].xaxis.set_major_locator(locmaj)
            ax[i,j].xaxis.set_minor_locator(locmin)
            ax[i,j].xaxis.set_major_locator(locmaj)
            ax[i,j].xaxis.set_minor_locator(locmin)
            ax[i,j].set_ylim(1,1e7)
            ax[i,j].set_xlim(xmin,20)
            
            ax[i,j].set_xscale('log')
            ax[i,j].set_yscale('log')
            k += 1
#            ax[i,j].yaxis.set_major_locator(locmaj)
#            ax[i,j].yaxis.set_minor_locator(locmin)
#            ax[i,j].yaxis.set_major_locator(locmaj)
#            ax[i,j].yaxis.set_minor_locator(locmin)


fig.suptitle(r'\textbf{Distribution for Incoherent Scattering (proton) off $^{40}$Ar}', fontsize=50)
fig.supxlabel(r'\textbf{Momentum Transfer} $Q$ (GeV)', fontsize=60)
# Save #
fig.savefig("../plots/tau_q2dist_incoh_p.png", dpi=400, bbox_inches='tight')
