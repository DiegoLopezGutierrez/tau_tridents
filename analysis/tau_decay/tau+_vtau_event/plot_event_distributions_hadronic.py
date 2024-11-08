import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.patheffects as pe
import csv

STYLE_DIR = '../../../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

DIST_DIR_1TAU = '../../../csv/distributions/vmu_to_vtau_tau+_mu-'

RMiss_background_filename = '../../../csv/events/AlexSousa_RMiss_Background.csv'
energy_RMiss = []
background_RMiss = []

RMiss_NuWro_background_filename = DIST_DIR_1TAU+'/RT-thresholds-no-single.txt'
background_RMiss_NuWro = []

write_dist = False

############################
#### Read Distributions ####
############################

def read_four_vectors_incoh(file_path, nucleon):
    print("Reading file: ", file_path)
    if nucleon == 'proton':
        nucleon_pid = '2212'
    if nucleon == 'neutron':
        nucleon_pid = '2112'
    with open(file_path, 'r') as txtfile:
        nu_tau_array = []
        muon_array = []
        pVis_array = []
        nucleon_array = []
        pMiss_array = []
        i = 1
        for line in txtfile:
            data = line.split(',')
            if data[0] == '16':  # Read and store nu_tau
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                E  = float(data[5])
                nu_tau = np.array([px, py, pz, E])
                nu_tau_array.append(nu_tau)
            if data[0] == '13':  # Read and store muon
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                E  = float(data[5])
                muon = np.array([px, py, pz, E])
                muon_array.append(muon)
            if data[0] == '999':  # Read and store pMiss
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                E  = float(data[5])
                pVis = np.array([px, py, pz, E])
                assert E <= 33, "Unphysical energy in line {}: {}".format(i, E)
                pVis_array.append(pVis)
            if data[0] == nucleon_pid:  # Read and store antimuon
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                E  = float(data[5])
                nuc = np.array([px, py, pz, E])
                nucleon_array.append(nuc)

        if nucleon == 'proton':
            pMiss_array = [pVis + nucleon + muon for pVis, nucleon, muon in zip(pVis_array, nucleon_array, muon_array)]
        if nucleon == 'neutron':
            pMiss_array = [pVis + muon for pVis, muon in zip(pVis_array, muon_array)]
        return nu_tau_array, muon_array, pMiss_array, nucleon_array

def read_four_vectors_coh(file_path):
    with open(file_path, 'r') as txtfile:
        nu_tau_array = []
        muon_array = []
        nubar_tau_array = []
        for line in txtfile:
            data = line.split(',')
            if data[0] == '16':  # Read and store nu_tau
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                E  = float(data[5])
                nu_tau = np.array([px, py, pz, E])
                nu_tau_array.append(nu_tau)
            if data[0] == '13':  # Read and store muon
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                E  = float(data[5])
                muon = np.array([px, py, pz, E])
                muon_array.append(muon)
            if data[0] == '-16':  # Read and store nu_tau_bar
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                E  = float(data[5])
                nubar_tau = np.array([px, py, pz, E])
                nubar_tau_array.append(nubar_tau)
        
        pMiss_array = [nu_tau + nubar_tau for nu_tau, nubar_tau in zip(nu_tau_array, nubar_tau_array)]
        return nu_tau_array, muon_array, nubar_tau_array, pMiss_array

def calculate_ptMiss(pMiss_array, max_val):
    ptMiss_array = []
    for i,pMiss in enumerate(pMiss_array):
        pxMiss = pMiss[0]
        pyMiss = pMiss[1]
        ptMiss = np.sqrt(pxMiss**2 + pyMiss**2)
        if ptMiss >= max_val:
            print("Unphysical ptMiss for Event {}: {}".format(i+1, ptMiss))
            print("\tpMiss = {}".format(pMiss))
        ptMiss_array.append(ptMiss)
    return ptMiss_array

def calculate_pzMiss(pMiss_array):
    pzMiss_array = []
    for pMiss in pMiss_array:
        pzMiss = pMiss[2]
        pzMiss_array.append(pzMiss)
    return pzMiss_array

def calculate_ptMuon(muon_array):
    ptMuon_array = []
    for muon in muon_array:
        px = muon[0]
        py = muon[1]
        ptMuon = np.sqrt(px**2 + py**2)
        ptMuon_array.append(ptMuon)
    return ptMuon_array

def calculate_RMiss(ptMiss_array, ptMuon_array):
    RMiss_array = []
    for ptMiss, ptMuon in zip(ptMiss_array, ptMuon_array):
        RMiss = ptMiss / (ptMiss + ptMuon)
        RMiss_array.append(RMiss)
    return RMiss_array

def calculate_ThetaMiss(ptMiss_array, pzMiss_array):
    ThetaMiss_array = []
    for ptMiss, pzMiss in zip(ptMiss_array, pzMiss_array):
        ThetaMiss = math.acos(pzMiss / np.sqrt(ptMiss**2 + pzMiss**2))
        ThetaMiss_array.append(math.degrees(ThetaMiss))
    return ThetaMiss_array

############################
###### Distributions #######
############################

nu_tau_p_33GeV, muon_p_33GeV, pMiss_p_33GeV, nucleon_p_33GeV = read_four_vectors_incoh(DIST_DIR_1TAU+'/nucleon/proton/tau+_vtau_events/tau_hadronic_decayed_distribution_p_33GeV.txt', 'proton')
nu_tau_n_33GeV, muon_n_33GeV, pMiss_n_33GeV, nucleon_n_33GeV = read_four_vectors_incoh(DIST_DIR_1TAU+'/nucleon/neutron/tau+_vtau_events/tau_hadronic_decayed_distribution_n_33GeV.txt', 'neutron')
nu_tau_Ar_47GeV, muon_Ar_47GeV, nubar_tau_Ar_47GeV, pMiss_Ar_47GeV = read_four_vectors_coh(DIST_DIR_1TAU+'/coherent/argon/tau+_vtau_events/tau_hadronic_decayed_distribution_Ar_47GeV.txt')

print(len(pMiss_p_33GeV))
print(len(pMiss_n_33GeV))
print(len(pMiss_Ar_47GeV))

ptMiss_p_33GeV = calculate_ptMiss(pMiss_p_33GeV, 33)
ptMiss_n_33GeV = calculate_ptMiss(pMiss_n_33GeV, 33)
ptMiss_Ar_47GeV = calculate_ptMiss(pMiss_Ar_47GeV, 47)

pzMiss_p_33GeV = calculate_pzMiss(pMiss_p_33GeV)
pzMiss_n_33GeV = calculate_pzMiss(pMiss_n_33GeV)
pzMiss_Ar_47GeV = calculate_pzMiss(pMiss_Ar_47GeV)

ptMuon_p_33GeV = calculate_ptMuon(muon_p_33GeV)
ptMuon_n_33GeV = calculate_ptMuon(muon_n_33GeV)
ptMuon_Ar_47GeV = calculate_ptMuon(muon_Ar_47GeV)

RMiss_p_33GeV = calculate_RMiss(ptMiss_p_33GeV, ptMuon_p_33GeV)
RMiss_n_33GeV = calculate_RMiss(ptMiss_n_33GeV, ptMuon_n_33GeV)
RMiss_incoh_33GeV = RMiss_p_33GeV + RMiss_n_33GeV
RMiss_Ar_47GeV = calculate_RMiss(ptMiss_Ar_47GeV, ptMuon_Ar_47GeV)

ThetaMiss_p_33GeV = calculate_ThetaMiss(ptMiss_p_33GeV, pzMiss_p_33GeV)
ThetaMiss_n_33GeV = calculate_ThetaMiss(ptMiss_n_33GeV, pzMiss_n_33GeV)
ThetaMiss_Ar_47GeV = calculate_ThetaMiss(ptMiss_Ar_47GeV, pzMiss_Ar_47GeV)

if write_dist:
    with open(DIST_DIR_1TAU+'/nucleon/proton/ptDist_1tau_nucleon_p_33GeV.txt','w',newline='') as outfile:
        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for pt_miss, pt_muon, r_miss in zip(ptMiss_p_33GeV, ptMuon_p_33GeV, RMiss_p_33GeV):
            writer.writerow([pt_miss, pt_muon, r_miss])
    with open(DIST_DIR_1TAU+'/nucleon/neutron/ptDist_1tau_nucleon_n_33GeV.txt','w',newline='') as outfile:
        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for pt_miss, pt_muon, r_miss in zip(ptMiss_n_33GeV, ptMuon_n_33GeV, RMiss_n_33GeV):
            writer.writerow([pt_miss, pt_muon, r_miss])
    with open(DIST_DIR_1TAU+'/coherent/argon/ptDist_1tau_Ar_47GeV.txt','w',newline='') as outfile:
        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for pt_miss, pt_muon, r_miss in zip(ptMiss_Ar_47GeV, ptMuon_Ar_47GeV, RMiss_Ar_47GeV):
            writer.writerow([pt_miss, pt_muon, r_miss])

with open(RMiss_background_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        energy_RMiss.append(float(row[0]))
        background_RMiss.append(float(row[1]))

with open(RMiss_NuWro_background_filename,'r') as txtfile:
    for line in txtfile:
        background_RMiss_NuWro.append(float(line))

##################
#### Binning ####
##################
NBINS = 40

# vmu -> vtau tau+ mu- #
print("vmu -> vtau tau+ mu-")
print("---------- proton - 33 GeV ----------")
ptMiss_p_33GeV_max = max(ptMiss_p_33GeV)
ptMiss_p_33GeV_min = min(ptMiss_p_33GeV)
ptMiss_p_33GeV_bins, ptMiss_p_33GeV_step = np.linspace(ptMiss_p_33GeV_min, ptMiss_p_33GeV_max, num=NBINS, retstep=True)
print("ptMiss:")
print("\tmin: ", ptMiss_p_33GeV_min)
print("\tmax: ", ptMiss_p_33GeV_max)
ptMuon_p_33GeV_max = max(ptMuon_p_33GeV)
ptMuon_p_33GeV_min = min(ptMuon_p_33GeV)
ptMuon_p_33GeV_bins, ptMuon_p_33GeV_step = np.linspace(ptMuon_p_33GeV_min, ptMuon_p_33GeV_max, num=NBINS, retstep=True)
print("ptMuon:")
print("\tmin: ", ptMuon_p_33GeV_min)
print("\tmax: ", ptMuon_p_33GeV_max)
RMiss_p_33GeV_max = max(RMiss_p_33GeV)
RMiss_p_33GeV_min = min(RMiss_p_33GeV)
RMiss_p_33GeV_bins, RMiss_p_33GeV_step = np.linspace(0, 1, num=NBINS, retstep=True)
print("bin width: ", RMiss_p_33GeV_step)
print("RMiss:")
print("\tmin: ", RMiss_p_33GeV_min)
print("\tmax: ", RMiss_p_33GeV_max)
ThetaMiss_p_33GeV_max = max(ThetaMiss_p_33GeV)
ThetaMiss_p_33GeV_min = min(ThetaMiss_p_33GeV)
ThetaMiss_p_33GeV_bins, ThetaMiss_p_33GeV_step = np.linspace(0, 1, num=NBINS, retstep=True)
print("bin width: ", ThetaMiss_p_33GeV_step)
print("ThetaMiss:")
print("\tmin: ", ThetaMiss_p_33GeV_min)
print("\tmax: ", ThetaMiss_p_33GeV_max)

print("---------- neutron - 33 GeV ----------")
ptMiss_n_33GeV_max = max(ptMiss_n_33GeV)
ptMiss_n_33GeV_min = min(ptMiss_n_33GeV)
ptMiss_n_33GeV_bins, ptMiss_n_33GeV_step = np.linspace(ptMiss_n_33GeV_min, ptMiss_n_33GeV_max, num=NBINS, retstep=True)
print("ptMiss:")
print("\tmin: ", ptMiss_n_33GeV_min)
print("\tmax: ", ptMiss_n_33GeV_max)
ptMuon_n_33GeV_max = max(ptMuon_n_33GeV)
ptMuon_n_33GeV_min = min(ptMuon_n_33GeV)
ptMuon_n_33GeV_bins, ptMuon_n_33GeV_step = np.linspace(ptMuon_n_33GeV_min, ptMuon_n_33GeV_max, num=NBINS, retstep=True)
print("ptMuon:")
print("\tmin: ", ptMuon_n_33GeV_min)
print("\tmax: ", ptMuon_n_33GeV_max)
RMiss_n_33GeV_max = max(RMiss_n_33GeV)
RMiss_n_33GeV_min = min(RMiss_n_33GeV)
RMiss_n_33GeV_bins, RMiss_n_33GeV_step = np.linspace(0, 1, num=NBINS, retstep=True)
print("bin width: ", RMiss_n_33GeV_step)
print("RMiss:")
print("\tmin: ", RMiss_n_33GeV_min)
print("\tmax: ", RMiss_n_33GeV_max)
ThetaMiss_n_33GeV_max = max(ThetaMiss_n_33GeV)
ThetaMiss_n_33GeV_min = min(ThetaMiss_n_33GeV)
ThetaMiss_n_33GeV_bins, ThetaMiss_n_33GeV_step = np.linspace(0, 1, num=NBINS, retstep=True)
print("bin width: ", ThetaMiss_n_33GeV_step)
print("ThetaMiss:")
print("\tmin: ", ThetaMiss_n_33GeV_min)
print("\tmax: ", ThetaMiss_n_33GeV_max)

print("---------- argon - 47 GeV ----------")
ptMiss_Ar_47GeV_max = max(ptMiss_Ar_47GeV)
ptMiss_Ar_47GeV_min = min(ptMiss_Ar_47GeV)
ptMiss_Ar_47GeV_bins, ptMiss_Ar_47GeV_step = np.linspace(ptMiss_Ar_47GeV_min, ptMiss_Ar_47GeV_max, num=NBINS, retstep=True)
print("ptMiss:")
print("\tmin: ", ptMiss_Ar_47GeV_min)
print("\tmax: ", ptMiss_Ar_47GeV_max)
ptMuon_Ar_47GeV_max = max(ptMuon_Ar_47GeV)
ptMuon_Ar_47GeV_min = min(ptMuon_Ar_47GeV)
ptMuon_Ar_47GeV_bins, ptMuon_Ar_47GeV_step = np.linspace(ptMuon_Ar_47GeV_min, ptMuon_Ar_47GeV_max, num=NBINS, retstep=True)
print("ptMuon:")
print("\tmin: ", ptMuon_Ar_47GeV_min)
print("\tmax: ", ptMuon_Ar_47GeV_max)
RMiss_Ar_47GeV_max = max(RMiss_Ar_47GeV)
RMiss_Ar_47GeV_min = min(RMiss_Ar_47GeV)
RMiss_Ar_47GeV_bins, RMiss_Ar_47GeV_step = np.linspace(0, 1, num=NBINS, retstep=True)
print("bin width: ", RMiss_Ar_47GeV_step)
print("RMiss:")
print("\tmin: ", RMiss_Ar_47GeV_min)
print("\tmax: ", RMiss_Ar_47GeV_max)
ThetaMiss_Ar_47GeV_max = max(ThetaMiss_Ar_47GeV)
ThetaMiss_Ar_47GeV_min = min(ThetaMiss_Ar_47GeV)
ThetaMiss_Ar_47GeV_bins, ThetaMiss_Ar_47GeV_step = np.linspace(0, 1, num=NBINS, retstep=True)
print("bin width: ", ThetaMiss_Ar_47GeV_step)
print("ThetaMiss:")
print("\tmin: ", ThetaMiss_Ar_47GeV_min)
print("\tmax: ", ThetaMiss_Ar_47GeV_max)

### Weights ###
NEVENTS_p = len(ptMiss_p_33GeV)
NEVENTS_n = len(ptMiss_n_33GeV)
NEVENTS_Ar = len(ptMiss_Ar_47GeV)

print("Total events:")
print(f"\tproton: {NEVENTS_p}")
print(f"\tneutron: {NEVENTS_n}")
print(f"\targon: {NEVENTS_Ar}")

# proton - 33 GeV #
wts_ptMiss_p_33GeV = np.ones_like(ptMiss_p_33GeV) / (NEVENTS_p * ptMiss_p_33GeV_step)
wts_ptMuon_p_33GeV = np.ones_like(ptMuon_p_33GeV) / (NEVENTS_p * ptMuon_p_33GeV_step)
wts_RMiss_p_33GeV = np.ones_like(RMiss_p_33GeV) / (NEVENTS_p * RMiss_p_33GeV_step)
wts_ThetaMiss_p_33GeV = np.ones_like(ThetaMiss_p_33GeV) / (NEVENTS_p * ThetaMiss_p_33GeV_step)

# neutron - 33 GeV #
wts_ptMiss_n_33GeV = np.ones_like(ptMiss_n_33GeV) / (NEVENTS_n * ptMiss_n_33GeV_step)
wts_ptMuon_n_33GeV = np.ones_like(ptMuon_n_33GeV) / (NEVENTS_n * ptMuon_n_33GeV_step)
wts_RMiss_n_33GeV = np.ones_like(RMiss_n_33GeV) / (NEVENTS_n * RMiss_n_33GeV_step)
wts_ThetaMiss_n_33GeV = np.ones_like(ThetaMiss_n_33GeV) / (NEVENTS_n * ThetaMiss_n_33GeV_step)

# incoherent - 33 GeV #
wts_RMiss_incoh_33GeV = np.ones_like(RMiss_incoh_33GeV) / ((NEVENTS_p + NEVENTS_n) * RMiss_n_33GeV_step)

# argon - 47 GeV #
wts_ptMiss_Ar_47GeV = np.ones_like(ptMiss_Ar_47GeV) / (NEVENTS_Ar * ptMiss_Ar_47GeV_step)
wts_ptMuon_Ar_47GeV = np.ones_like(ptMuon_Ar_47GeV) / (NEVENTS_Ar * ptMuon_Ar_47GeV_step)
wts_RMiss_Ar_47GeV = np.ones_like(RMiss_Ar_47GeV) / (NEVENTS_Ar * RMiss_Ar_47GeV_step)
wts_ThetaMiss_Ar_47GeV = np.ones_like(ThetaMiss_Ar_47GeV) / (NEVENTS_Ar * ThetaMiss_Ar_47GeV_step)

# background - Pedro #
wts_RMiss_Pedro = np.ones_like(background_RMiss_NuWro) / (2922963 * RMiss_p_33GeV_step)

##################
#### Plotting ####
##################

fig6, (ax61, ax62) = plt.subplots(1, 2, sharey=True, figsize=(24,15), tight_layout=True, gridspec_kw={'wspace':0}) # RMiss incoh and coh

signal_color = '#2B4ACA'
background_color = '#E70C64'
background_color_NuWro = '#DE7D39'
opacity = 0.50

### RMiss ###
ax61.hist(RMiss_incoh_33GeV, density=False, bins=RMiss_p_33GeV_bins, weights=wts_RMiss_incoh_33GeV, color=signal_color, alpha=opacity, edgecolor='black', lw=0.5)
ax61.hist(RMiss_incoh_33GeV, density=False, bins=RMiss_p_33GeV_bins, weights=wts_RMiss_incoh_33GeV, histtype='step', color=signal_color, alpha=1, lw=2)
ax61.hist(background_RMiss_NuWro, density=False, bins=RMiss_p_33GeV_bins, weights=wts_RMiss_Pedro, color=background_color_NuWro,alpha=opacity, edgecolor='black',lw=0.5)
ax61.hist(background_RMiss_NuWro, density=False, bins=RMiss_p_33GeV_bins, weights=wts_RMiss_Pedro, histtype='step', color=background_color_NuWro, alpha=1, lw=2, hatch='o')
ax61.hist(energy_RMiss, density=False, bins=RMiss_p_33GeV_bins, weights=background_RMiss, color=background_color,alpha=opacity, edgecolor='black',lw=0.5)
ax61.hist(energy_RMiss, density=False, bins=RMiss_p_33GeV_bins, weights=background_RMiss, histtype='step',color=background_color,alpha=1,lw=2, hatch='/')
ax62.hist(RMiss_Ar_47GeV, density=False, bins=RMiss_Ar_47GeV_bins, weights=wts_RMiss_Ar_47GeV, color=signal_color, alpha=opacity, edgecolor='black', lw=0.5)
ax62.hist(RMiss_Ar_47GeV, density=False, bins=RMiss_Ar_47GeV_bins, weights=wts_RMiss_Ar_47GeV, histtype='step', color=signal_color, alpha=1, lw=2)
ax62.hist(background_RMiss_NuWro, density=False, bins=RMiss_p_33GeV_bins, weights=wts_RMiss_Pedro, color=background_color_NuWro,alpha=opacity, edgecolor='black',lw=0.5)
ax62.hist(background_RMiss_NuWro, density=False, bins=RMiss_p_33GeV_bins, weights=wts_RMiss_Pedro, histtype='step', color=background_color_NuWro, alpha=1, lw=2, hatch='o')
ax62.hist(energy_RMiss, density=False, bins=RMiss_Ar_47GeV_bins, weights=background_RMiss, color=background_color,alpha=opacity, edgecolor='black',lw=0.5)
ax62.hist(energy_RMiss, density=False, bins=RMiss_Ar_47GeV_bins, weights=background_RMiss, histtype='step',color=background_color,alpha=1,lw=2, hatch='/')

### Labels ###
ax61.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / d$R^T_\mathrm{Miss}$ / N)', fontsize=50)

ax61.text(0.98,0.70,r'\textbf{Signal}',ha='right',transform=ax61.transAxes, fontsize=40, color=signal_color)
ax61.text(0.98,0.65,r'\textbf{(This work)}',ha='right',transform=ax61.transAxes, fontsize=40, color=signal_color)
ax61.text(0.15,0.87,r'\textbf{Bkg (DUNE)}',ha='left',transform=ax61.transAxes, fontsize=40, color=background_color)
ax61.text(0.05,0.80,r'\textbf{Bkg (This work - NuWro)}',ha='left',transform=ax61.transAxes, fontsize=40, color=background_color_NuWro)

ax62.text(0.76,0.73,r'\textbf{Signal}',ha='right',transform=ax62.transAxes, fontsize=40, color=signal_color)
ax62.text(0.76,0.68,r'\textbf{(This work)}',ha='right',transform=ax62.transAxes, fontsize=40, color=signal_color)
ax62.text(0.15,0.87,r'\textbf{Bkg (DUNE)}',ha='left',transform=ax62.transAxes, fontsize=40, color=background_color)
ax62.text(0.05,0.80,r'\textbf{Bkg (This work - NuWro)}',ha='left',transform=ax62.transAxes, fontsize=40, color=background_color_NuWro)

ax61.set_title(r'\textbf{Incoherent - $\langle E_\nu \rangle = 33$ GeV}',fontsize=50)
ax62.set_title(r'\textbf{Coherent - $\langle E_\nu \rangle = 47$ GeV}',fontsize=50)
# Location of supxlabel; 15 is the figure height
y_xlabel = np.sqrt(15) / 100 * 1.1
y_title  = 0.96
fig6.supxlabel(r'\textbf{Missing Momentum Ratio} $R^T_{\mathrm{Miss}}$', fontsize=50, y=y_xlabel)
fig6.suptitle(r'$\tau$ \textbf{Hadronic Decay Channels}',fontsize=60, y=y_title)

ax61.set_xlim(0, 1)
ax62.set_xlim(0, 1)
ax61.set_ylim(0, 3.2)

# Tick labels #
xmajor_RMiss = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
xmajor_ThetaMiss = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]

ax61.set_xticks(xmajor_RMiss, labels=['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'])
ax61.tick_params(which='both',right=False, bottom=True)
ax61.tick_params(axis='x',labelsize=40)
ax61.tick_params(axis='y',labelsize=40)

ax62.set_xticks(xmajor_RMiss, labels=['', '0.2', '0.4', '0.6', '0.8', '1.0'])
ax62.tick_params(which='both',right=False, bottom=True)
ax62.tick_params(axis='x',labelsize=40)
ax62.tick_params(axis='y',labelsize=40)

# Save #
fig6.savefig("../../../plots/RMiss_1tau_hadronic_tau+_vtau_events.png", dpi=100, bbox_inches='tight')
