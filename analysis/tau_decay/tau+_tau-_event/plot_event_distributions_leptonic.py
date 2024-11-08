import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.patheffects as pe
import csv

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

DIST_DIR_1TAU = '../csv/distributions/vmu_to_vtau_tau+_mu-'

RMiss_background_filename = '../csv/events/AlexSousa_RMiss_Background.csv'
energy_RMiss = []
background_RMiss = []

RMiss_NuWro_background_filename = DIST_DIR_1TAU+'/RT-thresholds-no-single.txt'
background_RMiss_NuWro = []

write_dist = False

############################
#### Read Distributions ####
############################

def read_four_vectors_incoh(file_path, nucleon):
    if nucleon == 'proton':
        nucleon_pid = '2212'
    if nucleon == 'neutron':
        nucleon_pid = '2112'
    with open(file_path, 'r') as txtfile:
        nu_tau_array = []
        muon_array = []
        pVis_array = []
        nucleon_array = []
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
                pVis_array.append(pVis)
            if data[0] == nucleon_pid:  # Read and store nucleon
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                E  = float(data[5])
                nuc = np.array([px, py, pz, E])
                nucleon_array.append(nuc)

        if nucleon == 'proton':
            print("proton pMiss")
            print("\tpVis_array: ", len(pVis_array))
            print("\tnucleon_array: ", len(nucleon_array))
            print("\tmuon_array: ", len(muon_array))
            pMiss_array = [pVis + nucleon + muon for pVis, nucleon, muon in zip(pVis_array, nucleon_array, muon_array)]
        if nucleon == 'neutron':
            print("neutron pMiss")
            print("\tpVis_array: ", len(pVis_array))
            print("\tmuon_array: ", len(muon_array))
            pMiss_array = [pVis + muon for pVis, muon in zip(pVis_array, muon_array)]
        return nu_tau_array, muon_array, pMiss_array, nucleon_array

def read_four_vectors_coh(file_path):
    with open(file_path, 'r') as txtfile:
        nu_tau_array = []
        muon_array = []
        neutrinos_array = []
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
            if data[0] == '999':  # Read and store nu_tau_bar
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                E  = float(data[5])
                neutrinos = np.array([px, py, pz, E])
                neutrinos_array.append(neutrinos)

        print("read coherent; nu_tau_array: ",len(nu_tau_array))
        print("read coherent; neutrinos_array: ",len(neutrinos_array))
        pMiss_array = [nu_tau + neutrinos for nu_tau, neutrinos in zip(nu_tau_array, neutrinos_array)]
        return nu_tau_array, muon_array, neutrinos_array, pMiss_array

def calculate_ptMiss(pMiss_array):
    ptMiss_array = []
    for pMiss in pMiss_array:
        pxMiss = pMiss[0]
        pyMiss = pMiss[1]
        ptMiss = np.sqrt(pxMiss**2 + pyMiss**2)
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

nu_tau_p_33GeV, muon_p_33GeV, pMiss_p_33GeV, nucleon_p_33GeV = read_four_vectors_incoh(DIST_DIR_1TAU+'/nucleon/proton/tau_leptonic_decayed_distribution_p_33GeV.txt', 'proton')
nu_tau_n_33GeV, muon_n_33GeV, pMiss_n_33GeV, nucleon_n_33GeV = read_four_vectors_incoh(DIST_DIR_1TAU+'/nucleon/neutron/tau_leptonic_decayed_distribution_n_33GeV.txt', 'neutron')
nu_tau_Ar_47GeV, muon_Ar_47GeV, neutrinos_Ar_47GeV, pMiss_Ar_47GeV = read_four_vectors_coh(DIST_DIR_1TAU+'/coherent/argon/tau_leptonic_decayed_distribution_Ar_47GeV.txt')

print(len(pMiss_p_33GeV))
print(len(pMiss_n_33GeV))
print(len(pMiss_Ar_47GeV))

ptMiss_p_33GeV = calculate_ptMiss(pMiss_p_33GeV)
ptMiss_n_33GeV = calculate_ptMiss(pMiss_n_33GeV)
ptMiss_Ar_47GeV = calculate_ptMiss(pMiss_Ar_47GeV)

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
#NEVENTS = 1e5
NEVENTS = len(ptMiss_p_33GeV)
print("Total events: ", NEVENTS)

# proton - 33 GeV #
wts_ptMiss_p_33GeV = np.ones_like(ptMiss_p_33GeV) / (NEVENTS * ptMiss_p_33GeV_step)
wts_ptMuon_p_33GeV = np.ones_like(ptMuon_p_33GeV) / (NEVENTS * ptMuon_p_33GeV_step)
wts_RMiss_p_33GeV = np.ones_like(RMiss_p_33GeV) / (NEVENTS * RMiss_p_33GeV_step)
wts_ThetaMiss_p_33GeV = np.ones_like(ThetaMiss_p_33GeV) / (NEVENTS * ThetaMiss_p_33GeV_step)

# neutron - 33 GeV #
wts_ptMiss_n_33GeV = np.ones_like(ptMiss_n_33GeV) / (NEVENTS * ptMiss_n_33GeV_step)
wts_ptMuon_n_33GeV = np.ones_like(ptMuon_n_33GeV) / (NEVENTS * ptMuon_n_33GeV_step)
wts_RMiss_n_33GeV = np.ones_like(RMiss_n_33GeV) / (NEVENTS * RMiss_n_33GeV_step)
wts_ThetaMiss_n_33GeV = np.ones_like(ThetaMiss_n_33GeV) / (NEVENTS * ThetaMiss_n_33GeV_step)

# incoherent - 33 GeV #
wts_RMiss_incoh_33GeV = np.ones_like(RMiss_incoh_33GeV) / (2 * NEVENTS * RMiss_n_33GeV_step)

# argon - 47 GeV #
wts_ptMiss_Ar_47GeV = np.ones_like(ptMiss_Ar_47GeV) / (NEVENTS * ptMiss_Ar_47GeV_step)
wts_ptMuon_Ar_47GeV = np.ones_like(ptMuon_Ar_47GeV) / (NEVENTS * ptMuon_Ar_47GeV_step)
wts_RMiss_Ar_47GeV = np.ones_like(RMiss_Ar_47GeV) / (NEVENTS * RMiss_Ar_47GeV_step)
wts_ThetaMiss_Ar_47GeV = np.ones_like(ThetaMiss_Ar_47GeV) / (NEVENTS * ThetaMiss_Ar_47GeV_step)

# background - Pedro #
wts_RMiss_Pedro = np.ones_like(background_RMiss_NuWro) / (2922963 * RMiss_p_33GeV_step)

##################
#### Plotting ####
##################

fig1, (ax11, ax12, ax13) = plt.subplots(1, 3, sharey=True, figsize=(36,15), tight_layout=True, gridspec_kw={'wspace':0}) # ptMiss
fig2, (ax21, ax22, ax23) = plt.subplots(1, 3, sharey=True, figsize=(36,15), tight_layout=True, gridspec_kw={'wspace':0}) # ptMuon
fig4, (ax41, ax42, ax43) = plt.subplots(1, 3, sharey=True, figsize=(36,15), tight_layout=True, gridspec_kw={'wspace':0}) # RMiss
fig5, (ax51, ax52, ax53) = plt.subplots(1, 3, sharey=True, figsize=(36,15), tight_layout=True, gridspec_kw={'wspace':0}) # ThetaMiss
fig6, (ax61, ax62) = plt.subplots(1, 2, sharey=True, figsize=(24,15), tight_layout=True, gridspec_kw={'wspace':0}) # RMiss incoh and coh

#signal_color = '#3E0077'
signal_color = '#2B4ACA'
#background_color = '#FFA500'
background_color = '#E70C64'
background_color_NuWro = '#DE7D39'
#signal_color = '#91B510'
#background_color = '#2CA5E3'
opacity = 0.50

### ptMiss ###
ax11.hist(ptMiss_p_33GeV, density=False, bins=ptMiss_p_33GeV_bins, weights=wts_ptMiss_p_33GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
ax11.hist(ptMiss_p_33GeV, density=False, bins=ptMiss_p_33GeV_bins, weights=wts_ptMiss_p_33GeV, histtype='step', color='black', alpha=1, lw=2)
ax12.hist(ptMiss_n_33GeV, density=False, bins=ptMiss_n_33GeV_bins, weights=wts_ptMiss_n_33GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
ax12.hist(ptMiss_n_33GeV, density=False, bins=ptMiss_n_33GeV_bins, weights=wts_ptMiss_n_33GeV, histtype='step', color='black', alpha=1, lw=2)
ax13.hist(ptMiss_Ar_47GeV, density=False, bins=ptMiss_Ar_47GeV_bins, weights=wts_ptMiss_Ar_47GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
ax13.hist(ptMiss_Ar_47GeV, density=False, bins=ptMiss_Ar_47GeV_bins, weights=wts_ptMiss_Ar_47GeV, histtype='step', color='black', alpha=1, lw=2)

### ptMuon ###
ax21.hist(ptMuon_p_33GeV, density=False, bins=ptMuon_p_33GeV_bins, weights=wts_ptMuon_p_33GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
ax21.hist(ptMuon_p_33GeV, density=False, bins=ptMuon_p_33GeV_bins, weights=wts_ptMuon_p_33GeV, histtype='step', color='black', alpha=1, lw=2)
ax22.hist(ptMuon_n_33GeV, density=False, bins=ptMuon_n_33GeV_bins, weights=wts_ptMuon_n_33GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
ax22.hist(ptMuon_n_33GeV, density=False, bins=ptMuon_n_33GeV_bins, weights=wts_ptMuon_n_33GeV, histtype='step', color='black', alpha=1, lw=2)
ax23.hist(ptMuon_Ar_47GeV, density=False, bins=ptMuon_Ar_47GeV_bins, weights=wts_ptMuon_Ar_47GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
ax23.hist(ptMuon_Ar_47GeV, density=False, bins=ptMuon_Ar_47GeV_bins, weights=wts_ptMuon_Ar_47GeV, histtype='step', color='black', alpha=1, lw=2)

### RMiss ###
ax41.hist(RMiss_p_33GeV, density=False, bins=RMiss_p_33GeV_bins, weights=wts_RMiss_p_33GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
ax41.hist(RMiss_p_33GeV, density=False, bins=RMiss_p_33GeV_bins, weights=wts_RMiss_p_33GeV, histtype='step', color='black', alpha=1, lw=2)
#ax41.hist(energy_RMiss, density=False, bins=RMiss_p_33GeV_bins, weights=background_RMiss, color='firebrick',alpha=0.3, edgecolor='black',lw=0.5)
#ax41.hist(energy_RMiss, density=False, bins=RMiss_p_33GeV_bins, weights=background_RMiss, histtype='step',color='black',alpha=1,lw=2)
ax42.hist(RMiss_n_33GeV, density=False, bins=RMiss_n_33GeV_bins, weights=wts_RMiss_n_33GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
ax42.hist(RMiss_n_33GeV, density=False, bins=RMiss_n_33GeV_bins, weights=wts_RMiss_n_33GeV, histtype='step', color='black', alpha=1, lw=2)
#ax42.hist(energy_RMiss, density=False, bins=RMiss_p_33GeV_bins, weights=background_RMiss, color='firebrick',alpha=0.3, edgecolor='black',lw=0.5)
#ax42.hist(energy_RMiss, density=False, bins=RMiss_p_33GeV_bins, weights=background_RMiss, histtype='step',color='black',alpha=1,lw=2)
ax43.hist(RMiss_Ar_47GeV, density=False, bins=RMiss_Ar_47GeV_bins, weights=wts_RMiss_Ar_47GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
ax43.hist(RMiss_Ar_47GeV, density=False, bins=RMiss_Ar_47GeV_bins, weights=wts_RMiss_Ar_47GeV, histtype='step', color='black', alpha=1, lw=2)
#ax43.hist(energy_RMiss, density=False, bins=RMiss_Ar_47GeV_bins, weights=background_RMiss, color='firebrick',alpha=0.3, edgecolor='black',lw=0.5)
#ax43.hist(energy_RMiss, density=False, bins=RMiss_Ar_47GeV_bins, weights=background_RMiss, histtype='step',color='black',alpha=1,lw=2)

### ThetaMiss ###
ax51.hist(ThetaMiss_p_33GeV, density=False, bins=ThetaMiss_p_33GeV_bins, weights=wts_ThetaMiss_p_33GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
ax51.hist(ThetaMiss_p_33GeV, density=False, bins=ThetaMiss_p_33GeV_bins, weights=wts_ThetaMiss_p_33GeV, histtype='step', color='black', alpha=1, lw=2)
ax52.hist(ThetaMiss_n_33GeV, density=False, bins=ThetaMiss_n_33GeV_bins, weights=wts_ThetaMiss_n_33GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
ax52.hist(ThetaMiss_n_33GeV, density=False, bins=ThetaMiss_n_33GeV_bins, weights=wts_ThetaMiss_n_33GeV, histtype='step', color='black', alpha=1, lw=2)
ax53.hist(ThetaMiss_Ar_47GeV, density=False, bins=ThetaMiss_Ar_47GeV_bins, weights=wts_ThetaMiss_Ar_47GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
ax53.hist(ThetaMiss_Ar_47GeV, density=False, bins=ThetaMiss_Ar_47GeV_bins, weights=wts_ThetaMiss_Ar_47GeV, histtype='step', color='black', alpha=1, lw=2)

### RMiss ###
ax61.hist(RMiss_incoh_33GeV, density=False, bins=RMiss_p_33GeV_bins, weights=wts_RMiss_incoh_33GeV, color=signal_color, alpha=opacity, edgecolor='black', lw=0.5)
ax61.hist(RMiss_incoh_33GeV, density=False, bins=RMiss_p_33GeV_bins, weights=wts_RMiss_incoh_33GeV, histtype='step', color=signal_color, alpha=1, lw=2)
ax61.hist(background_RMiss_NuWro, density=False, bins=RMiss_p_33GeV_bins, weights=wts_RMiss_Pedro, color=background_color_NuWro,alpha=opacity, edgecolor='black',lw=0.5)
ax61.hist(background_RMiss_NuWro, density=False, bins=RMiss_p_33GeV_bins, weights=wts_RMiss_Pedro, histtype='step', color=background_color_NuWro, alpha=1, lw=2, hatch='o')
ax61.hist(energy_RMiss, density=False, bins=RMiss_p_33GeV_bins, weights=background_RMiss, color=background_color,alpha=opacity, edgecolor='black',lw=0.5)
ax61.hist(energy_RMiss, density=False, bins=RMiss_p_33GeV_bins, weights=background_RMiss, histtype='step',color=background_color,alpha=1,lw=2, hatch='/')
ax62.hist(RMiss_Ar_47GeV, density=False, bins=RMiss_Ar_47GeV_bins, weights=wts_RMiss_Ar_47GeV, color=signal_color, alpha=opacity, edgecolor='black', lw=0.5)
ax62.hist(RMiss_Ar_47GeV, density=False, bins=RMiss_Ar_47GeV_bins, weights=wts_RMiss_Ar_47GeV, histtype='step', color=signal_color, alpha=1, lw=2)
ax62.hist(energy_RMiss, density=False, bins=RMiss_Ar_47GeV_bins, weights=background_RMiss, color=background_color,alpha=opacity, edgecolor='black',lw=0.5)
ax62.hist(energy_RMiss, density=False, bins=RMiss_Ar_47GeV_bins, weights=background_RMiss, histtype='step',color=background_color,alpha=1,lw=2, hatch='/')
ax62.hist(background_RMiss_NuWro, density=False, bins=RMiss_p_33GeV_bins, weights=wts_RMiss_Pedro, color=background_color_NuWro,alpha=opacity, edgecolor='black',lw=0.5)
ax62.hist(background_RMiss_NuWro, density=False, bins=RMiss_p_33GeV_bins, weights=wts_RMiss_Pedro, histtype='step', color=background_color_NuWro, alpha=1, lw=2, hatch='o')

### Labels ###
ax11.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / bin / N)', fontsize=40)
ax21.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / bin / N)', fontsize=40)
ax41.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / d$R^T_\mathrm{Miss}$ / N)', fontsize=40)
ax51.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / d$\theta_\mathrm{Miss}$ / N)', fontsize=40)
ax61.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / d$R^T_\mathrm{Miss}$ / N)', fontsize=50)

for ax in [ax11, ax12, ax13]:
    ax.set_xlabel(r'\textbf{Missing Transverse Momentum} $p^T_{\mathrm{Miss}}$ (GeV) ', fontsize=40)

for ax in [ax21, ax22, ax23]:
    ax.set_xlabel(r'\textbf{Muon Transverse Momentum} $p^T_\mu$ (GeV)', fontsize=40)

for ax in [ax41, ax42, ax43]:
    ax.set_xlabel(r'\textbf{Missing Momentum Ratio} $R^T_{\mathrm{Miss}}$', fontsize=40)

for ax in [ax51, ax52, ax53]:
    ax.set_xlabel(r'\textbf{Missing Polar Angle} $\theta_{\mathrm{Miss}}$', fontsize=40)

for ax in [ax11, ax21, ax51]:
    ax.text(0.90,0.90,r'\textbf{Incoherent p - 33 GeV}',ha='right',transform=ax.transAxes, fontsize=40)

for ax in [ax12, ax22, ax52]:
    ax.text(0.90,0.90,r'\textbf{Incoherent n - 33 GeV}',ha='right',transform=ax.transAxes, fontsize=40)

for ax in [ax13, ax23, ax53]:
    ax.text(0.90,0.90,r'\textbf{Coherent Ar - 47 GeV}',ha='right',transform=ax.transAxes, fontsize=40)

ax61.text(0.98,0.70,r'\textbf{Signal}',ha='right',transform=ax61.transAxes, fontsize=40, color=signal_color)
ax61.text(0.98,0.65,r'\textbf{(This work)}',ha='right',transform=ax61.transAxes, fontsize=40, color=signal_color)
ax61.text(0.15,0.90,r'\textbf{Bkg (ARRS - DUNE)}',ha='left',transform=ax61.transAxes, fontsize=40, color=background_color)
ax61.text(0.05,0.83,r'\textbf{Bkg (This work - NuWro)}',ha='left',transform=ax61.transAxes, fontsize=40, color=background_color_NuWro)

ax62.text(0.72,0.76,r'\textbf{Signal}',ha='right',transform=ax62.transAxes, fontsize=40, color=signal_color)
ax62.text(0.72,0.71,r'\textbf{(This work)}',ha='right',transform=ax62.transAxes, fontsize=40, color=signal_color)
ax62.text(0.15,0.90,r'\textbf{Bkg (ARRS - DUNE)}',ha='left',transform=ax62.transAxes, fontsize=40, color=background_color)
ax62.text(0.05,0.83,r'\textbf{Bkg (This work - NuWro)}',ha='left',transform=ax62.transAxes, fontsize=40, color=background_color_NuWro)

ax61.set_title(r'\textbf{Incoherent - $\langle E_\nu \rangle = 33$ GeV}',fontsize=50)
ax62.set_title(r'\textbf{Coherent - $\langle E_\nu \rangle = 47$ GeV}',fontsize=50)
# Location of supxlabel; 15 is the figure height
y_xlabel = np.sqrt(15) / 100 * 1.1
#y_title  = np.sqrt(15) / 100 * 9.0
y_title = 0.96
fig6.supxlabel(r'\textbf{Missing Momentum Ratio} $R^T_{\mathrm{Miss}}$', fontsize=50, y=y_xlabel)
fig6.suptitle(r'$\tau$ \textbf{Leptonic Decay Channels}', fontsize=60, y=y_title)

ax41.text(0.10,0.90,r'\textbf{Incoherent p - 33 GeV}',ha='left',transform=ax41.transAxes,fontsize=40)
ax42.text(0.10,0.90,r'\textbf{Incoherent n - 33 GeV}',ha='left',transform=ax42.transAxes,fontsize=40)
ax43.text(0.10,0.90,r'\textbf{Coherent Ar - 47 GeV}',ha='left',transform=ax43.transAxes,fontsize=40)

### Limits ###
# ptMiss #
ax11.set_xlim(ptMiss_p_33GeV_min, ptMiss_p_33GeV_max)
ax12.set_xlim(ptMiss_n_33GeV_min, ptMiss_n_33GeV_max)
ax13.set_xlim(ptMiss_Ar_47GeV_min, ptMiss_Ar_47GeV_max)

# ptMuon #
ax21.set_xlim(ptMuon_p_33GeV_min, ptMuon_p_33GeV_max)
ax22.set_xlim(ptMuon_n_33GeV_min, ptMuon_n_33GeV_max)
ax23.set_xlim(ptMuon_Ar_47GeV_min, ptMuon_Ar_47GeV_max)

# RMiss #
ax41.set_xlim(0, 1)
ax42.set_xlim(0, 1)
ax43.set_xlim(0, 1)

ax61.set_xlim(0, 1)
ax62.set_xlim(0, 1)
ax61.set_ylim(0, 3.0)

# Tick labels #
xmajor_RMiss = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
xmajor_ThetaMiss = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]

ax41.set_xticks(xmajor_RMiss, labels=['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'])
ax41.tick_params(which='both',right=False, bottom=True)
ax42.set_xticks(xmajor_RMiss, labels=['', '0.2', '0.4', '0.6', '0.8', '1.0'])
ax42.tick_params(which='both',right=False, bottom=True)
ax43.set_xticks(xmajor_RMiss, labels=['', '0.2', '0.4', '0.6', '0.8', '1.0'])
ax43.tick_params(which='both',right=False, bottom=True)

ax51.set_xticks(xmajor_ThetaMiss, labels=['0', '20', '40', '60', '80', '100', '120', '140', '160', '180'])
ax51.tick_params(which='both',right=False, bottom=True)
ax52.set_xticks(xmajor_ThetaMiss, labels=['', '20', '40', '60', '80', '100', '120', '140', '160', '180'])
ax52.tick_params(which='both',right=False, bottom=True)
ax53.set_xticks(xmajor_ThetaMiss, labels=['', '20', '40', '60', '80', '100', '120', '140', '160', '180'])
ax53.tick_params(which='both',right=False, bottom=True)

ax61.set_xticks(xmajor_RMiss, labels=['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'])
ax61.tick_params(which='both',right=False, bottom=True)
ax61.tick_params(axis='x',labelsize=40)
ax61.tick_params(axis='y',labelsize=40)

ax62.set_xticks(xmajor_RMiss, labels=['', '0.2', '0.4', '0.6', '0.8', '1.0'])
ax62.tick_params(which='both',right=False, bottom=True)
ax62.tick_params(axis='x',labelsize=40)
ax62.tick_params(axis='y',labelsize=40)

# Save #
fig1.savefig("../plots/ptMiss_1tau_leptonic.png", dpi=100, bbox_inches='tight')
fig2.savefig("../plots/ptMuon_1tau_leptonic.png", dpi=100, bbox_inches='tight')
fig4.savefig("../plots/RMiss_1tau_leptonic.png", dpi=100, bbox_inches='tight')
fig5.savefig("../plots/ThetaMiss_1tau_leptonic.png", dpi=100, bbox_inches='tight')
fig6.savefig("../plots/RMiss_1tau_leptonic.png", dpi=100, bbox_inches='tight')
