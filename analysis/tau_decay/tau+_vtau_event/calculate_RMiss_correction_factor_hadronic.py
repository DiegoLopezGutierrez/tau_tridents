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
background_value_RMiss = []
background_weight_RMiss = []

RMiss_NuWro_background_filename = DIST_DIR_1TAU+'/RT-thresholds-no-single.txt'
background_value_RMiss_NuWro = []

write_dist = False

Ev_array = [5,8,11,14,17,20,24,27,30,33,36,39,42,45,49,52,55,58,61,64,67,71,74,77,80,83,86,89,92,96,99,102,105,108,111,114,118,121,124]

R_cut = float(input("Cut value for RMiss: "))

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
                assert E <= float(Ev), "Unphysical energy in line {}: {}".format(i, E)
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
    print("Reading file: ", file_path)
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

def calculate_RMiss(ptMiss_array, ptMuon_array, cut):
    RMiss_array = []
    ntotal = 0
    npass = 0
    for ptMiss, ptMuon in zip(ptMiss_array, ptMuon_array):
        RMiss = ptMiss / (ptMiss + ptMuon)
        if RMiss >= cut:
            npass += 1
        RMiss_array.append(RMiss)
        ntotal += 1
    return RMiss_array, npass, ntotal

############################
###### Distributions #######
############################

#correction_factor_incoh_array = []
correction_factor_p_array = []
correction_factor_n_array = []
correction_factor_Ar_array = []

for Ev in Ev_array:
    nu_tau_p, muon_p, pMiss_p, nucleon_p = read_four_vectors_incoh(DIST_DIR_1TAU+f'/nucleon/proton/tau+_vtau_events/hadronic/tau_hadronic_decayed_distribution_p_{Ev}GeV.txt', 'proton')
    nu_tau_n, muon_n, pMiss_n, nucleon_n = read_four_vectors_incoh(DIST_DIR_1TAU+f'/nucleon/neutron/tau+_vtau_events/hadronic/tau_hadronic_decayed_distribution_n_{Ev}GeV.txt', 'neutron')
    nu_tau_Ar, muon_Ar, nubar_tau_Ar, pMiss_Ar = read_four_vectors_coh(DIST_DIR_1TAU+f'/coherent/argon/tau+_vtau_events/hadronic/tau_hadronic_decayed_distribution_Ar_{Ev}GeV.txt')

    ptMiss_p = calculate_ptMiss(pMiss_p, float(Ev))
    ptMiss_n = calculate_ptMiss(pMiss_n, float(Ev))
    ptMiss_Ar = calculate_ptMiss(pMiss_Ar, float(Ev))

    pzMiss_p = calculate_pzMiss(pMiss_p)
    pzMiss_n = calculate_pzMiss(pMiss_n)
    pzMiss_Ar = calculate_pzMiss(pMiss_Ar)

    ptMuon_p = calculate_ptMuon(muon_p)
    ptMuon_n = calculate_ptMuon(muon_n)
    ptMuon_Ar = calculate_ptMuon(muon_Ar)

    RMiss_p, npass_p, ntotal_p = calculate_RMiss(ptMiss_p, ptMuon_p, R_cut)
    RMiss_n, npass_n, ntotal_n = calculate_RMiss(ptMiss_n, ptMuon_n, R_cut)
#    RMiss_incoh = RMiss_p + RMiss_n
    RMiss_Ar, npass_Ar, ntotal_Ar = calculate_RMiss(ptMiss_Ar, ptMuon_Ar, R_cut)

#    correction_factor_incoh = (npass_p + npass_n) / (ntotal_p + ntotal_n)
    correction_factor_p  = npass_p / ntotal_p
    correction_factor_n  = npass_n / ntotal_n
    correction_factor_Ar = npass_Ar / ntotal_Ar

#    correction_factor_incoh_array.append(correction_factor_incoh)
    correction_factor_p_array.append(correction_factor_p)
    correction_factor_n_array.append(correction_factor_n)
    correction_factor_Ar_array.append(correction_factor_Ar)

with open('RMiss_correction_factors_hadronic.txt','w',newline='') as outfile:
    writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(["Neutrino Energy [GeV]", "Proton Correction Factor", "Neutron Correction Factor", "Argon Correction Factor"])
    for Ev, correction_factor_p, correction_factor_n, correction_factor_Ar in zip(Ev_array, correction_factor_p_array, correction_factor_n_array, correction_factor_Ar_array):
        writer.writerow([Ev, correction_factor_p, correction_factor_n, correction_factor_Ar])

## TO-DO: Consider calculating the correction factor for background events as well.
#with open(RMiss_background_filename,'r') as csvfile:
#    data = csv.reader(csvfile, delimiter=',')
#    for row in data:
#        background_value_RMiss.append(float(row[0]))     # This is the actual RMiss value (i.e. x-value)
#        background_weight_RMiss.append(float(row[1]))    # Since the histogram from Alex Sousa is digitized, these are the weights (i.e. y-value) of the histogram.
#
#with open(RMiss_NuWro_background_filename,'r') as txtfile:
#    for line in txtfile:
#        background_value_RMiss_NuWro.append(float(line))  # Pedro's NuWro simulation is just the actual RMiss values (i.e. x-values) that need to be binned and histogramed8

##################
#### Plotting ####
##################

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, figsize=(30,15), tight_layout=True, gridspec_kw={'wspace':0}) # RMiss incoh and coh

coh_color = '#2B4ACA'
p_color = '#E70C64'
n_color = '#DE7D39'

### RMiss ###
ax1.plot(Ev_array, correction_factor_p_array, '-', color=p_color, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax2.plot(Ev_array, correction_factor_n_array, '-', color=n_color, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax3.plot(Ev_array, correction_factor_Ar_array, '-', color=coh_color, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])

### Labels ###
ax1.set_ylabel(r'Correction Factor $\varepsilon_{R^T_\mathrm{Miss}}$', fontsize=50)

ax1.text(0.90,0.90,r'\textbf{Proton}',ha='right',transform=ax1.transAxes, fontsize=40, color=p_color)
ax2.text(0.90,0.90,r'\textbf{Neutron}',ha='right',transform=ax2.transAxes, fontsize=40, color=n_color)
ax3.text(0.90,0.90,r'\textbf{Argon}',ha='right',transform=ax3.transAxes, fontsize=40, color=coh_color)

# Location of supxlabel; 15 is the figure height
y_xlabel = np.sqrt(15) / 100 * 1.1
y_title  = 0.96
fig.supxlabel(r'\textbf{Neutrino Energy} $E_\nu$', fontsize=50, y=y_xlabel)
fig.suptitle(r'$\tau$ \textbf{Hadronic Decay Channels}, $R^T_\mathrm{Miss} \geq$ '+str(R_cut),fontsize=60, y=y_title)

fig.savefig("correction_factor_RMiss.png", dpi=100)
