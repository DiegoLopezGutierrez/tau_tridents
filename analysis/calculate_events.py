import numpy as np
import csv
from scipy.integrate import simpson, trapezoid
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe
import matplotlib as mpl

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

#######################
###### Constants ######
#######################

atomic_mass_unit = 1.6605e-27 # kg
MD = 1000                     # 1 tonne in kg
mmu = 0.105                   # mass of muon in GeV
GFermi = 1.1663787e-5         # Fermi constant in GeV^-2
aEM = 1/137.036               # Fine structure constant
sW2 = 0.23119                 # sin(theta_Weinberg)^2

me = 5.1099e-4 # GeV
mmu = 0.1056   # GeV
mtau = 1.7769  # GeV

## DUNE ##
N_POT = 1.1e21               # Number of POTs
MAr = 39.95*atomic_mass_unit # Mass of argon in atomic mass units
A_Ar = 40
Z_Ar = 18

## FASERnu ##
L_Run = 150                  # Run luminosity at FASERv; 150 fb^-1
L_Run2 = 3000                # Run luminosity at FASERv2; 3 ab^-1
MW = 183.84*atomic_mass_unit # Mass of tungsten in atomic mass units
A_W = 184
Z_W = 74
M_FASERv = 1.1               # Older papers report a FASERv mass of 1.2 tonnes but the current (03/2024) value is 1.1 instead.
M_FASERv2 = 20               # The Forward Physics Facility at the High Luminosity LHC (Feng et al. 2023) reports 20 tonnes for FASERv2.
RHO_W = 19300                # Density of tungsten in kg/m^3
Length_FASERv = 0.90         # FASERv will have a length of 0.90 m, which together with a cross-sectional area of 25 x 25 cm^2 leads to a total tungsten mass of 1.1 tonnes.

## MINOS ## Note: MINOS+ never ran in antineutrino mode.
N_POT_MINOS_neutrino = 10.56             # Total neutrino mode exposure of the MINOS detector throughout its lifetime is 10.56e20. Flux file is normalized to 1e20 POT already.
N_POT_MINOS_antineutrino = 3.36          # Total antineutrino mode exposure of the MINOS detector throughout its lifetime is 3.36e20. Flux file is normalized to 1e20 POT already.
N_POT_MINOSPlus_neutrino = 9.69          # Total neutrino mode exposure of the MINOS+ detector throughout its lifetime is 9.69e20. Flux file is normalized to 1e20 POT already.
MFe = 55.85*atomic_mass_unit             # Mass of iron in atomic mass units
A_Fe = 56
Z_Fe = 26
M_MINOS = 28.6               # MINOS and MINOS+ had a total detector mass of 980t of iron. Ballett et al. assume a fiducial volume of 28.6t for their analysis, which is standard.

## T2K -- INGRID ##
N_POT_INGRID = 3.9e1        # Total neutrino/antineutrino mode exposure of the INGRID detector throughout its lifetime is 3.9e21. Flux file is normalized to 1e20 POT already.
N_POT_INGRID_phase2 = 1e2   # Total neutrino/antineutrino mode exposure of the INGRID detector phase 2 throughout its lifetime is 1.0e22. Flux file is normalized to 1e20 POT already.
M_INGRID = 99.4

## SHiP ##
N_POT_SHiP = 2e20   # Number of POTs used by Magill and Plestid. Number given by SHiP people seems to be 5e13 instead.
M_SHiP = 8          # Total detector target mass of about 8 tonnes per arXiv:2112.01487.

## SBND ##
N_POT_SBND = 6.6    # Number of POTs used by Ballett et al. Flux file is normalized to 1e20 POT already.
M_SBND = 112        # Total detector mass of 112 tonnes of argon per their webiste. (Also used by Ballett).

## Cross Section Uncertainties ##
# Components # (see Altmannshofer et al. for details)
aEM = 1/137  # low q^2 usually, so using zero momentum value of fine structure
sigma_highQED_Ar = Z_Ar*aEM/(4*np.pi) + 0.02  # higher order QED corrections should give roughly Z*aEM/(4*pi) which is about 1% for Ar. Add 2% to be conservative.
sigma_highQED_W  = Z_W*aEM/(4*np.pi) + 0.02 # higher order QED corrections should give roughly Z*aEM/(4*pi) which is about 4% for W. Add 2% to be conservative.
sigma_highQED_Fe  = Z_Fe*aEM/(4*np.pi) + 0.02 # higher order QED corrections should give roughly Z*aEM/(4*pi) which is about 4% for W. Add 2% to be conservative.
sigma_highQED_p = 1*aEM/(4*np.pi) + 0.02
sigma_highQED_n = 0

sigma_form_factors_coh = 0.01 # 1% uncertainty on coherent cross sections due to form factors
sigma_form_factors_incoh = 0.03 # 3% uncertainty on incoherent cross sections due to form factors

sigma_highEW = 0.05 # 5% uncertainty due to higher order weak corrections

sigma_nuclear_modeling = 0.30 # for incoherent scattering, largest uncertainty due to other nuclear effects besides the included Pauli blocking. 30% to be conservative.

# Total #
sigma_total_coh_Ar = np.sqrt(sigma_highQED_Ar**2 + sigma_form_factors_coh**2 + sigma_highEW**2)
sigma_total_coh_W  = np.sqrt(sigma_highQED_W**2 + sigma_form_factors_coh**2 + sigma_highEW**2)
sigma_total_coh_Fe  = np.sqrt(sigma_highQED_Fe**2 + sigma_form_factors_coh**2 + sigma_highEW**2)

sigma_total_incoh_Ar = np.sqrt(sigma_highQED_Ar**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)
sigma_total_incoh_W  = np.sqrt(sigma_highQED_W**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)
sigma_total_incoh_Fe  = np.sqrt(sigma_highQED_Fe**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)

sigma_total_p = np.sqrt(sigma_highQED_p**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)
sigma_total_n = np.sqrt(sigma_highQED_n**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)

########################
###### Initialize ######
########################

### Directories ###
FLUX_DIR = '../csv/fluxes'
CROSS_SECTION_DIR = '../csv/cross_sections'
CORRECTION_FACTOR_DIR = './tau_decay/tau+_vtau_event'

### Fluxes ###
DUNE_neutrino_filename = FLUX_DIR + '/DUNE/histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_LBNEND_globes_flux.txt'
DUNE_antineutrino_filename = FLUX_DIR + '/DUNE/histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_antineutrino_LBNEND_globes_flux.txt'
DUNE_tau_opt_neutrino_filename = FLUX_DIR + '/DUNE/histos_g4lbne_v3r5p4_QGSP_BERT_TauOptimized_neutrino_LBNEND_globes_flux.txt'
DUNE_tau_opt_antineutrino_filename = FLUX_DIR + '/DUNE/histos_g4lbne_v3r5p4_QGSP_BERT_TauOptimized_antineutrino_LBNEND_globes_flux.txt'

FASERnu_filename = FLUX_DIR + '/FASERnu/vmu/FASERvmu.csv'
FASERnubar_filename = FLUX_DIR + '/FASERnu/vmubar/FASERvmubar.csv'
FASERnu2_filename = FLUX_DIR + '/FASERnu/FASERnu2/FASERv2_vmu_vmubar_flux.csv'

MINOS_neutrino_vmu_filename = FLUX_DIR + '/MINOS/MINOS_ND/neutrino_mode/MINOS_ND_neutrino_vmu_flux.csv'
MINOS_neutrino_vmubar_filename = FLUX_DIR + '/MINOS/MINOS_ND/neutrino_mode/MINOS_ND_neutrino_vmubar_flux.csv'
MINOS_antineutrino_vmu_filename = FLUX_DIR + '/MINOS/MINOS_ND/antineutrino_mode/MINOS_ND_antineutrino_vmu_flux.csv'
MINOS_antineutrino_vmubar_filename = FLUX_DIR + '/MINOS/MINOS_ND/antineutrino_mode/MINOS_ND_antineutrino_vmubar_flux.csv'

MINOSPlus_neutrino_vmu_filename = FLUX_DIR + '/MINOS/MINOS+_ND/neutrino_mode/MINOS+_ND_neutrino_vmu_flux.csv'
MINOSPlus_neutrino_vmubar_filename = FLUX_DIR + '/MINOS/MINOS+_ND/neutrino_mode/MINOS+_ND_neutrino_vmubar_flux.csv'
MINOSPlus_antineutrino_vmu_filename = FLUX_DIR + '/MINOS/MINOS+_ND/antineutrino_mode/MINOS+_ND_antineutrino_vmu_flux.csv'
MINOSPlus_antineutrino_vmubar_filename = FLUX_DIR + '/MINOS/MINOS+_ND/antineutrino_mode/MINOS+_ND_antineutrino_vmubar_flux.csv'

T2K_INGRID_filename = FLUX_DIR + '/T2K/INGRID/T2K_INGRID_neutrino_flux.csv'

SHiP_filename = FLUX_DIR + '/SHiP/SHiP.csv'

SBND_vmu_filename = FLUX_DIR + '/SBND/SBND_vmu_flux.csv'
SBND_vmubar_filename = FLUX_DIR + '/SBND/SBND_vmubar_flux.csv'
SBND_ve_filename = FLUX_DIR + '/SBND/SBND_ve_flux.csv'
SBND_vebar_filename = FLUX_DIR + '/SBND/SBND_vebar_flux.csv'

# DUNE Standard #
energy_DUNE = []
flux_DUNE_neutrino_vmu = []
flux_DUNE_neutrino_ve = []
flux_DUNE_neutrino_vmubar = []
flux_DUNE_neutrino_vebar = []

flux_DUNE_antineutrino_vmu = []
flux_DUNE_antineutrino_ve = []
flux_DUNE_antineutrino_vmubar = []
flux_DUNE_antineutrino_vebar = []

# DUNE Tau Optimized #
energy_DUNE_tau_opt = []
flux_DUNE_tau_opt_neutrino_vmu = []
flux_DUNE_tau_opt_neutrino_ve = []
flux_DUNE_tau_opt_neutrino_vmubar = []
flux_DUNE_tau_opt_neutrino_vebar = []

flux_DUNE_tau_opt_antineutrino_vmu = []
flux_DUNE_tau_opt_antineutrino_ve = []
flux_DUNE_tau_opt_antineutrino_vmubar = []
flux_DUNE_tau_opt_antineutrino_vebar = []

# FASERv #
energy_FASERvmu = []
bins_FASERvmu = []
flux_FASERvmu = []
flux_FASERvmubar = []

# FASERv2 #
energy_FASERv2 = []
flux_FASERv2_vmu_vmubar = []

# MINOS #
energy_MINOS = []
flux_MINOS_neutrino_vmu = []
flux_MINOS_neutrino_vmubar = []

flux_MINOS_antineutrino_vmu = []
flux_MINOS_antineutrino_vmubar = []

# MINOS+ #
energy_MINOSPlus = []
flux_MINOSPlus_neutrino_vmu = []
flux_MINOSPlus_neutrino_vmubar = []

flux_MINOSPlus_antineutrino_vmu = []
flux_MINOSPlus_antineutrino_vmubar = []

# T2K - INGRID #
energy_INGRID = []
flux_INGRID_neutrino_vmu = []
flux_INGRID_neutrino_vmubar = []
flux_INGRID_neutrino_ve = []
flux_INGRID_neutrino_vebar = []

# SHiP #
energy_SHiP = []
flux_SHiP_vmu = []
flux_SHiP_vmubar = []
flux_SHiP_ve = []
flux_SHiP_vebar = []

# SBND #
energy_SBND_vmu = []
energy_SBND_vmubar = []
energy_SBND_ve = []
energy_SBND_vebar = []
flux_SBND_vmu = []
flux_SBND_vmubar = []
flux_SBND_ve = []
flux_SBND_vebar = []

# Get bin edges. This is based on a visual inspection of the DUNE plot.
with open(FLUX_DIR + '/FASERnu/FASERvmu_bins.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        bins_FASERvmu.append(float(row[0]))

bins_FASERvmu.append(8150.00)
bins_FASERvmu = np.array(bins_FASERvmu)


### Cross Sections ###
## Filenames ##
# incoming vmu / vmubar #
xsec_1tau_coh_Ar_filename = CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/coherent/argon/vmu_to_vtau_tau+_mu-_coh_Ar_xsec.csv'
xsec_1tau_coh_W_filename = CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/coherent/tungsten/vmu_to_vtau_tau+_mu-_coh_W_xsec.csv'
xsec_1tau_coh_Fe_filename = CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/coherent/iron/vmu_to_vtau_tau+_mu-_coh_Fe_xsec.csv'
xsec_1tau_p_filename = CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv'
xsec_1tau_n_filename = CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/neutron/vmu_to_vtau_tau+_mu-_nucleon_n_xsec.csv'

xsec_2tau_coh_Ar_filename = CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/coherent/argon/vmu_to_vmu_tau+_tau-_coh_Ar_xsec.csv'
xsec_2tau_coh_W_filename = CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/coherent/tungsten/vmu_to_vmu_tau+_tau-_coh_W_xsec.csv'
xsec_2tau_p_filename = CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv'
xsec_2tau_n_filename = CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/neutron/vmu_to_vmu_tau+_tau-_nucleon_n_xsec.csv'

xsec_2mu_coh_Ar_filename  = CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/argon/vmu_to_vmu_mu+_mu-_coh_Ar_xsec.csv'
xsec_2mu_coh_W_filename = CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/tungsten/vmu_to_vmu_mu+_mu-_coh_W_xsec.csv'
xsec_2mu_coh_Fe_filename = CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/iron/vmu_to_vmu_mu+_mu-_coh_Fe_xsec.csv'
xsec_2mu_p_filename  = CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/proton/vmu_to_vmu_mu+_mu-_nucleon_p_xsec.csv'
xsec_2mu_n_filename  = CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/neutron/vmu_to_vmu_mu+_mu-_nucleon_n_xsec.csv'

xsec_1e1mu_coh_Ar_filename = CROSS_SECTION_DIR + '/vmu_to_ve_e+_mu-_xsec/coherent/argon/vmu_to_ve_e+_mu-_coh_Ar_xsec.csv'
xsec_1e1mu_p_filename = CROSS_SECTION_DIR + '/vmu_to_ve_e+_mu-_xsec/nucleon/proton/vmu_to_ve_e+_mu-_nucleon_p_xsec.csv'
xsec_1e1mu_n_filename = CROSS_SECTION_DIR + '/vmu_to_ve_e+_mu-_xsec/nucleon/neutron/vmu_to_ve_e+_mu-_nucleon_n_xsec.csv'

xsec_vmuCC_filename = CROSS_SECTION_DIR + '/vmuCC/vmuCC_xsec_perE_Formaggio.csv'

xsec_DIS_vmuCC_filename = CROSS_SECTION_DIR + '/DIS/vmu_DIS_CC_tungsten_xsec_perE.csv'
xsec_DIS_vmubarCC_filename = CROSS_SECTION_DIR + '/DIS/vmubar_DIS_CC_tungsten_xsec_perE.csv'
xsec_DIS_vtauCC_filename = CROSS_SECTION_DIR + '/DIS/vtau_DIS_CC_tungsten_xsec_perE.csv'
xsec_DIS_vtaubarCC_filename = CROSS_SECTION_DIR + '/DIS/vtaubar_DIS_CC_tungsten_xsec_perE.csv'

# incoming ve / vebar #
xsec_ve1tau_coh_Ar_filename = CROSS_SECTION_DIR + '/ve_to_vtau_tau+_e-_xsec/coherent/argon/ve_to_vtau_tau+_e-_coh_Ar_xsec.csv'
xsec_ve1tau_p_filename = CROSS_SECTION_DIR + '/ve_to_vtau_tau+_e-_xsec/nucleon/proton/ve_to_vtau_tau+_e-_nucleon_p_xsec.csv'
xsec_ve1tau_n_filename = CROSS_SECTION_DIR + '/ve_to_vtau_tau+_e-_xsec/nucleon/neutron/ve_to_vtau_tau+_e-_nucleon_n_xsec.csv'

# Coherent ; Argon #
energy_1tau_coh_Ar = []
xsec_1tau_coh_Ar = []
delta_1tau_coh_Ar = []

energy_2tau_coh_Ar = []
xsec_2tau_coh_Ar = []
delta_2tau_coh_Ar = []

energy_2mu_coh_Ar = []
xsec_2mu_coh_Ar = []
delta_2mu_coh_Ar = []

energy_ve1tau_coh_Ar = []
xsec_ve1tau_coh_Ar = []
delta_ve1tau_coh_Ar = []

energy_1e1mu_coh_Ar = []
xsec_1e1mu_coh_Ar = []
delta_1e1mu_coh_Ar = []

# Coherent ; Tungsten #
energy_1tau_coh_W = []
xsec_1tau_coh_W = []
delta_1tau_coh_W = []

energy_2tau_coh_W = []
xsec_2tau_coh_W = []
delta_2tau_coh_W = []

energy_2mu_coh_W = []
xsec_2mu_coh_W = []
delta_2mu_coh_W = []

# Coherent ; Iron #
energy_1tau_coh_Fe = []
xsec_1tau_coh_Fe = []
delta_1tau_coh_Fe = []

energy_2mu_coh_Fe = []
xsec_2mu_coh_Fe = []
delta_2mu_coh_Fe = []

# Incoherent ; proton ; Argon #
energy_1tau_p_Ar = []
xsec_1tau_p_Ar = []
delta_1tau_p_Ar = []

energy_2tau_p_Ar = []
xsec_2tau_p_Ar = []
delta_2tau_p_Ar = []

energy_2mu_p_Ar = []
xsec_2mu_p_Ar = []
delta_2mu_p_Ar = []

energy_ve1tau_p_Ar = []
xsec_ve1tau_p_Ar = []
delta_ve1tau_p_Ar = []

energy_1e1mu_p_Ar = []
xsec_1e1mu_p_Ar = []
delta_1e1mu_p_Ar = []

# Incoherent ; neutron ; Argon #
energy_1tau_n_Ar = []
xsec_1tau_n_Ar = []
delta_1tau_n_Ar = []

energy_2tau_n_Ar = []
xsec_2tau_n_Ar = []
delta_2tau_n_Ar = []

energy_2mu_n_Ar = []
xsec_2mu_n_Ar = []
delta_2mu_n_Ar = []

energy_ve1tau_n_Ar = []
xsec_ve1tau_n_Ar = []
delta_ve1tau_n_Ar = []

energy_1e1mu_n_Ar = []
xsec_1e1mu_n_Ar = []
delta_1e1mu_n_Ar = []

# Incoherent ; proton ; Tungsten #
energy_1tau_p_W = []
xsec_1tau_p_W = []
delta_1tau_p_W = []

energy_2tau_p_W = []
xsec_2tau_p_W = []
delta_2tau_p_W = []

energy_2mu_p_W = []
xsec_2mu_p_W = []
delta_2mu_p_W = []

# Incoherent ; neutron ; Tungsten #
energy_1tau_n_W = []
xsec_1tau_n_W = []
delta_1tau_n_W = []

energy_2tau_n_W = []
xsec_2tau_n_W = []
delta_2tau_n_W = []

energy_2mu_n_W = []
xsec_2mu_n_W = []
delta_2mu_n_W = []

# Incoherent ; proton ; Iron #
energy_1tau_p_Fe = []
xsec_1tau_p_Fe = []
delta_1tau_p_Fe = []

energy_2tau_p_Fe = []
xsec_2tau_p_Fe = []
delta_2tau_p_Fe = []

energy_2mu_p_Fe = []
xsec_2mu_p_Fe = []
delta_2mu_p_Fe = []

# Incoherent ; neutron ; Iron #
energy_1tau_n_Fe = []
xsec_1tau_n_Fe = []
delta_1tau_n_Fe = []

energy_2tau_n_Fe = []
xsec_2tau_n_Fe = []
delta_2tau_n_Fe = []

energy_2mu_n_Fe = []
xsec_2mu_n_Fe = []
delta_2mu_n_Fe = []

# vmuCC Formaggion and Zeller #
energy_vmuCC = []
xsec_vmuCC = []

# DIS vmuCC #
energy_DIS_vmuCC_tungsten = []
xsec_DIS_vmuCC_tungsten = []

# DIS vmubarCC #
energy_DIS_vmubarCC_tungsten = []
xsec_DIS_vmubarCC_tungsten = []

# DIS vtauCC #
energy_DIS_vtauCC_tungsten = []
xsec_DIS_vtauCC_tungsten = []

# DIS vtaubarCC #
energy_DIS_vtaubarCC_tungsten = []
xsec_DIS_vtaubarCC_tungsten = []

########################
###### Load Files ######
########################

##### Fluxes #####
### DUNE ; Standard Mode Flux ; Neutrino Mode ###
with open(DUNE_neutrino_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ' ')
    for row in data:
        energy_DUNE.append(float(row[0]))
        flux_DUNE_neutrino_vmu.append(float(row[2])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        flux_DUNE_neutrino_ve.append(float(row[1])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        flux_DUNE_neutrino_vmubar.append(float(row[5])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        flux_DUNE_neutrino_vebar.append(float(row[4])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]

### DUNE ; Standard Mode Flux ; Antineutrino Mode ###
with open(DUNE_antineutrino_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ' ')
    for row in data:
        flux_DUNE_antineutrino_vmu.append(float(row[2])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        flux_DUNE_antineutrino_ve.append(float(row[1])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        flux_DUNE_antineutrino_vmubar.append(float(row[5])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        flux_DUNE_antineutrino_vebar.append(float(row[4])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]

### DUNE ; Tau-optimized Flux ; Neutrino Mode ###
with open(DUNE_tau_opt_neutrino_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ' ')
    for row in data:
        energy_DUNE_tau_opt.append(float(row[0]))
        flux_DUNE_tau_opt_neutrino_vmu.append(float(row[2])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        flux_DUNE_tau_opt_neutrino_ve.append(float(row[1])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        flux_DUNE_tau_opt_neutrino_vmubar.append(float(row[5])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        flux_DUNE_tau_opt_neutrino_vebar.append(float(row[4])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]

### DUNE ; Tau-optimized Flux ; Antineutrino Mode ###
with open(DUNE_tau_opt_antineutrino_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ' ')
    for row in data:
        flux_DUNE_tau_opt_antineutrino_vmu.append(float(row[2])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        flux_DUNE_tau_opt_antineutrino_ve.append(float(row[1])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        flux_DUNE_tau_opt_antineutrino_vmubar.append(float(row[5])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        flux_DUNE_tau_opt_antineutrino_vebar.append(float(row[4])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]

### FASERvmu ; vmu flux ###
with open(FASERnu_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        energy = float(row[0]) # Energy in GeV
        flux = float(row[1]) # Flux has already been normalized by the 25 x 25 cm^2 area and 150 fb^-1; it is now in units of [m^-2 GeV^-1 fb].
        energy_FASERvmu.append(energy)
        flux_FASERvmu.append(flux)

### FASERvmu ; vmubar flux ###
with open(FASERnubar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        flux = float(row[1]) # Flux has already been normalized by the 25 x 25 cm^2 area and 150 fb^-1; it is now in units of [m^-2 GeV^-1 fb].
        flux_FASERvmubar.append(flux)

### FASERv2 ; vmu + vmubar flux ###
with open(FASERnu2_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        energy = float(row[0]) # Energy in GeV
        flux = float(row[1]) # Flux has not been normalized by the 0.5 x 0.5 m^2 area and 3 ab^-1; it is in units of [1].
        energy_FASERv2.append(energy)
        flux_FASERv2_vmu_vmubar.append(flux / (0.25 * 3000))

### MINOS ; Neutrino Mode ; vmu flux ###
with open(MINOS_neutrino_vmu_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        energy = float(row[0]) # Energy in GeV
        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
        energy_MINOS.append(energy)
        flux_MINOS_neutrino_vmu.append(flux)

### MINOS ; Neutrino Mode ; vmubar flux ###
with open(MINOS_neutrino_vmubar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
        flux_MINOS_neutrino_vmubar.append(flux)

### MINOS ; Antineutrino Mode ; vmu flux ###
with open(MINOS_antineutrino_vmu_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        energy = float(row[0]) # Energy in GeV
        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
        flux_MINOS_antineutrino_vmu.append(flux)

### MINOS ; Antineutrino Mode ; vmubar flux ###
with open(MINOS_antineutrino_vmubar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
        flux_MINOS_antineutrino_vmubar.append(flux)

### MINOS+ ; Neutrino Mode ; vmu flux ###
with open(MINOSPlus_neutrino_vmu_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        energy = float(row[0]) # Energy in GeV
        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
        energy_MINOSPlus.append(energy)
        flux_MINOSPlus_neutrino_vmu.append(flux)

### MINOS+ ; Neutrino Mode ; vmubar flux ###
with open(MINOSPlus_neutrino_vmubar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
        flux_MINOSPlus_neutrino_vmubar.append(flux)

### MINOS+ ; Antineutrino Mode ; vmu flux ###
with open(MINOSPlus_antineutrino_vmu_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
        flux_MINOSPlus_antineutrino_vmu.append(flux)

### MINOS+ ; Antineutrino Mode ; vmubar flux ###
with open(MINOSPlus_antineutrino_vmubar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
        flux_MINOSPlus_antineutrino_vmubar.append(flux)

### T2K ; Neutrino Mode ; flux ###
with open(T2K_INGRID_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        energy = float(row[0])
        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1].
        energy_INGRID.append(energy)
        flux_INGRID_neutrino_vmu.append(flux * 92.5/100.0)  # Flux composition of vmu is 92.5%. See Ballett et al.
        flux_INGRID_neutrino_vmubar.append(flux * 5.8/100.0) # Flux composition of vmu is 5.8%. See Ballett et al.
        flux_INGRID_neutrino_ve.append(flux * 1.5/100.0) # Flux composition of vmu is 1.5%. See Ballett et al.
        flux_INGRID_neutrino_vebar.append(flux * 0.2/100.0) # Flux composition of vmu is 0.2%. See Ballett et al.

### SHiP ###
with open(SHiP_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ',')
    for row in data:
        energy_SHiP.append(float(row[0]))
        flux_SHiP_vmu.append(float(row[1])) # Histogram has units of [m^-2 POT^-1]
        flux_SHiP_ve.append(float(row[3])) # Histogram has units of [m^-2 POT^-1]
        flux_SHiP_vmubar.append(float(row[2])) # Histogram has units of [m^-2 POT^-1]
        flux_SHiP_vebar.append(float(row[4])) # Histogram has units of [m^-2 POT^-1]

### SBND ; vmu flux ###
with open(SBND_vmu_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ',')
    for row in data:
        energy_SBND_vmu.append(float(row[0]))
        flux_SBND_vmu.append(float(row[1]))

### SBND ; vmubar flux ###
with open(SBND_vmubar_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ',')
    for row in data:
        energy_SBND_vmubar.append(float(row[0]))
        flux_SBND_vmubar.append(float(row[1]))

### SBND ; ve flux ###
with open(SBND_ve_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ',')
    for row in data:
        energy_SBND_ve.append(float(row[0]))
        flux_SBND_ve.append(float(row[1]))

### SBND ; vebar flux ###
with open(SBND_vebar_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ',')
    for row in data:
        energy_SBND_vebar.append(float(row[0]))
        flux_SBND_vebar.append(float(row[1]))

### Integrated Fluxes ###
# DUNE Standard #
DUNE_neutrino_vmu_integrated_flux = simpson(flux_DUNE_neutrino_vmu, x=energy_DUNE)                                                 # Integrated flux [Nv / m^2 POT]
DUNE_neutrino_vmubar_integrated_flux = simpson(flux_DUNE_neutrino_vmubar, x=energy_DUNE)                                           # Integrated flux [Nv / m^2 POT]
DUNE_neutrino_ve_integrated_flux = simpson(flux_DUNE_neutrino_ve, x=energy_DUNE)                                                   # Integrated flux [Nv / m^2 POT]
DUNE_neutrino_vebar_integrated_flux = simpson(flux_DUNE_neutrino_vebar, x=energy_DUNE)                                             # Integrated flux [Nv / m^2 POT]

DUNE_antineutrino_vmu_integrated_flux = simpson(flux_DUNE_antineutrino_vmu, x=energy_DUNE)                                         # Integrated flux [Nv / m^2 POT]
DUNE_antineutrino_vmubar_integrated_flux = simpson(flux_DUNE_antineutrino_vmubar, x=energy_DUNE)                                   # Integrated flux [Nv / m^2 POT]
DUNE_antineutrino_ve_integrated_flux = simpson(flux_DUNE_antineutrino_ve, x=energy_DUNE)                                           # Integrated flux [Nv / m^2 POT]
DUNE_antineutrino_vebar_integrated_flux = simpson(flux_DUNE_antineutrino_vebar, x=energy_DUNE)                                     # Integrated flux [Nv / m^2 POT]

# DUNE Tau Optimized #
DUNE_tau_opt_neutrino_vmu_integrated_flux = simpson(flux_DUNE_tau_opt_neutrino_vmu, x=energy_DUNE_tau_opt)                         # Integrated flux [Nv / m^2 POT]
DUNE_tau_opt_neutrino_vmubar_integrated_flux = simpson(flux_DUNE_tau_opt_neutrino_vmubar, x=energy_DUNE_tau_opt)                   # Integrated flux [Nv / m^2 POT]
DUNE_tau_opt_neutrino_ve_integrated_flux = simpson(flux_DUNE_tau_opt_neutrino_ve, x=energy_DUNE_tau_opt)                           # Integrated flux [Nv / m^2 POT]
DUNE_tau_opt_neutrino_vebar_integrated_flux = simpson(flux_DUNE_tau_opt_neutrino_vebar, x=energy_DUNE_tau_opt)                     # Integrated flux [Nv / m^2 POT]

DUNE_tau_opt_antineutrino_vmu_integrated_flux = simpson(flux_DUNE_tau_opt_antineutrino_vmu, x=energy_DUNE_tau_opt)                 # Integrated flux [Nv / m^2 POT]
DUNE_tau_opt_antineutrino_vmubar_integrated_flux = simpson(flux_DUNE_tau_opt_antineutrino_vmubar, x=energy_DUNE_tau_opt)           # Integrated flux [Nv / m^2 POT]
DUNE_tau_opt_antineutrino_ve_integrated_flux = simpson(flux_DUNE_tau_opt_antineutrino_ve, x=energy_DUNE_tau_opt)                   # Integrated flux [Nv / m^2 POT]
DUNE_tau_opt_antineutrino_vebar_integrated_flux = simpson(flux_DUNE_tau_opt_antineutrino_vebar, x=energy_DUNE_tau_opt)             # Integrated flux [Nv / m^2 POT]

# FASER #
#FASER_integrated_flux = simpson(np.ones(len(flux_FASERvmu)), x=list(reversed(flux_FASERvmu)))       # Integrated flux [Nv / m^2 fb^-1]
FASER_integrated_flux = sum(flux_FASERvmu)
FASER_vmubar_integrated_flux = sum(flux_FASERvmubar)
FASERv2_integrated_flux = sum(flux_FASERv2_vmu_vmubar)

# MINOS #
MINOS_neutrino_vmu_integrated_flux = simpson(flux_MINOS_neutrino_vmu, x=energy_MINOS)/1e20                                         # Integrated flux [Nv / m^2 POT]
MINOS_neutrino_vmubar_integrated_flux = simpson(flux_MINOS_neutrino_vmubar, x=energy_MINOS)/1e20                                   # Integrated flux [Nv / m^2 POT]
MINOS_antineutrino_vmu_integrated_flux = simpson(flux_MINOS_antineutrino_vmu, x=energy_MINOS)/1e20                                 # Integrated flux [Nv / m^2 POT]
MINOS_antineutrino_vmubar_integrated_flux = simpson(flux_MINOS_antineutrino_vmubar, x=energy_MINOS)/1e20                           # Integrated flux [Nv / m^2 POT]

# MINOS+ #
MINOSPlus_neutrino_vmu_integrated_flux = simpson(flux_MINOSPlus_neutrino_vmu, x=energy_MINOSPlus)/1e20                             # Integrated flux [Nv / m^2 POT]
MINOSPlus_neutrino_vmubar_integrated_flux = simpson(flux_MINOSPlus_neutrino_vmubar, x=energy_MINOSPlus)/1e20                       # Integrated flux [Nv / m^2 POT]

# T2K - INGRID #
INGRID_neutrino_vmu_integrated_flux = simpson(flux_INGRID_neutrino_vmu, x=energy_INGRID)/1e20                                      # Integrated flux [Nv / m^2 POT]
INGRID_neutrino_vmubar_integrated_flux = simpson(flux_INGRID_neutrino_vmubar, x=energy_INGRID)/1e20                                # Integrated flux [Nv / m^2 POT]
INGRID_neutrino_ve_integrated_flux = simpson(flux_INGRID_neutrino_ve, x=energy_INGRID)/1e20                                        # Integrated flux [Nv / m^2 POT]
INGRID_neutrino_vebar_integrated_flux = simpson(flux_INGRID_neutrino_vebar, x=energy_INGRID)/1e20                                  # Integrated flux [Nv / m^2 POT]

# SHiP #
SHiP_vmu_integrated_flux = sum(flux_SHiP_vmu)
SHiP_vmubar_integrated_flux = sum(flux_SHiP_vmubar)
SHiP_ve_integrated_flux = sum(flux_SHiP_ve)
SHiP_vebar_integrated_flux = sum(flux_SHiP_vebar)

# SBND #
SBND_vmu_integrated_flux = simpson(flux_SBND_vmu, x=energy_SBND_vmu)/1e20                                                          # Integrated flux [Nv / m^2 POT]
SBND_vmubar_integrated_flux = simpson(flux_SBND_vmu, x=energy_SBND_vmu)/1e20                                                       # Integrated flux [Nv / m^2 POT]
SBND_ve_integrated_flux = simpson(flux_SBND_ve, x=energy_SBND_ve)/1e20                                                             # Integrated flux [Nv / m^2 POT]
SBND_vebar_integrated_flux = simpson(flux_SBND_vebar, x=energy_SBND_vebar)/1e20                                                    # Integrated flux [Nv / m^2 POT]

##### Cross Sections #####
### Coherent ###
# vmu -> vtau tau+ mu- ; coherent ; Argon #
with open(xsec_1tau_coh_Ar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_coh_Ar.append(float(row[0]))
        xsec_1tau_coh_Ar.append(float(row[1]) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) * 1e-43
        delta_1tau_coh_Ar.append(delta)

# vmu -> vmu tau+ tau- ; coherent ; Argon #
with open(xsec_2tau_coh_Ar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_coh_Ar.append(float(row[0]))
        xsec_2tau_coh_Ar.append(float(row[1]) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) * 1e-43
        delta_2tau_coh_Ar.append(delta)

# vmu -> vmu mu+ mu- ; coherent ; Argon #
with open(xsec_2mu_coh_Ar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_coh_Ar.append(float(row[0]))
        xsec_2mu_coh_Ar.append(float(row[1]) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) * 1e-43
        delta_2mu_coh_Ar.append(delta)

# ve -> vtau tau+ e- ; coherent ; Argon #
with open(xsec_ve1tau_coh_Ar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_ve1tau_coh_Ar.append(float(row[0]))
        xsec_ve1tau_coh_Ar.append(float(row[1]) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) * 1e-43
        delta_ve1tau_coh_Ar.append(delta)

# vmu -> ve e+ mu- ; coherent ; Argon #
with open(xsec_1e1mu_coh_Ar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1e1mu_coh_Ar.append(float(row[0]))
        xsec_1e1mu_coh_Ar.append(float(row[1]) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) * 1e-43
        delta_1e1mu_coh_Ar.append(delta)

# vmu -> vtau tau+ mu- ; coherent ; Tungsten #
with open(xsec_1tau_coh_W_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_coh_W.append(float(row[0]))
        xsec_1tau_coh_W.append(float(row[1]) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_W * float(row[1])
        delta = (numerical_delta + physical_delta) * 1e-43
        delta_1tau_coh_W.append(delta)

# vmu -> vmu tau+ tau- ; coherent ; Tungsten #
with open(xsec_2tau_coh_W_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_coh_W.append(float(row[0]))
        xsec_2tau_coh_W.append(float(row[1]) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_W * float(row[1])
        delta = (numerical_delta + physical_delta) * 1e-43
        delta_2tau_coh_W.append(delta)

# vmu -> vmu mu+ mu- ; coherent ; Tungsten #
with open(xsec_2mu_coh_W_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_coh_W.append(float(row[0]))
        xsec_2mu_coh_W.append(float(row[1]) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_W * float(row[1])
        delta = (numerical_delta + physical_delta) * 1e-43
        delta_2mu_coh_W.append(delta)

# vmu -> vtau tau+ mu- ; coherent ; Iron #
with open(xsec_1tau_coh_Fe_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_coh_Fe.append(float(row[0]))
        xsec_1tau_coh_Fe.append(float(row[1]) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) * 1e-43
        delta_1tau_coh_Fe.append(delta)

# vmu -> vmu mu+ mu- ; coherent ; Tungsten #
with open(xsec_2mu_coh_Fe_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_coh_Fe.append(float(row[0]))
        xsec_2mu_coh_Fe.append(float(row[1]) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_W * float(row[1])
        delta = (numerical_delta + physical_delta) * 1e-43
        delta_2mu_coh_Fe.append(delta)

### Incoherent ; proton ; Argon ###
# vmu -> vtau tau+ mu- #
with open(xsec_1tau_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_p_Ar.append(float(row[0]))
        xsec_1tau_p_Ar.append(float(row[1]) * Z_Ar * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) * Z_Ar * 1e-43
        delta_1tau_p_Ar.append(delta)

# vmu -> vmu tau+ tau- #
with open(xsec_2tau_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_p_Ar.append(float(row[0]))
        xsec_2tau_p_Ar.append(float(row[1]) * Z_Ar * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) * Z_Ar * 1e-43
        delta_2tau_p_Ar.append(delta)

# vmu -> vmu mu+ mu- #
with open(xsec_2mu_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_p_Ar.append(float(row[0]))
        xsec_2mu_p_Ar.append(float(row[1]) * Z_Ar * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) * Z_Ar * 1e-43
        delta_2mu_p_Ar.append(delta)

# ve -> vtau tau+ e- #
with open(xsec_ve1tau_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_ve1tau_p_Ar.append(float(row[0]))
        xsec_ve1tau_p_Ar.append(float(row[1]) * Z_Ar * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) * Z_Ar * 1e-43
        delta_ve1tau_p_Ar.append(delta)

# vmu -> ve e+ mu- #
with open(xsec_1e1mu_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1e1mu_p_Ar.append(float(row[0]))
        xsec_1e1mu_p_Ar.append(float(row[1]) * Z_Ar * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) * Z_Ar * 1e-43
        delta_1e1mu_p_Ar.append(delta)

### Incoherent ; proton ; Tungsten ###
# vmu -> vtau tau+ mu- #
with open(xsec_1tau_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_p_W.append(float(row[0]))
        xsec_1tau_p_W.append(float(row[1]) * Z_W * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_W * float(row[1])
        delta = (numerical_delta + physical_delta) * Z_W * 1e-43
        delta_1tau_p_W.append(delta)

# vmu -> vmu tau+ tau- #
with open(xsec_2tau_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_p_W.append(float(row[0]))
        xsec_2tau_p_W.append(float(row[1]) * Z_W * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_W * float(row[1])
        delta = (numerical_delta + physical_delta) * Z_W * 1e-43
        delta_2tau_p_W.append(delta)

# vmu -> vmu mu+ mu- #
with open(xsec_2mu_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_p_W.append(float(row[0]))
        xsec_2mu_p_W.append(float(row[1]) * Z_W * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_W * float(row[1])
        delta = (numerical_delta + physical_delta) * Z_W * 1e-43
        delta_2mu_p_W.append(delta)

### Incoherent ; proton ; Iron ###
# vmu -> vtau tau+ mu- #
with open(xsec_1tau_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_p_Fe.append(float(row[0]))
        xsec_1tau_p_Fe.append(float(row[1]) * Z_Fe * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) * Z_Fe * 1e-43
        delta_1tau_p_Fe.append(delta)

# vmu -> vmu mu+ mu- #
with open(xsec_2mu_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_p_Fe.append(float(row[0]))
        xsec_2mu_p_Fe.append(float(row[1]) * Z_Fe * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) * Z_Fe * 1e-43
        delta_2mu_p_Fe.append(delta)

### Incoherent ; neutron ; Argon ###
# vmu -> vtau tau+ mu- #
with open(xsec_1tau_n_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_n_Ar.append(float(row[0]))
        xsec_1tau_n_Ar.append(float(row[1]) * (A_Ar - Z_Ar) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) * (A_Ar - Z_Ar) * 1e-43
        delta_1tau_n_Ar.append(delta)

# vmu -> vmu tau+ tau- #
with open(xsec_2tau_n_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_n_Ar.append(float(row[0]))
        xsec_2tau_n_Ar.append(float(row[1]) * (A_Ar - Z_Ar) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) * (A_Ar - Z_Ar) * 1e-43
        delta_2tau_n_Ar.append(delta)

# vmu -> vmu mu+ mu- #
with open(xsec_2mu_n_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_n_Ar.append(float(row[0]))
        xsec_2mu_n_Ar.append(float(row[1]) * (A_Ar - Z_Ar) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) * (A_Ar - Z_Ar) * 1e-43
        delta_2mu_n_Ar.append(delta)

# ve -> vtau tau+ e- #
with open(xsec_ve1tau_n_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_ve1tau_n_Ar.append(float(row[0]))
        xsec_ve1tau_n_Ar.append(float(row[1]) * (A_Ar - Z_Ar) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) * (A_Ar - Z_Ar) * 1e-43
        delta_ve1tau_n_Ar.append(delta)

# vmu -> ve e+ mu- #
with open(xsec_1e1mu_n_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1e1mu_n_Ar.append(float(row[0]))
        xsec_1e1mu_n_Ar.append(float(row[1]) * (A_Ar - Z_Ar) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) * (A_Ar - Z_Ar) * 1e-43
        delta_1e1mu_n_Ar.append(delta)

### Incoherent ; neutron ; Tungsten ###
# vmu -> vtau tau+ mu- #
with open(xsec_1tau_n_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_n_W.append(float(row[0]))
        xsec_1tau_n_W.append(float(row[1]) * (A_W - Z_W) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_W * float(row[1])
        delta = (numerical_delta + physical_delta) * (A_W - Z_W) * 1e-43
        delta_1tau_n_W.append(delta)

# vmu -> vmu tau+ tau- #
with open(xsec_2tau_n_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_n_W.append(float(row[0]))
        xsec_2tau_n_W.append(float(row[1]) * (A_W - Z_W) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_W * float(row[1])
        delta = (numerical_delta + physical_delta) * (A_W - Z_W) * 1e-43
        delta_2tau_n_W.append(delta)

# vmu -> vmu mu+ mu- #
with open(xsec_2mu_n_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_n_W.append(float(row[0]))
        xsec_2mu_n_W.append(float(row[1]) * (A_W - Z_W) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_W * float(row[1])
        delta = (numerical_delta + physical_delta) * (A_W - Z_W) * 1e-43
        delta_2mu_n_W.append(delta)

### Incoherent ; neutron ; Iron ###
# vmu -> vtau tau+ mu- #
with open(xsec_1tau_n_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_n_Fe.append(float(row[0]))
        xsec_1tau_n_Fe.append(float(row[1]) * (A_Fe - Z_Fe) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) * (A_Fe - Z_Fe) * 1e-43
        delta_1tau_n_Fe.append(delta)

# vmu -> vmu mu+ mu- #
with open(xsec_2mu_n_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_n_Fe.append(float(row[0]))
        xsec_2mu_n_Fe.append(float(row[1]) * (A_Fe - Z_Fe) * 1e-43) # Convert to m^2.
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_W * float(row[1])
        delta = (numerical_delta + physical_delta) * (A_Fe - Z_Fe) * 1e-43
        delta_2mu_n_Fe.append(delta)

### vmu X -> mu- X' ; vmuCC ; Formaggio & Zeller ###
# This vmuCC cross section applies to Tungsten only. #
with open(xsec_vmuCC_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) * energy * 1e-42 * A_W # digitized plot is in [1e-38 cm^2 / GeV] convert to [m^2]; also, xsec given in as xsec/nucleon -> xsec/W

        energy_vmuCC.append(energy)
        xsec_vmuCC.append(xsec)

### vmu X -> mu- X' ; DIS vmuCC ; Detecting and Studying High-Energy Collider Neutrinos with FASER at the LHC ###
# This vmuCC cross section applies to Tungsten only. #
with open(xsec_DIS_vmuCC_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) * energy * 1e-4 # digitized plot is in [cm^2 / GeV] convert to [m^2]

        energy_DIS_vmuCC_tungsten.append(energy)
        xsec_DIS_vmuCC_tungsten.append(xsec)

### vmubar X -> mu+ X' ; DIS vmuCC ; Detecting and Studying High-Energy Collider Neutrinos with FASER at the LHC ###
# This vmubarCC cross section applies to Tungsten only. #
with open(xsec_DIS_vmubarCC_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) * energy * 1e-4 # digitized plot is in [cm^2 / GeV] convert to [m^2]

        energy_DIS_vmubarCC_tungsten.append(energy)
        xsec_DIS_vmubarCC_tungsten.append(xsec)

### vtau X -> tau- X' ; DIS vtauCC ; Detecting and Studying High-Energy Collider Neutrinos with FASER at the LHC ###
# This vtauCC cross section applies to Tungsten only. #
with open(xsec_DIS_vtauCC_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) * energy * 1e-4 # digitized plot is in [cm^2 / GeV] convert to [m^2]

        energy_DIS_vtauCC_tungsten.append(energy)
        xsec_DIS_vtauCC_tungsten.append(xsec)

### vtaubar X -> tau+ X' ; DIS vtaubarCC ; Detecting and Studying High-Energy Collider Neutrinos with FASER at the LHC ###
# This vtaubarCC cross section applies to Tungsten only. #
with open(xsec_DIS_vtaubarCC_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) * energy * 1e-4 # digitized plot is in [cm^2 / GeV] convert to [m^2]

        energy_DIS_vtaubarCC_tungsten.append(energy)
        xsec_DIS_vtaubarCC_tungsten.append(xsec)

####################
##### EPA xsec #####
####################

LambdaQCD = 217e-3 # LambdaQCD is about 217 MeV for 5 quark flavors.
inverseGeV2_to_cm2 = (0.197e-13)**2 # 1 GeV^-1 = 0.197e-15 m

def EPA_XSec_Coh_Heaviside_FormFactor(E, m1, m2, axial, vector, target='argon'):
    if target == 'argon':
        A = A_Ar
        Z = Z_Ar
    if target == 'iron':
        A = A_Fe
        Z = Z_Fe
    if target == 'tungsten':
        A = A_W
        Z = Z_W

    Qmax = LambdaQCD / A**(1./3.)
    smax = 2.*E*Qmax

    sigma = (1./2.)*(axial**2 + vector**2)*(2.*(Z*aEM*GFermi)**2)/(9.*np.pi**3)*smax*np.log(smax/((m1+m2)**2)) # xsec in natural units [1/GeV^2]
    if sigma > 0:
        return float(sigma*inverseGeV2_to_cm2 * 1e-4) # xsec in m^2 now.
    else:
        return float(0)

EPA_Heaviside_xsec_1tau_coh_Ar = []
for E in energy_1tau_coh_Ar:
    xsec = EPA_XSec_Coh_Heaviside_FormFactor(E,mmu,mtau,1.,1.,target='argon')
    EPA_Heaviside_xsec_1tau_coh_Ar.append(xsec)

EPA_Heaviside_xsec_2mu_coh_Ar = []
for E in energy_2mu_coh_Ar:
    xsec = EPA_XSec_Coh_Heaviside_FormFactor(E,mmu,mmu,1./2.,1./2.+2*sW2,target='argon')
    EPA_Heaviside_xsec_2mu_coh_Ar.append(xsec)

###################################
##### RMiss Correction Factor #####
###################################

correction_factor_RMiss_hadronic_filename = CORRECTION_FACTOR_DIR+'/RMiss_correction_factors_hadronic.txt'

# Correction factors for RMiss start at 5 GeV. At 2 GeV the threshold prevents any events from happening. 
# To account for any (unlikely) events showing up below 5 GeV and above the threshold, I am setting eps = 1.0 for Ev = 0 GeV
# which follows the trend by the other eps at different Ev.
correction_factor_RMiss_energy = [0.0]
correction_factor_RMiss_p = [1.0]
correction_factor_RMiss_n = [1.0]
correction_factor_RMiss_coh = [1.0]

with open(correction_factor_RMiss_hadronic_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    first_line = True
    for row in data:
        if first_line:
            first_line = False
            continue
        energy = float(row[0])
        correction_factor_p = float(row[1])
        correction_factor_n = float(row[2])
        correction_factor_coh = float(row[3])
        correction_factor_RMiss_energy.append(energy)
        correction_factor_RMiss_p.append(correction_factor_p)
        correction_factor_RMiss_n.append(correction_factor_n)
        correction_factor_RMiss_coh.append(correction_factor_coh)

####################################
###### Cross Section Matching ######
####################################

def Interpolation(E1, E2, xsec1, xsec2, E):
    """
    Calculate the xsec for a given energy E given two energy values E1 and E2 with xsec1 and xsec2 respectively. E1 <= E <= E2.
    """
    assert E1 < E2
    assert E >= E1
    assert E <= E2

    slope = (xsec2 - xsec1)/(E2-E1)
    interp = slope*(E-E1) + xsec1
    return interp

def MatchXSec(flux_energy, trid_energy, trid_xsec):
    """
    Interpolate the trident xsec to the flux_energies and return the resulting list.
    """
    flux_size = len(flux_energy)
    trid_size = len(trid_energy)
    xsec_matched = []

    for ii in range(flux_size):
        if flux_energy[ii] > trid_energy[-1]: # If the flux energy is larger than the maximum trident energy, set the matched xsec to 0.
            xsec_matched.append(0.0)
            continue
        elif flux_energy[ii] < trid_energy[0]: # If the flux energy range starts below the trident one, set the matched xsec to 0. Otherwise, these values will be skipped.
            xsec_matched.append(0.0)
            continue
        else:
            for jj in range(trid_size):
                if flux_energy[ii] == trid_energy[jj]: # In the unlikely event that the flux energy equals the trident energy, just output the cross section there.
                    xsec_matched.append(trid_xsec[jj])
                    continue
                if ((flux_energy[ii] > trid_energy[jj]) and (flux_energy[ii] < trid_energy[jj+1])): # Find the trident energy range for the flux energy. Interpolate there.
                    if (trid_xsec[jj] < 0) or (trid_xsec[jj+1] < 0): # Some xsec are negative in the 2tau proton case; ad hoc solution for now.
                        xsec_matched.append(0.0)
                    else:
                        xsec_interp = Interpolation(trid_energy[jj], trid_energy[jj+1], trid_xsec[jj], trid_xsec[jj+1], flux_energy[ii])
                        xsec_matched.append(xsec_interp)
    
    assert len(flux_energy) == len(xsec_matched)
    return xsec_matched

def Limits(xsec, delta):
    "Return arrays of xsec+delta and xsec-delta"
    upper_limit = []
    lower_limit = []
    for i, j in zip(xsec, delta):
        upper_limit.append(i+j)
        lower_limit.append(i-j)
    return lower_limit, upper_limit

def MatchEps(flux_energy, eps_energy, eps):
    """
    Interpolate the trident xsec to the flux_energies and return the resulting list.
    """
    flux_size = len(flux_energy)
    eps_size = len(eps_energy)
    eps_matched = []

    for ii in range(flux_size):
        if flux_energy[ii] > eps_energy[-1]: # If the flux energy is larger than the maximum correction factor energy, set the matched correction factor to 0.
            eps_matched.append(0.0)
            continue
        elif flux_energy[ii] < eps_energy[0]: # If the flux energy range starts below the eps one (0 GeV), set the matched eps to 1. Otherwise, these values will be skipped.
            eps_matched.append(1.0)
            continue
        else:
            for jj in range(eps_size):
                if flux_energy[ii] == eps_energy[jj]: # In the unlikely event that the flux energy equals the correction factor energy, just output the correction factor there.
                    eps_matched.append(eps[jj])
                    continue
                if ((flux_energy[ii] > eps_energy[jj]) and (flux_energy[ii] < eps_energy[jj+1])): # Find the correction factor energy range for the flux energy. Interpolate there.
                    eps_interp = Interpolation(eps_energy[jj], eps_energy[jj+1], eps[jj], eps[jj+1], flux_energy[ii])
                    eps_matched.append(eps_interp)
    
    assert flux_size == len(eps_matched)
    return eps_matched

### Upper and lower cross sections ###
xsec_1tau_coh_Ar_lower, xsec_1tau_coh_Ar_upper = Limits(xsec_1tau_coh_Ar, delta_1tau_coh_Ar)
xsec_1tau_p_Ar_lower, xsec_1tau_p_Ar_upper = Limits(xsec_1tau_p_Ar, delta_1tau_p_Ar)
xsec_1tau_n_Ar_lower, xsec_1tau_n_Ar_upper = Limits(xsec_1tau_n_Ar, delta_1tau_n_Ar)

xsec_2tau_coh_Ar_lower, xsec_2tau_coh_Ar_upper = Limits(xsec_2tau_coh_Ar, delta_2tau_coh_Ar)
xsec_2tau_p_Ar_lower, xsec_2tau_p_Ar_upper = Limits(xsec_2tau_p_Ar, delta_2tau_p_Ar)
xsec_2tau_n_Ar_lower, xsec_2tau_n_Ar_upper = Limits(xsec_2tau_n_Ar, delta_2tau_n_Ar)

xsec_2mu_coh_Ar_lower, xsec_2mu_coh_Ar_upper = Limits(xsec_2mu_coh_Ar, delta_2mu_coh_Ar)
xsec_2mu_p_Ar_lower, xsec_2mu_p_Ar_upper = Limits(xsec_2mu_p_Ar, delta_2mu_p_Ar)
xsec_2mu_n_Ar_lower, xsec_2mu_n_Ar_upper = Limits(xsec_2mu_n_Ar, delta_2mu_n_Ar)

xsec_ve1tau_coh_Ar_lower, xsec_ve1tau_coh_Ar_upper = Limits(xsec_ve1tau_coh_Ar, delta_ve1tau_coh_Ar)
xsec_ve1tau_p_Ar_lower, xsec_ve1tau_p_Ar_upper = Limits(xsec_ve1tau_p_Ar, delta_ve1tau_p_Ar)
xsec_ve1tau_n_Ar_lower, xsec_ve1tau_n_Ar_upper = Limits(xsec_ve1tau_n_Ar, delta_ve1tau_n_Ar)

xsec_1e1mu_coh_Ar_lower, xsec_1e1mu_coh_Ar_upper = Limits(xsec_1e1mu_coh_Ar, delta_1e1mu_coh_Ar)
xsec_1e1mu_p_Ar_lower, xsec_1e1mu_p_Ar_upper = Limits(xsec_1e1mu_p_Ar, delta_1e1mu_p_Ar)
xsec_1e1mu_n_Ar_lower, xsec_1e1mu_n_Ar_upper = Limits(xsec_1e1mu_n_Ar, delta_1e1mu_n_Ar)

xsec_1tau_coh_W_lower, xsec_1tau_coh_W_upper = Limits(xsec_1tau_coh_W, delta_1tau_coh_W)
xsec_1tau_p_W_lower, xsec_1tau_p_W_upper = Limits(xsec_1tau_p_W, delta_1tau_p_W)
xsec_1tau_n_W_lower, xsec_1tau_n_W_upper = Limits(xsec_1tau_n_W, delta_1tau_n_W)

xsec_2tau_coh_W_lower, xsec_2tau_coh_W_upper = Limits(xsec_2tau_coh_W, delta_2tau_coh_W)
xsec_2tau_p_W_lower, xsec_2tau_p_W_upper = Limits(xsec_2tau_p_W, delta_2tau_p_W)
xsec_2tau_n_W_lower, xsec_2tau_n_W_upper = Limits(xsec_2tau_n_W, delta_2tau_n_W)

xsec_2mu_coh_W_lower, xsec_2mu_coh_W_upper = Limits(xsec_2mu_coh_W, delta_2mu_coh_W)
xsec_2mu_p_W_lower, xsec_2mu_p_W_upper = Limits(xsec_2mu_p_W, delta_2mu_p_W)
xsec_2mu_n_W_lower, xsec_2mu_n_W_upper = Limits(xsec_2mu_n_W, delta_2mu_n_W)

xsec_1tau_coh_Fe_lower, xsec_1tau_coh_Fe_upper = Limits(xsec_1tau_coh_Fe, delta_1tau_coh_Fe)
xsec_1tau_p_Fe_lower, xsec_1tau_p_Fe_upper = Limits(xsec_1tau_p_Fe, delta_1tau_p_Fe)
xsec_1tau_n_Fe_lower, xsec_1tau_n_Fe_upper = Limits(xsec_1tau_n_Fe, delta_1tau_n_Fe)

xsec_2mu_coh_Fe_lower, xsec_2mu_coh_Fe_upper = Limits(xsec_2mu_coh_Fe, delta_2mu_coh_Fe)
xsec_2mu_p_Fe_lower, xsec_2mu_p_Fe_upper = Limits(xsec_2mu_p_Fe, delta_2mu_p_Fe)
xsec_2mu_n_Fe_lower, xsec_2mu_n_Fe_upper = Limits(xsec_2mu_n_Fe, delta_2mu_n_Fe)

### DUNE Standard Flux ; Matched XSec ###
DUNE_xsec_1tau_coh_Ar_matched = MatchXSec(energy_DUNE, energy_1tau_coh_Ar, xsec_1tau_coh_Ar)
DUNE_xsec_1tau_p_Ar_matched = MatchXSec(energy_DUNE, energy_1tau_p_Ar, xsec_1tau_p_Ar)
DUNE_xsec_1tau_n_Ar_matched = MatchXSec(energy_DUNE, energy_1tau_n_Ar, xsec_1tau_n_Ar)

DUNE_xsec_2tau_coh_Ar_matched = MatchXSec(energy_DUNE, energy_2tau_coh_Ar, xsec_2tau_coh_Ar)
DUNE_xsec_2tau_p_Ar_matched = MatchXSec(energy_DUNE, energy_2tau_p_Ar, xsec_2tau_p_Ar)
DUNE_xsec_2tau_n_Ar_matched = MatchXSec(energy_DUNE, energy_2tau_n_Ar, xsec_2tau_n_Ar)

DUNE_xsec_2mu_coh_Ar_matched = MatchXSec(energy_DUNE, energy_2mu_coh_Ar, xsec_2mu_coh_Ar)
DUNE_xsec_2mu_p_Ar_matched = MatchXSec(energy_DUNE, energy_2mu_p_Ar, xsec_2mu_p_Ar)
DUNE_xsec_2mu_n_Ar_matched = MatchXSec(energy_DUNE, energy_2mu_n_Ar, xsec_2mu_n_Ar)

DUNE_xsec_ve1tau_coh_Ar_matched = MatchXSec(energy_DUNE, energy_ve1tau_coh_Ar, xsec_ve1tau_coh_Ar)
DUNE_xsec_ve1tau_p_Ar_matched = MatchXSec(energy_DUNE, energy_ve1tau_p_Ar, xsec_ve1tau_p_Ar)
DUNE_xsec_ve1tau_n_Ar_matched = MatchXSec(energy_DUNE, energy_ve1tau_n_Ar, xsec_ve1tau_n_Ar)

# EPA
DUNE_EPA_Heaviside_xsec_1tau_coh_Ar_matched = MatchXSec(energy_DUNE, energy_1tau_coh_Ar, EPA_Heaviside_xsec_1tau_coh_Ar)
DUNE_EPA_Heaviside_xsec_2mu_coh_Ar_matched = MatchXSec(energy_DUNE, energy_2mu_coh_Ar, EPA_Heaviside_xsec_2mu_coh_Ar)


# Upper and Lower limits #
DUNE_xsec_1tau_coh_Ar_matched_upper = MatchXSec(energy_DUNE, energy_1tau_coh_Ar, xsec_1tau_coh_Ar_upper)
DUNE_xsec_1tau_p_Ar_matched_upper = MatchXSec(energy_DUNE, energy_1tau_p_Ar, xsec_1tau_p_Ar_upper)
DUNE_xsec_1tau_n_Ar_matched_upper = MatchXSec(energy_DUNE, energy_1tau_n_Ar, xsec_1tau_n_Ar_upper)

DUNE_xsec_1tau_coh_Ar_matched_lower = MatchXSec(energy_DUNE, energy_1tau_coh_Ar, xsec_1tau_coh_Ar_lower)
DUNE_xsec_1tau_p_Ar_matched_lower = MatchXSec(energy_DUNE, energy_1tau_p_Ar, xsec_1tau_p_Ar_lower)
DUNE_xsec_1tau_n_Ar_matched_lower = MatchXSec(energy_DUNE, energy_1tau_n_Ar, xsec_1tau_n_Ar_lower)

DUNE_xsec_2tau_coh_Ar_matched_upper = MatchXSec(energy_DUNE, energy_2tau_coh_Ar, xsec_2tau_coh_Ar_upper)
DUNE_xsec_2tau_p_Ar_matched_upper = MatchXSec(energy_DUNE, energy_2tau_p_Ar, xsec_2tau_p_Ar_upper)
DUNE_xsec_2tau_n_Ar_matched_upper = MatchXSec(energy_DUNE, energy_2tau_n_Ar, xsec_2tau_n_Ar_upper)

DUNE_xsec_2tau_coh_Ar_matched_lower = MatchXSec(energy_DUNE, energy_2tau_coh_Ar, xsec_2tau_coh_Ar_lower)
DUNE_xsec_2tau_p_Ar_matched_lower = MatchXSec(energy_DUNE, energy_2tau_p_Ar, xsec_2tau_p_Ar_lower)
DUNE_xsec_2tau_n_Ar_matched_lower = MatchXSec(energy_DUNE, energy_2tau_n_Ar, xsec_2tau_n_Ar_lower)

DUNE_xsec_2mu_coh_Ar_matched_upper = MatchXSec(energy_DUNE, energy_2mu_coh_Ar, xsec_2mu_coh_Ar_upper)
DUNE_xsec_2mu_p_Ar_matched_upper = MatchXSec(energy_DUNE, energy_2mu_p_Ar, xsec_2mu_p_Ar_upper)
DUNE_xsec_2mu_n_Ar_matched_upper = MatchXSec(energy_DUNE, energy_2mu_n_Ar, xsec_2mu_n_Ar_upper)

DUNE_xsec_2mu_coh_Ar_matched_lower = MatchXSec(energy_DUNE, energy_2mu_coh_Ar, xsec_2mu_coh_Ar_lower)
DUNE_xsec_2mu_p_Ar_matched_lower = MatchXSec(energy_DUNE, energy_2mu_p_Ar, xsec_2mu_p_Ar_lower)
DUNE_xsec_2mu_n_Ar_matched_lower = MatchXSec(energy_DUNE, energy_2mu_n_Ar, xsec_2mu_n_Ar_lower)

DUNE_xsec_ve1tau_coh_Ar_matched_upper = MatchXSec(energy_DUNE, energy_ve1tau_coh_Ar, xsec_ve1tau_coh_Ar_upper)
DUNE_xsec_ve1tau_p_Ar_matched_upper = MatchXSec(energy_DUNE, energy_ve1tau_p_Ar, xsec_ve1tau_p_Ar_upper)
DUNE_xsec_ve1tau_n_Ar_matched_upper = MatchXSec(energy_DUNE, energy_ve1tau_n_Ar, xsec_ve1tau_n_Ar_upper)

DUNE_xsec_ve1tau_coh_Ar_matched_lower = MatchXSec(energy_DUNE, energy_ve1tau_coh_Ar, xsec_ve1tau_coh_Ar_lower)
DUNE_xsec_ve1tau_p_Ar_matched_lower = MatchXSec(energy_DUNE, energy_ve1tau_p_Ar, xsec_ve1tau_p_Ar_lower)
DUNE_xsec_ve1tau_n_Ar_matched_lower = MatchXSec(energy_DUNE, energy_ve1tau_n_Ar, xsec_ve1tau_n_Ar_lower)

### DUNE Tau-optimized Flux ; Matched XSec ###
DUNE_tau_opt_xsec_1tau_coh_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_1tau_coh_Ar, xsec_1tau_coh_Ar)
DUNE_tau_opt_xsec_1tau_p_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_1tau_p_Ar, xsec_1tau_p_Ar)
DUNE_tau_opt_xsec_1tau_n_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_1tau_n_Ar, xsec_1tau_n_Ar)

DUNE_tau_opt_xsec_2tau_coh_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_2tau_coh_Ar, xsec_2tau_coh_Ar)
DUNE_tau_opt_xsec_2tau_p_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_2tau_p_Ar, xsec_2tau_p_Ar)
DUNE_tau_opt_xsec_2tau_n_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_2tau_n_Ar, xsec_2tau_n_Ar)

DUNE_tau_opt_xsec_2mu_coh_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_2mu_coh_Ar, xsec_2mu_coh_Ar)
DUNE_tau_opt_xsec_2mu_p_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_2mu_p_Ar, xsec_2mu_p_Ar)
DUNE_tau_opt_xsec_2mu_n_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_2mu_n_Ar, xsec_2mu_n_Ar)

DUNE_tau_opt_xsec_ve1tau_coh_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_ve1tau_coh_Ar, xsec_ve1tau_coh_Ar)
DUNE_tau_opt_xsec_ve1tau_p_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_ve1tau_p_Ar, xsec_ve1tau_p_Ar)
DUNE_tau_opt_xsec_ve1tau_n_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_ve1tau_n_Ar, xsec_ve1tau_n_Ar)

# EPA
DUNE_tau_opt_EPA_Heaviside_xsec_1tau_coh_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_1tau_coh_Ar, EPA_Heaviside_xsec_1tau_coh_Ar)
DUNE_tau_opt_EPA_Heaviside_xsec_2mu_coh_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_2mu_coh_Ar, EPA_Heaviside_xsec_2mu_coh_Ar)

# Upper and Lower limits #
DUNE_tau_opt_xsec_1tau_coh_Ar_matched_upper = MatchXSec(energy_DUNE_tau_opt, energy_1tau_coh_Ar, xsec_1tau_coh_Ar_upper)
DUNE_tau_opt_xsec_1tau_p_Ar_matched_upper = MatchXSec(energy_DUNE_tau_opt, energy_1tau_p_Ar, xsec_1tau_p_Ar_upper)
DUNE_tau_opt_xsec_1tau_n_Ar_matched_upper = MatchXSec(energy_DUNE_tau_opt, energy_1tau_n_Ar, xsec_1tau_n_Ar_upper)

DUNE_tau_opt_xsec_1tau_coh_Ar_matched_lower = MatchXSec(energy_DUNE_tau_opt, energy_1tau_coh_Ar, xsec_1tau_coh_Ar_lower)
DUNE_tau_opt_xsec_1tau_p_Ar_matched_lower = MatchXSec(energy_DUNE_tau_opt, energy_1tau_p_Ar, xsec_1tau_p_Ar_lower)
DUNE_tau_opt_xsec_1tau_n_Ar_matched_lower = MatchXSec(energy_DUNE_tau_opt, energy_1tau_n_Ar, xsec_1tau_n_Ar_lower)

DUNE_tau_opt_xsec_2tau_coh_Ar_matched_upper = MatchXSec(energy_DUNE_tau_opt, energy_2tau_coh_Ar, xsec_2tau_coh_Ar_upper)
DUNE_tau_opt_xsec_2tau_p_Ar_matched_upper = MatchXSec(energy_DUNE_tau_opt, energy_2tau_p_Ar, xsec_2tau_p_Ar_upper)
DUNE_tau_opt_xsec_2tau_n_Ar_matched_upper = MatchXSec(energy_DUNE_tau_opt, energy_2tau_n_Ar, xsec_2tau_n_Ar_upper)

DUNE_tau_opt_xsec_2tau_coh_Ar_matched_lower = MatchXSec(energy_DUNE_tau_opt, energy_2tau_coh_Ar, xsec_2tau_coh_Ar_lower)
DUNE_tau_opt_xsec_2tau_p_Ar_matched_lower = MatchXSec(energy_DUNE_tau_opt, energy_2tau_p_Ar, xsec_2tau_p_Ar_lower)
DUNE_tau_opt_xsec_2tau_n_Ar_matched_lower = MatchXSec(energy_DUNE_tau_opt, energy_2tau_n_Ar, xsec_2tau_n_Ar_lower)

DUNE_tau_opt_xsec_2mu_coh_Ar_matched_upper = MatchXSec(energy_DUNE_tau_opt, energy_2mu_coh_Ar, xsec_2mu_coh_Ar_upper)
DUNE_tau_opt_xsec_2mu_p_Ar_matched_upper = MatchXSec(energy_DUNE_tau_opt, energy_2mu_p_Ar, xsec_2mu_p_Ar_upper)
DUNE_tau_opt_xsec_2mu_n_Ar_matched_upper = MatchXSec(energy_DUNE_tau_opt, energy_2mu_n_Ar, xsec_2mu_n_Ar_upper)

DUNE_tau_opt_xsec_2mu_coh_Ar_matched_lower = MatchXSec(energy_DUNE_tau_opt, energy_2mu_coh_Ar, xsec_2mu_coh_Ar_lower)
DUNE_tau_opt_xsec_2mu_p_Ar_matched_lower = MatchXSec(energy_DUNE_tau_opt, energy_2mu_p_Ar, xsec_2mu_p_Ar_lower)
DUNE_tau_opt_xsec_2mu_n_Ar_matched_lower = MatchXSec(energy_DUNE_tau_opt, energy_2mu_n_Ar, xsec_2mu_n_Ar_lower)

DUNE_tau_opt_xsec_ve1tau_coh_Ar_matched_upper = MatchXSec(energy_DUNE_tau_opt, energy_ve1tau_coh_Ar, xsec_ve1tau_coh_Ar_upper)
DUNE_tau_opt_xsec_ve1tau_p_Ar_matched_upper = MatchXSec(energy_DUNE_tau_opt, energy_ve1tau_p_Ar, xsec_ve1tau_p_Ar_upper)
DUNE_tau_opt_xsec_ve1tau_n_Ar_matched_upper = MatchXSec(energy_DUNE_tau_opt, energy_ve1tau_n_Ar, xsec_ve1tau_n_Ar_upper)

DUNE_tau_opt_xsec_ve1tau_coh_Ar_matched_lower = MatchXSec(energy_DUNE_tau_opt, energy_ve1tau_coh_Ar, xsec_ve1tau_coh_Ar_lower)
DUNE_tau_opt_xsec_ve1tau_p_Ar_matched_lower = MatchXSec(energy_DUNE_tau_opt, energy_ve1tau_p_Ar, xsec_ve1tau_p_Ar_lower)
DUNE_tau_opt_xsec_ve1tau_n_Ar_matched_lower = MatchXSec(energy_DUNE_tau_opt, energy_ve1tau_n_Ar, xsec_ve1tau_n_Ar_lower)

### FASERnu vmu Flux ; Matched XSec ###
FASERvmu_xsec_1tau_coh_W_matched = MatchXSec(energy_FASERvmu, energy_1tau_coh_W, xsec_1tau_coh_W)
FASERvmu_xsec_1tau_p_W_matched = MatchXSec(energy_FASERvmu, energy_1tau_p_W, xsec_1tau_p_W)
FASERvmu_xsec_1tau_n_W_matched = MatchXSec(energy_FASERvmu, energy_1tau_n_W, xsec_1tau_n_W)

FASERvmu_xsec_2tau_coh_W_matched = MatchXSec(energy_FASERvmu, energy_2tau_coh_W, xsec_2tau_coh_W)
FASERvmu_xsec_2tau_p_W_matched = MatchXSec(energy_FASERvmu, energy_2tau_p_W, xsec_2tau_p_W)
FASERvmu_xsec_2tau_n_W_matched = MatchXSec(energy_FASERvmu, energy_2tau_n_W, xsec_2tau_n_W)

FASERvmu_xsec_2mu_coh_W_matched = MatchXSec(energy_FASERvmu, energy_2mu_coh_W, xsec_2mu_coh_W)
FASERvmu_xsec_2mu_p_W_matched = MatchXSec(energy_FASERvmu, energy_2mu_p_W, xsec_2mu_p_W)
FASERvmu_xsec_2mu_n_W_matched = MatchXSec(energy_FASERvmu, energy_2mu_n_W, xsec_2mu_n_W)

# Cross check with vmu CC #
FASERvmu_xsec_vmuCC_matched = MatchXSec(energy_FASERvmu, energy_vmuCC, xsec_vmuCC)

FASERvmu_xsec_DIS_vmuCC_matched = MatchXSec(energy_FASERvmu, energy_DIS_vmuCC_tungsten, xsec_DIS_vmuCC_tungsten)
FASERvmu_xsec_DIS_vmubarCC_matched = MatchXSec(energy_FASERvmu, energy_DIS_vmubarCC_tungsten, xsec_DIS_vmubarCC_tungsten)

# Upper and Lower limits #
FASERvmu_xsec_1tau_coh_W_matched_upper = MatchXSec(energy_FASERvmu, energy_1tau_coh_W, xsec_1tau_coh_W_upper)
FASERvmu_xsec_1tau_p_W_matched_upper = MatchXSec(energy_FASERvmu, energy_1tau_p_W, xsec_1tau_p_W_upper)
FASERvmu_xsec_1tau_n_W_matched_upper = MatchXSec(energy_FASERvmu, energy_1tau_n_W, xsec_1tau_n_W_upper)

FASERvmu_xsec_1tau_coh_W_matched_lower = MatchXSec(energy_FASERvmu, energy_1tau_coh_W, xsec_1tau_coh_W_lower)
FASERvmu_xsec_1tau_p_W_matched_lower = MatchXSec(energy_FASERvmu, energy_1tau_p_W, xsec_1tau_p_W_lower)
FASERvmu_xsec_1tau_n_W_matched_lower = MatchXSec(energy_FASERvmu, energy_1tau_n_W, xsec_1tau_n_W_lower)

FASERvmu_xsec_2tau_coh_W_matched_upper = MatchXSec(energy_FASERvmu, energy_2tau_coh_W, xsec_2tau_coh_W_upper)
FASERvmu_xsec_2tau_p_W_matched_upper = MatchXSec(energy_FASERvmu, energy_2tau_p_W, xsec_2tau_p_W_upper)
FASERvmu_xsec_2tau_n_W_matched_upper = MatchXSec(energy_FASERvmu, energy_2tau_n_W, xsec_2tau_n_W_upper)

FASERvmu_xsec_2tau_coh_W_matched_lower = MatchXSec(energy_FASERvmu, energy_2tau_coh_W, xsec_2tau_coh_W_lower)
FASERvmu_xsec_2tau_p_W_matched_lower = MatchXSec(energy_FASERvmu, energy_2tau_p_W, xsec_2tau_p_W_lower)
FASERvmu_xsec_2tau_n_W_matched_lower = MatchXSec(energy_FASERvmu, energy_2tau_n_W, xsec_2tau_n_W_lower)

FASERvmu_xsec_2mu_coh_W_matched_upper = MatchXSec(energy_FASERvmu, energy_2mu_coh_W, xsec_2mu_coh_W_upper)
FASERvmu_xsec_2mu_p_W_matched_upper = MatchXSec(energy_FASERvmu, energy_2mu_p_W, xsec_2mu_p_W_upper)
FASERvmu_xsec_2mu_n_W_matched_upper = MatchXSec(energy_FASERvmu, energy_2mu_n_W, xsec_2mu_n_W_upper)

FASERvmu_xsec_2mu_coh_W_matched_lower = MatchXSec(energy_FASERvmu, energy_2mu_coh_W, xsec_2mu_coh_W_lower)
FASERvmu_xsec_2mu_p_W_matched_lower = MatchXSec(energy_FASERvmu, energy_2mu_p_W, xsec_2mu_p_W_lower)
FASERvmu_xsec_2mu_n_W_matched_lower = MatchXSec(energy_FASERvmu, energy_2mu_n_W, xsec_2mu_n_W_lower)

### FASERnu2 vmu + vmubar Flux ; Matched XSec ###
FASERv2_xsec_1tau_coh_W_matched = MatchXSec(energy_FASERv2, energy_1tau_coh_W, xsec_1tau_coh_W)
FASERv2_xsec_1tau_p_W_matched = MatchXSec(energy_FASERv2, energy_1tau_p_W, xsec_1tau_p_W)
FASERv2_xsec_1tau_n_W_matched = MatchXSec(energy_FASERv2, energy_1tau_n_W, xsec_1tau_n_W)

FASERv2_xsec_2tau_coh_W_matched = MatchXSec(energy_FASERv2, energy_2tau_coh_W, xsec_2tau_coh_W)
FASERv2_xsec_2tau_p_W_matched = MatchXSec(energy_FASERv2, energy_2tau_p_W, xsec_2tau_p_W)
FASERv2_xsec_2tau_n_W_matched = MatchXSec(energy_FASERv2, energy_2tau_n_W, xsec_2tau_n_W)

FASERv2_xsec_2mu_coh_W_matched = MatchXSec(energy_FASERv2, energy_2mu_coh_W, xsec_2mu_coh_W)
FASERv2_xsec_2mu_p_W_matched = MatchXSec(energy_FASERv2, energy_2mu_p_W, xsec_2mu_p_W)
FASERv2_xsec_2mu_n_W_matched = MatchXSec(energy_FASERv2, energy_2mu_n_W, xsec_2mu_n_W)

# Cross check with vmu CC #
FASERv2_xsec_vmuCC_matched = MatchXSec(energy_FASERv2, energy_vmuCC, xsec_vmuCC)

# Upper and Lower limits #
FASERv2_xsec_1tau_coh_W_matched_upper = MatchXSec(energy_FASERv2, energy_1tau_coh_W, xsec_1tau_coh_W_upper)
FASERv2_xsec_1tau_p_W_matched_upper = MatchXSec(energy_FASERv2, energy_1tau_p_W, xsec_1tau_p_W_upper)
FASERv2_xsec_1tau_n_W_matched_upper = MatchXSec(energy_FASERv2, energy_1tau_n_W, xsec_1tau_n_W_upper)

FASERv2_xsec_1tau_coh_W_matched_lower = MatchXSec(energy_FASERv2, energy_1tau_coh_W, xsec_1tau_coh_W_lower)
FASERv2_xsec_1tau_p_W_matched_lower = MatchXSec(energy_FASERv2, energy_1tau_p_W, xsec_1tau_p_W_lower)
FASERv2_xsec_1tau_n_W_matched_lower = MatchXSec(energy_FASERv2, energy_1tau_n_W, xsec_1tau_n_W_lower)

FASERv2_xsec_2tau_coh_W_matched_upper = MatchXSec(energy_FASERv2, energy_2tau_coh_W, xsec_2tau_coh_W_upper)
FASERv2_xsec_2tau_p_W_matched_upper = MatchXSec(energy_FASERv2, energy_2tau_p_W, xsec_2tau_p_W_upper)
FASERv2_xsec_2tau_n_W_matched_upper = MatchXSec(energy_FASERv2, energy_2tau_n_W, xsec_2tau_n_W_upper)

FASERv2_xsec_2tau_coh_W_matched_lower = MatchXSec(energy_FASERv2, energy_2tau_coh_W, xsec_2tau_coh_W_lower)
FASERv2_xsec_2tau_p_W_matched_lower = MatchXSec(energy_FASERv2, energy_2tau_p_W, xsec_2tau_p_W_lower)
FASERv2_xsec_2tau_n_W_matched_lower = MatchXSec(energy_FASERv2, energy_2tau_n_W, xsec_2tau_n_W_lower)

FASERv2_xsec_2mu_coh_W_matched_upper = MatchXSec(energy_FASERv2, energy_2mu_coh_W, xsec_2mu_coh_W_upper)
FASERv2_xsec_2mu_p_W_matched_upper = MatchXSec(energy_FASERv2, energy_2mu_p_W, xsec_2mu_p_W_upper)
FASERv2_xsec_2mu_n_W_matched_upper = MatchXSec(energy_FASERv2, energy_2mu_n_W, xsec_2mu_n_W_upper)

FASERv2_xsec_2mu_coh_W_matched_lower = MatchXSec(energy_FASERv2, energy_2mu_coh_W, xsec_2mu_coh_W_lower)
FASERv2_xsec_2mu_p_W_matched_lower = MatchXSec(energy_FASERv2, energy_2mu_p_W, xsec_2mu_p_W_lower)
FASERv2_xsec_2mu_n_W_matched_lower = MatchXSec(energy_FASERv2, energy_2mu_n_W, xsec_2mu_n_W_lower)

### MINOS Matched XSec ###
MINOS_xsec_vmuCC_matched = MatchXSec(energy_MINOS, energy_vmuCC, np.multiply(xsec_vmuCC, A_Fe/A_W))

MINOS_xsec_1tau_coh_Fe_matched = MatchXSec(energy_MINOS, energy_1tau_coh_Fe, xsec_1tau_coh_Fe)
MINOS_xsec_1tau_p_Fe_matched = MatchXSec(energy_MINOS, energy_1tau_p_Fe, xsec_1tau_p_Fe)
MINOS_xsec_1tau_n_Fe_matched = MatchXSec(energy_MINOS, energy_1tau_n_Fe, xsec_1tau_n_Fe)

MINOS_xsec_2mu_coh_Fe_matched = MatchXSec(energy_MINOS, energy_2mu_coh_Fe, xsec_2mu_coh_Fe)
MINOS_xsec_2mu_p_Fe_matched = MatchXSec(energy_MINOS, energy_2mu_p_Fe, xsec_2mu_p_Fe)
MINOS_xsec_2mu_n_Fe_matched = MatchXSec(energy_MINOS, energy_2mu_n_Fe, xsec_2mu_n_Fe)

# Upper and Lower limits #
MINOS_xsec_1tau_coh_Fe_matched_upper = MatchXSec(energy_MINOS, energy_1tau_coh_Fe, xsec_1tau_coh_Fe_upper)
MINOS_xsec_1tau_p_Fe_matched_upper = MatchXSec(energy_MINOS, energy_1tau_p_Fe, xsec_1tau_p_Fe_upper)
MINOS_xsec_1tau_n_Fe_matched_upper = MatchXSec(energy_MINOS, energy_1tau_n_Fe, xsec_1tau_n_Fe_upper)

MINOS_xsec_1tau_coh_Fe_matched_lower = MatchXSec(energy_MINOS, energy_1tau_coh_Fe, xsec_1tau_coh_Fe_lower)
MINOS_xsec_1tau_p_Fe_matched_lower = MatchXSec(energy_MINOS, energy_1tau_p_Fe, xsec_1tau_p_Fe_lower)
MINOS_xsec_1tau_n_Fe_matched_lower = MatchXSec(energy_MINOS, energy_1tau_n_Fe, xsec_1tau_n_Fe_lower)

MINOS_xsec_2mu_coh_Fe_matched_upper = MatchXSec(energy_MINOS, energy_2mu_coh_Fe, xsec_2mu_coh_Fe_upper)
MINOS_xsec_2mu_p_Fe_matched_upper = MatchXSec(energy_MINOS, energy_2mu_p_Fe, xsec_2mu_p_Fe_upper)
MINOS_xsec_2mu_n_Fe_matched_upper = MatchXSec(energy_MINOS, energy_2mu_n_Fe, xsec_2mu_n_Fe_upper)

MINOS_xsec_2mu_coh_Fe_matched_lower = MatchXSec(energy_MINOS, energy_2mu_coh_Fe, xsec_2mu_coh_Fe_lower)
MINOS_xsec_2mu_p_Fe_matched_lower = MatchXSec(energy_MINOS, energy_2mu_p_Fe, xsec_2mu_p_Fe_lower)
MINOS_xsec_2mu_n_Fe_matched_lower = MatchXSec(energy_MINOS, energy_2mu_n_Fe, xsec_2mu_n_Fe_lower)

### MINOSPlus Matched XSec ###
MINOSPlus_xsec_vmuCC_matched = MatchXSec(energy_MINOSPlus, energy_vmuCC, np.multiply(xsec_vmuCC, A_Fe/A_W))

MINOSPlus_xsec_1tau_coh_Fe_matched = MatchXSec(energy_MINOSPlus, energy_1tau_coh_Fe, xsec_1tau_coh_Fe)
MINOSPlus_xsec_1tau_p_Fe_matched = MatchXSec(energy_MINOSPlus, energy_1tau_p_Fe, xsec_1tau_p_Fe)
MINOSPlus_xsec_1tau_n_Fe_matched = MatchXSec(energy_MINOSPlus, energy_1tau_n_Fe, xsec_1tau_n_Fe)

MINOSPlus_xsec_2mu_coh_Fe_matched = MatchXSec(energy_MINOSPlus, energy_2mu_coh_Fe, xsec_2mu_coh_Fe)
MINOSPlus_xsec_2mu_p_Fe_matched = MatchXSec(energy_MINOSPlus, energy_2mu_p_Fe, xsec_2mu_p_Fe)
MINOSPlus_xsec_2mu_n_Fe_matched = MatchXSec(energy_MINOSPlus, energy_2mu_n_Fe, xsec_2mu_n_Fe)

# Upper and Lower limits #
MINOSPlus_xsec_1tau_coh_Fe_matched_upper = MatchXSec(energy_MINOSPlus, energy_1tau_coh_Fe, xsec_1tau_coh_Fe_upper)
MINOSPlus_xsec_1tau_p_Fe_matched_upper = MatchXSec(energy_MINOSPlus, energy_1tau_p_Fe, xsec_1tau_p_Fe_upper)
MINOSPlus_xsec_1tau_n_Fe_matched_upper = MatchXSec(energy_MINOSPlus, energy_1tau_n_Fe, xsec_1tau_n_Fe_upper)

MINOSPlus_xsec_1tau_coh_Fe_matched_lower = MatchXSec(energy_MINOSPlus, energy_1tau_coh_Fe, xsec_1tau_coh_Fe_lower)
MINOSPlus_xsec_1tau_p_Fe_matched_lower = MatchXSec(energy_MINOSPlus, energy_1tau_p_Fe, xsec_1tau_p_Fe_lower)
MINOSPlus_xsec_1tau_n_Fe_matched_lower = MatchXSec(energy_MINOSPlus, energy_1tau_n_Fe, xsec_1tau_n_Fe_lower)

MINOSPlus_xsec_2mu_coh_Fe_matched_upper = MatchXSec(energy_MINOSPlus, energy_2mu_coh_Fe, xsec_2mu_coh_Fe_upper)
MINOSPlus_xsec_2mu_p_Fe_matched_upper = MatchXSec(energy_MINOSPlus, energy_2mu_p_Fe, xsec_2mu_p_Fe_upper)
MINOSPlus_xsec_2mu_n_Fe_matched_upper = MatchXSec(energy_MINOSPlus, energy_2mu_n_Fe, xsec_2mu_n_Fe_upper)

MINOSPlus_xsec_2mu_coh_Fe_matched_lower = MatchXSec(energy_MINOSPlus, energy_2mu_coh_Fe, xsec_2mu_coh_Fe_lower)
MINOSPlus_xsec_2mu_p_Fe_matched_lower = MatchXSec(energy_MINOSPlus, energy_2mu_p_Fe, xsec_2mu_p_Fe_lower)
MINOSPlus_xsec_2mu_n_Fe_matched_lower = MatchXSec(energy_MINOSPlus, energy_2mu_n_Fe, xsec_2mu_n_Fe_lower)

### T2K - INGRID Matched XSec ###
INGRID_xsec_1tau_coh_Fe_matched = MatchXSec(energy_INGRID, energy_1tau_coh_Fe, xsec_1tau_coh_Fe)
INGRID_xsec_1tau_p_Fe_matched = MatchXSec(energy_INGRID, energy_1tau_p_Fe, xsec_1tau_p_Fe)
INGRID_xsec_1tau_n_Fe_matched = MatchXSec(energy_INGRID, energy_1tau_n_Fe, xsec_1tau_n_Fe)

INGRID_xsec_2mu_coh_Fe_matched = MatchXSec(energy_INGRID, energy_2mu_coh_Fe, xsec_2mu_coh_Fe)
INGRID_xsec_2mu_p_Fe_matched = MatchXSec(energy_INGRID, energy_2mu_p_Fe, xsec_2mu_p_Fe)
INGRID_xsec_2mu_n_Fe_matched = MatchXSec(energy_INGRID, energy_2mu_n_Fe, xsec_2mu_n_Fe)

# Upper and Lower limits #
INGRID_xsec_1tau_coh_Fe_matched_upper = MatchXSec(energy_INGRID, energy_1tau_coh_Fe, xsec_1tau_coh_Fe_upper)
INGRID_xsec_1tau_p_Fe_matched_upper = MatchXSec(energy_INGRID, energy_1tau_p_Fe, xsec_1tau_p_Fe_upper)
INGRID_xsec_1tau_n_Fe_matched_upper = MatchXSec(energy_INGRID, energy_1tau_n_Fe, xsec_1tau_n_Fe_upper)

INGRID_xsec_1tau_coh_Fe_matched_lower = MatchXSec(energy_INGRID, energy_1tau_coh_Fe, xsec_1tau_coh_Fe_lower)
INGRID_xsec_1tau_p_Fe_matched_lower = MatchXSec(energy_INGRID, energy_1tau_p_Fe, xsec_1tau_p_Fe_lower)
INGRID_xsec_1tau_n_Fe_matched_lower = MatchXSec(energy_INGRID, energy_1tau_n_Fe, xsec_1tau_n_Fe_lower)

INGRID_xsec_2mu_coh_Fe_matched_upper = MatchXSec(energy_INGRID, energy_2mu_coh_Fe, xsec_2mu_coh_Fe_upper)
INGRID_xsec_2mu_p_Fe_matched_upper = MatchXSec(energy_INGRID, energy_2mu_p_Fe, xsec_2mu_p_Fe_upper)
INGRID_xsec_2mu_n_Fe_matched_upper = MatchXSec(energy_INGRID, energy_2mu_n_Fe, xsec_2mu_n_Fe_upper)

INGRID_xsec_2mu_coh_Fe_matched_lower = MatchXSec(energy_INGRID, energy_2mu_coh_Fe, xsec_2mu_coh_Fe_lower)
INGRID_xsec_2mu_p_Fe_matched_lower = MatchXSec(energy_INGRID, energy_2mu_p_Fe, xsec_2mu_p_Fe_lower)
INGRID_xsec_2mu_n_Fe_matched_lower = MatchXSec(energy_INGRID, energy_2mu_n_Fe, xsec_2mu_n_Fe_lower)

### SHiP ; Matched XSec ###
SHiP_xsec_1tau_coh_W_matched = MatchXSec(energy_SHiP, energy_1tau_coh_W, xsec_1tau_coh_W)
SHiP_xsec_1tau_p_W_matched = MatchXSec(energy_SHiP, energy_1tau_p_W, xsec_1tau_p_W)
SHiP_xsec_1tau_n_W_matched = MatchXSec(energy_SHiP, energy_1tau_n_W, xsec_1tau_n_W)

SHiP_xsec_2tau_coh_W_matched = MatchXSec(energy_SHiP, energy_2tau_coh_W, xsec_2tau_coh_W)
SHiP_xsec_2tau_p_W_matched = MatchXSec(energy_SHiP, energy_2tau_p_W, xsec_2tau_p_W)
SHiP_xsec_2tau_n_W_matched = MatchXSec(energy_SHiP, energy_2tau_n_W, xsec_2tau_n_W)

SHiP_xsec_2mu_coh_W_matched = MatchXSec(energy_SHiP, energy_2mu_coh_W, xsec_2mu_coh_W)
SHiP_xsec_2mu_p_W_matched = MatchXSec(energy_SHiP, energy_2mu_p_W, xsec_2mu_p_W)
SHiP_xsec_2mu_n_W_matched = MatchXSec(energy_SHiP, energy_2mu_n_W, xsec_2mu_n_W)

# Upper and Lower limits #
SHiP_xsec_1tau_coh_W_matched_upper = MatchXSec(energy_SHiP, energy_1tau_coh_W, xsec_1tau_coh_W_upper)
SHiP_xsec_1tau_p_W_matched_upper = MatchXSec(energy_SHiP, energy_1tau_p_W, xsec_1tau_p_W_upper)
SHiP_xsec_1tau_n_W_matched_upper = MatchXSec(energy_SHiP, energy_1tau_n_W, xsec_1tau_n_W_upper)

SHiP_xsec_1tau_coh_W_matched_lower = MatchXSec(energy_SHiP, energy_1tau_coh_W, xsec_1tau_coh_W_lower)
SHiP_xsec_1tau_p_W_matched_lower = MatchXSec(energy_SHiP, energy_1tau_p_W, xsec_1tau_p_W_lower)
SHiP_xsec_1tau_n_W_matched_lower = MatchXSec(energy_SHiP, energy_1tau_n_W, xsec_1tau_n_W_lower)

SHiP_xsec_2tau_coh_W_matched_upper = MatchXSec(energy_SHiP, energy_2tau_coh_W, xsec_2tau_coh_W_upper)
SHiP_xsec_2tau_p_W_matched_upper = MatchXSec(energy_SHiP, energy_2tau_p_W, xsec_2tau_p_W_upper)
SHiP_xsec_2tau_n_W_matched_upper = MatchXSec(energy_SHiP, energy_2tau_n_W, xsec_2tau_n_W_upper)

SHiP_xsec_2tau_coh_W_matched_lower = MatchXSec(energy_SHiP, energy_2tau_coh_W, xsec_2tau_coh_W_lower)
SHiP_xsec_2tau_p_W_matched_lower = MatchXSec(energy_SHiP, energy_2tau_p_W, xsec_2tau_p_W_lower)
SHiP_xsec_2tau_n_W_matched_lower = MatchXSec(energy_SHiP, energy_2tau_n_W, xsec_2tau_n_W_lower)

SHiP_xsec_2mu_coh_W_matched_upper = MatchXSec(energy_SHiP, energy_2mu_coh_W, xsec_2mu_coh_W_upper)
SHiP_xsec_2mu_p_W_matched_upper = MatchXSec(energy_SHiP, energy_2mu_p_W, xsec_2mu_p_W_upper)
SHiP_xsec_2mu_n_W_matched_upper = MatchXSec(energy_SHiP, energy_2mu_n_W, xsec_2mu_n_W_upper)

SHiP_xsec_2mu_coh_W_matched_lower = MatchXSec(energy_SHiP, energy_2mu_coh_W, xsec_2mu_coh_W_lower)
SHiP_xsec_2mu_p_W_matched_lower = MatchXSec(energy_SHiP, energy_2mu_p_W, xsec_2mu_p_W_lower)
SHiP_xsec_2mu_n_W_matched_lower = MatchXSec(energy_SHiP, energy_2mu_n_W, xsec_2mu_n_W_lower)

### SBND ; Matched XSec ###
SBND_xsec_1e1mu_coh_Ar_matched = MatchXSec(energy_SBND_vmu, energy_1e1mu_coh_Ar, xsec_1e1mu_coh_Ar)
SBND_xsec_1e1mu_p_Ar_matched = MatchXSec(energy_SBND_vmu, energy_1e1mu_p_Ar, xsec_1e1mu_p_Ar)
SBND_xsec_1e1mu_n_Ar_matched = MatchXSec(energy_SBND_vmu, energy_1e1mu_n_Ar, xsec_1e1mu_n_Ar)

SBND_xsec_2mu_coh_Ar_matched = MatchXSec(energy_SBND_vmu, energy_2mu_coh_Ar, xsec_2mu_coh_Ar)
SBND_xsec_2mu_p_Ar_matched = MatchXSec(energy_SBND_vmu, energy_2mu_p_Ar, xsec_2mu_p_Ar)
SBND_xsec_2mu_n_Ar_matched = MatchXSec(energy_SBND_vmu, energy_2mu_n_Ar, xsec_2mu_n_Ar)

# Upper and Lower limits #
SBND_xsec_1e1mu_coh_Ar_matched_upper = MatchXSec(energy_SBND_vmu, energy_1e1mu_coh_Ar, xsec_1e1mu_coh_Ar_upper)
SBND_xsec_1e1mu_p_Ar_matched_upper = MatchXSec(energy_SBND_vmu, energy_1e1mu_p_Ar, xsec_1e1mu_p_Ar_upper)
SBND_xsec_1e1mu_n_Ar_matched_upper = MatchXSec(energy_SBND_vmu, energy_1e1mu_n_Ar, xsec_1e1mu_n_Ar_upper)

SBND_xsec_1e1mu_coh_Ar_matched_lower = MatchXSec(energy_SBND_vmu, energy_1e1mu_coh_Ar, xsec_1e1mu_coh_Ar_lower)
SBND_xsec_1e1mu_p_Ar_matched_lower = MatchXSec(energy_SBND_vmu, energy_1e1mu_p_Ar, xsec_1e1mu_p_Ar_lower)
SBND_xsec_1e1mu_n_Ar_matched_lower = MatchXSec(energy_SBND_vmu, energy_1e1mu_n_Ar, xsec_1e1mu_n_Ar_lower)

SBND_xsec_2mu_coh_Ar_matched_upper = MatchXSec(energy_SBND_vmu, energy_2mu_coh_Ar, xsec_2mu_coh_Ar_upper)
SBND_xsec_2mu_p_Ar_matched_upper = MatchXSec(energy_SBND_vmu, energy_2mu_p_Ar, xsec_2mu_p_Ar_upper)
SBND_xsec_2mu_n_Ar_matched_upper = MatchXSec(energy_SBND_vmu, energy_2mu_n_Ar, xsec_2mu_n_Ar_upper)

SBND_xsec_2mu_coh_Ar_matched_lower = MatchXSec(energy_SBND_vmu, energy_2mu_coh_Ar, xsec_2mu_coh_Ar_lower)
SBND_xsec_2mu_p_Ar_matched_lower = MatchXSec(energy_SBND_vmu, energy_2mu_p_Ar, xsec_2mu_p_Ar_lower)
SBND_xsec_2mu_n_Ar_matched_lower = MatchXSec(energy_SBND_vmu, energy_2mu_n_Ar, xsec_2mu_n_Ar_lower)

###### Correction Factors ######
# RMiss Correction Factor
DUNE_correction_factor_RMiss_1tau_coh_Ar_matched = MatchEps(energy_DUNE, correction_factor_RMiss_energy, correction_factor_RMiss_coh)
DUNE_correction_factor_RMiss_1tau_p_Ar_matched = MatchEps(energy_DUNE, correction_factor_RMiss_energy, correction_factor_RMiss_p)
DUNE_correction_factor_RMiss_1tau_n_Ar_matched = MatchEps(energy_DUNE, correction_factor_RMiss_energy, correction_factor_RMiss_n)

DUNE_tau_opt_correction_factor_RMiss_1tau_coh_Ar_matched = MatchEps(energy_DUNE_tau_opt, correction_factor_RMiss_energy, correction_factor_RMiss_coh)
DUNE_tau_opt_correction_factor_RMiss_1tau_p_Ar_matched = MatchEps(energy_DUNE_tau_opt, correction_factor_RMiss_energy, correction_factor_RMiss_p)
DUNE_tau_opt_correction_factor_RMiss_1tau_n_Ar_matched = MatchEps(energy_DUNE_tau_opt, correction_factor_RMiss_energy, correction_factor_RMiss_n)

##########################################
###### Number of Events Calculation ######
##########################################

def CalculateEvents(flux, flux_energy, xsec_matched, NTONNES=1, NYEAR=1, normalized=False, Phi=1, detector='DUNE', correction_factor=None):
    """
    Calculate the expected number of events for a given neutrino flux and trident process.

    N = XSecCon * MD * N_POT * NTONNES * NYEAR / MAr                                (DUNE)

    N = XSecCon * (MD * NTONNES / MW) * L_Run                                      (FASER)

    where N is the total number of expected events, XSecCon is the flux-convoluted trident 
    cross section

    XSecCon = int dPhi/dE xsec(E) dE,

    MD is the detector mass of 1 tonne, N_POT is the protons-on-target 1.1e21, L_Run is the run luminosity 
    in fb^-1, NTONNES is the number of tonnes in the detector (default is 1), NYEAR is the number of years 
    of running (default is 1) and MAr (MW) is the mass of argon (tungsten) in kg. The flag normalized 
    indicates if the provided flux is already normalized (default is False). If True, the number of events
    is calculated by

    N = Phi * XSecCon * MD * N_POT * NTONNES * NYEAR / MAr                         (DUNE)

    where Phi is the energy-integrated flux

    Phi = int dPhi/dE dE,

    given in this case as a user input (Phi=1 by default), and the trident cross section is
    now convoluted with the provided normalized flux

    XSecCon = int 1 / Phi dPhi/dE xsec(E) dE.

    The function also outputs the convoluted cross section normalized by the integrated flux in fb.
    """
    if detector in ['DUNE','MINOS_neutrino','MINOS_antineutrino','MINOS+','INGRID','INGRIDPhase2','SBND']:
        total_flux = simpson(flux, x=flux_energy)                    # [Nv / m^2 POT]
        integrand = np.multiply(flux, xsec_matched)                  # [Nv / GeV POT]
        integral = simpson(integrand, x=flux_energy)                 # [Nv / POT]
        XSecCon = integral / total_flux                              # [m^2]
        if detector == 'DUNE':
            N = integral * (MD * NTONNES / MAr) * (N_POT * NYEAR)                 # [1]
            if correction_factor != None:
                N_bin = [MD * NTONNES / MAr * N_POT * NYEAR * eps * f * sigma * 0.25 for eps, f, sigma in zip(correction_factor, flux, xsec_matched)]
                N = sum(N_bin)
                XSecCon = N / (MD * NTONNES / MAr * N_POT * NYEAR) / total_flux
                return N, XSecCon*1e43, N_bin
            else:
                N_bin = [MD * NTONNES / MAr * N_POT * NYEAR * f * sigma * 0.25 for f, sigma in zip(flux, xsec_matched)] # Events per bin; DUNE histograms have fixed 0.25 GeV energy bin width.            
                return N, XSecCon*1e43, N_bin
        elif detector == 'MINOS_neutrino':
            N = integral * (MD * NTONNES / MFe) * (N_POT_MINOS_neutrino)          # [1]
        elif detector == 'MINOS_antineutrino':
            N = integral * (MD * NTONNES / MFe) * (N_POT_MINOS_antineutrino)      # [1]
        elif detector == 'MINOS+':
            N = integral * (MD * NTONNES / MFe) * (N_POT_MINOSPlus_neutrino)      # [1]
        elif detector == 'INGRID':
            N = integral * (MD * NTONNES / MFe) * (N_POT_INGRID)                  # [1]
        elif detector == 'SBND':
            N = integral * (MD * NTONNES / MAr) * (N_POT_SBND)                    # [1]
        else:
            N = integral * (MD * NTONNES / MFe) * (N_POT_INGRID_phase2)           # [1]
        return N, XSecCon*1e43                                       # [1], [fb]
    if detector in ['FASER','FASER2','SHiP']:
        # Integration variable, flux, is in decreasing order so the integral will run from right to left and give negative. Reverse both lists to get positive, correct result.
        reversed_xsec_matched = list(reversed(xsec_matched))
        reversed_flux = list(reversed(flux))

        total_flux = sum(flux)                                      # [Nv / m^2 fb^-1] or [Nv / m^2 POT]
        integral = simpson(reversed_xsec_matched, x=reversed_flux)  # [Nv / fb^-1] or [Nv / POT]
        XSecCon = integral / total_flux                             # [m^2]
        if detector == 'FASER':
            N = integral * (MD * M_FASERv / MW) * L_Run             # [1]
        if detector == 'FASER2':
            N = integral * (MD * M_FASERv2 / MW) * L_Run2           # [1]
        if detector == 'SHiP':
            N = integral * (MD * M_SHiP / MW) * N_POT_SHiP          # [1]
        return N, XSecCon*1e43                                      # [1], [fb]

def CalculateEventsFASERv(flux, flux_energy, xsec_matched, NTONNES=1, NYEAR=1, detector='FASER'):
    flux = [f*150*25*25*1e-4 for f in flux]  # FASERv flux was normalized by 25x25 cm^2 and 150 fb^-1. We want the time-integrated neutrino flux for FASERv; that is just the v's.
    interaction_prob = [sigma * RHO_W * Length_FASERv / MW for sigma in xsec_matched]
    N_bin = [nu*prob for nu, prob in zip(flux, interaction_prob)]
    return sum(N_bin), 1

### DUNE Standard Flux Events ###
## Coherent ; Argon ##
# Neutrino Mode #
N_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, XSecCon_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, Nbin_1tau_coh_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_coh_Ar_DUNE_neutrino_vmu_67_3, XSecCon_2tau_coh_Ar_DUNE_neutrino_vmu_67_3, Nbin_2tau_coh_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, XSecCon_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, Nbin_2mu_coh_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_coh_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_1tau_coh_Ar_DUNE_neutrino_vmubar_67_3, Nbin_1tau_coh_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_coh_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_2tau_coh_Ar_DUNE_neutrino_vmubar_67_3, Nbin_2tau_coh_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_coh_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_2mu_coh_Ar_DUNE_neutrino_vmubar_67_3, Nbin_2mu_coh_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_coh_Ar_DUNE_neutrino_ve_67_3, XSecCon_1tau_coh_Ar_DUNE_neutrino_ve_67_3, Nbin_1tau_coh_Ar_DUNE_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_neutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_coh_Ar_DUNE_neutrino_vebar_67_3, XSecCon_1tau_coh_Ar_DUNE_neutrino_vebar_67_3, Nbin_1tau_coh_Ar_DUNE_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_neutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)

# EPA
N_EPA_Heaviside_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, XSecCon_EPA_Heaviside_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, Nbin_EPA_Heaviside_1tau_coh_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_EPA_Heaviside_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_EPA_Heaviside_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, XSecCon_EPA_Heaviside_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, Nbin_EPA_Heaviside_2mu_coh_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_EPA_Heaviside_xsec_2mu_coh_Ar_matched, NTONNES=67, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, upper_Nbin_1tau_coh_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_coh_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_neutrino_vmu_67_3, upper_Nbin_2tau_coh_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, upper_Nbin_2mu_coh_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, lower_Nbin_1tau_coh_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_coh_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_neutrino_vmu_67_3, lower_Nbin_2tau_coh_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, lower_Nbin_2mu_coh_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_coh_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_neutrino_vmubar_67_3, upper_Nbin_1tau_coh_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_coh_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_neutrino_vmubar_67_3, upper_Nbin_2tau_coh_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_coh_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_neutrino_vmubar_67_3, upper_Nbin_2mu_coh_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_coh_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_neutrino_vmubar_67_3, lower_Nbin_1tau_coh_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_coh_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_neutrino_vmubar_67_3, lower_Nbin_2tau_coh_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_coh_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_neutrino_vmubar_67_3, lower_Nbin_2mu_coh_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_coh_Ar_DUNE_neutrino_ve_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_neutrino_ve_67_3, upper_Nbin_1tau_coh_Ar_DUNE_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_neutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_coh_Ar_DUNE_neutrino_ve_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_neutrino_ve_67_3, lower_Nbin_1tau_coh_Ar_DUNE_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_neutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_coh_Ar_DUNE_neutrino_vebar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_neutrino_vebar_67_3, upper_Nbin_1tau_coh_Ar_DUNE_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_neutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_coh_Ar_DUNE_neutrino_vebar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_neutrino_vebar_67_3, lower_Nbin_1tau_coh_Ar_DUNE_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_neutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

# Antineutrino Mode #
N_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3, Nbin_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_coh_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_2tau_coh_Ar_DUNE_antineutrino_vmu_67_3, Nbin_2tau_coh_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_coh_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_2mu_coh_Ar_DUNE_antineutrino_vmu_67_3, Nbin_2mu_coh_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, Nbin_1tau_coh_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_2tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, Nbin_2tau_coh_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_coh_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_2mu_coh_Ar_DUNE_antineutrino_vmubar_67_3, Nbin_2mu_coh_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_coh_Ar_DUNE_antineutrino_ve_67_3, XSecCon_1tau_coh_Ar_DUNE_antineutrino_ve_67_3, Nbin_1tau_coh_Ar_DUNE_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_antineutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_coh_Ar_DUNE_antineutrino_vebar_67_3, XSecCon_1tau_coh_Ar_DUNE_antineutrino_vebar_67_3, Nbin_1tau_coh_Ar_DUNE_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3, upper_Nbin_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_coh_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_antineutrino_vmu_67_3, upper_Nbin_2tau_coh_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_coh_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_antineutrino_vmu_67_3, upper_Nbin_2mu_coh_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3, lower_Nbin_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_coh_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_antineutrino_vmu_67_3, lower_Nbin_2tau_coh_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_coh_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_antineutrino_vmu_67_3, lower_Nbin_2mu_coh_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, upper_Nbin_1tau_coh_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, upper_Nbin_2tau_coh_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_coh_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_antineutrino_vmubar_67_3, upper_Nbin_2mu_coh_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, lower_Nbin_1tau_coh_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, lower_Nbin_2tau_coh_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_coh_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_antineutrino_vmubar_67_3, lower_Nbin_2mu_coh_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_coh_Ar_DUNE_antineutrino_ve_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_antineutrino_ve_67_3, upper_Nbin_1tau_coh_Ar_DUNE_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_antineutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_coh_Ar_DUNE_antineutrino_ve_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_antineutrino_ve_67_3, lower_Nbin_1tau_coh_Ar_DUNE_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_antineutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_coh_Ar_DUNE_antineutrino_vebar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vebar_67_3, upper_Nbin_1tau_coh_Ar_DUNE_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_coh_Ar_DUNE_antineutrino_vebar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vebar_67_3, lower_Nbin_1tau_coh_Ar_DUNE_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

## Incoherent ; proton ; Argon ##
# Neutrino Mode #
N_1tau_p_Ar_DUNE_neutrino_vmu_67_3, XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3, Nbin_1tau_p_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_p_Ar_DUNE_neutrino_vmu_67_3, XSecCon_2tau_p_Ar_DUNE_neutrino_vmu_67_3, Nbin_2tau_p_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_p_Ar_DUNE_neutrino_vmu_67_3, XSecCon_2mu_p_Ar_DUNE_neutrino_vmu_67_3, Nbin_2mu_p_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_p_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_1tau_p_Ar_DUNE_neutrino_vmubar_67_3, Nbin_1tau_p_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_p_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_2tau_p_Ar_DUNE_neutrino_vmubar_67_3, Nbin_2tau_p_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_p_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_2mu_p_Ar_DUNE_neutrino_vmubar_67_3, Nbin_2mu_p_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_p_Ar_DUNE_neutrino_ve_67_3, XSecCon_1tau_p_Ar_DUNE_neutrino_ve_67_3, Nbin_1tau_p_Ar_DUNE_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_neutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_p_Ar_DUNE_neutrino_vebar_67_3, XSecCon_1tau_p_Ar_DUNE_neutrino_vebar_67_3, Nbin_1tau_p_Ar_DUNE_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_neutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_p_Ar_matched, NTONNES=67, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_p_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3, upper_Nbin_1tau_p_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_p_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_2tau_p_Ar_DUNE_neutrino_vmu_67_3, upper_Nbin_2tau_p_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_p_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_2mu_p_Ar_DUNE_neutrino_vmu_67_3, upper_Nbin_2mu_p_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_p_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3, lower_Nbin_1tau_p_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_p_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_2tau_p_Ar_DUNE_neutrino_vmu_67_3, lower_Nbin_2tau_p_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_p_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_2mu_p_Ar_DUNE_neutrino_vmu_67_3, lower_Nbin_2mu_p_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_p_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_neutrino_vmubar_67_3, upper_Nbin_1tau_p_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_p_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2tau_p_Ar_DUNE_neutrino_vmubar_67_3, upper_Nbin_2tau_p_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_p_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2mu_p_Ar_DUNE_neutrino_vmubar_67_3, upper_Nbin_2mu_p_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_p_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_neutrino_vmubar_67_3, lower_Nbin_1tau_p_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_p_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2tau_p_Ar_DUNE_neutrino_vmubar_67_3, lower_Nbin_2tau_p_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_p_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2mu_p_Ar_DUNE_neutrino_vmubar_67_3, lower_Nbin_2mu_p_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_p_Ar_DUNE_neutrino_ve_67_3, upper_XSecCon_1tau_p_Ar_DUNE_neutrino_ve_67_3, upper_Nbin_1tau_p_Ar_DUNE_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_neutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_p_Ar_DUNE_neutrino_ve_67_3, lower_XSecCon_1tau_p_Ar_DUNE_neutrino_ve_67_3, lower_Nbin_1tau_p_Ar_DUNE_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_neutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
upper_N_1tau_p_Ar_DUNE_neutrino_vebar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_neutrino_vebar_67_3, upper_Nbin_1tau_p_Ar_DUNE_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_neutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_p_Ar_DUNE_neutrino_vebar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_neutrino_vebar_67_3, lower_Nbin_1tau_p_Ar_DUNE_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_neutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)

# Antineutrino Mode #
N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3, Nbin_1tau_p_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_p_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_2tau_p_Ar_DUNE_antineutrino_vmu_67_3, Nbin_2tau_p_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_p_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_2mu_p_Ar_DUNE_antineutrino_vmu_67_3, Nbin_2mu_p_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3, Nbin_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3, Nbin_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3, Nbin_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_p_Ar_DUNE_antineutrino_ve_67_3, XSecCon_1tau_p_Ar_DUNE_antineutrino_ve_67_3, Nbin_1tau_p_Ar_DUNE_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_antineutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_p_Ar_DUNE_antineutrino_vebar_67_3, XSecCon_1tau_p_Ar_DUNE_antineutrino_vebar_67_3, Nbin_1tau_p_Ar_DUNE_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_p_Ar_matched, NTONNES=67, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3, upper_Nbin_1tau_p_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_p_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2tau_p_Ar_DUNE_antineutrino_vmu_67_3, upper_Nbin_2tau_p_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_p_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2mu_p_Ar_DUNE_antineutrino_vmu_67_3, upper_Nbin_2mu_p_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3, lower_Nbin_1tau_p_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_p_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2tau_p_Ar_DUNE_antineutrino_vmu_67_3, lower_Nbin_2tau_p_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_p_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2mu_p_Ar_DUNE_antineutrino_vmu_67_3, lower_Nbin_2mu_p_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3, upper_Nbin_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3, upper_Nbin_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3, upper_Nbin_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3, lower_Nbin_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3, lower_Nbin_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3, lower_Nbin_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_p_Ar_DUNE_antineutrino_ve_67_3, upper_XSecCon_1tau_p_Ar_DUNE_antineutrino_ve_67_3, upper_Nbin_1tau_p_Ar_DUNE_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_antineutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_p_Ar_DUNE_antineutrino_ve_67_3, lower_XSecCon_1tau_p_Ar_DUNE_antineutrino_ve_67_3, lower_Nbin_1tau_p_Ar_DUNE_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_antineutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
upper_N_1tau_p_Ar_DUNE_antineutrino_vebar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_antineutrino_vebar_67_3, upper_Nbin_1tau_p_Ar_DUNE_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_p_Ar_DUNE_antineutrino_vebar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_antineutrino_vebar_67_3, lower_Nbin_1tau_p_Ar_DUNE_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)

## Incoherent ; neutron ; Argon ##
# Neutrino Mode #
N_1tau_n_Ar_DUNE_neutrino_vmu_67_3, XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3, Nbin_1tau_n_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_n_Ar_DUNE_neutrino_vmu_67_3, XSecCon_2tau_n_Ar_DUNE_neutrino_vmu_67_3, Nbin_2tau_n_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_n_Ar_DUNE_neutrino_vmu_67_3, XSecCon_2mu_n_Ar_DUNE_neutrino_vmu_67_3, Nbin_2mu_n_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_n_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_1tau_n_Ar_DUNE_neutrino_vmubar_67_3, Nbin_1tau_n_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_n_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_2tau_n_Ar_DUNE_neutrino_vmubar_67_3, Nbin_2tau_n_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_n_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_2mu_n_Ar_DUNE_neutrino_vmubar_67_3, Nbin_2mu_n_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_n_Ar_DUNE_neutrino_ve_67_3, XSecCon_1tau_n_Ar_DUNE_neutrino_ve_67_3, Nbin_1tau_n_Ar_DUNE_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_neutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_n_Ar_DUNE_neutrino_vebar_67_3, XSecCon_1tau_n_Ar_DUNE_neutrino_vebar_67_3, Nbin_1tau_n_Ar_DUNE_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_neutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_n_Ar_matched, NTONNES=67, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_n_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3, upper_Nbin_1tau_n_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_n_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_2tau_n_Ar_DUNE_neutrino_vmu_67_3, upper_Nbin_2tau_n_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_n_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_2mu_n_Ar_DUNE_neutrino_vmu_67_3, upper_Nbin_2mu_n_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_n_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3, lower_Nbin_1tau_n_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_n_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_2tau_n_Ar_DUNE_neutrino_vmu_67_3, lower_Nbin_2tau_n_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_n_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_2mu_n_Ar_DUNE_neutrino_vmu_67_3, lower_Nbin_2mu_n_Ar_DUNE_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_n_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_neutrino_vmubar_67_3, upper_Nbin_1tau_n_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_n_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2tau_n_Ar_DUNE_neutrino_vmubar_67_3, upper_Nbin_2tau_n_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_n_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2mu_n_Ar_DUNE_neutrino_vmubar_67_3, upper_Nbin_2mu_n_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_n_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_neutrino_vmubar_67_3, lower_Nbin_1tau_n_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_n_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2tau_n_Ar_DUNE_neutrino_vmubar_67_3, lower_Nbin_2tau_n_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_n_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2mu_n_Ar_DUNE_neutrino_vmubar_67_3, lower_Nbin_2mu_n_Ar_DUNE_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_neutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_n_Ar_DUNE_neutrino_ve_67_3, upper_XSecCon_1tau_n_Ar_DUNE_neutrino_ve_67_3, upper_Nbin_1tau_n_Ar_DUNE_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_neutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_n_Ar_DUNE_neutrino_ve_67_3, lower_XSecCon_1tau_n_Ar_DUNE_neutrino_ve_67_3, lower_Nbin_1tau_n_Ar_DUNE_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_neutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_n_Ar_DUNE_neutrino_vebar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_neutrino_vebar_67_3, upper_Nbin_1tau_n_Ar_DUNE_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_neutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_n_Ar_DUNE_neutrino_vebar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_neutrino_vebar_67_3, lower_Nbin_1tau_n_Ar_DUNE_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_neutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

# Antineutrino Mode #
N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3, Nbin_1tau_n_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_n_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_2tau_n_Ar_DUNE_antineutrino_vmu_67_3, Nbin_2tau_n_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_n_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_2mu_n_Ar_DUNE_antineutrino_vmu_67_3, Nbin_2mu_n_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3, Nbin_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3, Nbin_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3, Nbin_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_n_Ar_DUNE_antineutrino_ve_67_3, XSecCon_1tau_n_Ar_DUNE_antineutrino_ve_67_3, Nbin_1tau_n_Ar_DUNE_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_antineutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_n_Ar_DUNE_antineutrino_vebar_67_3, XSecCon_1tau_n_Ar_DUNE_antineutrino_vebar_67_3, Nbin_1tau_n_Ar_DUNE_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_n_Ar_matched, NTONNES=67, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3, upper_Nbin_1tau_n_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_n_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2tau_n_Ar_DUNE_antineutrino_vmu_67_3, upper_Nbin_2tau_n_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_n_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2mu_n_Ar_DUNE_antineutrino_vmu_67_3, upper_Nbin_2mu_n_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3, lower_Nbin_1tau_n_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_n_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2tau_n_Ar_DUNE_antineutrino_vmu_67_3, lower_Nbin_2tau_n_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_n_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2mu_n_Ar_DUNE_antineutrino_vmu_67_3, lower_Nbin_2mu_n_Ar_DUNE_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3, upper_Nbin_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3, upper_Nbin_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3, upper_Nbin_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3, lower_Nbin_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3, lower_Nbin_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3, lower_Nbin_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vmubar, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_n_Ar_DUNE_antineutrino_ve_67_3, upper_XSecCon_1tau_n_Ar_DUNE_antineutrino_ve_67_3, upper_Nbin_1tau_n_Ar_DUNE_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_antineutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_n_Ar_DUNE_antineutrino_ve_67_3, lower_XSecCon_1tau_n_Ar_DUNE_antineutrino_ve_67_3, lower_Nbin_1tau_n_Ar_DUNE_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_antineutrino_ve, energy_DUNE, DUNE_xsec_ve1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_n_Ar_DUNE_antineutrino_vebar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_antineutrino_vebar_67_3, upper_Nbin_1tau_n_Ar_DUNE_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_n_Ar_DUNE_antineutrino_vebar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_antineutrino_vebar_67_3, lower_Nbin_1tau_n_Ar_DUNE_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_antineutrino_vebar, energy_DUNE, DUNE_xsec_ve1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

## Incoherent ; proton + neutron ; Argon ##
# Neutrino Mode #
N_1tau_incoh_DUNE_neutrino_vmu_67_3, XSecCon_1tau_incoh_DUNE_neutrino_vmu_67_3, Nbin_1tau_incoh_DUNE_neutrino_vmu_67_3 = (N_1tau_p_Ar_DUNE_neutrino_vmu_67_3 + N_1tau_n_Ar_DUNE_neutrino_vmu_67_3), (XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3 + XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_neutrino_vmu_67_3, Nbin_1tau_n_Ar_DUNE_neutrino_vmu_67_3)]
N_2tau_incoh_DUNE_neutrino_vmu_67_3, XSecCon_2tau_incoh_DUNE_neutrino_vmu_67_3, Nbin_2tau_incoh_DUNE_neutrino_vmu_67_3 = (N_2tau_p_Ar_DUNE_neutrino_vmu_67_3 + N_2tau_n_Ar_DUNE_neutrino_vmu_67_3), (XSecCon_2tau_p_Ar_DUNE_neutrino_vmu_67_3 + XSecCon_2tau_n_Ar_DUNE_neutrino_vmu_67_3), [sum(x) for x in zip(Nbin_2tau_p_Ar_DUNE_neutrino_vmu_67_3, Nbin_2tau_n_Ar_DUNE_neutrino_vmu_67_3)]
N_2mu_incoh_DUNE_neutrino_vmu_67_3, XSecCon_2mu_incoh_DUNE_neutrino_vmu_67_3, Nbin_2mu_incoh_DUNE_neutrino_vmu_67_3 = (N_2mu_p_Ar_DUNE_neutrino_vmu_67_3 + N_2mu_n_Ar_DUNE_neutrino_vmu_67_3), (XSecCon_2mu_p_Ar_DUNE_neutrino_vmu_67_3 + XSecCon_2mu_n_Ar_DUNE_neutrino_vmu_67_3), [sum(x) for x in zip(Nbin_2mu_p_Ar_DUNE_neutrino_vmu_67_3, Nbin_2mu_n_Ar_DUNE_neutrino_vmu_67_3)]

N_1tau_incoh_DUNE_neutrino_vmubar_67_3, XSecCon_1tau_incoh_DUNE_neutrino_vmubar_67_3, Nbin_1tau_incoh_DUNE_neutrino_vmubar_67_3 = (N_1tau_p_Ar_DUNE_neutrino_vmubar_67_3 + N_1tau_n_Ar_DUNE_neutrino_vmubar_67_3), (XSecCon_1tau_p_Ar_DUNE_neutrino_vmubar_67_3 + XSecCon_1tau_n_Ar_DUNE_neutrino_vmubar_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_neutrino_vmubar_67_3, Nbin_1tau_n_Ar_DUNE_neutrino_vmubar_67_3)]
N_2tau_incoh_DUNE_neutrino_vmubar_67_3, XSecCon_2tau_incoh_DUNE_neutrino_vmubar_67_3, Nbin_2tau_incoh_DUNE_neutrino_vmubar_67_3 = (N_2tau_p_Ar_DUNE_neutrino_vmubar_67_3 + N_2tau_n_Ar_DUNE_neutrino_vmubar_67_3), (XSecCon_2tau_p_Ar_DUNE_neutrino_vmubar_67_3 + XSecCon_2tau_n_Ar_DUNE_neutrino_vmubar_67_3), [sum(x) for x in zip(Nbin_2tau_p_Ar_DUNE_neutrino_vmubar_67_3, Nbin_2tau_n_Ar_DUNE_neutrino_vmubar_67_3)]
N_2mu_incoh_DUNE_neutrino_vmubar_67_3, XSecCon_2mu_incoh_DUNE_neutrino_vmubar_67_3, Nbin_2mu_incoh_DUNE_neutrino_vmubar_67_3 = (N_2mu_p_Ar_DUNE_neutrino_vmubar_67_3 + N_2mu_n_Ar_DUNE_neutrino_vmubar_67_3), (XSecCon_2mu_p_Ar_DUNE_neutrino_vmubar_67_3 + XSecCon_2mu_n_Ar_DUNE_neutrino_vmubar_67_3), [sum(x) for x in zip(Nbin_2mu_p_Ar_DUNE_neutrino_vmubar_67_3, Nbin_2mu_n_Ar_DUNE_neutrino_vmubar_67_3)]

N_1tau_incoh_DUNE_neutrino_ve_67_3, XSecCon_1tau_incoh_DUNE_neutrino_ve_67_3, Nbin_1tau_incoh_DUNE_neutrino_ve_67_3 = (N_1tau_p_Ar_DUNE_neutrino_ve_67_3 + N_1tau_n_Ar_DUNE_neutrino_ve_67_3), (XSecCon_1tau_p_Ar_DUNE_neutrino_ve_67_3 + XSecCon_1tau_n_Ar_DUNE_neutrino_ve_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_neutrino_ve_67_3, Nbin_1tau_n_Ar_DUNE_neutrino_ve_67_3)]
N_1tau_incoh_DUNE_neutrino_vebar_67_3, XSecCon_1tau_incoh_DUNE_neutrino_vebar_67_3, Nbin_1tau_incoh_DUNE_neutrino_vebar_67_3 = (N_1tau_p_Ar_DUNE_neutrino_vebar_67_3 + N_1tau_n_Ar_DUNE_neutrino_vebar_67_3), (XSecCon_1tau_p_Ar_DUNE_neutrino_vebar_67_3 + XSecCon_1tau_n_Ar_DUNE_neutrino_vebar_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_neutrino_vebar_67_3, Nbin_1tau_n_Ar_DUNE_neutrino_vebar_67_3)]

# Upper and Lower limits #
upper_N_1tau_incoh_DUNE_neutrino_vmu_67_3, upper_XSecCon_1tau_incoh_DUNE_neutrino_vmu_67_3, upper_Nbin_1tau_incoh_DUNE_neutrino_vmu_67_3 = (upper_N_1tau_p_Ar_DUNE_neutrino_vmu_67_3 + upper_N_1tau_n_Ar_DUNE_neutrino_vmu_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_neutrino_vmu_67_3, upper_Nbin_1tau_n_Ar_DUNE_neutrino_vmu_67_3)]
upper_N_2tau_incoh_DUNE_neutrino_vmu_67_3, upper_XSecCon_2tau_incoh_DUNE_neutrino_vmu_67_3, upper_Nbin_2tau_incoh_DUNE_neutrino_vmu_67_3 = (upper_N_2tau_p_Ar_DUNE_neutrino_vmu_67_3 + upper_N_2tau_n_Ar_DUNE_neutrino_vmu_67_3), (upper_XSecCon_2tau_p_Ar_DUNE_neutrino_vmu_67_3 + upper_XSecCon_2tau_n_Ar_DUNE_neutrino_vmu_67_3), [sum(x) for x in zip(upper_Nbin_2tau_p_Ar_DUNE_neutrino_vmu_67_3, upper_Nbin_2tau_n_Ar_DUNE_neutrino_vmu_67_3)]
upper_N_2mu_incoh_DUNE_neutrino_vmu_67_3, upper_XSecCon_2mu_incoh_DUNE_neutrino_vmu_67_3, upper_Nbin_2mu_incoh_DUNE_neutrino_vmu_67_3 = (upper_N_2mu_p_Ar_DUNE_neutrino_vmu_67_3 + upper_N_2mu_n_Ar_DUNE_neutrino_vmu_67_3), (upper_XSecCon_2mu_p_Ar_DUNE_neutrino_vmu_67_3 + upper_XSecCon_2mu_n_Ar_DUNE_neutrino_vmu_67_3), [sum(x) for x in zip(upper_Nbin_2mu_p_Ar_DUNE_neutrino_vmu_67_3, upper_Nbin_2mu_n_Ar_DUNE_neutrino_vmu_67_3)]

lower_N_1tau_incoh_DUNE_neutrino_vmu_67_3, lower_XSecCon_1tau_incoh_DUNE_neutrino_vmu_67_3, lower_Nbin_1tau_incoh_DUNE_neutrino_vmu_67_3 = (lower_N_1tau_p_Ar_DUNE_neutrino_vmu_67_3 + lower_N_1tau_n_Ar_DUNE_neutrino_vmu_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_neutrino_vmu_67_3, lower_Nbin_1tau_n_Ar_DUNE_neutrino_vmu_67_3)]
lower_N_2tau_incoh_DUNE_neutrino_vmu_67_3, lower_XSecCon_2tau_incoh_DUNE_neutrino_vmu_67_3, lower_Nbin_2tau_incoh_DUNE_neutrino_vmu_67_3 = (lower_N_2tau_p_Ar_DUNE_neutrino_vmu_67_3 + lower_N_2tau_n_Ar_DUNE_neutrino_vmu_67_3), (lower_XSecCon_2tau_p_Ar_DUNE_neutrino_vmu_67_3 + lower_XSecCon_2tau_n_Ar_DUNE_neutrino_vmu_67_3), [sum(x) for x in zip(lower_Nbin_2tau_p_Ar_DUNE_neutrino_vmu_67_3, lower_Nbin_2tau_n_Ar_DUNE_neutrino_vmu_67_3)]
lower_N_2mu_incoh_DUNE_neutrino_vmu_67_3, lower_XSecCon_2mu_incoh_DUNE_neutrino_vmu_67_3, lower_Nbin_2mu_incoh_DUNE_neutrino_vmu_67_3 = (lower_N_2mu_p_Ar_DUNE_neutrino_vmu_67_3 + lower_N_2mu_n_Ar_DUNE_neutrino_vmu_67_3), (lower_XSecCon_2mu_p_Ar_DUNE_neutrino_vmu_67_3 + lower_XSecCon_2mu_n_Ar_DUNE_neutrino_vmu_67_3), [sum(x) for x in zip(lower_Nbin_2mu_p_Ar_DUNE_neutrino_vmu_67_3, lower_Nbin_2mu_n_Ar_DUNE_neutrino_vmu_67_3)]

upper_N_1tau_incoh_DUNE_neutrino_vmubar_67_3, upper_XSecCon_1tau_incoh_DUNE_neutrino_vmubar_67_3, upper_Nbin_1tau_incoh_DUNE_neutrino_vmubar_67_3 = (upper_N_1tau_p_Ar_DUNE_neutrino_vmubar_67_3 + upper_N_1tau_n_Ar_DUNE_neutrino_vmubar_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_neutrino_vmubar_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_neutrino_vmubar_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_neutrino_vmubar_67_3, upper_Nbin_1tau_n_Ar_DUNE_neutrino_vmubar_67_3)]
upper_N_2tau_incoh_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2tau_incoh_DUNE_neutrino_vmubar_67_3, upper_Nbin_2tau_incoh_DUNE_neutrino_vmubar_67_3 = (upper_N_2tau_p_Ar_DUNE_neutrino_vmubar_67_3 + upper_N_2tau_n_Ar_DUNE_neutrino_vmubar_67_3), (upper_XSecCon_2tau_p_Ar_DUNE_neutrino_vmubar_67_3 + upper_XSecCon_2tau_n_Ar_DUNE_neutrino_vmubar_67_3), [sum(x) for x in zip(upper_Nbin_2tau_p_Ar_DUNE_neutrino_vmubar_67_3, upper_Nbin_2tau_n_Ar_DUNE_neutrino_vmubar_67_3)]
upper_N_2mu_incoh_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2mu_incoh_DUNE_neutrino_vmubar_67_3, upper_Nbin_2mu_incoh_DUNE_neutrino_vmubar_67_3 = (upper_N_2mu_p_Ar_DUNE_neutrino_vmubar_67_3 + upper_N_2mu_n_Ar_DUNE_neutrino_vmubar_67_3), (upper_XSecCon_2mu_p_Ar_DUNE_neutrino_vmubar_67_3 + upper_XSecCon_2mu_n_Ar_DUNE_neutrino_vmubar_67_3), [sum(x) for x in zip(upper_Nbin_2mu_p_Ar_DUNE_neutrino_vmubar_67_3, upper_Nbin_2mu_n_Ar_DUNE_neutrino_vmubar_67_3)]

lower_N_1tau_incoh_DUNE_neutrino_vmubar_67_3, lower_XSecCon_1tau_incoh_DUNE_neutrino_vmubar_67_3, lower_Nbin_1tau_incoh_DUNE_neutrino_vmubar_67_3 = (lower_N_1tau_p_Ar_DUNE_neutrino_vmubar_67_3 + lower_N_1tau_n_Ar_DUNE_neutrino_vmubar_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_neutrino_vmubar_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_neutrino_vmubar_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_neutrino_vmubar_67_3, lower_Nbin_1tau_n_Ar_DUNE_neutrino_vmubar_67_3)]
lower_N_2tau_incoh_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2tau_incoh_DUNE_neutrino_vmubar_67_3, lower_Nbin_2tau_incoh_DUNE_neutrino_vmubar_67_3 = (lower_N_2tau_p_Ar_DUNE_neutrino_vmubar_67_3 + lower_N_2tau_n_Ar_DUNE_neutrino_vmubar_67_3), (lower_XSecCon_2tau_p_Ar_DUNE_neutrino_vmubar_67_3 + lower_XSecCon_2tau_n_Ar_DUNE_neutrino_vmubar_67_3), [sum(x) for x in zip(lower_Nbin_2tau_p_Ar_DUNE_neutrino_vmubar_67_3, lower_Nbin_2tau_n_Ar_DUNE_neutrino_vmubar_67_3)]
lower_N_2mu_incoh_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2mu_incoh_DUNE_neutrino_vmubar_67_3, lower_Nbin_2mu_incoh_DUNE_neutrino_vmubar_67_3 = (lower_N_2mu_p_Ar_DUNE_neutrino_vmubar_67_3 + lower_N_2mu_n_Ar_DUNE_neutrino_vmubar_67_3), (lower_XSecCon_2mu_p_Ar_DUNE_neutrino_vmubar_67_3 + lower_XSecCon_2mu_n_Ar_DUNE_neutrino_vmubar_67_3), [sum(x) for x in zip(lower_Nbin_2mu_p_Ar_DUNE_neutrino_vmubar_67_3, lower_Nbin_2mu_n_Ar_DUNE_neutrino_vmubar_67_3)]

upper_N_1tau_incoh_DUNE_neutrino_ve_67_3, upper_XSecCon_1tau_incoh_DUNE_neutrino_ve_67_3, upper_Nbin_1tau_incoh_DUNE_neutrino_ve_67_3 = (upper_N_1tau_p_Ar_DUNE_neutrino_ve_67_3 + upper_N_1tau_n_Ar_DUNE_neutrino_ve_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_neutrino_ve_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_neutrino_ve_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_neutrino_ve_67_3, upper_Nbin_1tau_n_Ar_DUNE_neutrino_ve_67_3)]
lower_N_1tau_incoh_DUNE_neutrino_ve_67_3, lower_XSecCon_1tau_incoh_DUNE_neutrino_ve_67_3, lower_Nbin_1tau_incoh_DUNE_neutrino_ve_67_3 = (lower_N_1tau_p_Ar_DUNE_neutrino_ve_67_3 + lower_N_1tau_n_Ar_DUNE_neutrino_ve_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_neutrino_ve_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_neutrino_ve_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_neutrino_ve_67_3, lower_Nbin_1tau_n_Ar_DUNE_neutrino_ve_67_3)]

upper_N_1tau_incoh_DUNE_neutrino_vebar_67_3, upper_XSecCon_1tau_incoh_DUNE_neutrino_vebar_67_3, upper_Nbin_1tau_incoh_DUNE_neutrino_vebar_67_3 = (upper_N_1tau_p_Ar_DUNE_neutrino_vebar_67_3 + upper_N_1tau_n_Ar_DUNE_neutrino_vebar_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_neutrino_vebar_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_neutrino_vebar_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_neutrino_vebar_67_3, upper_Nbin_1tau_n_Ar_DUNE_neutrino_vebar_67_3)]
lower_N_1tau_incoh_DUNE_neutrino_vebar_67_3, lower_XSecCon_1tau_incoh_DUNE_neutrino_vebar_67_3, lower_Nbin_1tau_incoh_DUNE_neutrino_vebar_67_3 = (lower_N_1tau_p_Ar_DUNE_neutrino_vebar_67_3 + lower_N_1tau_n_Ar_DUNE_neutrino_vebar_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_neutrino_vebar_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_neutrino_vebar_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_neutrino_vebar_67_3, lower_Nbin_1tau_n_Ar_DUNE_neutrino_vebar_67_3)]

# Antineutrino Mode #
N_1tau_incoh_DUNE_antineutrino_vmu_67_3, XSecCon_1tau_incoh_DUNE_antineutrino_vmu_67_3, Nbin_1tau_incoh_DUNE_antineutrino_vmu_67_3 = (N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3 + N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3), (XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3 + XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_antineutrino_vmu_67_3, Nbin_1tau_n_Ar_DUNE_antineutrino_vmu_67_3)]
N_2tau_incoh_DUNE_antineutrino_vmu_67_3, XSecCon_2tau_incoh_DUNE_antineutrino_vmu_67_3, Nbin_2tau_incoh_DUNE_antineutrino_vmu_67_3 = (N_2tau_p_Ar_DUNE_antineutrino_vmu_67_3 + N_2tau_n_Ar_DUNE_antineutrino_vmu_67_3), (XSecCon_2tau_p_Ar_DUNE_antineutrino_vmu_67_3 + XSecCon_2tau_n_Ar_DUNE_antineutrino_vmu_67_3), [sum(x) for x in zip(Nbin_2tau_p_Ar_DUNE_antineutrino_vmu_67_3, Nbin_2tau_n_Ar_DUNE_antineutrino_vmu_67_3)]
N_2mu_incoh_DUNE_antineutrino_vmu_67_3, XSecCon_2mu_incoh_DUNE_antineutrino_vmu_67_3, Nbin_2mu_incoh_DUNE_antineutrino_vmu_67_3 = (N_2mu_p_Ar_DUNE_antineutrino_vmu_67_3 + N_2mu_n_Ar_DUNE_antineutrino_vmu_67_3), (XSecCon_2mu_p_Ar_DUNE_antineutrino_vmu_67_3 + XSecCon_2mu_n_Ar_DUNE_antineutrino_vmu_67_3), [sum(x) for x in zip(Nbin_2mu_p_Ar_DUNE_antineutrino_vmu_67_3, Nbin_2mu_n_Ar_DUNE_antineutrino_vmu_67_3)]

N_1tau_incoh_DUNE_antineutrino_vmubar_67_3, XSecCon_1tau_incoh_DUNE_antineutrino_vmubar_67_3, Nbin_1tau_incoh_DUNE_antineutrino_vmubar_67_3 = (N_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3 + N_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3), (XSecCon_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3 + XSecCon_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3, Nbin_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3)]
N_2tau_incoh_DUNE_antineutrino_vmubar_67_3, XSecCon_2tau_incoh_DUNE_antineutrino_vmubar_67_3, Nbin_2tau_incoh_DUNE_antineutrino_vmubar_67_3 = (N_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3 + N_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3), (XSecCon_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3 + XSecCon_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3), [sum(x) for x in zip(Nbin_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3, Nbin_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3)]
N_2mu_incoh_DUNE_antineutrino_vmubar_67_3, XSecCon_2mu_incoh_DUNE_antineutrino_vmubar_67_3, Nbin_2mu_incoh_DUNE_antineutrino_vmubar_67_3 = (N_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3 + N_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3), (XSecCon_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3 + XSecCon_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3), [sum(x) for x in zip(Nbin_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3, Nbin_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3)]

N_1tau_incoh_DUNE_antineutrino_ve_67_3, XSecCon_1tau_incoh_DUNE_antineutrino_ve_67_3, Nbin_1tau_incoh_DUNE_antineutrino_ve_67_3 = (N_1tau_p_Ar_DUNE_antineutrino_ve_67_3 + N_1tau_n_Ar_DUNE_antineutrino_ve_67_3), (XSecCon_1tau_p_Ar_DUNE_antineutrino_ve_67_3 + XSecCon_1tau_n_Ar_DUNE_antineutrino_ve_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_antineutrino_ve_67_3, Nbin_1tau_n_Ar_DUNE_antineutrino_ve_67_3)]
N_1tau_incoh_DUNE_antineutrino_vebar_67_3, XSecCon_1tau_incoh_DUNE_antineutrino_vebar_67_3, Nbin_1tau_incoh_DUNE_antineutrino_vebar_67_3 = (N_1tau_p_Ar_DUNE_antineutrino_vebar_67_3 + N_1tau_n_Ar_DUNE_antineutrino_vebar_67_3), (XSecCon_1tau_p_Ar_DUNE_antineutrino_vebar_67_3 + XSecCon_1tau_n_Ar_DUNE_antineutrino_vebar_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_antineutrino_vebar_67_3, Nbin_1tau_n_Ar_DUNE_antineutrino_vebar_67_3)]

# Upper and Lower limits #
upper_N_1tau_incoh_DUNE_antineutrino_vmu_67_3, upper_XSecCon_1tau_incoh_DUNE_antineutrino_vmu_67_3, upper_Nbin_1tau_incoh_DUNE_antineutrino_vmu_67_3 = (upper_N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3 + upper_N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_antineutrino_vmu_67_3, upper_Nbin_1tau_n_Ar_DUNE_antineutrino_vmu_67_3)]
upper_N_2tau_incoh_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2tau_incoh_DUNE_antineutrino_vmu_67_3, upper_Nbin_2tau_incoh_DUNE_antineutrino_vmu_67_3 = (upper_N_2tau_p_Ar_DUNE_antineutrino_vmu_67_3 + upper_N_2tau_n_Ar_DUNE_antineutrino_vmu_67_3), (upper_XSecCon_2tau_p_Ar_DUNE_antineutrino_vmu_67_3 + upper_XSecCon_2tau_n_Ar_DUNE_antineutrino_vmu_67_3), [sum(x) for x in zip(upper_Nbin_2tau_p_Ar_DUNE_antineutrino_vmu_67_3, upper_Nbin_2tau_n_Ar_DUNE_antineutrino_vmu_67_3)]
upper_N_2mu_incoh_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2mu_incoh_DUNE_antineutrino_vmu_67_3, upper_Nbin_2mu_incoh_DUNE_antineutrino_vmu_67_3 = (upper_N_2mu_p_Ar_DUNE_antineutrino_vmu_67_3 + upper_N_2mu_n_Ar_DUNE_antineutrino_vmu_67_3), (upper_XSecCon_2mu_p_Ar_DUNE_antineutrino_vmu_67_3 + upper_XSecCon_2mu_n_Ar_DUNE_antineutrino_vmu_67_3), [sum(x) for x in zip(upper_Nbin_2mu_p_Ar_DUNE_antineutrino_vmu_67_3, upper_Nbin_2mu_n_Ar_DUNE_antineutrino_vmu_67_3)]

lower_N_1tau_incoh_DUNE_antineutrino_vmu_67_3, lower_XSecCon_1tau_incoh_DUNE_antineutrino_vmu_67_3, lower_Nbin_1tau_incoh_DUNE_antineutrino_vmu_67_3 = (lower_N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3 + lower_N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_antineutrino_vmu_67_3, lower_Nbin_1tau_n_Ar_DUNE_antineutrino_vmu_67_3)]
lower_N_2tau_incoh_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2tau_incoh_DUNE_antineutrino_vmu_67_3, lower_Nbin_2tau_incoh_DUNE_antineutrino_vmu_67_3 = (lower_N_2tau_p_Ar_DUNE_antineutrino_vmu_67_3 + lower_N_2tau_n_Ar_DUNE_antineutrino_vmu_67_3), (lower_XSecCon_2tau_p_Ar_DUNE_antineutrino_vmu_67_3 + lower_XSecCon_2tau_n_Ar_DUNE_antineutrino_vmu_67_3), [sum(x) for x in zip(lower_Nbin_2tau_p_Ar_DUNE_antineutrino_vmu_67_3, lower_Nbin_2tau_n_Ar_DUNE_antineutrino_vmu_67_3)]
lower_N_2mu_incoh_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2mu_incoh_DUNE_antineutrino_vmu_67_3, lower_Nbin_2mu_incoh_DUNE_antineutrino_vmu_67_3 = (lower_N_2mu_p_Ar_DUNE_antineutrino_vmu_67_3 + lower_N_2mu_n_Ar_DUNE_antineutrino_vmu_67_3), (lower_XSecCon_2mu_p_Ar_DUNE_antineutrino_vmu_67_3 + lower_XSecCon_2mu_n_Ar_DUNE_antineutrino_vmu_67_3), [sum(x) for x in zip(lower_Nbin_2mu_p_Ar_DUNE_antineutrino_vmu_67_3, lower_Nbin_2mu_n_Ar_DUNE_antineutrino_vmu_67_3)]

upper_N_1tau_incoh_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_1tau_incoh_DUNE_antineutrino_vmubar_67_3, upper_Nbin_1tau_incoh_DUNE_antineutrino_vmubar_67_3 = (upper_N_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3 + upper_N_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3, upper_Nbin_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3)]
upper_N_2tau_incoh_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2tau_incoh_DUNE_antineutrino_vmubar_67_3, upper_Nbin_2tau_incoh_DUNE_antineutrino_vmubar_67_3 = (upper_N_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3 + upper_N_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3), (upper_XSecCon_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3 + upper_XSecCon_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3), [sum(x) for x in zip(upper_Nbin_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3, upper_Nbin_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3)]
upper_N_2mu_incoh_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2mu_incoh_DUNE_antineutrino_vmubar_67_3, upper_Nbin_2mu_incoh_DUNE_antineutrino_vmubar_67_3 = (upper_N_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3 + upper_N_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3), (upper_XSecCon_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3 + upper_XSecCon_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3), [sum(x) for x in zip(upper_Nbin_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3, upper_Nbin_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3)]

lower_N_1tau_incoh_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_1tau_incoh_DUNE_antineutrino_vmubar_67_3, lower_Nbin_1tau_incoh_DUNE_antineutrino_vmubar_67_3 = (lower_N_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3 + lower_N_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3, lower_Nbin_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3)]
lower_N_2tau_incoh_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2tau_incoh_DUNE_antineutrino_vmubar_67_3, lower_Nbin_2tau_incoh_DUNE_antineutrino_vmubar_67_3 = (lower_N_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3 + lower_N_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3), (lower_XSecCon_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3 + lower_XSecCon_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3), [sum(x) for x in zip(lower_Nbin_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3, lower_Nbin_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3)]
lower_N_2mu_incoh_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2mu_incoh_DUNE_antineutrino_vmubar_67_3, lower_Nbin_2mu_incoh_DUNE_antineutrino_vmubar_67_3 = (lower_N_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3 + lower_N_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3), (lower_XSecCon_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3 + lower_XSecCon_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3), [sum(x) for x in zip(lower_Nbin_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3, lower_Nbin_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3)]

upper_N_1tau_incoh_DUNE_antineutrino_ve_67_3, upper_XSecCon_1tau_incoh_DUNE_antineutrino_ve_67_3, upper_Nbin_1tau_incoh_DUNE_antineutrino_ve_67_3 = (upper_N_1tau_p_Ar_DUNE_antineutrino_ve_67_3 + upper_N_1tau_n_Ar_DUNE_antineutrino_ve_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_antineutrino_ve_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_antineutrino_ve_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_antineutrino_ve_67_3, upper_Nbin_1tau_n_Ar_DUNE_antineutrino_ve_67_3)]
lower_N_1tau_incoh_DUNE_antineutrino_ve_67_3, lower_XSecCon_1tau_incoh_DUNE_antineutrino_ve_67_3, lower_Nbin_1tau_incoh_DUNE_antineutrino_ve_67_3 = (lower_N_1tau_p_Ar_DUNE_antineutrino_ve_67_3 + lower_N_1tau_n_Ar_DUNE_antineutrino_ve_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_antineutrino_ve_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_antineutrino_ve_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_antineutrino_ve_67_3, lower_Nbin_1tau_n_Ar_DUNE_antineutrino_ve_67_3)]

upper_N_1tau_incoh_DUNE_antineutrino_vebar_67_3, upper_XSecCon_1tau_incoh_DUNE_antineutrino_vebar_67_3, upper_Nbin_1tau_incoh_DUNE_antineutrino_vebar_67_3 = (upper_N_1tau_p_Ar_DUNE_antineutrino_vebar_67_3 + upper_N_1tau_n_Ar_DUNE_antineutrino_vebar_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_antineutrino_vebar_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_antineutrino_vebar_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_antineutrino_vebar_67_3, upper_Nbin_1tau_n_Ar_DUNE_antineutrino_vebar_67_3)]
lower_N_1tau_incoh_DUNE_antineutrino_vebar_67_3, lower_XSecCon_1tau_incoh_DUNE_antineutrino_vebar_67_3, lower_Nbin_1tau_incoh_DUNE_antineutrino_vebar_67_3 = (lower_N_1tau_p_Ar_DUNE_antineutrino_vebar_67_3 + lower_N_1tau_n_Ar_DUNE_antineutrino_vebar_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_antineutrino_vebar_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_antineutrino_vebar_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_antineutrino_vebar_67_3, lower_Nbin_1tau_n_Ar_DUNE_antineutrino_vebar_67_3)]


### DUNE Tau-Optimized ND Flux Events ###
## Coherent ; Argon ##
# Neutrino Mode #
N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, Nbin_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, Nbin_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_coh_Ar_DUNE_tau_opt_neutrino_ve_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_ve_67_3, Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vebar_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vebar_67_3, Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)

# EPA
N_EPA_Heaviside_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_EPA_Heaviside_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_EPA_Heaviside_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_EPA_Heaviside_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_EPA_Heaviside_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_EPA_Heaviside_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_EPA_Heaviside_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_EPA_Heaviside_xsec_2mu_coh_Ar_matched, NTONNES=67, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_Nbin_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_Nbin_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_Nbin_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_Nbin_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_Nbin_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_Nbin_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_Nbin_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_Nbin_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_ve_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_ve_67_3, upper_Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_ve_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_ve_67_3, lower_Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vebar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vebar_67_3, upper_Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vebar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vebar_67_3, lower_Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

# Antineutrino Mode #
N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, Nbin_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, Nbin_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, Nbin_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, Nbin_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_ve_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_ve_67_3, Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_Nbin_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_Nbin_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_Nbin_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_Nbin_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_Nbin_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_Nbin_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_Nbin_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_Nbin_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_ve_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_ve_67_3, upper_Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_ve_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_ve_67_3, lower_Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, upper_Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, lower_Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)

## Incoherent ; proton ; Argon ##
# Neutrino Mode #
N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, Nbin_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, Nbin_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3, Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3, Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_p_Ar_matched, NTONNES=67, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_Nbin_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_Nbin_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_Nbin_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_Nbin_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_Nbin_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_Nbin_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_Nbin_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_Nbin_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3, upper_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3, lower_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3, upper_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3, lower_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)

# Antineutrino Mode #
N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, Nbin_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, Nbin_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, Nbin_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, Nbin_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3, Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_p_Ar_matched, NTONNES=67, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_Nbin_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_Nbin_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_Nbin_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_Nbin_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_Nbin_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_Nbin_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_Nbin_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_Nbin_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3, upper_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3, lower_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, upper_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, lower_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)

## Incoherent ; neutron ; Argon ##
# Neutrino Mode #
N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, Nbin_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, Nbin_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_n_Ar_matched, NTONNES=67, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_Nbin_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_Nbin_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_Nbin_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_Nbin_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_Nbin_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_Nbin_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_Nbin_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_Nbin_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_neutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

# Antineutrino Mode #
N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, Nbin_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, Nbin_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, Nbin_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, Nbin_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched, NTONNES=67, NYEAR=3)

N_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_n_Ar_matched, NTONNES=67, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_Nbin_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_Nbin_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_Nbin_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_Nbin_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_Nbin_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_Nbin_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched_upper, NTONNES=67, NYEAR=3)

lower_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_Nbin_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_Nbin_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmubar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_ve, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

upper_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
lower_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3 = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vebar, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_ve1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)

## Incoherent ; proton + neutron ; Argon ##
# Neutrino Mode #
N_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3 = (N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), (XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3)]
N_2tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_2tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3 = (N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), (XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), [sum(x) for x in zip(Nbin_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3)]
N_2mu_incoh_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2mu_incoh_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_2mu_incoh_DUNE_tau_opt_neutrino_vmu_67_3 = (N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), (XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), [sum(x) for x in zip(Nbin_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, Nbin_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3)]

N_1tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3 = (N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), (XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)]
N_2tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, Nbin_2tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3 = (N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), (XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), [sum(x) for x in zip(Nbin_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, Nbin_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)]
N_2mu_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2mu_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, Nbin_2mu_incoh_DUNE_tau_opt_neutrino_vmubar_67_3 = (N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), (XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), [sum(x) for x in zip(Nbin_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, Nbin_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)]

N_1tau_incoh_DUNE_tau_opt_neutrino_ve_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_ve_67_3, Nbin_1tau_incoh_DUNE_tau_opt_neutrino_ve_67_3 = (N_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3 + N_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3), (XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3 + XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3)]

N_1tau_incoh_DUNE_tau_opt_neutrino_vebar_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vebar_67_3, Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vebar_67_3 = (N_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3 + N_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3), (XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3 + XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3)]

# Upper and Lower limits #
upper_N_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, upper_Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3 = (upper_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + upper_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3)]
upper_N_2tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, upper_Nbin_2tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3 = (upper_N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + upper_N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), (upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), [sum(x) for x in zip(upper_Nbin_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_Nbin_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3)]
upper_N_2mu_incoh_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2mu_incoh_DUNE_tau_opt_neutrino_vmu_67_3, upper_Nbin_2mu_incoh_DUNE_tau_opt_neutrino_vmu_67_3 = (upper_N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + upper_N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), (upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), [sum(x) for x in zip(upper_Nbin_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_Nbin_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3)]

lower_N_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, lower_Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3 = (lower_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + lower_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3)]
lower_N_2tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, lower_Nbin_2tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3 = (lower_N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + lower_N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), (lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), [sum(x) for x in zip(lower_Nbin_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_Nbin_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3)]
lower_N_2mu_incoh_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2mu_incoh_DUNE_tau_opt_neutrino_vmu_67_3, lower_Nbin_2mu_incoh_DUNE_tau_opt_neutrino_vmu_67_3 = (lower_N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + lower_N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), (lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3 + lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3), [sum(x) for x in zip(lower_Nbin_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_Nbin_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3)]

upper_N_1tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, upper_Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3 = (upper_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + upper_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)]
upper_N_2tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, upper_Nbin_2tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3 = (upper_N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + upper_N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), (upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), [sum(x) for x in zip(upper_Nbin_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_Nbin_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)]
upper_N_2mu_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2mu_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, upper_Nbin_2mu_incoh_DUNE_tau_opt_neutrino_vmubar_67_3 = (upper_N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + upper_N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), (upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), [sum(x) for x in zip(upper_Nbin_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_Nbin_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)]

lower_N_1tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, lower_Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3 = (lower_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + lower_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)]
lower_N_2tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, lower_Nbin_2tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3 = (lower_N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + lower_N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), (lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), [sum(x) for x in zip(lower_Nbin_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_Nbin_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)]
lower_N_2mu_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2mu_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, lower_Nbin_2mu_incoh_DUNE_tau_opt_neutrino_vmubar_67_3 = (lower_N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + lower_N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), (lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3 + lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3), [sum(x) for x in zip(lower_Nbin_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_Nbin_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)]

upper_N_1tau_incoh_DUNE_tau_opt_neutrino_ve_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_ve_67_3, upper_Nbin_1tau_incoh_DUNE_tau_opt_neutrino_ve_67_3 = (upper_N_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3 + upper_N_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3)]
lower_N_1tau_incoh_DUNE_tau_opt_neutrino_ve_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_ve_67_3, lower_Nbin_1tau_incoh_DUNE_tau_opt_neutrino_ve_67_3 = (lower_N_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3 + lower_N_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3)]

upper_N_1tau_incoh_DUNE_tau_opt_neutrino_vebar_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vebar_67_3, upper_Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vebar_67_3 = (upper_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3 + upper_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3)]
lower_N_1tau_incoh_DUNE_tau_opt_neutrino_vebar_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vebar_67_3, lower_Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vebar_67_3 = (lower_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3 + lower_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3)]

# Antineutrino Mode #
N_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3 = (N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), (XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)]
N_2tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, Nbin_2tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3 = (N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), (XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), [sum(x) for x in zip(Nbin_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, Nbin_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)]
N_2mu_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2mu_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, Nbin_2mu_incoh_DUNE_tau_opt_antineutrino_vmu_67_3 = (N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), (XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), [sum(x) for x in zip(Nbin_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, Nbin_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)]

N_1tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3 = (N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), (XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)]
N_2tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, Nbin_2tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3 = (N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), (XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), [sum(x) for x in zip(Nbin_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, Nbin_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)]
N_2mu_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2mu_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, Nbin_2mu_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3 = (N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), (XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), [sum(x) for x in zip(Nbin_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, Nbin_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)]

N_1tau_incoh_DUNE_tau_opt_antineutrino_ve_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_ve_67_3, Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_ve_67_3 = (N_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3 + N_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3), (XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3 + XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3)]

N_1tau_incoh_DUNE_tau_opt_antineutrino_vebar_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vebar_67_3, Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vebar_67_3 = (N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3 + N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3), (XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3 + XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3)]

# Upper and Lower limits #
upper_N_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, upper_Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3 = (upper_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + upper_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)]
upper_N_2tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, upper_Nbin_2tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3 = (upper_N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + upper_N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), (upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), [sum(x) for x in zip(upper_Nbin_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_Nbin_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)]
upper_N_2mu_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2mu_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, upper_Nbin_2mu_incoh_DUNE_tau_opt_antineutrino_vmu_67_3 = (upper_N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + upper_N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), (upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), [sum(x) for x in zip(upper_Nbin_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_Nbin_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)]

lower_N_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, lower_Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3 = (lower_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + lower_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)]
lower_N_2tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, lower_Nbin_2tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3 = (lower_N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + lower_N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), (lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), [sum(x) for x in zip(lower_Nbin_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_Nbin_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)]
lower_N_2mu_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2mu_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, lower_Nbin_2mu_incoh_DUNE_tau_opt_antineutrino_vmu_67_3 = (lower_N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + lower_N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), (lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3 + lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3), [sum(x) for x in zip(lower_Nbin_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_Nbin_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)]

upper_N_1tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3 = (upper_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + upper_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)]
upper_N_2tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_Nbin_2tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3 = (upper_N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + upper_N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), (upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), [sum(x) for x in zip(upper_Nbin_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_Nbin_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)]
upper_N_2mu_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2mu_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_Nbin_2mu_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3 = (upper_N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + upper_N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), (upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), [sum(x) for x in zip(upper_Nbin_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_Nbin_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)]

lower_N_1tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3 = (lower_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + lower_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)]
lower_N_2tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_Nbin_2tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3 = (lower_N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + lower_N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), (lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), [sum(x) for x in zip(lower_Nbin_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_Nbin_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)]
lower_N_2mu_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2mu_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_Nbin_2mu_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3 = (lower_N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + lower_N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), (lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3 + lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3), [sum(x) for x in zip(lower_Nbin_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_Nbin_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)]

upper_N_1tau_incoh_DUNE_tau_opt_antineutrino_ve_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_ve_67_3, upper_Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_ve_67_3 = (upper_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3 + upper_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3)]
lower_N_1tau_incoh_DUNE_tau_opt_antineutrino_ve_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_ve_67_3, lower_Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_ve_67_3 = (lower_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3 + lower_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3)]

upper_N_1tau_incoh_DUNE_tau_opt_antineutrino_vebar_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vebar_67_3, upper_Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vebar_67_3 = (upper_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3 + upper_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3)]
lower_N_1tau_incoh_DUNE_tau_opt_antineutrino_vebar_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vebar_67_3, lower_Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vebar_67_3 = (lower_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3 + lower_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3)]

### FASER vmu Flux Events ###
## Coherent ; Tungsten ##
# vmu Flux #
N_1tau_coh_W_FASERvmu_1p2_3, XSecCon_1tau_coh_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_coh_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')
N_2tau_coh_W_FASERvmu_1p2_3, XSecCon_2tau_coh_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_coh_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')
N_2mu_coh_W_FASERvmu_1p2_3, XSecCon_2mu_coh_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_coh_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

# Upper and Lower limits #
upper_N_1tau_coh_W_FASERvmu_1p2_3, upper_XSecCon_1tau_coh_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_coh_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')
upper_N_2tau_coh_W_FASERvmu_1p2_3, upper_XSecCon_2tau_coh_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_coh_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')
upper_N_2mu_coh_W_FASERvmu_1p2_3, upper_XSecCon_2mu_coh_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_coh_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')

lower_N_1tau_coh_W_FASERvmu_1p2_3, lower_XSecCon_1tau_coh_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_coh_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')
lower_N_2tau_coh_W_FASERvmu_1p2_3, lower_XSecCon_2tau_coh_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_coh_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')
lower_N_2mu_coh_W_FASERvmu_1p2_3, lower_XSecCon_2mu_coh_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_coh_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')

# vmubar Flux #
N_1tau_coh_W_FASERvmubar_1p2_3, XSecCon_1tau_coh_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_1tau_coh_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')
N_2tau_coh_W_FASERvmubar_1p2_3, XSecCon_2tau_coh_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2tau_coh_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')
N_2mu_coh_W_FASERvmubar_1p2_3, XSecCon_2mu_coh_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2mu_coh_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

# Upper and Lower limits #
upper_N_1tau_coh_W_FASERvmubar_1p2_3, upper_XSecCon_1tau_coh_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_1tau_coh_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')
upper_N_2tau_coh_W_FASERvmubar_1p2_3, upper_XSecCon_2tau_coh_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2tau_coh_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')
upper_N_2mu_coh_W_FASERvmubar_1p2_3, upper_XSecCon_2mu_coh_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2mu_coh_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')

lower_N_1tau_coh_W_FASERvmubar_1p2_3, lower_XSecCon_1tau_coh_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_1tau_coh_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')
lower_N_2tau_coh_W_FASERvmubar_1p2_3, lower_XSecCon_2tau_coh_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2tau_coh_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')
lower_N_2mu_coh_W_FASERvmubar_1p2_3, lower_XSecCon_2mu_coh_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2mu_coh_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')

## Incoherent ; proton ; Tungsten ##
# vmu Flux #
N_1tau_p_FASERvmu_1p2_3, XSecCon_1tau_p_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_p_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')
N_2tau_p_FASERvmu_1p2_3, XSecCon_2tau_p_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_p_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')
N_2mu_p_FASERvmu_1p2_3, XSecCon_2mu_p_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_p_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

# Upper and Lower limits #
upper_N_1tau_p_W_FASERvmu_1p2_3, upper_XSecCon_1tau_p_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_p_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')
upper_N_2tau_p_W_FASERvmu_1p2_3, upper_XSecCon_2tau_p_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_p_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')
upper_N_2mu_p_W_FASERvmu_1p2_3, upper_XSecCon_2mu_p_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_p_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')

lower_N_1tau_p_W_FASERvmu_1p2_3, lower_XSecCon_1tau_p_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_p_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')
lower_N_2tau_p_W_FASERvmu_1p2_3, lower_XSecCon_2tau_p_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_p_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')
lower_N_2mu_p_W_FASERvmu_1p2_3, lower_XSecCon_2mu_p_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_p_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')

# vmubar Flux #
N_1tau_p_FASERvmubar_1p2_3, XSecCon_1tau_p_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_1tau_p_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')
N_2tau_p_FASERvmubar_1p2_3, XSecCon_2tau_p_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2tau_p_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')
N_2mu_p_FASERvmubar_1p2_3, XSecCon_2mu_p_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2mu_p_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

# Upper and Lower limits #
upper_N_1tau_p_W_FASERvmubar_1p2_3, upper_XSecCon_1tau_p_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_1tau_p_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')
upper_N_2tau_p_W_FASERvmubar_1p2_3, upper_XSecCon_2tau_p_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2tau_p_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')
upper_N_2mu_p_W_FASERvmubar_1p2_3, upper_XSecCon_2mu_p_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2mu_p_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')

lower_N_1tau_p_W_FASERvmubar_1p2_3, lower_XSecCon_1tau_p_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_1tau_p_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')
lower_N_2tau_p_W_FASERvmubar_1p2_3, lower_XSecCon_2tau_p_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2tau_p_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')
lower_N_2mu_p_W_FASERvmubar_1p2_3, lower_XSecCon_2mu_p_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2mu_p_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')

## Incoherent ; neutron ; Tungsten ##
# vmu Flux #
N_1tau_n_FASERvmu_1p2_3, XSecCon_1tau_n_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_n_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')
N_2tau_n_FASERvmu_1p2_3, XSecCon_2tau_n_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_n_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')
N_2mu_n_FASERvmu_1p2_3, XSecCon_2mu_n_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_n_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

# Upper and Lower limits #
upper_N_1tau_n_W_FASERvmu_1p2_3, upper_XSecCon_1tau_n_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_n_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')
upper_N_2tau_n_W_FASERvmu_1p2_3, upper_XSecCon_2tau_n_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_n_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')
upper_N_2mu_n_W_FASERvmu_1p2_3, upper_XSecCon_2mu_n_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_n_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')

lower_N_1tau_n_W_FASERvmu_1p2_3, lower_XSecCon_1tau_n_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_n_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')
lower_N_2tau_n_W_FASERvmu_1p2_3, lower_XSecCon_2tau_n_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_n_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')
lower_N_2mu_n_W_FASERvmu_1p2_3, lower_XSecCon_2mu_n_W_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_n_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')

# vmubar Flux #
N_1tau_n_FASERvmubar_1p2_3, XSecCon_1tau_n_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_1tau_n_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')
N_2tau_n_FASERvmubar_1p2_3, XSecCon_2tau_n_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2tau_n_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')
N_2mu_n_FASERvmubar_1p2_3, XSecCon_2mu_n_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2mu_n_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

# Upper and Lower limits #
upper_N_1tau_n_W_FASERvmubar_1p2_3, upper_XSecCon_1tau_n_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_1tau_n_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')
upper_N_2tau_n_W_FASERvmubar_1p2_3, upper_XSecCon_2tau_n_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2tau_n_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')
upper_N_2mu_n_W_FASERvmubar_1p2_3, upper_XSecCon_2mu_n_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2mu_n_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')

lower_N_1tau_n_W_FASERvmubar_1p2_3, lower_XSecCon_1tau_n_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_1tau_n_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')
lower_N_2tau_n_W_FASERvmubar_1p2_3, lower_XSecCon_2tau_n_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2tau_n_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')
lower_N_2mu_n_W_FASERvmubar_1p2_3, lower_XSecCon_2mu_n_W_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_2mu_n_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')

## Incoherent ; proton + neutron ; Tungsten ##
# vmu Flux #
N_1tau_incoh_FASERvmu_1p2_3, XSecCon_1tau_incoh_FASERvmu_1p2_3 = (N_1tau_p_FASERvmu_1p2_3 + N_1tau_n_FASERvmu_1p2_3), (XSecCon_1tau_p_FASERvmu_1p2_3 + XSecCon_1tau_n_FASERvmu_1p2_3)
N_2tau_incoh_FASERvmu_1p2_3, XSecCon_2tau_incoh_FASERvmu_1p2_3 = (N_2tau_p_FASERvmu_1p2_3 + N_2tau_n_FASERvmu_1p2_3), (XSecCon_2tau_p_FASERvmu_1p2_3 + XSecCon_2tau_n_FASERvmu_1p2_3)
N_2mu_incoh_FASERvmu_1p2_3, XSecCon_2mu_incoh_FASERvmu_1p2_3 = (N_2mu_p_FASERvmu_1p2_3 + N_2mu_n_FASERvmu_1p2_3), (XSecCon_2mu_p_FASERvmu_1p2_3 + XSecCon_2mu_n_FASERvmu_1p2_3)

# Upper and Lower limits #
upper_N_1tau_incoh_FASERvmu_1p2_3, upper_XSecCon_1tau_incoh_FASERvmu_1p2_3 = (upper_N_1tau_p_W_FASERvmu_1p2_3 + upper_N_1tau_n_W_FASERvmu_1p2_3), (upper_XSecCon_1tau_p_W_FASERvmu_1p2_3 + upper_XSecCon_1tau_n_W_FASERvmu_1p2_3)
upper_N_2tau_incoh_FASERvmu_1p2_3, upper_XSecCon_2tau_incoh_FASERvmu_1p2_3 = (upper_N_2tau_p_W_FASERvmu_1p2_3 + upper_N_2tau_n_W_FASERvmu_1p2_3), (upper_XSecCon_2tau_p_W_FASERvmu_1p2_3 + upper_XSecCon_2tau_n_W_FASERvmu_1p2_3)
upper_N_2mu_incoh_FASERvmu_1p2_3, upper_XSecCon_2mu_incoh_FASERvmu_1p2_3 = (upper_N_2mu_p_W_FASERvmu_1p2_3 + upper_N_2mu_n_W_FASERvmu_1p2_3), (upper_XSecCon_2mu_p_W_FASERvmu_1p2_3 + upper_XSecCon_2mu_n_W_FASERvmu_1p2_3)

lower_N_1tau_incoh_FASERvmu_1p2_3, lower_XSecCon_1tau_incoh_FASERvmu_1p2_3 = (lower_N_1tau_p_W_FASERvmu_1p2_3 + lower_N_1tau_n_W_FASERvmu_1p2_3), (lower_XSecCon_1tau_p_W_FASERvmu_1p2_3 + lower_XSecCon_1tau_n_W_FASERvmu_1p2_3)
lower_N_2tau_incoh_FASERvmu_1p2_3, lower_XSecCon_2tau_incoh_FASERvmu_1p2_3 = (lower_N_2tau_p_W_FASERvmu_1p2_3 + lower_N_2tau_n_W_FASERvmu_1p2_3), (lower_XSecCon_2tau_p_W_FASERvmu_1p2_3 + lower_XSecCon_2tau_n_W_FASERvmu_1p2_3)
lower_N_2mu_incoh_FASERvmu_1p2_3, lower_XSecCon_2mu_incoh_FASERvmu_1p2_3 = (lower_N_2mu_p_W_FASERvmu_1p2_3 + lower_N_2mu_n_W_FASERvmu_1p2_3), (lower_XSecCon_2mu_p_W_FASERvmu_1p2_3 + lower_XSecCon_2mu_n_W_FASERvmu_1p2_3)

# vmubar Flux #
N_1tau_incoh_FASERvmubar_1p2_3, XSecCon_1tau_incoh_FASERvmubar_1p2_3 = (N_1tau_p_FASERvmubar_1p2_3 + N_1tau_n_FASERvmubar_1p2_3), (XSecCon_1tau_p_FASERvmubar_1p2_3 + XSecCon_1tau_n_FASERvmubar_1p2_3)
N_2tau_incoh_FASERvmubar_1p2_3, XSecCon_2tau_incoh_FASERvmubar_1p2_3 = (N_2tau_p_FASERvmubar_1p2_3 + N_2tau_n_FASERvmubar_1p2_3), (XSecCon_2tau_p_FASERvmubar_1p2_3 + XSecCon_2tau_n_FASERvmubar_1p2_3)
N_2mu_incoh_FASERvmubar_1p2_3, XSecCon_2mu_incoh_FASERvmubar_1p2_3 = (N_2mu_p_FASERvmubar_1p2_3 + N_2mu_n_FASERvmubar_1p2_3), (XSecCon_2mu_p_FASERvmubar_1p2_3 + XSecCon_2mu_n_FASERvmubar_1p2_3)

# Upper and Lower limits #
upper_N_1tau_incoh_FASERvmubar_1p2_3, upper_XSecCon_1tau_incoh_FASERvmubar_1p2_3 = (upper_N_1tau_p_W_FASERvmubar_1p2_3 + upper_N_1tau_n_W_FASERvmubar_1p2_3), (upper_XSecCon_1tau_p_W_FASERvmubar_1p2_3 + upper_XSecCon_1tau_n_W_FASERvmubar_1p2_3)
upper_N_2tau_incoh_FASERvmubar_1p2_3, upper_XSecCon_2tau_incoh_FASERvmubar_1p2_3 = (upper_N_2tau_p_W_FASERvmubar_1p2_3 + upper_N_2tau_n_W_FASERvmubar_1p2_3), (upper_XSecCon_2tau_p_W_FASERvmubar_1p2_3 + upper_XSecCon_2tau_n_W_FASERvmubar_1p2_3)
upper_N_2mu_incoh_FASERvmubar_1p2_3, upper_XSecCon_2mu_incoh_FASERvmubar_1p2_3 = (upper_N_2mu_p_W_FASERvmubar_1p2_3 + upper_N_2mu_n_W_FASERvmubar_1p2_3), (upper_XSecCon_2mu_p_W_FASERvmubar_1p2_3 + upper_XSecCon_2mu_n_W_FASERvmubar_1p2_3)

lower_N_1tau_incoh_FASERvmubar_1p2_3, lower_XSecCon_1tau_incoh_FASERvmubar_1p2_3 = (lower_N_1tau_p_W_FASERvmubar_1p2_3 + lower_N_1tau_n_W_FASERvmubar_1p2_3), (lower_XSecCon_1tau_p_W_FASERvmubar_1p2_3 + lower_XSecCon_1tau_n_W_FASERvmubar_1p2_3)
lower_N_2tau_incoh_FASERvmubar_1p2_3, lower_XSecCon_2tau_incoh_FASERvmubar_1p2_3 = (lower_N_2tau_p_W_FASERvmubar_1p2_3 + lower_N_2tau_n_W_FASERvmubar_1p2_3), (lower_XSecCon_2tau_p_W_FASERvmubar_1p2_3 + lower_XSecCon_2tau_n_W_FASERvmubar_1p2_3)
lower_N_2mu_incoh_FASERvmubar_1p2_3, lower_XSecCon_2mu_incoh_FASERvmubar_1p2_3 = (lower_N_2mu_p_W_FASERvmubar_1p2_3 + lower_N_2mu_n_W_FASERvmubar_1p2_3), (lower_XSecCon_2mu_p_W_FASERvmubar_1p2_3 + lower_XSecCon_2mu_n_W_FASERvmubar_1p2_3)

## vmuCC cross check ##
# vmu Flux #
N_vmuCC_FASERvmu_1p2_3, XSecCon_vmuCC_FASERvmu_1p2_3 = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_vmuCC_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

# vmubar Flux #
N_vmuCC_FASERvmubar_1p2_3, XSecCon_vmuCC_FASERvmubar_1p2_3 = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_vmuCC_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

N_DIS_vmuCC_FASERvmu, _ = CalculateEventsFASERv(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_DIS_vmuCC_matched)
N_DIS_vmubarCC_FASERvmubar, _ = CalculateEventsFASERv(flux_FASERvmubar, energy_FASERvmu, FASERvmu_xsec_DIS_vmubarCC_matched)

### FASERv2 vmu + vmubar Flux Events ###
## Coherent ; Tungsten ##
N_1tau_coh_W_FASER2vmu, XSecCon_1tau_coh_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_1tau_coh_W_matched, detector='FASER2')
N_2tau_coh_W_FASER2vmu, XSecCon_2tau_coh_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2tau_coh_W_matched, detector='FASER2')
N_2mu_coh_W_FASER2vmu, XSecCon_2mu_coh_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2mu_coh_W_matched, detector='FASER2')

# Upper and Lower limits #
upper_N_1tau_coh_W_FASER2vmu, upper_XSecCon_1tau_coh_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_1tau_coh_W_matched_upper, detector='FASER2')
upper_N_2tau_coh_W_FASER2vmu, upper_XSecCon_2tau_coh_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2tau_coh_W_matched_upper, detector='FASER2')
upper_N_2mu_coh_W_FASER2vmu, upper_XSecCon_2mu_coh_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2mu_coh_W_matched_upper, detector='FASER2')

lower_N_1tau_coh_W_FASER2vmu, lower_XSecCon_1tau_coh_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_1tau_coh_W_matched_lower, detector='FASER2')
lower_N_2tau_coh_W_FASER2vmu, lower_XSecCon_2tau_coh_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2tau_coh_W_matched_lower, detector='FASER2')
lower_N_2mu_coh_W_FASER2vmu, lower_XSecCon_2mu_coh_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2mu_coh_W_matched_lower, detector='FASER2')

## Incoherent ; proton ; Tungsten ##
# vmu Flux #
N_1tau_p_FASER2vmu, XSecCon_1tau_p_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_1tau_p_W_matched, detector='FASER2')
N_2tau_p_FASER2vmu, XSecCon_2tau_p_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2tau_p_W_matched, detector='FASER2')
N_2mu_p_FASER2vmu, XSecCon_2mu_p_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2mu_p_W_matched, detector='FASER2')

# Upper and Lower limits #
upper_N_1tau_p_W_FASER2vmu, upper_XSecCon_1tau_p_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_1tau_p_W_matched_upper, detector='FASER2')
upper_N_2tau_p_W_FASER2vmu, upper_XSecCon_2tau_p_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2tau_p_W_matched_upper, detector='FASER2')
upper_N_2mu_p_W_FASER2vmu, upper_XSecCon_2mu_p_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2mu_p_W_matched_upper, detector='FASER2')

lower_N_1tau_p_W_FASER2vmu, lower_XSecCon_1tau_p_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_1tau_p_W_matched_lower, detector='FASER2')
lower_N_2tau_p_W_FASER2vmu, lower_XSecCon_2tau_p_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2tau_p_W_matched_lower, detector='FASER2')
lower_N_2mu_p_W_FASER2vmu, lower_XSecCon_2mu_p_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2mu_p_W_matched_lower, detector='FASER2')

## Incoherent ; neutron ; Tungsten ##
# vmu Flux #
N_1tau_n_FASER2vmu, XSecCon_1tau_n_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_1tau_n_W_matched, detector='FASER2')
N_2tau_n_FASER2vmu, XSecCon_2tau_n_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2tau_n_W_matched, detector='FASER2')
N_2mu_n_FASER2vmu, XSecCon_2mu_n_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2mu_n_W_matched, detector='FASER2')

# Upper and Lower limits #
upper_N_1tau_n_W_FASER2vmu, upper_XSecCon_1tau_n_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_1tau_n_W_matched_upper, detector='FASER2')
upper_N_2tau_n_W_FASER2vmu, upper_XSecCon_2tau_n_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2tau_n_W_matched_upper, detector='FASER2')
upper_N_2mu_n_W_FASER2vmu, upper_XSecCon_2mu_n_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2mu_n_W_matched_upper, detector='FASER2')

lower_N_1tau_n_W_FASER2vmu, lower_XSecCon_1tau_n_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_1tau_n_W_matched_lower, detector='FASER2')
lower_N_2tau_n_W_FASER2vmu, lower_XSecCon_2tau_n_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2tau_n_W_matched_lower, detector='FASER2')
lower_N_2mu_n_W_FASER2vmu, lower_XSecCon_2mu_n_W_FASER2vmu = CalculateEvents(flux_FASERv2_vmu_vmubar, energy_FASERv2, FASERv2_xsec_2mu_n_W_matched_lower, detector='FASER2')

## Incoherent ; proton + neutron ; Tungsten ##
# vmu Flux #
N_1tau_incoh_FASER2vmu, XSecCon_1tau_incoh_FASER2vmu = (N_1tau_p_FASER2vmu + N_1tau_n_FASER2vmu), (XSecCon_1tau_p_FASER2vmu + XSecCon_1tau_n_FASER2vmu)
N_2tau_incoh_FASER2vmu, XSecCon_2tau_incoh_FASER2vmu = (N_2tau_p_FASER2vmu + N_2tau_n_FASER2vmu), (XSecCon_2tau_p_FASER2vmu + XSecCon_2tau_n_FASER2vmu)
N_2mu_incoh_FASER2vmu, XSecCon_2mu_incoh_FASER2vmu = (N_2mu_p_FASER2vmu + N_2mu_n_FASER2vmu), (XSecCon_2mu_p_FASER2vmu + XSecCon_2mu_n_FASER2vmu)

# Upper and Lower limits #
upper_N_1tau_incoh_FASER2vmu, upper_XSecCon_1tau_incoh_FASER2vmu = (upper_N_1tau_p_W_FASER2vmu + upper_N_1tau_n_W_FASER2vmu), (upper_XSecCon_1tau_p_W_FASER2vmu + upper_XSecCon_1tau_n_W_FASER2vmu)
upper_N_2tau_incoh_FASER2vmu, upper_XSecCon_2tau_incoh_FASER2vmu = (upper_N_2tau_p_W_FASER2vmu + upper_N_2tau_n_W_FASER2vmu), (upper_XSecCon_2tau_p_W_FASER2vmu + upper_XSecCon_2tau_n_W_FASER2vmu)
upper_N_2mu_incoh_FASER2vmu, upper_XSecCon_2mu_incoh_FASER2vmu = (upper_N_2mu_p_W_FASER2vmu + upper_N_2mu_n_W_FASER2vmu), (upper_XSecCon_2mu_p_W_FASER2vmu + upper_XSecCon_2mu_n_W_FASER2vmu)

lower_N_1tau_incoh_FASER2vmu, lower_XSecCon_1tau_incoh_FASER2vmu = (lower_N_1tau_p_W_FASER2vmu + lower_N_1tau_n_W_FASER2vmu), (lower_XSecCon_1tau_p_W_FASER2vmu + lower_XSecCon_1tau_n_W_FASER2vmu)
lower_N_2tau_incoh_FASER2vmu, lower_XSecCon_2tau_incoh_FASER2vmu = (lower_N_2tau_p_W_FASER2vmu + lower_N_2tau_n_W_FASER2vmu), (lower_XSecCon_2tau_p_W_FASER2vmu + lower_XSecCon_2tau_n_W_FASER2vmu)
lower_N_2mu_incoh_FASER2vmu, lower_XSecCon_2mu_incoh_FASER2vmu = (lower_N_2mu_p_W_FASER2vmu + lower_N_2mu_n_W_FASER2vmu), (lower_XSecCon_2mu_p_W_FASER2vmu + lower_XSecCon_2mu_n_W_FASER2vmu)

### MINOS Events ###
## Coherent ; Iron ##
# Neutrino Mode #
N_1tau_coh_Fe_MINOS_neutrino_vmu, XSecCon_1tau_coh_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_1tau_coh_Fe_matched, NTONNES=M_MINOS, detector='MINOS_neutrino')
N_2mu_coh_Fe_MINOS_neutrino_vmu, XSecCon_2mu_coh_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_2mu_coh_Fe_matched, NTONNES=M_MINOS, detector='MINOS_neutrino')

N_1tau_coh_Fe_MINOS_neutrino_vmubar, XSecCon_1tau_coh_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_coh_Fe_matched, NTONNES=M_MINOS, detector='MINOS_neutrino')
N_2mu_coh_Fe_MINOS_neutrino_vmubar, XSecCon_2mu_coh_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_coh_Fe_matched, NTONNES=M_MINOS, detector='MINOS_neutrino')

# Upper and Lower limits #
upper_N_1tau_coh_Fe_MINOS_neutrino_vmu, upper_XSecCon_1tau_coh_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_1tau_coh_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_neutrino')
upper_N_2mu_coh_Fe_MINOS_neutrino_vmu, upper_XSecCon_2mu_coh_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_2mu_coh_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_neutrino')

upper_N_1tau_coh_Fe_MINOS_neutrino_vmubar, upper_XSecCon_1tau_coh_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_coh_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_neutrino')
upper_N_2mu_coh_Fe_MINOS_neutrino_vmubar, upper_XSecCon_2mu_coh_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_coh_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_neutrino')

lower_N_1tau_coh_Fe_MINOS_neutrino_vmu, lower_XSecCon_1tau_coh_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_1tau_coh_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_neutrino')
lower_N_2mu_coh_Fe_MINOS_neutrino_vmu, lower_XSecCon_2mu_coh_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_2mu_coh_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_neutrino')

lower_N_1tau_coh_Fe_MINOS_neutrino_vmubar, lower_XSecCon_1tau_coh_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_coh_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_neutrino')
lower_N_2mu_coh_Fe_MINOS_neutrino_vmubar, lower_XSecCon_2mu_coh_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_coh_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_neutrino')

# Antineutrino Mode #
N_1tau_coh_Fe_MINOS_antineutrino_vmu, XSecCon_1tau_coh_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_1tau_coh_Fe_matched, NTONNES=M_MINOS, detector='MINOS_antineutrino')
N_2mu_coh_Fe_MINOS_antineutrino_vmu, XSecCon_2mu_coh_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_2mu_coh_Fe_matched, NTONNES=M_MINOS, detector='MINOS_antineutrino')

N_1tau_coh_Fe_MINOS_antineutrino_vmubar, XSecCon_1tau_coh_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_coh_Fe_matched, NTONNES=M_MINOS, detector='MINOS_antineutrino')
N_2mu_coh_Fe_MINOS_antineutrino_vmubar, XSecCon_2mu_coh_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_coh_Fe_matched, NTONNES=M_MINOS, detector='MINOS_antineutrino')

# Upper and Lower limits #
upper_N_1tau_coh_Fe_MINOS_antineutrino_vmu, upper_XSecCon_1tau_coh_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_1tau_coh_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_antineutrino')
upper_N_2mu_coh_Fe_MINOS_antineutrino_vmu, upper_XSecCon_2mu_coh_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_2mu_coh_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_antineutrino')

upper_N_1tau_coh_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_1tau_coh_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_coh_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_antineutrino')
upper_N_2mu_coh_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_2mu_coh_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_coh_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_antineutrino')

lower_N_1tau_coh_Fe_MINOS_antineutrino_vmu, lower_XSecCon_1tau_coh_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_1tau_coh_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_antineutrino')
lower_N_2mu_coh_Fe_MINOS_antineutrino_vmu, lower_XSecCon_2mu_coh_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_2mu_coh_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_antineutrino')

lower_N_1tau_coh_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_1tau_coh_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_coh_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_antineutrino')
lower_N_2mu_coh_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_2mu_coh_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_coh_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_antineutrino')

## Incoherent ; proton ; Iron ##
# Neutrino Mode #
N_1tau_p_Fe_MINOS_neutrino_vmu, XSecCon_1tau_p_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_1tau_p_Fe_matched, NTONNES=M_MINOS, detector='MINOS_neutrino')
N_2mu_p_Fe_MINOS_neutrino_vmu, XSecCon_2mu_p_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_2mu_p_Fe_matched, NTONNES=M_MINOS, detector='MINOS_neutrino')

N_1tau_p_Fe_MINOS_neutrino_vmubar, XSecCon_1tau_p_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_p_Fe_matched, NTONNES=M_MINOS, detector='MINOS_neutrino')
N_2mu_p_Fe_MINOS_neutrino_vmubar, XSecCon_2mu_p_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_p_Fe_matched, NTONNES=M_MINOS, detector='MINOS_neutrino')

# Upper and Lower limits #
upper_N_1tau_p_Fe_MINOS_neutrino_vmu, upper_XSecCon_1tau_p_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_1tau_p_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_neutrino')
upper_N_2mu_p_Fe_MINOS_neutrino_vmu, upper_XSecCon_2mu_p_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_2mu_p_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_neutrino')

upper_N_1tau_p_Fe_MINOS_neutrino_vmubar, upper_XSecCon_1tau_p_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_p_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_neutrino')
upper_N_2mu_p_Fe_MINOS_neutrino_vmubar, upper_XSecCon_2mu_p_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_p_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_neutrino')

lower_N_1tau_p_Fe_MINOS_neutrino_vmu, lower_XSecCon_1tau_p_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_1tau_p_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_neutrino')
lower_N_2mu_p_Fe_MINOS_neutrino_vmu, lower_XSecCon_2mu_p_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_2mu_p_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_neutrino')

lower_N_1tau_p_Fe_MINOS_neutrino_vmubar, lower_XSecCon_1tau_p_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_p_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_neutrino')
lower_N_2mu_p_Fe_MINOS_neutrino_vmubar, lower_XSecCon_2mu_p_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_p_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_neutrino')

# Antineutrino Mode #
N_1tau_p_Fe_MINOS_antineutrino_vmu, XSecCon_1tau_p_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_1tau_p_Fe_matched, NTONNES=M_MINOS, detector='MINOS_antineutrino')
N_2mu_p_Fe_MINOS_antineutrino_vmu, XSecCon_2mu_p_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_2mu_p_Fe_matched, NTONNES=M_MINOS, detector='MINOS_antineutrino')

N_1tau_p_Fe_MINOS_antineutrino_vmubar, XSecCon_1tau_p_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_p_Fe_matched, NTONNES=M_MINOS, detector='MINOS_antineutrino')
N_2mu_p_Fe_MINOS_antineutrino_vmubar, XSecCon_2mu_p_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_p_Fe_matched, NTONNES=M_MINOS, detector='MINOS_antineutrino')

# Upper and Lower limits #
upper_N_1tau_p_Fe_MINOS_antineutrino_vmu, upper_XSecCon_1tau_p_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_1tau_p_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_antineutrino')
upper_N_2mu_p_Fe_MINOS_antineutrino_vmu, upper_XSecCon_2mu_p_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_2mu_p_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_antineutrino')

upper_N_1tau_p_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_1tau_p_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_p_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_antineutrino')
upper_N_2mu_p_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_2mu_p_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_p_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_antineutrino')

lower_N_1tau_p_Fe_MINOS_antineutrino_vmu, lower_XSecCon_1tau_p_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_1tau_p_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_antineutrino')
lower_N_2mu_p_Fe_MINOS_antineutrino_vmu, lower_XSecCon_2mu_p_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_2mu_p_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_antineutrino')

lower_N_1tau_p_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_1tau_p_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_p_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_antineutrino')
lower_N_2mu_p_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_2mu_p_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_p_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_antineutrino')

## Incoherent ; neutron ; Iron ##
# Neutrino Mode #
N_1tau_n_Fe_MINOS_neutrino_vmu, XSecCon_1tau_n_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_1tau_n_Fe_matched, NTONNES=M_MINOS, detector='MINOS_neutrino')
N_2mu_n_Fe_MINOS_neutrino_vmu, XSecCon_2mu_n_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_2mu_n_Fe_matched, NTONNES=M_MINOS, detector='MINOS_neutrino')

N_1tau_n_Fe_MINOS_neutrino_vmubar, XSecCon_1tau_n_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_n_Fe_matched, NTONNES=M_MINOS, detector='MINOS_neutrino')
N_2mu_n_Fe_MINOS_neutrino_vmubar, XSecCon_2mu_n_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_n_Fe_matched, NTONNES=M_MINOS, detector='MINOS_neutrino')

# Upper and Lower limits #
upper_N_1tau_n_Fe_MINOS_neutrino_vmu, upper_XSecCon_1tau_n_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_1tau_n_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_neutrino')
upper_N_2mu_n_Fe_MINOS_neutrino_vmu, upper_XSecCon_2mu_n_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_2mu_n_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_neutrino')

upper_N_1tau_n_Fe_MINOS_neutrino_vmubar, upper_XSecCon_1tau_n_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_n_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_neutrino')
upper_N_2mu_n_Fe_MINOS_neutrino_vmubar, upper_XSecCon_2mu_n_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_n_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_neutrino')

lower_N_1tau_n_Fe_MINOS_neutrino_vmu, lower_XSecCon_1tau_n_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_1tau_n_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_neutrino')
lower_N_2mu_n_Fe_MINOS_neutrino_vmu, lower_XSecCon_2mu_n_Fe_MINOS_neutrino_vmu = CalculateEvents(flux_MINOS_neutrino_vmu, energy_MINOS, MINOS_xsec_2mu_n_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_neutrino')

lower_N_1tau_n_Fe_MINOS_neutrino_vmubar, lower_XSecCon_1tau_n_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_n_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_neutrino')
lower_N_2mu_n_Fe_MINOS_neutrino_vmubar, lower_XSecCon_2mu_n_Fe_MINOS_neutrino_vmubar = CalculateEvents(flux_MINOS_neutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_n_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_neutrino')

# Antineutrino Mode #
N_1tau_n_Fe_MINOS_antineutrino_vmu, XSecCon_1tau_n_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_1tau_n_Fe_matched, NTONNES=M_MINOS, detector='MINOS_antineutrino')
N_2mu_n_Fe_MINOS_antineutrino_vmu, XSecCon_2mu_n_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_2mu_n_Fe_matched, NTONNES=M_MINOS, detector='MINOS_antineutrino')

N_1tau_n_Fe_MINOS_antineutrino_vmubar, XSecCon_1tau_n_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_n_Fe_matched, NTONNES=M_MINOS, detector='MINOS_antineutrino')
N_2mu_n_Fe_MINOS_antineutrino_vmubar, XSecCon_2mu_n_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_n_Fe_matched, NTONNES=M_MINOS, detector='MINOS_antineutrino')

# Upper and Lower limits #
upper_N_1tau_n_Fe_MINOS_antineutrino_vmu, upper_XSecCon_1tau_n_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_1tau_n_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_antineutrino')
upper_N_2mu_n_Fe_MINOS_antineutrino_vmu, upper_XSecCon_2mu_n_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_2mu_n_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_antineutrino')

upper_N_1tau_n_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_1tau_n_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_n_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_antineutrino')
upper_N_2mu_n_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_2mu_n_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_n_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS_antineutrino')

lower_N_1tau_n_Fe_MINOS_antineutrino_vmu, lower_XSecCon_1tau_n_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_1tau_n_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_antineutrino')
lower_N_2mu_n_Fe_MINOS_antineutrino_vmu, lower_XSecCon_2mu_n_Fe_MINOS_antineutrino_vmu = CalculateEvents(flux_MINOS_antineutrino_vmu, energy_MINOS, MINOS_xsec_2mu_n_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_antineutrino')

lower_N_1tau_n_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_1tau_n_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_1tau_n_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_antineutrino')
lower_N_2mu_n_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_2mu_n_Fe_MINOS_antineutrino_vmubar = CalculateEvents(flux_MINOS_antineutrino_vmubar, energy_MINOS, MINOS_xsec_2mu_n_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS_antineutrino')

## Incoherent ; proton + neutron ; Iron
# Neutrino Mode #
N_1tau_incoh_Fe_MINOS_neutrino_vmu, XSecCon_1tau_incoh_Fe_MINOS_neutrino_vmu = (N_1tau_p_Fe_MINOS_neutrino_vmu + N_1tau_n_Fe_MINOS_neutrino_vmu), (XSecCon_1tau_p_Fe_MINOS_neutrino_vmu + XSecCon_1tau_n_Fe_MINOS_neutrino_vmu)
N_2mu_incoh_Fe_MINOS_neutrino_vmu, XSecCon_2mu_incoh_Fe_MINOS_neutrino_vmu = (N_2mu_p_Fe_MINOS_neutrino_vmu + N_2mu_n_Fe_MINOS_neutrino_vmu), (XSecCon_2mu_p_Fe_MINOS_neutrino_vmu + XSecCon_2mu_n_Fe_MINOS_neutrino_vmu)

N_1tau_incoh_Fe_MINOS_neutrino_vmubar, XSecCon_1tau_incoh_Fe_MINOS_neutrino_vmubar = (N_1tau_p_Fe_MINOS_neutrino_vmubar + N_1tau_n_Fe_MINOS_neutrino_vmubar), (XSecCon_1tau_p_Fe_MINOS_neutrino_vmubar + XSecCon_1tau_n_Fe_MINOS_neutrino_vmubar)
N_2mu_incoh_Fe_MINOS_neutrino_vmubar, XSecCon_2mu_incoh_Fe_MINOS_neutrino_vmubar = (N_2mu_p_Fe_MINOS_neutrino_vmubar + N_2mu_n_Fe_MINOS_neutrino_vmubar), (XSecCon_2mu_p_Fe_MINOS_neutrino_vmubar + XSecCon_2mu_n_Fe_MINOS_neutrino_vmubar)

# Upper and Lower limits #
upper_N_1tau_incoh_Fe_MINOS_neutrino_vmu, upper_XSecCon_1tau_incoh_Fe_MINOS_neutrino_vmu = (upper_N_1tau_p_Fe_MINOS_neutrino_vmu + upper_N_1tau_n_Fe_MINOS_neutrino_vmu), (upper_XSecCon_1tau_p_Fe_MINOS_neutrino_vmu + upper_XSecCon_1tau_n_Fe_MINOS_neutrino_vmu)
upper_N_2mu_incoh_Fe_MINOS_neutrino_vmu, upper_XSecCon_2mu_incoh_Fe_MINOS_neutrino_vmu = (upper_N_2mu_p_Fe_MINOS_neutrino_vmu + upper_N_2mu_n_Fe_MINOS_neutrino_vmu), (upper_XSecCon_2mu_p_Fe_MINOS_neutrino_vmu + upper_XSecCon_2mu_n_Fe_MINOS_neutrino_vmu)

upper_N_1tau_incoh_Fe_MINOS_neutrino_vmubar, upper_XSecCon_1tau_incoh_Fe_MINOS_neutrino_vmubar = (upper_N_1tau_p_Fe_MINOS_neutrino_vmubar + upper_N_1tau_n_Fe_MINOS_neutrino_vmubar), (upper_XSecCon_1tau_p_Fe_MINOS_neutrino_vmubar + upper_XSecCon_1tau_n_Fe_MINOS_neutrino_vmubar)
upper_N_2mu_incoh_Fe_MINOS_neutrino_vmubar, upper_XSecCon_2mu_incoh_Fe_MINOS_neutrino_vmubar = (upper_N_2mu_p_Fe_MINOS_neutrino_vmubar + upper_N_2mu_n_Fe_MINOS_neutrino_vmubar), (upper_XSecCon_2mu_p_Fe_MINOS_neutrino_vmubar + upper_XSecCon_2mu_n_Fe_MINOS_neutrino_vmubar)

lower_N_1tau_incoh_Fe_MINOS_neutrino_vmu, lower_XSecCon_1tau_incoh_Fe_MINOS_neutrino_vmu = (lower_N_1tau_p_Fe_MINOS_neutrino_vmu + lower_N_1tau_n_Fe_MINOS_neutrino_vmu), (lower_XSecCon_1tau_p_Fe_MINOS_neutrino_vmu + lower_XSecCon_1tau_n_Fe_MINOS_neutrino_vmu)
lower_N_2mu_incoh_Fe_MINOS_neutrino_vmu, lower_XSecCon_2mu_incoh_Fe_MINOS_neutrino_vmu = (lower_N_2mu_p_Fe_MINOS_neutrino_vmu + lower_N_2mu_n_Fe_MINOS_neutrino_vmu), (lower_XSecCon_2mu_p_Fe_MINOS_neutrino_vmu + lower_XSecCon_2mu_n_Fe_MINOS_neutrino_vmu)

lower_N_1tau_incoh_Fe_MINOS_neutrino_vmubar, lower_XSecCon_1tau_incoh_Fe_MINOS_neutrino_vmubar = (lower_N_1tau_p_Fe_MINOS_neutrino_vmubar + lower_N_1tau_n_Fe_MINOS_neutrino_vmubar), (lower_XSecCon_1tau_p_Fe_MINOS_neutrino_vmubar + lower_XSecCon_1tau_n_Fe_MINOS_neutrino_vmubar)
lower_N_2mu_incoh_Fe_MINOS_neutrino_vmubar, lower_XSecCon_2mu_incoh_Fe_MINOS_neutrino_vmubar = (lower_N_2mu_p_Fe_MINOS_neutrino_vmubar + lower_N_2mu_n_Fe_MINOS_neutrino_vmubar), (lower_XSecCon_2mu_p_Fe_MINOS_neutrino_vmubar + lower_XSecCon_2mu_n_Fe_MINOS_neutrino_vmubar)

# Antineutrino Mode #
N_1tau_incoh_Fe_MINOS_antineutrino_vmu, XSecCon_1tau_incoh_Fe_MINOS_antineutrino_vmu = (N_1tau_p_Fe_MINOS_antineutrino_vmu + N_1tau_n_Fe_MINOS_antineutrino_vmu), (XSecCon_1tau_p_Fe_MINOS_antineutrino_vmu + XSecCon_1tau_n_Fe_MINOS_antineutrino_vmu)
N_2mu_incoh_Fe_MINOS_antineutrino_vmu, XSecCon_2mu_incoh_Fe_MINOS_antineutrino_vmu = (N_2mu_p_Fe_MINOS_antineutrino_vmu + N_2mu_n_Fe_MINOS_antineutrino_vmu), (XSecCon_2mu_p_Fe_MINOS_antineutrino_vmu + XSecCon_2mu_n_Fe_MINOS_antineutrino_vmu)

N_1tau_incoh_Fe_MINOS_antineutrino_vmubar, XSecCon_1tau_incoh_Fe_MINOS_antineutrino_vmubar = (N_1tau_p_Fe_MINOS_antineutrino_vmubar + N_1tau_n_Fe_MINOS_antineutrino_vmubar), (XSecCon_1tau_p_Fe_MINOS_antineutrino_vmubar + XSecCon_1tau_n_Fe_MINOS_antineutrino_vmubar)
N_2mu_incoh_Fe_MINOS_antineutrino_vmubar, XSecCon_2mu_incoh_Fe_MINOS_antineutrino_vmubar = (N_2mu_p_Fe_MINOS_antineutrino_vmubar + N_2mu_n_Fe_MINOS_antineutrino_vmubar), (XSecCon_2mu_p_Fe_MINOS_antineutrino_vmubar + XSecCon_2mu_n_Fe_MINOS_antineutrino_vmubar)

# Upper and Lower limits #
upper_N_1tau_incoh_Fe_MINOS_antineutrino_vmu, upper_XSecCon_1tau_incoh_Fe_MINOS_antineutrino_vmu = (upper_N_1tau_p_Fe_MINOS_antineutrino_vmu + upper_N_1tau_n_Fe_MINOS_antineutrino_vmu), (upper_XSecCon_1tau_p_Fe_MINOS_antineutrino_vmu + upper_XSecCon_1tau_n_Fe_MINOS_antineutrino_vmu)
upper_N_2mu_incoh_Fe_MINOS_antineutrino_vmu, upper_XSecCon_2mu_incoh_Fe_MINOS_antineutrino_vmu = (upper_N_2mu_p_Fe_MINOS_antineutrino_vmu + upper_N_2mu_n_Fe_MINOS_antineutrino_vmu), (upper_XSecCon_2mu_p_Fe_MINOS_antineutrino_vmu + upper_XSecCon_2mu_n_Fe_MINOS_antineutrino_vmu)

upper_N_1tau_incoh_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_1tau_incoh_Fe_MINOS_antineutrino_vmubar = (upper_N_1tau_p_Fe_MINOS_antineutrino_vmubar + upper_N_1tau_n_Fe_MINOS_antineutrino_vmubar), (upper_XSecCon_1tau_p_Fe_MINOS_antineutrino_vmubar + upper_XSecCon_1tau_n_Fe_MINOS_antineutrino_vmubar)
upper_N_2mu_incoh_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_2mu_incoh_Fe_MINOS_antineutrino_vmubar = (upper_N_2mu_p_Fe_MINOS_antineutrino_vmubar + upper_N_2mu_n_Fe_MINOS_antineutrino_vmubar), (upper_XSecCon_2mu_p_Fe_MINOS_antineutrino_vmubar + upper_XSecCon_2mu_n_Fe_MINOS_antineutrino_vmubar)

lower_N_1tau_incoh_Fe_MINOS_antineutrino_vmu, lower_XSecCon_1tau_incoh_Fe_MINOS_antineutrino_vmu = (lower_N_1tau_p_Fe_MINOS_antineutrino_vmu + lower_N_1tau_n_Fe_MINOS_antineutrino_vmu), (lower_XSecCon_1tau_p_Fe_MINOS_antineutrino_vmu + lower_XSecCon_1tau_n_Fe_MINOS_antineutrino_vmu)
lower_N_2mu_incoh_Fe_MINOS_antineutrino_vmu, lower_XSecCon_2mu_incoh_Fe_MINOS_antineutrino_vmu = (lower_N_2mu_p_Fe_MINOS_antineutrino_vmu + lower_N_2mu_n_Fe_MINOS_antineutrino_vmu), (lower_XSecCon_2mu_p_Fe_MINOS_antineutrino_vmu + lower_XSecCon_2mu_n_Fe_MINOS_antineutrino_vmu)

lower_N_1tau_incoh_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_1tau_incoh_Fe_MINOS_antineutrino_vmubar = (lower_N_1tau_p_Fe_MINOS_antineutrino_vmubar + lower_N_1tau_n_Fe_MINOS_antineutrino_vmubar), (lower_XSecCon_1tau_p_Fe_MINOS_antineutrino_vmubar + lower_XSecCon_1tau_n_Fe_MINOS_antineutrino_vmubar)
lower_N_2mu_incoh_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_2mu_incoh_Fe_MINOS_antineutrino_vmubar = (lower_N_2mu_p_Fe_MINOS_antineutrino_vmubar + lower_N_2mu_n_Fe_MINOS_antineutrino_vmubar), (lower_XSecCon_2mu_p_Fe_MINOS_antineutrino_vmubar + lower_XSecCon_2mu_n_Fe_MINOS_antineutrino_vmubar)

### MINOS+ Events ###
## Coherent ; Iron ##
# Neutrino Mode #
N_1tau_coh_Fe_MINOSPlus_neutrino_vmu, XSecCon_1tau_coh_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_1tau_coh_Fe_matched, NTONNES=M_MINOS, detector='MINOS+')
N_2mu_coh_Fe_MINOSPlus_neutrino_vmu, XSecCon_2mu_coh_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_2mu_coh_Fe_matched, NTONNES=M_MINOS, detector='MINOS+')

N_1tau_coh_Fe_MINOSPlus_neutrino_vmubar, XSecCon_1tau_coh_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_1tau_coh_Fe_matched, NTONNES=M_MINOS, detector='MINOS+')
N_2mu_coh_Fe_MINOSPlus_neutrino_vmubar, XSecCon_2mu_coh_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_2mu_coh_Fe_matched, NTONNES=M_MINOS, detector='MINOS+')

# Upper and Lower limits #
upper_N_1tau_coh_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_1tau_coh_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_1tau_coh_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS+')
upper_N_2mu_coh_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_2mu_coh_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_2mu_coh_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS+')

upper_N_1tau_coh_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_1tau_coh_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_1tau_coh_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS+')
upper_N_2mu_coh_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_2mu_coh_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_2mu_coh_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS+')

lower_N_1tau_coh_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_1tau_coh_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_1tau_coh_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS+')
lower_N_2mu_coh_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_2mu_coh_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_2mu_coh_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS+')

lower_N_1tau_coh_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_1tau_coh_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_1tau_coh_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS+')
lower_N_2mu_coh_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_2mu_coh_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_2mu_coh_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS+')

## Incoherent ; proton ; Iron ##
# Neutrino Mode #
N_1tau_p_Fe_MINOSPlus_neutrino_vmu, XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_1tau_p_Fe_matched, NTONNES=M_MINOS, detector='MINOS+')
N_2mu_p_Fe_MINOSPlus_neutrino_vmu, XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_2mu_p_Fe_matched, NTONNES=M_MINOS, detector='MINOS+')

N_1tau_p_Fe_MINOSPlus_neutrino_vmubar, XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_1tau_p_Fe_matched, NTONNES=M_MINOS, detector='MINOS+')
N_2mu_p_Fe_MINOSPlus_neutrino_vmubar, XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_2mu_p_Fe_matched, NTONNES=M_MINOS, detector='MINOS+')

# Upper and Lower limits #
upper_N_1tau_p_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_1tau_p_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS+')
upper_N_2mu_p_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_2mu_p_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS+')

upper_N_1tau_p_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_1tau_p_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS+')
upper_N_2mu_p_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_2mu_p_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS+')

lower_N_1tau_p_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_1tau_p_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS+')
lower_N_2mu_p_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_2mu_p_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS+')

lower_N_1tau_p_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_1tau_p_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS+')
lower_N_2mu_p_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_2mu_p_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS+')

## Incoherent ; neutron ; Iron ##
# Neutrino Mode #
N_1tau_n_Fe_MINOSPlus_neutrino_vmu, XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_1tau_n_Fe_matched, NTONNES=M_MINOS, detector='MINOS+')
N_2mu_n_Fe_MINOSPlus_neutrino_vmu, XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_2mu_n_Fe_matched, NTONNES=M_MINOS, detector='MINOS+')

N_1tau_n_Fe_MINOSPlus_neutrino_vmubar, XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_1tau_n_Fe_matched, NTONNES=M_MINOS, detector='MINOS+')
N_2mu_n_Fe_MINOSPlus_neutrino_vmubar, XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_2mu_n_Fe_matched, NTONNES=M_MINOS, detector='MINOS+')

# Upper and Lower limits #
upper_N_1tau_n_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_1tau_n_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS+')
upper_N_2mu_n_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_2mu_n_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS+')

upper_N_1tau_n_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_1tau_n_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS+')
upper_N_2mu_n_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_2mu_n_Fe_matched_upper, NTONNES=M_MINOS, detector='MINOS+')

lower_N_1tau_n_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_1tau_n_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS+')
lower_N_2mu_n_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmu = CalculateEvents(flux_MINOSPlus_neutrino_vmu, energy_MINOSPlus, MINOSPlus_xsec_2mu_n_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS+')

lower_N_1tau_n_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_1tau_n_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS+')
lower_N_2mu_n_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmubar = CalculateEvents(flux_MINOSPlus_neutrino_vmubar, energy_MINOSPlus, MINOSPlus_xsec_2mu_n_Fe_matched_lower, NTONNES=M_MINOS, detector='MINOS+')

## Incoherent ; proton + neutron ; Iron
# Neutrino Mode #
N_1tau_incoh_Fe_MINOSPlus_neutrino_vmu, XSecCon_1tau_incoh_Fe_MINOSPlus_neutrino_vmu = (N_1tau_p_Fe_MINOSPlus_neutrino_vmu + N_1tau_n_Fe_MINOSPlus_neutrino_vmu), (XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmu + XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmu)
N_2mu_incoh_Fe_MINOSPlus_neutrino_vmu, XSecCon_2mu_incoh_Fe_MINOSPlus_neutrino_vmu = (N_2mu_p_Fe_MINOSPlus_neutrino_vmu + N_2mu_n_Fe_MINOSPlus_neutrino_vmu), (XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmu + XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmu)

N_1tau_incoh_Fe_MINOSPlus_neutrino_vmubar, XSecCon_1tau_incoh_Fe_MINOSPlus_neutrino_vmubar = (N_1tau_p_Fe_MINOSPlus_neutrino_vmubar + N_1tau_n_Fe_MINOSPlus_neutrino_vmubar), (XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmubar + XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmubar)
N_2mu_incoh_Fe_MINOSPlus_neutrino_vmubar, XSecCon_2mu_incoh_Fe_MINOSPlus_neutrino_vmubar = (N_2mu_p_Fe_MINOSPlus_neutrino_vmubar + N_2mu_n_Fe_MINOSPlus_neutrino_vmubar), (XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmubar + XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmubar)

# Upper and Lower limits #
upper_N_1tau_incoh_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_1tau_incoh_Fe_MINOSPlus_neutrino_vmu = (upper_N_1tau_p_Fe_MINOSPlus_neutrino_vmu + upper_N_1tau_n_Fe_MINOSPlus_neutrino_vmu), (upper_XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmu + upper_XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmu)
upper_N_2mu_incoh_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_2mu_incoh_Fe_MINOSPlus_neutrino_vmu = (upper_N_2mu_p_Fe_MINOSPlus_neutrino_vmu + upper_N_2mu_n_Fe_MINOSPlus_neutrino_vmu), (upper_XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmu + upper_XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmu)

upper_N_1tau_incoh_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_1tau_incoh_Fe_MINOSPlus_neutrino_vmubar = (upper_N_1tau_p_Fe_MINOSPlus_neutrino_vmubar + upper_N_1tau_n_Fe_MINOSPlus_neutrino_vmubar), (upper_XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmubar + upper_XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmubar)
upper_N_2mu_incoh_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_2mu_incoh_Fe_MINOSPlus_neutrino_vmubar = (upper_N_2mu_p_Fe_MINOSPlus_neutrino_vmubar + upper_N_2mu_n_Fe_MINOSPlus_neutrino_vmubar), (upper_XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmubar + upper_XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmubar)

lower_N_1tau_incoh_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_1tau_incoh_Fe_MINOSPlus_neutrino_vmu = (lower_N_1tau_p_Fe_MINOSPlus_neutrino_vmu + lower_N_1tau_n_Fe_MINOSPlus_neutrino_vmu), (lower_XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmu + lower_XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmu)
lower_N_2mu_incoh_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_2mu_incoh_Fe_MINOSPlus_neutrino_vmu = (lower_N_2mu_p_Fe_MINOSPlus_neutrino_vmu + lower_N_2mu_n_Fe_MINOSPlus_neutrino_vmu), (lower_XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmu + lower_XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmu)

lower_N_1tau_incoh_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_1tau_incoh_Fe_MINOSPlus_neutrino_vmubar = (lower_N_1tau_p_Fe_MINOSPlus_neutrino_vmubar + lower_N_1tau_n_Fe_MINOSPlus_neutrino_vmubar), (lower_XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmubar + lower_XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmubar)
lower_N_2mu_incoh_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_2mu_incoh_Fe_MINOSPlus_neutrino_vmubar = (lower_N_2mu_p_Fe_MINOSPlus_neutrino_vmubar + lower_N_2mu_n_Fe_MINOSPlus_neutrino_vmubar), (lower_XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmubar + lower_XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmubar)

### T2K - INGRID Events ###
## Coherent ; Iron ##
# Neutrino Mode #
N_1tau_coh_Fe_INGRID_neutrino_vmu, XSecCon_1tau_coh_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_coh_Fe_matched, NTONNES=M_INGRID, detector='INGRID')
N_2mu_coh_Fe_INGRID_neutrino_vmu, XSecCon_2mu_coh_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_coh_Fe_matched, NTONNES=M_INGRID, detector='INGRID')

N_1tau_coh_Fe_INGRID_neutrino_vmubar, XSecCon_1tau_coh_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_coh_Fe_matched, NTONNES=M_INGRID, detector='INGRID')
N_2mu_coh_Fe_INGRID_neutrino_vmubar, XSecCon_2mu_coh_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_coh_Fe_matched, NTONNES=M_INGRID, detector='INGRID')

# Upper and Lower limits #
upper_N_1tau_coh_Fe_INGRID_neutrino_vmu, upper_XSecCon_1tau_coh_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_coh_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRID')
upper_N_2mu_coh_Fe_INGRID_neutrino_vmu, upper_XSecCon_2mu_coh_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_coh_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRID')

upper_N_1tau_coh_Fe_INGRID_neutrino_vmubar, upper_XSecCon_1tau_coh_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_coh_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRID')
upper_N_2mu_coh_Fe_INGRID_neutrino_vmubar, upper_XSecCon_2mu_coh_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_coh_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRID')

lower_N_1tau_coh_Fe_INGRID_neutrino_vmu, lower_XSecCon_1tau_coh_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_coh_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRID')
lower_N_2mu_coh_Fe_INGRID_neutrino_vmu, lower_XSecCon_2mu_coh_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_coh_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRID')

lower_N_1tau_coh_Fe_INGRID_neutrino_vmubar, lower_XSecCon_1tau_coh_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_coh_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRID')
lower_N_2mu_coh_Fe_INGRID_neutrino_vmubar, lower_XSecCon_2mu_coh_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_coh_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRID')

## Incoherent ; proton ; Iron ##
# Neutrino Mode #
N_1tau_p_Fe_INGRID_neutrino_vmu, XSecCon_1tau_p_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_p_Fe_matched, NTONNES=M_INGRID, detector='INGRID')
N_2mu_p_Fe_INGRID_neutrino_vmu, XSecCon_2mu_p_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_p_Fe_matched, NTONNES=M_INGRID, detector='INGRID')

N_1tau_p_Fe_INGRID_neutrino_vmubar, XSecCon_1tau_p_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_p_Fe_matched, NTONNES=M_INGRID, detector='INGRID')
N_2mu_p_Fe_INGRID_neutrino_vmubar, XSecCon_2mu_p_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_p_Fe_matched, NTONNES=M_INGRID, detector='INGRID')

# Upper and Lower limits #
upper_N_1tau_p_Fe_INGRID_neutrino_vmu, upper_XSecCon_1tau_p_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_p_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRID')
upper_N_2mu_p_Fe_INGRID_neutrino_vmu, upper_XSecCon_2mu_p_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_p_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRID')

upper_N_1tau_p_Fe_INGRID_neutrino_vmubar, upper_XSecCon_1tau_p_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_p_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRID')
upper_N_2mu_p_Fe_INGRID_neutrino_vmubar, upper_XSecCon_2mu_p_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_p_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRID')

lower_N_1tau_p_Fe_INGRID_neutrino_vmu, lower_XSecCon_1tau_p_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_p_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRID')
lower_N_2mu_p_Fe_INGRID_neutrino_vmu, lower_XSecCon_2mu_p_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_p_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRID')

lower_N_1tau_p_Fe_INGRID_neutrino_vmubar, lower_XSecCon_1tau_p_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_p_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRID')
lower_N_2mu_p_Fe_INGRID_neutrino_vmubar, lower_XSecCon_2mu_p_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_p_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRID')

## Incoherent ; neutron ; Iron ##
# Neutrino Mode #
N_1tau_n_Fe_INGRID_neutrino_vmu, XSecCon_1tau_n_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_n_Fe_matched, NTONNES=M_INGRID, detector='INGRID')
N_2mu_n_Fe_INGRID_neutrino_vmu, XSecCon_2mu_n_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_n_Fe_matched, NTONNES=M_INGRID, detector='INGRID')

N_1tau_n_Fe_INGRID_neutrino_vmubar, XSecCon_1tau_n_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_n_Fe_matched, NTONNES=M_INGRID, detector='INGRID')
N_2mu_n_Fe_INGRID_neutrino_vmubar, XSecCon_2mu_n_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_n_Fe_matched, NTONNES=M_INGRID, detector='INGRID')

# Upper and Lower limits #
upper_N_1tau_n_Fe_INGRID_neutrino_vmu, upper_XSecCon_1tau_n_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_n_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRID')
upper_N_2mu_n_Fe_INGRID_neutrino_vmu, upper_XSecCon_2mu_n_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_n_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRID')

upper_N_1tau_n_Fe_INGRID_neutrino_vmubar, upper_XSecCon_1tau_n_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_n_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRID')
upper_N_2mu_n_Fe_INGRID_neutrino_vmubar, upper_XSecCon_2mu_n_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_n_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRID')

lower_N_1tau_n_Fe_INGRID_neutrino_vmu, lower_XSecCon_1tau_n_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_n_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRID')
lower_N_2mu_n_Fe_INGRID_neutrino_vmu, lower_XSecCon_2mu_n_Fe_INGRID_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_n_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRID')

lower_N_1tau_n_Fe_INGRID_neutrino_vmubar, lower_XSecCon_1tau_n_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_n_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRID')
lower_N_2mu_n_Fe_INGRID_neutrino_vmubar, lower_XSecCon_2mu_n_Fe_INGRID_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_n_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRID')

## Incoherent ; proton + neutron ; Iron
# Neutrino Mode #
N_1tau_incoh_Fe_INGRID_neutrino_vmu, XSecCon_1tau_incoh_Fe_INGRID_neutrino_vmu = (N_1tau_p_Fe_INGRID_neutrino_vmu + N_1tau_n_Fe_INGRID_neutrino_vmu), (XSecCon_1tau_p_Fe_INGRID_neutrino_vmu + XSecCon_1tau_n_Fe_INGRID_neutrino_vmu)
N_2mu_incoh_Fe_INGRID_neutrino_vmu, XSecCon_2mu_incoh_Fe_INGRID_neutrino_vmu = (N_2mu_p_Fe_INGRID_neutrino_vmu + N_2mu_n_Fe_INGRID_neutrino_vmu), (XSecCon_2mu_p_Fe_INGRID_neutrino_vmu + XSecCon_2mu_n_Fe_INGRID_neutrino_vmu)

N_1tau_incoh_Fe_INGRID_neutrino_vmubar, XSecCon_1tau_incoh_Fe_INGRID_neutrino_vmubar = (N_1tau_p_Fe_INGRID_neutrino_vmubar + N_1tau_n_Fe_INGRID_neutrino_vmubar), (XSecCon_1tau_p_Fe_INGRID_neutrino_vmubar + XSecCon_1tau_n_Fe_INGRID_neutrino_vmubar)
N_2mu_incoh_Fe_INGRID_neutrino_vmubar, XSecCon_2mu_incoh_Fe_INGRID_neutrino_vmubar = (N_2mu_p_Fe_INGRID_neutrino_vmubar + N_2mu_n_Fe_INGRID_neutrino_vmubar), (XSecCon_2mu_p_Fe_INGRID_neutrino_vmubar + XSecCon_2mu_n_Fe_INGRID_neutrino_vmubar)

# Upper and Lower limits #
upper_N_1tau_incoh_Fe_INGRID_neutrino_vmu, upper_XSecCon_1tau_incoh_Fe_INGRID_neutrino_vmu = (upper_N_1tau_p_Fe_INGRID_neutrino_vmu + upper_N_1tau_n_Fe_INGRID_neutrino_vmu), (upper_XSecCon_1tau_p_Fe_INGRID_neutrino_vmu + upper_XSecCon_1tau_n_Fe_INGRID_neutrino_vmu)
upper_N_2mu_incoh_Fe_INGRID_neutrino_vmu, upper_XSecCon_2mu_incoh_Fe_INGRID_neutrino_vmu = (upper_N_2mu_p_Fe_INGRID_neutrino_vmu + upper_N_2mu_n_Fe_INGRID_neutrino_vmu), (upper_XSecCon_2mu_p_Fe_INGRID_neutrino_vmu + upper_XSecCon_2mu_n_Fe_INGRID_neutrino_vmu)

upper_N_1tau_incoh_Fe_INGRID_neutrino_vmubar, upper_XSecCon_1tau_incoh_Fe_INGRID_neutrino_vmubar = (upper_N_1tau_p_Fe_INGRID_neutrino_vmubar + upper_N_1tau_n_Fe_INGRID_neutrino_vmubar), (upper_XSecCon_1tau_p_Fe_INGRID_neutrino_vmubar + upper_XSecCon_1tau_n_Fe_INGRID_neutrino_vmubar)
upper_N_2mu_incoh_Fe_INGRID_neutrino_vmubar, upper_XSecCon_2mu_incoh_Fe_INGRID_neutrino_vmubar = (upper_N_2mu_p_Fe_INGRID_neutrino_vmubar + upper_N_2mu_n_Fe_INGRID_neutrino_vmubar), (upper_XSecCon_2mu_p_Fe_INGRID_neutrino_vmubar + upper_XSecCon_2mu_n_Fe_INGRID_neutrino_vmubar)

lower_N_1tau_incoh_Fe_INGRID_neutrino_vmu, lower_XSecCon_1tau_incoh_Fe_INGRID_neutrino_vmu = (lower_N_1tau_p_Fe_INGRID_neutrino_vmu + lower_N_1tau_n_Fe_INGRID_neutrino_vmu), (lower_XSecCon_1tau_p_Fe_INGRID_neutrino_vmu + lower_XSecCon_1tau_n_Fe_INGRID_neutrino_vmu)
lower_N_2mu_incoh_Fe_INGRID_neutrino_vmu, lower_XSecCon_2mu_incoh_Fe_INGRID_neutrino_vmu = (lower_N_2mu_p_Fe_INGRID_neutrino_vmu + lower_N_2mu_n_Fe_INGRID_neutrino_vmu), (lower_XSecCon_2mu_p_Fe_INGRID_neutrino_vmu + lower_XSecCon_2mu_n_Fe_INGRID_neutrino_vmu)

lower_N_1tau_incoh_Fe_INGRID_neutrino_vmubar, lower_XSecCon_1tau_incoh_Fe_INGRID_neutrino_vmubar = (lower_N_1tau_p_Fe_INGRID_neutrino_vmubar + lower_N_1tau_n_Fe_INGRID_neutrino_vmubar), (lower_XSecCon_1tau_p_Fe_INGRID_neutrino_vmubar + lower_XSecCon_1tau_n_Fe_INGRID_neutrino_vmubar)
lower_N_2mu_incoh_Fe_INGRID_neutrino_vmubar, lower_XSecCon_2mu_incoh_Fe_INGRID_neutrino_vmubar = (lower_N_2mu_p_Fe_INGRID_neutrino_vmubar + lower_N_2mu_n_Fe_INGRID_neutrino_vmubar), (lower_XSecCon_2mu_p_Fe_INGRID_neutrino_vmubar + lower_XSecCon_2mu_n_Fe_INGRID_neutrino_vmubar)

### T2K - INGRID Phase 2 Events ###
## Coherent ; Iron ##
# Neutrino Mode #
N_1tau_coh_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_1tau_coh_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_coh_Fe_matched, NTONNES=M_INGRID, detector='INGRIDPhase2')
N_2mu_coh_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_2mu_coh_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_coh_Fe_matched, NTONNES=M_INGRID, detector='INGRIDPhase2')

N_1tau_coh_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_1tau_coh_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_coh_Fe_matched, NTONNES=M_INGRID, detector='INGRIDPhase2')
N_2mu_coh_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_2mu_coh_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_coh_Fe_matched, NTONNES=M_INGRID, detector='INGRIDPhase2')

# Upper and Lower limits #
upper_N_1tau_coh_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_1tau_coh_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_coh_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRIDPhase2')
upper_N_2mu_coh_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_2mu_coh_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_coh_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRIDPhase2')

upper_N_1tau_coh_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_1tau_coh_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_coh_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRIDPhase2')
upper_N_2mu_coh_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_2mu_coh_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_coh_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRIDPhase2')

lower_N_1tau_coh_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_1tau_coh_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_coh_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRIDPhase2')
lower_N_2mu_coh_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_2mu_coh_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_coh_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRIDPhase2')

lower_N_1tau_coh_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_1tau_coh_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_coh_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRIDPhase2')
lower_N_2mu_coh_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_2mu_coh_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_coh_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRIDPhase2')

## Incoherent ; proton ; Iron ##
# Neutrino Mode #
N_1tau_p_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_p_Fe_matched, NTONNES=M_INGRID, detector='INGRIDPhase2')
N_2mu_p_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_p_Fe_matched, NTONNES=M_INGRID, detector='INGRIDPhase2')

N_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_p_Fe_matched, NTONNES=M_INGRID, detector='INGRIDPhase2')
N_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_p_Fe_matched, NTONNES=M_INGRID, detector='INGRIDPhase2')

# Upper and Lower limits #
upper_N_1tau_p_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_p_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRIDPhase2')
upper_N_2mu_p_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_p_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRIDPhase2')

upper_N_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_p_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRIDPhase2')
upper_N_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_p_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRIDPhase2')

lower_N_1tau_p_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_p_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRIDPhase2')
lower_N_2mu_p_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_p_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRIDPhase2')

lower_N_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_p_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRIDPhase2')
lower_N_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_p_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRIDPhase2')

## Incoherent ; neutron ; Iron ##
# Neutrino Mode #
N_1tau_n_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_n_Fe_matched, NTONNES=M_INGRID, detector='INGRIDPhase2')
N_2mu_n_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_n_Fe_matched, NTONNES=M_INGRID, detector='INGRIDPhase2')

N_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_n_Fe_matched, NTONNES=M_INGRID, detector='INGRIDPhase2')
N_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_n_Fe_matched, NTONNES=M_INGRID, detector='INGRIDPhase2')

# Upper and Lower limits #
upper_N_1tau_n_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_n_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRIDPhase2')
upper_N_2mu_n_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_n_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRIDPhase2')

upper_N_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_n_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRIDPhase2')
upper_N_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_n_Fe_matched_upper, NTONNES=M_INGRID, detector='INGRIDPhase2')

lower_N_1tau_n_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_1tau_n_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRIDPhase2')
lower_N_2mu_n_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmu = CalculateEvents(flux_INGRID_neutrino_vmu, energy_INGRID, INGRID_xsec_2mu_n_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRIDPhase2')

lower_N_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_1tau_n_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRIDPhase2')
lower_N_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar = CalculateEvents(flux_INGRID_neutrino_vmubar, energy_INGRID, INGRID_xsec_2mu_n_Fe_matched_lower, NTONNES=M_INGRID, detector='INGRIDPhase2')

## Incoherent ; proton + neutron ; Iron
# Neutrino Mode #
N_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmu = (N_1tau_p_Fe_INGRIDPhase2_neutrino_vmu + N_1tau_n_Fe_INGRIDPhase2_neutrino_vmu), (XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmu + XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmu)
N_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmu = (N_2mu_p_Fe_INGRIDPhase2_neutrino_vmu + N_2mu_n_Fe_INGRIDPhase2_neutrino_vmu), (XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmu + XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmu)

N_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmubar = (N_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar + N_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar), (XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar + XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar)
N_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmubar = (N_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar + N_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar), (XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar + XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar)

# Upper and Lower limits #
upper_N_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmu = (upper_N_1tau_p_Fe_INGRIDPhase2_neutrino_vmu + upper_N_1tau_n_Fe_INGRIDPhase2_neutrino_vmu), (upper_XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmu + upper_XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmu)
upper_N_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmu = (upper_N_2mu_p_Fe_INGRIDPhase2_neutrino_vmu + upper_N_2mu_n_Fe_INGRIDPhase2_neutrino_vmu), (upper_XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmu + upper_XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmu)

upper_N_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmubar = (upper_N_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar + upper_N_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar), (upper_XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar + upper_XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar)
upper_N_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmubar = (upper_N_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar + upper_N_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar), (upper_XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar + upper_XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar)

lower_N_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmu = (lower_N_1tau_p_Fe_INGRIDPhase2_neutrino_vmu + lower_N_1tau_n_Fe_INGRIDPhase2_neutrino_vmu), (lower_XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmu + lower_XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmu)
lower_N_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmu = (lower_N_2mu_p_Fe_INGRIDPhase2_neutrino_vmu + lower_N_2mu_n_Fe_INGRIDPhase2_neutrino_vmu), (lower_XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmu + lower_XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmu)

lower_N_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmubar = (lower_N_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar + lower_N_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar), (lower_XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar + lower_XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar)
lower_N_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmubar = (lower_N_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar + lower_N_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar), (lower_XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar + lower_XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar)

### SHiP Events ###
## Coherent ; Tungsten ##
N_1tau_coh_W_SHiP_vmu, XSecCon_1tau_coh_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_1tau_coh_W_matched, NTONNES=M_SHiP, detector='SHiP')
N_2tau_coh_W_SHiP_vmu, XSecCon_2tau_coh_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2tau_coh_W_matched, NTONNES=M_SHiP, detector='SHiP')
N_2mu_coh_W_SHiP_vmu, XSecCon_2mu_coh_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2mu_coh_W_matched, NTONNES=M_SHiP, detector='SHiP')

N_1tau_coh_W_SHiP_vmubar, XSecCon_1tau_coh_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_1tau_coh_W_matched, NTONNES=M_SHiP, detector='SHiP')
N_2tau_coh_W_SHiP_vmubar, XSecCon_2tau_coh_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2tau_coh_W_matched, NTONNES=M_SHiP, detector='SHiP')
N_2mu_coh_W_SHiP_vmubar, XSecCon_2mu_coh_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2mu_coh_W_matched, NTONNES=M_SHiP, detector='SHiP')

# Upper and Lower limits #
upper_N_1tau_coh_W_SHiP_vmu, upper_XSecCon_1tau_coh_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_1tau_coh_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')
upper_N_2tau_coh_W_SHiP_vmu, upper_XSecCon_2tau_coh_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2tau_coh_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')
upper_N_2mu_coh_W_SHiP_vmu, upper_XSecCon_2mu_coh_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2mu_coh_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')

lower_N_1tau_coh_W_SHiP_vmu, lower_XSecCon_1tau_coh_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_1tau_coh_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')
lower_N_2tau_coh_W_SHiP_vmu, lower_XSecCon_2tau_coh_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2tau_coh_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')
lower_N_2mu_coh_W_SHiP_vmu, lower_XSecCon_2mu_coh_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2mu_coh_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')

upper_N_1tau_coh_W_SHiP_vmubar, upper_XSecCon_1tau_coh_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_1tau_coh_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')
upper_N_2tau_coh_W_SHiP_vmubar, upper_XSecCon_2tau_coh_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2tau_coh_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')
upper_N_2mu_coh_W_SHiP_vmubar, upper_XSecCon_2mu_coh_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2mu_coh_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')

lower_N_1tau_coh_W_SHiP_vmubar, lower_XSecCon_1tau_coh_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_1tau_coh_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')
lower_N_2tau_coh_W_SHiP_vmubar, lower_XSecCon_2tau_coh_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2tau_coh_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')
lower_N_2mu_coh_W_SHiP_vmubar, lower_XSecCon_2mu_coh_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2mu_coh_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')

## Incoherent ; proton ; Tungsten ##
N_1tau_p_W_SHiP_vmu, XSecCon_1tau_p_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_1tau_p_W_matched, NTONNES=M_SHiP, detector='SHiP')
N_2tau_p_W_SHiP_vmu, XSecCon_2tau_p_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2tau_p_W_matched, NTONNES=M_SHiP, detector='SHiP')
N_2mu_p_W_SHiP_vmu, XSecCon_2mu_p_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2mu_p_W_matched, NTONNES=M_SHiP, detector='SHiP')

N_1tau_p_W_SHiP_vmubar, XSecCon_1tau_p_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_1tau_p_W_matched, NTONNES=M_SHiP, detector='SHiP')
N_2tau_p_W_SHiP_vmubar, XSecCon_2tau_p_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2tau_p_W_matched, NTONNES=M_SHiP, detector='SHiP')
N_2mu_p_W_SHiP_vmubar, XSecCon_2mu_p_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2mu_p_W_matched, NTONNES=M_SHiP, detector='SHiP')

# Upper and Lower limits #
upper_N_1tau_p_W_SHiP_vmu, upper_XSecCon_1tau_p_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_1tau_p_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')
upper_N_2tau_p_W_SHiP_vmu, upper_XSecCon_2tau_p_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2tau_p_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')
upper_N_2mu_p_W_SHiP_vmu, upper_XSecCon_2mu_p_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2mu_p_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')

lower_N_1tau_p_W_SHiP_vmu, lower_XSecCon_1tau_p_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_1tau_p_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')
lower_N_2tau_p_W_SHiP_vmu, lower_XSecCon_2tau_p_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2tau_p_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')
lower_N_2mu_p_W_SHiP_vmu, lower_XSecCon_2mu_p_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2mu_p_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')

upper_N_1tau_p_W_SHiP_vmubar, upper_XSecCon_1tau_p_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_1tau_p_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')
upper_N_2tau_p_W_SHiP_vmubar, upper_XSecCon_2tau_p_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2tau_p_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')
upper_N_2mu_p_W_SHiP_vmubar, upper_XSecCon_2mu_p_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2mu_p_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')

lower_N_1tau_p_W_SHiP_vmubar, lower_XSecCon_1tau_p_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_1tau_p_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')
lower_N_2tau_p_W_SHiP_vmubar, lower_XSecCon_2tau_p_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2tau_p_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')
lower_N_2mu_p_W_SHiP_vmubar, lower_XSecCon_2mu_p_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2mu_p_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')

## Incoherent ; neutron ; Tungsten ##
N_1tau_n_W_SHiP_vmu, XSecCon_1tau_n_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_1tau_n_W_matched, NTONNES=M_SHiP, detector='SHiP')
N_2tau_n_W_SHiP_vmu, XSecCon_2tau_n_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2tau_n_W_matched, NTONNES=M_SHiP, detector='SHiP')
N_2mu_n_W_SHiP_vmu, XSecCon_2mu_n_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2mu_n_W_matched, NTONNES=M_SHiP, detector='SHiP')

N_1tau_n_W_SHiP_vmubar, XSecCon_1tau_n_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_1tau_n_W_matched, NTONNES=M_SHiP, detector='SHiP')
N_2tau_n_W_SHiP_vmubar, XSecCon_2tau_n_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2tau_n_W_matched, NTONNES=M_SHiP, detector='SHiP')
N_2mu_n_W_SHiP_vmubar, XSecCon_2mu_n_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2mu_n_W_matched, NTONNES=M_SHiP, detector='SHiP')

# Upper and Lower limits #
upper_N_1tau_n_W_SHiP_vmu, upper_XSecCon_1tau_n_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_1tau_n_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')
upper_N_2tau_n_W_SHiP_vmu, upper_XSecCon_2tau_n_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2tau_n_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')
upper_N_2mu_n_W_SHiP_vmu, upper_XSecCon_2mu_n_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2mu_n_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')

lower_N_1tau_n_W_SHiP_vmu, lower_XSecCon_1tau_n_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_1tau_n_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')
lower_N_2tau_n_W_SHiP_vmu, lower_XSecCon_2tau_n_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2tau_n_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')
lower_N_2mu_n_W_SHiP_vmu, lower_XSecCon_2mu_n_W_SHiP_vmu = CalculateEvents(flux_SHiP_vmu, energy_SHiP, SHiP_xsec_2mu_n_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')

upper_N_1tau_n_W_SHiP_vmubar, upper_XSecCon_1tau_n_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_1tau_n_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')
upper_N_2tau_n_W_SHiP_vmubar, upper_XSecCon_2tau_n_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2tau_n_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')
upper_N_2mu_n_W_SHiP_vmubar, upper_XSecCon_2mu_n_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2mu_n_W_matched_upper, NTONNES=M_SHiP, detector='SHiP')

lower_N_1tau_n_W_SHiP_vmubar, lower_XSecCon_1tau_n_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_1tau_n_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')
lower_N_2tau_n_W_SHiP_vmubar, lower_XSecCon_2tau_n_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2tau_n_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')
lower_N_2mu_n_W_SHiP_vmubar, lower_XSecCon_2mu_n_W_SHiP_vmubar = CalculateEvents(flux_SHiP_vmubar, energy_SHiP, SHiP_xsec_2mu_n_W_matched_lower, NTONNES=M_SHiP, detector='SHiP')

## Incoherent ; proton + neutron ; Tungsten ##
N_1tau_incoh_SHiP_vmu, XSecCon_1tau_incoh_SHiP_vmu = (N_1tau_p_W_SHiP_vmu + N_1tau_n_W_SHiP_vmu), (XSecCon_1tau_p_W_SHiP_vmu + XSecCon_1tau_n_W_SHiP_vmu)
N_2tau_incoh_SHiP_vmu, XSecCon_2tau_incoh_SHiP_vmu = (N_2tau_p_W_SHiP_vmu + N_2tau_n_W_SHiP_vmu), (XSecCon_2tau_p_W_SHiP_vmu + XSecCon_2tau_n_W_SHiP_vmu)
N_2mu_incoh_SHiP_vmu, XSecCon_2mu_incoh_SHiP_vmu = (N_2mu_p_W_SHiP_vmu + N_2mu_n_W_SHiP_vmu), (XSecCon_2mu_p_W_SHiP_vmu + XSecCon_2mu_n_W_SHiP_vmu)

N_1tau_incoh_SHiP_vmubar, XSecCon_1tau_incoh_SHiP_vmubar = (N_1tau_p_W_SHiP_vmubar + N_1tau_n_W_SHiP_vmubar), (XSecCon_1tau_p_W_SHiP_vmubar + XSecCon_1tau_n_W_SHiP_vmubar)
N_2tau_incoh_SHiP_vmubar, XSecCon_2tau_incoh_SHiP_vmubar = (N_2tau_p_W_SHiP_vmubar + N_2tau_n_W_SHiP_vmubar), (XSecCon_2tau_p_W_SHiP_vmubar + XSecCon_2tau_n_W_SHiP_vmubar)
N_2mu_incoh_SHiP_vmubar, XSecCon_2mu_incoh_SHiP_vmubar = (N_2mu_p_W_SHiP_vmubar + N_2mu_n_W_SHiP_vmubar), (XSecCon_2mu_p_W_SHiP_vmubar + XSecCon_2mu_n_W_SHiP_vmubar)

# Upper and Lower limits #
upper_N_1tau_incoh_SHiP_vmu, upper_XSecCon_1tau_incoh_SHiP_vmu = (upper_N_1tau_p_W_SHiP_vmu + upper_N_1tau_n_W_SHiP_vmu), (upper_XSecCon_1tau_p_W_SHiP_vmu + upper_XSecCon_1tau_n_W_SHiP_vmu)
upper_N_2tau_incoh_SHiP_vmu, upper_XSecCon_2tau_incoh_SHiP_vmu = (upper_N_2tau_p_W_SHiP_vmu + upper_N_2tau_n_W_SHiP_vmu), (upper_XSecCon_2tau_p_W_SHiP_vmu + upper_XSecCon_2tau_n_W_SHiP_vmu)
upper_N_2mu_incoh_SHiP_vmu, upper_XSecCon_2mu_incoh_SHiP_vmu = (upper_N_2mu_p_W_SHiP_vmu + upper_N_2mu_n_W_SHiP_vmu), (upper_XSecCon_2mu_p_W_SHiP_vmu + upper_XSecCon_2mu_n_W_SHiP_vmu)

lower_N_1tau_incoh_SHiP_vmu, lower_XSecCon_1tau_incoh_SHiP_vmu = (lower_N_1tau_p_W_SHiP_vmu + lower_N_1tau_n_W_SHiP_vmu), (lower_XSecCon_1tau_p_W_SHiP_vmu + lower_XSecCon_1tau_n_W_SHiP_vmu)
lower_N_2tau_incoh_SHiP_vmu, lower_XSecCon_2tau_incoh_SHiP_vmu = (lower_N_2tau_p_W_SHiP_vmu + lower_N_2tau_n_W_SHiP_vmu), (lower_XSecCon_2tau_p_W_SHiP_vmu + lower_XSecCon_2tau_n_W_SHiP_vmu)
lower_N_2mu_incoh_SHiP_vmu, lower_XSecCon_2mu_incoh_SHiP_vmu = (lower_N_2mu_p_W_SHiP_vmu + lower_N_2mu_n_W_SHiP_vmu), (lower_XSecCon_2mu_p_W_SHiP_vmu + lower_XSecCon_2mu_n_W_SHiP_vmu)

upper_N_1tau_incoh_SHiP_vmubar, upper_XSecCon_1tau_incoh_SHiP_vmubar = (upper_N_1tau_p_W_SHiP_vmubar + upper_N_1tau_n_W_SHiP_vmubar), (upper_XSecCon_1tau_p_W_SHiP_vmubar + upper_XSecCon_1tau_n_W_SHiP_vmubar)
upper_N_2tau_incoh_SHiP_vmubar, upper_XSecCon_2tau_incoh_SHiP_vmubar = (upper_N_2tau_p_W_SHiP_vmubar + upper_N_2tau_n_W_SHiP_vmubar), (upper_XSecCon_2tau_p_W_SHiP_vmubar + upper_XSecCon_2tau_n_W_SHiP_vmubar)
upper_N_2mu_incoh_SHiP_vmubar, upper_XSecCon_2mu_incoh_SHiP_vmubar = (upper_N_2mu_p_W_SHiP_vmubar + upper_N_2mu_n_W_SHiP_vmubar), (upper_XSecCon_2mu_p_W_SHiP_vmubar + upper_XSecCon_2mu_n_W_SHiP_vmubar)

lower_N_1tau_incoh_SHiP_vmubar, lower_XSecCon_1tau_incoh_SHiP_vmubar = (lower_N_1tau_p_W_SHiP_vmubar + lower_N_1tau_n_W_SHiP_vmubar), (lower_XSecCon_1tau_p_W_SHiP_vmubar + lower_XSecCon_1tau_n_W_SHiP_vmubar)
lower_N_2tau_incoh_SHiP_vmubar, lower_XSecCon_2tau_incoh_SHiP_vmubar = (lower_N_2tau_p_W_SHiP_vmubar + lower_N_2tau_n_W_SHiP_vmubar), (lower_XSecCon_2tau_p_W_SHiP_vmubar + lower_XSecCon_2tau_n_W_SHiP_vmubar)
lower_N_2mu_incoh_SHiP_vmubar, lower_XSecCon_2mu_incoh_SHiP_vmubar = (lower_N_2mu_p_W_SHiP_vmubar + lower_N_2mu_n_W_SHiP_vmubar), (lower_XSecCon_2mu_p_W_SHiP_vmubar + lower_XSecCon_2mu_n_W_SHiP_vmubar)

### SBND Events ###
## Coherent ; Argon ##
N_1e1mu_coh_Ar_SBND_vmu, XSecCon_1e1mu_coh_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_1e1mu_coh_Ar_matched, NTONNES=M_SBND, detector='SBND')
#N_1e1mu_coh_Ar_SBND_vmubar, XSecCon_1e1mu_coh_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_1e1mu_coh_Ar_matched, NTONNES=M_SBND, detector='SBND')

N_2mu_coh_Ar_SBND_vmu, XSecCon_2mu_coh_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_2mu_coh_Ar_matched, NTONNES=M_SBND, detector='SBND')
#N_2mu_coh_Ar_SBND_vmubar, XSecCon_2mu_coh_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_2mu_coh_Ar_matched, NTONNES=M_SBND, detector='SBND')

# Upper and Lower limits #
upper_N_1e1mu_coh_Ar_SBND_vmu, upper_XSecCon_1e1mu_coh_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_1e1mu_coh_Ar_matched_upper, NTONNES=M_SBND, detector='SBND')
lower_N_1e1mu_coh_Ar_SBND_vmu, lower_XSecCon_1e1mu_coh_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_1e1mu_coh_Ar_matched_lower, NTONNES=M_SBND, detector='SBND')

#upper_N_1e1mu_coh_Ar_SBND_vmubar, upper_XSecCon_1e1mu_coh_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_1e1mu_coh_Ar_matched_upper, NTONNES=M_SBND, detector='SBND')
#lower_N_1e1mu_coh_Ar_SBND_vmubar, lower_XSecCon_1e1mu_coh_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_1e1mu_coh_Ar_matched_lower, NTONNES=M_SBND, detector='SBND')

upper_N_2mu_coh_Ar_SBND_vmu, upper_XSecCon_2mu_coh_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_2mu_coh_Ar_matched_upper, NTONNES=M_SBND, detector='SBND')
lower_N_2mu_coh_Ar_SBND_vmu, lower_XSecCon_2mu_coh_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_2mu_coh_Ar_matched_lower, NTONNES=M_SBND, detector='SBND')

#upper_N_2mu_coh_Ar_SBND_vmubar, upper_XSecCon_2mu_coh_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_2mu_coh_Ar_matched_upper, NTONNES=M_SBND, detector='SBND')
#lower_N_2mu_coh_Ar_SBND_vmubar, lower_XSecCon_2mu_coh_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_2mu_coh_Ar_matched_lower, NTONNES=M_SBND, detector='SBND')

## Incoherent ; proton ; Argon ##
N_1e1mu_p_Ar_SBND_vmu, XSecCon_1e1mu_p_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_1e1mu_p_Ar_matched, NTONNES=M_SBND, detector='SBND')
#N_1e1mu_p_Ar_SBND_vmubar, XSecCon_1e1mu_p_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_1e1mu_p_Ar_matched, NTONNES=M_SBND, detector='SBND')

N_2mu_p_Ar_SBND_vmu, XSecCon_2mu_p_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_2mu_p_Ar_matched, NTONNES=M_SBND, detector='SBND')
#N_2mu_p_Ar_SBND_vmubar, XSecCon_2mu_p_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_2mu_p_Ar_matched, NTONNES=M_SBND, detector='SBND')

# Upper and Lower limits #
upper_N_1e1mu_p_Ar_SBND_vmu, upper_XSecCon_1e1mu_p_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_1e1mu_p_Ar_matched_upper, NTONNES=M_SBND, detector='SBND')
lower_N_1e1mu_p_Ar_SBND_vmu, lower_XSecCon_1e1mu_p_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_1e1mu_p_Ar_matched_lower, NTONNES=M_SBND, detector='SBND')

#upper_N_1e1mu_p_Ar_SBND_vmubar, upper_XSecCon_1e1mu_p_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_1e1mu_p_Ar_matched_upper, NTONNES=M_SBND, detector='SBND')
#lower_N_1e1mu_p_Ar_SBND_vmubar, lower_XSecCon_1e1mu_p_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_1e1mu_p_Ar_matched_lower, NTONNES=M_SBND, detector='SBND')

upper_N_2mu_p_Ar_SBND_vmu, upper_XSecCon_2mu_p_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_2mu_p_Ar_matched_upper, NTONNES=M_SBND, detector='SBND')
lower_N_2mu_p_Ar_SBND_vmu, lower_XSecCon_2mu_p_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_2mu_p_Ar_matched_lower, NTONNES=M_SBND, detector='SBND')

#upper_N_2mu_p_Ar_SBND_vmubar, upper_XSecCon_2mu_p_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_2mu_p_Ar_matched_upper, NTONNES=M_SBND, detector='SBND')
#lower_N_2mu_p_Ar_SBND_vmubar, lower_XSecCon_2mu_p_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_2mu_p_Ar_matched_lower, NTONNES=M_SBND, detector='SBND')

## Incoherent ; neutron ; Argon ##
N_1e1mu_n_Ar_SBND_vmu, XSecCon_1e1mu_n_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_1e1mu_n_Ar_matched, NTONNES=M_SBND, detector='SBND')
#N_1e1mu_n_Ar_SBND_vmubar, XSecCon_1e1mu_n_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_1e1mu_n_Ar_matched, NTONNES=M_SBND, detector='SBND')

N_2mu_n_Ar_SBND_vmu, XSecCon_2mu_n_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_2mu_n_Ar_matched, NTONNES=M_SBND, detector='SBND')
#N_2mu_n_Ar_SBND_vmubar, XSecCon_2mu_n_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_2mu_n_Ar_matched, NTONNES=M_SBND, detector='SBND')

# Upper and Lower limits #
upper_N_1e1mu_n_Ar_SBND_vmu, upper_XSecCon_1e1mu_n_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_1e1mu_n_Ar_matched_upper, NTONNES=M_SBND, detector='SBND')
lower_N_1e1mu_n_Ar_SBND_vmu, lower_XSecCon_1e1mu_n_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_1e1mu_n_Ar_matched_lower, NTONNES=M_SBND, detector='SBND')

#upper_N_1e1mu_n_Ar_SBND_vmubar, upper_XSecCon_1e1mu_n_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_1e1mu_n_Ar_matched_upper, NTONNES=M_SBND, detector='SBND')
#lower_N_1e1mu_n_Ar_SBND_vmubar, lower_XSecCon_1e1mu_n_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_1e1mu_n_Ar_matched_lower, NTONNES=M_SBND, detector='SBND')

upper_N_2mu_n_Ar_SBND_vmu, upper_XSecCon_2mu_n_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_2mu_n_Ar_matched_upper, NTONNES=M_SBND, detector='SBND')
lower_N_2mu_n_Ar_SBND_vmu, lower_XSecCon_2mu_n_Ar_SBND_vmu = CalculateEvents(flux_SBND_vmu, energy_SBND_vmu, SBND_xsec_2mu_n_Ar_matched_lower, NTONNES=M_SBND, detector='SBND')

#upper_N_2mu_n_Ar_SBND_vmubar, upper_XSecCon_2mu_n_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_2mu_n_Ar_matched_upper, NTONNES=M_SBND, detector='SBND')
#lower_N_2mu_n_Ar_SBND_vmubar, lower_XSecCon_2mu_n_Ar_SBND_vmubar = CalculateEvents(flux_SBND_vmubar, energy_SBND_vmubar, SBND_xsec_2mu_n_Ar_matched_lower, NTONNES=M_SBND, detector='SBND')

## Incoherent ; proton + neutron ; Argon ##
N_1e1mu_incoh_SBND_vmu, XSecCon_1e1mu_incoh_SBND_vmu = (N_1e1mu_p_Ar_SBND_vmu + N_1e1mu_n_Ar_SBND_vmu), (XSecCon_1e1mu_p_Ar_SBND_vmu + XSecCon_1e1mu_n_Ar_SBND_vmu)
#N_1e1mu_incoh_SBND_vmubar, XSecCon_1e1mu_incoh_SBND_vmubar = (N_1e1mu_p_Ar_SBND_vmubar + N_1e1mu_n_Ar_SBND_vmubar), (XSecCon_1e1mu_p_Ar_SBND_vmubar + XSecCon_1e1mu_n_Ar_SBND_vmubar)

N_2mu_incoh_SBND_vmu, XSecCon_2mu_incoh_SBND_vmu = (N_2mu_p_Ar_SBND_vmu + N_2mu_n_Ar_SBND_vmu), (XSecCon_2mu_p_Ar_SBND_vmu + XSecCon_2mu_n_Ar_SBND_vmu)
#N_2mu_incoh_SBND_vmubar, XSecCon_2mu_incoh_SBND_vmubar = (N_2mu_p_Ar_SBND_vmubar + N_2mu_n_Ar_SBND_vmubar), (XSecCon_2mu_p_Ar_SBND_vmubar + XSecCon_2mu_n_Ar_SBND_vmubar)

# Upper and Lower limits #
upper_N_2mu_incoh_SBND_vmu, upper_XSecCon_2mu_incoh_SBND_vmu = (upper_N_2mu_p_Ar_SBND_vmu + upper_N_2mu_n_Ar_SBND_vmu), (upper_XSecCon_2mu_p_Ar_SBND_vmu + upper_XSecCon_2mu_n_Ar_SBND_vmu)
lower_N_2mu_incoh_SBND_vmu, lower_XSecCon_2mu_incoh_SBND_vmu = (lower_N_2mu_p_Ar_SBND_vmu + lower_N_2mu_n_Ar_SBND_vmu), (lower_XSecCon_2mu_p_Ar_SBND_vmu + lower_XSecCon_2mu_n_Ar_SBND_vmu)

#upper_N_2mu_incoh_SBND_vmubar, upper_XSecCon_2mu_incoh_SBND_vmubar = (upper_N_2mu_p_Ar_SBND_vmubar + upper_N_2mu_n_Ar_SBND_vmubar), (upper_XSecCon_2mu_p_Ar_SBND_vmubar + upper_XSecCon_2mu_n_Ar_SBND_vmubar)
#lower_N_2mu_incoh_SBND_vmubar, lower_XSecCon_2mu_incoh_SBND_vmubar = (lower_N_2mu_p_Ar_SBND_vmubar + lower_N_2mu_n_Ar_SBND_vmubar), (lower_XSecCon_2mu_p_Ar_SBND_vmubar + lower_XSecCon_2mu_n_Ar_SBND_vmubar)

upper_N_1e1mu_incoh_SBND_vmu, upper_XSecCon_1e1mu_incoh_SBND_vmu = (upper_N_1e1mu_p_Ar_SBND_vmu + upper_N_1e1mu_n_Ar_SBND_vmu), (upper_XSecCon_1e1mu_p_Ar_SBND_vmu + upper_XSecCon_1e1mu_n_Ar_SBND_vmu)
lower_N_1e1mu_incoh_SBND_vmu, lower_XSecCon_1e1mu_incoh_SBND_vmu = (lower_N_1e1mu_p_Ar_SBND_vmu + lower_N_1e1mu_n_Ar_SBND_vmu), (lower_XSecCon_1e1mu_p_Ar_SBND_vmu + lower_XSecCon_1e1mu_n_Ar_SBND_vmu)

#upper_N_1e1mu_incoh_SBND_vmubar, upper_XSecCon_1e1mu_incoh_SBND_vmubar = (upper_N_1e1mu_p_Ar_SBND_vmubar + upper_N_1e1mu_n_Ar_SBND_vmubar), (upper_XSecCon_1e1mu_p_Ar_SBND_vmubar + upper_XSecCon_1e1mu_n_Ar_SBND_vmubar)
#lower_N_1e1mu_incoh_SBND_vmubar, lower_XSecCon_1e1mu_incoh_SBND_vmubar = (lower_N_1e1mu_p_Ar_SBND_vmubar + lower_N_1e1mu_n_Ar_SBND_vmubar), (lower_XSecCon_1e1mu_p_Ar_SBND_vmubar + lower_XSecCon_1e1mu_n_Ar_SBND_vmubar)


#########################
### Correction Factor ###
#########################

## DUNE ; Neutrino Mode ##
N_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_coh_Ar_matched)

N_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_p_Ar_matched)
N_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_n_Ar_matched)

N_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss = (N_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss + N_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss), (XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss + XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss)]

# Upper and lower limits #
upper_N_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_coh_Ar_matched)
lower_N_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_coh_Ar_matched)

upper_N_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_p_Ar_matched)
lower_N_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_p_Ar_matched)

upper_N_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_n_Ar_matched)
lower_N_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_neutrino_vmu, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_n_Ar_matched)

upper_N_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss = (upper_N_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss + upper_N_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss), (upper_XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss + upper_XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss)]
lower_N_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss = (lower_N_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss + lower_N_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss), (lower_XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss + lower_XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss)]

## DUNE ; Antineutrino Mode ##
N_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_coh_Ar_matched)

N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_p_Ar_matched)
N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_n_Ar_matched)

N_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss = (N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss + N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss), (XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss + XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss)]

# Upper and lower limits #
upper_N_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_coh_Ar_matched)
lower_N_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_coh_Ar_matched)

upper_N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_p_Ar_matched)
lower_N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_p_Ar_matched)

upper_N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_n_Ar_matched)
lower_N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_antineutrino_vmu, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3, correction_factor=DUNE_correction_factor_RMiss_1tau_n_Ar_matched)

upper_N_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss = (upper_N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss + upper_N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss), (upper_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss + upper_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss)]
lower_N_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss = (lower_N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss + lower_N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss), (lower_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss + lower_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss)]

## DUNE Tau Optimized ; Neutrino Mode ##
N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_coh_Ar_matched)

N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_p_Ar_matched)
N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_n_Ar_matched)

N_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss = (N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss + N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss), (XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss + XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss)]

# Upper and lower limits #
upper_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_coh_Ar_matched)
lower_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_coh_Ar_matched)

upper_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_p_Ar_matched)
lower_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_p_Ar_matched)

upper_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_n_Ar_matched)
lower_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_neutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_n_Ar_matched)

upper_N_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss = (upper_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss + upper_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss), (upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss + upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss)]
lower_N_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss = (lower_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss + lower_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss), (lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss + lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss)]

## DUNE Tau Optimized ; Antineutrino Mode ##
N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_coh_Ar_matched)

N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_p_Ar_matched)
N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_n_Ar_matched)

N_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss = (N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss + N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss), (XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss + XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss), [sum(x) for x in zip(Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss)]

# Upper and lower limits #
upper_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_coh_Ar_matched)
lower_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_coh_Ar_matched)

upper_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_p_Ar_matched)
lower_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_p_Ar_matched)

upper_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_n_Ar_matched)
lower_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss = CalculateEvents(flux_DUNE_tau_opt_antineutrino_vmu, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3, correction_factor=DUNE_tau_opt_correction_factor_RMiss_1tau_n_Ar_matched)

upper_N_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss = (upper_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss + upper_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss), (upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss + upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss), [sum(x) for x in zip(upper_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss)]
lower_N_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss = (lower_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss + lower_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss), (lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss + lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss), [sum(x) for x in zip(lower_Nbin_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_Nbin_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss)]

############################
###### Event Printout ######
############################

def PrintOutEvent(filename, fluxname, tonnes, years, events, xsec_con, upper_events=0, lower_events=0, upper_xsec_con=0, lower_xsec_con=0, detector='DUNE'):
    if detector == 'DUNE':
        print("\t{} tonne Ar, {} year, {} flux: {:.3e} + {:.3e} - {:.3e} ; XSecCon: {:.3e} + {:.3e} - {:.3e}".format(tonnes, years, fluxname, events, upper_events-events, events-lower_events, xsec_con, upper_xsec_con-xsec_con, xsec_con-lower_xsec_con), file=filename)
    if detector == 'FASER':
        print("\t{} tonne W, {} fb^-1 run L, {} flux: {:.3e} + {:.3e} - {:.3e} ; XSecCon: {:.3e} + {:.3e} - {:.3e}".format(M_FASERv, L_Run, fluxname, events, upper_events-events, events-lower_events, xsec_con, upper_xsec_con-xsec_con, xsec_con-lower_xsec_con), file=filename)
    if detector == 'FASER2':
        print("\t{} tonne W, {} fb^-1 run L, {} flux: {:.3e} + {:.3e} - {:.3e} ; XSecCon: {:.3e} + {:.3e} - {:.3e}".format(M_FASERv2, L_Run2, fluxname, events, upper_events-events, events-lower_events, xsec_con, upper_xsec_con-xsec_con, xsec_con-lower_xsec_con), file=filename)
    if detector == 'MINOS_neutrino':
        print("\t{} tonnes Fe, {}e20 POT, {} flux: {:.3e} + {:.3e} - {:.3e} ; XSecCon: {:.3e} + {:.3e} - {:.3e}".format(tonnes, N_POT_MINOS_neutrino, fluxname, events, upper_events-events, events-lower_events, xsec_con, upper_xsec_con-xsec_con, xsec_con-lower_xsec_con), file=filename)
    if detector == 'MINOS_antineutrino':
        print("\t{} tonnes Fe, {}e20 POT, {} flux: {:.3e} + {:.3e} - {:.3e} ; XSecCon: {:.3e} + {:.3e} - {:.3e}".format(tonnes, N_POT_MINOS_antineutrino, fluxname, events, upper_events-events, events-lower_events, xsec_con, upper_xsec_con-xsec_con, xsec_con-lower_xsec_con), file=filename)
    if detector == 'MINOS+':
        print("\t{} tonnes Fe, {}e20 POT, {} flux: {:.3e} + {:.3e} - {:.3e} ; XSecCon: {:.3e} + {:.3e} - {:.3e}".format(tonnes, N_POT_MINOSPlus_neutrino, fluxname, events, upper_events-events, events-lower_events, xsec_con, upper_xsec_con-xsec_con, xsec_con-lower_xsec_con), file=filename)
    if detector == 'INGRID':
        print("\t{} tonnes Fe, {}e20 POT, {} flux: {:.3e} + {:.3e} - {:.3e} ; XSecCon: {:.3e} + {:.3e} - {:.3e}".format(M_INGRID, N_POT_INGRID, fluxname, events, upper_events-events, events-lower_events, xsec_con, upper_xsec_con-xsec_con, xsec_con-lower_xsec_con), file=filename)
    if detector == 'INGRID2':
        print("\t{} tonnes Fe, {}e20 POT, {} flux: {:.3e} + {:.3e} - {:.3e} ; XSecCon: {:.3e} + {:.3e} - {:.3e}".format(M_INGRID, N_POT_INGRID_phase2, fluxname, events, upper_events-events, events-lower_events, xsec_con, upper_xsec_con-xsec_con, xsec_con-lower_xsec_con), file=filename)
    if detector == 'SHiP':
        print("\t{} tonne W, {} POT, {} flux: {:.3e} + {:.3e} - {:.3e} ; XSecCon: {:.3e} + {:.3e} - {:.3e}".format(M_SHiP, N_POT_SHiP, fluxname, events, upper_events-events, events-lower_events, xsec_con, upper_xsec_con-xsec_con, xsec_con-lower_xsec_con), file=filename)
    if detector == 'SBND':
        print("\t{} tonne Ar, {} POT, {} flux: {:.3e} + {:.3e} - {:.3e} ; XSecCon: {:.3e} + {:.3e} - {:.3e}".format(M_SBND, N_POT_SBND, fluxname, events, upper_events-events, events-lower_events, xsec_con, upper_xsec_con-xsec_con, xsec_con-lower_xsec_con), file=filename)



def WriteOutFile(filename):
    with open(filename,'w') as textfile:
        print("Events expected at DUNE ND or FASERv for:", file=textfile)
        print("The following integrated fluxes were used:", file=textfile)
        print("%%%%%%%%%% DUNE -- NEUTRINO MODE %%%%%%%%%%", file=textfile)
        print("-- vmu --", file=textfile)
        print("\tDUNE Standard Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_neutrino_vmu_integrated_flux), file=textfile)
        print("\tDUNE Tau-optimized Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_tau_opt_neutrino_vmu_integrated_flux), file=textfile)
        print("-- vmubar --", file=textfile)
        print("\tDUNE Standard Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_neutrino_vmubar_integrated_flux), file=textfile)
        print("\tDUNE Tau-optimized Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_tau_opt_neutrino_vmubar_integrated_flux), file=textfile)
        print("-- ve --", file=textfile)
        print("\tDUNE Standard Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_neutrino_ve_integrated_flux), file=textfile)
        print("\tDUNE Tau-optimized Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_tau_opt_neutrino_ve_integrated_flux), file=textfile)
        print("-- vebar --", file=textfile)
        print("\tDUNE Standard Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_neutrino_vebar_integrated_flux), file=textfile)
        print("\tDUNE Tau-optimized Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_tau_opt_neutrino_vebar_integrated_flux), file=textfile)
        print("%%%%%%%% DUNE -- ANTINEUTRINO MODE %%%%%%%%", file=textfile)
        print("-- vmu --", file=textfile)
        print("\tDUNE Standard Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_antineutrino_vmu_integrated_flux), file=textfile)
        print("\tDUNE Tau-optimized Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_tau_opt_antineutrino_vmu_integrated_flux), file=textfile)
        print("-- vmubar --", file=textfile)
        print("\tDUNE Standard Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_antineutrino_vmubar_integrated_flux), file=textfile)
        print("\tDUNE Tau-optimized Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_tau_opt_antineutrino_vmubar_integrated_flux), file=textfile)
        print("-- ve --", file=textfile)
        print("\tDUNE Standard Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_antineutrino_ve_integrated_flux), file=textfile)
        print("\tDUNE Tau-optimized Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_tau_opt_antineutrino_ve_integrated_flux), file=textfile)
        print("-- vebar --", file=textfile)
        print("\tDUNE Standard Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_antineutrino_vebar_integrated_flux), file=textfile)
        print("\tDUNE Tau-optimized Flux: {:.3e} [neutrinos / m^2 POT]".format(DUNE_tau_opt_antineutrino_vebar_integrated_flux), file=textfile)
        print("%%%%%%%%%%%%%%%%%% FASER %%%%%%%%%%%%%%%%%%", file=textfile)
        print("\tFASER vmu Flux: {:.3e} [neutrinos / m^2 fb^-1]".format(FASER_integrated_flux), file=textfile)
        print("\tFASER vmubar Flux: {:.3e} [neutrinos / m^2 fb^-1]".format(FASER_vmubar_integrated_flux), file=textfile)
        print("\tFASERv vmu + vmubar Flux: {:.3e} [neutrinos / m^2 fb^-1]".format(FASER_integrated_flux + FASER_vmubar_integrated_flux), file=textfile)
        print("\tFASERv2 vmu + vmubar Flux: {:.3e} [neutrinos / m^2 fb^-1]".format(FASERv2_integrated_flux), file=textfile)
        print("%%%%%%%%% MINOS -- NEUTRINO MODE %%%%%%%%%%", file=textfile)
        print("-- vmu --", file=textfile)
        print("\tMINOS: {:.3e} [neutrinos / m^2 POT]".format(MINOS_neutrino_vmu_integrated_flux), file=textfile)
        print("\tMINOS+: {:.3e} [neutrinos / m^2 POT]".format(MINOSPlus_neutrino_vmu_integrated_flux), file=textfile)
        print("-- vmubar --", file=textfile)
        print("\tMINOS: {:.3e} [neutrinos / m^2 POT]".format(MINOS_neutrino_vmubar_integrated_flux), file=textfile)
        print("\tMINOS+: {:.3e} [neutrinos / m^2 POT]".format(MINOSPlus_neutrino_vmubar_integrated_flux), file=textfile)
        print("%%%%%%% MINOS -- ANTINEUTRINO MODE %%%%%%%%", file=textfile)
        print("-- vmu --", file=textfile)
        print("\tMINOS: {:.3e} [neutrinos / m^2 POT]".format(MINOS_antineutrino_vmu_integrated_flux), file=textfile)
        print("-- vmubar --", file=textfile)
        print("\tMINOS: {:.3e} [neutrinos / m^2 POT]".format(MINOS_antineutrino_vmubar_integrated_flux), file=textfile)
        print("%%%%%%%%% T2K -- INGRID %%%%%%%%%%", file=textfile)
        print("-- vmu --", file=textfile)
        print("\tINGRID: {:.3e} [neutrinos / m^2 POT]".format(INGRID_neutrino_vmu_integrated_flux), file=textfile)
        print("-- vmubar --", file=textfile)
        print("\tINGRID: {:.3e} [neutrinos / m^2 POT]".format(INGRID_neutrino_vmubar_integrated_flux), file=textfile)
        print("-- ve --", file=textfile)
        print("\tINGRID: {:.3e} [neutrinos / m^2 POT]".format(INGRID_neutrino_ve_integrated_flux), file=textfile)
        print("-- vebar --", file=textfile)
        print("\tINGRID: {:.3e} [neutrinos / m^2 POT]".format(INGRID_neutrino_vebar_integrated_flux), file=textfile)
        print("%%%%%%%%%%%%%%%%%% SHiP %%%%%%%%%%%%%%%%%%", file=textfile)
        print("\tSHiP vmu Flux: {:.3e} [neutrinos / m^2 POT]".format(SHiP_vmu_integrated_flux), file=textfile)
        print("\tSHiP vmubar Flux: {:.3e} [neutrinos / m^2 POT]".format(SHiP_vmubar_integrated_flux), file=textfile)
        print("\tSHiP ve Flux: {:.3e} [neutrinos / m^2 POT]".format(SHiP_ve_integrated_flux), file=textfile)
        print("\tSHiP vebar Flux: {:.3e} [neutrinos / m^2 POT]".format(SHiP_vebar_integrated_flux), file=textfile)
        print("%%%%%%%%%%%%%%%%%% SBND %%%%%%%%%%%%%%%%%%", file=textfile)
        print("\tSBND vmu Flux: {:.3e} [neutrinos / m^2 POT]".format(SBND_vmu_integrated_flux), file=textfile)
        print("\tSBND vmubar Flux: {:.3e} [neutrinos / m^2 POT]".format(SBND_vmubar_integrated_flux), file=textfile)
        print("\tSBND ve Flux: {:.3e} [neutrinos / m^2 POT]".format(SBND_ve_integrated_flux), file=textfile)
        print("\tSBND vebar Flux: {:.3e} [neutrinos / m^2 POT]".format(SBND_vebar_integrated_flux), file=textfile)
        print("\n", file=textfile)

        print("%%%%%%%%%%%%%%%%%% DUNE -- NEUTRINO MODE %%%%%%%%%%%%%%%%%%", file=textfile)
        print("---------------------- vmu ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, XSecCon_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, upper_N_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, lower_N_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE EPA", 67, 3, N_EPA_Heaviside_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, XSecCon_EPA_Heaviside_1tau_coh_Ar_DUNE_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE EPA tau-opt", 67, 3, N_EPA_Heaviside_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_EPA_Heaviside_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3)
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_p_Ar_DUNE_neutrino_vmu_67_3, XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3, upper_N_1tau_p_Ar_DUNE_neutrino_vmu_67_3, lower_N_1tau_p_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3)
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_n_Ar_DUNE_neutrino_vmu_67_3, XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3, upper_N_1tau_n_Ar_DUNE_neutrino_vmu_67_3, lower_N_1tau_n_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3)
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_incoh_DUNE_neutrino_vmu_67_3, XSecCon_1tau_incoh_DUNE_neutrino_vmu_67_3, upper_N_1tau_incoh_DUNE_neutrino_vmu_67_3, lower_N_1tau_incoh_DUNE_neutrino_vmu_67_3, upper_XSecCon_1tau_incoh_DUNE_neutrino_vmu_67_3, lower_XSecCon_1tau_incoh_DUNE_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, upper_N_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, lower_N_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3)
        print("\n", file=textfile)

        print("2tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_coh_Ar_DUNE_neutrino_vmu_67_3, XSecCon_2tau_coh_Ar_DUNE_neutrino_vmu_67_3, upper_N_2tau_coh_Ar_DUNE_neutrino_vmu_67_3, lower_N_2tau_coh_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_N_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_N_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3)
        print("2tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_p_Ar_DUNE_neutrino_vmu_67_3, XSecCon_2tau_p_Ar_DUNE_neutrino_vmu_67_3, upper_N_2tau_p_Ar_DUNE_neutrino_vmu_67_3, lower_N_2tau_p_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_2tau_p_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_2tau_p_Ar_DUNE_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3)
        print("2tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_n_Ar_DUNE_neutrino_vmu_67_3, XSecCon_2tau_n_Ar_DUNE_neutrino_vmu_67_3, upper_N_2tau_n_Ar_DUNE_neutrino_vmu_67_3, lower_N_2tau_n_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_2tau_n_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_2tau_n_Ar_DUNE_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3)
        print("2tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_incoh_DUNE_neutrino_vmu_67_3, XSecCon_2tau_incoh_DUNE_neutrino_vmu_67_3, upper_N_2tau_incoh_DUNE_neutrino_vmu_67_3, lower_N_2tau_incoh_DUNE_neutrino_vmu_67_3, upper_XSecCon_2tau_incoh_DUNE_neutrino_vmu_67_3, lower_XSecCon_2tau_incoh_DUNE_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, upper_N_2tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, lower_N_2tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3)
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, XSecCon_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, upper_N_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, lower_N_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_N_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_N_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE", 67, 3, sum(Nbin_2mu_coh_Ar_DUNE_neutrino_vmu_67_3), XSecCon_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, sum(upper_Nbin_2mu_coh_Ar_DUNE_neutrino_vmu_67_3), sum(lower_Nbin_2mu_coh_Ar_DUNE_neutrino_vmu_67_3), upper_XSecCon_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE EPA", 67, 3, N_EPA_Heaviside_2mu_coh_Ar_DUNE_neutrino_vmu_67_3, XSecCon_EPA_Heaviside_2mu_coh_Ar_DUNE_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE EPA tau-opt", 67, 3, N_EPA_Heaviside_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_EPA_Heaviside_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3)
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_p_Ar_DUNE_neutrino_vmu_67_3, XSecCon_2mu_p_Ar_DUNE_neutrino_vmu_67_3, upper_N_2mu_p_Ar_DUNE_neutrino_vmu_67_3, lower_N_2mu_p_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_2mu_p_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_2mu_p_Ar_DUNE_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3)
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_n_Ar_DUNE_neutrino_vmu_67_3, XSecCon_2mu_n_Ar_DUNE_neutrino_vmu_67_3, upper_N_2mu_n_Ar_DUNE_neutrino_vmu_67_3, lower_N_2mu_n_Ar_DUNE_neutrino_vmu_67_3, upper_XSecCon_2mu_n_Ar_DUNE_neutrino_vmu_67_3, lower_XSecCon_2mu_n_Ar_DUNE_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3)
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_incoh_DUNE_neutrino_vmu_67_3, XSecCon_2mu_incoh_DUNE_neutrino_vmu_67_3, upper_N_2mu_incoh_DUNE_neutrino_vmu_67_3, lower_N_2mu_incoh_DUNE_neutrino_vmu_67_3, upper_XSecCon_2mu_incoh_DUNE_neutrino_vmu_67_3, lower_XSecCon_2mu_incoh_DUNE_neutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_incoh_DUNE_tau_opt_neutrino_vmu_67_3, XSecCon_2mu_incoh_DUNE_tau_opt_neutrino_vmu_67_3, upper_N_2mu_incoh_DUNE_tau_opt_neutrino_vmu_67_3, lower_N_2mu_incoh_DUNE_tau_opt_neutrino_vmu_67_3, upper_XSecCon_2mu_incoh_DUNE_tau_opt_neutrino_vmu_67_3, lower_XSecCon_2mu_incoh_DUNE_tau_opt_neutrino_vmu_67_3)
        print("\n", file=textfile)

        print("---------------------- vmubar ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_coh_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_1tau_coh_Ar_DUNE_neutrino_vmubar_67_3, upper_N_1tau_coh_Ar_DUNE_neutrino_vmubar_67_3, lower_N_1tau_coh_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_neutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_p_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_1tau_p_Ar_DUNE_neutrino_vmubar_67_3, upper_N_1tau_p_Ar_DUNE_neutrino_vmubar_67_3, lower_N_1tau_p_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_neutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_n_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_1tau_n_Ar_DUNE_neutrino_vmubar_67_3, upper_N_1tau_n_Ar_DUNE_neutrino_vmubar_67_3, lower_N_1tau_n_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_neutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_incoh_DUNE_neutrino_vmubar_67_3, XSecCon_1tau_incoh_DUNE_neutrino_vmubar_67_3, upper_N_1tau_incoh_DUNE_neutrino_vmubar_67_3, lower_N_1tau_incoh_DUNE_neutrino_vmubar_67_3, upper_XSecCon_1tau_incoh_DUNE_neutrino_vmubar_67_3, lower_XSecCon_1tau_incoh_DUNE_neutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, upper_N_1tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, lower_N_1tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3)
        print("\n", file=textfile)

        print("2tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_coh_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_2tau_coh_Ar_DUNE_neutrino_vmubar_67_3, upper_N_2tau_coh_Ar_DUNE_neutrino_vmubar_67_3, lower_N_2tau_coh_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_neutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_N_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_N_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)
        print("2tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_p_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_2tau_p_Ar_DUNE_neutrino_vmubar_67_3, upper_N_2tau_p_Ar_DUNE_neutrino_vmubar_67_3, lower_N_2tau_p_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2tau_p_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2tau_p_Ar_DUNE_neutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_N_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)
        print("2tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_n_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_2tau_n_Ar_DUNE_neutrino_vmubar_67_3, upper_N_2tau_n_Ar_DUNE_neutrino_vmubar_67_3, lower_N_2tau_n_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2tau_n_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2tau_n_Ar_DUNE_neutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_N_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)
        print("2tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_incoh_DUNE_neutrino_vmubar_67_3, XSecCon_2tau_incoh_DUNE_neutrino_vmubar_67_3, upper_N_2tau_incoh_DUNE_neutrino_vmubar_67_3, lower_N_2tau_incoh_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2tau_incoh_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2tau_incoh_DUNE_neutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, upper_N_2tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, lower_N_2tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2tau_incoh_DUNE_tau_opt_neutrino_vmubar_67_3)
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_coh_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_2mu_coh_Ar_DUNE_neutrino_vmubar_67_3, upper_N_2mu_coh_Ar_DUNE_neutrino_vmubar_67_3, lower_N_2mu_coh_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_neutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_N_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_N_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_p_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_2mu_p_Ar_DUNE_neutrino_vmubar_67_3, upper_N_2mu_p_Ar_DUNE_neutrino_vmubar_67_3, lower_N_2mu_p_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2mu_p_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2mu_p_Ar_DUNE_neutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_N_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_n_Ar_DUNE_neutrino_vmubar_67_3, XSecCon_2mu_n_Ar_DUNE_neutrino_vmubar_67_3, upper_N_2mu_n_Ar_DUNE_neutrino_vmubar_67_3, lower_N_2mu_n_Ar_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2mu_n_Ar_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2mu_n_Ar_DUNE_neutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_N_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_neutrino_vmubar_67_3)
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_incoh_DUNE_neutrino_vmubar_67_3, XSecCon_2mu_incoh_DUNE_neutrino_vmubar_67_3, upper_N_2mu_incoh_DUNE_neutrino_vmubar_67_3, lower_N_2mu_incoh_DUNE_neutrino_vmubar_67_3, upper_XSecCon_2mu_incoh_DUNE_neutrino_vmubar_67_3, lower_XSecCon_2mu_incoh_DUNE_neutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, XSecCon_2mu_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, upper_N_2mu_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, lower_N_2mu_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, upper_XSecCon_2mu_incoh_DUNE_tau_opt_neutrino_vmubar_67_3, lower_XSecCon_2mu_incoh_DUNE_tau_opt_neutrino_vmubar_67_3)
        print("\n", file=textfile)

        print("---------------------- ve ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_coh_Ar_DUNE_neutrino_ve_67_3, XSecCon_1tau_coh_Ar_DUNE_neutrino_ve_67_3, upper_N_1tau_coh_Ar_DUNE_neutrino_ve_67_3, lower_N_1tau_coh_Ar_DUNE_neutrino_ve_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_neutrino_ve_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_neutrino_ve_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_coh_Ar_DUNE_tau_opt_neutrino_ve_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_ve_67_3, upper_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_ve_67_3, lower_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_ve_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_ve_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_ve_67_3)
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_p_Ar_DUNE_neutrino_ve_67_3, XSecCon_1tau_p_Ar_DUNE_neutrino_ve_67_3, upper_N_1tau_p_Ar_DUNE_neutrino_ve_67_3, lower_N_1tau_p_Ar_DUNE_neutrino_ve_67_3, upper_XSecCon_1tau_p_Ar_DUNE_neutrino_ve_67_3, lower_XSecCon_1tau_p_Ar_DUNE_neutrino_ve_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3, upper_N_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3, lower_N_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_ve_67_3)
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_n_Ar_DUNE_neutrino_ve_67_3, XSecCon_1tau_n_Ar_DUNE_neutrino_ve_67_3, upper_N_1tau_n_Ar_DUNE_neutrino_ve_67_3, lower_N_1tau_n_Ar_DUNE_neutrino_ve_67_3, upper_XSecCon_1tau_n_Ar_DUNE_neutrino_ve_67_3, lower_XSecCon_1tau_n_Ar_DUNE_neutrino_ve_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3, upper_N_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3, lower_N_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_ve_67_3)
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_incoh_DUNE_neutrino_ve_67_3, XSecCon_1tau_incoh_DUNE_neutrino_ve_67_3, upper_N_1tau_incoh_DUNE_neutrino_ve_67_3, lower_N_1tau_incoh_DUNE_neutrino_ve_67_3, upper_XSecCon_1tau_incoh_DUNE_neutrino_ve_67_3, lower_XSecCon_1tau_incoh_DUNE_neutrino_ve_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_incoh_DUNE_tau_opt_neutrino_ve_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_ve_67_3, upper_N_1tau_incoh_DUNE_tau_opt_neutrino_ve_67_3, lower_N_1tau_incoh_DUNE_tau_opt_neutrino_ve_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_ve_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_ve_67_3)
        print("\n", file=textfile)

        print("---------------------- vebar ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_coh_Ar_DUNE_neutrino_vebar_67_3, XSecCon_1tau_coh_Ar_DUNE_neutrino_vebar_67_3, upper_N_1tau_coh_Ar_DUNE_neutrino_vebar_67_3, lower_N_1tau_coh_Ar_DUNE_neutrino_vebar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_neutrino_vebar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_neutrino_vebar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vebar_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vebar_67_3, upper_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vebar_67_3, lower_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vebar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vebar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vebar_67_3)
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_p_Ar_DUNE_neutrino_vebar_67_3, XSecCon_1tau_p_Ar_DUNE_neutrino_vebar_67_3, upper_N_1tau_p_Ar_DUNE_neutrino_vebar_67_3, lower_N_1tau_p_Ar_DUNE_neutrino_vebar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_neutrino_vebar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_neutrino_vebar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3, upper_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3, lower_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vebar_67_3)
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_n_Ar_DUNE_neutrino_vebar_67_3, XSecCon_1tau_n_Ar_DUNE_neutrino_vebar_67_3, upper_N_1tau_n_Ar_DUNE_neutrino_vebar_67_3, lower_N_1tau_n_Ar_DUNE_neutrino_vebar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_neutrino_vebar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_neutrino_vebar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3, upper_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3, lower_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vebar_67_3)
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_incoh_DUNE_neutrino_vebar_67_3, XSecCon_1tau_incoh_DUNE_neutrino_vebar_67_3, upper_N_1tau_incoh_DUNE_neutrino_vebar_67_3, lower_N_1tau_incoh_DUNE_neutrino_vebar_67_3, upper_XSecCon_1tau_incoh_DUNE_neutrino_vebar_67_3, lower_XSecCon_1tau_incoh_DUNE_neutrino_vebar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_incoh_DUNE_tau_opt_neutrino_vebar_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vebar_67_3, upper_N_1tau_incoh_DUNE_tau_opt_neutrino_vebar_67_3, lower_N_1tau_incoh_DUNE_tau_opt_neutrino_vebar_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vebar_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vebar_67_3)
        print("\n", file=textfile)

        print("%%%%%%%%%%%%%%%%% DUNE -- ANTINEUTRINO MODE %%%%%%%%%%%%%%%%%", file=textfile)
        print("---------------------- vmu ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3, upper_N_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3, lower_N_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3, upper_N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3, lower_N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3, upper_N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3, lower_N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_incoh_DUNE_antineutrino_vmu_67_3, XSecCon_1tau_incoh_DUNE_antineutrino_vmu_67_3, upper_N_1tau_incoh_DUNE_antineutrino_vmu_67_3, lower_N_1tau_incoh_DUNE_antineutrino_vmu_67_3, upper_XSecCon_1tau_incoh_DUNE_antineutrino_vmu_67_3, lower_XSecCon_1tau_incoh_DUNE_antineutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, upper_N_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, lower_N_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3)
        print("\n", file=textfile)

        print("2tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_coh_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_2tau_coh_Ar_DUNE_antineutrino_vmu_67_3, upper_N_2tau_coh_Ar_DUNE_antineutrino_vmu_67_3, lower_N_2tau_coh_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_antineutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_N_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_N_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)
        print("2tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_p_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_2tau_p_Ar_DUNE_antineutrino_vmu_67_3, upper_N_2tau_p_Ar_DUNE_antineutrino_vmu_67_3, lower_N_2tau_p_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2tau_p_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2tau_p_Ar_DUNE_antineutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)
        print("2tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_n_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_2tau_n_Ar_DUNE_antineutrino_vmu_67_3, upper_N_2tau_n_Ar_DUNE_antineutrino_vmu_67_3, lower_N_2tau_n_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2tau_n_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2tau_n_Ar_DUNE_antineutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)
        print("2tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_incoh_DUNE_antineutrino_vmu_67_3, XSecCon_2tau_incoh_DUNE_antineutrino_vmu_67_3, upper_N_2tau_incoh_DUNE_antineutrino_vmu_67_3, lower_N_2tau_incoh_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2tau_incoh_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2tau_incoh_DUNE_antineutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, upper_N_2tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, lower_N_2tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3)
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_coh_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_2mu_coh_Ar_DUNE_antineutrino_vmu_67_3, upper_N_2mu_coh_Ar_DUNE_antineutrino_vmu_67_3, lower_N_2mu_coh_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_antineutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_N_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_N_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_p_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_2mu_p_Ar_DUNE_antineutrino_vmu_67_3, upper_N_2mu_p_Ar_DUNE_antineutrino_vmu_67_3, lower_N_2mu_p_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2mu_p_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2mu_p_Ar_DUNE_antineutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_n_Ar_DUNE_antineutrino_vmu_67_3, XSecCon_2mu_n_Ar_DUNE_antineutrino_vmu_67_3, upper_N_2mu_n_Ar_DUNE_antineutrino_vmu_67_3, lower_N_2mu_n_Ar_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2mu_n_Ar_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2mu_n_Ar_DUNE_antineutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_incoh_DUNE_antineutrino_vmu_67_3, XSecCon_2mu_incoh_DUNE_antineutrino_vmu_67_3, upper_N_2mu_incoh_DUNE_antineutrino_vmu_67_3, lower_N_2mu_incoh_DUNE_antineutrino_vmu_67_3, upper_XSecCon_2mu_incoh_DUNE_antineutrino_vmu_67_3, lower_XSecCon_2mu_incoh_DUNE_antineutrino_vmu_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, XSecCon_2mu_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, upper_N_2mu_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, lower_N_2mu_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, upper_XSecCon_2mu_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, lower_XSecCon_2mu_incoh_DUNE_tau_opt_antineutrino_vmu_67_3)
        print("\n", file=textfile)

        print("vmuCC Cross Check", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_vmuCC_FASERvmu_1p2_3, XSecCon_vmuCC_FASERvmu_1p2_3, detector='FASER')
        print("\n", file=textfile)

        print("---------------------- vmubar ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, upper_N_1tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, lower_N_1tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3, upper_N_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3, lower_N_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3, upper_N_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3, lower_N_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_incoh_DUNE_antineutrino_vmubar_67_3, XSecCon_1tau_incoh_DUNE_antineutrino_vmubar_67_3, upper_N_1tau_incoh_DUNE_antineutrino_vmubar_67_3, lower_N_1tau_incoh_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_1tau_incoh_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_1tau_incoh_DUNE_antineutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_N_1tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_N_1tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3)
        print("\n", file=textfile)

        print("2tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_2tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, upper_N_2tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, lower_N_2tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_antineutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_N_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_N_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)
        print("2tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3, upper_N_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3, lower_N_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2tau_p_Ar_DUNE_antineutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_N_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)
        print("2tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3, upper_N_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3, lower_N_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2tau_n_Ar_DUNE_antineutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_N_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)
        print("2tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_incoh_DUNE_antineutrino_vmubar_67_3, XSecCon_2tau_incoh_DUNE_antineutrino_vmubar_67_3, upper_N_2tau_incoh_DUNE_antineutrino_vmubar_67_3, lower_N_2tau_incoh_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2tau_incoh_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2tau_incoh_DUNE_antineutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_N_2tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_N_2tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2tau_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3)
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_coh_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_2mu_coh_Ar_DUNE_antineutrino_vmubar_67_3, upper_N_2mu_coh_Ar_DUNE_antineutrino_vmubar_67_3, lower_N_2mu_coh_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_antineutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_N_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_N_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3, upper_N_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3, lower_N_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2mu_p_Ar_DUNE_antineutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_N_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3, XSecCon_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3, upper_N_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3, lower_N_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2mu_n_Ar_DUNE_antineutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_N_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_antineutrino_vmubar_67_3)
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_incoh_DUNE_antineutrino_vmubar_67_3, XSecCon_2mu_incoh_DUNE_antineutrino_vmubar_67_3, upper_N_2mu_incoh_DUNE_antineutrino_vmubar_67_3, lower_N_2mu_incoh_DUNE_antineutrino_vmubar_67_3, upper_XSecCon_2mu_incoh_DUNE_antineutrino_vmubar_67_3, lower_XSecCon_2mu_incoh_DUNE_antineutrino_vmubar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, XSecCon_2mu_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_N_2mu_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_N_2mu_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, upper_XSecCon_2mu_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3, lower_XSecCon_2mu_incoh_DUNE_tau_opt_antineutrino_vmubar_67_3)
        print("\n", file=textfile)

        print("---------------------- ve ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_coh_Ar_DUNE_antineutrino_ve_67_3, XSecCon_1tau_coh_Ar_DUNE_antineutrino_ve_67_3, upper_N_1tau_coh_Ar_DUNE_antineutrino_ve_67_3, lower_N_1tau_coh_Ar_DUNE_antineutrino_ve_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_antineutrino_ve_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_antineutrino_ve_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_ve_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_ve_67_3, upper_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_ve_67_3, lower_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_ve_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_ve_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_ve_67_3)
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_p_Ar_DUNE_antineutrino_ve_67_3, XSecCon_1tau_p_Ar_DUNE_antineutrino_ve_67_3, upper_N_1tau_p_Ar_DUNE_antineutrino_ve_67_3, lower_N_1tau_p_Ar_DUNE_antineutrino_ve_67_3, upper_XSecCon_1tau_p_Ar_DUNE_antineutrino_ve_67_3, lower_XSecCon_1tau_p_Ar_DUNE_antineutrino_ve_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3, upper_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3, lower_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_ve_67_3)
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_n_Ar_DUNE_antineutrino_ve_67_3, XSecCon_1tau_n_Ar_DUNE_antineutrino_ve_67_3, upper_N_1tau_n_Ar_DUNE_antineutrino_ve_67_3, lower_N_1tau_n_Ar_DUNE_antineutrino_ve_67_3, upper_XSecCon_1tau_n_Ar_DUNE_antineutrino_ve_67_3, lower_XSecCon_1tau_n_Ar_DUNE_antineutrino_ve_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3, upper_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3, lower_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_ve_67_3)
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_incoh_DUNE_antineutrino_ve_67_3, XSecCon_1tau_incoh_DUNE_antineutrino_ve_67_3, upper_N_1tau_incoh_DUNE_antineutrino_ve_67_3, lower_N_1tau_incoh_DUNE_antineutrino_ve_67_3, upper_XSecCon_1tau_incoh_DUNE_antineutrino_ve_67_3, lower_XSecCon_1tau_incoh_DUNE_antineutrino_ve_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_incoh_DUNE_tau_opt_antineutrino_ve_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_ve_67_3, upper_N_1tau_incoh_DUNE_tau_opt_antineutrino_ve_67_3, lower_N_1tau_incoh_DUNE_tau_opt_antineutrino_ve_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_ve_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_ve_67_3)
        print("\n", file=textfile)

        print("---------------------- vebar ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_coh_Ar_DUNE_antineutrino_vebar_67_3, XSecCon_1tau_coh_Ar_DUNE_antineutrino_vebar_67_3, upper_N_1tau_coh_Ar_DUNE_antineutrino_vebar_67_3, lower_N_1tau_coh_Ar_DUNE_antineutrino_vebar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vebar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vebar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, upper_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, lower_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vebar_67_3)
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_p_Ar_DUNE_antineutrino_vebar_67_3, XSecCon_1tau_p_Ar_DUNE_antineutrino_vebar_67_3, upper_N_1tau_p_Ar_DUNE_antineutrino_vebar_67_3, lower_N_1tau_p_Ar_DUNE_antineutrino_vebar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_antineutrino_vebar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_antineutrino_vebar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, upper_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, lower_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vebar_67_3)
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_n_Ar_DUNE_antineutrino_vebar_67_3, XSecCon_1tau_n_Ar_DUNE_antineutrino_vebar_67_3, upper_N_1tau_n_Ar_DUNE_antineutrino_vebar_67_3, lower_N_1tau_n_Ar_DUNE_antineutrino_vebar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_antineutrino_vebar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_antineutrino_vebar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, upper_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, lower_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vebar_67_3)
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_incoh_DUNE_antineutrino_vebar_67_3, XSecCon_1tau_incoh_DUNE_antineutrino_vebar_67_3, upper_N_1tau_incoh_DUNE_antineutrino_vebar_67_3, lower_N_1tau_incoh_DUNE_antineutrino_vebar_67_3, upper_XSecCon_1tau_incoh_DUNE_antineutrino_vebar_67_3, lower_XSecCon_1tau_incoh_DUNE_antineutrino_vebar_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_incoh_DUNE_tau_opt_antineutrino_vebar_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vebar_67_3, upper_N_1tau_incoh_DUNE_tau_opt_antineutrino_vebar_67_3, lower_N_1tau_incoh_DUNE_tau_opt_antineutrino_vebar_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vebar_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vebar_67_3)
        print("\n", file=textfile)

        print("%%%%%%%%%%%%%%%%%% DUNE -- NEUTRINO MODE %%%%%%%%%%%%%%%%%%", file=textfile)
        print("%%%%%%%%%%%%%%%%% with correction factors %%%%%%%%%%%%%%%%%", file=textfile)
        print("---------------------- vmu ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss)
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_p_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss)
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_n_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss)
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss)
        print("\n", file=textfile)

        print("%%%%%%%%%%%%%%%%% DUNE -- ANTINEUTRINO MODE %%%%%%%%%%%%%%%%%", file=textfile)
        print("%%%%%%%%%%%%%%%%% with correction factors %%%%%%%%%%%%%%%%%", file=textfile)
        print("---------------------- vmu ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss)
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_p_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss)
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_n_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss)
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_N_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_N_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, upper_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, lower_XSecCon_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss)
        print("\n", file=textfile)

        print("%%%%%%%%%%%%%%%%% FASER %%%%%%%%%%%%%%%%%", file=textfile)
        print("---------------------- vmu ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_coh_W_FASERvmu_1p2_3, XSecCon_1tau_coh_W_FASERvmu_1p2_3, upper_N_1tau_coh_W_FASERvmu_1p2_3, lower_N_1tau_coh_W_FASERvmu_1p2_3, upper_XSecCon_1tau_coh_W_FASERvmu_1p2_3, lower_XSecCon_1tau_coh_W_FASERvmu_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_coh_W_FASER2vmu, XSecCon_1tau_coh_W_FASER2vmu, upper_N_1tau_coh_W_FASER2vmu, lower_N_1tau_coh_W_FASER2vmu, upper_XSecCon_1tau_coh_W_FASER2vmu, lower_XSecCon_1tau_coh_W_FASER2vmu, detector='FASER2')
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_p_FASERvmu_1p2_3, XSecCon_1tau_p_FASERvmu_1p2_3, upper_N_1tau_p_W_FASERvmu_1p2_3, lower_N_1tau_p_W_FASERvmu_1p2_3, upper_XSecCon_1tau_p_W_FASERvmu_1p2_3, lower_XSecCon_1tau_p_W_FASERvmu_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_p_FASER2vmu, XSecCon_1tau_p_FASER2vmu, upper_N_1tau_p_W_FASER2vmu, lower_N_1tau_p_W_FASER2vmu, upper_XSecCon_1tau_p_W_FASER2vmu, lower_XSecCon_1tau_p_W_FASER2vmu, detector='FASER2')
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_n_FASERvmu_1p2_3, XSecCon_1tau_n_FASERvmu_1p2_3, upper_N_1tau_n_W_FASERvmu_1p2_3, lower_N_1tau_n_W_FASERvmu_1p2_3, upper_XSecCon_1tau_n_W_FASERvmu_1p2_3, lower_XSecCon_1tau_n_W_FASERvmu_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_n_FASER2vmu, XSecCon_1tau_n_FASER2vmu, upper_N_1tau_n_W_FASER2vmu, lower_N_1tau_n_W_FASER2vmu, upper_XSecCon_1tau_n_W_FASER2vmu, lower_XSecCon_1tau_n_W_FASER2vmu, detector='FASER2')
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_incoh_FASERvmu_1p2_3, XSecCon_1tau_incoh_FASERvmu_1p2_3, upper_N_1tau_incoh_FASERvmu_1p2_3, lower_N_1tau_incoh_FASERvmu_1p2_3, upper_XSecCon_1tau_incoh_FASERvmu_1p2_3, lower_XSecCon_1tau_incoh_FASERvmu_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_incoh_FASER2vmu, XSecCon_1tau_incoh_FASER2vmu, upper_N_1tau_incoh_FASER2vmu, lower_N_1tau_incoh_FASER2vmu, upper_XSecCon_1tau_incoh_FASER2vmu, lower_XSecCon_1tau_incoh_FASER2vmu, detector='FASER2')
        print("\n", file=textfile)


        print("2tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_coh_W_FASERvmu_1p2_3, XSecCon_2tau_coh_W_FASERvmu_1p2_3, upper_N_2tau_coh_W_FASERvmu_1p2_3, lower_N_2tau_coh_W_FASERvmu_1p2_3, upper_XSecCon_2tau_coh_W_FASERvmu_1p2_3, lower_XSecCon_2tau_coh_W_FASERvmu_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_coh_W_FASER2vmu, XSecCon_2tau_coh_W_FASER2vmu, upper_N_2tau_coh_W_FASER2vmu, lower_N_2tau_coh_W_FASER2vmu, upper_XSecCon_2tau_coh_W_FASER2vmu, lower_XSecCon_2tau_coh_W_FASER2vmu, detector='FASER2')
        print("2tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_p_FASERvmu_1p2_3, XSecCon_2tau_p_FASERvmu_1p2_3, upper_N_2tau_p_W_FASERvmu_1p2_3, lower_N_2tau_p_W_FASERvmu_1p2_3, upper_XSecCon_2tau_p_W_FASERvmu_1p2_3, lower_XSecCon_2tau_p_W_FASERvmu_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_p_FASER2vmu, XSecCon_2tau_p_FASER2vmu, upper_N_2tau_p_W_FASER2vmu, lower_N_2tau_p_W_FASER2vmu, upper_XSecCon_2tau_p_W_FASER2vmu, lower_XSecCon_2tau_p_W_FASER2vmu, detector='FASER2')
        print("2tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_n_FASERvmu_1p2_3, XSecCon_2tau_n_FASERvmu_1p2_3, upper_N_2tau_n_W_FASERvmu_1p2_3, lower_N_2tau_n_W_FASERvmu_1p2_3, upper_XSecCon_2tau_n_W_FASERvmu_1p2_3, lower_XSecCon_2tau_n_W_FASERvmu_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_n_FASER2vmu, XSecCon_2tau_n_FASER2vmu, upper_N_2tau_n_W_FASER2vmu, lower_N_2tau_n_W_FASER2vmu, upper_XSecCon_2tau_n_W_FASER2vmu, lower_XSecCon_2tau_n_W_FASER2vmu, detector='FASER2')
        print("2tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_incoh_FASERvmu_1p2_3, XSecCon_2tau_incoh_FASERvmu_1p2_3, upper_N_2tau_incoh_FASERvmu_1p2_3, lower_N_2tau_incoh_FASERvmu_1p2_3, upper_XSecCon_2tau_incoh_FASERvmu_1p2_3, lower_XSecCon_2tau_incoh_FASERvmu_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_incoh_FASER2vmu, XSecCon_2tau_incoh_FASER2vmu, upper_N_2tau_incoh_FASER2vmu, lower_N_2tau_incoh_FASER2vmu, upper_XSecCon_2tau_incoh_FASER2vmu, lower_XSecCon_2tau_incoh_FASER2vmu, detector='FASER2')
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_coh_W_FASERvmu_1p2_3, XSecCon_2mu_coh_W_FASERvmu_1p2_3, upper_N_2mu_coh_W_FASERvmu_1p2_3, lower_N_2mu_coh_W_FASERvmu_1p2_3, upper_XSecCon_2mu_coh_W_FASERvmu_1p2_3, lower_XSecCon_2mu_coh_W_FASERvmu_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_coh_W_FASER2vmu, XSecCon_2mu_coh_W_FASER2vmu, upper_N_2mu_coh_W_FASER2vmu, lower_N_2mu_coh_W_FASER2vmu, upper_XSecCon_2mu_coh_W_FASER2vmu, lower_XSecCon_2mu_coh_W_FASER2vmu, detector='FASER2')
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_p_FASERvmu_1p2_3, XSecCon_2mu_p_FASERvmu_1p2_3, upper_N_2mu_p_W_FASERvmu_1p2_3, lower_N_2mu_p_W_FASERvmu_1p2_3, upper_XSecCon_2mu_p_W_FASERvmu_1p2_3, lower_XSecCon_2mu_p_W_FASERvmu_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_p_FASER2vmu, XSecCon_2mu_p_FASER2vmu, upper_N_2mu_p_W_FASER2vmu, lower_N_2mu_p_W_FASER2vmu, upper_XSecCon_2mu_p_W_FASER2vmu, lower_XSecCon_2mu_p_W_FASER2vmu, detector='FASER2')
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_n_FASERvmu_1p2_3, XSecCon_2mu_n_FASERvmu_1p2_3, upper_N_2mu_n_W_FASERvmu_1p2_3, lower_N_2mu_n_W_FASERvmu_1p2_3, upper_XSecCon_2mu_n_W_FASERvmu_1p2_3, lower_XSecCon_2mu_n_W_FASERvmu_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_n_FASER2vmu, XSecCon_2mu_n_FASER2vmu, upper_N_2mu_n_W_FASER2vmu, lower_N_2mu_n_W_FASER2vmu, upper_XSecCon_2mu_n_W_FASER2vmu, lower_XSecCon_2mu_n_W_FASER2vmu, detector='FASER2')
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_incoh_FASERvmu_1p2_3, XSecCon_2mu_incoh_FASERvmu_1p2_3, upper_N_2mu_incoh_FASERvmu_1p2_3, lower_N_2mu_incoh_FASERvmu_1p2_3, upper_XSecCon_2mu_incoh_FASERvmu_1p2_3, lower_XSecCon_2mu_incoh_FASERvmu_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_incoh_FASER2vmu, XSecCon_2mu_incoh_FASER2vmu, upper_N_2mu_incoh_FASER2vmu, lower_N_2mu_incoh_FASER2vmu, upper_XSecCon_2mu_incoh_FASER2vmu, lower_XSecCon_2mu_incoh_FASER2vmu, detector='FASER2')
        print("\n", file=textfile)

        print("vmuCC Cross Check", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_vmuCC_FASERvmu_1p2_3, XSecCon_vmuCC_FASERvmu_1p2_3, detector='FASER')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_DIS_vmuCC_FASERvmu, 1, detector='FASER')
        print("\n", file=textfile)

        print("---------------------- vmubar ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_coh_W_FASERvmubar_1p2_3, XSecCon_1tau_coh_W_FASERvmubar_1p2_3, upper_N_1tau_coh_W_FASERvmubar_1p2_3, lower_N_1tau_coh_W_FASERvmubar_1p2_3, upper_XSecCon_1tau_coh_W_FASERvmubar_1p2_3, lower_XSecCon_1tau_coh_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_coh_W_FASER2vmubar, XSecCon_1tau_coh_W_FASER2vmubar, upper_N_1tau_coh_W_FASER2vmubar, lower_N_1tau_coh_W_FASER2vmubar, upper_XSecCon_1tau_coh_W_FASER2vmubar, lower_XSecCon_1tau_coh_W_FASER2vmubar, detector='FASER2')
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_p_FASERvmubar_1p2_3, XSecCon_1tau_p_FASERvmubar_1p2_3, upper_N_1tau_p_W_FASERvmubar_1p2_3, lower_N_1tau_p_W_FASERvmubar_1p2_3, upper_XSecCon_1tau_p_W_FASERvmubar_1p2_3, lower_XSecCon_1tau_p_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_p_FASER2vmubar, XSecCon_1tau_p_FASER2vmubar, upper_N_1tau_p_W_FASER2vmubar, lower_N_1tau_p_W_FASER2vmubar, upper_XSecCon_1tau_p_W_FASER2vmubar, lower_XSecCon_1tau_p_W_FASER2vmubar, detector='FASER2')
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_n_FASERvmubar_1p2_3, XSecCon_1tau_n_FASERvmubar_1p2_3, upper_N_1tau_n_W_FASERvmubar_1p2_3, lower_N_1tau_n_W_FASERvmubar_1p2_3, upper_XSecCon_1tau_n_W_FASERvmubar_1p2_3, lower_XSecCon_1tau_n_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_n_FASER2vmubar, XSecCon_1tau_n_FASER2vmubar, upper_N_1tau_n_W_FASER2vmubar, lower_N_1tau_n_W_FASER2vmubar, upper_XSecCon_1tau_n_W_FASER2vmubar, lower_XSecCon_1tau_n_W_FASER2vmubar, detector='FASER2')
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_incoh_FASERvmubar_1p2_3, XSecCon_1tau_incoh_FASERvmubar_1p2_3, upper_N_1tau_incoh_FASERvmubar_1p2_3, lower_N_1tau_incoh_FASERvmubar_1p2_3, upper_XSecCon_1tau_incoh_FASERvmubar_1p2_3, lower_XSecCon_1tau_incoh_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_incoh_FASER2vmubar, XSecCon_1tau_incoh_FASER2vmubar, upper_N_1tau_incoh_FASER2vmubar, lower_N_1tau_incoh_FASER2vmubar, upper_XSecCon_1tau_incoh_FASER2vmubar, lower_XSecCon_1tau_incoh_FASER2vmubar, detector='FASER2')
        print("\n", file=textfile)

        print("2tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_coh_W_FASERvmubar_1p2_3, XSecCon_2tau_coh_W_FASERvmubar_1p2_3, upper_N_2tau_coh_W_FASERvmubar_1p2_3, lower_N_2tau_coh_W_FASERvmubar_1p2_3, upper_XSecCon_2tau_coh_W_FASERvmubar_1p2_3, lower_XSecCon_2tau_coh_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_coh_W_FASER2vmubar, XSecCon_2tau_coh_W_FASER2vmubar, upper_N_2tau_coh_W_FASER2vmubar, lower_N_2tau_coh_W_FASER2vmubar, upper_XSecCon_2tau_coh_W_FASER2vmubar, lower_XSecCon_2tau_coh_W_FASER2vmubar, detector='FASER2')
        print("2tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_p_FASERvmubar_1p2_3, XSecCon_2tau_p_FASERvmubar_1p2_3, upper_N_2tau_p_W_FASERvmubar_1p2_3, lower_N_2tau_p_W_FASERvmubar_1p2_3, upper_XSecCon_2tau_p_W_FASERvmubar_1p2_3, lower_XSecCon_2tau_p_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_p_FASER2vmubar, XSecCon_2tau_p_FASER2vmubar, upper_N_2tau_p_W_FASER2vmubar, lower_N_2tau_p_W_FASER2vmubar, upper_XSecCon_2tau_p_W_FASER2vmubar, lower_XSecCon_2tau_p_W_FASER2vmubar, detector='FASER2')
        print("2tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_n_FASERvmubar_1p2_3, XSecCon_2tau_n_FASERvmubar_1p2_3, upper_N_2tau_n_W_FASERvmubar_1p2_3, lower_N_2tau_n_W_FASERvmubar_1p2_3, upper_XSecCon_2tau_n_W_FASERvmubar_1p2_3, lower_XSecCon_2tau_n_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_n_FASER2vmubar, XSecCon_2tau_n_FASER2vmubar, upper_N_2tau_n_W_FASER2vmubar, lower_N_2tau_n_W_FASER2vmubar, upper_XSecCon_2tau_n_W_FASER2vmubar, lower_XSecCon_2tau_n_W_FASER2vmubar, detector='FASER2')
        print("2tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_incoh_FASERvmubar_1p2_3, XSecCon_2tau_incoh_FASERvmubar_1p2_3, upper_N_2tau_incoh_FASERvmubar_1p2_3, lower_N_2tau_incoh_FASERvmubar_1p2_3, upper_XSecCon_2tau_incoh_FASERvmubar_1p2_3, lower_XSecCon_2tau_incoh_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_incoh_FASER2vmubar, XSecCon_2tau_incoh_FASER2vmubar, upper_N_2tau_incoh_FASER2vmubar, lower_N_2tau_incoh_FASER2vmubar, upper_XSecCon_2tau_incoh_FASER2vmubar, lower_XSecCon_2tau_incoh_FASER2vmubar, detector='FASER2')
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_coh_W_FASERvmubar_1p2_3, XSecCon_2mu_coh_W_FASERvmubar_1p2_3, upper_N_2mu_coh_W_FASERvmubar_1p2_3, lower_N_2mu_coh_W_FASERvmubar_1p2_3, upper_XSecCon_2mu_coh_W_FASERvmubar_1p2_3, lower_XSecCon_2mu_coh_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_coh_W_FASER2vmubar, XSecCon_2mu_coh_W_FASER2vmubar, upper_N_2mu_coh_W_FASER2vmubar, lower_N_2mu_coh_W_FASER2vmubar, upper_XSecCon_2mu_coh_W_FASER2vmubar, lower_XSecCon_2mu_coh_W_FASER2vmubar, detector='FASER2')
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_p_FASERvmubar_1p2_3, XSecCon_2mu_p_FASERvmubar_1p2_3, upper_N_2mu_p_W_FASERvmubar_1p2_3, lower_N_2mu_p_W_FASERvmubar_1p2_3, upper_XSecCon_2mu_p_W_FASERvmubar_1p2_3, lower_XSecCon_2mu_p_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_p_FASER2vmubar, XSecCon_2mu_p_FASER2vmubar, upper_N_2mu_p_W_FASER2vmubar, lower_N_2mu_p_W_FASER2vmubar, upper_XSecCon_2mu_p_W_FASER2vmubar, lower_XSecCon_2mu_p_W_FASER2vmubar, detector='FASER2')
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_n_FASERvmubar_1p2_3, XSecCon_2mu_n_FASERvmubar_1p2_3, upper_N_2mu_n_W_FASERvmubar_1p2_3, lower_N_2mu_n_W_FASERvmubar_1p2_3, upper_XSecCon_2mu_n_W_FASERvmubar_1p2_3, lower_XSecCon_2mu_n_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_n_FASER2vmubar, XSecCon_2mu_n_FASER2vmubar, upper_N_2mu_n_W_FASER2vmubar, lower_N_2mu_n_W_FASER2vmubar, upper_XSecCon_2mu_n_W_FASER2vmubar, lower_XSecCon_2mu_n_W_FASER2vmubar, detector='FASER2')
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_incoh_FASERvmubar_1p2_3, XSecCon_2mu_incoh_FASERvmubar_1p2_3, upper_N_2mu_incoh_FASERvmubar_1p2_3, lower_N_2mu_incoh_FASERvmubar_1p2_3, upper_XSecCon_2mu_incoh_FASERvmubar_1p2_3, lower_XSecCon_2mu_incoh_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_incoh_FASER2vmubar, XSecCon_2mu_incoh_FASER2vmubar, upper_N_2mu_incoh_FASER2vmubar, lower_N_2mu_incoh_FASER2vmubar, upper_XSecCon_2mu_incoh_FASER2vmubar, lower_XSecCon_2mu_incoh_FASER2vmubar, detector='FASER2')
        print("\n", file=textfile)

        print("vmuCC Cross Check", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_vmuCC_FASERvmubar_1p2_3, XSecCon_vmuCC_FASERvmubar_1p2_3, detector='FASER')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_DIS_vmubarCC_FASERvmubar, 1, detector='FASER')
        print("\n", file=textfile)

        print("---------------------- vmu + vmubar ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_coh_W_FASERvmu_1p2_3 + N_1tau_coh_W_FASERvmubar_1p2_3, XSecCon_1tau_coh_W_FASERvmu_1p2_3 + XSecCon_1tau_coh_W_FASERvmubar_1p2_3, upper_N_1tau_coh_W_FASERvmu_1p2_3 + upper_N_1tau_coh_W_FASERvmubar_1p2_3, lower_N_1tau_coh_W_FASERvmu_1p2_3 + lower_N_1tau_coh_W_FASERvmubar_1p2_3, upper_XSecCon_1tau_coh_W_FASERvmu_1p2_3 + upper_XSecCon_1tau_coh_W_FASERvmubar_1p2_3, lower_XSecCon_1tau_coh_W_FASERvmu_1p2_3 + lower_XSecCon_1tau_coh_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_coh_W_FASER2vmu + N_1tau_coh_W_FASER2vmubar, XSecCon_1tau_coh_W_FASER2vmu + XSecCon_1tau_coh_W_FASER2vmubar, upper_N_1tau_coh_W_FASER2vmu + upper_N_1tau_coh_W_FASER2vmubar, lower_N_1tau_coh_W_FASER2vmu + lower_N_1tau_coh_W_FASER2vmubar, upper_XSecCon_1tau_coh_W_FASER2vmu + upper_XSecCon_1tau_coh_W_FASER2vmubar, lower_XSecCon_1tau_coh_W_FASER2vmu + lower_XSecCon_1tau_coh_W_FASER2vmubar, detector='FASER2')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_coh_W_FASER2vmu, XSecCon_1tau_coh_W_FASER2vmu, upper_N_1tau_coh_W_FASER2vmu, lower_N_1tau_coh_W_FASER2vmu, upper_XSecCon_1tau_coh_W_FASER2vmu, lower_XSecCon_1tau_coh_W_FASER2vmu, detector='FASER2')
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_p_FASERvmu_1p2_3 + N_1tau_p_FASERvmubar_1p2_3, XSecCon_1tau_p_FASERvmu_1p2_3 + XSecCon_1tau_p_FASERvmubar_1p2_3, upper_N_1tau_p_W_FASERvmu_1p2_3 + upper_N_1tau_p_W_FASERvmubar_1p2_3, lower_N_1tau_p_W_FASERvmu_1p2_3 + lower_N_1tau_p_W_FASERvmubar_1p2_3, upper_XSecCon_1tau_p_W_FASERvmu_1p2_3 + upper_XSecCon_1tau_p_W_FASERvmubar_1p2_3, lower_XSecCon_1tau_p_W_FASERvmu_1p2_3 + lower_XSecCon_1tau_p_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_p_FASER2vmu + N_1tau_p_FASER2vmubar, XSecCon_1tau_p_FASER2vmu + XSecCon_1tau_p_FASER2vmubar, upper_N_1tau_p_W_FASER2vmu + upper_N_1tau_p_W_FASER2vmubar, lower_N_1tau_p_W_FASER2vmu + lower_N_1tau_p_W_FASER2vmubar, upper_XSecCon_1tau_p_W_FASER2vmu + upper_XSecCon_1tau_p_W_FASER2vmubar, lower_XSecCon_1tau_p_W_FASER2vmu + lower_XSecCon_1tau_p_W_FASER2vmubar, detector='FASER2')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_p_FASER2vmu, XSecCon_1tau_p_FASER2vmu, upper_N_1tau_p_W_FASER2vmu, lower_N_1tau_p_W_FASER2vmu, upper_XSecCon_1tau_p_W_FASER2vmu, lower_XSecCon_1tau_p_W_FASER2vmu, detector='FASER2')
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_n_FASERvmu_1p2_3 + N_1tau_n_FASERvmubar_1p2_3, XSecCon_1tau_n_FASERvmu_1p2_3 + XSecCon_1tau_n_FASERvmubar_1p2_3, upper_N_1tau_n_W_FASERvmu_1p2_3 + upper_N_1tau_n_W_FASERvmubar_1p2_3, lower_N_1tau_n_W_FASERvmu_1p2_3 + lower_N_1tau_n_W_FASERvmubar_1p2_3, upper_XSecCon_1tau_n_W_FASERvmu_1p2_3 + upper_XSecCon_1tau_n_W_FASERvmubar_1p2_3, lower_XSecCon_1tau_n_W_FASERvmu_1p2_3 + lower_XSecCon_1tau_n_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_n_FASER2vmu + N_1tau_n_FASER2vmubar, XSecCon_1tau_n_FASER2vmu + XSecCon_1tau_n_FASER2vmubar, upper_N_1tau_n_W_FASER2vmu + upper_N_1tau_n_W_FASER2vmubar, lower_N_1tau_n_W_FASER2vmu + lower_N_1tau_n_W_FASER2vmubar, upper_XSecCon_1tau_n_W_FASER2vmu + upper_XSecCon_1tau_n_W_FASER2vmubar, lower_XSecCon_1tau_n_W_FASER2vmu + lower_XSecCon_1tau_n_W_FASER2vmubar, detector='FASER2')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_n_FASER2vmu, XSecCon_1tau_n_FASER2vmu, upper_N_1tau_n_W_FASER2vmu, lower_N_1tau_n_W_FASER2vmu, upper_XSecCon_1tau_n_W_FASER2vmu, lower_XSecCon_1tau_n_W_FASER2vmu, detector='FASER2')
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_incoh_FASERvmu_1p2_3 + N_1tau_incoh_FASERvmubar_1p2_3, XSecCon_1tau_incoh_FASERvmu_1p2_3 + XSecCon_1tau_incoh_FASERvmubar_1p2_3, upper_N_1tau_incoh_FASERvmu_1p2_3 + upper_N_1tau_incoh_FASERvmubar_1p2_3, lower_N_1tau_incoh_FASERvmu_1p2_3 + lower_N_1tau_incoh_FASERvmubar_1p2_3, upper_XSecCon_1tau_incoh_FASERvmu_1p2_3 + upper_XSecCon_1tau_incoh_FASERvmubar_1p2_3, lower_XSecCon_1tau_incoh_FASERvmu_1p2_3 + lower_XSecCon_1tau_incoh_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_incoh_FASER2vmu + N_1tau_incoh_FASER2vmubar, XSecCon_1tau_incoh_FASER2vmu + XSecCon_1tau_incoh_FASER2vmubar, upper_N_1tau_incoh_FASER2vmu + upper_N_1tau_incoh_FASER2vmubar, lower_N_1tau_incoh_FASER2vmu + lower_N_1tau_incoh_FASER2vmubar, upper_XSecCon_1tau_incoh_FASER2vmu + upper_XSecCon_1tau_incoh_FASER2vmubar, lower_XSecCon_1tau_incoh_FASER2vmu + lower_XSecCon_1tau_incoh_FASER2vmubar, detector='FASER2')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_incoh_FASER2vmu, XSecCon_1tau_incoh_FASER2vmu, upper_N_1tau_incoh_FASER2vmu, lower_N_1tau_incoh_FASER2vmu, upper_XSecCon_1tau_incoh_FASER2vmu, lower_XSecCon_1tau_incoh_FASER2vmu, detector='FASER2')
        print("\n", file=textfile)

        print("2tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_coh_W_FASERvmu_1p2_3 + N_2tau_coh_W_FASERvmubar_1p2_3, XSecCon_2tau_coh_W_FASERvmu_1p2_3 + XSecCon_2tau_coh_W_FASERvmubar_1p2_3, upper_N_2tau_coh_W_FASERvmu_1p2_3 + upper_N_2tau_coh_W_FASERvmubar_1p2_3, lower_N_2tau_coh_W_FASERvmu_1p2_3 + lower_N_2tau_coh_W_FASERvmubar_1p2_3, upper_XSecCon_2tau_coh_W_FASERvmu_1p2_3 + upper_XSecCon_2tau_coh_W_FASERvmubar_1p2_3, lower_XSecCon_2tau_coh_W_FASERvmu_1p2_3 + lower_XSecCon_2tau_coh_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_coh_W_FASER2vmu + N_2tau_coh_W_FASER2vmubar, XSecCon_2tau_coh_W_FASER2vmu + XSecCon_2tau_coh_W_FASER2vmubar, upper_N_2tau_coh_W_FASER2vmu + upper_N_2tau_coh_W_FASER2vmubar, lower_N_2tau_coh_W_FASER2vmu + lower_N_2tau_coh_W_FASER2vmubar, upper_XSecCon_2tau_coh_W_FASER2vmu + upper_XSecCon_2tau_coh_W_FASER2vmubar, lower_XSecCon_2tau_coh_W_FASER2vmu + lower_XSecCon_2tau_coh_W_FASER2vmubar, detector='FASER2')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_coh_W_FASER2vmu, XSecCon_2tau_coh_W_FASER2vmu, upper_N_2tau_coh_W_FASER2vmu, lower_N_2tau_coh_W_FASER2vmu, upper_XSecCon_2tau_coh_W_FASER2vmu, lower_XSecCon_2tau_coh_W_FASER2vmu, detector='FASER2')
        print("2tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_p_FASERvmu_1p2_3 + N_2tau_p_FASERvmubar_1p2_3, XSecCon_2tau_p_FASERvmu_1p2_3 + XSecCon_2tau_p_FASERvmubar_1p2_3, upper_N_2tau_p_W_FASERvmu_1p2_3 + upper_N_2tau_p_W_FASERvmubar_1p2_3, lower_N_2tau_p_W_FASERvmu_1p2_3 + lower_N_2tau_p_W_FASERvmubar_1p2_3, upper_XSecCon_2tau_p_W_FASERvmu_1p2_3 + upper_XSecCon_2tau_p_W_FASERvmubar_1p2_3, lower_XSecCon_2tau_p_W_FASERvmu_1p2_3 + lower_XSecCon_2tau_p_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_p_FASER2vmu + N_2tau_p_FASER2vmubar, XSecCon_2tau_p_FASER2vmu + XSecCon_2tau_p_FASER2vmubar, upper_N_2tau_p_W_FASER2vmu + upper_N_2tau_p_W_FASER2vmubar, lower_N_2tau_p_W_FASER2vmu + lower_N_2tau_p_W_FASER2vmubar, upper_XSecCon_2tau_p_W_FASER2vmu + upper_XSecCon_2tau_p_W_FASER2vmubar, lower_XSecCon_2tau_p_W_FASER2vmu + lower_XSecCon_2tau_p_W_FASER2vmubar, detector='FASER2')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_p_FASER2vmu, XSecCon_2tau_p_FASER2vmu, upper_N_2tau_p_W_FASER2vmu, lower_N_2tau_p_W_FASER2vmu, upper_XSecCon_2tau_p_W_FASER2vmu, lower_XSecCon_2tau_p_W_FASER2vmu, detector='FASER2')
        print("2tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_n_FASERvmu_1p2_3 + N_2tau_n_FASERvmubar_1p2_3, XSecCon_2tau_n_FASERvmu_1p2_3 + XSecCon_2tau_n_FASERvmubar_1p2_3, upper_N_2tau_n_W_FASERvmu_1p2_3 + upper_N_2tau_n_W_FASERvmubar_1p2_3, lower_N_2tau_n_W_FASERvmu_1p2_3 + lower_N_2tau_n_W_FASERvmubar_1p2_3, upper_XSecCon_2tau_n_W_FASERvmu_1p2_3 + upper_XSecCon_2tau_n_W_FASERvmubar_1p2_3, lower_XSecCon_2tau_n_W_FASERvmu_1p2_3 + lower_XSecCon_2tau_n_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_n_FASER2vmu + N_2tau_n_FASER2vmubar, XSecCon_2tau_n_FASER2vmu + XSecCon_2tau_n_FASER2vmubar, upper_N_2tau_n_W_FASER2vmu + upper_N_2tau_n_W_FASER2vmubar, lower_N_2tau_n_W_FASER2vmu + lower_N_2tau_n_W_FASER2vmubar, upper_XSecCon_2tau_n_W_FASER2vmu + upper_XSecCon_2tau_n_W_FASER2vmubar, lower_XSecCon_2tau_n_W_FASER2vmu + lower_XSecCon_2tau_n_W_FASER2vmubar, detector='FASER2')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_n_FASER2vmu, XSecCon_2tau_n_FASER2vmu, upper_N_2tau_n_W_FASER2vmu, lower_N_2tau_n_W_FASER2vmu, upper_XSecCon_2tau_n_W_FASER2vmu, lower_XSecCon_2tau_n_W_FASER2vmu, detector='FASER2')
        print("2tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_incoh_FASERvmu_1p2_3 + N_2tau_incoh_FASERvmubar_1p2_3, XSecCon_2tau_incoh_FASERvmu_1p2_3 + XSecCon_2tau_incoh_FASERvmubar_1p2_3, upper_N_2tau_incoh_FASERvmu_1p2_3 + upper_N_2tau_incoh_FASERvmubar_1p2_3, lower_N_2tau_incoh_FASERvmu_1p2_3 + lower_N_2tau_incoh_FASERvmubar_1p2_3, upper_XSecCon_2tau_incoh_FASERvmu_1p2_3 + upper_XSecCon_2tau_incoh_FASERvmubar_1p2_3, lower_XSecCon_2tau_incoh_FASERvmu_1p2_3 + lower_XSecCon_2tau_incoh_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_incoh_FASER2vmu + N_2tau_incoh_FASER2vmubar, XSecCon_2tau_incoh_FASER2vmu + XSecCon_2tau_incoh_FASER2vmubar, upper_N_2tau_incoh_FASER2vmu + upper_N_2tau_incoh_FASER2vmubar, lower_N_2tau_incoh_FASER2vmu + lower_N_2tau_incoh_FASER2vmubar, upper_XSecCon_2tau_incoh_FASER2vmu + upper_XSecCon_2tau_incoh_FASER2vmubar, lower_XSecCon_2tau_incoh_FASER2vmu + lower_XSecCon_2tau_incoh_FASER2vmubar, detector='FASER2')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_incoh_FASER2vmu, XSecCon_2tau_incoh_FASER2vmu, upper_N_2tau_incoh_FASER2vmu, lower_N_2tau_incoh_FASER2vmu, upper_XSecCon_2tau_incoh_FASER2vmu, lower_XSecCon_2tau_incoh_FASER2vmu, detector='FASER2')
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_coh_W_FASERvmu_1p2_3 +  N_2mu_coh_W_FASERvmubar_1p2_3, XSecCon_2mu_coh_W_FASERvmu_1p2_3 + XSecCon_2mu_coh_W_FASERvmubar_1p2_3, upper_N_2mu_coh_W_FASERvmu_1p2_3 + upper_N_2mu_coh_W_FASERvmubar_1p2_3, lower_N_2mu_coh_W_FASERvmu_1p2_3 + lower_N_2mu_coh_W_FASERvmubar_1p2_3, upper_XSecCon_2mu_coh_W_FASERvmu_1p2_3 + upper_XSecCon_2mu_coh_W_FASERvmubar_1p2_3, lower_XSecCon_2mu_coh_W_FASERvmu_1p2_3 + lower_XSecCon_2mu_coh_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_coh_W_FASER2vmu +  N_2mu_coh_W_FASER2vmubar, XSecCon_2mu_coh_W_FASER2vmu + XSecCon_2mu_coh_W_FASER2vmubar, upper_N_2mu_coh_W_FASER2vmu + upper_N_2mu_coh_W_FASER2vmubar, lower_N_2mu_coh_W_FASER2vmu + lower_N_2mu_coh_W_FASER2vmubar, upper_XSecCon_2mu_coh_W_FASER2vmu + upper_XSecCon_2mu_coh_W_FASER2vmubar, lower_XSecCon_2mu_coh_W_FASER2vmu + lower_XSecCon_2mu_coh_W_FASER2vmubar, detector='FASER2')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_coh_W_FASER2vmu, XSecCon_2mu_coh_W_FASER2vmu, upper_N_2mu_coh_W_FASER2vmu, lower_N_2mu_coh_W_FASER2vmu, upper_XSecCon_2mu_coh_W_FASER2vmu, lower_XSecCon_2mu_coh_W_FASER2vmu, detector='FASER2')
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_p_FASERvmu_1p2_3 + N_2mu_p_FASERvmubar_1p2_3, XSecCon_2mu_p_FASERvmu_1p2_3 + XSecCon_2mu_p_FASERvmubar_1p2_3, upper_N_2mu_p_W_FASERvmu_1p2_3 + upper_N_2mu_p_W_FASERvmubar_1p2_3, lower_N_2mu_p_W_FASERvmu_1p2_3 + lower_N_2mu_p_W_FASERvmubar_1p2_3, upper_XSecCon_2mu_p_W_FASERvmu_1p2_3 + upper_XSecCon_2mu_p_W_FASERvmubar_1p2_3, lower_XSecCon_2mu_p_W_FASERvmu_1p2_3 + lower_XSecCon_2mu_p_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_p_FASER2vmu + N_2mu_p_FASER2vmubar, XSecCon_2mu_p_FASER2vmu + XSecCon_2mu_p_FASER2vmubar, upper_N_2mu_p_W_FASER2vmu + upper_N_2mu_p_W_FASER2vmubar, lower_N_2mu_p_W_FASER2vmu + lower_N_2mu_p_W_FASER2vmubar, upper_XSecCon_2mu_p_W_FASER2vmu + upper_XSecCon_2mu_p_W_FASER2vmubar, lower_XSecCon_2mu_p_W_FASER2vmu + lower_XSecCon_2mu_p_W_FASER2vmubar, detector='FASER2')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_p_FASER2vmu, XSecCon_2mu_p_FASER2vmu, upper_N_2mu_p_W_FASER2vmu, lower_N_2mu_p_W_FASER2vmu, upper_XSecCon_2mu_p_W_FASER2vmu, lower_XSecCon_2mu_p_W_FASER2vmu, detector='FASER2')
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_n_FASERvmu_1p2_3 + N_2mu_n_FASERvmubar_1p2_3, XSecCon_2mu_n_FASERvmu_1p2_3 + XSecCon_2mu_n_FASERvmubar_1p2_3, upper_N_2mu_n_W_FASERvmu_1p2_3 + upper_N_2mu_n_W_FASERvmubar_1p2_3, lower_N_2mu_n_W_FASERvmu_1p2_3 + lower_N_2mu_n_W_FASERvmubar_1p2_3, upper_XSecCon_2mu_n_W_FASERvmu_1p2_3 + upper_XSecCon_2mu_n_W_FASERvmubar_1p2_3, lower_XSecCon_2mu_n_W_FASERvmu_1p2_3 + lower_XSecCon_2mu_n_W_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_n_FASER2vmu + N_2mu_n_FASER2vmubar, XSecCon_2mu_n_FASER2vmu + XSecCon_2mu_n_FASER2vmubar, upper_N_2mu_n_W_FASER2vmu + upper_N_2mu_n_W_FASER2vmubar, lower_N_2mu_n_W_FASER2vmu + lower_N_2mu_n_W_FASER2vmubar, upper_XSecCon_2mu_n_W_FASER2vmu + upper_XSecCon_2mu_n_W_FASER2vmubar, lower_XSecCon_2mu_n_W_FASER2vmu + lower_XSecCon_2mu_n_W_FASER2vmubar, detector='FASER2')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_n_FASER2vmu, XSecCon_2mu_n_FASER2vmu, upper_N_2mu_n_W_FASER2vmu, lower_N_2mu_n_W_FASER2vmu, upper_XSecCon_2mu_n_W_FASER2vmu, lower_XSecCon_2mu_n_W_FASER2vmu, detector='FASER2')
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_incoh_FASERvmu_1p2_3 + N_2mu_incoh_FASERvmubar_1p2_3, XSecCon_2mu_incoh_FASERvmu_1p2_3 + XSecCon_2mu_incoh_FASERvmubar_1p2_3, upper_N_2mu_incoh_FASERvmu_1p2_3 + upper_N_2mu_incoh_FASERvmubar_1p2_3, lower_N_2mu_incoh_FASERvmu_1p2_3 + lower_N_2mu_incoh_FASERvmubar_1p2_3, upper_XSecCon_2mu_incoh_FASERvmu_1p2_3 + upper_XSecCon_2mu_incoh_FASERvmubar_1p2_3, lower_XSecCon_2mu_incoh_FASERvmu_1p2_3 + lower_XSecCon_2mu_incoh_FASERvmubar_1p2_3, detector='FASER')
#        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_incoh_FASER2vmu + N_2mu_incoh_FASER2vmubar, XSecCon_2mu_incoh_FASER2vmu + XSecCon_2mu_incoh_FASER2vmubar, upper_N_2mu_incoh_FASER2vmu + upper_N_2mu_incoh_FASER2vmubar, lower_N_2mu_incoh_FASER2vmu + lower_N_2mu_incoh_FASER2vmubar, upper_XSecCon_2mu_incoh_FASER2vmu + upper_XSecCon_2mu_incoh_FASER2vmubar, lower_XSecCon_2mu_incoh_FASER2vmu + lower_XSecCon_2mu_incoh_FASER2vmubar, detector='FASER2')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_incoh_FASER2vmu, XSecCon_2mu_incoh_FASER2vmu, upper_N_2mu_incoh_FASER2vmu, lower_N_2mu_incoh_FASER2vmu, upper_XSecCon_2mu_incoh_FASER2vmu, lower_XSecCon_2mu_incoh_FASER2vmu, detector='FASER2')
        print("\n", file=textfile)

        print("vmuCC Cross Check", file=textfile)
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_vmuCC_FASERvmu_1p2_3 + N_vmuCC_FASERvmubar_1p2_3, XSecCon_vmuCC_FASERvmu_1p2_3 + XSecCon_vmuCC_FASERvmubar_1p2_3, detector='FASER')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_DIS_vmuCC_FASERvmu + N_DIS_vmubarCC_FASERvmubar, 1, detector='FASER')
        print("\n", file=textfile)

        print("%%%%%%%%%%%% MINOS and MINOS+ -- NEUTRINO MODE %%%%%%%%%%%%", file=textfile)
        print("---------------------- vmu ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_coh_Fe_MINOS_neutrino_vmu, XSecCon_1tau_coh_Fe_MINOS_neutrino_vmu, upper_N_1tau_coh_Fe_MINOS_neutrino_vmu, lower_N_1tau_coh_Fe_MINOS_neutrino_vmu, upper_XSecCon_1tau_coh_Fe_MINOS_neutrino_vmu, lower_XSecCon_1tau_coh_Fe_MINOS_neutrino_vmu, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_1tau_coh_Fe_MINOSPlus_neutrino_vmu, XSecCon_1tau_coh_Fe_MINOSPlus_neutrino_vmu, upper_N_1tau_coh_Fe_MINOSPlus_neutrino_vmu, lower_N_1tau_coh_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_1tau_coh_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_1tau_coh_Fe_MINOSPlus_neutrino_vmu, detector='MINOS+')
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_p_Fe_MINOS_neutrino_vmu, XSecCon_1tau_p_Fe_MINOS_neutrino_vmu, upper_N_1tau_p_Fe_MINOS_neutrino_vmu, lower_N_1tau_p_Fe_MINOS_neutrino_vmu, upper_XSecCon_1tau_p_Fe_MINOS_neutrino_vmu, lower_XSecCon_1tau_p_Fe_MINOS_neutrino_vmu, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_1tau_p_Fe_MINOSPlus_neutrino_vmu, XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmu, upper_N_1tau_p_Fe_MINOSPlus_neutrino_vmu, lower_N_1tau_p_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmu, detector='MINOS+')
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_n_Fe_MINOS_neutrino_vmu, XSecCon_1tau_n_Fe_MINOS_neutrino_vmu, upper_N_1tau_n_Fe_MINOS_neutrino_vmu, lower_N_1tau_n_Fe_MINOS_neutrino_vmu, upper_XSecCon_1tau_n_Fe_MINOS_neutrino_vmu, lower_XSecCon_1tau_n_Fe_MINOS_neutrino_vmu, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_1tau_n_Fe_MINOSPlus_neutrino_vmu, XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmu, upper_N_1tau_n_Fe_MINOSPlus_neutrino_vmu, lower_N_1tau_n_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmu, detector='MINOS+')
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_incoh_Fe_MINOS_neutrino_vmu, XSecCon_1tau_incoh_Fe_MINOS_neutrino_vmu, upper_N_1tau_incoh_Fe_MINOS_neutrino_vmu, lower_N_1tau_incoh_Fe_MINOS_neutrino_vmu, upper_XSecCon_1tau_incoh_Fe_MINOS_neutrino_vmu, lower_XSecCon_1tau_incoh_Fe_MINOS_neutrino_vmu, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_1tau_incoh_Fe_MINOSPlus_neutrino_vmu, XSecCon_1tau_incoh_Fe_MINOSPlus_neutrino_vmu, upper_N_1tau_incoh_Fe_MINOSPlus_neutrino_vmu, lower_N_1tau_incoh_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_1tau_incoh_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_1tau_incoh_Fe_MINOSPlus_neutrino_vmu, detector='MINOS+')
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_coh_Fe_MINOS_neutrino_vmu, XSecCon_2mu_coh_Fe_MINOS_neutrino_vmu, upper_N_2mu_coh_Fe_MINOS_neutrino_vmu, lower_N_2mu_coh_Fe_MINOS_neutrino_vmu, upper_XSecCon_2mu_coh_Fe_MINOS_neutrino_vmu, lower_XSecCon_2mu_coh_Fe_MINOS_neutrino_vmu, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_2mu_coh_Fe_MINOSPlus_neutrino_vmu, XSecCon_2mu_coh_Fe_MINOSPlus_neutrino_vmu, upper_N_2mu_coh_Fe_MINOSPlus_neutrino_vmu, lower_N_2mu_coh_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_2mu_coh_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_2mu_coh_Fe_MINOSPlus_neutrino_vmu, detector='MINOS+')
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_p_Fe_MINOS_neutrino_vmu, XSecCon_2mu_p_Fe_MINOS_neutrino_vmu, upper_N_2mu_p_Fe_MINOS_neutrino_vmu, lower_N_2mu_p_Fe_MINOS_neutrino_vmu, upper_XSecCon_2mu_p_Fe_MINOS_neutrino_vmu, lower_XSecCon_2mu_p_Fe_MINOS_neutrino_vmu, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_2mu_p_Fe_MINOSPlus_neutrino_vmu, XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmu, upper_N_2mu_p_Fe_MINOSPlus_neutrino_vmu, lower_N_2mu_p_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmu, detector='MINOS+')
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_n_Fe_MINOS_neutrino_vmu, XSecCon_2mu_n_Fe_MINOS_neutrino_vmu, upper_N_2mu_n_Fe_MINOS_neutrino_vmu, lower_N_2mu_n_Fe_MINOS_neutrino_vmu, upper_XSecCon_2mu_n_Fe_MINOS_neutrino_vmu, lower_XSecCon_2mu_n_Fe_MINOS_neutrino_vmu, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_2mu_n_Fe_MINOSPlus_neutrino_vmu, XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmu, upper_N_2mu_n_Fe_MINOSPlus_neutrino_vmu, lower_N_2mu_n_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmu, detector='MINOS+')
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_incoh_Fe_MINOS_neutrino_vmu, XSecCon_2mu_incoh_Fe_MINOS_neutrino_vmu, upper_N_2mu_incoh_Fe_MINOS_neutrino_vmu, lower_N_2mu_incoh_Fe_MINOS_neutrino_vmu, upper_XSecCon_2mu_incoh_Fe_MINOS_neutrino_vmu, lower_XSecCon_2mu_incoh_Fe_MINOS_neutrino_vmu, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_2mu_incoh_Fe_MINOSPlus_neutrino_vmu, XSecCon_2mu_incoh_Fe_MINOSPlus_neutrino_vmu, upper_N_2mu_incoh_Fe_MINOSPlus_neutrino_vmu, lower_N_2mu_incoh_Fe_MINOSPlus_neutrino_vmu, upper_XSecCon_2mu_incoh_Fe_MINOSPlus_neutrino_vmu, lower_XSecCon_2mu_incoh_Fe_MINOSPlus_neutrino_vmu, detector='MINOS+')
        print("\n", file=textfile)

        print("---------------------- vmubar ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_coh_Fe_MINOS_neutrino_vmubar, XSecCon_1tau_coh_Fe_MINOS_neutrino_vmubar, upper_N_1tau_coh_Fe_MINOS_neutrino_vmubar, lower_N_1tau_coh_Fe_MINOS_neutrino_vmubar, upper_XSecCon_1tau_coh_Fe_MINOS_neutrino_vmubar, lower_XSecCon_1tau_coh_Fe_MINOS_neutrino_vmubar, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_1tau_coh_Fe_MINOSPlus_neutrino_vmubar, XSecCon_1tau_coh_Fe_MINOSPlus_neutrino_vmubar, upper_N_1tau_coh_Fe_MINOSPlus_neutrino_vmubar, lower_N_1tau_coh_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_1tau_coh_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_1tau_coh_Fe_MINOSPlus_neutrino_vmubar, detector='MINOS+')
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_p_Fe_MINOS_neutrino_vmubar, XSecCon_1tau_p_Fe_MINOS_neutrino_vmubar, upper_N_1tau_p_Fe_MINOS_neutrino_vmubar, lower_N_1tau_p_Fe_MINOS_neutrino_vmubar, upper_XSecCon_1tau_p_Fe_MINOS_neutrino_vmubar, lower_XSecCon_1tau_p_Fe_MINOS_neutrino_vmubar, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_1tau_p_Fe_MINOSPlus_neutrino_vmubar, XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmubar, upper_N_1tau_p_Fe_MINOSPlus_neutrino_vmubar, lower_N_1tau_p_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_1tau_p_Fe_MINOSPlus_neutrino_vmubar, detector='MINOS+')
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_n_Fe_MINOS_neutrino_vmubar, XSecCon_1tau_n_Fe_MINOS_neutrino_vmubar, upper_N_1tau_n_Fe_MINOS_neutrino_vmubar, lower_N_1tau_n_Fe_MINOS_neutrino_vmubar, upper_XSecCon_1tau_n_Fe_MINOS_neutrino_vmubar, lower_XSecCon_1tau_n_Fe_MINOS_neutrino_vmubar, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_1tau_n_Fe_MINOSPlus_neutrino_vmubar, XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmubar, upper_N_1tau_n_Fe_MINOSPlus_neutrino_vmubar, lower_N_1tau_n_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_1tau_n_Fe_MINOSPlus_neutrino_vmubar, detector='MINOS+')
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_incoh_Fe_MINOS_neutrino_vmubar, XSecCon_1tau_incoh_Fe_MINOS_neutrino_vmubar, upper_N_1tau_incoh_Fe_MINOS_neutrino_vmubar, lower_N_1tau_incoh_Fe_MINOS_neutrino_vmubar, upper_XSecCon_1tau_incoh_Fe_MINOS_neutrino_vmubar, lower_XSecCon_1tau_incoh_Fe_MINOS_neutrino_vmubar, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_1tau_incoh_Fe_MINOSPlus_neutrino_vmubar, XSecCon_1tau_incoh_Fe_MINOSPlus_neutrino_vmubar, upper_N_1tau_incoh_Fe_MINOSPlus_neutrino_vmubar, lower_N_1tau_incoh_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_1tau_incoh_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_1tau_incoh_Fe_MINOSPlus_neutrino_vmubar, detector='MINOS+')
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_coh_Fe_MINOS_neutrino_vmubar, XSecCon_2mu_coh_Fe_MINOS_neutrino_vmubar, upper_N_2mu_coh_Fe_MINOS_neutrino_vmubar, lower_N_2mu_coh_Fe_MINOS_neutrino_vmubar, upper_XSecCon_2mu_coh_Fe_MINOS_neutrino_vmubar, lower_XSecCon_2mu_coh_Fe_MINOS_neutrino_vmubar, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_2mu_coh_Fe_MINOSPlus_neutrino_vmubar, XSecCon_2mu_coh_Fe_MINOSPlus_neutrino_vmubar, upper_N_2mu_coh_Fe_MINOSPlus_neutrino_vmubar, lower_N_2mu_coh_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_2mu_coh_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_2mu_coh_Fe_MINOSPlus_neutrino_vmubar, detector='MINOS+')
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_p_Fe_MINOS_neutrino_vmubar, XSecCon_2mu_p_Fe_MINOS_neutrino_vmubar, upper_N_2mu_p_Fe_MINOS_neutrino_vmubar, lower_N_2mu_p_Fe_MINOS_neutrino_vmubar, upper_XSecCon_2mu_p_Fe_MINOS_neutrino_vmubar, lower_XSecCon_2mu_p_Fe_MINOS_neutrino_vmubar, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_2mu_p_Fe_MINOSPlus_neutrino_vmubar, XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmubar, upper_N_2mu_p_Fe_MINOSPlus_neutrino_vmubar, lower_N_2mu_p_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_2mu_p_Fe_MINOSPlus_neutrino_vmubar, detector='MINOS+')
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_n_Fe_MINOS_neutrino_vmubar, XSecCon_2mu_n_Fe_MINOS_neutrino_vmubar, upper_N_2mu_n_Fe_MINOS_neutrino_vmubar, lower_N_2mu_n_Fe_MINOS_neutrino_vmubar, upper_XSecCon_2mu_n_Fe_MINOS_neutrino_vmubar, lower_XSecCon_2mu_n_Fe_MINOS_neutrino_vmubar, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_2mu_n_Fe_MINOSPlus_neutrino_vmubar, XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmubar, upper_N_2mu_n_Fe_MINOSPlus_neutrino_vmubar, lower_N_2mu_n_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_2mu_n_Fe_MINOSPlus_neutrino_vmubar, detector='MINOS+')
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_incoh_Fe_MINOS_neutrino_vmubar, XSecCon_2mu_incoh_Fe_MINOS_neutrino_vmubar, upper_N_2mu_incoh_Fe_MINOS_neutrino_vmubar, lower_N_2mu_incoh_Fe_MINOS_neutrino_vmubar, upper_XSecCon_2mu_incoh_Fe_MINOS_neutrino_vmubar, lower_XSecCon_2mu_incoh_Fe_MINOS_neutrino_vmubar, detector='MINOS_neutrino')
        PrintOutEvent(textfile, "MINOS+", M_MINOS, 1, N_2mu_incoh_Fe_MINOSPlus_neutrino_vmubar, XSecCon_2mu_incoh_Fe_MINOSPlus_neutrino_vmubar, upper_N_2mu_incoh_Fe_MINOSPlus_neutrino_vmubar, lower_N_2mu_incoh_Fe_MINOSPlus_neutrino_vmubar, upper_XSecCon_2mu_incoh_Fe_MINOSPlus_neutrino_vmubar, lower_XSecCon_2mu_incoh_Fe_MINOSPlus_neutrino_vmubar, detector='MINOS+')
        print("\n", file=textfile)

        print("%%%%%%%%%% MINOS and MINOS+ -- ANTINEUTRINO MODE %%%%%%%%%%", file=textfile)
        print("---------------------- vmu ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_coh_Fe_MINOS_antineutrino_vmu, XSecCon_1tau_coh_Fe_MINOS_antineutrino_vmu, upper_N_1tau_coh_Fe_MINOS_antineutrino_vmu, lower_N_1tau_coh_Fe_MINOS_antineutrino_vmu, upper_XSecCon_1tau_coh_Fe_MINOS_antineutrino_vmu, lower_XSecCon_1tau_coh_Fe_MINOS_antineutrino_vmu, detector='MINOS_antineutrino')
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_p_Fe_MINOS_antineutrino_vmu, XSecCon_1tau_p_Fe_MINOS_antineutrino_vmu, upper_N_1tau_p_Fe_MINOS_antineutrino_vmu, lower_N_1tau_p_Fe_MINOS_antineutrino_vmu, upper_XSecCon_1tau_p_Fe_MINOS_antineutrino_vmu, lower_XSecCon_1tau_p_Fe_MINOS_antineutrino_vmu, detector='MINOS_antineutrino')
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_n_Fe_MINOS_antineutrino_vmu, XSecCon_1tau_n_Fe_MINOS_antineutrino_vmu, upper_N_1tau_n_Fe_MINOS_antineutrino_vmu, lower_N_1tau_n_Fe_MINOS_antineutrino_vmu, upper_XSecCon_1tau_n_Fe_MINOS_antineutrino_vmu, lower_XSecCon_1tau_n_Fe_MINOS_antineutrino_vmu, detector='MINOS_antineutrino')
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_incoh_Fe_MINOS_antineutrino_vmu, XSecCon_1tau_incoh_Fe_MINOS_antineutrino_vmu, upper_N_1tau_incoh_Fe_MINOS_antineutrino_vmu, lower_N_1tau_incoh_Fe_MINOS_antineutrino_vmu, upper_XSecCon_1tau_incoh_Fe_MINOS_antineutrino_vmu, lower_XSecCon_1tau_incoh_Fe_MINOS_antineutrino_vmu, detector='MINOS_antineutrino')
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_coh_Fe_MINOS_antineutrino_vmu, XSecCon_2mu_coh_Fe_MINOS_antineutrino_vmu, upper_N_2mu_coh_Fe_MINOS_antineutrino_vmu, lower_N_2mu_coh_Fe_MINOS_antineutrino_vmu, upper_XSecCon_2mu_coh_Fe_MINOS_antineutrino_vmu, lower_XSecCon_2mu_coh_Fe_MINOS_antineutrino_vmu, detector='MINOS_antineutrino')
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_p_Fe_MINOS_antineutrino_vmu, XSecCon_2mu_p_Fe_MINOS_antineutrino_vmu, upper_N_2mu_p_Fe_MINOS_antineutrino_vmu, lower_N_2mu_p_Fe_MINOS_antineutrino_vmu, upper_XSecCon_2mu_p_Fe_MINOS_antineutrino_vmu, lower_XSecCon_2mu_p_Fe_MINOS_antineutrino_vmu, detector='MINOS_antineutrino')
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_n_Fe_MINOS_antineutrino_vmu, XSecCon_2mu_n_Fe_MINOS_antineutrino_vmu, upper_N_2mu_n_Fe_MINOS_antineutrino_vmu, lower_N_2mu_n_Fe_MINOS_antineutrino_vmu, upper_XSecCon_2mu_n_Fe_MINOS_antineutrino_vmu, lower_XSecCon_2mu_n_Fe_MINOS_antineutrino_vmu, detector='MINOS_antineutrino')
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_incoh_Fe_MINOS_antineutrino_vmu, XSecCon_2mu_incoh_Fe_MINOS_antineutrino_vmu, upper_N_2mu_incoh_Fe_MINOS_antineutrino_vmu, lower_N_2mu_incoh_Fe_MINOS_antineutrino_vmu, upper_XSecCon_2mu_incoh_Fe_MINOS_antineutrino_vmu, lower_XSecCon_2mu_incoh_Fe_MINOS_antineutrino_vmu, detector='MINOS_antineutrino')
        print("\n", file=textfile)

        print("---------------------- vmubar ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_coh_Fe_MINOS_antineutrino_vmubar, XSecCon_1tau_coh_Fe_MINOS_antineutrino_vmubar, upper_N_1tau_coh_Fe_MINOS_antineutrino_vmubar, lower_N_1tau_coh_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_1tau_coh_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_1tau_coh_Fe_MINOS_antineutrino_vmubar, detector='MINOS_antineutrino')
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_p_Fe_MINOS_antineutrino_vmubar, XSecCon_1tau_p_Fe_MINOS_antineutrino_vmubar, upper_N_1tau_p_Fe_MINOS_antineutrino_vmubar, lower_N_1tau_p_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_1tau_p_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_1tau_p_Fe_MINOS_antineutrino_vmubar, detector='MINOS_antineutrino')
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_n_Fe_MINOS_antineutrino_vmubar, XSecCon_1tau_n_Fe_MINOS_antineutrino_vmubar, upper_N_1tau_n_Fe_MINOS_antineutrino_vmubar, lower_N_1tau_n_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_1tau_n_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_1tau_n_Fe_MINOS_antineutrino_vmubar, detector='MINOS_antineutrino')
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_1tau_incoh_Fe_MINOS_antineutrino_vmubar, XSecCon_1tau_incoh_Fe_MINOS_antineutrino_vmubar, upper_N_1tau_incoh_Fe_MINOS_antineutrino_vmubar, lower_N_1tau_incoh_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_1tau_incoh_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_1tau_incoh_Fe_MINOS_antineutrino_vmubar, detector='MINOS_antineutrino')
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_coh_Fe_MINOS_antineutrino_vmubar, XSecCon_2mu_coh_Fe_MINOS_antineutrino_vmubar, upper_N_2mu_coh_Fe_MINOS_antineutrino_vmubar, lower_N_2mu_coh_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_2mu_coh_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_2mu_coh_Fe_MINOS_antineutrino_vmubar, detector='MINOS_antineutrino')
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_p_Fe_MINOS_antineutrino_vmubar, XSecCon_2mu_p_Fe_MINOS_antineutrino_vmubar, upper_N_2mu_p_Fe_MINOS_antineutrino_vmubar, lower_N_2mu_p_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_2mu_p_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_2mu_p_Fe_MINOS_antineutrino_vmubar, detector='MINOS_antineutrino')
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_n_Fe_MINOS_antineutrino_vmubar, XSecCon_2mu_n_Fe_MINOS_antineutrino_vmubar, upper_N_2mu_n_Fe_MINOS_antineutrino_vmubar, lower_N_2mu_n_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_2mu_n_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_2mu_n_Fe_MINOS_antineutrino_vmubar, detector='MINOS_antineutrino')
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "MINOS", M_MINOS, 1, N_2mu_incoh_Fe_MINOS_antineutrino_vmubar, XSecCon_2mu_incoh_Fe_MINOS_antineutrino_vmubar, upper_N_2mu_incoh_Fe_MINOS_antineutrino_vmubar, lower_N_2mu_incoh_Fe_MINOS_antineutrino_vmubar, upper_XSecCon_2mu_incoh_Fe_MINOS_antineutrino_vmubar, lower_XSecCon_2mu_incoh_Fe_MINOS_antineutrino_vmubar, detector='MINOS_antineutrino')
        print("\n", file=textfile)

        print("%%%%%%%%%%%% T2K -- INGRID %%%%%%%%%%%%", file=textfile)
        print("---------------------- vmu ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_1tau_coh_Fe_INGRID_neutrino_vmu, XSecCon_1tau_coh_Fe_INGRID_neutrino_vmu, upper_N_1tau_coh_Fe_INGRID_neutrino_vmu, lower_N_1tau_coh_Fe_INGRID_neutrino_vmu, upper_XSecCon_1tau_coh_Fe_INGRID_neutrino_vmu, lower_XSecCon_1tau_coh_Fe_INGRID_neutrino_vmu, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_1tau_coh_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_1tau_coh_Fe_INGRIDPhase2_neutrino_vmu, upper_N_1tau_coh_Fe_INGRIDPhase2_neutrino_vmu, lower_N_1tau_coh_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_1tau_coh_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_1tau_coh_Fe_INGRIDPhase2_neutrino_vmu, detector='INGRID2')
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_1tau_p_Fe_INGRID_neutrino_vmu, XSecCon_1tau_p_Fe_INGRID_neutrino_vmu, upper_N_1tau_p_Fe_INGRID_neutrino_vmu, lower_N_1tau_p_Fe_INGRID_neutrino_vmu, upper_XSecCon_1tau_p_Fe_INGRID_neutrino_vmu, lower_XSecCon_1tau_p_Fe_INGRID_neutrino_vmu, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_1tau_p_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmu, upper_N_1tau_p_Fe_INGRIDPhase2_neutrino_vmu, lower_N_1tau_p_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmu, detector='INGRID2')
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_1tau_n_Fe_INGRID_neutrino_vmu, XSecCon_1tau_n_Fe_INGRID_neutrino_vmu, upper_N_1tau_n_Fe_INGRID_neutrino_vmu, lower_N_1tau_n_Fe_INGRID_neutrino_vmu, upper_XSecCon_1tau_n_Fe_INGRID_neutrino_vmu, lower_XSecCon_1tau_n_Fe_INGRID_neutrino_vmu, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_1tau_n_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmu, upper_N_1tau_n_Fe_INGRIDPhase2_neutrino_vmu, lower_N_1tau_n_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmu, detector='INGRID2')
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_1tau_incoh_Fe_INGRID_neutrino_vmu, XSecCon_1tau_incoh_Fe_INGRID_neutrino_vmu, upper_N_1tau_incoh_Fe_INGRID_neutrino_vmu, lower_N_1tau_incoh_Fe_INGRID_neutrino_vmu, upper_XSecCon_1tau_incoh_Fe_INGRID_neutrino_vmu, lower_XSecCon_1tau_incoh_Fe_INGRID_neutrino_vmu, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmu, upper_N_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmu, lower_N_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmu, detector='INGRID2')
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_2mu_coh_Fe_INGRID_neutrino_vmu, XSecCon_2mu_coh_Fe_INGRID_neutrino_vmu, upper_N_2mu_coh_Fe_INGRID_neutrino_vmu, lower_N_2mu_coh_Fe_INGRID_neutrino_vmu, upper_XSecCon_2mu_coh_Fe_INGRID_neutrino_vmu, lower_XSecCon_2mu_coh_Fe_INGRID_neutrino_vmu, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_2mu_coh_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_2mu_coh_Fe_INGRIDPhase2_neutrino_vmu, upper_N_2mu_coh_Fe_INGRIDPhase2_neutrino_vmu, lower_N_2mu_coh_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_2mu_coh_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_2mu_coh_Fe_INGRIDPhase2_neutrino_vmu, detector='INGRID2')
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_2mu_p_Fe_INGRID_neutrino_vmu, XSecCon_2mu_p_Fe_INGRID_neutrino_vmu, upper_N_2mu_p_Fe_INGRID_neutrino_vmu, lower_N_2mu_p_Fe_INGRID_neutrino_vmu, upper_XSecCon_2mu_p_Fe_INGRID_neutrino_vmu, lower_XSecCon_2mu_p_Fe_INGRID_neutrino_vmu, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_2mu_p_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmu, upper_N_2mu_p_Fe_INGRIDPhase2_neutrino_vmu, lower_N_2mu_p_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmu, detector='INGRID2')
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_2mu_n_Fe_INGRID_neutrino_vmu, XSecCon_2mu_n_Fe_INGRID_neutrino_vmu, upper_N_2mu_n_Fe_INGRID_neutrino_vmu, lower_N_2mu_n_Fe_INGRID_neutrino_vmu, upper_XSecCon_2mu_n_Fe_INGRID_neutrino_vmu, lower_XSecCon_2mu_n_Fe_INGRID_neutrino_vmu, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_2mu_n_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmu, upper_N_2mu_n_Fe_INGRIDPhase2_neutrino_vmu, lower_N_2mu_n_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmu, detector='INGRID2')
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_2mu_incoh_Fe_INGRID_neutrino_vmu, XSecCon_2mu_incoh_Fe_INGRID_neutrino_vmu, upper_N_2mu_incoh_Fe_INGRID_neutrino_vmu, lower_N_2mu_incoh_Fe_INGRID_neutrino_vmu, upper_XSecCon_2mu_incoh_Fe_INGRID_neutrino_vmu, lower_XSecCon_2mu_incoh_Fe_INGRID_neutrino_vmu, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmu, XSecCon_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmu, upper_N_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmu, lower_N_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmu, upper_XSecCon_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmu, lower_XSecCon_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmu, detector='INGRID2')
        print("\n", file=textfile)

        print("---------------------- vmubar ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_1tau_coh_Fe_INGRID_neutrino_vmubar, XSecCon_1tau_coh_Fe_INGRID_neutrino_vmubar, upper_N_1tau_coh_Fe_INGRID_neutrino_vmubar, lower_N_1tau_coh_Fe_INGRID_neutrino_vmubar, upper_XSecCon_1tau_coh_Fe_INGRID_neutrino_vmubar, lower_XSecCon_1tau_coh_Fe_INGRID_neutrino_vmubar, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_1tau_coh_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_1tau_coh_Fe_INGRIDPhase2_neutrino_vmubar, upper_N_1tau_coh_Fe_INGRIDPhase2_neutrino_vmubar, lower_N_1tau_coh_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_1tau_coh_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_1tau_coh_Fe_INGRIDPhase2_neutrino_vmubar, detector='INGRID2')
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_1tau_p_Fe_INGRID_neutrino_vmubar, XSecCon_1tau_p_Fe_INGRID_neutrino_vmubar, upper_N_1tau_p_Fe_INGRID_neutrino_vmubar, lower_N_1tau_p_Fe_INGRID_neutrino_vmubar, upper_XSecCon_1tau_p_Fe_INGRID_neutrino_vmubar, lower_XSecCon_1tau_p_Fe_INGRID_neutrino_vmubar, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar, upper_N_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar, lower_N_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_1tau_p_Fe_INGRIDPhase2_neutrino_vmubar, detector='INGRID2')
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_1tau_n_Fe_INGRID_neutrino_vmubar, XSecCon_1tau_n_Fe_INGRID_neutrino_vmubar, upper_N_1tau_n_Fe_INGRID_neutrino_vmubar, lower_N_1tau_n_Fe_INGRID_neutrino_vmubar, upper_XSecCon_1tau_n_Fe_INGRID_neutrino_vmubar, lower_XSecCon_1tau_n_Fe_INGRID_neutrino_vmubar, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar, upper_N_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar, lower_N_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_1tau_n_Fe_INGRIDPhase2_neutrino_vmubar, detector='INGRID2')
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_1tau_incoh_Fe_INGRID_neutrino_vmubar, XSecCon_1tau_incoh_Fe_INGRID_neutrino_vmubar, upper_N_1tau_incoh_Fe_INGRID_neutrino_vmubar, lower_N_1tau_incoh_Fe_INGRID_neutrino_vmubar, upper_XSecCon_1tau_incoh_Fe_INGRID_neutrino_vmubar, lower_XSecCon_1tau_incoh_Fe_INGRID_neutrino_vmubar, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmubar, upper_N_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmubar, lower_N_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_1tau_incoh_Fe_INGRIDPhase2_neutrino_vmubar, detector='INGRID2')
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_2mu_coh_Fe_INGRID_neutrino_vmubar, XSecCon_2mu_coh_Fe_INGRID_neutrino_vmubar, upper_N_2mu_coh_Fe_INGRID_neutrino_vmubar, lower_N_2mu_coh_Fe_INGRID_neutrino_vmubar, upper_XSecCon_2mu_coh_Fe_INGRID_neutrino_vmubar, lower_XSecCon_2mu_coh_Fe_INGRID_neutrino_vmubar, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_2mu_coh_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_2mu_coh_Fe_INGRIDPhase2_neutrino_vmubar, upper_N_2mu_coh_Fe_INGRIDPhase2_neutrino_vmubar, lower_N_2mu_coh_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_2mu_coh_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_2mu_coh_Fe_INGRIDPhase2_neutrino_vmubar, detector='INGRID2')
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_2mu_p_Fe_INGRID_neutrino_vmubar, XSecCon_2mu_p_Fe_INGRID_neutrino_vmubar, upper_N_2mu_p_Fe_INGRID_neutrino_vmubar, lower_N_2mu_p_Fe_INGRID_neutrino_vmubar, upper_XSecCon_2mu_p_Fe_INGRID_neutrino_vmubar, lower_XSecCon_2mu_p_Fe_INGRID_neutrino_vmubar, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar, upper_N_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar, lower_N_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_2mu_p_Fe_INGRIDPhase2_neutrino_vmubar, detector='INGRID2')
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_2mu_n_Fe_INGRID_neutrino_vmubar, XSecCon_2mu_n_Fe_INGRID_neutrino_vmubar, upper_N_2mu_n_Fe_INGRID_neutrino_vmubar, lower_N_2mu_n_Fe_INGRID_neutrino_vmubar, upper_XSecCon_2mu_n_Fe_INGRID_neutrino_vmubar, lower_XSecCon_2mu_n_Fe_INGRID_neutrino_vmubar, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar, upper_N_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar, lower_N_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_2mu_n_Fe_INGRIDPhase2_neutrino_vmubar, detector='INGRID2')
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "INGRID", M_INGRID, 1, N_2mu_incoh_Fe_INGRID_neutrino_vmubar, XSecCon_2mu_incoh_Fe_INGRID_neutrino_vmubar, upper_N_2mu_incoh_Fe_INGRID_neutrino_vmubar, lower_N_2mu_incoh_Fe_INGRID_neutrino_vmubar, upper_XSecCon_2mu_incoh_Fe_INGRID_neutrino_vmubar, lower_XSecCon_2mu_incoh_Fe_INGRID_neutrino_vmubar, detector='INGRID')
        PrintOutEvent(textfile, "INGRID Phase 2", M_INGRID, 1, N_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmubar, XSecCon_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmubar, upper_N_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmubar, lower_N_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmubar, upper_XSecCon_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmubar, lower_XSecCon_2mu_incoh_Fe_INGRIDPhase2_neutrino_vmubar, detector='INGRID2')
        print("\n", file=textfile)

        print("%%%%%%%%%%%%%%%%%% SHiP %%%%%%%%%%%%%%%%%%", file=textfile)
        print("---------------------- vmu ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_1tau_coh_W_SHiP_vmu, XSecCon_1tau_coh_W_SHiP_vmu, upper_N_1tau_coh_W_SHiP_vmu, lower_N_1tau_coh_W_SHiP_vmu, upper_XSecCon_1tau_coh_W_SHiP_vmu, lower_XSecCon_1tau_coh_W_SHiP_vmu, detector='SHiP')
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_1tau_p_W_SHiP_vmu, XSecCon_1tau_p_W_SHiP_vmu, upper_N_1tau_p_W_SHiP_vmu, lower_N_1tau_p_W_SHiP_vmu, upper_XSecCon_1tau_p_W_SHiP_vmu, lower_XSecCon_1tau_p_W_SHiP_vmu, detector='SHiP')
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_1tau_n_W_SHiP_vmu, XSecCon_1tau_n_W_SHiP_vmu, upper_N_1tau_n_W_SHiP_vmu, lower_N_1tau_n_W_SHiP_vmu, upper_XSecCon_1tau_n_W_SHiP_vmu, lower_XSecCon_1tau_n_W_SHiP_vmu, detector='SHiP')
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_1tau_incoh_SHiP_vmu, XSecCon_1tau_incoh_SHiP_vmu, upper_N_1tau_incoh_SHiP_vmu, lower_N_1tau_incoh_SHiP_vmu, upper_XSecCon_1tau_incoh_SHiP_vmu, lower_XSecCon_1tau_incoh_SHiP_vmu, detector='SHiP')
        print("\n", file=textfile)

        print("2tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2tau_coh_W_SHiP_vmu, XSecCon_2tau_coh_W_SHiP_vmu, upper_N_2tau_coh_W_SHiP_vmu, lower_N_2tau_coh_W_SHiP_vmu, upper_XSecCon_2tau_coh_W_SHiP_vmu, lower_XSecCon_2tau_coh_W_SHiP_vmu, detector='SHiP')
        print("2tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2tau_p_W_SHiP_vmu, XSecCon_2tau_p_W_SHiP_vmu, upper_N_2tau_p_W_SHiP_vmu, lower_N_2tau_p_W_SHiP_vmu, upper_XSecCon_2tau_p_W_SHiP_vmu, lower_XSecCon_2tau_p_W_SHiP_vmu, detector='SHiP')
        print("2tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2tau_n_W_SHiP_vmu, XSecCon_2tau_n_W_SHiP_vmu, upper_N_2tau_n_W_SHiP_vmu, lower_N_2tau_n_W_SHiP_vmu, upper_XSecCon_2tau_n_W_SHiP_vmu, lower_XSecCon_2tau_n_W_SHiP_vmu, detector='SHiP')
        print("2tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2tau_incoh_SHiP_vmu, XSecCon_2tau_incoh_SHiP_vmu, upper_N_2tau_incoh_SHiP_vmu, lower_N_2tau_incoh_SHiP_vmu, upper_XSecCon_2tau_incoh_SHiP_vmu, lower_XSecCon_2tau_incoh_SHiP_vmu, detector='SHiP')
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2mu_coh_W_SHiP_vmu, XSecCon_2mu_coh_W_SHiP_vmu, upper_N_2mu_coh_W_SHiP_vmu, lower_N_2mu_coh_W_SHiP_vmu, upper_XSecCon_2mu_coh_W_SHiP_vmu, lower_XSecCon_2mu_coh_W_SHiP_vmu, detector='SHiP')
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2mu_p_W_SHiP_vmu, XSecCon_2mu_p_W_SHiP_vmu, upper_N_2mu_p_W_SHiP_vmu, lower_N_2mu_p_W_SHiP_vmu, upper_XSecCon_2mu_p_W_SHiP_vmu, lower_XSecCon_2mu_p_W_SHiP_vmu, detector='SHiP')
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2mu_n_W_SHiP_vmu, XSecCon_2mu_n_W_SHiP_vmu, upper_N_2mu_n_W_SHiP_vmu, lower_N_2mu_n_W_SHiP_vmu, upper_XSecCon_2mu_n_W_SHiP_vmu, lower_XSecCon_2mu_n_W_SHiP_vmu, detector='SHiP')
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2mu_incoh_SHiP_vmu, XSecCon_2mu_incoh_SHiP_vmu, upper_N_2mu_incoh_SHiP_vmu, lower_N_2mu_incoh_SHiP_vmu, upper_XSecCon_2mu_incoh_SHiP_vmu, lower_XSecCon_2mu_incoh_SHiP_vmu, detector='SHiP')
        print("\n", file=textfile)

        print("---------------------- vmubar ----------------------", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_1tau_coh_W_SHiP_vmubar, XSecCon_1tau_coh_W_SHiP_vmubar, upper_N_1tau_coh_W_SHiP_vmubar, lower_N_1tau_coh_W_SHiP_vmubar, upper_XSecCon_1tau_coh_W_SHiP_vmubar, lower_XSecCon_1tau_coh_W_SHiP_vmubar, detector='SHiP')
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_1tau_p_W_SHiP_vmubar, XSecCon_1tau_p_W_SHiP_vmubar, upper_N_1tau_p_W_SHiP_vmubar, lower_N_1tau_p_W_SHiP_vmubar, upper_XSecCon_1tau_p_W_SHiP_vmubar, lower_XSecCon_1tau_p_W_SHiP_vmubar, detector='SHiP')
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_1tau_n_W_SHiP_vmubar, XSecCon_1tau_n_W_SHiP_vmubar, upper_N_1tau_n_W_SHiP_vmubar, lower_N_1tau_n_W_SHiP_vmubar, upper_XSecCon_1tau_n_W_SHiP_vmubar, lower_XSecCon_1tau_n_W_SHiP_vmubar, detector='SHiP')
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_1tau_incoh_SHiP_vmubar, XSecCon_1tau_incoh_SHiP_vmubar, upper_N_1tau_incoh_SHiP_vmubar, lower_N_1tau_incoh_SHiP_vmubar, upper_XSecCon_1tau_incoh_SHiP_vmubar, lower_XSecCon_1tau_incoh_SHiP_vmubar, detector='SHiP')
        print("\n", file=textfile)

        print("2tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2tau_coh_W_SHiP_vmubar, XSecCon_2tau_coh_W_SHiP_vmubar, upper_N_2tau_coh_W_SHiP_vmubar, lower_N_2tau_coh_W_SHiP_vmubar, upper_XSecCon_2tau_coh_W_SHiP_vmubar, lower_XSecCon_2tau_coh_W_SHiP_vmubar, detector='SHiP')
        print("2tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2tau_p_W_SHiP_vmubar, XSecCon_2tau_p_W_SHiP_vmubar, upper_N_2tau_p_W_SHiP_vmubar, lower_N_2tau_p_W_SHiP_vmubar, upper_XSecCon_2tau_p_W_SHiP_vmubar, lower_XSecCon_2tau_p_W_SHiP_vmubar, detector='SHiP')
        print("2tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2tau_n_W_SHiP_vmubar, XSecCon_2tau_n_W_SHiP_vmubar, upper_N_2tau_n_W_SHiP_vmubar, lower_N_2tau_n_W_SHiP_vmubar, upper_XSecCon_2tau_n_W_SHiP_vmubar, lower_XSecCon_2tau_n_W_SHiP_vmubar, detector='SHiP')
        print("2tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2tau_incoh_SHiP_vmubar, XSecCon_2tau_incoh_SHiP_vmubar, upper_N_2tau_incoh_SHiP_vmubar, lower_N_2tau_incoh_SHiP_vmubar, upper_XSecCon_2tau_incoh_SHiP_vmubar, lower_XSecCon_2tau_incoh_SHiP_vmubar, detector='SHiP')
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2mu_coh_W_SHiP_vmubar, XSecCon_2mu_coh_W_SHiP_vmubar, upper_N_2mu_coh_W_SHiP_vmubar, lower_N_2mu_coh_W_SHiP_vmubar, upper_XSecCon_2mu_coh_W_SHiP_vmubar, lower_XSecCon_2mu_coh_W_SHiP_vmubar, detector='SHiP')
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2mu_p_W_SHiP_vmubar, XSecCon_2mu_p_W_SHiP_vmubar, upper_N_2mu_p_W_SHiP_vmubar, lower_N_2mu_p_W_SHiP_vmubar, upper_XSecCon_2mu_p_W_SHiP_vmubar, lower_XSecCon_2mu_p_W_SHiP_vmubar, detector='SHiP')
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2mu_n_W_SHiP_vmubar, XSecCon_2mu_n_W_SHiP_vmubar, upper_N_2mu_n_W_SHiP_vmubar, lower_N_2mu_n_W_SHiP_vmubar, upper_XSecCon_2mu_n_W_SHiP_vmubar, lower_XSecCon_2mu_n_W_SHiP_vmubar, detector='SHiP')
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "SHiP", 67, 3, N_2mu_incoh_SHiP_vmubar, XSecCon_2mu_incoh_SHiP_vmubar, upper_N_2mu_incoh_SHiP_vmubar, lower_N_2mu_incoh_SHiP_vmubar, upper_XSecCon_2mu_incoh_SHiP_vmubar, lower_XSecCon_2mu_incoh_SHiP_vmubar, detector='SHiP')
        print("\n", file=textfile)

        print("%%%%%%%%%%%%%%%%%% SBND %%%%%%%%%%%%%%%%%%", file=textfile)
        print("---------------------- vmu ----------------------", file=textfile)
        print("1e1mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "SBND", 67, 3, N_1e1mu_coh_Ar_SBND_vmu, XSecCon_1e1mu_coh_Ar_SBND_vmu, upper_N_1e1mu_coh_Ar_SBND_vmu, lower_N_1e1mu_coh_Ar_SBND_vmu, upper_XSecCon_1e1mu_coh_Ar_SBND_vmu, lower_XSecCon_1e1mu_coh_Ar_SBND_vmu, detector='SBND')
        print("1e1mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "SBND", 67, 3, N_1e1mu_p_Ar_SBND_vmu, XSecCon_1e1mu_p_Ar_SBND_vmu, upper_N_1e1mu_p_Ar_SBND_vmu, lower_N_1e1mu_p_Ar_SBND_vmu, upper_XSecCon_1e1mu_p_Ar_SBND_vmu, lower_XSecCon_1e1mu_p_Ar_SBND_vmu, detector='SBND')
        print("1e1mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "SBND", 67, 3, N_1e1mu_n_Ar_SBND_vmu, XSecCon_1e1mu_n_Ar_SBND_vmu, upper_N_1e1mu_n_Ar_SBND_vmu, lower_N_1e1mu_n_Ar_SBND_vmu, upper_XSecCon_1e1mu_n_Ar_SBND_vmu, lower_XSecCon_1e1mu_n_Ar_SBND_vmu, detector='SBND')
        print("1e1mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "SBND", 67, 3, N_1e1mu_incoh_SBND_vmu, XSecCon_1e1mu_incoh_SBND_vmu, upper_N_1e1mu_incoh_SBND_vmu, lower_N_1e1mu_incoh_SBND_vmu, upper_XSecCon_1e1mu_incoh_SBND_vmu, lower_XSecCon_1e1mu_incoh_SBND_vmu, detector='SBND')
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "SBND", 67, 3, N_2mu_coh_Ar_SBND_vmu, XSecCon_2mu_coh_Ar_SBND_vmu, upper_N_2mu_coh_Ar_SBND_vmu, lower_N_2mu_coh_Ar_SBND_vmu, upper_XSecCon_2mu_coh_Ar_SBND_vmu, lower_XSecCon_2mu_coh_Ar_SBND_vmu, detector='SBND')
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "SBND", 67, 3, N_2mu_p_Ar_SBND_vmu, XSecCon_2mu_p_Ar_SBND_vmu, upper_N_2mu_p_Ar_SBND_vmu, lower_N_2mu_p_Ar_SBND_vmu, upper_XSecCon_2mu_p_Ar_SBND_vmu, lower_XSecCon_2mu_p_Ar_SBND_vmu, detector='SBND')
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "SBND", 67, 3, N_2mu_n_Ar_SBND_vmu, XSecCon_2mu_n_Ar_SBND_vmu, upper_N_2mu_n_Ar_SBND_vmu, lower_N_2mu_n_Ar_SBND_vmu, upper_XSecCon_2mu_n_Ar_SBND_vmu, lower_XSecCon_2mu_n_Ar_SBND_vmu, detector='SBND')
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "SBND", 67, 3, N_2mu_incoh_SBND_vmu, XSecCon_2mu_incoh_SBND_vmu, upper_N_2mu_incoh_SBND_vmu, lower_N_2mu_incoh_SBND_vmu, upper_XSecCon_2mu_incoh_SBND_vmu, lower_XSecCon_2mu_incoh_SBND_vmu, detector='SBND')
        print("\n", file=textfile)

#        print("---------------------- vmubar ----------------------", file=textfile)
#        print("2mu Coherent:", file=textfile)
#        PrintOutEvent(textfile, "SHiP", 67, 3, N_2mu_coh_W_SHiP_vmubar, XSecCon_2mu_coh_W_SHiP_vmubar, upper_N_2mu_coh_W_SHiP_vmubar, lower_N_2mu_coh_W_SHiP_vmubar, upper_XSecCon_2mu_coh_W_SHiP_vmubar, lower_XSecCon_2mu_coh_W_SHiP_vmubar, detector='SHiP')
#        print("2mu Incoherent ; proton:", file=textfile)
#        PrintOutEvent(textfile, "SHiP", 67, 3, N_2mu_p_W_SHiP_vmubar, XSecCon_2mu_p_W_SHiP_vmubar, upper_N_2mu_p_W_SHiP_vmubar, lower_N_2mu_p_W_SHiP_vmubar, upper_XSecCon_2mu_p_W_SHiP_vmubar, lower_XSecCon_2mu_p_W_SHiP_vmubar, detector='SHiP')
#        print("2mu Incoherent ; neutron:", file=textfile)
#        PrintOutEvent(textfile, "SHiP", 67, 3, N_2mu_n_W_SHiP_vmubar, XSecCon_2mu_n_W_SHiP_vmubar, upper_N_2mu_n_W_SHiP_vmubar, lower_N_2mu_n_W_SHiP_vmubar, upper_XSecCon_2mu_n_W_SHiP_vmubar, lower_XSecCon_2mu_n_W_SHiP_vmubar, detector='SHiP')
#        print("2mu Incoherent ; proton + neutron:", file=textfile)
#        PrintOutEvent(textfile, "SHiP", 67, 3, N_2mu_incoh_SHiP_vmubar, XSecCon_2mu_incoh_SHiP_vmubar, upper_N_2mu_incoh_SHiP_vmubar, lower_N_2mu_incoh_SHiP_vmubar, upper_XSecCon_2mu_incoh_SHiP_vmubar, lower_XSecCon_2mu_incoh_SHiP_vmubar, detector='SHiP')
#        print("\n", file=textfile)


WriteOutFile("number_of_events.txt")

### Plotting ###

fig1, ax1 = plt.subplots(1, 1, figsize=(15, 10), tight_layout=True) # Neutrino Mode -- Coherent
fig2, ax2 = plt.subplots(1, 1, figsize=(15, 10), tight_layout=True) # Neutrino Mode -- Incoherent
fig3, ax3 = plt.subplots(1, 1, figsize=(15, 10), tight_layout=True) # Antineutrino Mode -- Coherent
fig4, ax4 = plt.subplots(1, 1, figsize=(15, 10), tight_layout=True) # Antineutrino Mode -- Incoherent

fig5, ax5 = plt.subplots(1, 1, figsize=(15, 10), tight_layout=True) # EPA Neutrino Mode -- Coherent

# with correction factors #
fig6, ax6 = plt.subplots(1, 1, figsize=(15, 10), tight_layout=True) # Neutrino Mode -- Coherent
fig7, ax7 = plt.subplots(1, 1, figsize=(15, 10), tight_layout=True) # Neutrino Mode -- Incoherent
fig8, ax8 = plt.subplots(1, 1, figsize=(15, 10), tight_layout=True) # Antineutrino Mode -- Coherent
fig9, ax9 = plt.subplots(1, 1, figsize=(15, 10), tight_layout=True) # Antineutrino Mode -- Incoherent

mean_1tau_coh_DUNE_neutrino = np.sum(np.multiply(energy_DUNE, Nbin_1tau_coh_Ar_DUNE_neutrino_vmu_67_3)) / np.sum(Nbin_1tau_coh_Ar_DUNE_neutrino_vmu_67_3)
mean_1tau_incoh_DUNE_neutrino = np.sum(np.multiply(energy_DUNE, Nbin_1tau_incoh_DUNE_neutrino_vmu_67_3)) / np.sum(Nbin_1tau_incoh_DUNE_neutrino_vmu_67_3)
mean_1tau_coh_DUNE_antineutrino = np.sum(np.multiply(energy_DUNE, Nbin_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3)) / np.sum(Nbin_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3)
mean_1tau_incoh_DUNE_antineutrino = np.sum(np.multiply(energy_DUNE, Nbin_1tau_incoh_DUNE_antineutrino_vmu_67_3)) / np.sum(Nbin_1tau_incoh_DUNE_antineutrino_vmu_67_3)

mean_1tau_coh_DUNE_tau_opt_neutrino = np.sum(np.multiply(energy_DUNE_tau_opt, Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3)) / np.sum(Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3)
mean_1tau_incoh_DUNE_tau_opt_neutrino = np.sum(np.multiply(energy_DUNE_tau_opt, Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3)) / np.sum(Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3)
mean_1tau_coh_DUNE_tau_opt_antineutrino = np.sum(np.multiply(energy_DUNE_tau_opt, Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)) / np.sum(Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3)
mean_1tau_incoh_DUNE_tau_opt_antineutrino = np.sum(np.multiply(energy_DUNE_tau_opt, Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3)) / np.sum(Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3)

mean_EPA_Heaviside_1tau_coh_DUNE_neutrino = np.sum(np.multiply(energy_DUNE, Nbin_EPA_Heaviside_1tau_coh_Ar_DUNE_neutrino_vmu_67_3)) / np.sum(Nbin_EPA_Heaviside_1tau_coh_Ar_DUNE_neutrino_vmu_67_3)

mean_EPA_Heaviside_1tau_coh_DUNE_tau_opt_neutrino = np.sum(np.multiply(energy_DUNE_tau_opt, Nbin_EPA_Heaviside_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3)) / np.sum(Nbin_EPA_Heaviside_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3)


mean_1tau_coh_DUNE_neutrino_correction_factor_RMiss = np.sum(np.multiply(energy_DUNE, Nbin_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss)) / np.sum(Nbin_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss)
mean_1tau_incoh_DUNE_neutrino_correction_factor_RMiss = np.sum(np.multiply(energy_DUNE, Nbin_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss)) / np.sum(Nbin_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss)
mean_1tau_coh_DUNE_antineutrino_correction_factor_RMiss = np.sum(np.multiply(energy_DUNE, Nbin_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss)) / np.sum(Nbin_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss)
mean_1tau_incoh_DUNE_antineutrino_correction_factor_RMiss = np.sum(np.multiply(energy_DUNE, Nbin_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss)) / np.sum(Nbin_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss)

mean_1tau_coh_DUNE_tau_opt_neutrino_correction_factor_RMiss = np.sum(np.multiply(energy_DUNE_tau_opt, Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss)) / np.sum(Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss)
mean_1tau_incoh_DUNE_tau_opt_neutrino_correction_factor_RMiss = np.sum(np.multiply(energy_DUNE_tau_opt, Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss)) / np.sum(Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss)
mean_1tau_coh_DUNE_tau_opt_antineutrino_correction_factor_RMiss = np.sum(np.multiply(energy_DUNE_tau_opt, Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss)) / np.sum(Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss)
mean_1tau_incoh_DUNE_tau_opt_antineutrino_correction_factor_RMiss = np.sum(np.multiply(energy_DUNE_tau_opt, Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss)) / np.sum(Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss)


standard_color = '#1E3282'
tau_opt_color = '#B22222'

opacity = 0.5

ax1.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, color=standard_color, alpha=opacity, edgecolor='black', lw=0.5, zorder=2)
ax1.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, histtype='step', color=standard_color, alpha=1, lw=2, zorder=2, hatch='/')
ax1.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, color=tau_opt_color, alpha=opacity, edgecolor='black', lw=0.5)
ax1.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, histtype='step', color=tau_opt_color, alpha=1, lw=2)
ax1.text(0.98,0.90,r'\textbf{Standard} $\langle E_\nu \rangle = $'+ f' {mean_1tau_coh_DUNE_neutrino:.1f}' + r' GeV',ha='right',transform=ax1.transAxes, fontsize=30, color=standard_color)
ax1.text(0.98,0.80,r'\textbf{$\tau$-Optimized} $\langle E_\nu \rangle = $'+ f' {mean_1tau_coh_DUNE_tau_opt_neutrino:.1f}' + r' GeV',ha='right',transform=ax1.transAxes, fontsize=30, color=tau_opt_color)

ax1.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=standard_color, alpha=1, lw=3)
ax1.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=tau_opt_color, alpha=1, lw=3)

ax2.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_neutrino_vmu_67_3, color=standard_color, alpha=opacity, edgecolor='black', lw=0.5, zorder=2)
ax2.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_neutrino_vmu_67_3, histtype='step', color=standard_color, alpha=1, lw=2, zorder=2, hatch='/')
ax2.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, color=tau_opt_color, alpha=opacity, edgecolor='black', lw=0.5)
ax2.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3, histtype='step', color=tau_opt_color, alpha=1, lw=2)
ax2.text(0.98,0.90,r'\textbf{Standard} $\langle E_\nu \rangle = $'+ f' {mean_1tau_incoh_DUNE_neutrino:.1f}' + r' GeV',ha='right',transform=ax2.transAxes, fontsize=30, color=standard_color)
ax2.text(0.98,0.80,r'\textbf{$\tau$-Optimized} $\langle E_\nu \rangle = $'+ f' {mean_1tau_incoh_DUNE_tau_opt_neutrino:.1f}' + r' GeV',ha='right',transform=ax2.transAxes, fontsize=30, color=tau_opt_color)

ax2.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=standard_color, alpha=1, lw=3)
ax2.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=tau_opt_color, alpha=1, lw=3)

ax3.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3, color=standard_color, alpha=opacity, edgecolor='black', lw=0.5, zorder=2)
ax3.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3, histtype='step', color=standard_color, alpha=1, lw=2, zorder=2, hatch='/')
ax3.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, color=tau_opt_color, alpha=opacity, edgecolor='black', lw=0.5)
ax3.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3, histtype='step', color=tau_opt_color, alpha=1, lw=2)
ax3.text(0.02,0.90,r'\textbf{Standard} $\langle E_\nu \rangle = $'+ f' {mean_1tau_coh_DUNE_antineutrino:.1f}' + r' GeV',ha='left',transform=ax3.transAxes, fontsize=29, color=standard_color)
ax3.text(0.02,0.80,r'\textbf{$\tau$-Optimized} $\langle E_\nu \rangle = $'+ f' {mean_1tau_coh_DUNE_tau_opt_antineutrino:.1f}' + r' GeV',ha='left',transform=ax3.transAxes, fontsize=29, color=tau_opt_color)

ax3.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=standard_color, alpha=1, lw=3)
ax3.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=tau_opt_color, alpha=1, lw=3)

ax4.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_antineutrino_vmu_67_3, color=standard_color, alpha=opacity, edgecolor='black', lw=0.5, zorder=2)
ax4.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_antineutrino_vmu_67_3, histtype='step', color=standard_color, alpha=1, lw=2, zorder=2, hatch='/')
ax4.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, color=tau_opt_color, alpha=opacity, edgecolor='black', lw=0.5)
ax4.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3, histtype='step', color=tau_opt_color, alpha=1, lw=2)
ax4.text(0.98,0.90,r'\textbf{Standard} $\langle E_\nu \rangle = $'+ f' {mean_1tau_incoh_DUNE_antineutrino:.1f}' + r' GeV',ha='right',transform=ax4.transAxes, fontsize=30, color=standard_color)
ax4.text(0.98,0.80,r'\textbf{$\tau$-Optimized} $\langle E_\nu \rangle = $'+ f' {mean_1tau_incoh_DUNE_tau_opt_antineutrino:.1f}' + r' GeV',ha='right',transform=ax4.transAxes, fontsize=30, color=tau_opt_color)

ax4.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=standard_color, alpha=1, lw=3)
ax4.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=tau_opt_color, alpha=1, lw=3)

ax5.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_EPA_Heaviside_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, color=standard_color, alpha=opacity, edgecolor='black', lw=0.5, zorder=2)
ax5.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_EPA_Heaviside_1tau_coh_Ar_DUNE_neutrino_vmu_67_3, histtype='step', color=standard_color, alpha=1, lw=2, zorder=2, hatch='/')
ax5.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_EPA_Heaviside_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, color=tau_opt_color, alpha=opacity, edgecolor='black', lw=0.5)
ax5.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_EPA_Heaviside_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3, histtype='step', color=tau_opt_color, alpha=1, lw=2)
ax5.text(0.98,0.90,r'\textbf{Standard} $\langle E_\nu \rangle = $'+ f' {mean_EPA_Heaviside_1tau_coh_DUNE_neutrino:.1f}' + r' GeV',ha='right',transform=ax1.transAxes, fontsize=30, color=standard_color)
ax5.text(0.98,0.80,r'\textbf{$\tau$-Optimized} $\langle E_\nu \rangle = $'+ f' {mean_EPA_Heaviside_1tau_coh_DUNE_tau_opt_neutrino:.1f}' + r' GeV',ha='right',transform=ax1.transAxes, fontsize=30, color=tau_opt_color)

# with correction factors #
ax6.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, color=standard_color, alpha=opacity, edgecolor='black', lw=0.5, zorder=2)
ax6.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=standard_color, alpha=1, lw=2, zorder=2, hatch='/')
ax6.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, color=tau_opt_color, alpha=opacity, edgecolor='black', lw=0.5)
ax6.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=tau_opt_color, alpha=1, lw=2)
ax6.text(0.98,0.90,r'\textbf{Standard} $\langle E_\nu \rangle = $'+ f' {mean_1tau_coh_DUNE_neutrino_correction_factor_RMiss:.1f}' + r' GeV',ha='right',transform=ax6.transAxes, fontsize=30, color=standard_color)
ax6.text(0.98,0.80,r'\textbf{$\tau$-Optimized} $\langle E_\nu \rangle = $'+ f' {mean_1tau_coh_DUNE_tau_opt_neutrino_correction_factor_RMiss:.1f}' + r' GeV',ha='right',transform=ax6.transAxes, fontsize=30, color=tau_opt_color)

ax7.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, color=standard_color, alpha=opacity, edgecolor='black', lw=0.5, zorder=2)
ax7.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_neutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=standard_color, alpha=1, lw=2, zorder=2, hatch='/')
ax7.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, color=tau_opt_color, alpha=opacity, edgecolor='black', lw=0.5)
ax7.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_tau_opt_neutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=tau_opt_color, alpha=1, lw=2)
ax7.text(0.98,0.90,r'\textbf{Standard} $\langle E_\nu \rangle = $'+ f' {mean_1tau_incoh_DUNE_neutrino_correction_factor_RMiss:.1f}' + r' GeV',ha='right',transform=ax7.transAxes, fontsize=30, color=standard_color)
ax7.text(0.98,0.80,r'\textbf{$\tau$-Optimized} $\langle E_\nu \rangle = $'+ f' {mean_1tau_incoh_DUNE_tau_opt_neutrino_correction_factor_RMiss:.1f}' + r' GeV',ha='right',transform=ax7.transAxes, fontsize=30, color=tau_opt_color)

ax8.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, color=standard_color, alpha=opacity, edgecolor='black', lw=0.5, zorder=2)
ax8.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=standard_color, alpha=1, lw=2, zorder=2, hatch='/')
ax8.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, color=tau_opt_color, alpha=opacity, edgecolor='black', lw=0.5)
ax8.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_coh_Ar_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=tau_opt_color, alpha=1, lw=2)
ax8.text(0.02,0.90,r'\textbf{Standard} $\langle E_\nu \rangle = $'+ f' {mean_1tau_coh_DUNE_antineutrino_correction_factor_RMiss:.1f}' + r' GeV',ha='left',transform=ax8.transAxes, fontsize=29, color=standard_color)
ax8.text(0.02,0.80,r'\textbf{$\tau$-Optimized} $\langle E_\nu \rangle = $'+ f' {mean_1tau_coh_DUNE_tau_opt_antineutrino_correction_factor_RMiss:.1f}' + r' GeV',ha='left',transform=ax8.transAxes, fontsize=29, color=tau_opt_color)

ax9.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, color=standard_color, alpha=opacity, edgecolor='black', lw=0.5, zorder=2)
ax9.hist(energy_DUNE, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_antineutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=standard_color, alpha=1, lw=2, zorder=2, hatch='/')
ax9.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, color=tau_opt_color, alpha=opacity, edgecolor='black', lw=0.5)
ax9.hist(energy_DUNE_tau_opt, density=False, bins=40, range=(0.0,125.375), weights=Nbin_1tau_incoh_DUNE_tau_opt_antineutrino_vmu_67_3_correction_factor_RMiss, histtype='step', color=tau_opt_color, alpha=1, lw=2)
ax9.text(0.98,0.90,r'\textbf{Standard} $\langle E_\nu \rangle = $'+ f' {mean_1tau_incoh_DUNE_antineutrino_correction_factor_RMiss:.1f}' + r' GeV',ha='right',transform=ax9.transAxes, fontsize=30, color=standard_color)
ax9.text(0.98,0.80,r'\textbf{$\tau$-Optimized} $\langle E_\nu \rangle = $'+ f' {mean_1tau_incoh_DUNE_tau_opt_antineutrino_correction_factor_RMiss:.1f}' + r' GeV',ha='right',transform=ax9.transAxes, fontsize=30, color=tau_opt_color)

for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]:
    ax.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^- (N_\mathrm{events} / \mathrm{bin})$', fontsize=40)
    ax.set_xlabel(r'\textbf{Neutrino Energy} $E_\nu$ (GeV) ', fontsize=40)
    ax.set_xlim(0.0, 100.0)

ax1.set_title(r"\textbf{DUNE $\nu$ Mode - Coherent Scattering}", fontsize=40)
ax2.set_title(r"\textbf{DUNE $\nu$ Mode - Incoherent Scattering}", fontsize=40)
ax3.set_title(r"\textbf{DUNE $\bar{\nu}$ Mode - Coherent Scattering}", fontsize=40)
ax4.set_title(r"\textbf{DUNE $\bar{\nu}$ Mode - Incoherent Scattering}", fontsize=40)

ax5.set_title(r"\textbf{DUNE $\nu$ Mode - Coherent Scattering}", fontsize=40)

ax6.set_title(r"\textbf{DUNE $\nu$ Mode - Coherent Scattering}", fontsize=40)
ax7.set_title(r"\textbf{DUNE $\nu$ Mode - Incoherent Scattering}", fontsize=40)
ax8.set_title(r"\textbf{DUNE $\bar{\nu}$ Mode - Coherent Scattering}", fontsize=40)
ax9.set_title(r"\textbf{DUNE $\bar{\nu}$ Mode - Incoherent Scattering}", fontsize=40)

fig1.savefig("../plots/Nbin_1tau_coh_DUNE_neutrino.png", dpi=200, bbox_inches='tight')
fig2.savefig("../plots/Nbin_1tau_incoh_DUNE_neutrino.png", dpi=200, bbox_inches='tight')
fig3.savefig("../plots/Nbin_1tau_coh_DUNE_antineutrino.png", dpi=200, bbox_inches='tight')
fig4.savefig("../plots/Nbin_1tau_incoh_DUNE_antineutrino.png", dpi=200, bbox_inches='tight')

fig5.savefig("../plots/Nbin_EPA_Heaviside_1tau_coh_DUNE_neutrino.png", dpi=200, bbox_inches='tight')

fig6.savefig("../plots/Nbin_1tau_coh_DUNE_neutrino_correction_factor_RMiss.png", dpi=200, bbox_inches='tight')
fig7.savefig("../plots/Nbin_1tau_incoh_DUNE_neutrino_correction_factor_RMiss.png", dpi=200, bbox_inches='tight')
fig8.savefig("../plots/Nbin_1tau_coh_DUNE_antineutrino_correction_factor_RMiss.png", dpi=200, bbox_inches='tight')
fig9.savefig("../plots/Nbin_1tau_incoh_DUNE_antineutrino_correction_factor_RMiss.png", dpi=200, bbox_inches='tight')
