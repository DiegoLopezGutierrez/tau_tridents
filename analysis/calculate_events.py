import numpy as np
import csv
from scipy.integrate import simpson
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

## DUNE ##
N_POT = 1.1e21               # Number of POTs
MAr = 39.95*atomic_mass_unit # Mass of argon in atomic mass units
A_Ar = 40
Z_Ar = 18

#Phi_Alt = 1.04e-3

## FASERnu ##
L_Run = 150                  # Run luminosity at FASERv; 150 fb^-1
MW = 183.84*atomic_mass_unit # Mass of tungsten in atomic mass units
A_W = 184
Z_W = 74

## Cross Section Uncertainties ##
# Components # (see Altmannshofer et al. for details)
aEM = 1/137  # low q^2 usually, so using zero momentum value of fine structure
sigma_highQED_Ar = Z_Ar*aEM/(4*np.pi) + 0.02  # higher order QED corrections should give roughly Z*aEM/(4*pi) which is about 1% for Ar. Add 2% to be conservative.
sigma_highQED_W  = Z_W*aEM/(4*np.pi) + 0.02 # higher order QED corrections should give roughly Z*aEM/(4*pi) which is about 4% for W. Add 2% to be conservative.
sigma_highQED_p = 1*aEM/(4*np.pi) + 0.02
sigma_highQED_n = 0

sigma_form_factors_coh = 0.01 # 1% uncertainty on coherent cross sections due to form factors
sigma_form_factors_incoh = 0.03 # 3% uncertainty on incoherent cross sections due to form factors

sigma_highEW = 0.05 # 5% uncertainty due to higher order weak corrections

sigma_nuclear_modeling = 0.30 # for incoherent scattering, largest uncertainty due to other nuclear effects besides the included Pauli blocking. 30% to be conservative.

# Total #
sigma_total_coh_Ar = np.sqrt(sigma_highQED_Ar**2 + sigma_form_factors_coh**2 + sigma_highEW**2)
sigma_total_coh_W  = np.sqrt(sigma_highQED_W**2 + sigma_form_factors_coh**2 + sigma_highEW**2)
print(sigma_total_coh_Ar, sigma_total_coh_W)

sigma_total_incoh_Ar = np.sqrt(sigma_highQED_Ar**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)
sigma_total_incoh_W  = np.sqrt(sigma_highQED_W**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)
print(sigma_total_incoh_Ar, sigma_total_incoh_W)

sigma_total_p = np.sqrt(sigma_highQED_p**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)
sigma_total_n = np.sqrt(sigma_highQED_n**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)
print(sigma_total_p, sigma_total_n)

########################
###### Initialize ######
########################

### Directories ###
FLUX_DIR = '../csv/fluxes'
CROSS_SECTION_DIR = '../csv/cross_sections'

### Fluxes ###
DUNE_filename = FLUX_DIR + '/DUNE/histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_LBNEND_globes_flux.txt'
DUNE_tau_opt_numu_flux = FLUX_DIR + '/DUNE/histos_g4lbne_v3r5p4_QGSP_BERT_TauOptimized_neutrino_LBNEND_globes_flux.txt'
FASERnu_filename = FLUX_DIR + '/FASERnu/vmu/FASERvmu.csv'

energy_DUNE = []
flux_DUNE = []

energy_DUNE_tau_opt = []
flux_DUNE_tau_opt = []

energy_FASERvmu = []
flux_FASERvmu = []
bins_FASERvmu = []

# Get bin edges. This is based on a visual inspection of the DUNE plot.
with open(FLUX_DIR + '/FASERnu/FASERvmu_bins.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        bins_FASERvmu.append(float(row[0]))

bins_FASERvmu.append(8150.00)
bins_FASERvmu = np.array(bins_FASERvmu)


### Cross Sections ###
# Filenames #
xsec_1tau_coh_Ar_filename = CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/coherent/argon/vmu_to_vtau_tau+_mu-_coh_Ar_xsec.csv'
xsec_1tau_coh_W_filename = CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/coherent/tungsten/vmu_to_vtau_tau+_mu-_coh_W_xsec.csv'
xsec_1tau_p_filename = CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv'
xsec_1tau_n_filename = CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/new_neutron/vmu_to_vtau_tau+_mu-_nucleon_n_xsec.csv'

xsec_2tau_coh_Ar_filename = CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/coherent/argon/vmu_to_vmu_tau+_tau-_coh_Ar_xsec.csv'
xsec_2tau_coh_W_filename = CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/coherent/tungsten/vmu_to_vmu_tau+_tau-_coh_W_xsec.csv'
xsec_2tau_p_filename = CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv'
xsec_2tau_n_filename = CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/new_neutron/vmu_to_vmu_tau+_tau-_nucleon_n_xsec.csv'

xsec_2mu_coh_Ar_filename  = CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/argon/vmu_to_vmu_mu+_mu-_coh_Ar_xsec.csv'
xsec_2mu_p_filename  = CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/proton/vmu_to_vmu_mu+_mu-_nucleon_p_xsec.csv'
xsec_2mu_n_filename  = CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/neutron/vmu_to_vmu_mu+_mu-_nucleon_n_xsec.csv'

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

# Coherent ; Tungsten #
energy_1tau_coh_W = []
xsec_1tau_coh_W = []
delta_1tau_coh_W = []

energy_2tau_coh_W = []
xsec_2tau_coh_W = []
delta_2tau_coh_W = []

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

########################
###### Load Files ######
########################

##### Fluxes #####
### DUNE ; Standard Mode Flux ###
with open(DUNE_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ' ')
    for row in data:
        energy_DUNE.append(float(row[0]))
        flux_DUNE.append(float(row[2])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]

### DUNE ; Tau-optimized Flux ###
with open(DUNE_tau_opt_numu_flux,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ' ')
    for row in data:
        energy_DUNE_tau_opt.append(float(row[0]))
        flux_DUNE_tau_opt.append(float(row[2])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]

### FASERvmu ; vmu flux ###
with open(FASERnu_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        energy = float(row[0]) # Energy in GeV
        flux = float(row[1]) # Flux has already been normalized by the 25 x 25 cm^2 area and 150 fb^-1; it is now in units of [m^-2 GeV^-1 fb].
        energy_FASERvmu.append(energy)
        flux_FASERvmu.append(flux)

### Integrated Fluxes ###
DUNE_integrated_flux = simpson(flux_DUNE, x=energy_DUNE)                            # Integrated flux [Nv / m^2 POT yr]
DUNE_tau_opt_integrated_flux = simpson(flux_DUNE_tau_opt, x=energy_DUNE_tau_opt)    # Integrated flux [Nv / m^2 POT yr]
#Alt120_integrated_flux = Phi_Alt
#Alt_integrated_flux = Phi_Alt
FASER_integrated_flux = simpson(flux_FASERvmu, x=energy_FASERvmu)                   # Integrated flux [Nv / m^2 fb^-1]

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

xsec_1tau_coh_W_lower, xsec_1tau_coh_W_upper = Limits(xsec_1tau_coh_W, delta_1tau_coh_W)
xsec_1tau_p_W_lower, xsec_1tau_p_W_upper = Limits(xsec_1tau_p_W, delta_1tau_p_W)
xsec_1tau_n_W_lower, xsec_1tau_n_W_upper = Limits(xsec_1tau_n_W, delta_1tau_n_W)

xsec_2tau_coh_W_lower, xsec_2tau_coh_W_upper = Limits(xsec_2tau_coh_W, delta_2tau_coh_W)
xsec_2tau_p_W_lower, xsec_2tau_p_W_upper = Limits(xsec_2tau_p_W, delta_2tau_p_W)
xsec_2tau_n_W_lower, xsec_2tau_n_W_upper = Limits(xsec_2tau_n_W, delta_2tau_n_W)

xsec_2mu_p_W_lower, xsec_2mu_p_W_upper = Limits(xsec_2mu_p_W, delta_2mu_p_W)
xsec_2mu_n_W_lower, xsec_2mu_n_W_upper = Limits(xsec_2mu_n_W, delta_2mu_n_W)

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

### FASERnu vmu Flux ; Matched XSec ###
FASERvmu_xsec_1tau_coh_W_matched = MatchXSec(energy_FASERvmu, energy_1tau_coh_W, xsec_1tau_coh_W)
FASERvmu_xsec_1tau_p_W_matched = MatchXSec(energy_FASERvmu, energy_1tau_p_W, xsec_1tau_p_W)
FASERvmu_xsec_1tau_n_W_matched = MatchXSec(energy_FASERvmu, energy_1tau_n_W, xsec_1tau_n_W)

FASERvmu_xsec_2tau_coh_W_matched = MatchXSec(energy_FASERvmu, energy_2tau_coh_W, xsec_2tau_coh_W)
FASERvmu_xsec_2tau_p_W_matched = MatchXSec(energy_FASERvmu, energy_2tau_p_W, xsec_2tau_p_W)
FASERvmu_xsec_2tau_n_W_matched = MatchXSec(energy_FASERvmu, energy_2tau_n_W, xsec_2tau_n_W)

FASERvmu_xsec_2mu_p_W_matched = MatchXSec(energy_FASERvmu, energy_2mu_p_W, xsec_2mu_p_W)
FASERvmu_xsec_2mu_n_W_matched = MatchXSec(energy_FASERvmu, energy_2mu_n_W, xsec_2mu_n_W)

#print('XSec matched: ',len(FASERvmu_xsec_1tau_coh_W_matched), FASERvmu_xsec_1tau_coh_W_matched)
#print('Flux energy: ',len(energy_FASERvmu), energy_FASERvmu)
#print('Flux: ',len(flux_FASERvmu), flux_FASERvmu)

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

FASERvmu_xsec_2mu_p_W_matched_upper = MatchXSec(energy_FASERvmu, energy_2mu_p_W, xsec_2mu_p_W_upper)
FASERvmu_xsec_2mu_n_W_matched_upper = MatchXSec(energy_FASERvmu, energy_2mu_n_W, xsec_2mu_n_W_upper)

FASERvmu_xsec_2mu_p_W_matched_lower = MatchXSec(energy_FASERvmu, energy_2mu_p_W, xsec_2mu_p_W_lower)
FASERvmu_xsec_2mu_n_W_matched_lower = MatchXSec(energy_FASERvmu, energy_2mu_n_W, xsec_2mu_n_W_lower)

##########################################
###### Number of Events Calculation ######
##########################################

def CalculateEvents(flux, flux_energy, xsec_matched, NTONNES=1, NYEAR=1, normalized=False, Phi=1, detector='DUNE'):
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
    integrand = np.multiply(flux, xsec_matched)
    XSecCon = simpson(integrand, x=flux_energy)
    
    if not normalized:
        Phi = simpson(flux, x=flux_energy)
        if detector == 'DUNE':
            N = XSecCon * (MD * NTONNES / MAr) * (N_POT * NYEAR)
        if detector == 'FASER':
            N = XSecCon * (MD * NTONNES / MW) * L_Run
        return N, XSecCon*1e43 / Phi 
    else:
        if detector == 'DUNE':
            N = Phi * XSecCon * MD * N_POT * NTONNES * NYEAR / MAr
        return N, XSecCon*1e43

### DUNE Standard Flux Events ###
## Coherent ; Argon ##
N_1tau_coh_Ar_DUNE_1_1, XSecCon_1tau_coh_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched)
N_1tau_coh_Ar_DUNE_67_3, XSecCon_1tau_coh_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_coh_Ar_DUNE_147_3, XSecCon_1tau_coh_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched, NTONNES=147, NYEAR=3)

N_2tau_coh_Ar_DUNE_1_1, XSecCon_2tau_coh_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched)
N_2tau_coh_Ar_DUNE_67_3, XSecCon_2tau_coh_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_coh_Ar_DUNE_147_3, XSecCon_2tau_coh_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched, NTONNES=147, NYEAR=3)

N_2mu_coh_Ar_DUNE_1_1, XSecCon_2mu_coh_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched)
N_2mu_coh_Ar_DUNE_67_3, XSecCon_2mu_coh_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_coh_Ar_DUNE_147_3, XSecCon_2mu_coh_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched, NTONNES=147, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_coh_Ar_DUNE_1_1, upper_XSecCon_1tau_coh_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_upper)
upper_N_1tau_coh_Ar_DUNE_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_1tau_coh_Ar_DUNE_147_3, upper_XSecCon_1tau_coh_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_upper, NTONNES=147, NYEAR=3)

upper_N_2tau_coh_Ar_DUNE_1_1, upper_XSecCon_2tau_coh_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched_upper)
upper_N_2tau_coh_Ar_DUNE_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_coh_Ar_DUNE_147_3, upper_XSecCon_2tau_coh_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched_upper, NTONNES=147, NYEAR=3)

upper_N_2mu_coh_Ar_DUNE_1_1, upper_XSecCon_2mu_coh_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched_upper)
upper_N_2mu_coh_Ar_DUNE_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_coh_Ar_DUNE_147_3, upper_XSecCon_2mu_coh_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched_upper, NTONNES=147, NYEAR=3)

lower_N_1tau_coh_Ar_DUNE_1_1, lower_XSecCon_1tau_coh_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_lower)
lower_N_1tau_coh_Ar_DUNE_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_1tau_coh_Ar_DUNE_147_3, lower_XSecCon_1tau_coh_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched_lower, NTONNES=147, NYEAR=3)

lower_N_2tau_coh_Ar_DUNE_1_1, lower_XSecCon_2tau_coh_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched_lower)
lower_N_2tau_coh_Ar_DUNE_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_coh_Ar_DUNE_147_3, lower_XSecCon_2tau_coh_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched_lower, NTONNES=147, NYEAR=3)

lower_N_2mu_coh_Ar_DUNE_1_1, lower_XSecCon_2mu_coh_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched_lower)
lower_N_2mu_coh_Ar_DUNE_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_coh_Ar_DUNE_147_3, lower_XSecCon_2mu_coh_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched_lower, NTONNES=147, NYEAR=3)

## Incoherent ; proton ; Argon ##
N_1tau_p_Ar_DUNE_1_1, XSecCon_1tau_p_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched)
N_1tau_p_Ar_DUNE_67_3, XSecCon_1tau_p_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_p_Ar_DUNE_147_3, XSecCon_1tau_p_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched, NTONNES=147, NYEAR=3)

N_2tau_p_Ar_DUNE_1_1, XSecCon_2tau_p_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched)
N_2tau_p_Ar_DUNE_67_3, XSecCon_2tau_p_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_p_Ar_DUNE_147_3, XSecCon_2tau_p_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched, NTONNES=147, NYEAR=3)

N_2mu_p_Ar_DUNE_1_1, XSecCon_2mu_p_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched)
N_2mu_p_Ar_DUNE_67_3, XSecCon_2mu_p_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_p_Ar_DUNE_147_3, XSecCon_2mu_p_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched, NTONNES=147, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_p_Ar_DUNE_1_1, upper_XSecCon_1tau_p_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_upper)
upper_N_1tau_p_Ar_DUNE_67_3, upper_XSecCon_1tau_p_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_1tau_p_Ar_DUNE_147_3, upper_XSecCon_1tau_p_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_upper, NTONNES=147, NYEAR=3)

upper_N_2tau_p_Ar_DUNE_1_1, upper_XSecCon_2tau_p_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched_upper)
upper_N_2tau_p_Ar_DUNE_67_3, upper_XSecCon_2tau_p_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_p_Ar_DUNE_147_3, upper_XSecCon_2tau_p_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched_upper, NTONNES=147, NYEAR=3)

upper_N_2mu_p_Ar_DUNE_1_1, upper_XSecCon_2mu_p_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched_upper)
upper_N_2mu_p_Ar_DUNE_67_3, upper_XSecCon_2mu_p_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_p_Ar_DUNE_147_3, upper_XSecCon_2mu_p_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched_upper, NTONNES=147, NYEAR=3)

lower_N_1tau_p_Ar_DUNE_1_1, lower_XSecCon_1tau_p_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_lower)
lower_N_1tau_p_Ar_DUNE_67_3, lower_XSecCon_1tau_p_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_1tau_p_Ar_DUNE_147_3, lower_XSecCon_1tau_p_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_p_Ar_matched_lower, NTONNES=147, NYEAR=3)

lower_N_2tau_p_Ar_DUNE_1_1, lower_XSecCon_2tau_p_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched_lower)
lower_N_2tau_p_Ar_DUNE_67_3, lower_XSecCon_2tau_p_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_p_Ar_DUNE_147_3, lower_XSecCon_2tau_p_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_p_Ar_matched_lower, NTONNES=147, NYEAR=3)

lower_N_2mu_p_Ar_DUNE_1_1, lower_XSecCon_2mu_p_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched_lower)
lower_N_2mu_p_Ar_DUNE_67_3, lower_XSecCon_2mu_p_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_p_Ar_DUNE_147_3, lower_XSecCon_2mu_p_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_p_Ar_matched_lower, NTONNES=147, NYEAR=3)

## Incoherent ; neutron ; Argon ##
N_1tau_n_Ar_DUNE_1_1, XSecCon_1tau_n_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched)
N_1tau_n_Ar_DUNE_67_3, XSecCon_1tau_n_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_n_Ar_DUNE_147_3, XSecCon_1tau_n_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched, NTONNES=147, NYEAR=3)

N_2tau_n_Ar_DUNE_1_1, XSecCon_2tau_n_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched)
N_2tau_n_Ar_DUNE_67_3, XSecCon_2tau_n_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_n_Ar_DUNE_147_3, XSecCon_2tau_n_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched, NTONNES=147, NYEAR=3)

N_2mu_n_Ar_DUNE_1_1, XSecCon_2mu_n_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched)
N_2mu_n_Ar_DUNE_67_3, XSecCon_2mu_n_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_n_Ar_DUNE_147_3, XSecCon_2mu_n_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched, NTONNES=147, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_n_Ar_DUNE_1_1, upper_XSecCon_1tau_n_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_upper)
upper_N_1tau_n_Ar_DUNE_67_3, upper_XSecCon_1tau_n_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_1tau_n_Ar_DUNE_147_3, upper_XSecCon_1tau_n_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_upper, NTONNES=147, NYEAR=3)

upper_N_2tau_n_Ar_DUNE_1_1, upper_XSecCon_2tau_n_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched_upper)
upper_N_2tau_n_Ar_DUNE_67_3, upper_XSecCon_2tau_n_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_n_Ar_DUNE_147_3, upper_XSecCon_2tau_n_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched_upper, NTONNES=147, NYEAR=3)

upper_N_2mu_n_Ar_DUNE_1_1, upper_XSecCon_2mu_n_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched_upper)
upper_N_2mu_n_Ar_DUNE_67_3, upper_XSecCon_2mu_n_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_n_Ar_DUNE_147_3, upper_XSecCon_2mu_n_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched_upper, NTONNES=147, NYEAR=3)

lower_N_1tau_n_Ar_DUNE_1_1, lower_XSecCon_1tau_n_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_lower)
lower_N_1tau_n_Ar_DUNE_67_3, lower_XSecCon_1tau_n_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_1tau_n_Ar_DUNE_147_3, lower_XSecCon_1tau_n_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_n_Ar_matched_lower, NTONNES=147, NYEAR=3)

lower_N_2tau_n_Ar_DUNE_1_1, lower_XSecCon_2tau_n_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched_lower)
lower_N_2tau_n_Ar_DUNE_67_3, lower_XSecCon_2tau_n_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_n_Ar_DUNE_147_3, lower_XSecCon_2tau_n_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_n_Ar_matched_lower, NTONNES=147, NYEAR=3)

lower_N_2mu_n_Ar_DUNE_1_1, lower_XSecCon_2mu_n_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched_lower)
lower_N_2mu_n_Ar_DUNE_67_3, lower_XSecCon_2mu_n_Ar_DUNE_67_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_n_Ar_DUNE_147_3, lower_XSecCon_2mu_n_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_n_Ar_matched_lower, NTONNES=147, NYEAR=3)

## Incoherent ; proton + neutron ; Argon ##
N_1tau_incoh_DUNE_1_1, XSecCon_1tau_incoh_DUNE_1_1 = (N_1tau_p_Ar_DUNE_1_1 + N_1tau_n_Ar_DUNE_1_1), (XSecCon_1tau_p_Ar_DUNE_1_1 + XSecCon_1tau_n_Ar_DUNE_1_1)
N_1tau_incoh_DUNE_67_3, XSecCon_1tau_incoh_DUNE_67_3 = (N_1tau_p_Ar_DUNE_67_3 + N_1tau_n_Ar_DUNE_67_3), (XSecCon_1tau_p_Ar_DUNE_67_3 + XSecCon_1tau_n_Ar_DUNE_67_3)
N_1tau_incoh_DUNE_147_3, XSecCon_1tau_incoh_DUNE_147_3 = (N_1tau_p_Ar_DUNE_147_3 + N_1tau_n_Ar_DUNE_147_3), (XSecCon_1tau_p_Ar_DUNE_147_3 + XSecCon_1tau_n_Ar_DUNE_147_3) 

N_2tau_incoh_DUNE_1_1, XSecCon_2tau_incoh_DUNE_1_1 = (N_2tau_p_Ar_DUNE_1_1 + N_2tau_n_Ar_DUNE_1_1), (XSecCon_2tau_p_Ar_DUNE_1_1 + XSecCon_2tau_n_Ar_DUNE_1_1)
N_2tau_incoh_DUNE_67_3, XSecCon_2tau_incoh_DUNE_67_3 = (N_2tau_p_Ar_DUNE_67_3 + N_2tau_n_Ar_DUNE_67_3), (XSecCon_2tau_p_Ar_DUNE_67_3 + XSecCon_2tau_n_Ar_DUNE_67_3)
N_2tau_incoh_DUNE_147_3, XSecCon_2tau_incoh_DUNE_147_3 = (N_2tau_p_Ar_DUNE_147_3 + N_2tau_n_Ar_DUNE_147_3), (XSecCon_2tau_p_Ar_DUNE_147_3 + XSecCon_2tau_n_Ar_DUNE_147_3) 

N_2mu_incoh_DUNE_1_1, XSecCon_2mu_incoh_DUNE_1_1 = (N_2mu_p_Ar_DUNE_1_1 + N_2mu_n_Ar_DUNE_1_1), (XSecCon_2mu_p_Ar_DUNE_1_1 + XSecCon_2mu_n_Ar_DUNE_1_1)
N_2mu_incoh_DUNE_67_3, XSecCon_2mu_incoh_DUNE_67_3 = (N_2mu_p_Ar_DUNE_67_3 + N_2mu_n_Ar_DUNE_67_3), (XSecCon_2mu_p_Ar_DUNE_67_3 + XSecCon_2mu_n_Ar_DUNE_67_3)
N_2mu_incoh_DUNE_147_3, XSecCon_2mu_incoh_DUNE_147_3 = (N_2mu_p_Ar_DUNE_147_3 + N_2mu_n_Ar_DUNE_147_3), (XSecCon_2mu_p_Ar_DUNE_147_3 + XSecCon_2mu_n_Ar_DUNE_147_3) 

# Upper and Lower limits #
upper_N_1tau_incoh_DUNE_1_1, upper_XSecCon_1tau_incoh_DUNE_1_1 = (upper_N_1tau_p_Ar_DUNE_1_1 + upper_N_1tau_n_Ar_DUNE_1_1), (upper_XSecCon_1tau_p_Ar_DUNE_1_1 + upper_XSecCon_1tau_n_Ar_DUNE_1_1)
upper_N_1tau_incoh_DUNE_67_3, upper_XSecCon_1tau_incoh_DUNE_67_3 = (upper_N_1tau_p_Ar_DUNE_67_3 + upper_N_1tau_n_Ar_DUNE_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_67_3)
upper_N_1tau_incoh_DUNE_147_3, upper_XSecCon_1tau_incoh_DUNE_147_3 = (upper_N_1tau_p_Ar_DUNE_147_3 + upper_N_1tau_n_Ar_DUNE_147_3), (upper_XSecCon_1tau_p_Ar_DUNE_147_3 + upper_XSecCon_1tau_n_Ar_DUNE_147_3) 

upper_N_2tau_incoh_DUNE_1_1, upper_XSecCon_2tau_incoh_DUNE_1_1 = (upper_N_2tau_p_Ar_DUNE_1_1 + upper_N_2tau_n_Ar_DUNE_1_1), (upper_XSecCon_2tau_p_Ar_DUNE_1_1 + upper_XSecCon_2tau_n_Ar_DUNE_1_1)
upper_N_2tau_incoh_DUNE_67_3, upper_XSecCon_2tau_incoh_DUNE_67_3 = (upper_N_2tau_p_Ar_DUNE_67_3 + upper_N_2tau_n_Ar_DUNE_67_3), (upper_XSecCon_2tau_p_Ar_DUNE_67_3 + upper_XSecCon_2tau_n_Ar_DUNE_67_3)
upper_N_2tau_incoh_DUNE_147_3, upper_XSecCon_2tau_incoh_DUNE_147_3 = (upper_N_2tau_p_Ar_DUNE_147_3 + upper_N_2tau_n_Ar_DUNE_147_3), (upper_XSecCon_2tau_p_Ar_DUNE_147_3 + upper_XSecCon_2tau_n_Ar_DUNE_147_3) 

upper_N_2mu_incoh_DUNE_1_1, upper_XSecCon_2mu_incoh_DUNE_1_1 = (upper_N_2mu_p_Ar_DUNE_1_1 + upper_N_2mu_n_Ar_DUNE_1_1), (upper_XSecCon_2mu_p_Ar_DUNE_1_1 + upper_XSecCon_2mu_n_Ar_DUNE_1_1)
upper_N_2mu_incoh_DUNE_67_3, upper_XSecCon_2mu_incoh_DUNE_67_3 = (upper_N_2mu_p_Ar_DUNE_67_3 + upper_N_2mu_n_Ar_DUNE_67_3), (upper_XSecCon_2mu_p_Ar_DUNE_67_3 + upper_XSecCon_2mu_n_Ar_DUNE_67_3)
upper_N_2mu_incoh_DUNE_147_3, upper_XSecCon_2mu_incoh_DUNE_147_3 = (upper_N_2mu_p_Ar_DUNE_147_3 + upper_N_2mu_n_Ar_DUNE_147_3), (upper_XSecCon_2mu_p_Ar_DUNE_147_3 + upper_XSecCon_2mu_n_Ar_DUNE_147_3) 

lower_N_1tau_incoh_DUNE_1_1, lower_XSecCon_1tau_incoh_DUNE_1_1 = (lower_N_1tau_p_Ar_DUNE_1_1 + lower_N_1tau_n_Ar_DUNE_1_1), (lower_XSecCon_1tau_p_Ar_DUNE_1_1 + lower_XSecCon_1tau_n_Ar_DUNE_1_1)
lower_N_1tau_incoh_DUNE_67_3, lower_XSecCon_1tau_incoh_DUNE_67_3 = (lower_N_1tau_p_Ar_DUNE_67_3 + lower_N_1tau_n_Ar_DUNE_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_67_3)
lower_N_1tau_incoh_DUNE_147_3, lower_XSecCon_1tau_incoh_DUNE_147_3 = (lower_N_1tau_p_Ar_DUNE_147_3 + lower_N_1tau_n_Ar_DUNE_147_3), (lower_XSecCon_1tau_p_Ar_DUNE_147_3 + lower_XSecCon_1tau_n_Ar_DUNE_147_3) 

lower_N_2tau_incoh_DUNE_1_1, lower_XSecCon_2tau_incoh_DUNE_1_1 = (lower_N_2tau_p_Ar_DUNE_1_1 + lower_N_2tau_n_Ar_DUNE_1_1), (lower_XSecCon_2tau_p_Ar_DUNE_1_1 + lower_XSecCon_2tau_n_Ar_DUNE_1_1)
lower_N_2tau_incoh_DUNE_67_3, lower_XSecCon_2tau_incoh_DUNE_67_3 = (lower_N_2tau_p_Ar_DUNE_67_3 + lower_N_2tau_n_Ar_DUNE_67_3), (lower_XSecCon_2tau_p_Ar_DUNE_67_3 + lower_XSecCon_2tau_n_Ar_DUNE_67_3)
lower_N_2tau_incoh_DUNE_147_3, lower_XSecCon_2tau_incoh_DUNE_147_3 = (lower_N_2tau_p_Ar_DUNE_147_3 + lower_N_2tau_n_Ar_DUNE_147_3), (lower_XSecCon_2tau_p_Ar_DUNE_147_3 + lower_XSecCon_2tau_n_Ar_DUNE_147_3) 

lower_N_2mu_incoh_DUNE_1_1, lower_XSecCon_2mu_incoh_DUNE_1_1 = (lower_N_2mu_p_Ar_DUNE_1_1 + lower_N_2mu_n_Ar_DUNE_1_1), (lower_XSecCon_2mu_p_Ar_DUNE_1_1 + lower_XSecCon_2mu_n_Ar_DUNE_1_1)
lower_N_2mu_incoh_DUNE_67_3, lower_XSecCon_2mu_incoh_DUNE_67_3 = (lower_N_2mu_p_Ar_DUNE_67_3 + lower_N_2mu_n_Ar_DUNE_67_3), (lower_XSecCon_2mu_p_Ar_DUNE_67_3 + lower_XSecCon_2mu_n_Ar_DUNE_67_3)
lower_N_2mu_incoh_DUNE_147_3, lower_XSecCon_2mu_incoh_DUNE_147_3 = (lower_N_2mu_p_Ar_DUNE_147_3 + lower_N_2mu_n_Ar_DUNE_147_3), (lower_XSecCon_2mu_p_Ar_DUNE_147_3 + lower_XSecCon_2mu_n_Ar_DUNE_147_3) 


### DUNE Tau-Optimized ND Flux Events ###
## Coherent ; Argon ##
N_1tau_coh_Ar_DUNE_tau_opt_1_1, XSecCon_1tau_coh_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched)
N_1tau_coh_Ar_DUNE_tau_opt_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_coh_Ar_DUNE_tau_opt_147_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched, NTONNES=147, NYEAR=3)

N_2tau_coh_Ar_DUNE_tau_opt_1_1, XSecCon_2tau_coh_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched)
N_2tau_coh_Ar_DUNE_tau_opt_67_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_coh_Ar_DUNE_tau_opt_147_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched, NTONNES=147, NYEAR=3)

N_2mu_coh_Ar_DUNE_tau_opt_1_1, XSecCon_2mu_coh_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched)
N_2mu_coh_Ar_DUNE_tau_opt_67_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_coh_Ar_DUNE_tau_opt_147_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched, NTONNES=147, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_coh_Ar_DUNE_tau_opt_1_1, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_upper)
upper_N_1tau_coh_Ar_DUNE_tau_opt_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_1tau_coh_Ar_DUNE_tau_opt_147_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_upper, NTONNES=147, NYEAR=3)

upper_N_2tau_coh_Ar_DUNE_tau_opt_1_1, upper_XSecCon_2tau_coh_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched_upper)
upper_N_2tau_coh_Ar_DUNE_tau_opt_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_coh_Ar_DUNE_tau_opt_147_3, upper_XSecCon_2tau_coh_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched_upper, NTONNES=147, NYEAR=3)

upper_N_2mu_coh_Ar_DUNE_tau_opt_1_1, upper_XSecCon_2mu_coh_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched_upper)
upper_N_2mu_coh_Ar_DUNE_tau_opt_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_coh_Ar_DUNE_tau_opt_147_3, upper_XSecCon_2mu_coh_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched_upper, NTONNES=147, NYEAR=3)

lower_N_1tau_coh_Ar_DUNE_tau_opt_1_1, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_lower)
lower_N_1tau_coh_Ar_DUNE_tau_opt_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_1tau_coh_Ar_DUNE_tau_opt_147_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched_lower, NTONNES=147, NYEAR=3)

lower_N_2tau_coh_Ar_DUNE_tau_opt_1_1, lower_XSecCon_2tau_coh_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched_lower)
lower_N_2tau_coh_Ar_DUNE_tau_opt_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_coh_Ar_DUNE_tau_opt_147_3, lower_XSecCon_2tau_coh_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched_lower, NTONNES=147, NYEAR=3)

lower_N_2mu_coh_Ar_DUNE_tau_opt_1_1, lower_XSecCon_2mu_coh_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched_lower)
lower_N_2mu_coh_Ar_DUNE_tau_opt_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_coh_Ar_DUNE_tau_opt_147_3, lower_XSecCon_2mu_coh_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched_lower, NTONNES=147, NYEAR=3)

## Incoherent ; proton ; Argon ##
N_1tau_p_Ar_DUNE_tau_opt_1_1, XSecCon_1tau_p_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched)
N_1tau_p_Ar_DUNE_tau_opt_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_p_Ar_DUNE_tau_opt_147_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched, NTONNES=147, NYEAR=3)

N_2tau_p_Ar_DUNE_tau_opt_1_1, XSecCon_2tau_p_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched)
N_2tau_p_Ar_DUNE_tau_opt_67_3, XSecCon_2tau_p_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_p_Ar_DUNE_tau_opt_147_3, XSecCon_2tau_p_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched, NTONNES=147, NYEAR=3)

N_2mu_p_Ar_DUNE_tau_opt_1_1, XSecCon_2mu_p_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched)
N_2mu_p_Ar_DUNE_tau_opt_67_3, XSecCon_2mu_p_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_p_Ar_DUNE_tau_opt_147_3, XSecCon_2mu_p_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched, NTONNES=147, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_p_Ar_DUNE_tau_opt_1_1, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_upper)
upper_N_1tau_p_Ar_DUNE_tau_opt_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_1tau_p_Ar_DUNE_tau_opt_147_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_upper, NTONNES=147, NYEAR=3)

upper_N_2tau_p_Ar_DUNE_tau_opt_1_1, upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched_upper)
upper_N_2tau_p_Ar_DUNE_tau_opt_67_3, upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_p_Ar_DUNE_tau_opt_147_3, upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched_upper, NTONNES=147, NYEAR=3)

upper_N_2mu_p_Ar_DUNE_tau_opt_1_1, upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched_upper)
upper_N_2mu_p_Ar_DUNE_tau_opt_67_3, upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_p_Ar_DUNE_tau_opt_147_3, upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched_upper, NTONNES=147, NYEAR=3)

lower_N_1tau_p_Ar_DUNE_tau_opt_1_1, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_lower)
lower_N_1tau_p_Ar_DUNE_tau_opt_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_1tau_p_Ar_DUNE_tau_opt_147_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_Ar_matched_lower, NTONNES=147, NYEAR=3)

lower_N_2tau_p_Ar_DUNE_tau_opt_1_1, lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched_lower)
lower_N_2tau_p_Ar_DUNE_tau_opt_67_3, lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_p_Ar_DUNE_tau_opt_147_3, lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_Ar_matched_lower, NTONNES=147, NYEAR=3)

lower_N_2mu_p_Ar_DUNE_tau_opt_1_1, lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched_lower)
lower_N_2mu_p_Ar_DUNE_tau_opt_67_3, lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_p_Ar_DUNE_tau_opt_147_3, lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_Ar_matched_lower, NTONNES=147, NYEAR=3)

## Incoherent ; neutron ; Argon ##
N_1tau_n_Ar_DUNE_tau_opt_1_1, XSecCon_1tau_n_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched)
N_1tau_n_Ar_DUNE_tau_opt_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_1tau_n_Ar_DUNE_tau_opt_147_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched, NTONNES=147, NYEAR=3)

N_2tau_n_Ar_DUNE_tau_opt_1_1, XSecCon_2tau_n_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched)
N_2tau_n_Ar_DUNE_tau_opt_67_3, XSecCon_2tau_n_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2tau_n_Ar_DUNE_tau_opt_147_3, XSecCon_2tau_n_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched, NTONNES=147, NYEAR=3)

N_2mu_n_Ar_DUNE_tau_opt_1_1, XSecCon_2mu_n_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched)
N_2mu_n_Ar_DUNE_tau_opt_67_3, XSecCon_2mu_n_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched, NTONNES=67, NYEAR=3)
N_2mu_n_Ar_DUNE_tau_opt_147_3, XSecCon_2mu_n_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched, NTONNES=147, NYEAR=3)

# Upper and Lower limits #
upper_N_1tau_n_Ar_DUNE_tau_opt_1_1, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_upper)
upper_N_1tau_n_Ar_DUNE_tau_opt_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_1tau_n_Ar_DUNE_tau_opt_147_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_upper, NTONNES=147, NYEAR=3)

upper_N_2tau_n_Ar_DUNE_tau_opt_1_1, upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched_upper)
upper_N_2tau_n_Ar_DUNE_tau_opt_67_3, upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2tau_n_Ar_DUNE_tau_opt_147_3, upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched_upper, NTONNES=147, NYEAR=3)

upper_N_2mu_n_Ar_DUNE_tau_opt_1_1, upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched_upper)
upper_N_2mu_n_Ar_DUNE_tau_opt_67_3, upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched_upper, NTONNES=67, NYEAR=3)
upper_N_2mu_n_Ar_DUNE_tau_opt_147_3, upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched_upper, NTONNES=147, NYEAR=3)

lower_N_1tau_n_Ar_DUNE_tau_opt_1_1, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_lower)
lower_N_1tau_n_Ar_DUNE_tau_opt_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_1tau_n_Ar_DUNE_tau_opt_147_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_Ar_matched_lower, NTONNES=147, NYEAR=3)

lower_N_2tau_n_Ar_DUNE_tau_opt_1_1, lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched_lower)
lower_N_2tau_n_Ar_DUNE_tau_opt_67_3, lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2tau_n_Ar_DUNE_tau_opt_147_3, lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_Ar_matched_lower, NTONNES=147, NYEAR=3)

lower_N_2mu_n_Ar_DUNE_tau_opt_1_1, lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched_lower)
lower_N_2mu_n_Ar_DUNE_tau_opt_67_3, lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_67_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched_lower, NTONNES=67, NYEAR=3)
lower_N_2mu_n_Ar_DUNE_tau_opt_147_3, lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_Ar_matched_lower, NTONNES=147, NYEAR=3)


## Incoherent ; proton + neutron ; Argon ##
N_1tau_incoh_DUNE_tau_opt_1_1, XSecCon_1tau_incoh_DUNE_tau_opt_1_1 = (N_1tau_p_Ar_DUNE_tau_opt_1_1 + N_1tau_n_Ar_DUNE_tau_opt_1_1), (XSecCon_1tau_p_Ar_DUNE_tau_opt_1_1 + XSecCon_1tau_n_Ar_DUNE_tau_opt_1_1)
N_1tau_incoh_DUNE_tau_opt_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_67_3 = (N_1tau_p_Ar_DUNE_tau_opt_67_3 + N_1tau_n_Ar_DUNE_tau_opt_67_3), (XSecCon_1tau_p_Ar_DUNE_tau_opt_67_3 + XSecCon_1tau_n_Ar_DUNE_tau_opt_67_3)
N_1tau_incoh_DUNE_tau_opt_147_3, XSecCon_1tau_incoh_DUNE_tau_opt_147_3 = (N_1tau_p_Ar_DUNE_tau_opt_147_3 + N_1tau_n_Ar_DUNE_tau_opt_147_3), (XSecCon_1tau_p_Ar_DUNE_tau_opt_147_3 + XSecCon_1tau_n_Ar_DUNE_tau_opt_147_3) 

N_2tau_incoh_DUNE_tau_opt_1_1, XSecCon_2tau_incoh_DUNE_tau_opt_1_1 = (N_2tau_p_Ar_DUNE_tau_opt_1_1 + N_2tau_n_Ar_DUNE_tau_opt_1_1), (XSecCon_2tau_p_Ar_DUNE_tau_opt_1_1 + XSecCon_2tau_n_Ar_DUNE_tau_opt_1_1)
N_2tau_incoh_DUNE_tau_opt_67_3, XSecCon_2tau_incoh_DUNE_tau_opt_67_3 = (N_2tau_p_Ar_DUNE_tau_opt_67_3 + N_2tau_n_Ar_DUNE_tau_opt_67_3), (XSecCon_2tau_p_Ar_DUNE_tau_opt_67_3 + XSecCon_2tau_n_Ar_DUNE_tau_opt_67_3)
N_2tau_incoh_DUNE_tau_opt_147_3, XSecCon_2tau_incoh_DUNE_tau_opt_147_3 = (N_2tau_p_Ar_DUNE_tau_opt_147_3 + N_2tau_n_Ar_DUNE_tau_opt_147_3), (XSecCon_2tau_p_Ar_DUNE_tau_opt_147_3 + XSecCon_2tau_n_Ar_DUNE_tau_opt_147_3) 

N_2mu_incoh_DUNE_tau_opt_1_1, XSecCon_2mu_incoh_DUNE_tau_opt_1_1 = (N_2mu_p_Ar_DUNE_tau_opt_1_1 + N_2mu_n_Ar_DUNE_tau_opt_1_1), (XSecCon_2mu_p_Ar_DUNE_tau_opt_1_1 + XSecCon_2mu_n_Ar_DUNE_tau_opt_1_1)
N_2mu_incoh_DUNE_tau_opt_67_3, XSecCon_2mu_incoh_DUNE_tau_opt_67_3 = (N_2mu_p_Ar_DUNE_tau_opt_67_3 + N_2mu_n_Ar_DUNE_tau_opt_67_3), (XSecCon_2mu_p_Ar_DUNE_tau_opt_67_3 + XSecCon_2mu_n_Ar_DUNE_tau_opt_67_3)
N_2mu_incoh_DUNE_tau_opt_147_3, XSecCon_2mu_incoh_DUNE_tau_opt_147_3 = (N_2mu_p_Ar_DUNE_tau_opt_147_3 + N_2mu_n_Ar_DUNE_tau_opt_147_3), (XSecCon_2mu_p_Ar_DUNE_tau_opt_147_3 + XSecCon_2mu_n_Ar_DUNE_tau_opt_147_3) 

# Upper and Lower limits #
upper_N_1tau_incoh_DUNE_tau_opt_1_1, upper_XSecCon_1tau_incoh_DUNE_tau_opt_1_1 = (upper_N_1tau_p_Ar_DUNE_tau_opt_1_1 + upper_N_1tau_n_Ar_DUNE_tau_opt_1_1), (upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_1_1 + upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_1_1)
upper_N_1tau_incoh_DUNE_tau_opt_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_67_3 = (upper_N_1tau_p_Ar_DUNE_tau_opt_67_3 + upper_N_1tau_n_Ar_DUNE_tau_opt_67_3), (upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_67_3 + upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_67_3)
upper_N_1tau_incoh_DUNE_tau_opt_147_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_147_3 = (upper_N_1tau_p_Ar_DUNE_tau_opt_147_3 + upper_N_1tau_n_Ar_DUNE_tau_opt_147_3), (upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_147_3 + upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_147_3) 

upper_N_2tau_incoh_DUNE_tau_opt_1_1, upper_XSecCon_2tau_incoh_DUNE_tau_opt_1_1 = (upper_N_2tau_p_Ar_DUNE_tau_opt_1_1 + upper_N_2tau_n_Ar_DUNE_tau_opt_1_1), (upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_1_1 + upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_1_1)
upper_N_2tau_incoh_DUNE_tau_opt_67_3, upper_XSecCon_2tau_incoh_DUNE_tau_opt_67_3 = (upper_N_2tau_p_Ar_DUNE_tau_opt_67_3 + upper_N_2tau_n_Ar_DUNE_tau_opt_67_3), (upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_67_3 + upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_67_3)
upper_N_2tau_incoh_DUNE_tau_opt_147_3, upper_XSecCon_2tau_incoh_DUNE_tau_opt_147_3 = (upper_N_2tau_p_Ar_DUNE_tau_opt_147_3 + upper_N_2tau_n_Ar_DUNE_tau_opt_147_3), (upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_147_3 + upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_147_3) 

upper_N_2mu_incoh_DUNE_tau_opt_1_1, upper_XSecCon_2mu_incoh_DUNE_tau_opt_1_1 = (upper_N_2mu_p_Ar_DUNE_tau_opt_1_1 + upper_N_2mu_n_Ar_DUNE_tau_opt_1_1), (upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_1_1 + upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_1_1)
upper_N_2mu_incoh_DUNE_tau_opt_67_3, upper_XSecCon_2mu_incoh_DUNE_tau_opt_67_3 = (upper_N_2mu_p_Ar_DUNE_tau_opt_67_3 + upper_N_2mu_n_Ar_DUNE_tau_opt_67_3), (upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_67_3 + upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_67_3)
upper_N_2mu_incoh_DUNE_tau_opt_147_3, upper_XSecCon_2mu_incoh_DUNE_tau_opt_147_3 = (upper_N_2mu_p_Ar_DUNE_tau_opt_147_3 + upper_N_2mu_n_Ar_DUNE_tau_opt_147_3), (upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_147_3 + upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_147_3) 

lower_N_1tau_incoh_DUNE_tau_opt_1_1, lower_XSecCon_1tau_incoh_DUNE_tau_opt_1_1 = (lower_N_1tau_p_Ar_DUNE_tau_opt_1_1 + lower_N_1tau_n_Ar_DUNE_tau_opt_1_1), (lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_1_1 + lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_1_1)
lower_N_1tau_incoh_DUNE_tau_opt_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_67_3 = (lower_N_1tau_p_Ar_DUNE_tau_opt_67_3 + lower_N_1tau_n_Ar_DUNE_tau_opt_67_3), (lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_67_3 + lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_67_3)
lower_N_1tau_incoh_DUNE_tau_opt_147_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_147_3 = (lower_N_1tau_p_Ar_DUNE_tau_opt_147_3 + lower_N_1tau_n_Ar_DUNE_tau_opt_147_3), (lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_147_3 + lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_147_3) 

lower_N_2tau_incoh_DUNE_tau_opt_1_1, lower_XSecCon_2tau_incoh_DUNE_tau_opt_1_1 = (lower_N_2tau_p_Ar_DUNE_tau_opt_1_1 + lower_N_2tau_n_Ar_DUNE_tau_opt_1_1), (lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_1_1 + lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_1_1)
lower_N_2tau_incoh_DUNE_tau_opt_67_3, lower_XSecCon_2tau_incoh_DUNE_tau_opt_67_3 = (lower_N_2tau_p_Ar_DUNE_tau_opt_67_3 + lower_N_2tau_n_Ar_DUNE_tau_opt_67_3), (lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_67_3 + lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_67_3)
lower_N_2tau_incoh_DUNE_tau_opt_147_3, lower_XSecCon_2tau_incoh_DUNE_tau_opt_147_3 = (lower_N_2tau_p_Ar_DUNE_tau_opt_147_3 + lower_N_2tau_n_Ar_DUNE_tau_opt_147_3), (lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_147_3 + lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_147_3) 

lower_N_2mu_incoh_DUNE_tau_opt_1_1, lower_XSecCon_2mu_incoh_DUNE_tau_opt_1_1 = (lower_N_2mu_p_Ar_DUNE_tau_opt_1_1 + lower_N_2mu_n_Ar_DUNE_tau_opt_1_1), (lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_1_1 + lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_1_1)
lower_N_2mu_incoh_DUNE_tau_opt_67_3, lower_XSecCon_2mu_incoh_DUNE_tau_opt_67_3 = (lower_N_2mu_p_Ar_DUNE_tau_opt_67_3 + lower_N_2mu_n_Ar_DUNE_tau_opt_67_3), (lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_67_3 + lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_67_3)
lower_N_2mu_incoh_DUNE_tau_opt_147_3, lower_XSecCon_2mu_incoh_DUNE_tau_opt_147_3 = (lower_N_2mu_p_Ar_DUNE_tau_opt_147_3 + lower_N_2mu_n_Ar_DUNE_tau_opt_147_3), (lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_147_3 + lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_147_3) 


### FASER vmu Flux Events ###
## Coherent ; Tungsten ##
N_1tau_coh_W_FASERvmu_1_1, XSecCon_1tau_coh_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_coh_W_matched, detector='FASER')
N_1tau_coh_W_FASERvmu_1p2_3, XSecCon_1tau_coh_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_coh_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

N_2tau_coh_W_FASERvmu_1_1, XSecCon_2tau_coh_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_coh_W_matched, detector='FASER')
N_2tau_coh_W_FASERvmu_1p2_3, XSecCon_2tau_coh_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_coh_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

# Upper and Lower limits #
upper_N_1tau_coh_W_FASERvmu_1_1, upper_XSecCon_1tau_coh_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_coh_W_matched_upper, detector='FASER')
upper_N_1tau_coh_W_FASERvmu_1p2_3, upper_XSecCon_1tau_coh_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_coh_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')

upper_N_2tau_coh_W_FASERvmu_1_1, upper_XSecCon_2tau_coh_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_coh_W_matched_upper, detector='FASER')
upper_N_2tau_coh_W_FASERvmu_1p2_3, upper_XSecCon_2tau_coh_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_coh_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')

#upper_N_2mu_coh_W_FASERvmu_1_1, upper_XSecCon_2mu_coh_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_coh_W_matched_upper, detector='FASER')

lower_N_1tau_coh_W_FASERvmu_1_1, lower_XSecCon_1tau_coh_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_coh_W_matched_lower, detector='FASER')
lower_N_1tau_coh_W_FASERvmu_1p2_3, lower_XSecCon_1tau_coh_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_coh_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')

lower_N_2tau_coh_W_FASERvmu_1_1, lower_XSecCon_2tau_coh_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_coh_W_matched_lower, detector='FASER')
lower_N_2tau_coh_W_FASERvmu_1p2_3, lower_XSecCon_2tau_coh_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_coh_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')

#lower_N_2mu_coh_W_FASERvmu_1_1, lower_XSecCon_2mu_coh_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_coh_W_matched_lower, detector='FASER')

## Incoherent ; proton ; Tungsten ##
N_1tau_p_FASERvmu_1_1, XSecCon_1tau_p_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_p_W_matched, detector='FASER')
N_1tau_p_FASERvmu_1p2_3, XSecCon_1tau_p_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_p_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

N_2tau_p_FASERvmu_1_1, XSecCon_2tau_p_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_p_W_matched, detector='FASER')
N_2tau_p_FASERvmu_1p2_3, XSecCon_2tau_p_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_p_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

N_2mu_p_FASERvmu_1_1, XSecCon_2mu_p_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_p_W_matched, detector='FASER')
N_2mu_p_FASERvmu_1p2_3, XSecCon_2mu_p_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_p_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

# Upper and Lower limits #
upper_N_1tau_p_W_FASERvmu_1_1, upper_XSecCon_1tau_p_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_p_W_matched_upper, detector='FASER')
upper_N_1tau_p_W_FASERvmu_1p2_3, upper_XSecCon_1tau_p_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_p_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')

upper_N_2tau_p_W_FASERvmu_1_1, upper_XSecCon_2tau_p_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_p_W_matched_upper, detector='FASER')
upper_N_2tau_p_W_FASERvmu_1p2_3, upper_XSecCon_2tau_p_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_p_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')

upper_N_2mu_p_W_FASERvmu_1_1, upper_XSecCon_2mu_p_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_p_W_matched_upper, detector='FASER')
upper_N_2mu_p_W_FASERvmu_1p2_3, upper_XSecCon_2mu_p_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_p_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')

lower_N_1tau_p_W_FASERvmu_1_1, lower_XSecCon_1tau_p_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_p_W_matched_lower, detector='FASER')
lower_N_1tau_p_W_FASERvmu_1p2_3, lower_XSecCon_1tau_p_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_p_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')

lower_N_2tau_p_W_FASERvmu_1_1, lower_XSecCon_2tau_p_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_p_W_matched_lower, detector='FASER')
lower_N_2tau_p_W_FASERvmu_1p2_3, lower_XSecCon_2tau_p_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_p_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')

lower_N_2mu_p_W_FASERvmu_1_1, lower_XSecCon_2mu_p_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_p_W_matched_lower, detector='FASER')
lower_N_2mu_p_W_FASERvmu_1p2_3, lower_XSecCon_2mu_p_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_p_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')

## Incoherent ; neutron ; Tungsten ##
N_1tau_n_FASERvmu_1_1, XSecCon_1tau_n_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_n_W_matched, detector='FASER')
N_1tau_n_FASERvmu_1p2_3, XSecCon_1tau_n_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_n_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

N_2tau_n_FASERvmu_1_1, XSecCon_2tau_n_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_n_W_matched, detector='FASER')
N_2tau_n_FASERvmu_1p2_3, XSecCon_2tau_n_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_n_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

N_2mu_n_FASERvmu_1_1, XSecCon_2mu_n_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_n_W_matched, detector='FASER')
N_2mu_n_FASERvmu_1p2_3, XSecCon_2mu_n_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_n_W_matched, NTONNES=1.2, NYEAR=3, detector='FASER')

# Upper and Lower limits #
upper_N_1tau_n_W_FASERvmu_1_1, upper_XSecCon_1tau_n_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_n_W_matched_upper, detector='FASER')
upper_N_1tau_n_W_FASERvmu_1p2_3, upper_XSecCon_1tau_n_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_n_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')

upper_N_2tau_n_W_FASERvmu_1_1, upper_XSecCon_2tau_n_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_n_W_matched_upper, detector='FASER')
upper_N_2tau_n_W_FASERvmu_1p2_3, upper_XSecCon_2tau_n_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_n_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')

upper_N_2mu_n_W_FASERvmu_1_1, upper_XSecCon_2mu_n_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_n_W_matched_upper, detector='FASER')
upper_N_2mu_n_W_FASERvmu_1p2_3, upper_XSecCon_2mu_n_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_n_W_matched_upper, NTONNES=1.2, NYEAR=3, detector='FASER')

lower_N_1tau_n_W_FASERvmu_1_1, lower_XSecCon_1tau_n_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_n_W_matched_lower, detector='FASER')
lower_N_1tau_n_W_FASERvmu_1p2_3, lower_XSecCon_1tau_n_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_1tau_n_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')

lower_N_2tau_n_W_FASERvmu_1_1, lower_XSecCon_2tau_n_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_n_W_matched_lower, detector='FASER')
lower_N_2tau_n_W_FASERvmu_1p2_3, lower_XSecCon_2tau_n_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2tau_n_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')

lower_N_2mu_n_W_FASERvmu_1_1, lower_XSecCon_2mu_n_W_FASERvmu_1_1 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_n_W_matched_lower, detector='FASER')
lower_N_2mu_n_W_FASERvmu_1p2_3, lower_XSecCon_2mu_n_W_FASERvmu_1p2_3 = CalculateEvents(flux_FASERvmu, energy_FASERvmu, FASERvmu_xsec_2mu_n_W_matched_lower, NTONNES=1.2, NYEAR=3, detector='FASER')

# Incoherent ; proton + neutron ; Tungsten #
N_1tau_incoh_FASERvmu_1_1, XSecCon_1tau_incoh_FASERvmu_1_1 = (N_1tau_p_FASERvmu_1_1 + N_1tau_n_FASERvmu_1_1), (XSecCon_1tau_p_FASERvmu_1_1 + XSecCon_1tau_n_FASERvmu_1_1)
N_1tau_incoh_FASERvmu_1p2_3, XSecCon_1tau_incoh_FASERvmu_1p2_3 = (N_1tau_p_FASERvmu_1p2_3 + N_1tau_n_FASERvmu_1p2_3), (XSecCon_1tau_p_FASERvmu_1p2_3 + XSecCon_1tau_n_FASERvmu_1p2_3)

N_2tau_incoh_FASERvmu_1_1, XSecCon_2tau_incoh_FASERvmu_1_1 = (N_2tau_p_FASERvmu_1_1 + N_2tau_n_FASERvmu_1_1), (XSecCon_2tau_p_FASERvmu_1_1 + XSecCon_2tau_n_FASERvmu_1_1)
N_2tau_incoh_FASERvmu_1p2_3, XSecCon_2tau_incoh_FASERvmu_1p2_3 = (N_2tau_p_FASERvmu_1p2_3 + N_2tau_n_FASERvmu_1p2_3), (XSecCon_2tau_p_FASERvmu_1p2_3 + XSecCon_2tau_n_FASERvmu_1p2_3)

N_2mu_incoh_FASERvmu_1_1, XSecCon_2mu_incoh_FASERvmu_1_1 = (N_2mu_p_FASERvmu_1_1 + N_2mu_n_FASERvmu_1_1), (XSecCon_2mu_p_FASERvmu_1_1 + XSecCon_2mu_n_FASERvmu_1_1)
N_2mu_incoh_FASERvmu_1p2_3, XSecCon_2mu_incoh_FASERvmu_1p2_3 = (N_2mu_p_FASERvmu_1p2_3 + N_2mu_n_FASERvmu_1p2_3), (XSecCon_2mu_p_FASERvmu_1p2_3 + XSecCon_2mu_n_FASERvmu_1p2_3)

# Upper and Lower limits #
upper_N_1tau_incoh_FASERvmu_1_1, upper_XSecCon_1tau_incoh_FASERvmu_1_1 = (upper_N_1tau_p_W_FASERvmu_1_1 + upper_N_1tau_n_W_FASERvmu_1_1), (upper_XSecCon_1tau_p_W_FASERvmu_1_1 + upper_XSecCon_1tau_n_W_FASERvmu_1_1)
upper_N_1tau_incoh_FASERvmu_1p2_3, upper_XSecCon_1tau_incoh_FASERvmu_1p2_3 = (upper_N_1tau_p_W_FASERvmu_1p2_3 + upper_N_1tau_n_W_FASERvmu_1p2_3), (upper_XSecCon_1tau_p_W_FASERvmu_1p2_3 + upper_XSecCon_1tau_n_W_FASERvmu_1p2_3)

upper_N_2tau_incoh_FASERvmu_1_1, upper_XSecCon_2tau_incoh_FASERvmu_1_1 = (upper_N_2tau_p_W_FASERvmu_1_1 + upper_N_2tau_n_W_FASERvmu_1_1), (upper_XSecCon_2tau_p_W_FASERvmu_1_1 + upper_XSecCon_2tau_n_W_FASERvmu_1_1)
upper_N_2tau_incoh_FASERvmu_1p2_3, upper_XSecCon_2tau_incoh_FASERvmu_1p2_3 = (upper_N_2tau_p_W_FASERvmu_1p2_3 + upper_N_2tau_n_W_FASERvmu_1p2_3), (upper_XSecCon_2tau_p_W_FASERvmu_1p2_3 + upper_XSecCon_2tau_n_W_FASERvmu_1p2_3)

upper_N_2mu_incoh_FASERvmu_1_1, upper_XSecCon_2mu_incoh_FASERvmu_1_1 = (upper_N_2mu_p_W_FASERvmu_1_1 + upper_N_2mu_n_W_FASERvmu_1_1), (upper_XSecCon_2mu_p_W_FASERvmu_1_1 + upper_XSecCon_2mu_n_W_FASERvmu_1_1)
upper_N_2mu_incoh_FASERvmu_1p2_3, upper_XSecCon_2mu_incoh_FASERvmu_1p2_3 = (upper_N_2mu_p_W_FASERvmu_1p2_3 + upper_N_2mu_n_W_FASERvmu_1p2_3), (upper_XSecCon_2mu_p_W_FASERvmu_1p2_3 + upper_XSecCon_2mu_n_W_FASERvmu_1p2_3)

lower_N_1tau_incoh_FASERvmu_1_1, lower_XSecCon_1tau_incoh_FASERvmu_1_1 = (lower_N_1tau_p_W_FASERvmu_1_1 + lower_N_1tau_n_W_FASERvmu_1_1), (lower_XSecCon_1tau_p_W_FASERvmu_1_1 + lower_XSecCon_1tau_n_W_FASERvmu_1_1)
lower_N_1tau_incoh_FASERvmu_1p2_3, lower_XSecCon_1tau_incoh_FASERvmu_1p2_3 = (lower_N_1tau_p_W_FASERvmu_1p2_3 + lower_N_1tau_n_W_FASERvmu_1p2_3), (lower_XSecCon_1tau_p_W_FASERvmu_1p2_3 + lower_XSecCon_1tau_n_W_FASERvmu_1p2_3)

lower_N_2tau_incoh_FASERvmu_1_1, lower_XSecCon_2tau_incoh_FASERvmu_1_1 = (lower_N_2tau_p_W_FASERvmu_1_1 + lower_N_2tau_n_W_FASERvmu_1_1), (lower_XSecCon_2tau_p_W_FASERvmu_1_1 + lower_XSecCon_2tau_n_W_FASERvmu_1_1)
lower_N_2tau_incoh_FASERvmu_1p2_3, lower_XSecCon_2tau_incoh_FASERvmu_1p2_3 = (lower_N_2tau_p_W_FASERvmu_1p2_3 + lower_N_2tau_n_W_FASERvmu_1p2_3), (lower_XSecCon_2tau_p_W_FASERvmu_1p2_3 + lower_XSecCon_2tau_n_W_FASERvmu_1p2_3)

lower_N_2mu_incoh_FASERvmu_1_1, lower_XSecCon_2mu_incoh_FASERvmu_1_1 = (lower_N_2mu_p_W_FASERvmu_1_1 + lower_N_2mu_n_W_FASERvmu_1_1), (lower_XSecCon_2mu_p_W_FASERvmu_1_1 + lower_XSecCon_2mu_n_W_FASERvmu_1_1)
lower_N_2mu_incoh_FASERvmu_1p2_3, lower_XSecCon_2mu_incoh_FASERvmu_1p2_3 = (lower_N_2mu_p_W_FASERvmu_1p2_3 + lower_N_2mu_n_W_FASERvmu_1p2_3), (lower_XSecCon_2mu_p_W_FASERvmu_1p2_3 + lower_XSecCon_2mu_n_W_FASERvmu_1p2_3)

############################
###### Event Printout ######
############################

def PrintOutEvent(filename, fluxname, tonnes, years, events, xsec_con, upper_events=0, lower_events=0, upper_xsec_con=0, lower_xsec_con=0, detector='DUNE'):
    if detector == 'DUNE':
        print("\t{} tonne Ar, {} year, {} flux: {:.3e} + {:.3e} - {:.3e} ; XSecCon: {:.3e} + {:.3e} - {:.3e}".format(tonnes, years, fluxname, events, upper_events-events, events-lower_events, xsec_con, upper_xsec_con-xsec_con, xsec_con-lower_xsec_con), file=filename)
    if detector == 'FASER':
        print("\t{} tonne W, {} fb^-1 run L, {} flux: {:.3e} + {:.3e} - {:.3e} ; XSecCon: {:.3e} + {:.3e} - {:.3e}".format(tonnes, L_Run, fluxname, events, upper_events-events, events-lower_events, xsec_con, upper_xsec_con-xsec_con, xsec_con-lower_xsec_con), file=filename)

def WriteOutFile(filename):
    with open(filename,'w') as textfile:
        print("Events expected at DUNE ND or FASERv for:", file=textfile)
        print("The following integrated fluxes were used:", file=textfile)
        print("\tDUNE Standard Flux: {:.3e} [m^-2 POT^-1 yr^-1]".format(DUNE_integrated_flux), file=textfile)
        print("\tDUNE Tau-optimized Flux: {:.3e} [m^-2 POT^-1 yr^-1]".format(DUNE_tau_opt_integrated_flux), file=textfile)
#        print("\tAlt. Digitized Flux: {:.3e}".format(Alt_integrated_flux), file=textfile)
#        print("\tAlt. 120 GeV Flux: {:.3e}".format(Alt120_integrated_flux), file=textfile)
        print("\tFASER vmu Flux: {:.3e} [fb m^-2]".format(FASER_integrated_flux), file=textfile)
        print("\n", file=textfile)

        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_1tau_coh_Ar_DUNE_1_1, XSecCon_1tau_coh_Ar_DUNE_1_1, upper_N_1tau_coh_Ar_DUNE_1_1, lower_N_1tau_coh_Ar_DUNE_1_1, upper_XSecCon_1tau_coh_Ar_DUNE_1_1, lower_XSecCon_1tau_coh_Ar_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_coh_Ar_DUNE_67_3, XSecCon_1tau_coh_Ar_DUNE_67_3, upper_N_1tau_coh_Ar_DUNE_67_3, lower_N_1tau_coh_Ar_DUNE_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_67_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_1tau_coh_Ar_DUNE_147_3, XSecCon_1tau_coh_Ar_DUNE_147_3, upper_N_1tau_coh_Ar_DUNE_147_3, lower_N_1tau_coh_Ar_DUNE_147_3, upper_XSecCon_1tau_coh_Ar_DUNE_147_3, lower_XSecCon_1tau_coh_Ar_DUNE_147_3)
#        PrintOutEvent(textfile, "Alt120", 1, 1, N_1tau_coh_Ar_Alt120_1_1, XSecCon_1tau_coh_Ar_Alt120_1_1)
#        PrintOutEvent(textfile, "Alt120", 67, 3, N_1tau_coh_Ar_Alt120_67_3, XSecCon_1tau_coh_Ar_Alt120_67_3)
#        PrintOutEvent(textfile, "Alt120", 147, 3, N_1tau_coh_Ar_Alt120_147_3, XSecCon_1tau_coh_Ar_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_1tau_coh_Ar_DUNE_tau_opt_1_1, XSecCon_1tau_coh_Ar_DUNE_tau_opt_1_1, upper_N_1tau_coh_Ar_DUNE_tau_opt_1_1, lower_N_1tau_coh_Ar_DUNE_tau_opt_1_1, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_1_1, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_coh_Ar_DUNE_tau_opt_67_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_67_3, upper_N_1tau_coh_Ar_DUNE_tau_opt_67_3, lower_N_1tau_coh_Ar_DUNE_tau_opt_67_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_67_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_1tau_coh_Ar_DUNE_tau_opt_147_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_147_3, upper_N_1tau_coh_Ar_DUNE_tau_opt_147_3, lower_N_1tau_coh_Ar_DUNE_tau_opt_147_3, upper_XSecCon_1tau_coh_Ar_DUNE_tau_opt_147_3, lower_XSecCon_1tau_coh_Ar_DUNE_tau_opt_147_3)
        PrintOutEvent(textfile, "FASERvmu", 1, 1, N_1tau_coh_W_FASERvmu_1_1, XSecCon_1tau_coh_W_FASERvmu_1_1, upper_N_1tau_coh_W_FASERvmu_1_1, lower_N_1tau_coh_W_FASERvmu_1_1, upper_XSecCon_1tau_coh_W_FASERvmu_1_1, lower_XSecCon_1tau_coh_W_FASERvmu_1_1, detector='FASER')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_coh_W_FASERvmu_1p2_3, XSecCon_1tau_coh_W_FASERvmu_1p2_3, upper_N_1tau_coh_W_FASERvmu_1p2_3, lower_N_1tau_coh_W_FASERvmu_1p2_3, upper_XSecCon_1tau_coh_W_FASERvmu_1p2_3, lower_XSecCon_1tau_coh_W_FASERvmu_1p2_3, detector='FASER')
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_1tau_p_Ar_DUNE_1_1, XSecCon_1tau_p_Ar_DUNE_1_1, upper_N_1tau_p_Ar_DUNE_1_1, lower_N_1tau_p_Ar_DUNE_1_1, upper_XSecCon_1tau_p_Ar_DUNE_1_1, lower_XSecCon_1tau_p_Ar_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_p_Ar_DUNE_67_3, XSecCon_1tau_p_Ar_DUNE_67_3, upper_N_1tau_p_Ar_DUNE_67_3, lower_N_1tau_p_Ar_DUNE_67_3, upper_XSecCon_1tau_p_Ar_DUNE_67_3, lower_XSecCon_1tau_p_Ar_DUNE_67_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_1tau_p_Ar_DUNE_147_3, XSecCon_1tau_p_Ar_DUNE_147_3, upper_N_1tau_p_Ar_DUNE_147_3, lower_N_1tau_p_Ar_DUNE_147_3, upper_XSecCon_1tau_p_Ar_DUNE_147_3, lower_XSecCon_1tau_p_Ar_DUNE_147_3)
#        PrintOutEvent(textfile, "Alt120", 1, 1, N_1tau_p_Alt120_1_1, XSecCon_1tau_p_Alt120_1_1)
#        PrintOutEvent(textfile, "Alt120", 67, 3, N_1tau_p_Alt120_67_3, XSecCon_1tau_p_Alt120_67_3)
#        PrintOutEvent(textfile, "Alt120", 147, 3, N_1tau_p_Alt120_147_3, XSecCon_1tau_p_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_1tau_p_Ar_DUNE_tau_opt_1_1, XSecCon_1tau_p_Ar_DUNE_tau_opt_1_1, upper_N_1tau_p_Ar_DUNE_tau_opt_1_1, lower_N_1tau_p_Ar_DUNE_tau_opt_1_1, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_1_1, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_p_Ar_DUNE_tau_opt_67_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_67_3, upper_N_1tau_p_Ar_DUNE_tau_opt_67_3, lower_N_1tau_p_Ar_DUNE_tau_opt_67_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_67_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_1tau_p_Ar_DUNE_tau_opt_147_3, XSecCon_1tau_p_Ar_DUNE_tau_opt_147_3, upper_N_1tau_p_Ar_DUNE_tau_opt_147_3, lower_N_1tau_p_Ar_DUNE_tau_opt_147_3, upper_XSecCon_1tau_p_Ar_DUNE_tau_opt_147_3, lower_XSecCon_1tau_p_Ar_DUNE_tau_opt_147_3)
        PrintOutEvent(textfile, "FASERvmu", 1, 1, N_1tau_p_FASERvmu_1_1, XSecCon_1tau_p_FASERvmu_1_1, upper_N_1tau_p_W_FASERvmu_1_1, lower_N_1tau_p_W_FASERvmu_1_1, upper_XSecCon_1tau_p_W_FASERvmu_1_1, lower_XSecCon_1tau_p_W_FASERvmu_1_1, detector='FASER')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_p_FASERvmu_1p2_3, XSecCon_1tau_p_FASERvmu_1p2_3, upper_N_1tau_p_W_FASERvmu_1p2_3, lower_N_1tau_p_W_FASERvmu_1p2_3, upper_XSecCon_1tau_p_W_FASERvmu_1p2_3, lower_XSecCon_1tau_p_W_FASERvmu_1p2_3, detector='FASER')
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_1tau_n_Ar_DUNE_1_1, XSecCon_1tau_n_Ar_DUNE_1_1, upper_N_1tau_n_Ar_DUNE_1_1, lower_N_1tau_n_Ar_DUNE_1_1, upper_XSecCon_1tau_n_Ar_DUNE_1_1, lower_XSecCon_1tau_n_Ar_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_n_Ar_DUNE_67_3, XSecCon_1tau_n_Ar_DUNE_67_3, upper_N_1tau_n_Ar_DUNE_67_3, lower_N_1tau_n_Ar_DUNE_67_3, upper_XSecCon_1tau_n_Ar_DUNE_67_3, lower_XSecCon_1tau_n_Ar_DUNE_67_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_1tau_n_Ar_DUNE_147_3, XSecCon_1tau_n_Ar_DUNE_147_3, upper_N_1tau_n_Ar_DUNE_147_3, lower_N_1tau_n_Ar_DUNE_147_3, upper_XSecCon_1tau_n_Ar_DUNE_147_3, lower_XSecCon_1tau_n_Ar_DUNE_147_3)
#        PrintOutEvent(textfile, "Alt120", 1, 1, N_1tau_n_Alt120_1_1, XSecCon_1tau_n_Alt120_1_1)
#        PrintOutEvent(textfile, "Alt120", 67, 3, N_1tau_n_Alt120_67_3, XSecCon_1tau_n_Alt120_67_3)
#        PrintOutEvent(textfile, "Alt120", 147, 3, N_1tau_n_Alt120_147_3, XSecCon_1tau_n_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_1tau_n_Ar_DUNE_tau_opt_1_1, XSecCon_1tau_n_Ar_DUNE_tau_opt_1_1, upper_N_1tau_n_Ar_DUNE_tau_opt_1_1, lower_N_1tau_n_Ar_DUNE_tau_opt_1_1, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_1_1, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_n_Ar_DUNE_tau_opt_67_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_67_3, upper_N_1tau_n_Ar_DUNE_tau_opt_67_3, lower_N_1tau_n_Ar_DUNE_tau_opt_67_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_67_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_1tau_n_Ar_DUNE_tau_opt_147_3, XSecCon_1tau_n_Ar_DUNE_tau_opt_147_3, upper_N_1tau_n_Ar_DUNE_tau_opt_147_3, lower_N_1tau_n_Ar_DUNE_tau_opt_147_3, upper_XSecCon_1tau_n_Ar_DUNE_tau_opt_147_3, lower_XSecCon_1tau_n_Ar_DUNE_tau_opt_147_3)
        PrintOutEvent(textfile, "FASERvmu", 1, 1, N_1tau_n_FASERvmu_1_1, XSecCon_1tau_n_FASERvmu_1_1, upper_N_1tau_n_W_FASERvmu_1_1, lower_N_1tau_n_W_FASERvmu_1_1, upper_XSecCon_1tau_n_W_FASERvmu_1_1, lower_XSecCon_1tau_n_W_FASERvmu_1_1, detector='FASER')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_n_FASERvmu_1p2_3, XSecCon_1tau_n_FASERvmu_1p2_3, upper_N_1tau_n_W_FASERvmu_1p2_3, lower_N_1tau_n_W_FASERvmu_1p2_3, upper_XSecCon_1tau_n_W_FASERvmu_1p2_3, lower_XSecCon_1tau_n_W_FASERvmu_1p2_3, detector='FASER')
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_1tau_incoh_DUNE_1_1, XSecCon_1tau_incoh_DUNE_1_1, upper_N_1tau_incoh_DUNE_1_1, lower_N_1tau_incoh_DUNE_1_1, upper_XSecCon_1tau_incoh_DUNE_1_1, lower_XSecCon_1tau_incoh_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_1tau_incoh_DUNE_67_3, XSecCon_1tau_incoh_DUNE_67_3, upper_N_1tau_incoh_DUNE_67_3, lower_N_1tau_incoh_DUNE_67_3, upper_XSecCon_1tau_incoh_DUNE_67_3, lower_XSecCon_1tau_incoh_DUNE_67_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_1tau_incoh_DUNE_147_3, XSecCon_1tau_incoh_DUNE_147_3, upper_N_1tau_incoh_DUNE_147_3, lower_N_1tau_incoh_DUNE_147_3, upper_XSecCon_1tau_incoh_DUNE_147_3, lower_XSecCon_1tau_incoh_DUNE_147_3)
#        PrintOutEvent(textfile, "Alt120", 1, 1, N_1tau_incoh_Alt120_1_1, XSecCon_1tau_incoh_Alt120_1_1)
#        PrintOutEvent(textfile, "Alt120", 67, 3, N_1tau_incoh_Alt120_67_3, XSecCon_1tau_incoh_Alt120_67_3)
#        PrintOutEvent(textfile, "Alt120", 147, 3, N_1tau_incoh_Alt120_147_3, XSecCon_1tau_incoh_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_1tau_incoh_DUNE_tau_opt_1_1, XSecCon_1tau_incoh_DUNE_tau_opt_1_1, upper_N_1tau_incoh_DUNE_tau_opt_1_1, lower_N_1tau_incoh_DUNE_tau_opt_1_1, upper_XSecCon_1tau_incoh_DUNE_tau_opt_1_1, lower_XSecCon_1tau_incoh_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_1tau_incoh_DUNE_tau_opt_67_3, XSecCon_1tau_incoh_DUNE_tau_opt_67_3, upper_N_1tau_incoh_DUNE_tau_opt_67_3, lower_N_1tau_incoh_DUNE_tau_opt_67_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_67_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_1tau_incoh_DUNE_tau_opt_147_3, XSecCon_1tau_incoh_DUNE_tau_opt_147_3, upper_N_1tau_incoh_DUNE_tau_opt_147_3, lower_N_1tau_incoh_DUNE_tau_opt_147_3, upper_XSecCon_1tau_incoh_DUNE_tau_opt_147_3, lower_XSecCon_1tau_incoh_DUNE_tau_opt_147_3)
        PrintOutEvent(textfile, "FASERvmu", 1, 1, N_1tau_incoh_FASERvmu_1_1, XSecCon_1tau_incoh_FASERvmu_1_1, upper_N_1tau_incoh_FASERvmu_1_1, lower_N_1tau_incoh_FASERvmu_1_1, upper_XSecCon_1tau_incoh_FASERvmu_1_1, lower_XSecCon_1tau_incoh_FASERvmu_1_1, detector='FASER')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_1tau_incoh_FASERvmu_1p2_3, XSecCon_1tau_incoh_FASERvmu_1p2_3, upper_N_1tau_incoh_FASERvmu_1p2_3, lower_N_1tau_incoh_FASERvmu_1p2_3, upper_XSecCon_1tau_incoh_FASERvmu_1p2_3, lower_XSecCon_1tau_incoh_FASERvmu_1p2_3, detector='FASER')
        print("\n", file=textfile)


        print("2tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2tau_coh_Ar_DUNE_1_1, XSecCon_2tau_coh_Ar_DUNE_1_1, upper_N_2tau_coh_Ar_DUNE_1_1, lower_N_2tau_coh_Ar_DUNE_1_1, upper_XSecCon_2tau_coh_Ar_DUNE_1_1, lower_XSecCon_2tau_coh_Ar_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_coh_Ar_DUNE_67_3, XSecCon_2tau_coh_Ar_DUNE_67_3, upper_N_2tau_coh_Ar_DUNE_67_3, lower_N_2tau_coh_Ar_DUNE_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_67_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2tau_coh_Ar_DUNE_147_3, XSecCon_2tau_coh_Ar_DUNE_147_3, upper_N_2tau_coh_Ar_DUNE_147_3, lower_N_2tau_coh_Ar_DUNE_147_3, upper_XSecCon_2tau_coh_Ar_DUNE_147_3, lower_XSecCon_2tau_coh_Ar_DUNE_147_3)
#        PrintOutEvent(textfile, "Alt120", 1, 1, N_2tau_coh_Ar_Alt120_1_1, XSecCon_2tau_coh_Ar_Alt120_1_1)
#        PrintOutEvent(textfile, "Alt120", 67, 3, N_2tau_coh_Ar_Alt120_67_3, XSecCon_2tau_coh_Ar_Alt120_67_3)
#        PrintOutEvent(textfile, "Alt120", 147, 3, N_2tau_coh_Ar_Alt120_147_3, XSecCon_2tau_coh_Ar_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2tau_coh_Ar_DUNE_tau_opt_1_1, XSecCon_2tau_coh_Ar_DUNE_tau_opt_1_1, upper_N_2tau_coh_Ar_DUNE_tau_opt_1_1, lower_N_2tau_coh_Ar_DUNE_tau_opt_1_1, upper_XSecCon_2tau_coh_Ar_DUNE_tau_opt_1_1, lower_XSecCon_2tau_coh_Ar_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_coh_Ar_DUNE_tau_opt_67_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_67_3, upper_N_2tau_coh_Ar_DUNE_tau_opt_67_3, lower_N_2tau_coh_Ar_DUNE_tau_opt_67_3, upper_XSecCon_2tau_coh_Ar_DUNE_tau_opt_67_3, lower_XSecCon_2tau_coh_Ar_DUNE_tau_opt_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2tau_coh_Ar_DUNE_tau_opt_147_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_147_3, upper_N_2tau_coh_Ar_DUNE_tau_opt_147_3, lower_N_2tau_coh_Ar_DUNE_tau_opt_147_3, upper_XSecCon_2tau_coh_Ar_DUNE_tau_opt_147_3, lower_XSecCon_2tau_coh_Ar_DUNE_tau_opt_147_3)
        PrintOutEvent(textfile, "FASERvmu", 1, 1, N_2tau_coh_W_FASERvmu_1_1, XSecCon_2tau_coh_W_FASERvmu_1_1, upper_N_2tau_coh_W_FASERvmu_1_1, lower_N_2tau_coh_W_FASERvmu_1_1, upper_XSecCon_2tau_coh_W_FASERvmu_1_1, lower_XSecCon_2tau_coh_W_FASERvmu_1_1, detector='FASER')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_coh_W_FASERvmu_1p2_3, XSecCon_2tau_coh_W_FASERvmu_1p2_3, upper_N_2tau_coh_W_FASERvmu_1p2_3, lower_N_2tau_coh_W_FASERvmu_1p2_3, upper_XSecCon_2tau_coh_W_FASERvmu_1p2_3, lower_XSecCon_2tau_coh_W_FASERvmu_1p2_3, detector='FASER')
        print("2tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2tau_p_Ar_DUNE_1_1, XSecCon_2tau_p_Ar_DUNE_1_1, upper_N_2tau_p_Ar_DUNE_1_1, lower_N_2tau_p_Ar_DUNE_1_1, upper_XSecCon_2tau_p_Ar_DUNE_1_1, lower_XSecCon_2tau_p_Ar_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_p_Ar_DUNE_67_3, XSecCon_2tau_p_Ar_DUNE_67_3, upper_N_2tau_p_Ar_DUNE_67_3, lower_N_2tau_p_Ar_DUNE_67_3, upper_XSecCon_2tau_p_Ar_DUNE_67_3, lower_XSecCon_2tau_p_Ar_DUNE_67_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2tau_p_Ar_DUNE_147_3, XSecCon_2tau_p_Ar_DUNE_147_3, upper_N_2tau_p_Ar_DUNE_147_3, lower_N_2tau_p_Ar_DUNE_147_3, upper_XSecCon_2tau_p_Ar_DUNE_147_3, lower_XSecCon_2tau_p_Ar_DUNE_147_3)
#        PrintOutEvent(textfile, "Alt120", 1, 1, N_2tau_p_Alt120_1_1, XSecCon_2tau_p_Alt120_1_1)
#        PrintOutEvent(textfile, "Alt120", 67, 3, N_2tau_p_Alt120_67_3, XSecCon_2tau_p_Alt120_67_3)
#        PrintOutEvent(textfile, "Alt120", 147, 3, N_2tau_p_Alt120_147_3, XSecCon_2tau_p_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2tau_p_Ar_DUNE_tau_opt_1_1, XSecCon_2tau_p_Ar_DUNE_tau_opt_1_1, upper_N_2tau_p_Ar_DUNE_tau_opt_1_1, lower_N_2tau_p_Ar_DUNE_tau_opt_1_1, upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_1_1, lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_p_Ar_DUNE_tau_opt_67_3, XSecCon_2tau_p_Ar_DUNE_tau_opt_67_3, upper_N_2tau_p_Ar_DUNE_tau_opt_67_3, lower_N_2tau_p_Ar_DUNE_tau_opt_67_3, upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_67_3, lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2tau_p_Ar_DUNE_tau_opt_147_3, XSecCon_2tau_p_Ar_DUNE_tau_opt_147_3, upper_N_2tau_p_Ar_DUNE_tau_opt_147_3, lower_N_2tau_p_Ar_DUNE_tau_opt_147_3, upper_XSecCon_2tau_p_Ar_DUNE_tau_opt_147_3, lower_XSecCon_2tau_p_Ar_DUNE_tau_opt_147_3)
        PrintOutEvent(textfile, "FASERvmu", 1, 1, N_2tau_p_FASERvmu_1_1, XSecCon_2tau_p_FASERvmu_1_1, upper_N_2tau_p_W_FASERvmu_1_1, lower_N_2tau_p_W_FASERvmu_1_1, upper_XSecCon_2tau_p_W_FASERvmu_1_1, lower_XSecCon_2tau_p_W_FASERvmu_1_1, detector='FASER')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_p_FASERvmu_1p2_3, XSecCon_2tau_p_FASERvmu_1p2_3, upper_N_2tau_p_W_FASERvmu_1p2_3, lower_N_2tau_p_W_FASERvmu_1p2_3, upper_XSecCon_2tau_p_W_FASERvmu_1p2_3, lower_XSecCon_2tau_p_W_FASERvmu_1p2_3, detector='FASER')
        print("2tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2tau_n_Ar_DUNE_1_1, XSecCon_2tau_n_Ar_DUNE_1_1, upper_N_2tau_n_Ar_DUNE_1_1, lower_N_2tau_n_Ar_DUNE_1_1, upper_XSecCon_2tau_n_Ar_DUNE_1_1, lower_XSecCon_2tau_n_Ar_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_n_Ar_DUNE_67_3, XSecCon_2tau_n_Ar_DUNE_67_3, upper_N_2tau_n_Ar_DUNE_67_3, lower_N_2tau_n_Ar_DUNE_67_3, upper_XSecCon_2tau_n_Ar_DUNE_67_3, lower_XSecCon_2tau_n_Ar_DUNE_67_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2tau_n_Ar_DUNE_147_3, XSecCon_2tau_n_Ar_DUNE_147_3, upper_N_2tau_n_Ar_DUNE_147_3, lower_N_2tau_n_Ar_DUNE_147_3, upper_XSecCon_2tau_n_Ar_DUNE_147_3, lower_XSecCon_2tau_n_Ar_DUNE_147_3)
#        PrintOutEvent(textfile, "Alt120", 1, 1, N_2tau_n_Alt120_1_1, XSecCon_2tau_n_Alt120_1_1)
#        PrintOutEvent(textfile, "Alt120", 67, 3, N_2tau_n_Alt120_67_3, XSecCon_2tau_n_Alt120_67_3)
#        PrintOutEvent(textfile, "Alt120", 147, 3, N_2tau_n_Alt120_147_3, XSecCon_2tau_n_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2tau_n_Ar_DUNE_tau_opt_1_1, XSecCon_2tau_n_Ar_DUNE_tau_opt_1_1, upper_N_2tau_n_Ar_DUNE_tau_opt_1_1, lower_N_2tau_n_Ar_DUNE_tau_opt_1_1, upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_1_1, lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_n_Ar_DUNE_tau_opt_67_3, XSecCon_2tau_n_Ar_DUNE_tau_opt_67_3, upper_N_2tau_n_Ar_DUNE_tau_opt_67_3, lower_N_2tau_n_Ar_DUNE_tau_opt_67_3, upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_67_3, lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2tau_n_Ar_DUNE_tau_opt_147_3, XSecCon_2tau_n_Ar_DUNE_tau_opt_147_3, upper_N_2tau_n_Ar_DUNE_tau_opt_147_3, lower_N_2tau_n_Ar_DUNE_tau_opt_147_3, upper_XSecCon_2tau_n_Ar_DUNE_tau_opt_147_3, lower_XSecCon_2tau_n_Ar_DUNE_tau_opt_147_3)
        PrintOutEvent(textfile, "FASERvmu", 1, 1, N_2tau_n_FASERvmu_1_1, XSecCon_2tau_n_FASERvmu_1_1, upper_N_2tau_n_W_FASERvmu_1_1, lower_N_2tau_n_W_FASERvmu_1_1, upper_XSecCon_2tau_n_W_FASERvmu_1_1, lower_XSecCon_2tau_n_W_FASERvmu_1_1, detector='FASER')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_n_FASERvmu_1p2_3, XSecCon_2tau_n_FASERvmu_1p2_3, upper_N_2tau_n_W_FASERvmu_1p2_3, lower_N_2tau_n_W_FASERvmu_1p2_3, upper_XSecCon_2tau_n_W_FASERvmu_1p2_3, lower_XSecCon_2tau_n_W_FASERvmu_1p2_3, detector='FASER')
        print("2tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2tau_incoh_DUNE_1_1, XSecCon_2tau_incoh_DUNE_1_1, upper_N_2tau_incoh_DUNE_1_1, lower_N_2tau_incoh_DUNE_1_1, upper_XSecCon_2tau_incoh_DUNE_1_1, lower_XSecCon_2tau_incoh_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2tau_incoh_DUNE_67_3, XSecCon_2tau_incoh_DUNE_67_3, upper_N_2tau_incoh_DUNE_67_3, lower_N_2tau_incoh_DUNE_67_3, upper_XSecCon_2tau_incoh_DUNE_67_3, lower_XSecCon_2tau_incoh_DUNE_67_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2tau_incoh_DUNE_147_3, XSecCon_2tau_incoh_DUNE_147_3, upper_N_2tau_incoh_DUNE_147_3, lower_N_2tau_incoh_DUNE_147_3, upper_XSecCon_2tau_incoh_DUNE_147_3, lower_XSecCon_2tau_incoh_DUNE_147_3)
#        PrintOutEvent(textfile, "Alt120", 1, 1, N_2tau_incoh_Alt120_1_1, XSecCon_2tau_incoh_Alt120_1_1)
#        PrintOutEvent(textfile, "Alt120", 67, 3, N_2tau_incoh_Alt120_67_3, XSecCon_2tau_incoh_Alt120_67_3)
#        PrintOutEvent(textfile, "Alt120", 147, 3, N_2tau_incoh_Alt120_147_3, XSecCon_2tau_incoh_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2tau_incoh_DUNE_tau_opt_1_1, XSecCon_2tau_incoh_DUNE_tau_opt_1_1, upper_N_2tau_incoh_DUNE_tau_opt_1_1, lower_N_2tau_incoh_DUNE_tau_opt_1_1, upper_XSecCon_2tau_incoh_DUNE_tau_opt_1_1, lower_XSecCon_2tau_incoh_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2tau_incoh_DUNE_tau_opt_67_3, XSecCon_2tau_incoh_DUNE_tau_opt_67_3, upper_N_2tau_incoh_DUNE_tau_opt_67_3, lower_N_2tau_incoh_DUNE_tau_opt_67_3, upper_XSecCon_2tau_incoh_DUNE_tau_opt_67_3, lower_XSecCon_2tau_incoh_DUNE_tau_opt_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2tau_incoh_DUNE_tau_opt_147_3, XSecCon_2tau_incoh_DUNE_tau_opt_147_3, upper_N_2tau_incoh_DUNE_tau_opt_147_3, lower_N_2tau_incoh_DUNE_tau_opt_147_3, upper_XSecCon_2tau_incoh_DUNE_tau_opt_147_3, lower_XSecCon_2tau_incoh_DUNE_tau_opt_147_3)
        PrintOutEvent(textfile, "FASERvmu", 1, 1, N_2tau_incoh_FASERvmu_1_1, XSecCon_2tau_incoh_FASERvmu_1_1, upper_N_2tau_incoh_FASERvmu_1_1, lower_N_2tau_incoh_FASERvmu_1_1, upper_XSecCon_2tau_incoh_FASERvmu_1_1, lower_XSecCon_2tau_incoh_FASERvmu_1_1, detector='FASER')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2tau_incoh_FASERvmu_1p2_3, XSecCon_2tau_incoh_FASERvmu_1p2_3, upper_N_2tau_incoh_FASERvmu_1p2_3, lower_N_2tau_incoh_FASERvmu_1p2_3, upper_XSecCon_2tau_incoh_FASERvmu_1p2_3, lower_XSecCon_2tau_incoh_FASERvmu_1p2_3, detector='FASER')
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2mu_coh_Ar_DUNE_1_1, XSecCon_2mu_coh_Ar_DUNE_1_1, upper_N_2mu_coh_Ar_DUNE_1_1, lower_N_2mu_coh_Ar_DUNE_1_1, upper_XSecCon_2mu_coh_Ar_DUNE_1_1, lower_XSecCon_2mu_coh_Ar_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_coh_Ar_DUNE_67_3, XSecCon_2mu_coh_Ar_DUNE_67_3, upper_N_2mu_coh_Ar_DUNE_67_3, lower_N_2mu_coh_Ar_DUNE_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_67_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2mu_coh_Ar_DUNE_147_3, XSecCon_2mu_coh_Ar_DUNE_147_3, upper_N_2mu_coh_Ar_DUNE_147_3, lower_N_2mu_coh_Ar_DUNE_147_3, upper_XSecCon_2mu_coh_Ar_DUNE_147_3, lower_XSecCon_2mu_coh_Ar_DUNE_147_3)
#        PrintOutEvent(textfile, "Alt120", 1, 1, N_2mu_coh_Ar_Alt120_1_1, XSecCon_2mu_coh_Ar_Alt120_1_1)
#        PrintOutEvent(textfile, "Alt120", 67, 3, N_2mu_coh_Ar_Alt120_67_3, XSecCon_2mu_coh_Ar_Alt120_67_3)
#        PrintOutEvent(textfile, "Alt120", 147, 3, N_2mu_coh_Ar_Alt120_147_3, XSecCon_2mu_coh_Ar_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2mu_coh_Ar_DUNE_tau_opt_1_1, XSecCon_2mu_coh_Ar_DUNE_tau_opt_1_1, upper_N_2mu_coh_Ar_DUNE_tau_opt_1_1, lower_N_2mu_coh_Ar_DUNE_tau_opt_1_1, upper_XSecCon_2mu_coh_Ar_DUNE_tau_opt_1_1, lower_XSecCon_2mu_coh_Ar_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_coh_Ar_DUNE_tau_opt_67_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_67_3, upper_N_2mu_coh_Ar_DUNE_tau_opt_67_3, lower_N_2mu_coh_Ar_DUNE_tau_opt_67_3, upper_XSecCon_2mu_coh_Ar_DUNE_tau_opt_67_3, lower_XSecCon_2mu_coh_Ar_DUNE_tau_opt_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2mu_coh_Ar_DUNE_tau_opt_147_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_147_3, upper_N_2mu_coh_Ar_DUNE_tau_opt_147_3, lower_N_2mu_coh_Ar_DUNE_tau_opt_147_3, upper_XSecCon_2mu_coh_Ar_DUNE_tau_opt_147_3, lower_XSecCon_2mu_coh_Ar_DUNE_tau_opt_147_3)
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2mu_p_Ar_DUNE_1_1, XSecCon_2mu_p_Ar_DUNE_1_1, upper_N_2mu_p_Ar_DUNE_1_1, lower_N_2mu_p_Ar_DUNE_1_1, upper_XSecCon_2mu_p_Ar_DUNE_1_1, lower_XSecCon_2mu_p_Ar_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_p_Ar_DUNE_67_3, XSecCon_2mu_p_Ar_DUNE_67_3, upper_N_2mu_p_Ar_DUNE_67_3, lower_N_2mu_p_Ar_DUNE_67_3, upper_XSecCon_2mu_p_Ar_DUNE_67_3, lower_XSecCon_2mu_p_Ar_DUNE_67_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2mu_p_Ar_DUNE_147_3, XSecCon_2mu_p_Ar_DUNE_147_3, upper_N_2mu_p_Ar_DUNE_147_3, lower_N_2mu_p_Ar_DUNE_147_3, upper_XSecCon_2mu_p_Ar_DUNE_147_3, lower_XSecCon_2mu_p_Ar_DUNE_147_3)
#        PrintOutEvent(textfile, "Alt120", 1, 1, N_2mu_p_Alt120_1_1, XSecCon_2mu_p_Alt120_1_1)
#        PrintOutEvent(textfile, "Alt120", 67, 3, N_2mu_p_Alt120_67_3, XSecCon_2mu_p_Alt120_67_3)
#        PrintOutEvent(textfile, "Alt120", 147, 3, N_2mu_p_Alt120_147_3, XSecCon_2mu_p_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2mu_p_Ar_DUNE_tau_opt_1_1, XSecCon_2mu_p_Ar_DUNE_tau_opt_1_1, upper_N_2mu_p_Ar_DUNE_tau_opt_1_1, lower_N_2mu_p_Ar_DUNE_tau_opt_1_1, upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_1_1, lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_p_Ar_DUNE_tau_opt_67_3, XSecCon_2mu_p_Ar_DUNE_tau_opt_67_3, upper_N_2mu_p_Ar_DUNE_tau_opt_67_3, lower_N_2mu_p_Ar_DUNE_tau_opt_67_3, upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_67_3, lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2mu_p_Ar_DUNE_tau_opt_147_3, XSecCon_2mu_p_Ar_DUNE_tau_opt_147_3, upper_N_2mu_p_Ar_DUNE_tau_opt_147_3, lower_N_2mu_p_Ar_DUNE_tau_opt_147_3, upper_XSecCon_2mu_p_Ar_DUNE_tau_opt_147_3, lower_XSecCon_2mu_p_Ar_DUNE_tau_opt_147_3)
        PrintOutEvent(textfile, "FASERvmu", 1, 1, N_2mu_p_FASERvmu_1_1, XSecCon_2mu_p_FASERvmu_1_1, upper_N_2mu_p_W_FASERvmu_1_1, lower_N_2mu_p_W_FASERvmu_1_1, upper_XSecCon_2mu_p_W_FASERvmu_1_1, lower_XSecCon_2mu_p_W_FASERvmu_1_1, detector='FASER')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_p_FASERvmu_1p2_3, XSecCon_2mu_p_FASERvmu_1p2_3, upper_N_2mu_p_W_FASERvmu_1p2_3, lower_N_2mu_p_W_FASERvmu_1p2_3, upper_XSecCon_2mu_p_W_FASERvmu_1p2_3, lower_XSecCon_2mu_p_W_FASERvmu_1p2_3, detector='FASER')
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2mu_n_Ar_DUNE_1_1, XSecCon_2mu_n_Ar_DUNE_1_1, upper_N_2mu_n_Ar_DUNE_1_1, lower_N_2mu_n_Ar_DUNE_1_1, upper_XSecCon_2mu_n_Ar_DUNE_1_1, lower_XSecCon_2mu_n_Ar_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_n_Ar_DUNE_67_3, XSecCon_2mu_n_Ar_DUNE_67_3, upper_N_2mu_n_Ar_DUNE_67_3, lower_N_2mu_n_Ar_DUNE_67_3, upper_XSecCon_2mu_n_Ar_DUNE_67_3, lower_XSecCon_2mu_n_Ar_DUNE_67_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2mu_n_Ar_DUNE_147_3, XSecCon_2mu_n_Ar_DUNE_147_3, upper_N_2mu_n_Ar_DUNE_147_3, lower_N_2mu_n_Ar_DUNE_147_3, upper_XSecCon_2mu_n_Ar_DUNE_147_3, lower_XSecCon_2mu_n_Ar_DUNE_147_3)
#        PrintOutEvent(textfile, "Alt120", 1, 1, N_2mu_n_Alt120_1_1, XSecCon_2mu_n_Alt120_1_1)
#        PrintOutEvent(textfile, "Alt120", 67, 3, N_2mu_n_Alt120_67_3, XSecCon_2mu_n_Alt120_67_3)
#        PrintOutEvent(textfile, "Alt120", 147, 3, N_2mu_n_Alt120_147_3, XSecCon_2mu_n_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2mu_n_Ar_DUNE_tau_opt_1_1, XSecCon_2mu_n_Ar_DUNE_tau_opt_1_1, upper_N_2mu_n_Ar_DUNE_tau_opt_1_1, lower_N_2mu_n_Ar_DUNE_tau_opt_1_1, upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_1_1, lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_n_Ar_DUNE_tau_opt_67_3, XSecCon_2mu_n_Ar_DUNE_tau_opt_67_3, upper_N_2mu_n_Ar_DUNE_tau_opt_67_3, lower_N_2mu_n_Ar_DUNE_tau_opt_67_3, upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_67_3, lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2mu_n_Ar_DUNE_tau_opt_147_3, XSecCon_2mu_n_Ar_DUNE_tau_opt_147_3, upper_N_2mu_n_Ar_DUNE_tau_opt_147_3, lower_N_2mu_n_Ar_DUNE_tau_opt_147_3, upper_XSecCon_2mu_n_Ar_DUNE_tau_opt_147_3, lower_XSecCon_2mu_n_Ar_DUNE_tau_opt_147_3)
        PrintOutEvent(textfile, "FASERvmu", 1, 1, N_2mu_n_FASERvmu_1_1, XSecCon_2mu_n_FASERvmu_1_1, upper_N_2mu_n_W_FASERvmu_1_1, lower_N_2mu_n_W_FASERvmu_1_1, upper_XSecCon_2mu_n_W_FASERvmu_1_1, lower_XSecCon_2mu_n_W_FASERvmu_1_1, detector='FASER')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_n_FASERvmu_1p2_3, XSecCon_2mu_n_FASERvmu_1p2_3, upper_N_2mu_n_W_FASERvmu_1p2_3, lower_N_2mu_n_W_FASERvmu_1p2_3, upper_XSecCon_2mu_n_W_FASERvmu_1p2_3, lower_XSecCon_2mu_n_W_FASERvmu_1p2_3, detector='FASER')
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2mu_incoh_DUNE_1_1, XSecCon_2mu_incoh_DUNE_1_1, upper_N_2mu_incoh_DUNE_1_1, lower_N_2mu_incoh_DUNE_1_1, upper_XSecCon_2mu_incoh_DUNE_1_1, lower_XSecCon_2mu_incoh_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 67, 3, N_2mu_incoh_DUNE_67_3, XSecCon_2mu_incoh_DUNE_67_3, upper_N_2mu_incoh_DUNE_67_3, lower_N_2mu_incoh_DUNE_67_3, upper_XSecCon_2mu_incoh_DUNE_67_3, lower_XSecCon_2mu_incoh_DUNE_67_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2mu_incoh_DUNE_147_3, XSecCon_2mu_incoh_DUNE_147_3, upper_N_2mu_incoh_DUNE_147_3, lower_N_2mu_incoh_DUNE_147_3, upper_XSecCon_2mu_incoh_DUNE_147_3, lower_XSecCon_2mu_incoh_DUNE_147_3)
#        PrintOutEvent(textfile, "Alt120", 1, 1, N_2mu_incoh_Alt120_1_1, XSecCon_2mu_incoh_Alt120_1_1)
#        PrintOutEvent(textfile, "Alt120", 67, 3, N_2mu_incoh_Alt120_67_3, XSecCon_2mu_incoh_Alt120_67_3)
#        PrintOutEvent(textfile, "Alt120", 147, 3, N_2mu_incoh_Alt120_147_3, XSecCon_2mu_incoh_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2mu_incoh_DUNE_tau_opt_1_1, XSecCon_2mu_incoh_DUNE_tau_opt_1_1, upper_N_2mu_incoh_DUNE_tau_opt_1_1, lower_N_2mu_incoh_DUNE_tau_opt_1_1, upper_XSecCon_2mu_incoh_DUNE_tau_opt_1_1, lower_XSecCon_2mu_incoh_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 67, 3, N_2mu_incoh_DUNE_tau_opt_67_3, XSecCon_2mu_incoh_DUNE_tau_opt_67_3, upper_N_2mu_incoh_DUNE_tau_opt_67_3, lower_N_2mu_incoh_DUNE_tau_opt_67_3, upper_XSecCon_2mu_incoh_DUNE_tau_opt_67_3, lower_XSecCon_2mu_incoh_DUNE_tau_opt_67_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2mu_incoh_DUNE_tau_opt_147_3, XSecCon_2mu_incoh_DUNE_tau_opt_147_3, upper_N_2mu_incoh_DUNE_tau_opt_147_3, lower_N_2mu_incoh_DUNE_tau_opt_147_3, upper_XSecCon_2mu_incoh_DUNE_tau_opt_147_3, lower_XSecCon_2mu_incoh_DUNE_tau_opt_147_3)
        PrintOutEvent(textfile, "FASERvmu", 1, 1, N_2mu_incoh_FASERvmu_1_1, XSecCon_2mu_incoh_FASERvmu_1_1, upper_N_2mu_incoh_FASERvmu_1_1, lower_N_2mu_incoh_FASERvmu_1_1, upper_XSecCon_2mu_incoh_FASERvmu_1_1, lower_XSecCon_2mu_incoh_FASERvmu_1_1, detector='FASER')
        PrintOutEvent(textfile, "FASERvmu", 1.2, 3, N_2mu_incoh_FASERvmu_1p2_3, XSecCon_2mu_incoh_FASERvmu_1p2_3, upper_N_2mu_incoh_FASERvmu_1p2_3, lower_N_2mu_incoh_FASERvmu_1p2_3, upper_XSecCon_2mu_incoh_FASERvmu_1p2_3, lower_XSecCon_2mu_incoh_FASERvmu_1p2_3, detector='FASER')
        print("\n", file=textfile)


WriteOutFile("number_of_events.txt")

################################
###### Debugging Plotting ######
################################

#print("DUNE integrated flux: ", DUNE_integrated_flux)
#print("DUNE tau-opt integrated flux: ", DUNE_tau_opt_integrated_flux)
#print("Altmannshofer integrated flux: ", Phi_Alt)
#
#print("Length of DUNE energy flux array: ", len(energy_DUNE))
#print("Length of DUNE 2mu coh Ar matched xsec array: ", len(DUNE_xsec_2mu_coh_Ar_matched))
#print("Length of DUNE 2mu p matched xsec array: ", len(DUNE_xsec_2mu_p_matched))
#
#print("Length of Alt. digitized energy flux array: ", len(energy_Alt))
#print("Length of Alt. digitized 2mu coh Ar matched xsec array: ", len(Alt_xsec_2mu_coh_Ar_matched))
#print("Length of Alt. digitized 2mu p matched xsec array: ", len(Alt_xsec_2mu_p_matched))
#
#print("Length of Alt. 120 energy flux array: ", len(energy_Alt120))
#print("Length of Alt. 120 2mu coh Ar matched xsec array: ", len(Alt120_xsec_2mu_coh_Ar_matched))
#print("Length of Alt. 120 2mu p matched xsec array: ", len(Alt120_xsec_2mu_p_matched))
#
#print("Length of DUNE tau-opt energy flux array: ", len(energy_DUNE_tau_opt))
#print("Length of DUNE tau-opt 2mu coh Ar matched xsec array: ", len(DUNE_tau_opt_xsec_2mu_coh_Ar_matched))
#print("Length of DUNE tau-opt 2mu p matched xsec array: ", len(DUNE_tau_opt_xsec_2mu_p_matched))
#
#flux_DUNE_norm = np.divide(flux_DUNE, DUNE_integrated_flux) # Normalize DUNE flux
#flux_DUNE_tau_opt_norm = np.divide(flux_DUNE_tau_opt, DUNE_tau_opt_integrated_flux) # Normalize DUNE tau-opt flux
#
#### Plot normalized fluxes ###
#fig1, ax1 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)
#
##ax1.hist(energy_DUNE, bins=bins_DUNE, weights=flux_DUNE_norm, histtype='stepfilled', label=r'DUNE', color='navy', alpha=0.5, lw=0.5)
##ax1.hist(energy_DUNE, bins=bins_DUNE, weights=flux_DUNE_norm, histtype='step', color='black', lw=2, alpha=1)
#
#ax1.hist(energy_Alt120, bins=bins_Alt120, weights=flux_Alt120, histtype='stepfilled', label=r'Alt. 120 GeV', color='orange', alpha=0.3, lw=0.5)
#ax1.hist(energy_Alt120, bins=bins_Alt120, weights=flux_Alt120, histtype='step', color='black',lw=2, alpha=1)
#
#ax1.plot(energy_DUNE, flux_DUNE_norm, '-', color='navy', label='DUNE', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
#ax1.plot(energy_Alt, flux_Alt, '-', color='orange', label='Alt. digitized', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
#ax1.plot(energy_DUNE_tau_opt, flux_DUNE_tau_opt_norm, '-', color='firebrick', label=r'DUNE $\tau-$opt', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
#
#ax1.set_xlabel('Energy (GeV)')
#ax1.set_ylabel(r'$\frac{1}{\Phi}\frac{\mathrm{d}\Phi}{\mathrm{d}E}$ (GeV$^{-1}$)')
#ax1.set_xscale('log')
#ax1.set_yscale('log')
#ax1.set_xlim(0.3, 100)
##ax1.set_ylim(5e-4, 0.500)
#ax1.set_ylim(1e-6, 0.500)
#ax1.legend(loc='upper right')
#ax1.set_title(r"$\nu_\mu$ Normalized Flux")
#
#fig1.savefig("../plots/fluxes.png", dpi=400)
#
#### Plot fluxes with cross sections ###
#fig2, ax2 = plt.subplots(2, 3, figsize=(45, 20), sharex=True, tight_layout=True)
#
#col_map = mpl.colormaps['Dark2'].resampled(4)
#rgb_col = np.linspace(0,1,num=4)
#c1 = col_map(rgb_col[0])
#c2 = col_map(rgb_col[1])
#c3 = col_map(rgb_col[2])
#c4 = col_map(rgb_col[3])
#
### Data ##
## 2mu Coherent ; Argon #
#mu2_coh_DUNE_integrand = np.multiply(flux_DUNE_norm, DUNE_xsec_2mu_coh_Ar_matched)
#mu2_coh_DUNE_tau_opt_integrand = np.multiply(flux_DUNE_tau_opt_norm, DUNE_tau_opt_xsec_2mu_coh_Ar_matched)
#
## 2mu Incoherent ; proton + neutron ; Argon #
#mu2_incoh_DUNE_integrand = np.multiply(flux_DUNE_norm, [sum(xsecs) for xsecs in zip(DUNE_xsec_2mu_p_matched, DUNE_xsec_2mu_n_matched)])
#mu2_incoh_DUNE_tau_opt_integrand = np.multiply(flux_DUNE_tau_opt_norm, [sum(xsecs) for xsecs in zip(DUNE_tau_opt_xsec_2mu_p_matched, DUNE_tau_opt_xsec_2mu_n_matched)])
#
## 1tau Coherent ; Argon #
#tau1_coh_DUNE_integrand = np.multiply(flux_DUNE_norm, DUNE_xsec_1tau_coh_Ar_matched)
#tau1_coh_DUNE_tau_opt_integrand = np.multiply(flux_DUNE_tau_opt_norm, DUNE_tau_opt_xsec_1tau_coh_Ar_matched)
#
## 1tau Incoherent ; proton + neutron ; Argon #
#tau1_incoh_DUNE_integrand = np.multiply(flux_DUNE_norm, [sum(xsecs) for xsecs in zip(DUNE_xsec_1tau_p_matched, DUNE_xsec_1tau_n_matched)])
#tau1_incoh_DUNE_tau_opt_integrand = np.multiply(flux_DUNE_tau_opt_norm, [sum(xsecs) for xsecs in zip(DUNE_tau_opt_xsec_1tau_p_matched, DUNE_tau_opt_xsec_1tau_n_matched)])
#
## 2tau Coherent ; Argon #
#tau2_coh_DUNE_integrand = np.multiply(flux_DUNE_norm, DUNE_xsec_2tau_coh_Ar_matched)
#tau2_coh_DUNE_tau_opt_integrand = np.multiply(flux_DUNE_tau_opt_norm, DUNE_tau_opt_xsec_2tau_coh_Ar_matched)
#
## 2tau Incoherent ; proton + neutron ; Argon #
#tau2_incoh_DUNE_integrand = np.multiply(flux_DUNE_norm, [sum(xsecs) for xsecs in zip(DUNE_xsec_2tau_p_matched, DUNE_xsec_2tau_n_matched)])
#tau2_incoh_DUNE_tau_opt_integrand = np.multiply(flux_DUNE_tau_opt_norm, [sum(xsecs) for xsecs in zip(DUNE_tau_opt_xsec_2tau_p_matched, DUNE_tau_opt_xsec_2tau_n_matched)])
#
### Plotting ##
## 2mu Coherent ; Argon #
#ax2[0,0].plot(energy_DUNE, mu2_coh_DUNE_integrand, c = c1, label = 'DUNE', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[0,0].plot(energy_DUNE_tau_opt, mu2_coh_DUNE_tau_opt_integrand, c = c2, label = r'DUNE $\tau-$opt.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[0,0].set_title(r'$\nu_\mu\mu^+\mu^-$ Coh.')
#
## 2mu Incoherent ; Argon #
#ax2[1,0].plot(energy_DUNE, mu2_incoh_DUNE_integrand, c = c3, label = 'DUNE', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[1,0].plot(energy_DUNE_tau_opt, mu2_incoh_DUNE_tau_opt_integrand, c = c4, label = r'DUNE $\tau-$opt.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[1,0].set_title(r'$\nu_\mu\mu^+\mu^-$ Incoh.')
#
## 1tau Coherent ; Argon #
#ax2[0,1].plot(energy_DUNE, tau1_coh_DUNE_integrand, c = c1, label = 'DUNE', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[0,1].plot(energy_DUNE_tau_opt, tau1_coh_DUNE_tau_opt_integrand, c = c2, label = r'DUNE $\tau-$opt.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[0,1].set_title(r'$\nu_\tau\tau^+\mu^-$ Coh.')
#
## 1tau Incoherent ; Argon #
#ax2[1,1].plot(energy_DUNE, tau1_incoh_DUNE_integrand, c = c3, label = 'DUNE', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[1,1].plot(energy_DUNE_tau_opt, tau1_incoh_DUNE_tau_opt_integrand, c = c4, label = r'DUNE $\tau-$opt.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[1,1].set_title(r'$\nu_\tau\tau^+\mu^-$ Incoh.')
#
## 2tau Coherent ; Argon #
#ax2[0,2].plot(energy_DUNE, tau2_coh_DUNE_integrand, c = c1, label = 'DUNE', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[0,2].plot(energy_DUNE_tau_opt, tau2_coh_DUNE_tau_opt_integrand, c = c2, label = r'DUNE $\tau-$opt.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[0,2].set_title(r'$\nu_\mu\tau^+\tau^-$ Coh.')
#
## 2tau Incoherent ; Argon #
#ax2[1,2].plot(energy_DUNE, tau2_incoh_DUNE_integrand, c = c3, label = 'DUNE', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[1,2].plot(energy_DUNE_tau_opt, tau2_incoh_DUNE_tau_opt_integrand, c = c4, label = r'DUNE $\tau-$opt.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[1,2].set_title(r'$\nu_\mu\tau^+\tau^-$ Incoh.')
#
#for ax1D in ax2:   # ax2 is a 2D array so ax will be a 1D array
#    for ax in ax1D:
#        ax.set_xlabel('Energy (GeV)')
#        ax.set_ylabel(r'Norm. Flux $\times$ Cross Section')
#        ax.set_xscale('log')
#        ax.set_yscale('log')
#        ax.grid()
#        ax.legend()
#
#fig2.savefig("eventCalc_integrands.png", dpi=400)
#fig2, ax2 = plt.subplots(2, 1, figsize=(15,20), tight_layout=True)
#
## 2mu Coherent ; Argon #
#ax2[0].plot(energy_Alt120, Alt120_xsec_2mu_coh_Ar_matched, '--', color='firebrick', label='Alt. 120 Gev', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[0].plot(energy_Alt, Alt_xsec_2mu_coh_Ar_matched, '--', color='c', label='Alt. digitized', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[0].plot(energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched, '--', color='goldenrod', label='DUNE', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[0].plot(energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched, '--', color='green', label=r'DUNE $\nu_\tau$', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[0].plot(energy_2mu_coh_Ar, xsec_2mu_coh_Ar, '--', color='navy', label='MC', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[0].set_title('Argon', fontsize=40)
##ax2[0].set_xlim(0.8, 25)
#ax2[0].set_xlim(0.8, 2)
#
## 2mu Nucleon ; proton #
#ax2[1].plot(energy_Alt120, Alt120_xsec_2mu_p_matched, '--', color='firebrick', label='Alt. 120 Gev', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[1].plot(energy_Alt, Alt_xsec_2mu_p_matched, '--', color='c', label='Alt. digitized', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[1].plot(energy_DUNE, DUNE_xsec_2mu_p_matched, '--', color='goldenrod', label='DUNE', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[1].plot(energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_matched, '--', color='green', label=r'DUNE $\nu_\tau$', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[1].plot(energy_2mu_p, xsec_2mu_p, '--', color='navy', label='MC', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2[1].set_title('Proton', fontsize=40)
##ax2[1].set_xlim(0.05, 25)
#ax2[1].set_xlim(0.05, 1)
#
#for i in [0,1]:
#    ax2[i].set_xlabel('Energy (GeV)')
#    ax2[i].set_ylabel(r'Cross Section (m$^2$)')
#    ax2[i].set_xscale('log')
#    ax2[i].set_yscale('log')
#    ax2[i].grid()
#    ax2[i].legend(loc='upper left', fontsize=25)
#
#fig2.suptitle(r'$\nu_\mu\to\nu_\mu \mu^+ \mu^-$', fontsize=50)
#fig2.savefig("eventCalc_xsec.png", dpi=400)
