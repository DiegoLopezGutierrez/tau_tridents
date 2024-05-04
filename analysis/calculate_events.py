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

N_POT = 1.1e21               # Number of POTs
MD = 1000                    # Mass of the DUNE ND detector in kg (1 tonne)
MAr = 39.95*atomic_mass_unit # Mass of argon in atomic mass units
A = 40
Z = 18

Phi_Alt = 1.04e-3

########################
###### Initialize ######
########################

### Directories ###
FLUX_DIR = '../csv/fluxes'
CROSS_SECTION_DIR = '../csv/cross_sections'

### Fluxes ###
#DUNE_filename = FLUX_DIR + '/DUNE/DUNE_diff_flux.csv'  # dPhi/dE
#DUNE_filename = FLUX_DIR + '/DUNE/DUNE_hist.csv'  # dPhi/dE
DUNE_filename = FLUX_DIR + '/DUNE/histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_LBNEND_globes_flux.txt'
Alt_filename = FLUX_DIR + '/Altmannshofer/vmu_normalized_flux_Altmannshofer_digitized.csv' # 1/Phi * dPhi/dE
DUNE_tau_opt_numu_flux = FLUX_DIR + '/DUNE/histos_g4lbne_v3r5p4_QGSP_BERT_TauOptimized_neutrino_LBNEND_globes_flux.txt'
Alt120_filename = FLUX_DIR + '/Altmannshofer/numu_flux_120.csv' # dPhi/Phi. Needs to divide by dE.

energy_DUNE = []
flux_DUNE = []
bins_DUNE = []

energy_Alt = []
flux_Alt = []

energy_Alt120 = []
flux_Alt120 = []
bins_Alt120 = []

energy_DUNE_tau_opt = []
flux_DUNE_tau_opt = []

### Cross Sections ###
# Filenames #
xsec_1tau_coh_Ar_filename = CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/coherent/argon/vmu_to_vtau_tau+_mu-_coh_Ar_xsec.csv'
xsec_1tau_p_filename = CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv'
xsec_1tau_n_filename = CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/neutron/vmu_to_vtau_tau+_mu-_nucleon_n_xsec.csv'

xsec_2tau_coh_Ar_filename = CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/coherent/argon/vmu_to_vmu_tau+_tau-_coh_Ar_xsec.csv'
xsec_2tau_p_filename = CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv'
xsec_2tau_n_filename = CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/neutron/vmu_to_vmu_tau+_tau-_nucleon_n_xsec.csv'

xsec_2mu_coh_Ar_filename  = CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/argon/vmu_to_vmu_mu+_mu-_coh_Ar_xsec.csv'
xsec_2mu_p_filename  = CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/proton/vmu_to_vmu_mu+_mu-_nucleon_p_xsec.csv'
xsec_2mu_n_filename  = CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/neutron/vmu_to_vmu_mu+_mu-_nucleon_n_xsec.csv'

# Coherent ; Argon #
energy_1tau_coh_Ar = []
xsec_1tau_coh_Ar = []

energy_2tau_coh_Ar = []
xsec_2tau_coh_Ar = []

energy_2mu_coh_Ar = []
xsec_2mu_coh_Ar = []

# Incoherent ; proton ; Argon #
energy_1tau_p = []
xsec_1tau_p = []

energy_2tau_p = []
xsec_2tau_p = []

energy_2mu_p = []
xsec_2mu_p = []

# Incoherent ; neutron ; Argon #
energy_1tau_n = []
xsec_1tau_n = []

energy_2tau_n = []
xsec_2tau_n = []

energy_2mu_n = []
xsec_2mu_n = []

########################
###### Load Files ######
########################

##### Fluxes #####
### DUNE ; Standard Mode Flux ###
#with open(DUNE_filename,'r') as csvfile:
#    data = csv.reader(csvfile, delimiter = ',')
#    for row in data:
#        energy_DUNE.append(float(row[0]))
#        flux_DUNE.append(float(row[1]) * 1e4)  DUNE differential flux is in [cm^-2 GeV^-1 POT^-1]; convert to [m^-2 GeV^-1 POT^-1]

### DUNE ; Standard Mode Flux ; Histogram version ###
#with open(DUNE_filename,'r') as csvfile:
#    data = csv.reader(csvfile, delimiter = ',')
#    for row in data:
#        low_bin = float(row[0])
#        high_bin = float(row[1])
#        flux = float(row[2])
#        delta = high_bin - low_bin
#        energy = delta/2 + low_bin
#        energy_DUNE.append(energy)
#        flux_DUNE.append(flux)
#        bins_DUNE.append(low_bin)
#
#bins_DUNE.append(120.0)

### DUNE ; Standard Mode Flux ###
with open(DUNE_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ' ')
    for row in data:
        energy_DUNE.append(float(row[0]))
        flux_DUNE.append(float(row[2])) # Histogram has units of [m^-2 GeV^-1 POT^-1]

### DUNE ; Altmannshofer digitized ###
with open(Alt_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_Alt.append(float(row[0]))
        flux_Alt.append(float(row[1]))

### DUNE ; Tau-optimized Flux ###
with open(DUNE_tau_opt_numu_flux,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ' ')
    for row in data:
        energy_DUNE_tau_opt.append(float(row[0]))
        flux_DUNE_tau_opt.append(float(row[2])) # Histogram has units of [m^-2 GeV^-1 POT^-1]

### DUNE ; Altmannshofer 120 GeV code ###
with open(Alt120_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        low_bin = float(row[0])
        high_bin = float(row[1])
        flux = float(row[2])
        delta = high_bin - low_bin
        energy = delta/2 + low_bin
        energy_Alt120.append(energy)
        flux_Alt120.append(flux / delta)   # To get the correct flux from the Altmannshofer code, we must divide the flux by the width of the energy bin. This value matches DUNE.
        bins_Alt120.append(low_bin)

bins_Alt120.append(68.0)

### Integrated Fluxes ###
DUNE_integrated_flux = simpson(flux_DUNE, x=energy_DUNE)
DUNE_tau_opt_integrated_flux = simpson(flux_DUNE_tau_opt, x=energy_DUNE_tau_opt)
Alt120_integrated_flux = Phi_Alt
Alt_integrated_flux = Phi_Alt

##### Cross Sections #####

### Coherent ###
# vmu -> vtau tau+ mu- ; coherent ; Argon #
with open(xsec_1tau_coh_Ar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_coh_Ar.append(float(row[0]))
        xsec_1tau_coh_Ar.append(float(row[1]) * 1e-43) # Convert to m^2.

# vmu -> vmu tau+ tau- ; coherent ; Argon #
with open(xsec_2tau_coh_Ar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_coh_Ar.append(float(row[0]))
        xsec_2tau_coh_Ar.append(float(row[1]) * 1e-43) # Convert to m^2.

# vmu -> vmu mu+ mu- ; coherent ; Argon #
with open(xsec_2mu_coh_Ar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_coh_Ar.append(float(row[0]))
        xsec_2mu_coh_Ar.append(float(row[1]) * 1e-43) # Convert to m^2.

### Incoherent ; proton ; Argon ###
# vmu -> vtau tau+ mu- #
with open(xsec_1tau_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_p.append(float(row[0]))
        xsec_1tau_p.append(float(row[1]) * Z * 1e-43) # Convert to m^2.

# vmu -> vmu tau+ tau- #
with open(xsec_2tau_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_p.append(float(row[0]))
        xsec_2tau_p.append(float(row[1]) * Z * 1e-43) # Convert to m^2.

# vmu -> vmu mu+ mu- #
with open(xsec_2mu_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_p.append(float(row[0]))
        xsec_2mu_p.append(float(row[1]) * Z * 1e-43) # Convert to m^2.

### Incoherent ; neutron ; Argon ###
# vmu -> vtau tau+ mu- #
with open(xsec_1tau_n_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_n.append(float(row[0]))
        xsec_1tau_n.append(float(row[1]) * (A - Z) * 1e-43) # Convert to m^2.

# vmu -> vmu tau+ tau- #
with open(xsec_2tau_n_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_n.append(float(row[0]))
        xsec_2tau_n.append(float(row[1]) * (A - Z) * 1e-43) # Convert to m^2.

# vmu -> vmu mu+ mu- #
with open(xsec_2mu_n_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_n.append(float(row[0]))
        xsec_2mu_n.append(float(row[1]) * (A - Z) * 1e-43) # Convert to m^2.

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
        for jj in range(trid_size):
            if flux_energy[ii] < trid_energy[0]: # If the flux energy range starts below the trident one, set the matched xsec to 0. Otherwise, these values will be skipped.r
                xsec_matched.append(0.0)
                break
            elif ((flux_energy[ii] >= trid_energy[jj]) and (flux_energy[ii] <= trid_energy[jj+1])): # Find the trident energy range for the flux energy. Interpolate there.
                if (trid_xsec[jj] < 0) or (trid_xsec[jj+1] < 0): # Some xsec are negative in the 2tau proton case; ad hoc solution for now.
                    xsec_matched.append(0.0)
                else:
                    xsec_interp = Interpolation(trid_energy[jj], trid_energy[jj+1], trid_xsec[jj], trid_xsec[jj+1], flux_energy[ii])
                    xsec_matched.append(xsec_interp)
    
    return xsec_matched

### DUNE Standard Flux ; Matched XSec ###
DUNE_xsec_1tau_coh_Ar_matched = MatchXSec(energy_DUNE, energy_1tau_coh_Ar, xsec_1tau_coh_Ar)
DUNE_xsec_1tau_p_matched = MatchXSec(energy_DUNE, energy_1tau_p, xsec_1tau_p)
DUNE_xsec_1tau_n_matched = MatchXSec(energy_DUNE, energy_1tau_n, xsec_1tau_n)

DUNE_xsec_2tau_coh_Ar_matched = MatchXSec(energy_DUNE, energy_2tau_coh_Ar, xsec_2tau_coh_Ar)
DUNE_xsec_2tau_p_matched = MatchXSec(energy_DUNE, energy_2tau_p, xsec_2tau_p)
DUNE_xsec_2tau_n_matched = MatchXSec(energy_DUNE, energy_2tau_n, xsec_2tau_n)

DUNE_xsec_2mu_coh_Ar_matched = MatchXSec(energy_DUNE, energy_2mu_coh_Ar, xsec_2mu_coh_Ar)
DUNE_xsec_2mu_p_matched = MatchXSec(energy_DUNE, energy_2mu_p, xsec_2mu_p)
DUNE_xsec_2mu_n_matched = MatchXSec(energy_DUNE, energy_2mu_n, xsec_2mu_n)

### Alt. Digitized Flux ; Matched XSec ###
Alt_xsec_1tau_coh_Ar_matched = MatchXSec(energy_Alt, energy_1tau_coh_Ar, xsec_1tau_coh_Ar)
Alt_xsec_1tau_p_matched = MatchXSec(energy_Alt, energy_1tau_p, xsec_1tau_p)
Alt_xsec_1tau_n_matched = MatchXSec(energy_Alt, energy_1tau_n, xsec_1tau_n)

Alt_xsec_2tau_coh_Ar_matched = MatchXSec(energy_Alt, energy_2tau_coh_Ar, xsec_2tau_coh_Ar)
Alt_xsec_2tau_p_matched = MatchXSec(energy_Alt, energy_2tau_p, xsec_2tau_p)
Alt_xsec_2tau_n_matched = MatchXSec(energy_Alt, energy_2tau_n, xsec_2tau_n)

Alt_xsec_2mu_coh_Ar_matched = MatchXSec(energy_Alt, energy_2mu_coh_Ar, xsec_2mu_coh_Ar)
Alt_xsec_2mu_p_matched = MatchXSec(energy_Alt, energy_2mu_p, xsec_2mu_p)
Alt_xsec_2mu_n_matched = MatchXSec(energy_Alt, energy_2mu_n, xsec_2mu_n)

### Alt. 120 GeV Flux ; Matched XSec ###
Alt120_xsec_2mu_coh_Ar_matched = MatchXSec(energy_Alt120, energy_2mu_coh_Ar, xsec_2mu_coh_Ar)
Alt120_xsec_2mu_p_matched = MatchXSec(energy_Alt120, energy_2mu_p, xsec_2mu_p)
Alt120_xsec_2mu_n_matched = MatchXSec(energy_Alt120, energy_2mu_n, xsec_2mu_n)

Alt120_xsec_2tau_coh_Ar_matched = MatchXSec(energy_Alt120, energy_2tau_coh_Ar, xsec_2tau_coh_Ar)
Alt120_xsec_2tau_p_matched = MatchXSec(energy_Alt120, energy_2tau_p, xsec_2tau_p)
Alt120_xsec_2tau_n_matched = MatchXSec(energy_Alt120, energy_2tau_n, xsec_2tau_n)

Alt120_xsec_1tau_coh_Ar_matched = MatchXSec(energy_Alt120, energy_1tau_coh_Ar, xsec_1tau_coh_Ar)
Alt120_xsec_1tau_p_matched = MatchXSec(energy_Alt120, energy_1tau_p, xsec_1tau_p)
Alt120_xsec_1tau_n_matched = MatchXSec(energy_Alt120, energy_1tau_n, xsec_1tau_n)

### DUNE Tau-optimized Flux ; Matched XSec ###
DUNE_tau_opt_xsec_1tau_coh_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_1tau_coh_Ar, xsec_1tau_coh_Ar)
DUNE_tau_opt_xsec_1tau_p_matched = MatchXSec(energy_DUNE_tau_opt, energy_1tau_p, xsec_1tau_p)
DUNE_tau_opt_xsec_1tau_n_matched = MatchXSec(energy_DUNE_tau_opt, energy_1tau_n, xsec_1tau_n)

DUNE_tau_opt_xsec_2tau_coh_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_2tau_coh_Ar, xsec_2tau_coh_Ar)
DUNE_tau_opt_xsec_2tau_p_matched = MatchXSec(energy_DUNE_tau_opt, energy_2tau_p, xsec_2tau_p)
DUNE_tau_opt_xsec_2tau_n_matched = MatchXSec(energy_DUNE_tau_opt, energy_2tau_n, xsec_2tau_n)

DUNE_tau_opt_xsec_2mu_coh_Ar_matched = MatchXSec(energy_DUNE_tau_opt, energy_2mu_coh_Ar, xsec_2mu_coh_Ar)
DUNE_tau_opt_xsec_2mu_p_matched = MatchXSec(energy_DUNE_tau_opt, energy_2mu_p, xsec_2mu_p)
DUNE_tau_opt_xsec_2mu_n_matched = MatchXSec(energy_DUNE_tau_opt, energy_2mu_n, xsec_2mu_n)

##########################################
###### Number of Events Calculation ######
##########################################

def CalculateEvents(flux, flux_energy, xsec_matched, NTONNES=1, NYEAR=1, normalized=False, Phi=1):
    """
    Calculate the expected number of events for a given neutrino flux and trident process.

    N = XSecCon * MD * N_POT * NTONNES * NYEAR / MAr

    where N is the total number of expected events, XSecCon is the flux-convoluted trident 
    cross section

    XSecCon = int dPhi/dE xsec(E) dE,

    MD is the detector mass of 1 tonne, N_POT is the protons-on-target 1.1e21, NTONNES is
    the number of tonnes in the detector (default is 1), NYEAR is the number of years of 
    running (default is 1) and MAr is the mass of argon in kg. The flag normalized indicates
    if the provided flux is already normalized (default is False). If True, the number of
    events is calculated by

    N = Phi * XSecCon * MD * N_POT * NTONNES * NYEAR / MAr

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
        N = XSecCon * MD * N_POT * NTONNES * NYEAR / MAr
        return N, XSecCon*1e43 / Phi 
    else:
        N = Phi * XSecCon * MD * N_POT * NTONNES * NYEAR / MAr
        return N, XSecCon*1e43

### DUNE Standard Flux Events ###
# Coherent ; Argon #
N_1tau_coh_Ar_DUNE_1_1, XSecCon_1tau_coh_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched)
N_1tau_coh_Ar_DUNE_50_3, XSecCon_1tau_coh_Ar_DUNE_50_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched, NTONNES=50, NYEAR=3)
N_1tau_coh_Ar_DUNE_147_3, XSecCon_1tau_coh_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched, NTONNES=147, NYEAR=3)

N_2tau_coh_Ar_DUNE_1_1, XSecCon_2tau_coh_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched)
N_2tau_coh_Ar_DUNE_50_3, XSecCon_2tau_coh_Ar_DUNE_50_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched, NTONNES=50, NYEAR=3)
N_2tau_coh_Ar_DUNE_147_3, XSecCon_2tau_coh_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched, NTONNES=147, NYEAR=3)

N_2mu_coh_Ar_DUNE_1_1, XSecCon_2mu_coh_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched)
N_2mu_coh_Ar_DUNE_50_3, XSecCon_2mu_coh_Ar_DUNE_50_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched, NTONNES=50, NYEAR=3)
N_2mu_coh_Ar_DUNE_147_3, XSecCon_2mu_coh_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched, NTONNES=147, NYEAR=3)

# Incoherent ; proton ; Argon #
N_1tau_p_DUNE_1_1, XSecCon_1tau_p_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_p_matched)
N_1tau_p_DUNE_50_3, XSecCon_1tau_p_DUNE_50_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_p_matched, NTONNES=50, NYEAR=3)
N_1tau_p_DUNE_147_3, XSecCon_1tau_p_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_p_matched, NTONNES=147, NYEAR=3)

N_2tau_p_DUNE_1_1, XSecCon_2tau_p_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_p_matched)
N_2tau_p_DUNE_50_3, XSecCon_2tau_p_DUNE_50_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_p_matched, NTONNES=50, NYEAR=3)
N_2tau_p_DUNE_147_3, XSecCon_2tau_p_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_p_matched, NTONNES=147, NYEAR=3)

N_2mu_p_DUNE_1_1, XSecCon_2mu_p_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_p_matched)
N_2mu_p_DUNE_50_3, XSecCon_2mu_p_DUNE_50_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_p_matched, NTONNES=50, NYEAR=3)
N_2mu_p_DUNE_147_3, XSecCon_2mu_p_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_p_matched, NTONNES=147, NYEAR=3)

# Incoherent ; neutron ; Argon #
N_1tau_n_DUNE_1_1, XSecCon_1tau_n_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_n_matched)
N_1tau_n_DUNE_50_3, XSecCon_1tau_n_DUNE_50_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_n_matched, NTONNES=50, NYEAR=3)
N_1tau_n_DUNE_147_3, XSecCon_1tau_n_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_n_matched, NTONNES=147, NYEAR=3)

N_2tau_n_DUNE_1_1, XSecCon_2tau_n_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_n_matched)
N_2tau_n_DUNE_50_3, XSecCon_2tau_n_DUNE_50_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_n_matched, NTONNES=50, NYEAR=3)
N_2tau_n_DUNE_147_3, XSecCon_2tau_n_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_n_matched, NTONNES=147, NYEAR=3)

N_2mu_n_DUNE_1_1, XSecCon_2mu_n_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_n_matched)
N_2mu_n_DUNE_50_3, XSecCon_2mu_n_DUNE_50_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_n_matched, NTONNES=50, NYEAR=3)
N_2mu_n_DUNE_147_3, XSecCon_2mu_n_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_n_matched, NTONNES=147, NYEAR=3)

# Incoherent ; proton + neutron ; Argon #
N_1tau_incoh_DUNE_1_1, XSecCon_1tau_incoh_DUNE_1_1 = (N_1tau_p_DUNE_1_1 + N_1tau_n_DUNE_1_1), (XSecCon_1tau_p_DUNE_1_1 + XSecCon_1tau_n_DUNE_1_1)
N_1tau_incoh_DUNE_50_3, XSecCon_1tau_incoh_DUNE_50_3 = (N_1tau_p_DUNE_50_3 + N_1tau_n_DUNE_50_3), (XSecCon_1tau_p_DUNE_50_3 + XSecCon_1tau_n_DUNE_50_3)
N_1tau_incoh_DUNE_147_3, XSecCon_1tau_incoh_DUNE_147_3 = (N_1tau_p_DUNE_147_3 + N_1tau_n_DUNE_147_3), (XSecCon_1tau_p_DUNE_147_3 + XSecCon_1tau_n_DUNE_147_3) 

N_2tau_incoh_DUNE_1_1, XSecCon_2tau_incoh_DUNE_1_1 = (N_2tau_p_DUNE_1_1 + N_2tau_n_DUNE_1_1), (XSecCon_2tau_p_DUNE_1_1 + XSecCon_2tau_n_DUNE_1_1)
N_2tau_incoh_DUNE_50_3, XSecCon_2tau_incoh_DUNE_50_3 = (N_2tau_p_DUNE_50_3 + N_2tau_n_DUNE_50_3), (XSecCon_2tau_p_DUNE_50_3 + XSecCon_2tau_n_DUNE_50_3)
N_2tau_incoh_DUNE_147_3, XSecCon_2tau_incoh_DUNE_147_3 = (N_2tau_p_DUNE_147_3 + N_2tau_n_DUNE_147_3), (XSecCon_2tau_p_DUNE_147_3 + XSecCon_2tau_n_DUNE_147_3) 

N_2mu_incoh_DUNE_1_1, XSecCon_2mu_incoh_DUNE_1_1 = (N_2mu_p_DUNE_1_1 + N_2mu_n_DUNE_1_1), (XSecCon_2mu_p_DUNE_1_1 + XSecCon_2mu_n_DUNE_1_1)
N_2mu_incoh_DUNE_50_3, XSecCon_2mu_incoh_DUNE_50_3 = (N_2mu_p_DUNE_50_3 + N_2mu_n_DUNE_50_3), (XSecCon_2mu_p_DUNE_50_3 + XSecCon_2mu_n_DUNE_50_3)
N_2mu_incoh_DUNE_147_3, XSecCon_2mu_incoh_DUNE_147_3 = (N_2mu_p_DUNE_147_3 + N_2mu_n_DUNE_147_3), (XSecCon_2mu_p_DUNE_147_3 + XSecCon_2mu_n_DUNE_147_3) 

### Altmannshofer Flux Events ###
# Coherent ; Argon #
N_1tau_coh_Ar_Alt_1_1, XSecCon_1tau_coh_Ar_Alt_1_1 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_1tau_coh_Ar_matched, normalized=True, Phi=Phi_Alt)
N_1tau_coh_Ar_Alt_50_3, XSecCon_1tau_coh_Ar_Alt_50_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_1tau_coh_Ar_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_1tau_coh_Ar_Alt_147_3, XSecCon_1tau_coh_Ar_Alt_147_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_1tau_coh_Ar_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2tau_coh_Ar_Alt_1_1, XSecCon_2tau_coh_Ar_Alt_1_1 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2tau_coh_Ar_matched, normalized=True, Phi=Phi_Alt)
N_2tau_coh_Ar_Alt_50_3, XSecCon_2tau_coh_Ar_Alt_50_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2tau_coh_Ar_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2tau_coh_Ar_Alt_147_3, XSecCon_2tau_coh_Ar_Alt_147_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2tau_coh_Ar_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2mu_coh_Ar_Alt_1_1, XSecCon_2mu_coh_Ar_Alt_1_1 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2mu_coh_Ar_matched, normalized=True, Phi=Phi_Alt)
N_2mu_coh_Ar_Alt_50_3, XSecCon_2mu_coh_Ar_Alt_50_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2mu_coh_Ar_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2mu_coh_Ar_Alt_147_3, XSecCon_2mu_coh_Ar_Alt_147_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2mu_coh_Ar_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

# Incoherent ; proton ; Argon #
N_1tau_p_Alt_1_1, XSecCon_1tau_p_Alt_1_1 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_1tau_p_matched, normalized=True, Phi=Phi_Alt)
N_1tau_p_Alt_50_3, XSecCon_1tau_p_Alt_50_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_1tau_p_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_1tau_p_Alt_147_3, XSecCon_1tau_p_Alt_147_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_1tau_p_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2tau_p_Alt_1_1, XSecCon_2tau_p_Alt_1_1 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2tau_p_matched, normalized=True, Phi=Phi_Alt)
N_2tau_p_Alt_50_3, XSecCon_2tau_p_Alt_50_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2tau_p_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2tau_p_Alt_147_3, XSecCon_2tau_p_Alt_147_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2tau_p_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2mu_p_Alt_1_1, XSecCon_2mu_p_Alt_1_1 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2mu_p_matched, normalized=True, Phi=Phi_Alt)
N_2mu_p_Alt_50_3, XSecCon_2mu_p_Alt_50_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2mu_p_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2mu_p_Alt_147_3, XSecCon_2mu_p_Alt_147_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2mu_p_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

# Incoherent ; neutron ; Argon
N_1tau_n_Alt_1_1, XSecCon_1tau_n_Alt_1_1 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_1tau_n_matched, normalized=True, Phi=Phi_Alt)
N_1tau_n_Alt_50_3, XSecCon_1tau_n_Alt_50_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_1tau_n_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_1tau_n_Alt_147_3, XSecCon_1tau_n_Alt_147_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_1tau_n_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2tau_n_Alt_1_1, XSecCon_2tau_n_Alt_1_1 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2tau_n_matched, normalized=True, Phi=Phi_Alt)
N_2tau_n_Alt_50_3, XSecCon_2tau_n_Alt_50_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2tau_n_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2tau_n_Alt_147_3, XSecCon_2tau_n_Alt_147_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2tau_n_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2mu_n_Alt_1_1, XSecCon_2mu_n_Alt_1_1 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2mu_n_matched, normalized=True, Phi=Phi_Alt)
N_2mu_n_Alt_50_3, XSecCon_2mu_n_Alt_50_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2mu_n_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2mu_n_Alt_147_3, XSecCon_2mu_n_Alt_147_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2mu_n_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

# Incoherent ; proton + neutron ; Argon #
N_1tau_incoh_Alt_1_1, XSecCon_1tau_incoh_Alt_1_1 = (N_1tau_p_Alt_1_1 + N_1tau_n_Alt_1_1), (XSecCon_1tau_p_Alt_1_1 + XSecCon_1tau_n_Alt_1_1)
N_1tau_incoh_Alt_50_3, XSecCon_1tau_incoh_Alt_50_3 = (N_1tau_p_Alt_50_3 + N_1tau_n_Alt_50_3), (XSecCon_1tau_p_Alt_50_3 + XSecCon_1tau_n_Alt_50_3)
N_1tau_incoh_Alt_147_3, XSecCon_1tau_incoh_Alt_147_3 = (N_1tau_p_Alt_147_3 + N_1tau_n_Alt_147_3), (XSecCon_1tau_p_Alt_147_3 + XSecCon_1tau_n_Alt_147_3) 

N_2tau_incoh_Alt_1_1, XSecCon_2tau_incoh_Alt_1_1 = (N_2tau_p_Alt_1_1 + N_2tau_n_Alt_1_1), (XSecCon_2tau_p_Alt_1_1 + XSecCon_2tau_n_Alt_1_1)
N_2tau_incoh_Alt_50_3, XSecCon_2tau_incoh_Alt_50_3 = (N_2tau_p_Alt_50_3 + N_2tau_n_Alt_50_3), (XSecCon_2tau_p_Alt_50_3 + XSecCon_2tau_n_Alt_50_3)
N_2tau_incoh_Alt_147_3, XSecCon_2tau_incoh_Alt_147_3 = (N_2tau_p_Alt_147_3 + N_2tau_n_Alt_147_3), (XSecCon_2tau_p_Alt_147_3 + XSecCon_2tau_n_Alt_147_3) 

N_2mu_incoh_Alt_1_1, XSecCon_2mu_incoh_Alt_1_1 = (N_2mu_p_Alt_1_1 + N_2mu_n_Alt_1_1), (XSecCon_2mu_p_Alt_1_1 + XSecCon_2mu_n_Alt_1_1)
N_2mu_incoh_Alt_50_3, XSecCon_2mu_incoh_Alt_50_3 = (N_2mu_p_Alt_50_3 + N_2mu_n_Alt_50_3), (XSecCon_2mu_p_Alt_50_3 + XSecCon_2mu_n_Alt_50_3)
N_2mu_incoh_Alt_147_3, XSecCon_2mu_incoh_Alt_147_3 = (N_2mu_p_Alt_147_3 + N_2mu_n_Alt_147_3), (XSecCon_2mu_p_Alt_147_3 + XSecCon_2mu_n_Alt_147_3) 


### Altmannshofer 120 GeV Flux Events (code) ###
# Coherent ; Argon #
N_1tau_coh_Ar_Alt120_1_1, XSecCon_1tau_coh_Ar_Alt120_1_1 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_1tau_coh_Ar_matched, normalized=True, Phi=Phi_Alt)
N_1tau_coh_Ar_Alt120_50_3, XSecCon_1tau_coh_Ar_Alt120_50_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_1tau_coh_Ar_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_1tau_coh_Ar_Alt120_147_3, XSecCon_1tau_coh_Ar_Alt120_147_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_1tau_coh_Ar_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2tau_coh_Ar_Alt120_1_1, XSecCon_2tau_coh_Ar_Alt120_1_1 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2tau_coh_Ar_matched, normalized=True, Phi=Phi_Alt)
N_2tau_coh_Ar_Alt120_50_3, XSecCon_2tau_coh_Ar_Alt120_50_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2tau_coh_Ar_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2tau_coh_Ar_Alt120_147_3, XSecCon_2tau_coh_Ar_Alt120_147_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2tau_coh_Ar_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2mu_coh_Ar_Alt120_1_1, XSecCon_2mu_coh_Ar_Alt120_1_1 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2mu_coh_Ar_matched, normalized=True, Phi=Phi_Alt)
N_2mu_coh_Ar_Alt120_50_3, XSecCon_2mu_coh_Ar_Alt120_50_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2mu_coh_Ar_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2mu_coh_Ar_Alt120_147_3, XSecCon_2mu_coh_Ar_Alt120_147_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2mu_coh_Ar_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

# Incoherent ; proton ; Argon #
N_1tau_p_Alt120_1_1, XSecCon_1tau_p_Alt120_1_1 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_1tau_p_matched, normalized=True, Phi=Phi_Alt)
N_1tau_p_Alt120_50_3, XSecCon_1tau_p_Alt120_50_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_1tau_p_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_1tau_p_Alt120_147_3, XSecCon_1tau_p_Alt120_147_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_1tau_p_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2tau_p_Alt120_1_1, XSecCon_2tau_p_Alt120_1_1 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2tau_p_matched, normalized=True, Phi=Phi_Alt)
N_2tau_p_Alt120_50_3, XSecCon_2tau_p_Alt120_50_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2tau_p_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2tau_p_Alt120_147_3, XSecCon_2tau_p_Alt120_147_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2tau_p_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2mu_p_Alt120_1_1, XSecCon_2mu_p_Alt120_1_1 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2mu_p_matched, normalized=True, Phi=Phi_Alt)
N_2mu_p_Alt120_50_3, XSecCon_2mu_p_Alt120_50_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2mu_p_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2mu_p_Alt120_147_3, XSecCon_2mu_p_Alt120_147_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2mu_p_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

# Incoherent ; neutron ; Argon #
N_1tau_n_Alt120_1_1, XSecCon_1tau_n_Alt120_1_1 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_1tau_n_matched, normalized=True, Phi=Phi_Alt)
N_1tau_n_Alt120_50_3, XSecCon_1tau_n_Alt120_50_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_1tau_n_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_1tau_n_Alt120_147_3, XSecCon_1tau_n_Alt120_147_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_1tau_n_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2tau_n_Alt120_1_1, XSecCon_2tau_n_Alt120_1_1 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2tau_n_matched, normalized=True, Phi=Phi_Alt)
N_2tau_n_Alt120_50_3, XSecCon_2tau_n_Alt120_50_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2tau_n_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2tau_n_Alt120_147_3, XSecCon_2tau_n_Alt120_147_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2tau_n_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2mu_n_Alt120_1_1, XSecCon_2mu_n_Alt120_1_1 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2mu_n_matched, normalized=True, Phi=Phi_Alt)
N_2mu_n_Alt120_50_3, XSecCon_2mu_n_Alt120_50_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2mu_n_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2mu_n_Alt120_147_3, XSecCon_2mu_n_Alt120_147_3 = CalculateEvents(flux_Alt120, energy_Alt120, Alt120_xsec_2mu_n_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

# Incoherent ; proton + neutron ; Argon #
N_1tau_incoh_Alt120_1_1, XSecCon_1tau_incoh_Alt120_1_1 = (N_1tau_p_Alt120_1_1 + N_1tau_n_Alt120_1_1), (XSecCon_1tau_p_Alt120_1_1 + XSecCon_1tau_n_Alt120_1_1)
N_1tau_incoh_Alt120_50_3, XSecCon_1tau_incoh_Alt120_50_3 = (N_1tau_p_Alt120_50_3 + N_1tau_n_Alt120_50_3), (XSecCon_1tau_p_Alt120_50_3 + XSecCon_1tau_n_Alt120_50_3)
N_1tau_incoh_Alt120_147_3, XSecCon_1tau_incoh_Alt120_147_3 = (N_1tau_p_Alt120_147_3 + N_1tau_n_Alt120_147_3), (XSecCon_1tau_p_Alt120_147_3 + XSecCon_1tau_n_Alt120_147_3) 

N_2tau_incoh_Alt120_1_1, XSecCon_2tau_incoh_Alt120_1_1 = (N_2tau_p_Alt120_1_1 + N_2tau_n_Alt120_1_1), (XSecCon_2tau_p_Alt120_1_1 + XSecCon_2tau_n_Alt120_1_1)
N_2tau_incoh_Alt120_50_3, XSecCon_2tau_incoh_Alt120_50_3 = (N_2tau_p_Alt120_50_3 + N_2tau_n_Alt120_50_3), (XSecCon_2tau_p_Alt120_50_3 + XSecCon_2tau_n_Alt120_50_3)
N_2tau_incoh_Alt120_147_3, XSecCon_2tau_incoh_Alt120_147_3 = (N_2tau_p_Alt120_147_3 + N_2tau_n_Alt120_147_3), (XSecCon_2tau_p_Alt120_147_3 + XSecCon_2tau_n_Alt120_147_3) 

N_2mu_incoh_Alt120_1_1, XSecCon_2mu_incoh_Alt120_1_1 = (N_2mu_p_Alt120_1_1 + N_2mu_n_Alt120_1_1), (XSecCon_2mu_p_Alt120_1_1 + XSecCon_2mu_n_Alt120_1_1)
N_2mu_incoh_Alt120_50_3, XSecCon_2mu_incoh_Alt120_50_3 = (N_2mu_p_Alt120_50_3 + N_2mu_n_Alt120_50_3), (XSecCon_2mu_p_Alt120_50_3 + XSecCon_2mu_n_Alt120_50_3)
N_2mu_incoh_Alt120_147_3, XSecCon_2mu_incoh_Alt120_147_3 = (N_2mu_p_Alt120_147_3 + N_2mu_n_Alt120_147_3), (XSecCon_2mu_p_Alt120_147_3 + XSecCon_2mu_n_Alt120_147_3) 


### DUNE Tau-Optimized ND Flux Events ###
# Coherent ; Argon #
N_1tau_coh_Ar_DUNE_tau_opt_1_1, XSecCon_1tau_coh_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched)
N_1tau_coh_Ar_DUNE_tau_opt_50_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_50_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched, NTONNES=50, NYEAR=3)
N_1tau_coh_Ar_DUNE_tau_opt_147_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_coh_Ar_matched, NTONNES=147, NYEAR=3)

N_2tau_coh_Ar_DUNE_tau_opt_1_1, XSecCon_2tau_coh_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched)
N_2tau_coh_Ar_DUNE_tau_opt_50_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_50_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched, NTONNES=50, NYEAR=3)
N_2tau_coh_Ar_DUNE_tau_opt_147_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_coh_Ar_matched, NTONNES=147, NYEAR=3)

N_2mu_coh_Ar_DUNE_tau_opt_1_1, XSecCon_2mu_coh_Ar_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched)
N_2mu_coh_Ar_DUNE_tau_opt_50_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_50_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched, NTONNES=50, NYEAR=3)
N_2mu_coh_Ar_DUNE_tau_opt_147_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_coh_Ar_matched, NTONNES=147, NYEAR=3)

# Incoherent ; proton ; Argon #
N_1tau_p_DUNE_tau_opt_1_1, XSecCon_1tau_p_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_matched)
N_1tau_p_DUNE_tau_opt_50_3, XSecCon_1tau_p_DUNE_tau_opt_50_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_matched, NTONNES=50, NYEAR=3)
N_1tau_p_DUNE_tau_opt_147_3, XSecCon_1tau_p_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_p_matched, NTONNES=147, NYEAR=3)

N_2tau_p_DUNE_tau_opt_1_1, XSecCon_2tau_p_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_matched)
N_2tau_p_DUNE_tau_opt_50_3, XSecCon_2tau_p_DUNE_tau_opt_50_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_matched, NTONNES=50, NYEAR=3)
N_2tau_p_DUNE_tau_opt_147_3, XSecCon_2tau_p_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_p_matched, NTONNES=147, NYEAR=3)

N_2mu_p_DUNE_tau_opt_1_1, XSecCon_2mu_p_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_matched)
N_2mu_p_DUNE_tau_opt_50_3, XSecCon_2mu_p_DUNE_tau_opt_50_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_matched, NTONNES=50, NYEAR=3)
N_2mu_p_DUNE_tau_opt_147_3, XSecCon_2mu_p_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_p_matched, NTONNES=147, NYEAR=3)

# Incoherent ; neutron ; Argon #
N_1tau_n_DUNE_tau_opt_1_1, XSecCon_1tau_n_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_matched)
N_1tau_n_DUNE_tau_opt_50_3, XSecCon_1tau_n_DUNE_tau_opt_50_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_matched, NTONNES=50, NYEAR=3)
N_1tau_n_DUNE_tau_opt_147_3, XSecCon_1tau_n_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_1tau_n_matched, NTONNES=147, NYEAR=3)

N_2tau_n_DUNE_tau_opt_1_1, XSecCon_2tau_n_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_matched)
N_2tau_n_DUNE_tau_opt_50_3, XSecCon_2tau_n_DUNE_tau_opt_50_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_matched, NTONNES=50, NYEAR=3)
N_2tau_n_DUNE_tau_opt_147_3, XSecCon_2tau_n_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2tau_n_matched, NTONNES=147, NYEAR=3)

N_2mu_n_DUNE_tau_opt_1_1, XSecCon_2mu_n_DUNE_tau_opt_1_1 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_matched)
N_2mu_n_DUNE_tau_opt_50_3, XSecCon_2mu_n_DUNE_tau_opt_50_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_matched, NTONNES=50, NYEAR=3)
N_2mu_n_DUNE_tau_opt_147_3, XSecCon_2mu_n_DUNE_tau_opt_147_3 = CalculateEvents(flux_DUNE_tau_opt, energy_DUNE_tau_opt, DUNE_tau_opt_xsec_2mu_n_matched, NTONNES=147, NYEAR=3)

# Incoherent ; proton + neutron ; Argon #
N_1tau_incoh_DUNE_tau_opt_1_1, XSecCon_1tau_incoh_DUNE_tau_opt_1_1 = (N_1tau_p_DUNE_tau_opt_1_1 + N_1tau_n_DUNE_tau_opt_1_1), (XSecCon_1tau_p_DUNE_tau_opt_1_1 + XSecCon_1tau_n_DUNE_tau_opt_1_1)
N_1tau_incoh_DUNE_tau_opt_50_3, XSecCon_1tau_incoh_DUNE_tau_opt_50_3 = (N_1tau_p_DUNE_tau_opt_50_3 + N_1tau_n_DUNE_tau_opt_50_3), (XSecCon_1tau_p_DUNE_tau_opt_50_3 + XSecCon_1tau_n_DUNE_tau_opt_50_3)
N_1tau_incoh_DUNE_tau_opt_147_3, XSecCon_1tau_incoh_DUNE_tau_opt_147_3 = (N_1tau_p_DUNE_tau_opt_147_3 + N_1tau_n_DUNE_tau_opt_147_3), (XSecCon_1tau_p_DUNE_tau_opt_147_3 + XSecCon_1tau_n_DUNE_tau_opt_147_3) 

N_2tau_incoh_DUNE_tau_opt_1_1, XSecCon_2tau_incoh_DUNE_tau_opt_1_1 = (N_2tau_p_DUNE_tau_opt_1_1 + N_2tau_n_DUNE_tau_opt_1_1), (XSecCon_2tau_p_DUNE_tau_opt_1_1 + XSecCon_2tau_n_DUNE_tau_opt_1_1)
N_2tau_incoh_DUNE_tau_opt_50_3, XSecCon_2tau_incoh_DUNE_tau_opt_50_3 = (N_2tau_p_DUNE_tau_opt_50_3 + N_2tau_n_DUNE_tau_opt_50_3), (XSecCon_2tau_p_DUNE_tau_opt_50_3 + XSecCon_2tau_n_DUNE_tau_opt_50_3)
N_2tau_incoh_DUNE_tau_opt_147_3, XSecCon_2tau_incoh_DUNE_tau_opt_147_3 = (N_2tau_p_DUNE_tau_opt_147_3 + N_2tau_n_DUNE_tau_opt_147_3), (XSecCon_2tau_p_DUNE_tau_opt_147_3 + XSecCon_2tau_n_DUNE_tau_opt_147_3) 

N_2mu_incoh_DUNE_tau_opt_1_1, XSecCon_2mu_incoh_DUNE_tau_opt_1_1 = (N_2mu_p_DUNE_tau_opt_1_1 + N_2mu_n_DUNE_tau_opt_1_1), (XSecCon_2mu_p_DUNE_tau_opt_1_1 + XSecCon_2mu_n_DUNE_tau_opt_1_1)
N_2mu_incoh_DUNE_tau_opt_50_3, XSecCon_2mu_incoh_DUNE_tau_opt_50_3 = (N_2mu_p_DUNE_tau_opt_50_3 + N_2mu_n_DUNE_tau_opt_50_3), (XSecCon_2mu_p_DUNE_tau_opt_50_3 + XSecCon_2mu_n_DUNE_tau_opt_50_3)
N_2mu_incoh_DUNE_tau_opt_147_3, XSecCon_2mu_incoh_DUNE_tau_opt_147_3 = (N_2mu_p_DUNE_tau_opt_147_3 + N_2mu_n_DUNE_tau_opt_147_3), (XSecCon_2mu_p_DUNE_tau_opt_147_3 + XSecCon_2mu_n_DUNE_tau_opt_147_3) 

############################
###### Event Printout ######
############################

def PrintOutEvent(filename, fluxname, tonnes, years, events, xsec_con):
    print("\t{} tonne Ar, {} year, {} flux: {:.3e} ; XSecCon: {:.3e}".format(tonnes, years, fluxname, events, xsec_con), file=filename)

def WriteOutFile(filename):
    with open(filename,'w') as textfile:
        print("Events expected at DUNE ND for:", file=textfile)
        print("The following integrated fluxes were used:", file=textfile)
        print("\tDUNE Standard Flux: {:.3e}".format(DUNE_integrated_flux), file=textfile)
        print("\tDUNE Tau-optimized Flux: {:.3e}".format(DUNE_tau_opt_integrated_flux), file=textfile)
        print("\tAlt. Digitized Flux: {:.3e}".format(Alt_integrated_flux), file=textfile)
        print("\tAlt. 120 GeV Flux: {:.3e}".format(Alt120_integrated_flux), file=textfile)
        print("\n", file=textfile)

        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_1tau_coh_Ar_DUNE_1_1, XSecCon_1tau_coh_Ar_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_1tau_coh_Ar_DUNE_50_3, XSecCon_1tau_coh_Ar_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_1tau_coh_Ar_DUNE_147_3, XSecCon_1tau_coh_Ar_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_1tau_coh_Ar_Alt_1_1, XSecCon_1tau_coh_Ar_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_1tau_coh_Ar_Alt_50_3, XSecCon_1tau_coh_Ar_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_1tau_coh_Ar_Alt_147_3, XSecCon_1tau_coh_Ar_Alt_147_3)
        PrintOutEvent(textfile, "Alt120", 1, 1, N_1tau_coh_Ar_Alt120_1_1, XSecCon_1tau_coh_Ar_Alt120_1_1)
        PrintOutEvent(textfile, "Alt120", 50, 3, N_1tau_coh_Ar_Alt120_50_3, XSecCon_1tau_coh_Ar_Alt120_50_3)
        PrintOutEvent(textfile, "Alt120", 147, 3, N_1tau_coh_Ar_Alt120_147_3, XSecCon_1tau_coh_Ar_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_1tau_coh_Ar_DUNE_tau_opt_1_1, XSecCon_1tau_coh_Ar_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 50, 3, N_1tau_coh_Ar_DUNE_tau_opt_50_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_50_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_1tau_coh_Ar_DUNE_tau_opt_147_3, XSecCon_1tau_coh_Ar_DUNE_tau_opt_147_3)
        print("1tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_1tau_p_DUNE_1_1, XSecCon_1tau_p_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_1tau_p_DUNE_50_3, XSecCon_1tau_p_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_1tau_p_DUNE_147_3, XSecCon_1tau_p_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_1tau_p_Alt_1_1, XSecCon_1tau_p_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_1tau_p_Alt_50_3, XSecCon_1tau_p_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_1tau_p_Alt_147_3, XSecCon_1tau_p_Alt_147_3)
        PrintOutEvent(textfile, "Alt120", 1, 1, N_1tau_p_Alt120_1_1, XSecCon_1tau_p_Alt120_1_1)
        PrintOutEvent(textfile, "Alt120", 50, 3, N_1tau_p_Alt120_50_3, XSecCon_1tau_p_Alt120_50_3)
        PrintOutEvent(textfile, "Alt120", 147, 3, N_1tau_p_Alt120_147_3, XSecCon_1tau_p_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_1tau_p_DUNE_tau_opt_1_1, XSecCon_1tau_p_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 50, 3, N_1tau_p_DUNE_tau_opt_50_3, XSecCon_1tau_p_DUNE_tau_opt_50_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_1tau_p_DUNE_tau_opt_147_3, XSecCon_1tau_p_DUNE_tau_opt_147_3)
        print("1tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_1tau_n_DUNE_1_1, XSecCon_1tau_n_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_1tau_n_DUNE_50_3, XSecCon_1tau_n_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_1tau_n_DUNE_147_3, XSecCon_1tau_n_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_1tau_n_Alt_1_1, XSecCon_1tau_n_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_1tau_n_Alt_50_3, XSecCon_1tau_n_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_1tau_n_Alt_147_3, XSecCon_1tau_n_Alt_147_3)
        PrintOutEvent(textfile, "Alt120", 1, 1, N_1tau_n_Alt120_1_1, XSecCon_1tau_n_Alt120_1_1)
        PrintOutEvent(textfile, "Alt120", 50, 3, N_1tau_n_Alt120_50_3, XSecCon_1tau_n_Alt120_50_3)
        PrintOutEvent(textfile, "Alt120", 147, 3, N_1tau_n_Alt120_147_3, XSecCon_1tau_n_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_1tau_n_DUNE_tau_opt_1_1, XSecCon_1tau_n_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 50, 3, N_1tau_n_DUNE_tau_opt_50_3, XSecCon_1tau_n_DUNE_tau_opt_50_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_1tau_n_DUNE_tau_opt_147_3, XSecCon_1tau_n_DUNE_tau_opt_147_3)
        print("1tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_1tau_incoh_DUNE_1_1, XSecCon_1tau_incoh_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_1tau_incoh_DUNE_50_3, XSecCon_1tau_incoh_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_1tau_incoh_DUNE_147_3, XSecCon_1tau_incoh_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_1tau_incoh_Alt_1_1, XSecCon_1tau_incoh_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_1tau_incoh_Alt_50_3, XSecCon_1tau_incoh_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_1tau_incoh_Alt_147_3, XSecCon_1tau_incoh_Alt_147_3)
        PrintOutEvent(textfile, "Alt120", 1, 1, N_1tau_incoh_Alt120_1_1, XSecCon_1tau_incoh_Alt120_1_1)
        PrintOutEvent(textfile, "Alt120", 50, 3, N_1tau_incoh_Alt120_50_3, XSecCon_1tau_incoh_Alt120_50_3)
        PrintOutEvent(textfile, "Alt120", 147, 3, N_1tau_incoh_Alt120_147_3, XSecCon_1tau_incoh_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_1tau_incoh_DUNE_tau_opt_1_1, XSecCon_1tau_incoh_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 50, 3, N_1tau_incoh_DUNE_tau_opt_50_3, XSecCon_1tau_incoh_DUNE_tau_opt_50_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_1tau_incoh_DUNE_tau_opt_147_3, XSecCon_1tau_incoh_DUNE_tau_opt_147_3)
        print("\n", file=textfile)


        print("2tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2tau_coh_Ar_DUNE_1_1, XSecCon_2tau_coh_Ar_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_2tau_coh_Ar_DUNE_50_3, XSecCon_2tau_coh_Ar_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2tau_coh_Ar_DUNE_147_3, XSecCon_2tau_coh_Ar_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_2tau_coh_Ar_Alt_1_1, XSecCon_2tau_coh_Ar_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_2tau_coh_Ar_Alt_50_3, XSecCon_2tau_coh_Ar_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_2tau_coh_Ar_Alt_147_3, XSecCon_2tau_coh_Ar_Alt_147_3)
        PrintOutEvent(textfile, "Alt120", 1, 1, N_2tau_coh_Ar_Alt120_1_1, XSecCon_2tau_coh_Ar_Alt120_1_1)
        PrintOutEvent(textfile, "Alt120", 50, 3, N_2tau_coh_Ar_Alt120_50_3, XSecCon_2tau_coh_Ar_Alt120_50_3)
        PrintOutEvent(textfile, "Alt120", 147, 3, N_2tau_coh_Ar_Alt120_147_3, XSecCon_2tau_coh_Ar_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2tau_coh_Ar_DUNE_tau_opt_1_1, XSecCon_2tau_coh_Ar_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 50, 3, N_2tau_coh_Ar_DUNE_tau_opt_50_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_50_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2tau_coh_Ar_DUNE_tau_opt_147_3, XSecCon_2tau_coh_Ar_DUNE_tau_opt_147_3)
        print("2tau Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2tau_p_DUNE_1_1, XSecCon_2tau_p_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_2tau_p_DUNE_50_3, XSecCon_2tau_p_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2tau_p_DUNE_147_3, XSecCon_2tau_p_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_2tau_p_Alt_1_1, XSecCon_2tau_p_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_2tau_p_Alt_50_3, XSecCon_2tau_p_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_2tau_p_Alt_147_3, XSecCon_2tau_p_Alt_147_3)
        PrintOutEvent(textfile, "Alt120", 1, 1, N_2tau_p_Alt120_1_1, XSecCon_2tau_p_Alt120_1_1)
        PrintOutEvent(textfile, "Alt120", 50, 3, N_2tau_p_Alt120_50_3, XSecCon_2tau_p_Alt120_50_3)
        PrintOutEvent(textfile, "Alt120", 147, 3, N_2tau_p_Alt120_147_3, XSecCon_2tau_p_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2tau_p_DUNE_tau_opt_1_1, XSecCon_2tau_p_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 50, 3, N_2tau_p_DUNE_tau_opt_50_3, XSecCon_2tau_p_DUNE_tau_opt_50_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2tau_p_DUNE_tau_opt_147_3, XSecCon_2tau_p_DUNE_tau_opt_147_3)
        print("2tau Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2tau_n_DUNE_1_1, XSecCon_2tau_n_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_2tau_n_DUNE_50_3, XSecCon_2tau_n_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2tau_n_DUNE_147_3, XSecCon_2tau_n_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_2tau_n_Alt_1_1, XSecCon_2tau_n_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_2tau_n_Alt_50_3, XSecCon_2tau_n_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_2tau_n_Alt_147_3, XSecCon_2tau_n_Alt_147_3)
        PrintOutEvent(textfile, "Alt120", 1, 1, N_2tau_n_Alt120_1_1, XSecCon_2tau_n_Alt120_1_1)
        PrintOutEvent(textfile, "Alt120", 50, 3, N_2tau_n_Alt120_50_3, XSecCon_2tau_n_Alt120_50_3)
        PrintOutEvent(textfile, "Alt120", 147, 3, N_2tau_n_Alt120_147_3, XSecCon_2tau_n_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2tau_n_DUNE_tau_opt_1_1, XSecCon_2tau_n_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 50, 3, N_2tau_n_DUNE_tau_opt_50_3, XSecCon_2tau_n_DUNE_tau_opt_50_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2tau_n_DUNE_tau_opt_147_3, XSecCon_2tau_n_DUNE_tau_opt_147_3)
        print("2tau Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2tau_incoh_DUNE_1_1, XSecCon_2tau_incoh_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_2tau_incoh_DUNE_50_3, XSecCon_2tau_incoh_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2tau_incoh_DUNE_147_3, XSecCon_2tau_incoh_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_2tau_incoh_Alt_1_1, XSecCon_2tau_incoh_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_2tau_incoh_Alt_50_3, XSecCon_2tau_incoh_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_2tau_incoh_Alt_147_3, XSecCon_2tau_incoh_Alt_147_3)
        PrintOutEvent(textfile, "Alt120", 1, 1, N_2tau_incoh_Alt120_1_1, XSecCon_2tau_incoh_Alt120_1_1)
        PrintOutEvent(textfile, "Alt120", 50, 3, N_2tau_incoh_Alt120_50_3, XSecCon_2tau_incoh_Alt120_50_3)
        PrintOutEvent(textfile, "Alt120", 147, 3, N_2tau_incoh_Alt120_147_3, XSecCon_2tau_incoh_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2tau_incoh_DUNE_tau_opt_1_1, XSecCon_2tau_incoh_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 50, 3, N_2tau_incoh_DUNE_tau_opt_50_3, XSecCon_2tau_incoh_DUNE_tau_opt_50_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2tau_incoh_DUNE_tau_opt_147_3, XSecCon_2tau_incoh_DUNE_tau_opt_147_3)
        print("\n", file=textfile)

        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2mu_coh_Ar_DUNE_1_1, XSecCon_2mu_coh_Ar_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_2mu_coh_Ar_DUNE_50_3, XSecCon_2mu_coh_Ar_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2mu_coh_Ar_DUNE_147_3, XSecCon_2mu_coh_Ar_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_2mu_coh_Ar_Alt_1_1, XSecCon_2mu_coh_Ar_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_2mu_coh_Ar_Alt_50_3, XSecCon_2mu_coh_Ar_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_2mu_coh_Ar_Alt_147_3, XSecCon_2mu_coh_Ar_Alt_147_3)
        PrintOutEvent(textfile, "Alt120", 1, 1, N_2mu_coh_Ar_Alt120_1_1, XSecCon_2mu_coh_Ar_Alt120_1_1)
        PrintOutEvent(textfile, "Alt120", 50, 3, N_2mu_coh_Ar_Alt120_50_3, XSecCon_2mu_coh_Ar_Alt120_50_3)
        PrintOutEvent(textfile, "Alt120", 147, 3, N_2mu_coh_Ar_Alt120_147_3, XSecCon_2mu_coh_Ar_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2mu_coh_Ar_DUNE_tau_opt_1_1, XSecCon_2mu_coh_Ar_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 50, 3, N_2mu_coh_Ar_DUNE_tau_opt_50_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_50_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2mu_coh_Ar_DUNE_tau_opt_147_3, XSecCon_2mu_coh_Ar_DUNE_tau_opt_147_3)
        print("2mu Incoherent ; proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2mu_p_DUNE_1_1, XSecCon_2mu_p_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_2mu_p_DUNE_50_3, XSecCon_2mu_p_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2mu_p_DUNE_147_3, XSecCon_2mu_p_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_2mu_p_Alt_1_1, XSecCon_2mu_p_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_2mu_p_Alt_50_3, XSecCon_2mu_p_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_2mu_p_Alt_147_3, XSecCon_2mu_p_Alt_147_3)
        PrintOutEvent(textfile, "Alt120", 1, 1, N_2mu_p_Alt120_1_1, XSecCon_2mu_p_Alt120_1_1)
        PrintOutEvent(textfile, "Alt120", 50, 3, N_2mu_p_Alt120_50_3, XSecCon_2mu_p_Alt120_50_3)
        PrintOutEvent(textfile, "Alt120", 147, 3, N_2mu_p_Alt120_147_3, XSecCon_2mu_p_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2mu_p_DUNE_tau_opt_1_1, XSecCon_2mu_p_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 50, 3, N_2mu_p_DUNE_tau_opt_50_3, XSecCon_2mu_p_DUNE_tau_opt_50_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2mu_p_DUNE_tau_opt_147_3, XSecCon_2mu_p_DUNE_tau_opt_147_3)
        print("2mu Incoherent ; neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2mu_n_DUNE_1_1, XSecCon_2mu_n_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_2mu_n_DUNE_50_3, XSecCon_2mu_n_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2mu_n_DUNE_147_3, XSecCon_2mu_n_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_2mu_n_Alt_1_1, XSecCon_2mu_n_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_2mu_n_Alt_50_3, XSecCon_2mu_n_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_2mu_n_Alt_147_3, XSecCon_2mu_n_Alt_147_3)
        PrintOutEvent(textfile, "Alt120", 1, 1, N_2mu_n_Alt120_1_1, XSecCon_2mu_n_Alt120_1_1)
        PrintOutEvent(textfile, "Alt120", 50, 3, N_2mu_n_Alt120_50_3, XSecCon_2mu_n_Alt120_50_3)
        PrintOutEvent(textfile, "Alt120", 147, 3, N_2mu_n_Alt120_147_3, XSecCon_2mu_n_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2mu_n_DUNE_tau_opt_1_1, XSecCon_2mu_n_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 50, 3, N_2mu_n_DUNE_tau_opt_50_3, XSecCon_2mu_n_DUNE_tau_opt_50_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2mu_n_DUNE_tau_opt_147_3, XSecCon_2mu_n_DUNE_tau_opt_147_3)
        print("2mu Incoherent ; proton + neutron:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2mu_incoh_DUNE_1_1, XSecCon_2mu_incoh_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_2mu_incoh_DUNE_50_3, XSecCon_2mu_incoh_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2mu_incoh_DUNE_147_3, XSecCon_2mu_incoh_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_2mu_incoh_Alt_1_1, XSecCon_2mu_incoh_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_2mu_incoh_Alt_50_3, XSecCon_2mu_incoh_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_2mu_incoh_Alt_147_3, XSecCon_2mu_incoh_Alt_147_3)
        PrintOutEvent(textfile, "Alt120", 1, 1, N_2mu_incoh_Alt120_1_1, XSecCon_2mu_incoh_Alt120_1_1)
        PrintOutEvent(textfile, "Alt120", 50, 3, N_2mu_incoh_Alt120_50_3, XSecCon_2mu_incoh_Alt120_50_3)
        PrintOutEvent(textfile, "Alt120", 147, 3, N_2mu_incoh_Alt120_147_3, XSecCon_2mu_incoh_Alt120_147_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 1, 1, N_2mu_incoh_DUNE_tau_opt_1_1, XSecCon_2mu_incoh_DUNE_tau_opt_1_1)
        PrintOutEvent(textfile, "DUNE tau-opt", 50, 3, N_2mu_incoh_DUNE_tau_opt_50_3, XSecCon_2mu_incoh_DUNE_tau_opt_50_3)
        PrintOutEvent(textfile, "DUNE tau-opt", 147, 3, N_2mu_incoh_DUNE_tau_opt_147_3, XSecCon_2mu_incoh_DUNE_tau_opt_147_3)
        print("\n", file=textfile)


WriteOutFile("number_of_events.txt")

################################
###### Debugging Plotting ######
################################

print("DUNE integrated flux: ", DUNE_integrated_flux)
print("DUNE tau-opt integrated flux: ", DUNE_tau_opt_integrated_flux)
print("Altmannshofer integrated flux: ", Phi_Alt)

print("Length of DUNE energy flux array: ", len(energy_DUNE))
print("Length of DUNE 2mu coh Ar matched xsec array: ", len(DUNE_xsec_2mu_coh_Ar_matched))
print("Length of DUNE 2mu p matched xsec array: ", len(DUNE_xsec_2mu_p_matched))

print("Length of Alt. digitized energy flux array: ", len(energy_Alt))
print("Length of Alt. digitized 2mu coh Ar matched xsec array: ", len(Alt_xsec_2mu_coh_Ar_matched))
print("Length of Alt. digitized 2mu p matched xsec array: ", len(Alt_xsec_2mu_p_matched))

print("Length of Alt. 120 energy flux array: ", len(energy_Alt120))
print("Length of Alt. 120 2mu coh Ar matched xsec array: ", len(Alt120_xsec_2mu_coh_Ar_matched))
print("Length of Alt. 120 2mu p matched xsec array: ", len(Alt120_xsec_2mu_p_matched))

print("Length of DUNE tau-opt energy flux array: ", len(energy_DUNE_tau_opt))
print("Length of DUNE tau-opt 2mu coh Ar matched xsec array: ", len(DUNE_tau_opt_xsec_2mu_coh_Ar_matched))
print("Length of DUNE tau-opt 2mu p matched xsec array: ", len(DUNE_tau_opt_xsec_2mu_p_matched))

flux_DUNE_norm = np.divide(flux_DUNE, DUNE_integrated_flux) # Normalize DUNE flux
flux_DUNE_tau_opt_norm = np.divide(flux_DUNE_tau_opt, DUNE_tau_opt_integrated_flux) # Normalize DUNE tau-opt flux

### Plot normalized fluxes ###
fig1, ax1 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)

#ax1.hist(energy_DUNE, bins=bins_DUNE, weights=flux_DUNE_norm, histtype='stepfilled', label=r'DUNE', color='navy', alpha=0.5, lw=0.5)
#ax1.hist(energy_DUNE, bins=bins_DUNE, weights=flux_DUNE_norm, histtype='step', color='black', lw=2, alpha=1)

ax1.hist(energy_Alt120, bins=bins_Alt120, weights=flux_Alt120, histtype='stepfilled', label=r'Alt. 120 GeV', color='orange', alpha=0.3, lw=0.5)
ax1.hist(energy_Alt120, bins=bins_Alt120, weights=flux_Alt120, histtype='step', color='black',lw=2, alpha=1)

ax1.plot(energy_DUNE, flux_DUNE_norm, '-', color='navy', label='DUNE', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax1.plot(energy_Alt, flux_Alt, '-', color='orange', label='Alt. digitized', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax1.plot(energy_DUNE_tau_opt, flux_DUNE_tau_opt_norm, '-', color='firebrick', label=r'DUNE $\tau-$opt', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])

ax1.set_xlabel('Energy (GeV)')
ax1.set_ylabel(r'$\frac{1}{\Phi}\frac{\mathrm{d}\Phi}{\mathrm{d}E}$ (GeV$^{-1}$)')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(0.3, 100)
#ax1.set_ylim(5e-4, 0.500)
ax1.set_ylim(1e-6, 0.500)
ax1.legend(loc='upper right')
ax1.set_title(r"$\nu_\mu$ Normalized Flux")

fig1.savefig("../plots/fluxes.png", dpi=400)

### Plot fluxes with cross sections ###
fig2, ax2 = plt.subplots(2, 3, figsize=(45, 20), sharex=True, tight_layout=True)

col_map = mpl.colormaps['Dark2'].resampled(4)
rgb_col = np.linspace(0,1,num=4)
c1 = col_map(rgb_col[0])
c2 = col_map(rgb_col[1])
c3 = col_map(rgb_col[2])
c4 = col_map(rgb_col[3])

## Data ##
# 2mu Coherent ; Argon #
mu2_coh_DUNE_integrand = np.multiply(flux_DUNE_norm, DUNE_xsec_2mu_coh_Ar_matched)
mu2_coh_DUNE_tau_opt_integrand = np.multiply(flux_DUNE_tau_opt_norm, DUNE_tau_opt_xsec_2mu_coh_Ar_matched)

# 2mu Incoherent ; proton + neutron ; Argon #
mu2_incoh_DUNE_integrand = np.multiply(flux_DUNE_norm, [sum(xsecs) for xsecs in zip(DUNE_xsec_2mu_p_matched, DUNE_xsec_2mu_n_matched)])
mu2_incoh_DUNE_tau_opt_integrand = np.multiply(flux_DUNE_tau_opt_norm, [sum(xsecs) for xsecs in zip(DUNE_tau_opt_xsec_2mu_p_matched, DUNE_tau_opt_xsec_2mu_n_matched)])

# 1tau Coherent ; Argon #
tau1_coh_DUNE_integrand = np.multiply(flux_DUNE_norm, DUNE_xsec_1tau_coh_Ar_matched)
tau1_coh_DUNE_tau_opt_integrand = np.multiply(flux_DUNE_tau_opt_norm, DUNE_tau_opt_xsec_1tau_coh_Ar_matched)

# 1tau Incoherent ; proton + neutron ; Argon #
tau1_incoh_DUNE_integrand = np.multiply(flux_DUNE_norm, [sum(xsecs) for xsecs in zip(DUNE_xsec_1tau_p_matched, DUNE_xsec_1tau_n_matched)])
tau1_incoh_DUNE_tau_opt_integrand = np.multiply(flux_DUNE_tau_opt_norm, [sum(xsecs) for xsecs in zip(DUNE_tau_opt_xsec_1tau_p_matched, DUNE_tau_opt_xsec_1tau_n_matched)])

# 2tau Coherent ; Argon #
tau2_coh_DUNE_integrand = np.multiply(flux_DUNE_norm, DUNE_xsec_2tau_coh_Ar_matched)
tau2_coh_DUNE_tau_opt_integrand = np.multiply(flux_DUNE_tau_opt_norm, DUNE_tau_opt_xsec_2tau_coh_Ar_matched)

# 2tau Incoherent ; proton + neutron ; Argon #
tau2_incoh_DUNE_integrand = np.multiply(flux_DUNE_norm, [sum(xsecs) for xsecs in zip(DUNE_xsec_2tau_p_matched, DUNE_xsec_2tau_n_matched)])
tau2_incoh_DUNE_tau_opt_integrand = np.multiply(flux_DUNE_tau_opt_norm, [sum(xsecs) for xsecs in zip(DUNE_tau_opt_xsec_2tau_p_matched, DUNE_tau_opt_xsec_2tau_n_matched)])

## Plotting ##
# 2mu Coherent ; Argon #
ax2[0,0].plot(energy_DUNE, mu2_coh_DUNE_integrand, c = c1, label = 'DUNE', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2[0,0].plot(energy_DUNE_tau_opt, mu2_coh_DUNE_tau_opt_integrand, c = c2, label = r'DUNE $\tau-$opt.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2[0,0].set_title(r'$\nu_\mu\mu^+\mu^-$ Coh.')

# 2mu Incoherent ; Argon #
ax2[1,0].plot(energy_DUNE, mu2_incoh_DUNE_integrand, c = c3, label = 'DUNE', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2[1,0].plot(energy_DUNE_tau_opt, mu2_incoh_DUNE_tau_opt_integrand, c = c4, label = r'DUNE $\tau-$opt.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2[1,0].set_title(r'$\nu_\mu\mu^+\mu^-$ Incoh.')

# 1tau Coherent ; Argon #
ax2[0,1].plot(energy_DUNE, tau1_coh_DUNE_integrand, c = c1, label = 'DUNE', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2[0,1].plot(energy_DUNE_tau_opt, tau1_coh_DUNE_tau_opt_integrand, c = c2, label = r'DUNE $\tau-$opt.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2[0,1].set_title(r'$\nu_\tau\tau^+\mu^-$ Coh.')

# 1tau Incoherent ; Argon #
ax2[1,1].plot(energy_DUNE, tau1_incoh_DUNE_integrand, c = c3, label = 'DUNE', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2[1,1].plot(energy_DUNE_tau_opt, tau1_incoh_DUNE_tau_opt_integrand, c = c4, label = r'DUNE $\tau-$opt.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2[1,1].set_title(r'$\nu_\tau\tau^+\mu^-$ Incoh.')

# 2tau Coherent ; Argon #
ax2[0,2].plot(energy_DUNE, tau2_coh_DUNE_integrand, c = c1, label = 'DUNE', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2[0,2].plot(energy_DUNE_tau_opt, tau2_coh_DUNE_tau_opt_integrand, c = c2, label = r'DUNE $\tau-$opt.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2[0,2].set_title(r'$\nu_\mu\tau^+\tau^-$ Coh.')

# 2tau Incoherent ; Argon #
ax2[1,2].plot(energy_DUNE, tau2_incoh_DUNE_integrand, c = c3, label = 'DUNE', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2[1,2].plot(energy_DUNE_tau_opt, tau2_incoh_DUNE_tau_opt_integrand, c = c4, label = r'DUNE $\tau-$opt.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2[1,2].set_title(r'$\nu_\mu\tau^+\tau^-$ Incoh.')

for ax1D in ax2:   # ax2 is a 2D array so ax will be a 1D array
    for ax in ax1D:
        ax.set_xlabel('Energy (GeV)')
        ax.set_ylabel(r'Norm. Flux $\times$ Cross Section')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid()
        ax.legend()

fig2.savefig("eventCalc_integrands.png", dpi=400)
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
