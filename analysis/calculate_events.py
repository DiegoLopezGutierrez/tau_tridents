import numpy as np
import csv
from scipy.integrate import simpson

#######################
###### Constants ######
#######################

atomic_mass_unit = 1.6605e-27 # kg

N_POT = 1.1e21               # Number of POTs
MD = 1000                    # Mass of the DUNE ND detector in kg (1 tonne)
MAr = 39.95*atomic_mass_unit # Mass of argon in atomic mass units

N_tonne = 147 # Should be 50 now for liquid argon.
N_year = 3

Phi_Alt = 1.04e-3

########################
###### Initialize ######
########################

### Directories ###
FLUX_DIR = '../csv/fluxes'
CROSS_SECTION_DIR = '../csv/cross_sections'

### Fluxes ###
DUNE_filename = FLUX_DIR + '/DUNE/DUNE_diff_flux.csv'  # dPhi/dE
Alt_filename = FLUX_DIR + '/Altmannshofer/vmu_normalized_flux_Altmannshofer_digitized.csv' # 1/Phi * dPhi/dE
#DUNE_FD_tau-opt_flux = 'csv/DUNE_tau-opt_FD_flux.csv'

energy_DUNE = []
flux_DUNE = []

energy_Alt = []
flux_Alt = []

#energy_DUNE_tau-opt = []
#diff_flux_DUNE_tau-opt = []

### Cross Sections ###
# Filenames #
xsec_1tau_coh_Ar_filename = CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/coherent/argon/vmu_to_vtau_tau+_mu-_coh_Ar_xsec.csv'
xsec_1tau_p_filename = CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv'

xsec_2tau_coh_Ar_filename = CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/coherent/argon/vmu_to_vmu_tau+_tau-_coh_Ar_xsec.csv'
xsec_2tau_p_filename = CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv'

xsec_2mu_coh_Ar_filename  = CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/argon/vmu_to_vmu_mu+_mu-_coh_Ar_xsec.csv'
xsec_2mu_p_filename  = CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/proton/vmu_to_vmu_mu+_mu-_nucleon_p_xsec.csv'

xsec_CC_filename = CROSS_SECTION_DIR + '/vmuCC/vmuCC_xsec_perE_Formaggio.csv'

# Coherent ; Argon #
energy_1tau_coh_Ar = []
xsec_1tau_coh_Ar = []

energy_2tau_coh_Ar = []
xsec_2tau_coh_Ar = []

energy_2mu_coh_Ar = []
xsec_2mu_coh_Ar = []

# Nucleon ; p #
energy_1tau_p = []
xsec_1tau_p = []

energy_2tau_p = []
xsec_2tau_p = []

energy_2mu_p = []
xsec_2mu_p = []



########################
###### Load Files ######
########################

##### Fluxes #####
### DUNE ; Standard Mode Flux ###
with open(DUNE_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_DUNE.append(float(row[0]))
        flux_DUNE.append(float(row[1]) * 1e4) # DUNE differential flux is in [cm^-2 GeV^-1 POT^-1]; convert to [m^-2 GeV^-1 POT^-1]

### DUNE ; Altmannshofer digitized ###
with open(Alt_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_Alt.append(float(row[0]))
        flux_Alt.append(float(row[1]))

#with open(DUNE_FD_tau-opt_flux,'r') as csvfile:
#    data = csv.reader(csvfile, delimiter = ',')
#    for row in data:
#        energy_DUNE_tau-opt.append(float(row[0]))
#        diff_flux_DUNE_tau-opt.append(float(row[1]) * 1e-21)

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

### Nucleon ###
# vmu -> vtau tau+ mu- ; nucleon ; proton #
with open(xsec_1tau_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_p.append(float(row[0]))
        xsec_1tau_p.append(float(row[1]) * 1e-43) # Convert to m^2.

# vmu -> vmu tau+ tau- ; nucleon ; proton #
with open(xsec_2tau_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_p.append(float(row[0]))
        xsec_2tau_p.append(float(row[1]) * 1e-43) # Convert to m^2.

# vmu -> vmu mu+ mu- ; nucleon ; proton #
with open(xsec_2mu_p_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_p.append(float(row[0]))
        xsec_2mu_p.append(float(row[1]) * 1e-43) # Convert to m^2.


#with open(xsec_cc_filename,'r') as csvfile:
#    data = csv.reader(csvfile, delimiter = ',')
#    for row in data:
#        energy_cc.append(float(row[0]))

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

DUNE_xsec_1tau_coh_Ar_matched = MatchXSec(energy_DUNE, energy_1tau_coh_Ar, xsec_1tau_coh_Ar)
DUNE_xsec_1tau_p_matched = MatchXSec(energy_DUNE, energy_1tau_p, xsec_1tau_p)

DUNE_xsec_2tau_coh_Ar_matched = MatchXSec(energy_DUNE, energy_2tau_coh_Ar, xsec_2tau_coh_Ar)
DUNE_xsec_2tau_p_matched = MatchXSec(energy_DUNE, energy_2tau_p, xsec_2tau_p)

DUNE_xsec_2mu_coh_Ar_matched = MatchXSec(energy_DUNE, energy_2mu_coh_Ar, xsec_2mu_coh_Ar)
DUNE_xsec_2mu_p_matched = MatchXSec(energy_DUNE, energy_2mu_p, xsec_2mu_p)

Alt_xsec_1tau_coh_Ar_matched = MatchXSec(energy_Alt, energy_1tau_coh_Ar, xsec_1tau_coh_Ar)
Alt_xsec_1tau_p_matched = MatchXSec(energy_Alt, energy_1tau_p, xsec_1tau_p)

Alt_xsec_2tau_coh_Ar_matched = MatchXSec(energy_Alt, energy_2tau_coh_Ar, xsec_2tau_coh_Ar)
Alt_xsec_2tau_p_matched = MatchXSec(energy_Alt, energy_2tau_p, xsec_2tau_p)

Alt_xsec_2mu_coh_Ar_matched = MatchXSec(energy_Alt, energy_2mu_coh_Ar, xsec_2mu_coh_Ar)
Alt_xsec_2mu_p_matched = MatchXSec(energy_Alt, energy_2mu_p, xsec_2mu_p)

# Match xsec to energy given by DUNE flux.
#for i in range(len(energy_DUNE)):
#    for j in range(len(energy_1tau)):
#        if (energy_DUNE[i] < energy_1tau[0]): # Check that DUNE energy val is within the energy range of 1tau trident energies
#            xsec_1tau_matched.append(0.0) # and set xsec to 0; otherwise, these values will be skipped altogether
#            break
#        elif ((energy_DUNE[i] >= energy_1tau[j]) and (energy_DUNE[i] <= energy_1tau[j+1])): # If the DUNE energy val is between two energy vals in 1tau trident
#            xsec_1tau_matched.append(xsec_1tau[j]) # Add the 1tau trident cross section of the lower of the two energy vals to the xsec array matched at the DUNE energy val.
#
#    for j in range(len(energy_2tau)):
#        if (energy_DUNE[i] < energy_2tau[0]):
#            xsec_2tau_matched.append(0.0)
#            break
#        elif ((energy_DUNE[i] >= energy_2tau[j]) and (energy_DUNE[i] <= energy_2tau[j+1])):
#            xsec_2tau_matched.append(xsec_2tau[j])
#
#    for j in range(len(energy_2mu)):
#        if (energy_DUNE[i] < energy_2mu[0]):
#            xsec_2mu_matched.append(0.0)
#            break
#        elif ((energy_DUNE[i] >= energy_2mu[j]) and (energy_DUNE[i] <= energy_2mu[j+1])):
#            xsec_2mu_matched.append(xsec_2mu[j])
#
#for i in range(len(energy_Alt)):
#    for j in range(len(energy_1tau)):
#        if (energy_Alt[i] < energy_1tau[0]):
#            xsec_1tau_matched_Alt.append(0.0)
#            break
#        elif ((energy_Alt[i] >= energy_1tau[j]) and (energy_Alt[i] <= energy_1tau[j+1])):
#            xsec_1tau_matched_Alt.append(xsec_1tau[j])
#
#    for j in range(len(energy_2tau)):
#        if (energy_Alt[i] < energy_2tau[0]):
#            xsec_2tau_matched_Alt.append(0.0)
#            break
#        elif ((energy_Alt[i] >= energy_2tau[j]) and (energy_Alt[i] <= energy_2tau[j+1])):
#            xsec_2tau_matched_Alt.append(xsec_2tau[j])
#
#    for j in range(len(energy_2mu)):
#        if (energy_Alt[i] < energy_2mu[0]):
#            xsec_2mu_matched_Alt.append(0.0)
#            break
#        elif ((energy_Alt[i] >= energy_2mu[j]) and (energy_Alt[i] <= energy_2mu[j+1])):
#            xsec_2mu_matched_Alt.append(xsec_2mu[j])

#print("Size of DUNE energy array: ", len(energy_DUNE))
#print("Size of Altmannshofer energy array: ", len(energy_Alt))

#print("Size of cross section array for 1 tau: ", len(xsec_1tau_matched))
#print("Size of cross section array for 2 tau: ", len(xsec_2tau_matched))
#print("Size of cross section array for 2 mu: ", len(xsec_2mu_matched))

#print("Size of cross section array for 1 tau (Alt): ", len(xsec_1tau_matched_Alt))
#print("Size of cross section array for 2 tau (Alt): ", len(xsec_2tau_matched_Alt))
#print("Size of cross section array for 2 mu (Alt): ", len(xsec_2mu_matched_Alt))

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

## DUNE Standard Flux Events ##
N_1tau_coh_Ar_DUNE_1_1, XSecCon_1tau_coh_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched)
N_1tau_coh_Ar_DUNE_50_3, XSecCon_1tau_coh_Ar_DUNE_50_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched, NTONNES=50, NYEAR=3)
N_1tau_coh_Ar_DUNE_147_3, XSecCon_1tau_coh_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_coh_Ar_matched, NTONNES=147, NYEAR=3)

N_2tau_coh_Ar_DUNE_1_1, XSecCon_2tau_coh_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched)
N_2tau_coh_Ar_DUNE_50_3, XSecCon_2tau_coh_Ar_DUNE_50_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched, NTONNES=50, NYEAR=3)
N_2tau_coh_Ar_DUNE_147_3, XSecCon_2tau_coh_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_coh_Ar_matched, NTONNES=147, NYEAR=3)

N_2mu_coh_Ar_DUNE_1_1, XSecCon_2mu_coh_Ar_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched)
N_2mu_coh_Ar_DUNE_50_3, XSecCon_2mu_coh_Ar_DUNE_50_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched, NTONNES=50, NYEAR=3)
N_2mu_coh_Ar_DUNE_147_3, XSecCon_2mu_coh_Ar_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_coh_Ar_matched, NTONNES=147, NYEAR=3)

N_1tau_p_DUNE_1_1, XSecCon_1tau_p_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_p_matched)
N_1tau_p_DUNE_50_3, XSecCon_1tau_p_DUNE_50_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_p_matched, NTONNES=50, NYEAR=3)
N_1tau_p_DUNE_147_3, XSecCon_1tau_p_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_1tau_p_matched, NTONNES=147, NYEAR=3)

N_2tau_p_DUNE_1_1, XSecCon_2tau_p_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_p_matched)
N_2tau_p_DUNE_50_3, XSecCon_2tau_p_DUNE_50_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_p_matched, NTONNES=50, NYEAR=3)
N_2tau_p_DUNE_147_3, XSecCon_2tau_p_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2tau_p_matched, NTONNES=147, NYEAR=3)

N_2mu_p_DUNE_1_1, XSecCon_2mu_p_DUNE_1_1 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_p_matched)
N_2mu_p_DUNE_50_3, XSecCon_2mu_p_DUNE_50_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_p_matched, NTONNES=50, NYEAR=3)
N_2mu_p_DUNE_147_3, XSecCon_2mu_p_DUNE_147_3 = CalculateEvents(flux_DUNE, energy_DUNE, DUNE_xsec_2mu_p_matched, NTONNES=147, NYEAR=3)

## Altmannshofer Flux Events ##
N_1tau_coh_Ar_Alt_1_1, XSecCon_1tau_coh_Ar_Alt_1_1 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_1tau_coh_Ar_matched, normalized=True, Phi=Phi_Alt)
N_1tau_coh_Ar_Alt_50_3, XSecCon_1tau_coh_Ar_Alt_50_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_1tau_coh_Ar_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_1tau_coh_Ar_Alt_147_3, XSecCon_1tau_coh_Ar_Alt_147_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_1tau_coh_Ar_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2tau_coh_Ar_Alt_1_1, XSecCon_2tau_coh_Ar_Alt_1_1 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2tau_coh_Ar_matched, normalized=True, Phi=Phi_Alt)
N_2tau_coh_Ar_Alt_50_3, XSecCon_2tau_coh_Ar_Alt_50_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2tau_coh_Ar_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2tau_coh_Ar_Alt_147_3, XSecCon_2tau_coh_Ar_Alt_147_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2tau_coh_Ar_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2mu_coh_Ar_Alt_1_1, XSecCon_2mu_coh_Ar_Alt_1_1 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2mu_coh_Ar_matched, normalized=True, Phi=Phi_Alt)
N_2mu_coh_Ar_Alt_50_3, XSecCon_2mu_coh_Ar_Alt_50_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2mu_coh_Ar_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2mu_coh_Ar_Alt_147_3, XSecCon_2mu_coh_Ar_Alt_147_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2mu_coh_Ar_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_1tau_p_Alt_1_1, XSecCon_1tau_p_Alt_1_1 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_1tau_p_matched, normalized=True, Phi=Phi_Alt)
N_1tau_p_Alt_50_3, XSecCon_1tau_p_Alt_50_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_1tau_p_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_1tau_p_Alt_147_3, XSecCon_1tau_p_Alt_147_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_1tau_p_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2tau_p_Alt_1_1, XSecCon_2tau_p_Alt_1_1 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2tau_p_matched, normalized=True, Phi=Phi_Alt)
N_2tau_p_Alt_50_3, XSecCon_2tau_p_Alt_50_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2tau_p_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2tau_p_Alt_147_3, XSecCon_2tau_p_Alt_147_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2tau_p_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

N_2mu_p_Alt_1_1, XSecCon_2mu_p_Alt_1_1 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2mu_p_matched, normalized=True, Phi=Phi_Alt)
N_2mu_p_Alt_50_3, XSecCon_2mu_p_Alt_50_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2mu_p_matched, NTONNES=50, NYEAR=3, normalized=True, Phi=Phi_Alt)
N_2mu_p_Alt_147_3, XSecCon_2mu_p_Alt_147_3 = CalculateEvents(flux_Alt, energy_Alt, Alt_xsec_2mu_p_matched, NTONNES=147, NYEAR=3, normalized=True, Phi=Phi_Alt)

def PrintOutEvent(filename, fluxname, tonnes, years, events, xsec_con):
    print("\t{} tonne Ar, {} year, {} flux: {:.3f} ; XSecCon: {}".format(tonnes, years, fluxname, events, xsec_con), file=filename)

def WriteOutFile(filename):
    with open(filename,'w') as textfile:
        print("Events expected at DUNE ND for:", file=textfile)
        print("1tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_1tau_coh_Ar_DUNE_1_1, XSecCon_1tau_coh_Ar_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_1tau_coh_Ar_DUNE_50_3, XSecCon_1tau_coh_Ar_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_1tau_coh_Ar_DUNE_147_3, XSecCon_1tau_coh_Ar_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_1tau_coh_Ar_Alt_1_1, XSecCon_1tau_coh_Ar_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_1tau_coh_Ar_Alt_50_3, XSecCon_1tau_coh_Ar_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_1tau_coh_Ar_Alt_147_3, XSecCon_1tau_coh_Ar_Alt_147_3)
        print("1tau Proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_1tau_p_DUNE_1_1, XSecCon_1tau_p_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_1tau_p_DUNE_50_3, XSecCon_1tau_p_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_1tau_p_DUNE_147_3, XSecCon_1tau_p_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_1tau_p_Alt_1_1, XSecCon_1tau_p_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_1tau_p_Alt_50_3, XSecCon_1tau_p_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_1tau_p_Alt_147_3, XSecCon_1tau_p_Alt_147_3)
        print("\n", file=textfile)
        print("2tau Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2tau_coh_Ar_DUNE_1_1, XSecCon_2tau_coh_Ar_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_2tau_coh_Ar_DUNE_50_3, XSecCon_2tau_coh_Ar_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2tau_coh_Ar_DUNE_147_3, XSecCon_2tau_coh_Ar_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_2tau_coh_Ar_Alt_1_1, XSecCon_2tau_coh_Ar_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_2tau_coh_Ar_Alt_50_3, XSecCon_2tau_coh_Ar_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_2tau_coh_Ar_Alt_147_3, XSecCon_2tau_coh_Ar_Alt_147_3)
        print("2tau Proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2tau_p_DUNE_1_1, XSecCon_2tau_p_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_2tau_p_DUNE_50_3, XSecCon_2tau_p_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2tau_p_DUNE_147_3, XSecCon_2tau_p_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_2tau_p_Alt_1_1, XSecCon_2tau_p_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_2tau_p_Alt_50_3, XSecCon_2tau_p_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_2tau_p_Alt_147_3, XSecCon_2tau_p_Alt_147_3)
        print("\n", file=textfile)
        print("2mu Coherent:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2mu_coh_Ar_DUNE_1_1, XSecCon_2mu_coh_Ar_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_2mu_coh_Ar_DUNE_50_3, XSecCon_2mu_coh_Ar_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2mu_coh_Ar_DUNE_147_3, XSecCon_2mu_coh_Ar_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_2mu_coh_Ar_Alt_1_1, XSecCon_2mu_coh_Ar_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_2mu_coh_Ar_Alt_50_3, XSecCon_2mu_coh_Ar_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_2mu_coh_Ar_Alt_147_3, XSecCon_2mu_coh_Ar_Alt_147_3)
        print("2mu Proton:", file=textfile)
        PrintOutEvent(textfile, "DUNE", 1, 1, N_2mu_p_DUNE_1_1, XSecCon_2mu_p_DUNE_1_1)
        PrintOutEvent(textfile, "DUNE", 50, 3, N_2mu_p_DUNE_50_3, XSecCon_2mu_p_DUNE_50_3)
        PrintOutEvent(textfile, "DUNE", 147, 3, N_2mu_p_DUNE_147_3, XSecCon_2mu_p_DUNE_147_3)
        PrintOutEvent(textfile, "Alt", 1, 1, N_2mu_p_Alt_1_1, XSecCon_2mu_p_Alt_1_1)
        PrintOutEvent(textfile, "Alt", 50, 3, N_2mu_p_Alt_50_3, XSecCon_2mu_p_Alt_50_3)
        PrintOutEvent(textfile, "Alt", 147, 3, N_2mu_p_Alt_147_3, XSecCon_2mu_p_Alt_147_3)
        print("\n", file=textfile)


WriteOutFile("number_of_events.txt")

#print("Events expected at DUNE ND for 1tau trident:")
#print("\t1 tonne Ar, 1 year, DUNE flux: ", N_1tau)
#print("\t147 tonne Ar, 3 year, DUNE flux: ", N_1tau * N_tonne * N_year)
#print("\t147 tonne Ar, 10 year, DUNE flux: ", N_1tau * N_tonne * 10)
#print("\t1 tonne Ar, 1 year, Altmannshofer flux: ", N_1tau_Alt)
#print("\t147 tonne Ar, 3 year, Altmannshofer flux: ", N_1tau_Alt * N_tonne * N_year)
#print("\t147 tonne Ar, 10 year, Altmannshofer flux: ", N_1tau_Alt * N_tonne * 10)
#
## Calculate the expected number of events
#Phi = simpson(diff_flux_DUNE, x=energy_DUNE)
#print("Integrated neutrino flux DUNE: ", Phi)
#print("Integrated neutrino flux Altmannshofer: ", Phi_Alt)
#
#sigma_1tau = np.multiply(diff_flux_DUNE, xsec_1tau_matched) / Phi
#sigma_1tau_Alt = np.multiply(norm_flux_Alt, xsec_1tau_matched_Alt)
##print("Size of convoluted sigma integrand for 1 tau: ", len(sigma_1tau))
##print("Size of convoluted sigma integrand for 1 tau (Alt): ", len(sigma_1tau_Alt))
#
#Sigma_1tau = simpson(sigma_1tau, x=energy_DUNE)
#Sigma_1tau_Alt = simpson(sigma_1tau_Alt, x=energy_Alt)
#
##print(sigma_1tau)
##print(Sigma_1tau)
#
#N_1tau = Phi * Sigma_1tau * MD * N_POT / MAr
#N_1tau_Alt = Phi_Alt * Sigma_1tau_Alt * MD * N_POT / MAr
#
#print("Events expected at DUNE ND for 1tau trident:")
#print("\t1 tonne Ar, 1 year, DUNE flux: ", N_1tau)
#print("\t147 tonne Ar, 3 year, DUNE flux: ", N_1tau * N_tonne * N_year)
#print("\t147 tonne Ar, 10 year, DUNE flux: ", N_1tau * N_tonne * 10)
#print("\t1 tonne Ar, 1 year, Altmannshofer flux: ", N_1tau_Alt)
#print("\t147 tonne Ar, 3 year, Altmannshofer flux: ", N_1tau_Alt * N_tonne * N_year)
#print("\t147 tonne Ar, 10 year, Altmannshofer flux: ", N_1tau_Alt * N_tonne * 10)
#
#sigma_2tau = np.multiply(diff_flux_DUNE, xsec_2tau_matched) / Phi
#sigma_2tau_Alt = np.multiply(norm_flux_Alt, xsec_2tau_matched_Alt)
##print("Size of convoluted sigma integrand for 2 tau: ", len(sigma_2tau))
##print("Size of convoluted sigma integrand for 2 tau (Alt): ", len(sigma_2tau_Alt))
#
#Sigma_2tau = simpson(sigma_2tau, x=energy_DUNE)
#Sigma_2tau_Alt = simpson(sigma_2tau_Alt, x=energy_Alt)
#
##print(sigma_2tau)
##print(Sigma_2tau)
#
#N_2tau = Phi * Sigma_2tau * MD * N_POT / MAr
#N_2tau_Alt = Phi_Alt * Sigma_2tau_Alt * MD * N_POT / MAr
#
#print("Events expected at DUNE ND for 2tau trident:")
#print("\t1 tonne Ar, 1 year, DUNE flux: ", N_2tau)
#print("\t147 tonne Ar, 3 year, DUNE flux: ", N_2tau * N_tonne * N_year)
#print("\t147 tonne Ar, 10 year, DUNE flux: ", N_2tau * N_tonne * 10)
#print("\t1 tonne Ar, 1 year, Altmannshofer flux: ", N_2tau_Alt)
#print("\t147 tonne Ar, 3 year, Altmannshofer flux: ", N_2tau_Alt * N_tonne * N_year)
#print("\t147 tonne Ar, 10 year, Altmannshofer flux: ", N_2tau_Alt * N_tonne * 10)
#
#
#sigma_2mu = np.multiply(diff_flux_DUNE, xsec_2mu_matched) / Phi
#sigma_2mu_Alt = np.multiply(norm_flux_Alt, xsec_2mu_matched_Alt)
##print("Size of convoluted sigma integrand for 2 mu: ", len(sigma_2mu))
##print("Size of convoluted sigma integrand for 2 mu (Alt): ", len(sigma_2mu_Alt))
#
#Sigma_2mu = simpson(sigma_2mu, x=energy_DUNE)
#Sigma_2mu_Alt = simpson(sigma_2mu_Alt, x=energy_Alt)
#
##print(sigma_2mu)
##print(Sigma_2mu)
#
#N_2mu = Phi * Sigma_2mu * MD * N_POT / MAr
#N_2mu_Alt = Phi_Alt * Sigma_2mu_Alt * MD * N_POT / MAr
#
#print("Events expected at DUNE ND for 2mu trident:")
#print("\t1 tonne Ar, 1 year, DUNE flux: ", N_2mu)
#print("\t147 tonne Ar, 3 year, DUNE flux: ", N_2mu * N_tonne * N_year)
#print("\t147 tonne Ar, 10 year, DUNE flux: ", N_2mu * N_tonne * 10)
#print("\t1 tonne Ar, 1 year, Altmannshofer flux: ", N_2mu_Alt)
#print("\t147 tonne Ar, 3 year, Altmannshofer flux: ", N_2mu_Alt * N_tonne * N_year)
#print("\t147 tonne Ar, 10 year, Altmannshofer flux: ", N_2mu_Alt * N_tonne * 10)
#
