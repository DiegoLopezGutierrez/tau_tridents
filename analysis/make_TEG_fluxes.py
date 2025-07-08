import numpy as np
import csv 
from scipy.integrate import simpson

#####################################################################################################################################
# The TEG generator can store neutrino fluxes to calculate convoluted cross sections.
# The format of the stored fluxes is as follows: [E_low, E_high, P_norm]
# Here, E_low [GeV] and E_high [GeV] denote the bin's lower and upper end.
# P_norm [1] is the normalized probability distribution. To go from a flux [m^-2 GeV^-1] to P_norm, you must:
#     i) flux [m^-2 GeV^-1] -> flux * bin_width [m^-2]
#     ii) flux * bin_width [m^-2] -> flux * bin_width [m^-2] / total_flux [m^-2] = P_norm [1]
# If wanting to "plot" P_norm, you must first divide by the bin_widths. This renormalizes the probabilities to
# account for the variable bin widths and gives you the actual probability distribution function. You'll end up with the normalized 
# flux [GeV^-1]. Naturally, to recover the flux, you must multiply this normalized flux [GeV^-1] by the total integrated flux [m^-2].
#
# For FASERv, we have the flux in Phi(E) [m^-2 fb]. P_norm would just be:
#     i) flux [m^2 fb] -> flux [m^2 fb] / total_flux [m^2 fb] = P_norm [1]
#####################################################################################################################################

#######################
###### Constants ######
#######################

# Integrated fluxes [neutrinos / m^2 POT]
DUNE_neutrino_vmu = 1.044e-03
DUNE_neutrino_vmubar = 7.493e-05
DUNE_neutrino_ve = 1.082e-05
DUNE_neutrino_vebar = 2.352e-06

DUNE_tau_opt_neutrino_vmu = 1.149e-03
DUNE_tau_opt_neutrino_vmubar = 5.103e-05
DUNE_tau_opt_neutrino_ve = 8.776e-06
DUNE_tau_opt_neutrino_vebar = 2.057e-06

DUNE_antineutrino_vmu = 8.997e-05
DUNE_antineutrino_vmubar = 9.389e-04
DUNE_antineutrino_ve = 2.943e-06
DUNE_antineutrino_vebar = 8.905e-06

DUNE_tau_opt_antineutrino_vmu = 6.309e-05
DUNE_tau_opt_antineutrino_vmubar = 9.263e-04
DUNE_tau_opt_antineutrino_ve = 2.502e-06
DUNE_tau_opt_antineutrino_vebar = 5.808e-06

SBND_neutrino_vmu = 1.611e-04
SBND_neutrino_vmubar = 1.171e-05
SBND_neutrino_ve = 9.600e-07
SBND_neutrino_vebar = 9.926e-08

# FASERv integrated flux [neutrinos / m^2 fb^-1]
FASERv_vmu = 1.122e+12
FASERv_vmubar = 8.874e+11

########################
###### Initialize ######
########################

### Directories ###
FLUX_DIR = '../csv/fluxes'

### Fluxes ###
DUNE_neutrino_filename = FLUX_DIR + '/DUNE/histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_LBNEND_globes_flux.txt'
DUNE_antineutrino_filename = FLUX_DIR + '/DUNE/histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_antineutrino_LBNEND_globes_flux.txt'
DUNE_tau_opt_neutrino_filename = FLUX_DIR + '/DUNE/histos_g4lbne_v3r5p4_QGSP_BERT_TauOptimized_neutrino_LBNEND_globes_flux.txt'
DUNE_tau_opt_antineutrino_filename = FLUX_DIR + '/DUNE/histos_g4lbne_v3r5p4_QGSP_BERT_TauOptimized_antineutrino_LBNEND_globes_flux.txt'

SBND_filename = FLUX_DIR + '/SBND/BNB_SBND_flux.csv'

#MINOS_neutrino_vmu_filename = FLUX_DIR + '/MINOS/MINOS_ND/neutrino_mode/MINOS_ND_neutrino_vmu_flux.csv'
#MINOS_neutrino_vmubar_filename = FLUX_DIR + '/MINOS/MINOS_ND/neutrino_mode/MINOS_ND_neutrino_vmubar_flux.csv'
#MINOS_antineutrino_vmu_filename = FLUX_DIR + '/MINOS/MINOS_ND/antineutrino_mode/MINOS_ND_antineutrino_vmu_flux.csv'
#MINOS_antineutrino_vmubar_filename = FLUX_DIR + '/MINOS/MINOS_ND/antineutrino_mode/MINOS_ND_antineutrino_vmubar_flux.csv'
#
#MINOSPlus_neutrino_vmu_filename = FLUX_DIR + '/MINOS/MINOS+_ND/neutrino_mode/MINOS+_ND_neutrino_vmu_flux.csv'
#MINOSPlus_neutrino_vmubar_filename = FLUX_DIR + '/MINOS/MINOS+_ND/neutrino_mode/MINOS+_ND_neutrino_vmubar_flux.csv'
#MINOSPlus_antineutrino_vmu_filename = FLUX_DIR + '/MINOS/MINOS+_ND/antineutrino_mode/MINOS+_ND_antineutrino_vmu_flux.csv'
#MINOSPlus_antineutrino_vmubar_filename = FLUX_DIR + '/MINOS/MINOS+_ND/antineutrino_mode/MINOS+_ND_antineutrino_vmubar_flux.csv'
#
#T2K_INGRID_filename = FLUX_DIR + '/T2K/INGRID/T2K_INGRID_neutrino_flux.csv'

FASERv_vmu_filename = FLUX_DIR + '/FASERnu/vmu/FASERvmu.csv'
FASERv_vmubar_filename = FLUX_DIR + '/FASERnu/vmubar/FASERvmubar.csv'

# DUNE Standard #
energy_DUNE_low = []
energy_DUNE_high = []
flux_DUNE_neutrino_vmu = []
PDF_DUNE_neutrino_vmu = []
PDF_DUNE_neutrino_ve = []
PDF_DUNE_neutrino_vmubar = []
PDF_DUNE_neutrino_vebar = []

PDF_DUNE_antineutrino_vmu = []
PDF_DUNE_antineutrino_ve = []
PDF_DUNE_antineutrino_vmubar = []
PDF_DUNE_antineutrino_vebar = []

# DUNE Tau Optimized #
energy_DUNE_tau_opt_low = []
energy_DUNE_tau_opt_high = []
PDF_DUNE_tau_opt_neutrino_vmu = []
PDF_DUNE_tau_opt_neutrino_ve = []
PDF_DUNE_tau_opt_neutrino_vmubar = []
PDF_DUNE_tau_opt_neutrino_vebar = []

PDF_DUNE_tau_opt_antineutrino_vmu = []
PDF_DUNE_tau_opt_antineutrino_ve = []
PDF_DUNE_tau_opt_antineutrino_vmubar = []
PDF_DUNE_tau_opt_antineutrino_vebar = []

# SBND #
energy_SBND_low = []
energy_SBND_high = []

PDF_SBND_neutrino_vmu = []
PDF_SBND_neutrino_vmubar = []
PDF_SBND_neutrino_ve = []
PDF_SBND_neutrino_vebar = []

# MINOS #
#energy_MINOS = []
#flux_MINOS_neutrino_vmu = []
#flux_MINOS_neutrino_vmubar = []
#
#flux_MINOS_antineutrino_vmu = []
#flux_MINOS_antineutrino_vmubar = []
#
## MINOS+ #
#energy_MINOSPlus = []
#flux_MINOSPlus_neutrino_vmu = []
#flux_MINOSPlus_neutrino_vmubar = []
#
#flux_MINOSPlus_antineutrino_vmu = []
#flux_MINOSPlus_antineutrino_vmubar = []
#
## T2K - INGRID #
#energy_INGRID = []
#flux_INGRID_neutrino_vmu = []
#flux_INGRID_neutrino_vmubar = []
#flux_INGRID_neutrino_ve = []
#flux_INGRID_neutrino_vebar = []
#flux_INGRID_total = []

# FASERv #
energy_FASERv_low = []
energy_FASERv_high = []

PDF_FASERv_vmu = []
PDF_FASERv_vmubar = []

########################
###### Load Files ######
########################

def flux_to_PDF(flux, bin_width, total_flux):
    return flux*bin_width/total_flux

##### Fluxes #####

# DUNE Std and DUNE Tau-Opt fluxes from Laura Fields have a fixed bin width of 0.250 GeV starting at 0 GeV up to 125.5 GeV
# The energy value from her files is the average between the lower and upper end of the bin.

### DUNE ; Standard Mode Flux ; Neutrino Mode ###
with open(DUNE_neutrino_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ' ')
    for row in data:
        energy = float(row[0])
        elow = energy - 0.125
        ehigh = energy + 0.125
        energy_DUNE_low.append(elow)
        energy_DUNE_high.append(ehigh)
        flux_DUNE_neutrino_vmu.append(float(row[2]))
        PDF_DUNE_neutrino_vmu.append(flux_to_PDF(float(row[2]), 0.250, DUNE_neutrino_vmu)) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        PDF_DUNE_neutrino_ve.append(flux_to_PDF(float(row[1]), 0.250, DUNE_neutrino_ve)) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        PDF_DUNE_neutrino_vmubar.append(flux_to_PDF(float(row[5]), 0.250, DUNE_neutrino_vmubar)) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        PDF_DUNE_neutrino_vebar.append(flux_to_PDF(float(row[4]), 0.250, DUNE_neutrino_vebar)) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]

### DUNE ; Standard Mode Flux ; Antineutrino Mode ###
with open(DUNE_antineutrino_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ' ')
    for row in data:
        PDF_DUNE_antineutrino_vmu.append(float(row[2])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        PDF_DUNE_antineutrino_ve.append(float(row[1])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        PDF_DUNE_antineutrino_vmubar.append(float(row[5])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        PDF_DUNE_antineutrino_vebar.append(float(row[4])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]

### DUNE ; Tau-optimized Flux ; Neutrino Mode ###
with open(DUNE_tau_opt_neutrino_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ' ')
    for row in data:
        energy = float(row[0])
        elow = energy - 0.125
        ehigh = energy + 0.125
        energy_DUNE_tau_opt_low.append(elow)
        energy_DUNE_tau_opt_high.append(ehigh)
        PDF_DUNE_tau_opt_neutrino_vmu.append(flux_to_PDF(float(row[2]), 0.250, DUNE_tau_opt_neutrino_vmu)) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        PDF_DUNE_tau_opt_neutrino_ve.append(flux_to_PDF(float(row[1]), 0.250, DUNE_tau_opt_neutrino_ve)) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        PDF_DUNE_tau_opt_neutrino_vmubar.append(flux_to_PDF(float(row[5]), 0.250, DUNE_tau_opt_neutrino_vmubar)) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        PDF_DUNE_tau_opt_neutrino_vebar.append(flux_to_PDF(float(row[4]), 0.250, DUNE_tau_opt_neutrino_vebar)) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]

### DUNE ; Tau-optimized Flux ; Antineutrino Mode ###
with open(DUNE_tau_opt_antineutrino_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ' ')
    for row in data:
        PDF_DUNE_tau_opt_antineutrino_vmu.append(float(row[2])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        PDF_DUNE_tau_opt_antineutrino_ve.append(float(row[1])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        PDF_DUNE_tau_opt_antineutrino_vmubar.append(float(row[5])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]
        PDF_DUNE_tau_opt_antineutrino_vebar.append(float(row[4])) # Histogram has units of [m^-2 GeV^-1 POT^-1 yr^-1]

### SBND ###
with open(SBND_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter=',')
    for row in data:
        elow = float(row[0])
        ehigh = float(row[1])
        bin_width = ehigh - elow
        energy_SBND_low.append(elow)
        energy_SBND_high.append(ehigh)
        PDF_SBND_neutrino_vmu.append(flux_to_PDF(float(row[2]) / 50.0 / 1e6 * 1e3, bin_width, SBND_neutrino_vmu))
        PDF_SBND_neutrino_ve.append(flux_to_PDF(float(row[4]) / 50.0 / 1e6 * 1e3, bin_width, SBND_neutrino_ve))
        PDF_SBND_neutrino_vmubar.append(flux_to_PDF(float(row[3]) / 50.0 / 1e6 * 1e3, bin_width, SBND_neutrino_vmubar))
        PDF_SBND_neutrino_vebar.append(flux_to_PDF(float(row[5]) / 50.0 / 1e6 * 1e3, bin_width, SBND_neutrino_vebar))

### FASERv ###
with open(FASERv_vmu_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter=',')
    for row in data:
        elow = float(row[0])
        ehigh = float(row[1])
        bin_width = ehigh - elow
        flux = float(row[3])
        energy_FASERv_low.append(elow)
        energy_FASERv_high.append(ehigh)
        PDF_FASERv_vmu.append(flux_to_PDF(flux, 1.0, FASERv_vmu))

### FASERv ###
with open(FASERv_vmubar_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter=',')
    for row in data:
        flux = float(row[3])
        PDF_FASERv_vmubar.append(flux_to_PDF(flux, 1.0, FASERv_vmubar))

### MINOS ; Neutrino Mode ; vmu flux ###
#with open(MINOS_neutrino_vmu_filename,'r') as csvfile:
#    data = csv.reader(csvfile, delimiter=',')
#    for row in data:
#        energy = float(row[0]) # Energy in GeV
#        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
#        energy_MINOS.append(energy)
#        flux_MINOS_neutrino_vmu.append(flux * 1e-20)
#
#### MINOS ; Neutrino Mode ; vmubar flux ###
#with open(MINOS_neutrino_vmubar_filename,'r') as csvfile:
#    data = csv.reader(csvfile, delimiter=',')
#    for row in data:
#        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
#        flux_MINOS_neutrino_vmubar.append(flux * 1e-20)
#
#### MINOS ; Antineutrino Mode ; vmu flux ###
#with open(MINOS_antineutrino_vmu_filename,'r') as csvfile:
#    data = csv.reader(csvfile, delimiter=',')
#    for row in data:
#        energy = float(row[0]) # Energy in GeV
#        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
#        flux_MINOS_antineutrino_vmu.append(flux * 1e-20)
#
#### MINOS ; Antineutrino Mode ; vmubar flux ###
#with open(MINOS_antineutrino_vmubar_filename,'r') as csvfile:
#    data = csv.reader(csvfile, delimiter=',')
#    for row in data:
#        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
#        flux_MINOS_antineutrino_vmubar.append(flux * 1e-20)
#
#### MINOS+ ; Neutrino Mode ; vmu flux ###
#with open(MINOSPlus_neutrino_vmu_filename,'r') as csvfile:
#    data = csv.reader(csvfile, delimiter=',')
#    for row in data:
#        energy = float(row[0]) # Energy in GeV
#        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
#        energy_MINOSPlus.append(energy)
#        flux_MINOSPlus_neutrino_vmu.append(flux * 1e-20)
#
#### MINOS+ ; Neutrino Mode ; vmubar flux ###
#with open(MINOSPlus_neutrino_vmubar_filename,'r') as csvfile:
#    data = csv.reader(csvfile, delimiter=',')
#    for row in data:
#        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
#        flux_MINOSPlus_neutrino_vmubar.append(flux * 1e-20)
#
#### T2K ; Neutrino Mode ; flux ###
#with open(T2K_INGRID_filename,'r') as csvfile:
#    data = csv.reader(csvfile, delimiter=',')
#    for row in data:
#        energy = float(row[0])
#        flux = float(row[1]) * 1e-20 # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1].
#        energy_INGRID.append(energy)
#        flux_INGRID_neutrino_vmu.append(flux * 92.5/100.0)  # Flux composition of vmu is 92.5%. See Ballett et al.
#        flux_INGRID_neutrino_vmubar.append(flux * 5.8/100.0) # Flux composition of vmu is 5.8%. See Ballett et al.
#        flux_INGRID_neutrino_ve.append(flux * 1.5/100.0) # Flux composition of vmu is 1.5%. See Ballett et al.
#        flux_INGRID_neutrino_vebar.append(flux * 0.2/100.0) # Flux composition of vmu is 0.2%. See Ballett et al.


########################
####### Exporting ######
########################

def export_TEGflux(filename,energy_low,energy_high,PDF):
    with open(filename,'w',newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        for elow, ehigh, pdf in zip(energy_low, energy_high, PDF):
            writer.writerow([elow, ehigh, pdf])

EXPORT_FLUX_DIR = '../csv/fluxes/TEG_fluxes'

export_TEGflux(EXPORT_FLUX_DIR+'/SBND_neutrino_vmu.csv',energy_SBND_low,energy_SBND_high,PDF_SBND_neutrino_vmu)
export_TEGflux(EXPORT_FLUX_DIR+'/SBND_neutrino_ve.csv',energy_SBND_low,energy_SBND_high,PDF_SBND_neutrino_ve)
export_TEGflux(EXPORT_FLUX_DIR+'/SBND_neutrino_vmubar.csv',energy_SBND_low,energy_SBND_high,PDF_SBND_neutrino_vmubar)
export_TEGflux(EXPORT_FLUX_DIR+'/SBND_neutrino_vebar.csv',energy_SBND_low,energy_SBND_high,PDF_SBND_neutrino_vebar)

export_TEGflux(EXPORT_FLUX_DIR+'/DUNE_neutrino_vmu.csv',energy_DUNE_low,energy_DUNE_high,PDF_DUNE_neutrino_vmu)
export_TEGflux(EXPORT_FLUX_DIR+'/DUNE_neutrino_ve.csv',energy_DUNE_low,energy_DUNE_high,PDF_DUNE_neutrino_ve)
export_TEGflux(EXPORT_FLUX_DIR+'/DUNE_neutrino_vmubar.csv',energy_DUNE_low,energy_DUNE_high,PDF_DUNE_neutrino_vmubar)
export_TEGflux(EXPORT_FLUX_DIR+'/DUNE_neutrino_vebar.csv',energy_DUNE_low,energy_DUNE_high,PDF_DUNE_neutrino_vebar)

export_TEGflux(EXPORT_FLUX_DIR+'/DUNE_tau_opt_neutrino_vmu.csv',energy_DUNE_tau_opt_low,energy_DUNE_tau_opt_high,PDF_DUNE_tau_opt_neutrino_vmu)
export_TEGflux(EXPORT_FLUX_DIR+'/DUNE_tau_opt_neutrino_ve.csv',energy_DUNE_tau_opt_low,energy_DUNE_tau_opt_high,PDF_DUNE_tau_opt_neutrino_ve)
export_TEGflux(EXPORT_FLUX_DIR+'/DUNE_tau_opt_neutrino_vmubar.csv',energy_DUNE_tau_opt_low,energy_DUNE_tau_opt_high,PDF_DUNE_tau_opt_neutrino_vmubar)
export_TEGflux(EXPORT_FLUX_DIR+'/DUNE_tau_opt_neutrino_vebar.csv',energy_DUNE_tau_opt_low,energy_DUNE_tau_opt_high,PDF_DUNE_tau_opt_neutrino_vebar)

export_TEGflux(EXPORT_FLUX_DIR+'/FASERv_vmu.csv',energy_FASERv_low,energy_FASERv_high,PDF_FASERv_vmu)
export_TEGflux(EXPORT_FLUX_DIR+'/FASERv_vmubar.csv',energy_FASERv_low,energy_FASERv_high,PDF_FASERv_vmubar)
