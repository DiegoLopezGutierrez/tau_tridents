import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe
import csv 
from scipy.integrate import simpson

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

#######################
###### Constants ######
#######################

atomic_mass_unit = 1.6605e-27 # kg

N_POT_DUNE = 1.1e21          # Number of POTs
L_FASER = 150                # FASERv luminosity of 150 fb^-1
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

# FASER #
energy_FASERvmu = []
bins_FASERvmu = []
flux_FASERvmu = []
flux_FASERvmubar = []
flux_FASERv_vmu_vmubar = []

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
flux_INGRID_total = []

# Get bin edges. This is based on a visual inspection of the DUNE plot.
with open(FLUX_DIR + '/FASERnu/FASERvmu_bins.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        bins_FASERvmu.append(float(row[0]))

bins_FASERvmu.append(8150.00)
bins_FASERvmu = np.array(bins_FASERvmu)

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
        flux = float(row[1]) # Flux has already been normalized by the 25 x 25 cm^2 area and 150 fb^-1; it is now in units of [m^-2 fb].
        energy_FASERvmu.append(energy)
        flux_FASERvmu.append(flux)

### FASERvmu ; vmubar flux ###
with open(FASERnubar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        flux = float(row[1]) # Flux has already been normalized by the 25 x 25 cm^2 area and 150 fb^-1; it is now in units of [m^-2 fb].
        flux_FASERvmubar.append(flux)

flux_FASERv_vmu_vmubar = [sum(flux) for flux in zip(flux_FASERvmu, flux_FASERvmubar)]

### FASERv2 ; vmu + vmubar flux ###
with open(FASERnu2_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        energy = float(row[0]) # Energy in GeV
        flux = float(row[1]) # Flux has not been normalized by the 0.5 x 0.5 m^2 area and 3 ab^-1; it is in units of [1].
        energy_FASERv2.append(energy)
        flux_FASERv2_vmu_vmubar.append(flux / (0.25 * 3000)) # Now flux is in units of [m^-2 fb]

### MINOS ; Neutrino Mode ; vmu flux ###
with open(MINOS_neutrino_vmu_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        energy = float(row[0]) # Energy in GeV
        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
        energy_MINOS.append(energy)
        flux_MINOS_neutrino_vmu.append(flux * 1e-20)

### MINOS ; Neutrino Mode ; vmubar flux ###
with open(MINOS_neutrino_vmubar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
        flux_MINOS_neutrino_vmubar.append(flux * 1e-20)

### MINOS ; Antineutrino Mode ; vmu flux ###
with open(MINOS_antineutrino_vmu_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        energy = float(row[0]) # Energy in GeV
        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
        flux_MINOS_antineutrino_vmu.append(flux * 1e-20)

### MINOS ; Antineutrino Mode ; vmubar flux ###
with open(MINOS_antineutrino_vmubar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
        flux_MINOS_antineutrino_vmubar.append(flux * 1e-20)

### MINOS+ ; Neutrino Mode ; vmu flux ###
with open(MINOSPlus_neutrino_vmu_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        energy = float(row[0]) # Energy in GeV
        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
        energy_MINOSPlus.append(energy)
        flux_MINOSPlus_neutrino_vmu.append(flux * 1e-20)

### MINOS+ ; Neutrino Mode ; vmubar flux ###
with open(MINOSPlus_neutrino_vmubar_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        flux = float(row[1]) # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1]
        flux_MINOSPlus_neutrino_vmubar.append(flux * 1e-20)

### T2K ; Neutrino Mode ; flux ###
with open(T2K_INGRID_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        energy = float(row[0])
        flux = float(row[1]) * 1e-20 # Histogram has units of [m^-2 GeV^-1 (1e20 POT)^-1].
        energy_INGRID.append(energy)
        flux_INGRID_neutrino_vmu.append(flux * 92.5/100.0)  # Flux composition of vmu is 92.5%. See Ballett et al.
        flux_INGRID_neutrino_vmubar.append(flux * 5.8/100.0) # Flux composition of vmu is 5.8%. See Ballett et al.
        flux_INGRID_neutrino_ve.append(flux * 1.5/100.0) # Flux composition of vmu is 1.5%. See Ballett et al.
        flux_INGRID_neutrino_vebar.append(flux * 0.2/100.0) # Flux composition of vmu is 0.2%. See Ballett et al.


########################
####### Plotting #######
########################

### Plot normalized fluxes ###
fig1, ax1 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  # DUNE and DUNE tau-optimized fluxes - Neutrino
fig2, ax2 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  # DUNE and DUNE tau-optimized fluxes - Antineutrino
fig3, ax3 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  # FASER flux
fig4, ax4 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  # MINOS and MINOS+ fluxes - Neutrino
fig5, ax5 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  # MINOS flux - Antineutrino
fig6, ax6 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  # T2K INGRID flux

## Colors ##
# https://davidmathlogic.com/colorblind/#%23FFA500-%238A2BE2-%23B22222-%231E3282-%2391B510-%232ca5e3-%23E70C64-%2327F596
cols = ['#FFA500',
        '#B22222',
        '#91B510',
        '#E70C64',
        '#8A2BE2',
        '#1E3282',
        '#2CA5E3',
        '#27F596']

col_vmu = '#FFA500'
col_vmubar = '#8A2BE2'
col_ve = '#B22222'
col_vebar = '#1E3282'

## DUNE ##
# Neutrino #
ax1.plot(energy_DUNE, flux_DUNE_neutrino_vmu, '-', color=col_vmu, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax1.plot(energy_DUNE, flux_DUNE_neutrino_vmubar, '-', color=col_vmubar, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax1.plot(energy_DUNE, flux_DUNE_neutrino_ve, '-', color=col_ve, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax1.plot(energy_DUNE, flux_DUNE_neutrino_vebar, '-', color=col_vebar, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax1.text(1,3.8e-4,r'\textbf{$\nu_\mu$ Std.}',color=col_vmu,rotation=0,fontsize=28)
ax1.text(2,1.9e-5,r'\textbf{$\bar{\nu}_\mu$ Std.}',color=col_vmubar,rotation=0,fontsize=28)
ax1.text(1,3.3e-6,r'\textbf{$\nu_e$ Std.}',color=col_ve,rotation=0,fontsize=28)
ax1.text(1,3.5e-7,r'\textbf{$\bar{\nu}_e$ Std.}',color=col_vebar,rotation=0,fontsize=28)

ax1.plot(energy_DUNE_tau_opt, flux_DUNE_tau_opt_neutrino_vmu, '--', color=col_vmu, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()], alpha=0.75)
ax1.plot(energy_DUNE_tau_opt, flux_DUNE_tau_opt_neutrino_vmubar, '--', color=col_vmubar, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax1.plot(energy_DUNE_tau_opt, flux_DUNE_tau_opt_neutrino_ve, '--', color=col_ve, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax1.plot(energy_DUNE_tau_opt, flux_DUNE_tau_opt_neutrino_vebar, '--', color=col_vebar, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax1.text(8,1.2e-4,r'\textbf{$\nu_\mu$ $\tau$-Opt.}',color=col_vmu,rotation=0,fontsize=28)
ax1.text(1,1e-5,r'\textbf{$\bar{\nu}_\mu$ $\tau$-Opt.}',color=col_vmubar,rotation=0,fontsize=28)
ax1.text(1,8.5e-7,r'\textbf{$\nu_e$ $\tau$-Opt.}',color=col_ve,rotation=0,fontsize=28)
ax1.text(1,1e-7,r'\textbf{$\bar{\nu}_e$ $\tau$-Opt.}',color=col_vebar,rotation=0,fontsize=28)

# Antineutrino #
ax2.plot(energy_DUNE, flux_DUNE_antineutrino_vmu, '-', color=col_vmu, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax2.plot(energy_DUNE, flux_DUNE_antineutrino_vmubar, '-', color=col_vmubar, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax2.plot(energy_DUNE, flux_DUNE_antineutrino_ve, '-', color=col_ve, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax2.plot(energy_DUNE, flux_DUNE_antineutrino_vebar, '-', color=col_vebar, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax2.text(1,3.8e-4,r'\textbf{$\bar{\nu}_\mu$ Std.}',color=col_vmubar,rotation=0,fontsize=28)
ax2.text(2,1.9e-5,r'\textbf{$\nu_\mu$ Std.}',color=col_vmu,rotation=0,fontsize=28)
ax2.text(1,3e-6,r'\textbf{$\bar{\nu}_e$ Std.}',color=col_vebar,rotation=0,fontsize=28)
ax2.text(2.8,3.3e-7,r'\textbf{$\nu_e$ Std.}',color=col_ve,rotation=0,fontsize=28)

ax2.plot(energy_DUNE_tau_opt, flux_DUNE_tau_opt_antineutrino_vmu, '--', color=col_vmu, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax2.plot(energy_DUNE_tau_opt, flux_DUNE_tau_opt_antineutrino_vmubar, '--', color=col_vmubar, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax2.plot(energy_DUNE_tau_opt, flux_DUNE_tau_opt_antineutrino_ve, '--', color=col_ve, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax2.plot(energy_DUNE_tau_opt, flux_DUNE_tau_opt_antineutrino_vebar, '--', color=col_vebar, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax2.text(8.5,8e-5,r'\textbf{$\bar{\nu}_\mu$ $\tau$-Opt.}',color=col_vmubar,rotation=0,fontsize=28)
ax2.text(0.3,0.8e-5,r'\textbf{$\nu_\mu$ $\tau$-Opt.}',color=col_vmu,rotation=0,fontsize=28)
ax2.text(1,7e-7,r'\textbf{$\bar{\nu}_e$ $\tau$-Opt.}',color=col_vebar,rotation=0,fontsize=28)
ax2.text(0.8,1.7e-7,r'\textbf{$\nu_e$ $\tau$-Opt.}',color=col_ve,rotation=0,fontsize=28)

## FASER ##
ax3.plot(energy_FASERvmu, flux_FASERv_vmu_vmubar, '-', color=col_vmu, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax3.text(20,1e10,r'\textbf{$\nu_\mu + \bar{\nu}_\mu$ FASER$\nu$}',color=col_vmu,rotation=0,fontsize=30)

## MINOS and MINOS+ ##
# Neutrino #
ax4.plot(energy_MINOS, flux_MINOS_neutrino_vmu, '-', color=col_vmu, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax4.plot(energy_MINOS, flux_MINOS_neutrino_vmubar, '-', color=col_vmubar, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax4.plot(energy_MINOSPlus, flux_MINOSPlus_neutrino_vmu, '--', color=col_vmu, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax4.plot(energy_MINOSPlus, flux_MINOSPlus_neutrino_vmubar, '--', color=col_vmubar, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax4.text(6,3.5e-6,r'\textbf{$\bar{\nu}_\mu$ MINOS}',color=col_vmubar,rotation=0,fontsize=30)
ax4.text(6.5,1e-5,r'\textbf{$\nu_\mu$ MINOS}',color=col_vmu,rotation=0,fontsize=30)
ax4.text(6,1e-6,r'\textbf{$\bar{\nu}_\mu$ MINOS+}',color=col_vmubar,rotation=0,fontsize=30)
ax4.text(11,1e-5,r'\textbf{$\nu_\mu$ MINOS+}',color=col_vmu,rotation=0,fontsize=30)

# Antineutrino #
ax5.plot(energy_MINOS, flux_MINOS_antineutrino_vmu, '-', color=col_vmu, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax5.plot(energy_MINOS, flux_MINOS_antineutrino_vmubar, '-', color=col_vmubar, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax5.text(5.5,1e-5,r'\textbf{$\bar{\nu}_\mu$ MINOS}',color=col_vmubar,rotation=0,fontsize=30)
ax5.text(10,4e-6,r'\textbf{$\nu_\mu$ MINOS}',color=col_vmu,rotation=0,fontsize=30)

## T2K INGRID ##
ax6.plot(energy_INGRID, flux_INGRID_neutrino_vmu, '-', color=col_vmu, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax6.plot(energy_INGRID, flux_INGRID_neutrino_vmubar, '-', color=col_vmubar, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax6.plot(energy_INGRID, flux_INGRID_neutrino_ve, '-', color=col_ve, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax6.plot(energy_INGRID, flux_INGRID_neutrino_vebar, '-', color=col_vebar, label='_hidden', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax6.text(2,1e-5,r'\textbf{$\bar{\nu}_\mu$}',color=col_vmubar,rotation=0,fontsize=35)
ax6.text(2,1.5e-4,r'\textbf{$\nu_\mu$}',color=col_vmu,rotation=0,fontsize=35)
ax6.text(2,3.5e-7,r'\textbf{$\bar{\nu}_e$}',color=col_vebar,rotation=0,fontsize=35)
ax6.text(2,2.5e-6,r'\textbf{$\nu_e$}',color=col_ve,rotation=0,fontsize=35)

## Styling ##
for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.set_xlabel(r'\textbf{Neutrino Energy} $E_\nu$ (GeV)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.legend(loc='upper right')
    ax.xaxis.grid(True, linestyle='--', which='major',color='grey', alpha=.45)
    ax.yaxis.grid(True, linestyle='--', which='major',color='grey', alpha=.45)

ax1.set_ylabel(r'$\frac{\mathrm{d}\Phi}{\mathrm{d}E_\nu}$ ($\nu$ / m$^2$ / GeV / POT)')
ax2.set_ylabel(r'$\frac{\mathrm{d}\Phi}{\mathrm{d}E_\nu}$ ($\nu$ / m$^2$ / GeV / POT)')
ax3.set_ylabel(r'$\Phi$ ($\nu$ / m$^2$ / fb$^{-1}$)')
ax4.set_ylabel(r'$\frac{\mathrm{d}\Phi}{\mathrm{d}E_\nu}$ ($\nu$ / m$^2$ / GeV / POT)')
ax5.set_ylabel(r'$\frac{\mathrm{d}\Phi}{\mathrm{d}E_\nu}$ ($\nu$ / m$^2$ / GeV / POT)')
ax6.set_ylabel(r'$\frac{\mathrm{d}\Phi}{\mathrm{d}E_\nu}$ ($\nu$ / m$^2$ / GeV / POT)')

ax1.set_title(r"\textbf{Flux at DUNE in $\nu$ Mode}")
ax2.set_title(r"\textbf{Flux at DUNE in $\bar{\nu}$ Mode}")
ax3.set_title(r"\textbf{Flux at FASER$\nu$}")
ax4.set_title(r"\textbf{Flux at MINOS and MINOS+ in $\nu$ Mode}")
ax5.set_title(r"\textbf{Flux at MINOS in $\bar{\nu}$ Mode}")
ax6.set_title(r"\textbf{Flux at T2K INGRID}")

ax1.set_ylim(4e-8,1e-3)
ax1.set_xlim(0.1, 100)
ax2.set_ylim(4e-8,1e-3)
ax2.set_xlim(0.1, 100)

ax4.set_xscale('linear')
ax5.set_xscale('linear')
ax6.set_xscale('linear')

fig1.savefig("../plots/flux_DUNE_neutrino.png", dpi=100)
fig2.savefig("../plots/flux_DUNE_antineutrino.png", dpi=100)
fig3.savefig("../plots/flux_FASER.png", dpi=100)
fig4.savefig("../plots/flux_MINOS_neutrino.png", dpi=100)
fig5.savefig("../plots/flux_MINOS_antineutrino.png", dpi=100)
fig6.savefig("../plots/flux_INGRID.png", dpi=100)
