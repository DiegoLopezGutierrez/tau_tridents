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

N_POT = 1.1e21               # Number of POTs
FASER_L = 150                # FASERv luminosity of 150 fb^-1
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
DUNE_filename = FLUX_DIR + '/DUNE/histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_LBNEND_globes_flux.txt'  # Nv / m^2 GeV POT yr
Alt_filename = FLUX_DIR + '/Altmannshofer/vmu_normalized_flux_Altmannshofer_digitized.csv' # 1/Phi * dPhi/dE
DUNE_tau_opt_numu_flux = FLUX_DIR + '/DUNE/histos_g4lbne_v3r5p4_QGSP_BERT_TauOptimized_neutrino_LBNEND_globes_flux.txt'  # Nv / m^2 GeV POT yr
Alt120_filename = FLUX_DIR + '/Altmannshofer/numu_flux_120.csv' # dPhi/Phi. Needs to divide by dE.
FASERvmu_filename = FLUX_DIR + '/FASERnu/vmu/FASERvmu.csv' # Nv / m^2 GeV fb^-1

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

########################
###### Load Files ######
########################

### DUNE ; Standard Mode Flux ###
with open(DUNE_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ' ')
    for row in data:
        energy_DUNE.append(float(row[0]))
        flux_DUNE.append(float(row[2])) # Histogram has units of [m^-2 GeV^-1 POT^-1]

### DUNE ; Tau-optimized Flux ###
with open(DUNE_tau_opt_numu_flux,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ' ')
    for row in data:
        energy_DUNE_tau_opt.append(float(row[0]))
        flux_DUNE_tau_opt.append(float(row[2])) # Histogram has units of [m^-2 GeV^-1 POT^-1]

### FASERvmu ; vmu flux ###
with open(FASERvmu_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        energy = float(row[0]) # Energy in GeV
        flux = float(row[1]) # Number of neutrinos [Nv / m^2 GeV fb^-1]
        energy_FASERvmu.append(energy)
        flux_FASERvmu.append(flux)

### Integrated Fluxes ###
DUNE_integrated_flux = simpson(flux_DUNE, x=energy_DUNE) # [Nv / m^2 POT yr]
DUNE_tau_opt_integrated_flux = simpson(flux_DUNE_tau_opt, x=energy_DUNE_tau_opt) # [Nv / m^2 POT yr]
FASER_integrated_flux = simpson(flux_FASERvmu, x=energy_FASERvmu) # [Nv / m^2 fb^-1]

flux_DUNE_norm = np.divide(flux_DUNE, DUNE_integrated_flux) # Normalize DUNE flux [1/GeV]
flux_DUNE_tau_opt_norm = np.divide(flux_DUNE_tau_opt, DUNE_tau_opt_integrated_flux) # Normalize DUNE tau-opt flux [1/GeV]
flux_FASERvmu_norm = np.divide(flux_FASERvmu, FASER_integrated_flux) # Normalize FASERv flux [1/GeV]

### Plot normalized fluxes ###
fig1, ax1 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)
fig2, ax2 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)

#ax1.hist(energy_DUNE, bins=bins_DUNE, weights=flux_DUNE_norm, histtype='stepfilled', label=r'DUNE', color='navy', alpha=0.5, lw=0.5)
#ax1.hist(energy_DUNE, bins=bins_DUNE, weights=flux_DUNE_norm, histtype='step', color='black', lw=2, alpha=1)

ax1.plot(energy_DUNE, flux_DUNE_norm, '-', color='orange', label='DUNE Std.', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax1.fill_between(energy_DUNE, 0, flux_DUNE_norm, color='orange', alpha=0.5, label='_hidden')
ax1.plot(energy_DUNE_tau_opt, flux_DUNE_tau_opt_norm, '-', color='firebrick', label=r'DUNE $\tau$-opt', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax1.fill_between(energy_DUNE_tau_opt, 0, flux_DUNE_tau_opt_norm, color='firebrick', alpha=0.5, label='_hidden')

ax1.plot(energy_FASERvmu, flux_FASERvmu_norm, '-', color='blueviolet', label=r'FASER$\nu$', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax1.fill_between(energy_FASERvmu, 0, flux_FASERvmu_norm, color='blueviolet', alpha=0.5, label='_hidden')

ax1.set_xlabel(r'\textbf{Neutrino Energy} $E_\nu$ (GeV)')
ax1.set_ylabel(r'$\frac{1}{\Phi}\frac{\mathrm{d}\Phi}{\mathrm{d}E_\nu}$ (GeV$^{-1}$)')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(0.3, 1e4)
#ax1.set_ylim(5e-4, 0.500)
ax1.set_ylim(1e-8, 0.500)
ax1.legend(loc='upper right')
ax1.set_title(r"$\nu_\mu$ \textbf{Normalized Flux}")

fig1.savefig("../plots/fluxes.png", dpi=400)
