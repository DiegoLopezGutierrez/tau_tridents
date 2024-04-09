import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe
import csv 
from scipy.integrate import simpson

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

FLUXES_DIR = '../csv/fluxes'

DUNE_filename = FLUXES_DIR+'/DUNE/DUNE_diff_flux.csv'
Alt_numu80_filename = FLUXES_DIR+'/Altmannshofer/numu_flux_80_ref.csv'
Alt_numu120_filename = FLUXES_DIR+'/Altmannshofer/numu_flux_120.csv'
Alt_numu_digitized_filename = FLUXES_DIR+'/Altmannshofer/vmu_normalized_flux_Altmannshofer_digitized.csv'

energy_DUNE = []
flux_DUNE = []

energy_Alt_numu80 = []
bins_Alt_numu80 = []
flux_Alt_numu80 = []

energy_Alt_numu120 = []
bins_Alt_numu120 = []
flux_Alt_numu120 = []

energy_Alt_numu_digitized = []
flux_Alt_numu_digitized = []

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

flux_DUNE_norm = np.divide(flux_DUNE, simpson(flux_DUNE, x=energy_DUNE))

### DUNE ; Altmannshofer numu 80 GeV ###
with open(Alt_numu80_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for i, row in enumerate(data):
        energy_Alt_numu80.append(float(row[0]))
        bins_Alt_numu80.append(float(row[0]))
        flux_Alt_numu80.append(float(row[2]))  # Unsure if this is normalized.
        #if i == len(data):
        #    bins_Alt_numu80.append(float(row[1])) # Add right bin edge for last line
    #final_row = data.readlines()[-1]
    #bins_Alt_numu80.append(float(final_row[1]))

bins_Alt_numu80.append(68) # ad hoc solution to add last bin edge.

### DUNE ; Altmannshofer numu 120 GeV ###
with open(Alt_numu120_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for i, row in enumerate(data):
        low_bin = float(row[0])
        high_bin = float(row[1])
        energy = (low_bin+high_bin)/2
        energy_Alt_numu120.append(energy)
        bins_Alt_numu120.append(float(row[0]))
        flux_Alt_numu120.append(float(row[2]) / (high_bin - low_bin))  # Unsure if this is normalized.
        #if i == len(data):
        #    bins_Alt_numu80.append(float(row[1])) # Add right bin edge for last line
    #final_row = data.readlines()[-1]
    #bins_Alt_numu80.append(float(final_row[1]))

bins_Alt_numu120.append(68.0) # ad hoc solution to add last bin edge.
integrated_flux = simpson(flux_Alt_numu120, x=energy_Alt_numu120)
flux_Alt_numu120 = np.divide(flux_Alt_numu120, integrated_flux)

### DUNE ; Altmannshofer digitized ###
with open(Alt_numu_digitized_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_Alt_numu_digitized.append(float(row[0]))
        flux_Alt_numu_digitized.append(float(row[1]))    # Normalized flux

# Get bin edges. This is based on a visual inspection of the DUNE plot.
bin_edges = []
with open('bin_edges_120.txt') as txt_file:
    for row in txt_file:
        bin_edges.append(float(row))

bin_edges = np.array(bin_edges)

print("Alt. 120 bins:\n\t", bins_Alt_numu120)
print("Size: ", len(bins_Alt_numu120))
print("Alt. 120 energy:\n\t", energy_Alt_numu120)
print("Size: ", len(energy_Alt_numu120))
print("Alt. 120 flux:\n\t", flux_Alt_numu120)
print("Size: ", len(flux_Alt_numu120))
print("Integrated flux: ", integrated_flux)

fig1, ax1 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Normalized fluxes

ax1.hist(energy_DUNE, bins=bin_edges, weights=flux_DUNE_norm, histtype='stepfilled', label=r'DUNE', color='navy', alpha=0.5, lw=0.5)
ax1.hist(energy_DUNE, bins=bin_edges, weights=flux_DUNE_norm, histtype='step', color='black',lw=2, alpha=1)
#ax1.hist(energy_Alt_numu80, bins=bins_Alt_numu80, weights=flux_Alt_numu80, histtype='stepfilled', label=r'Alt. $\nu_\mu$ 80', color='blueviolet', alpha=0.5)
ax1.hist(energy_Alt_numu120, bins=bins_Alt_numu120, weights=flux_Alt_numu120, histtype='stepfilled', label=r'Alt. $\nu_\mu$ 120', color='green', alpha=0.3)
ax1.hist(energy_Alt_numu120, bins=bins_Alt_numu120, weights=flux_Alt_numu120, histtype='step', color='black',lw=2, alpha=1)

ax1.plot(energy_Alt_numu_digitized, flux_Alt_numu_digitized, '-', color='orange', label='Alt. digitized', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()],)

ax1.set_xlabel('Energy (GeV)')
ax1.set_ylabel(r'$\frac{1}{\Phi}\frac{\mathrm{d}\Phi}{\mathrm{d}E}$ (GeV$^{-1}$)')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(0.3,20)
ax1.set_ylim(5e-4, 0.500)
ax1.legend(loc='upper right')

fig1.savefig("fluxes.png")
