import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe
import matplotlib as mpl
from scipy.integrate import simpson

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

FLUX_DIR = '../csv/fluxes/FASERnu/vmu/'

E_pion = []
E_kaon = []
E_hyperon = []
E_charm = []
E_bottom = []

E_toni = []
bins_toni = []

E_fast_vmu = []
E_fast_vmubar = []

flux_pion = []
flux_kaon = []
flux_hyperon = []
flux_charm = []
flux_bottom = []
flux_total = []

flux_toni = []
flux_fast_vmu = []
flux_fast_vmubar = []

pion_filename = FLUX_DIR+'FASERvmu_pion.csv'
kaon_filename = FLUX_DIR+'FASERvmu_kaon.csv'
hyperon_filename = FLUX_DIR+'FASERvmu_hyperon.csv'
charm_filename = FLUX_DIR+'FASERvmu_charm.csv'
bottom_filename = FLUX_DIR+'FASERvmu_bottom.csv'

normalization = 25*25*1e-4 # 25 x 25 cm^2 to m^2
normalization *= 150 # 150 fb^-1 luminosity
#normalization = 1

toni_filename = FLUX_DIR+'../FASERvmu_Toni.csv'

fast_vmu_filename = FLUX_DIR + '../FastNeutrinoFluxSimulation/FASER_14.txt'
fast_vmubar_filename = FLUX_DIR + '../FastNeutrinoFluxSimulation/FASER_-14.txt'

with open(pion_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',') 
    for row in data:
        E = float(row[0])
        N = float(row[1])  # N neutrinos going in a 25 cm x 25 cm cross sectional area
        N /= normalization
        E_pion.append(E)
        flux_pion.append(N)

with open(kaon_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',') 
    for row in data:
        E = float(row[0])
        N = float(row[1])  # N neutrinos going in a 25 cm x 25 cm cross sectional area
        N /= normalization
        E_kaon.append(E)
        flux_kaon.append(N)

with open(hyperon_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',') 
    for row in data:
        E = float(row[0])
        N = float(row[1])  # N neutrinos going in a 25 cm x 25 cm cross sectional area
        N /= normalization
        E_hyperon.append(E)
        flux_hyperon.append(N)

with open(charm_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',') 
    for row in data:
        E = float(row[0])
        N = float(row[1])  # N neutrinos going in a 25 cm x 25 cm cross sectional area
        N /= normalization
        E_charm.append(E)
        flux_charm.append(N)

with open(bottom_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',') 
    for row in data:
        E = float(row[0])
        N = float(row[1])  # N neutrinos going in a 25 cm x 25 cm cross sectional area
        N /= normalization
        E_bottom.append(E)
        flux_bottom.append(N)

bins_FASERvmu = []

# Get bin edges. This is based on a visual inspection of the DUNE plot.
with open('../csv/fluxes/FASERnu/FASERvmu_bins.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        bins_FASERvmu.append(float(row[0]))

bins_FASERvmu.append(8150.00)
bins_FASERvmu = np.array(bins_FASERvmu)

total_entries_toni = 0

with open(toni_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        low_E = float(row[0])
        high_E = float(row[1])
        N = float(row[2])
        total_entries_toni += N
        bins_toni.append(low_E)
        E_toni.append(low_E+(high_E-low_E)/2)
        flux_toni.append(N)

with open(fast_vmu_filename,'r') as txtfile:
    for line in txtfile:
        data = line.split()
        energy = float(data[0])
        flux = float(data[1])    # Flux is v's assuming 150 fb^-1 and 25x25 cm^2
        E_fast_vmu.append(energy)
        flux_fast_vmu.append(flux)

with open(fast_vmubar_filename,'r') as txtfile:
    for line in txtfile:
        data = line.split()
        energy = float(data[0])
        flux = float(data[1])    # Flux is v's assuming 150 fb^-1 and 25x25 cm^2
        E_fast_vmubar.append(energy)
        flux_fast_vmubar.append(flux)

bins_toni.append(10000.0)

fig1, ax1 = plt.subplots(1, 1, figsize=(17,13), tight_layout=True)
fig2, ax2 = plt.subplots(1, 1, figsize=(15,13), tight_layout=True)

col_map = mpl.colormaps['Dark2'].resampled(6)
rgb_col = np.linspace(0,1,6)

c_pion = col_map(rgb_col[0])
c_kaon = col_map(rgb_col[1])
c_hyperon = col_map(rgb_col[2])
c_charm = col_map(rgb_col[3])
c_bottom = col_map(rgb_col[4])
c_total = col_map(rgb_col[5])

hist_pion, bins_pion, patches_pion = ax1.hist(E_pion, bins=bins_FASERvmu, weights=flux_pion, histtype='step', color=c_pion,lw=4, alpha=1, label='pion')
hist_kaon, bins_kaon, patches_kaon = ax1.hist(E_kaon, bins=bins_FASERvmu, weights=flux_kaon, histtype='step', color=c_kaon,lw=4, alpha=1, label='kaon')
hist_hyperon, bins_hyperon, patches_hyperon = ax1.hist(E_hyperon, bins=bins_FASERvmu, weights=flux_hyperon, histtype='step', color=c_hyperon,lw=4, alpha=1, label='hyperon')
hist_charm, bins_charm, patches_charm = ax1.hist(E_charm, bins=bins_FASERvmu, weights=flux_charm, histtype='step', color=c_charm,lw=4, alpha=1, label='charm')
hist_bottom, bins_bottom, patches_bottom = ax1.hist(E_bottom, bins=bins_FASERvmu, weights=flux_bottom, histtype='step', color=c_bottom,lw=4, alpha=1, label='bottom')

hist_total = []

total_entries = 0

for e_bin, N_pion, N_kaon, N_hyperon, N_charm, N_bottom in zip(bins_pion, hist_pion, hist_kaon, hist_hyperon, hist_charm, hist_bottom):
    N_total = (N_pion + N_kaon + N_hyperon + N_charm + N_bottom)
    total_entries += N_total
    hist_total.append(N_total)

ax1.hist(bins_pion[:-1], bins=bins_pion, weights=hist_total, histtype='step', color=c_total,lw=5, alpha=1, label='total')

# Comparison with other FASER flux
hist_toni, bins_toni, patches_toni = ax2.hist(E_toni, bins=bins_toni, weights=flux_toni, histtype='step', color='red',lw=4, alpha=1, label='Toni')
ax2.plot(E_fast_vmu, flux_fast_vmu, '-', label=r'$\nu_\mu$', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax2.plot(E_fast_vmubar, flux_fast_vmubar, '-', label=r'$\bar{\nu}_\mu$', path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
ax2.hist(bins_pion[:-1], bins=bins_pion, weights=[total*normalization for total in hist_total], histtype='step', color=c_total, lw=5, alpha=1, label='total')

ax1.set_xlabel(r'\textbf{Neutrino Energy} (GeV)')
#ax1.set_ylabel(r'$\frac{\mathrm{d}\Phi}{\mathrm{d}E}$ ($N_\nu$ / m$^{2}$ GeV fb$^{-1}$)')
ax1.set_ylabel(r'$\mathrm{d}\Phi$ ($N_\nu$ / m$^2$ fb$^{-1}$)')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylim(1e6,1e12)
ax1.set_xlim(1e1,5e4)
ax1.legend(loc='upper right')
ax1.set_title(r"\textbf{FASER}$\nu$ $\nu_\mu$ \textbf{flux}", fontsize=40)

fig1.savefig("../plots/FASERnu_component_fluxes.png", dpi=400)

ax2.set_xlabel(r'\textbf{Neutrino Energy} (GeV)')
ax2.set_ylabel(r'\textbf{Neutrinos} (1 / bin)')
ax2.set_xscale('log')
ax2.set_yscale('log')
#ax2.set_ylim(1e6,1e12)
#ax2.set_xlim(1e1,5e4)
ax2.legend(loc='upper right')
ax2.set_title(r"\textbf{FASER}$\nu$ $\nu_\mu$ \textbf{flux - Toni Makela}", fontsize=40)

fig2.savefig("../plots/FASERnu_Toni.png", dpi=400)

### Export histogram ###
middle_energies = []
for i in range(len(bins_pion)-1):
    middle_point = (bins_pion[i+1]-bins_pion[i])/2 + bins_pion[i]
    middle_energies.append(middle_point)

print(len(bins_pion), len(middle_energies), len(hist_total))
print("Total entries: ", total_entries)
print("Total entries - Toni: ", total_entries_toni)

with open('FASERvmu.csv','w',newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for energy, flux in zip(middle_energies, hist_total):
        writer.writerow([energy,flux])
