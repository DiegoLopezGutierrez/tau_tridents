import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe
import matplotlib as mpl
from scipy.integrate import simpson

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

FLUX_DIR = '../csv/fluxes/FASERnu/vmubar/'

E_pion = []
E_kaon = []
E_hyperon = []
E_charm = []
E_bottom = []

flux_pion = []
flux_kaon = []
flux_hyperon = []
flux_charm = []
flux_bottom = []
flux_total = []

pion_filename = FLUX_DIR+'FASERvmubar_pion.csv'
kaon_filename = FLUX_DIR+'FASERvmubar_kaon.csv'
hyperon_filename = FLUX_DIR+'FASERvmubar_hyperon.csv'
charm_filename = FLUX_DIR+'FASERvmubar_charm.csv'
bottom_filename = FLUX_DIR+'FASERvmubar_bottom.csv'

normalization = 25*25*1e-4 # 25 x 25 cm^2 to m^2
normalization *= 150 # 150 fb^-1 luminosity

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

######## Plotting ########  
fig1, ax1 = plt.subplots(1, 1, figsize=(17,13), tight_layout=True)

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

ax1.set_xlabel(r'\textbf{Neutrino Energy} (GeV)')
ax1.set_ylabel(r'$\mathrm{d}\Phi$ ($N_\nu$ / m$^2$ fb$^{-1}$)')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylim(1e6,1e12)
ax1.set_xlim(1e1,5e4)
ax1.legend(loc='upper right')
ax1.set_title(r"\textbf{FASER}$\nu$ $\bar{\nu}_\mu$ \textbf{flux}", fontsize=40)

fig1.savefig("../plots/FASERnu_vmubar_component_fluxes.png", dpi=400)

### Export histogram ###
middle_energies = []
for i in range(len(bins_pion)-1):
    middle_point = (bins_pion[i+1]-bins_pion[i])/2 + bins_pion[i]
    middle_energies.append(middle_point)

print(len(bins_pion), len(middle_energies), len(hist_total))
print("Total entries: ", total_entries)

with open('FASERvmubar.csv','w',newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for energy, flux in zip(middle_energies, hist_total):
        writer.writerow([energy,flux])
