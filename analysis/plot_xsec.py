import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import csv 

plt.rcParams['text.usetex'] = True

CROSS_SECTION_DIR = '../csv/cross_sections'

energy_1tau = [] 
xsec_1tau = []
delta_1tau = []

energy_2tau_ar = []
xsec_2tau_ar = []
delta_2tau_ar = []

energy_2tau_p = []
xsec_2tau_p = []
delta_2tau_p = []

energy_2mu = []
xsec_2mu = []
delta_2mu = []

energy_2mu_dig = []
xsec_2mu_dig = []

energy_numuCC = []
xsec_numuCC = []

with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/coherent/argon/vmu_to_vtau_tau+_mu-_xsec.csv','r') as csvfile: 
    data = csv.reader(csvfile, delimiter = ',') 
    for row in data:
        energy_1tau.append(float(row[0])) 
        xsec_1tau.append(float(row[1]))
        delta_1tau.append(float(row[2]))

with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/coherent/argon/vmu_to_vmu_tau+_tau-_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_ar.append(float(row[0]))
        xsec_2tau_ar.append(float(row[1]))
        delta_2tau_ar.append(float(row[2]))

with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/incoherent/proton/vmu_p_to_vmu_tau+_tau-_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_p.append(float(row[0]))
        xsec_2tau_p.append(float(row[1]))
        delta_2tau_p.append(float(row[2]))

with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/argon/vmu_to_vmu_mu+_mu-_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu.append(float(row[0]))
        xsec_2mu.append(float(row[1]))
        delta_2mu.append(float(row[2]))

with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/argon/vmu_to_vmu_mu+_mu-_xsec_digitized.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_dig.append(float(row[0]))
        xsec_2mu_dig.append(float(row[1]))

with open(CROSS_SECTION_DIR + '/vmuCC/vmuCC_xsec_perE.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) * energy / 10 # digitized plot is in [1e-38 cm^2 / GeV]; 1 fb = 1e-39 cm^2
        energy_numuCC.append(energy)
        xsec_numuCC.append(xsec)


fig, ax = plt.subplots(1, 1, tight_layout=True)

ax.errorbar(energy_1tau, xsec_1tau, yerr=delta_1tau, color = 'firebrick', label = r"$\nu_\mu$ Ar $\to \nu_\tau \tau^+ \mu^-$ Ar")
ax.errorbar(energy_2tau_ar, xsec_2tau_ar, yerr=delta_2tau_ar, color = 'c', label = r"$\nu_\mu$ Ar $\to \nu_\mu \tau^+ \tau^-$ Ar")
ax.errorbar(energy_2tau_p, xsec_2tau_p, yerr=delta_2tau_p, color = 'c', linestyle = 'dashed', label = r'$\nu_\mu p \to \nu_\mu \tau^+ \tau^- p$')
ax.errorbar(energy_2mu, xsec_2mu, yerr=delta_2mu, color = 'olivedrab', label = r'$\nu_\mu$ Ar $\to \nu_\mu \mu^+ \mu^-$ Ar')

ax.scatter(energy_2mu_dig, xsec_2mu_dig, color = 'goldenrod', s=5, label = r'$\nu_\mu$ Ar $\to \nu_\mu \mu^+ \mu^-$ Ar')
ax.scatter(energy_numuCC, xsec_numuCC, color = 'blueviolet', s=5, label = r'$\nu_\mu N \to \mu^- X$')

ax.set_xlabel('Energy (GeV)') 
ax.set_ylabel(r'$\sigma$ (fb)') 
ax.set_ylim(1e-9, 1e2)
ax.legend() 
ax.set_yscale('log')
ax.set_xscale('log')

xmajor = [1, 5, 10, 50, 100]
ax.set_xticks(xmajor, labels=['1', '5', '10', '50', '100'])
#ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
#ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
ax.grid(which='major', axis='both')

fig.savefig("../plots/xsec.png", dpi=400, bbox_inches='tight') 
