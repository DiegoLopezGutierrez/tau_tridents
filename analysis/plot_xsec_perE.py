import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import csv 

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

CROSS_SECTION_DIR = '../csv/cross_sections'

###########################
###### Coherent xsec ######
###########################

energy_1tau_coh_Ar = [] 
xsec_1tau_coh_Ar = []
delta_1tau_coh_Ar = []

energy_2tau_coh_Ar = []
xsec_2tau_coh_Ar = []
delta_2tau_coh_Ar = []

energy_2mu_coh_Ar = []
xsec_2mu_coh_Ar = []
delta_2mu_coh_Ar = []

############################
##### Incoherent xsec ######
############################

### Proton ###
energy_2tau_p_Ar = []
xsec_2tau_p_Ar = []
delta_2tau_p_Ar = []

energy_1tau_p_Ar = []
xsec_1tau_p_Ar = []
delta_1tau_p_Ar = []

energy_2mu_p_Ar = []
xsec_2mu_p_Ar = []
delta_2mu_p_Ar = []

### Neutron ###
energy_2tau_n_Ar = []
xsec_2tau_n_Ar = []
delta_2tau_n_Ar = []

energy_1tau_n_Ar = []
xsec_1tau_n_Ar = []
delta_1tau_n_Ar = []

energy_2mu_n_Ar = []
xsec_2mu_n_Ar = []
delta_2mu_n_Ar = []

### Proton + Neutron ###
energy_2tau_incoh_Ar = []
xsec_2tau_incoh_Ar = []
delta_2tau_incoh_Ar = []

energy_1tau_incoh_Ar = []
xsec_1tau_incoh_Ar = []
delta_1tau_incoh_Ar = []

energy_2mu_incoh_Ar = []
xsec_2mu_incoh_Ar = []
delta_2mu_incoh_Ar = []

##########################
###### Nucleon xsec ######
##########################

### Proton ###
energy_1tau_p = []
xsec_1tau_p = []
delta_1tau_p = []

energy_2tau_p = []
xsec_2tau_p = []
delta_2tau_p = []

energy_2mu_p = []
xsec_2mu_p = []
delta_2mu_p = []

### Neutron ###
energy_2mu_n = []
xsec_2mu_n = []
delta_2mu_n = []

energy_1tau_n = []
xsec_1tau_n = []
delta_1tau_n = []

energy_2tau_n = []
xsec_2tau_n = []
delta_2tau_n = []

#########################
###### Other xsec #######
#########################

energy_2mu_coh_Ar_Alt = []
xsec_2mu_coh_Ar_Alt = []

energy_2mu_incoh_p_Ar_Alt = []
xsec_2mu_incoh_p_Ar_Alt = []

energy_2mu_incoh_n_Ar_Alt = []
xsec_2mu_incoh_n_Ar_Alt = []

energy_numuCC_FZ = []
xsec_numuCC_FZ = []


###############################
##### Load CSV xsec files #####
###############################
# Cross section / energy units will be in [1e-38 cm^2 / GeV]

##### Coherent #####
### vmu -> vtau tau+ mu- ; coherent ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/coherent/new_argon/vmu_to_vtau_tau+_mu-_coh_Ar_xsec.csv','r') as csvfile: 
    data = csv.reader(csvfile, delimiter = ',') 
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) / energy * 10
        delta = float(row[2]) / energy * 10

        energy_1tau_coh_Ar.append(energy) 
        xsec_1tau_coh_Ar.append(xsec)
        delta_1tau_coh_Ar.append(delta)

### vmu -> vmu tau+ tau- ; coherent ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/coherent/argon/vmu_to_vmu_tau+_tau-_coh_Ar_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) / energy * 10
        delta = float(row[2]) / energy * 10

        energy_2tau_coh_Ar.append(energy)
        xsec_2tau_coh_Ar.append(xsec)
        delta_2tau_coh_Ar.append(delta)

### vmu -> vmu mu+ mu- ; coherent ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/argon/vmu_to_vmu_mu+_mu-_coh_Ar_xsec_1em2_220.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) / energy * 10
        delta = float(row[2]) / energy * 10

        energy_2mu_coh_Ar.append(energy)
        xsec_2mu_coh_Ar.append(xsec)
        delta_2mu_coh_Ar.append(delta)

##### Incoherent #####
### vmu -> vmu tau+ tau- ; incoherent ; proton ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) * 18 / energy * 10
        delta = float(row[2]) * 18 / energy * 10

        energy_2tau_p_Ar.append(energy)
        xsec_2tau_p_Ar.append(xsec)
        delta_2tau_p_Ar.append(delta)

### vmu -> vtau tau+ mu- ; incoherent ; proton ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_p_Ar.append(float(row[0]))
        xsec_1tau_p_Ar.append(float(row[1]) * 18)
        delta_1tau_p_Ar.append(float(row[2]) * 18)

### vmu -> vmu mu+ mu- ; incoherent ; proton ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/proton/vmu_to_vmu_mu+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_p_Ar.append(float(row[0]))
        xsec_2mu_p_Ar.append(float(row[1]) * 18)
        delta_2mu_p_Ar.append(float(row[2]) * 18)

### vmu -> vmu mu+ mu- ; incoherent ; neutron ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/neutron/vmu_to_vmu_mu+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_n_Ar.append(float(row[0]))
        xsec_2mu_n_Ar.append(float(row[1]) * (40-18))
        delta_2mu_n_Ar.append(float(row[2]) * (40-18))

### vmu -> vmu tau+ tau- ; incoherent ; neutron ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/neutron/vmu_to_vmu_tau+_tau-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_n_Ar.append(float(row[0]))
        xsec_2tau_n_Ar.append(float(row[1]) * (40-18))
        delta_2tau_n_Ar.append(float(row[2]) * (40-18))

### vmu -> vtau tau+ mu- ; incoherent ; neutron ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/neutron/vmu_to_vtau_tau+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_n_Ar.append(float(row[0]))
        xsec_1tau_n_Ar.append(float(row[1]) * (40-18))
        delta_1tau_n_Ar.append(float(row[2]) * (40-18))

##### Coherent #####
### vmu -> vtau tau+ mu- ; coherent ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/coherent/argon/vmu_to_vtau_tau+_mu-_coh_Ar_xsec.csv','r') as csvfile: 
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) / energy * 10
        delta = float(row[2]) / energy * 10
        energy_1tau_coh_Ar.append(energy)
        xsec_1tau_coh_Ar.append(xsec)
        delta_1tau_coh_Ar.append(delta)

### vmu -> vmu tau+ tau- ; coherent ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/coherent/argon/vmu_to_vmu_tau+_tau-_coh_Ar_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) / energy * 10
        delta = float(row[2]) / energy * 10
        energy_2tau_coh_Ar.append(energy)
        xsec_2tau_coh_Ar.append(xsec)
        delta_2tau_coh_Ar.append(delta)

### vmu -> vmu mu+ mu- ; coherent ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/argon/vmu_to_vmu_mu+_mu-_coh_Ar_xsec_1em2_220.csv','r') as csvfile:  
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) / energy * 10
        delta = float(row[2]) / energy * 10
        energy_2mu_coh_Ar.append(energy)
        xsec_2mu_coh_Ar.append(xsec)
        delta_2mu_coh_Ar.append(delta)

##### Incoherent #####
### vmu -> vmu tau+ tau- ; incoherent ; proton ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) * 18 / energy * 10
        delta = float(row[2]) * 18 / energy * 10
        energy_2tau_p_Ar.append(energy)
        xsec_2tau_p_Ar.append(xsec)
        delta_2tau_p_Ar.append(delta)

### vmu -> vtau tau+ mu- ; incoherent ; proton ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) * 18 / energy * 10
        delta = float(row[2]) * 18 / energy * 10
        energy_1tau_p_Ar.append(energy)
        xsec_1tau_p_Ar.append(xsec)
        delta_1tau_p_Ar.append(delta)

### vmu -> vmu mu+ mu- ; incoherent ; proton ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/proton/vmu_to_vmu_mu+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) * 18 / energy * 10
        delta = float(row[2]) * 18 / energy * 10
        energy_2mu_p_Ar.append(energy)
        xsec_2mu_p_Ar.append(xsec)
        delta_2mu_p_Ar.append(delta)

### vmu -> vmu mu+ mu- ; incoherent ; neutron ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/neutron/vmu_to_vmu_mu+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) * (40-18) / energy * 10
        delta = float(row[2]) * (40-18) / energy * 10
        energy_2mu_n_Ar.append(energy)
        xsec_2mu_n_Ar.append(xsec)
        delta_2mu_n_Ar.append(delta)

### vmu -> vmu mu+ mu- ; incoherent ; total ; Argon ###
energy_2mu_incoh_Ar = energy_2mu_p_Ar
xsec_2mu_incoh_Ar   = [sum(i) for i in zip(xsec_2mu_n_Ar,xsec_2mu_p_Ar)]
delta_2mu_incoh_Ar  = [sum(i) for i in zip(delta_2mu_n_Ar, delta_2mu_p_Ar)]

##### Nucleon #####
### vmu -> vmu tau+ tau- ; nucleon ; proton ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) / energy * 10
        delta = float(row[2]) / energy * 10
        energy_2tau_p.append(energy)
        xsec_2tau_p.append(xsec)
        delta_2tau_p.append(delta)

### vmu -> vtau tau+ mu- ; nucleon ; proton ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) / energy * 10
        delta = float(row[2]) / energy * 10
        energy_1tau_p.append(energy)
        xsec_1tau_p.append(xsec)
        delta_1tau_p.append(delta)

### vmu -> vmu mu+ mu- ; nucleon ; proton ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/proton/vmu_to_vmu_mu+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) / energy * 10
        delta = float(row[2]) / energy * 10
        energy_2mu_p.append(energy)
        xsec_2mu_p.append(xsec)
        delta_2mu_p.append(delta)

### vmu -> vmu mu+ mu- ; nucleon ; neutron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/neutron/vmu_to_vmu_mu+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) / energy * 10
        delta = float(row[2]) / energy * 10
        energy_2mu_n.append(energy)
        xsec_2mu_n.append(xsec)
        delta_2mu_n.append(delta)

##### Digitized #####
### vmu -> vmu mu+ mu- ; coherent ; Argon ; Altmannshofer ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/argon/vmu_to_vmu_mu+_mu-_coh_Ar_xsec_Altmannshofer.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) / energy * 10
        energy_2mu_coh_Ar_Alt.append(energy)
        xsec_2mu_coh_Ar_Alt.append(xsec)

### vmu -> vmu mu+ mu- ; incoherent ; proton ; Argon ; Altmannshofer et al. ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/incoherent/proton/vmu_to_vmu_mu+_mu-_incoh_p_Ar_xsec_Altmannshofer.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) / energy * 10
        energy_2mu_incoh_p_Ar_Alt.append(energy)
        xsec_2mu_incoh_p_Ar_Alt.append(xsec)

### vmu -> vmu mu+ mu- ; incoherent ; neutron ; Argon ; Altmannshofer et al. ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/incoherent/neutron/vmu_to_vmu_mu+_mu-_incoh_n_Ar_xsec_Altmannshofer.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) / energy * 10
        energy_2mu_incoh_n_Ar_Alt.append(energy)
        xsec_2mu_incoh_n_Ar_Alt.append(xsec)

### vmu X -> mu- X' ; vmuCC ; Formaggio & Zeller
with open(CROSS_SECTION_DIR + '/vmuCC/vmuCC_xsec_perE_Formaggio.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_numuCC_FZ.append(float(row[0]))
        xsec_numuCC_FZ.append(float(row[1]))

####################
##### Plotting #####
####################

fig1, ax1 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)
fig2, ax2 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)

### vmu -> vtau tau+ mu- ###
ax1.errorbar(energy_1tau_coh_Ar, xsec_1tau_coh_Ar, yerr=delta_1tau_coh_Ar, color = 'firebrick', label = r"$\tau^+ \mu^-$ Ar")
ax1.errorbar(energy_1tau_p_Ar, xsec_1tau_p_Ar, yerr=delta_1tau_p_Ar, color = 'firebrick', linestyle = 'dashed', label = r'$\tau^+ \mu^- (p)$')
ax2.errorbar(energy_1tau_p, xsec_1tau_p, yerr=delta_1tau_p, color = 'firebrick', linestyle = 'dotted', label = r'$\tau^+ \mu^- (p)$')

### vmu -> vmu tau+ tau- ###
ax1.errorbar(energy_2tau_coh_Ar, xsec_2tau_coh_Ar, yerr=delta_2tau_coh_Ar, color = 'c', label = r"$\tau^+ \tau^-$ Ar")
ax1.errorbar(energy_2tau_p_Ar, xsec_2tau_p_Ar, yerr=delta_2tau_p_Ar, color = 'c', linestyle = 'dashed', label = r'$\tau^+ \tau^- (p)$')
ax2.errorbar(energy_2tau_p, xsec_2tau_p, yerr=delta_2tau_p, color = 'c', linestyle = 'dotted', label = r'$\tau^+ \tau^- (p)$')

### vmu -> vmu mu+ mu- ###
ax1.errorbar(energy_2mu_coh_Ar, xsec_2mu_coh_Ar, yerr=delta_2mu_coh_Ar, color = 'blueviolet', label = r'$\mu^+ \mu^-$ Ar')
ax1.errorbar(energy_2mu_p_Ar, xsec_2mu_p_Ar, yerr=delta_2mu_p_Ar, color = 'blueviolet', linestyle = 'dashed', label = r'$\mu^+ \mu^- (p)$')
ax1.errorbar(energy_2mu_n_Ar, xsec_2mu_n_Ar, yerr=delta_2mu_n_Ar, color = 'blueviolet', linestyle = 'dashdot', label = r'$\mu^+ \mu^- (n)$')
ax1.errorbar(energy_2mu_n_Ar, xsec_2mu_incoh_Ar, yerr=delta_2mu_incoh_Ar, color = 'blueviolet', linestyle = 'dotted', label = r'$\mu^+ \mu^- (p+n)$')
ax2.errorbar(energy_2mu_p, xsec_2mu_p, yerr=delta_2mu_p, color = 'blueviolet', linestyle = 'dotted', label = r'$\mu^+ \mu^- (p)$')
ax2.errorbar(energy_2mu_n, xsec_2mu_n, yerr=delta_2mu_n, color = 'blueviolet', linestyle = 'dashdot', label = r'$\mu^+ \mu^- (n)$')

### Digitized ###
ax1.scatter(energy_numuCC_FZ, xsec_numuCC_FZ, color = 'black', s=20, edgecolors='k', linewidths=1, alpha=0.75, label = r'$\nu_\mu$CC')
ax2.scatter(energy_numuCC_FZ, xsec_numuCC_FZ, color = 'black', s=20, edgecolors='k', linewidths=1, alpha=0.75, label = r'$\nu_\mu$CC')

ax1.set_xlabel('Energy (GeV)') 
ax1.set_ylabel(r'$\sigma$ / Energy $(10^{-38}\textrm{cm}^2 / \textrm{GeV})$') 
ax1.set_ylim(1e-9, 1.5)
ax1.set_xlim(0.1,600)
#ax1.legend(loc='lower right') 
ax1.legend(fontsize='25')
ax1.set_yscale('log')
ax1.set_xscale('log')

ax2.set_xlabel('Energy (GeV)') 
ax2.set_ylabel(r'$\sigma$ / Energy $(10^{-38}\textrm{cm}^2 / \textrm{GeV})$') 
ax2.set_ylim(1e-9, 1.5)
ax2.set_xlim(0.1,600)
ax2.legend(loc='lower right') 
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_title('Fermi gas model')

xmajor1 = [0.1, 1, 5, 10, 50, 100]
xmajor2 = [0.1, 1, 5, 10, 50, 100]
ax1.set_xticks(xmajor1, labels=['0.1', '1', '5', '10', '50', '100'])
ax1.grid(which='major', axis='both')
ax2.set_xticks(xmajor2, labels=['0.1', '1', '5', '10', '50', '100'])
ax2.grid(which='major', axis='both')

fig1.savefig("../plots/xsec_perE.png", dpi=400, bbox_inches='tight')
fig2.savefig("../plots/xsec_perE_nucleons.png", dpi=400, bbox_inches='tight') 
