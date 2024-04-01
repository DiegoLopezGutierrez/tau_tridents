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

energy_2tau_p_Ar = []
xsec_2tau_p_Ar = []
delta_2tau_p_Ar = []

energy_1tau_p_Ar = []
xsec_1tau_p_Ar = []
delta_1tau_p_Ar = []

energy_2mu_p_Ar = []
xsec_2mu_p_Ar = []
delta_2mu_p_Ar = []

energy_2mu_n_Ar = []
xsec_2mu_n_Ar = []
delta_2mu_n_Ar = []

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

##### Coherent #####
### vmu -> vtau tau+ mu- ; coherent ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/coherent/argon/vmu_to_vtau_tau+_mu-_coh_Ar_xsec.csv','r') as csvfile: 
    data = csv.reader(csvfile, delimiter = ',') 
    for row in data:
        energy_1tau_coh_Ar.append(float(row[0])) 
        xsec_1tau_coh_Ar.append(float(row[1]))
        delta_1tau_coh_Ar.append(float(row[2]))

### vmu -> vmu tau+ tau- ; coherent ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/coherent/argon/vmu_to_vmu_tau+_tau-_coh_Ar_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_coh_Ar.append(float(row[0]))
        xsec_2tau_coh_Ar.append(float(row[1]))
        delta_2tau_coh_Ar.append(float(row[2]))

### vmu -> vmu mu+ mu- ; coherent ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/argon/vmu_to_vmu_mu+_mu-_coh_Ar_xsec_1em2_220.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_coh_Ar.append(float(row[0]))
        xsec_2mu_coh_Ar.append(float(row[1]))
        delta_2mu_coh_Ar.append(float(row[2]))

##### Incoherent #####
### vmu -> vmu tau+ tau- ; incoherent ; proton ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_p_Ar.append(float(row[0]))
        xsec_2tau_p_Ar.append(float(row[1]) * 18)
        delta_2tau_p_Ar.append(float(row[2]) * 18)

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

### vmu -> vmu mu+ mu- ; incoherent ; total ; Argon ###
energy_2mu_incoh_Ar = energy_2mu_p_Ar
xsec_2mu_incoh_Ar   = [sum(i) for i in zip(xsec_2mu_n_Ar,xsec_2mu_p_Ar)]
delta_2mu_incoh_Ar  = [sum(i) for i in zip(delta_2mu_n_Ar, delta_2mu_p_Ar)]

##### Nucleon #####
### vmu -> vmu tau+ tau- ; nucleon ; proton ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_p.append(float(row[0]))
        xsec_2tau_p.append(float(row[1]))
        delta_2tau_p.append(float(row[2]))

### vmu -> vtau tau+ mu- ; nucleon ; proton ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_p.append(float(row[0]))
        xsec_1tau_p.append(float(row[1]))
        delta_1tau_p.append(float(row[2]))

### vmu -> vmu mu+ mu- ; nucleon ; proton ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/proton/vmu_to_vmu_mu+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_p.append(float(row[0]))
        xsec_2mu_p.append(float(row[1]))
        delta_2mu_p.append(float(row[2]))

### vmu -> vmu mu+ mu- ; nucleon ; neutron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/neutron/vmu_to_vmu_mu+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_n.append(float(row[0]))
        xsec_2mu_n.append(float(row[1]))
        delta_2mu_n.append(float(row[2]))

##### Digitized #####
### vmu -> vmu mu+ mu- ; coherent ; Argon ; Altmannshofer et al. ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/argon/vmu_to_vmu_mu+_mu-_coh_Ar_xsec_Altmannshofer.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_coh_Ar_Alt.append(float(row[0]))
        xsec_2mu_coh_Ar_Alt.append(float(row[1]))

### vmu -> vmu mu+ mu- ; incoherent ; proton ; Argon ; Altmannshofer et al. ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/incoherent/proton/vmu_to_vmu_mu+_mu-_incoh_p_Ar_xsec_Altmannshofer.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_incoh_p_Ar_Alt.append(float(row[0]))
        xsec_2mu_incoh_p_Ar_Alt.append(float(row[1]))

### vmu -> vmu mu+ mu- ; incoherent ; neutron ; Argon ; Altmannshofer et al. ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/incoherent/neutron/vmu_to_vmu_mu+_mu-_incoh_n_Ar_xsec_Altmannshofer.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_incoh_n_Ar_Alt.append(float(row[0]))
        xsec_2mu_incoh_n_Ar_Alt.append(float(row[1]))

### vmu X -> mu- X' ; vmuCC ; Formaggio & Zeller ###
with open(CROSS_SECTION_DIR + '/vmuCC/vmuCC_xsec_perE_Formaggio.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) * energy / 10 # digitized plot is in [1e-38 cm^2 / GeV]; 1 fb = 1e-39 cm^2
        energy_numuCC_FZ.append(energy)
        xsec_numuCC_FZ.append(xsec)

####################
#### Thresholds ####
####################

mmu = 0.105658
mtau = 1.778

Mproton = 0.938272
Mneutron = 0.939565

MArgon = 39.95*0.9315
MTungsten = 183.84*0.9315

def thresh_calc(lep1, lep2, target):
    return ((lep1 + lep2 + target)**2 - target**2)/(2*target)

##### Coherent ; Argon #####
thresh_2mu_coh_Ar = thresh_calc(mmu, mmu, MArgon)
thresh_1tau_coh_Ar = thresh_calc(mtau, mmu, MArgon)
thresh_2tau_coh_Ar = thresh_calc(mtau, mtau, MArgon)

##### Nucleon #####
thresh_2mu_p = thresh_calc(mmu, mmu, Mproton)
thresh_2mu_n = thresh_calc(mmu, mmu, Mneutron)
thresh_1tau_p = thresh_calc(mtau, mmu, Mproton)
thresh_2tau_p = thresh_calc(mtau, mtau, Mproton)

####################
##### Plotting #####
####################

fig1, ax1 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Validation cross sections
fig2, ax2 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Tau cross sections + vmuCC
fig3, ax3 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Nucleon cross sections

### vmu -> vtau tau+ mu- ###
# Plots #
ax2.errorbar(energy_1tau_coh_Ar, xsec_1tau_coh_Ar, yerr=delta_1tau_coh_Ar, color = 'firebrick', label = r"$\tau^+ \mu^-$")
ax2.errorbar(energy_1tau_p_Ar, xsec_1tau_p_Ar, yerr=delta_1tau_p_Ar, color = 'firebrick', linestyle = 'dashed', label = r'$\tau^+ \mu^- (p)$')
ax3.errorbar(energy_1tau_p, xsec_1tau_p, yerr=delta_1tau_p, color = 'firebrick', linestyle = 'dotted', label = r'$\tau^+ \mu^- (p)$')

# Thresholds #
ax2.axvline(thresh_1tau_coh_Ar,color='gray',linestyle='--')
ax2.text(thresh_1tau_coh_Ar*1.1, 2.5e-2, r'$1\tau$ Ar $=$ {:.1f} GeV'.format(thresh_1tau_coh_Ar), color='gray')
ax3.axvline(thresh_1tau_p,color='gray',linestyle='--')
ax3.text(thresh_1tau_p*1.1, 0.005, r'$1\tau (p) = ${:.1f} GeV'.format(thresh_1tau_p), color='gray')

### vmu -> vmu tau+ tau- ###
# Plots #
ax2.errorbar(energy_2tau_coh_Ar, xsec_2tau_coh_Ar, yerr=delta_2tau_coh_Ar, color = 'c', label = r"$\tau^+ \tau^-$")
ax2.errorbar(energy_2tau_p_Ar, xsec_2tau_p_Ar, yerr=delta_2tau_p_Ar, color = 'c', linestyle = 'dashed', label = r'$\tau^+ \tau^- (p)$')
ax3.errorbar(energy_2tau_p, xsec_2tau_p, yerr=delta_2tau_p, color = 'c', linestyle = 'dotted', label = r'$\tau^+ \tau^- (p)$')

# Thresholds #
ax2.axvline(thresh_2tau_coh_Ar,color='gray',linestyle='--')
ax2.text(thresh_2tau_coh_Ar*1.1, 1e-3, r'$2\tau$ Ar $=$ {:.1f} GeV'.format(thresh_2tau_coh_Ar), color='gray')
ax3.axvline(thresh_2tau_p,color='gray',linestyle='--')
ax3.text(thresh_2tau_p*1.1, 0.5, r'$2\tau (p) = ${:.1f} GeV'.format(thresh_2tau_p), color='gray')

### vmu -> vmu mu+ mu- ###
# Plots #
ax1.errorbar(energy_2mu_coh_Ar, xsec_2mu_coh_Ar, yerr=delta_2mu_coh_Ar, color = 'blueviolet', label = r'$\mu^+ \mu^-$')
ax1.errorbar(energy_2mu_p_Ar, xsec_2mu_p_Ar, yerr=delta_2mu_p_Ar, color = 'c', linestyle = 'dashed', label = r'$\mu^+ \mu^- (p)$')
ax1.errorbar(energy_2mu_n_Ar, xsec_2mu_n_Ar, yerr=delta_2mu_n_Ar, color = 'c', linestyle = 'dashdot', label = r'$\mu^+ \mu^- (n)$')
ax1.errorbar(energy_2mu_n_Ar, xsec_2mu_incoh_Ar, yerr=delta_2mu_incoh_Ar, color = 'c', linestyle = 'solid', label = r'$\mu^+ \mu^- (p+n)$')
ax3.errorbar(energy_2mu_p, xsec_2mu_p, yerr=delta_2mu_p, color = 'blueviolet', linestyle = 'dotted', label = r'$\mu^+ \mu^- (p)$')
ax3.errorbar(energy_2mu_n, xsec_2mu_n, yerr=delta_2mu_n, color = 'blueviolet', linestyle = 'dashdot', label = r'$\mu^+ \mu^- (n)$')

# Thresholds #
ax1.axvline(thresh_2mu_coh_Ar,color='gray',linestyle='--')
ax1.text(thresh_2mu_coh_Ar*1.1, 1, r'$2\mu$ Ar $=$ {:.1f} GeV'.format(thresh_2mu_coh_Ar), color='gray')
ax3.axvline(thresh_2mu_p,color='gray',linestyle='--')
ax3.text(thresh_2mu_p*1.1, 0.5, r'$2\mu (p) = ${:.1f} GeV'.format(thresh_2mu_p), color='gray')
ax3.axvline(thresh_2mu_n,color='gray',linestyle='--')
ax3.text(thresh_2mu_n*1.1, 0.005, r'$2\mu (n) =${:.1f} GeV'.format(thresh_2mu_n), color='gray')

### Digitized ###
ax1.scatter(energy_2mu_coh_Ar_Alt, xsec_2mu_coh_Ar_Alt, marker='o', color = 'black', s=20, edgecolors='k', linewidths=1, alpha=0.75) 
ax1.scatter(energy_2mu_incoh_p_Ar_Alt, xsec_2mu_incoh_p_Ar_Alt, color = 'black', s=20, marker='^', edgecolors='k', linewidths=1, alpha=0.75)
ax1.scatter(energy_2mu_incoh_n_Ar_Alt, xsec_2mu_incoh_n_Ar_Alt, color = 'black', s=20, marker='v', edgecolors='k', linewidths=1, alpha=0.75)

ax2.scatter(energy_numuCC_FZ, xsec_numuCC_FZ, color = 'black', s=20, edgecolors='k', linewidths=1, alpha=0.75, label = r'$\nu_\mu$CC')

### Plotting options ###
ax1.set_xlabel('Energy (GeV)') 
ax1.set_ylabel('Cross Section (fb)') 
ax1.set_ylim(1e-14, 1e2)
ax1.set_xlim(0.1,500)
ax1.legend(loc='lower right') 
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_title(r'$^{40}$Ar')

ax2.set_xlabel('Energy (GeV)') 
ax2.set_ylabel('Cross Section (fb)') 
ax2.set_ylim(1e-24, 1e2)
ax2.set_xlim(0.1,500)
ax2.legend(loc='lower right') 
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_title(r'$^{40}$Ar')

ax3.set_xlabel('Energy (GeV)') 
ax3.set_ylabel('Cross Section (fb)') 
ax3.set_ylim(1e-24, 1e1)
ax3.legend(loc='lower right') 
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.set_title('Fermi gas model')

xmajor = [0.1, 1, 5, 10, 50, 100]
ax1.set_xticks(xmajor, labels=['0.1', '1', '5', '10', '50', '100'])
ax1.grid(which='major', axis='both')
ax2.set_xticks(xmajor, labels=['0.1', '1', '5', '10', '50', '100'])
ax2.grid(which='major', axis='both')
ax3.set_xticks(xmajor, labels=['0.1', '1', '5', '10', '50', '100'])
ax3.grid(which='major', axis='both')

fig1.savefig("../plots/xsec_validation.png", dpi=400, bbox_inches='tight')
fig2.savefig("../plots/xsec_tau.png", dpi=400, bbox_inches='tight')
fig3.savefig("../plots/xsec_nucleons.png", dpi=400, bbox_inches='tight') 
