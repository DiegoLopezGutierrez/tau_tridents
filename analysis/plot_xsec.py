import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe
import csv 

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

CROSS_SECTION_DIR = '../csv/cross_sections'

Z_Ar = 18
A_Ar = 40

Z_W = 74
A_W = 184

###########################
###### Coherent xsec ######
###########################

### Argon ###
energy_1tau_coh_Ar = [] 
xsec_1tau_coh_Ar = []
delta_1tau_coh_Ar = []

energy_2tau_coh_Ar = []
xsec_2tau_coh_Ar = []
delta_2tau_coh_Ar = []

energy_2mu_coh_Ar = []
xsec_2mu_coh_Ar = []
delta_2mu_coh_Ar = []

### Tungsten ###
energy_1tau_coh_W = [] 
xsec_1tau_coh_W = []
delta_1tau_coh_W = []

energy_2tau_coh_W = []
xsec_2tau_coh_W = []
delta_2tau_coh_W = []

############################
##### Incoherent xsec ######
############################

##### Argon #####
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

##### Tungsten #####
### Proton ###
energy_2tau_p_W = []
xsec_2tau_p_W = []
delta_2tau_p_W = []

energy_1tau_p_W = []
xsec_1tau_p_W = []
delta_1tau_p_W = []

energy_2mu_p_W = []
xsec_2mu_p_W = []
delta_2mu_p_W = []

### Neutron ###
energy_2tau_n_W = []
xsec_2tau_n_W = []
delta_2tau_n_W = []

energy_1tau_n_W = []
xsec_1tau_n_W = []
delta_1tau_n_W = []

energy_2mu_n_W = []
xsec_2mu_n_W = []
delta_2mu_n_W = []

### Proton + Neutron ###
energy_2tau_incoh_W = []
xsec_2tau_incoh_W = []
delta_2tau_incoh_W = []

energy_1tau_incoh_W = []
xsec_1tau_incoh_W = []
delta_1tau_incoh_W = []

energy_2mu_incoh_W = []
xsec_2mu_incoh_W = []
delta_2mu_incoh_W = []

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

### vmu -> vtau tau+ mu- ; coherent ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/coherent/tungsten/vmu_to_vtau_tau+_mu-_coh_W_xsec.csv','r') as csvfile: 
    data = csv.reader(csvfile, delimiter = ',') 
    for row in data:
        energy_1tau_coh_W.append(float(row[0]))
        xsec_1tau_coh_W.append(float(row[1]))
        delta_1tau_coh_W.append(float(row[2]))

### vmu -> vmu tau+ tau- ; coherent ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/coherent/tungsten/vmu_to_vmu_tau+_tau-_coh_W_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_coh_W.append(float(row[0]))
        xsec_2tau_coh_W.append(float(row[1]))
        delta_2tau_coh_W.append(float(row[2]))

##### Incoherent #####
### vmu -> vmu tau+ tau- ; incoherent ; proton ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_p_Ar.append(float(row[0]))
        xsec_2tau_p_Ar.append(float(row[1]) * Z_Ar)
        delta_2tau_p_Ar.append(float(row[2]) * Z_Ar)

### vmu -> vtau tau+ mu- ; incoherent ; proton ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_p_Ar.append(float(row[0]))
        xsec_1tau_p_Ar.append(float(row[1]) * Z_Ar)
        delta_1tau_p_Ar.append(float(row[2]) * Z_Ar)

### vmu -> vmu mu+ mu- ; incoherent ; proton ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/proton/vmu_to_vmu_mu+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_p_Ar.append(float(row[0]))
        xsec_2mu_p_Ar.append(float(row[1]) * Z_Ar)
        delta_2mu_p_Ar.append(float(row[2]) * Z_Ar)

### vmu -> vmu mu+ mu- ; incoherent ; neutron ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/neutron/vmu_to_vmu_mu+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu_n_Ar.append(float(row[0]))
        xsec_2mu_n_Ar.append(float(row[1]) * (A_Ar-Z_Ar))
        delta_2mu_n_Ar.append(float(row[2]) * (A_Ar-Z_Ar))

### vmu -> vmu tau+ tau- ; incoherent ; neutron ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/new_neutron/vmu_to_vmu_tau+_tau-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_n_Ar.append(float(row[0]))
        xsec_2tau_n_Ar.append(float(row[1]) * (A_Ar-Z_Ar))
        delta_2tau_n_Ar.append(float(row[2]) * (A_Ar-Z_Ar))

### vmu -> vtau tau+ mu- ; incoherent ; neutron ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/new_neutron/vmu_to_vtau_tau+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_n_Ar.append(float(row[0]))
        xsec_1tau_n_Ar.append(float(row[1]) * (A_Ar-Z_Ar))
        delta_1tau_n_Ar.append(float(row[2]) * (A_Ar-Z_Ar))

### vmu -> vmu tau+ tau- ; incoherent ; proton ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_p_W.append(float(row[0]))
        xsec_2tau_p_W.append(float(row[1]) * Z_W)
        delta_2tau_p_W.append(float(row[2]) * Z_W)

### vmu -> vtau tau+ mu- ; incoherent ; proton ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_p_W.append(float(row[0]))
        xsec_1tau_p_W.append(float(row[1]) * Z_W)
        delta_1tau_p_W.append(float(row[2]) * Z_W)

### vmu -> vmu tau+ tau- ; incoherent ; neutron ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/new_neutron/vmu_to_vmu_tau+_tau-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_n_W.append(float(row[0]))
        xsec_2tau_n_W.append(float(row[1]) * (A_W-Z_W))
        delta_2tau_n_W.append(float(row[2]) * (A_W-Z_W))

### vmu -> vtau tau+ mu- ; incoherent ; neutron ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/new_neutron/vmu_to_vtau_tau+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_n_W.append(float(row[0]))
        xsec_1tau_n_W.append(float(row[1]) * (A_W-Z_W))
        delta_1tau_n_W.append(float(row[2]) * (A_W-Z_W))

### vmu -> vmu mu+ mu- ; incoherent ; total ; Argon ###
energy_2mu_incoh_Ar = energy_2mu_p_Ar
xsec_2mu_incoh_Ar   = [sum(i) for i in zip(xsec_2mu_n_Ar,xsec_2mu_p_Ar)]
delta_2mu_incoh_Ar  = [sum(i) for i in zip(delta_2mu_n_Ar, delta_2mu_p_Ar)]
print("Size of incoh 2mu: ", len(energy_2mu_incoh_Ar), len(xsec_2mu_incoh_Ar))

### vmu -> vtau tau+ mu- ; incoherent ; total ; Argon ###
energy_1tau_incoh_Ar = energy_1tau_n_Ar
xsec_1tau_incoh_Ar   = [sum(i) for i in zip(xsec_1tau_n_Ar,xsec_1tau_p_Ar)]
delta_1tau_incoh_Ar  = [sum(i) for i in zip(delta_1tau_n_Ar, delta_1tau_p_Ar)]
print("Size of incoh 1tau: ", len(energy_1tau_incoh_Ar), len(xsec_1tau_incoh_Ar))
print(energy_1tau_incoh_Ar[-1])

### vmu -> vmu tau+ tau- ; incoherent ; total ; Argon ###
energy_2tau_incoh_Ar = energy_2tau_n_Ar
xsec_2tau_incoh_Ar   = [sum(i) for i in zip(xsec_2tau_n_Ar,xsec_2tau_p_Ar)]
delta_2tau_incoh_Ar  = [sum(i) for i in zip(delta_2tau_n_Ar, delta_2tau_p_Ar)]
print("Size of incoh 2tau: ", len(energy_2tau_incoh_Ar), len(xsec_2tau_incoh_Ar))
print(energy_2tau_incoh_Ar[-1])

### vmu -> vtau tau+ mu- ; incoherent ; total ; Tungsten ###
energy_1tau_incoh_W = energy_1tau_n_W
xsec_1tau_incoh_W   = [sum(i) for i in zip(xsec_1tau_n_W,xsec_1tau_p_W)]
delta_1tau_incoh_W  = [sum(i) for i in zip(delta_1tau_n_W, delta_1tau_p_W)]
print("Size of incoh 1tau: ", len(energy_1tau_incoh_W), len(xsec_1tau_incoh_W))
print(energy_1tau_incoh_W[-1])

### vmu -> vmu tau+ tau- ; incoherent ; total ; Tungsten ###
energy_2tau_incoh_W = energy_2tau_n_W
xsec_2tau_incoh_W   = [sum(i) for i in zip(xsec_2tau_n_W,xsec_2tau_p_W)]
delta_2tau_incoh_W  = [sum(i) for i in zip(delta_2tau_n_W, delta_2tau_p_W)]
print("Size of incoh 2tau: ", len(energy_2tau_incoh_W), len(xsec_2tau_incoh_W))
print(energy_2tau_incoh_W[-1])

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

### vmu -> vmu tau+ tau- ; nucleon ; neutron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/new_neutron/vmu_to_vmu_tau+_tau-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau_n.append(float(row[0]))
        xsec_2tau_n.append(float(row[1]))
        delta_2tau_n.append(float(row[2]))

### vmu -> vtau tau+ mu- ; nucleon ; neutron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/new_neutron/vmu_to_vtau_tau+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau_n.append(float(row[0]))
        xsec_1tau_n.append(float(row[1]))
        delta_1tau_n.append(float(row[2]))

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
        xsec = float(row[1]) * energy / 10 * A_Ar # digitized plot is in [1e-38 cm^2 / GeV]; 1 fb = 1e-39 cm^2; also, xsec given in as xsec/nucleon -> xsec/Ar
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

###### Coherent ; Argon #####
#thresh_2mu_coh_Ar = thresh_calc(mmu, mmu, MArgon)
#thresh_1tau_coh_Ar = thresh_calc(mtau, mmu, MArgon)
#thresh_2tau_coh_Ar = thresh_calc(mtau, mtau, MArgon)
#
###### Nucleon #####
#thresh_2mu_p = thresh_calc(mmu, mmu, Mproton)
#thresh_2mu_n = thresh_calc(mmu, mmu, Mneutron)
#thresh_1tau_p = thresh_calc(mtau, mmu, Mproton)
#thresh_2tau_p = thresh_calc(mtau, mtau, Mproton)

####################
##### Plotting #####
####################

fig1, ax1 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Validation cross sections
fig2, ax2 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Tau cross sections + vmuCC
fig3, ax3 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Nucleon cross sections
fig4, ax4 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Tungsten tau cross sections

### vmu -> vtau tau+ mu- ###
# Argon #
ax2.errorbar(energy_1tau_coh_Ar, xsec_1tau_coh_Ar, yerr=delta_1tau_coh_Ar, color = 'firebrick', label = r"$\tau^+ \mu^-$ Coh.", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.errorbar(energy_1tau_incoh_Ar, xsec_1tau_incoh_Ar, yerr=delta_1tau_incoh_Ar, color = 'firebrick', linestyle = 'dashed', label = r'$\tau^+ \mu^-$ Incoh.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.text(50,3e-4,r'$\tau^+ \mu^-$ {\bf Coh.}',color='firebrick',rotation=20,fontsize=30)
ax2.text(13,2e-4,r'$\tau^+ \mu^-$ {\bf Incoh.}',color='firebrick',rotation=20,fontsize=30)

# Tungsten #
ax4.errorbar(energy_1tau_coh_W, xsec_1tau_coh_W, yerr=delta_1tau_coh_W, color = 'firebrick', label = r"$\tau^+ \mu^-$ Coh.", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax4.errorbar(energy_1tau_incoh_W, xsec_1tau_incoh_W, yerr=delta_1tau_incoh_W, color = 'firebrick', linestyle = 'dashed', label = r'$\tau^+ \mu^-$ Incoh.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax4.text(10,5e-7,r'$\tau^+ \mu^-$ {\bf Coh.}',color='firebrick',rotation=37,fontsize=30) 
ax4.text(10,8e-4,r'$\tau^+ \mu^-$ {\bf Incoh.}',color='firebrick',rotation=20,fontsize=30) 

# Nucleon #
ax3.errorbar(energy_1tau_p, xsec_1tau_p, yerr=delta_1tau_p, color = 'firebrick', linestyle = 'solid', label = r'$\tau^+ \mu^- (p)$', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.errorbar(energy_1tau_n, xsec_1tau_n, yerr=delta_1tau_n, color = 'firebrick', linestyle = 'dashed', label = r'$\tau^+ \mu^- (n)$', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.text(5,1e-7,r'$\tau^+ \mu^-$ {\bf p}',color='firebrick',rotation=65,fontsize=30)
ax3.text(300,2e-3,r'$\tau^+ \mu^-$ {\bf n}',color='firebrick',rotation=24,fontsize=30)

### vmu -> vmu tau+ tau- ###
# Argon #
ax2.errorbar(energy_2tau_coh_Ar, xsec_2tau_coh_Ar, yerr=delta_2tau_coh_Ar, color = 'c', label = r"$\tau^+ \tau^-$ Coh.", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.errorbar(energy_2tau_incoh_Ar, xsec_2tau_incoh_Ar, yerr=delta_2tau_incoh_Ar, color = 'c', linestyle = 'dashed', label = r'$\tau^+ \tau^-$ Incoh.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.text(50,2e-8,r'$\tau^+ \tau^-$ {\bf Coh.}',color='c',rotation=47,fontsize=30)
ax2.text(25,9e-7,r'$\tau^+ \tau^-$ {\bf Incoh.}',color='c',rotation=30,fontsize=30)

# Tungsten #
ax4.errorbar(energy_2tau_coh_W, xsec_2tau_coh_W, yerr=delta_2tau_coh_W, color = 'c', label = r"$\tau^+ \tau^-$ Coh.", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax4.errorbar(energy_2tau_incoh_W, xsec_2tau_incoh_W, yerr=delta_2tau_incoh_W, color = 'c', linestyle = 'dashed', label = r'$\tau^+ \tau^-$ Incoh.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])

ax4.text(50,1e-7,r'$\tau^+ \tau^-$ {\bf Coh.}',color='c',rotation=34,fontsize=30)
ax4.text(50,1e-4,r'$\tau^+ \tau^-$ {\bf Incoh.}',color='c',rotation=18,fontsize=30)

# Nucleon #
ax3.errorbar(energy_2tau_p, xsec_2tau_p, yerr=delta_2tau_p, color = 'c', linestyle = 'solid', label = r'$\tau^+ \tau^- (p)$', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.errorbar(energy_2tau_n, xsec_2tau_n, yerr=delta_2tau_n, color = 'c', linestyle = 'dashed', label = r'$\tau^+ \tau^- (n)$', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.text(17,1e-7,r'$\tau^+ \tau^-$ {\bf p}',color='c',rotation=62,fontsize=30)
ax3.text(200,1e-5,r'$\tau^+ \tau^-$ {\bf n}',color='c',rotation=34,fontsize=30)

### vmu -> vmu mu+ mu- ###
# Argon #
ax1.errorbar(energy_2mu_coh_Ar, xsec_2mu_coh_Ar, yerr=delta_2mu_coh_Ar, color = 'blueviolet', label = r'$\mu^+ \mu^-$ Coh.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax1.errorbar(energy_2mu_n_Ar, xsec_2mu_incoh_Ar, yerr=delta_2mu_incoh_Ar, color = 'blueviolet', linestyle = 'dashed', label = r'$\mu^+ \mu^-$ Incoh.', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])


# Nucleon #
ax3.errorbar(energy_2mu_p, xsec_2mu_p, yerr=delta_2mu_p, color = 'blueviolet', linestyle = 'solid', label = r'$\mu^+ \mu^- (p)$', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.errorbar(energy_2mu_n, xsec_2mu_n, yerr=delta_2mu_n, color = 'blueviolet', linestyle = 'dashed', label = r'$\mu^+ \mu^- (n)$', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.text(0.8,1e-6,r'$\mu^+ \mu^-$ {\bf p}',color='blueviolet',rotation=40,fontsize=30)
ax3.text(1.7,1e-6,r'$\mu^+ \mu^-$ {\bf n}',color='blueviolet',rotation=30,fontsize=30)

### Digitized ###
ax1.scatter(energy_2mu_coh_Ar_Alt, xsec_2mu_coh_Ar_Alt, marker='o', color = 'black', s=20, edgecolors='k', linewidths=1, alpha=0.75) 
ax1.scatter(energy_2mu_incoh_p_Ar_Alt, xsec_2mu_incoh_p_Ar_Alt, color = 'black', s=20, marker='^', edgecolors='k', linewidths=1, alpha=0.75)
ax1.scatter(energy_2mu_incoh_n_Ar_Alt, xsec_2mu_incoh_n_Ar_Alt, color = 'black', s=20, marker='v', edgecolors='k', linewidths=1, alpha=0.75)

ax2.scatter(energy_numuCC_FZ, xsec_numuCC_FZ, color = 'black', s=20, edgecolors='k', linewidths=1, alpha=0.75, label = r'$\nu_\mu$CC')
ax2.text(10,1,r'$\nu_\mu$ {\bf CC}',color='black',rotation=10,fontsize=30)

### Plotting options ###
ax1.set_xlabel(r'{\bf Neutrino Energy} (GeV)') 
ax1.set_ylabel(r'{\bf Cross Section / Ar} (fb)') 
ax1.set_ylim(1e-10, 1e2)
ax1.set_xlim(0.1,500)
ax1.legend(loc='lower right') 
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_title(r'{\bf Argon $^{40}$Ar}', fontsize=40)

ax2.set_xlabel(r'{\bf Neutrino Energy} (GeV)') 
ax2.set_ylabel(r'{\bf Cross Section / Ar} (fb)') 
ax2.set_ylim(1e-15, 1e2)
ax2.set_xlim(2,1100)
#ax2.legend(loc='lower right') 
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_title(r'{\bf Argon $^{40}$Ar}', fontsize=40)

ax3.set_xlabel(r'{\bf Neutrino Energy} (GeV)') 
ax3.set_ylabel(r'{\bf Cross Section / Nucleon} (fb)') 
ax3.set_ylim(1e-12, 1e0)
ax3.set_xlim(0.1,1100)
#ax3.legend(loc='lower right') 
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.set_title(r'{\bf Fermi gas model}', fontsize=40)

ax4.set_xlabel(r'{\bf Neutrino Energy} (GeV)') 
ax4.set_ylabel(r'{\bf Cross Section / W} (fb)') 
ax4.set_ylim(1e-15, 1e2)
ax4.set_xlim(2,1100)
#ax4.legend(loc='lower right') 
ax4.set_yscale('log')
ax4.set_xscale('log')
ax4.set_title(r'{\bf Tungsten $^{184}$W}', fontsize=40)

xmajor1 = [0.1, 1, 5, 10, 50, 100]
xmajor2 = [5, 10, 50, 100, 1000]
xmajor3 = [0.1, 1, 5, 10, 50, 100, 1000]
xmajor4 = [5, 10, 50, 100, 1000]

ax1.set_xticks(xmajor1, labels=['0.1', '1', '5', '10', '50', '100'])
ax1.grid(which='major', axis='both')
ax2.set_xticks(xmajor2, labels=['5', '10', '50', '100', '1000'])
ax2.grid(which='major', axis='both')
ax3.set_xticks(xmajor3, labels=['0.1', '1', '5', '10', '50', '100', '1000'])
ax3.grid(which='major', axis='both')
ax4.set_xticks(xmajor4, labels=['5', '10', '50', '100', '1000'])
ax4.grid(which='major', axis='both')

fig1.savefig("../plots/xsec_mu_argon.png", dpi=400, bbox_inches='tight')
fig2.savefig("../plots/xsec_tau_argon.png", dpi=400, bbox_inches='tight')
fig3.savefig("../plots/xsec_nucleons.png", dpi=400, bbox_inches='tight')
fig4.savefig("../plots/xsec_tau_tungsten.png", dpi=400, bbox_inches='tight')
