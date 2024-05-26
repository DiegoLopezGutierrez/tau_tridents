import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe
import csv 
import numpy as np

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

CROSS_SECTION_DIR = '../csv/cross_sections'

###########################
###### Constants ######
###########################

## Nuclear Parameters ##
Z_Ar = 18
A_Ar = 40

Z_W = 74
A_W = 184

## Flux Ranges ##
DUNE_std_low = 0.125
DUNE_std_high = 60.0

DUNE_tau_opt_low = 1

## Cross Section Uncertainties ##
# Components # (see Altmannshofer et al. for details)
aEM = 1/137  # low q^2 usually, so using zero momentum value of fine structure
sigma_highQED_Ar = Z_Ar*aEM/(4*np.pi) + 0.02  # higher order QED corrections should give roughly Z*aEM/(4*pi) which is about 1% for Ar. Add 2% to be conservative.
sigma_highQED_W  = Z_W*aEM/(4*np.pi) + 0.02 # higher order QED corrections should give roughly Z*aEM/(4*pi) which is about 4% for W. Add 2% to be conservative.
sigma_highQED_p = 1*aEM/(4*np.pi) + 0.02
sigma_highQED_n = 0

sigma_form_factors_coh = 0.01 # 1% uncertainty on coherent cross sections due to form factors
sigma_form_factors_incoh = 0.03 # 3% uncertainty on incoherent cross sections due to form factors

sigma_highEW = 0.05 # 5% uncertainty due to higher order weak corrections

sigma_nuclear_modeling = 0.30 # for incoherent scattering, largest uncertainty due to other nuclear effects besides the included Pauli blocking. 30% to be conservative.

# Total #
sigma_total_coh_Ar = np.sqrt(sigma_highQED_Ar**2 + sigma_form_factors_coh**2 + sigma_highEW**2)
sigma_total_coh_W  = np.sqrt(sigma_highQED_W**2 + sigma_form_factors_coh**2 + sigma_highEW**2)
print(sigma_total_coh_Ar, sigma_total_coh_W)

sigma_total_incoh_Ar = np.sqrt(sigma_highQED_Ar**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)
sigma_total_incoh_W  = np.sqrt(sigma_highQED_W**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)
print(sigma_total_incoh_Ar, sigma_total_incoh_W)

sigma_total_p = np.sqrt(sigma_highQED_p**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)
sigma_total_n = np.sqrt(sigma_highQED_n**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)
print(sigma_total_p, sigma_total_n)

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

energy_2mu_coh_W = []
xsec_2mu_coh_W = []
delta_2mu_coh_W = []

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

energy_2mu_coh_Ar_Ballett = []
xsec_2mu_coh_Ar_Ballett = []

energy_2tau_coh_Ar_Beacom = []
xsec_2tau_coh_Ar_Beacom = []

energy_2tau_incoh_Ar_Beacom = []
xsec_2tau_incoh_Ar_Beacom = []

energy_2mu_coh_Ar_Magill = []
xsec_2mu_coh_Ar_Magill = []

###############################
##### Load CSV xsec files #####
###############################

##### Conversion Factor #####
fb_to_cm2 = 1e-39 # 1 fb = 1e-39 cm^2

##### Coherent #####
### vmu -> vtau tau+ mu- ; coherent ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/coherent/argon/vmu_to_vtau_tau+_mu-_coh_Ar_xsec.csv','r') as csvfile: 
    data = csv.reader(csvfile, delimiter = ',') 
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_1tau_coh_Ar.append(energy) 
        xsec_1tau_coh_Ar.append(xsec_perE)
        delta_1tau_coh_Ar.append(delta)

### vmu -> vmu tau+ tau- ; coherent ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/coherent/argon/vmu_to_vmu_tau+_tau-_coh_Ar_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_2tau_coh_Ar.append(energy)
        xsec_2tau_coh_Ar.append(xsec_perE)
        delta_2tau_coh_Ar.append(delta)

### vmu -> vmu mu+ mu- ; coherent ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/argon/vmu_to_vmu_mu+_mu-_coh_Ar_xsec_1em2_220.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_2mu_coh_Ar.append(energy)
        xsec_2mu_coh_Ar.append(xsec_perE)
        delta_2mu_coh_Ar.append(delta)

### vmu -> vtau tau+ mu- ; coherent ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/coherent/tungsten/vmu_to_vtau_tau+_mu-_coh_W_xsec.csv','r') as csvfile: 
    data = csv.reader(csvfile, delimiter = ',') 
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_W * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_1tau_coh_W.append(energy)
        xsec_1tau_coh_W.append(xsec_perE)
        delta_1tau_coh_W.append(delta)

### vmu -> vmu tau+ tau- ; coherent ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/coherent/tungsten/vmu_to_vmu_tau+_tau-_coh_W_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_W * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_2tau_coh_W.append(energy)
        xsec_2tau_coh_W.append(xsec_perE)
        delta_2tau_coh_W.append(delta)

### vmu -> vmu mu+ mu- ; coherent ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/tungsten/vmu_to_vmu_mu+_mu-_coh_W_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_W * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_2mu_coh_W.append(energy)
        xsec_2mu_coh_W.append(xsec_perE)
        delta_2mu_coh_W.append(delta)

##### Incoherent #####
### vmu -> vmu tau+ tau- ; incoherent ; proton ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/new_proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * Z_Ar * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * Z_Ar * fb_to_cm2

        energy_2tau_p_Ar.append(energy)
        xsec_2tau_p_Ar.append(xsec_perE)
        delta_2tau_p_Ar.append(delta)

### vmu -> vtau tau+ mu- ; incoherent ; proton ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/new_proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * Z_Ar * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * Z_Ar * fb_to_cm2

        energy_1tau_p_Ar.append(energy)
        xsec_1tau_p_Ar.append(xsec_perE)
        delta_1tau_p_Ar.append(delta)

### vmu -> vmu mu+ mu- ; incoherent ; proton ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/new_proton/vmu_to_vmu_mu+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * Z_Ar * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * Z_Ar * fb_to_cm2

        energy_2mu_p_Ar.append(energy)
        xsec_2mu_p_Ar.append(xsec_perE)
        delta_2mu_p_Ar.append(delta)

### vmu -> vmu mu+ mu- ; incoherent ; neutron ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/new_neutron/vmu_to_vmu_mu+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * (A_Ar-Z_Ar) * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * (A_Ar-Z_Ar) * fb_to_cm2

        energy_2mu_n_Ar.append(energy)
        xsec_2mu_n_Ar.append(xsec_perE)
        delta_2mu_n_Ar.append(delta)

### vmu -> vmu tau+ tau- ; incoherent ; neutron ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/new_neutron/vmu_to_vmu_tau+_tau-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * (A_Ar-Z_Ar) * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * (A_Ar-Z_Ar) * fb_to_cm2

        energy_2tau_n_Ar.append(energy)
        xsec_2tau_n_Ar.append(xsec_perE)
        delta_2tau_n_Ar.append(delta)

### vmu -> vtau tau+ mu- ; incoherent ; neutron ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/new_neutron/vmu_to_vtau_tau+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * (A_Ar-Z_Ar) * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * (A_Ar-Z_Ar) * fb_to_cm2

        energy_1tau_n_Ar.append(energy)
        xsec_1tau_n_Ar.append(xsec_perE)
        delta_1tau_n_Ar.append(delta)

### vmu -> vmu tau+ tau- ; incoherent ; proton ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/new_proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * Z_W * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_W * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * Z_W * fb_to_cm2

        energy_2tau_p_W.append(energy)
        xsec_2tau_p_W.append(xsec_perE)
        delta_2tau_p_W.append(delta)

### vmu -> vtau tau+ mu- ; incoherent ; proton ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/new_proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * Z_W * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_W * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * Z_W * fb_to_cm2

        energy_1tau_p_W.append(energy)
        xsec_1tau_p_W.append(xsec_perE)
        delta_1tau_p_W.append(delta)

### vmu -> vmu tau+ tau- ; incoherent ; neutron ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/new_neutron/vmu_to_vmu_tau+_tau-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * (A_W-Z_W) * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_W * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * (A_W-Z_W) * fb_to_cm2

        energy_2tau_n_W.append(energy)
        xsec_2tau_n_W.append(xsec_perE)
        delta_2tau_n_W.append(delta)

### vmu -> vtau tau+ mu- ; incoherent ; neutron ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/new_neutron/vmu_to_vtau_tau+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * (A_W-Z_W) * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_W * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * (A_W-Z_W) * fb_to_cm2

        energy_1tau_n_W.append(energy)
        xsec_1tau_n_W.append(xsec_perE)
        delta_1tau_n_W.append(delta)

### vmu -> vmu mu+ mu- ; incoherent ; total ; Argon ###
energy_2mu_incoh_Ar = energy_2mu_p_Ar
xsec_2mu_incoh_Ar   = [sum(i) for i in zip(xsec_2mu_n_Ar,xsec_2mu_p_Ar)]
delta_2mu_incoh_Ar  = [np.sqrt(i**2 + j**2) for i,j in zip(delta_2mu_n_Ar, delta_2mu_p_Ar)]

### vmu -> vtau tau+ mu- ; incoherent ; total ; Argon ###
energy_1tau_incoh_Ar = energy_1tau_n_Ar
xsec_1tau_incoh_Ar   = [sum(i) for i in zip(xsec_1tau_n_Ar,xsec_1tau_p_Ar)]
delta_1tau_incoh_Ar  = [np.sqrt(i**2 + j**2) for i,j in zip(delta_1tau_n_Ar, delta_1tau_p_Ar)]

### vmu -> vmu tau+ tau- ; incoherent ; total ; Argon ###
energy_2tau_incoh_Ar = energy_2tau_n_Ar
xsec_2tau_incoh_Ar   = [sum(i) for i in zip(xsec_2tau_n_Ar,xsec_2tau_p_Ar)]
delta_2tau_incoh_Ar  = [np.sqrt(i**2 + j**2) for i,j in zip(delta_2tau_n_Ar, delta_2tau_p_Ar)]

### vmu -> vtau tau+ mu- ; incoherent ; total ; Tungsten ###
energy_1tau_incoh_W = energy_1tau_n_W
xsec_1tau_incoh_W   = [sum(i) for i in zip(xsec_1tau_n_W,xsec_1tau_p_W)]
delta_1tau_incoh_W  = [np.sqrt(i**2 + j**2) for i,j in zip(delta_1tau_n_W, delta_1tau_p_W)]

### vmu -> vmu tau+ tau- ; incoherent ; total ; Tungsten ###
energy_2tau_incoh_W = energy_2tau_n_W
xsec_2tau_incoh_W   = [sum(i) for i in zip(xsec_2tau_n_W,xsec_2tau_p_W)]
delta_2tau_incoh_W  = [np.sqrt(i**2 + j**2) for i,j in zip(delta_2tau_n_W, delta_2tau_p_W)]

##### Nucleon #####
### vmu -> vmu tau+ tau- ; nucleon ; proton ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/new_proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_p * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_2tau_p.append(energy)
        xsec_2tau_p.append(xsec_perE)
        delta_2tau_p.append(delta)

### vmu -> vtau tau+ mu- ; nucleon ; proton ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/new_proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_p * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_1tau_p.append(energy)
        xsec_1tau_p.append(xsec_perE)
        delta_1tau_p.append(delta)

### vmu -> vmu mu+ mu- ; nucleon ; proton ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/new_proton/vmu_to_vmu_mu+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_p * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_2mu_p.append(energy)
        xsec_2mu_p.append(xsec_perE)
        delta_2mu_p.append(delta)

### vmu -> vmu mu+ mu- ; nucleon ; neutron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/new_neutron/vmu_to_vmu_mu+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_n * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_2mu_n.append(energy)
        xsec_2mu_n.append(xsec_perE)
        delta_2mu_n.append(delta)

### vmu -> vmu tau+ tau- ; nucleon ; neutron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/new_neutron/vmu_to_vmu_tau+_tau-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_n * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_2tau_n.append(energy)
        xsec_2tau_n.append(xsec_perE)
        delta_2tau_n.append(delta)

### vmu -> vtau tau+ mu- ; nucleon ; neutron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/new_neutron/vmu_to_vtau_tau+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_n * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_1tau_n.append(energy)
        xsec_1tau_n.append(xsec_perE)
        delta_1tau_n.append(delta)

##### Digitized #####
### vmu -> vmu mu+ mu- ; coherent ; Argon ; Altmannshofer et al. ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/argon/vmu_to_vmu_mu+_mu-_coh_Ar_xsec_Altmannshofer.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2

        energy_2mu_coh_Ar_Alt.append(energy)
        xsec_2mu_coh_Ar_Alt.append(xsec_perE)

### vmu -> vmu mu+ mu- ; incoherent ; proton ; Argon ; Altmannshofer et al. ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/incoherent/proton/vmu_to_vmu_mu+_mu-_incoh_p_Ar_xsec_Altmannshofer.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2

        energy_2mu_incoh_p_Ar_Alt.append(energy)
        xsec_2mu_incoh_p_Ar_Alt.append(xsec_perE)

### vmu -> vmu mu+ mu- ; incoherent ; neutron ; Argon ; Altmannshofer et al. ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/incoherent/neutron/vmu_to_vmu_mu+_mu-_incoh_n_Ar_xsec_Altmannshofer.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2

        energy_2mu_incoh_n_Ar_Alt.append(energy)
        xsec_2mu_incoh_n_Ar_Alt.append(xsec_perE)

### vmu X -> mu- X' ; vmuCC ; Formaggio & Zeller ###
with open(CROSS_SECTION_DIR + '/vmuCC/vmuCC_xsec_perE_Formaggio.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) * 10 * A_Ar * fb_to_cm2 # digitized plot is in [1e-38 cm^2 / GeV]; also, xsec given in as xsec/nucleon -> xsec/Ar

        energy_numuCC_FZ.append(energy)
        xsec_numuCC_FZ.append(xsec)

### vmu -> vmu mu+ mu- ; coherent ; Argon ; Ballett et al. ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/argon/Ballett_coh_Ar.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) * (Z_Ar)**2 * 1e-44 / energy # Ballett et al. values are normalized to Z^2 1e-44 cm^2. 

        energy_2mu_coh_Ar_Ballett.append(energy)
        xsec_2mu_coh_Ar_Ballett.append(xsec_perE)

### vmu -> vmu mu+ mu- ; coherent ; Argon ; Zhou and Beacom ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/digitized/oxygen/Zhou+Beacom_2tau_coh_O.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) * (8/Z_Ar)**2  # Zhou+Beacom values are already normalized by Enu but are for oxygen. Use (Z_O/Z_Ar)^2 for rough comparison.

        energy_2tau_coh_Ar_Beacom.append(energy)
        xsec_2tau_coh_Ar_Beacom.append(xsec_perE)

### vmu -> vmu mu+ mu- ; incoherent ; Argon ; Zhou and Beacom ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/digitized/oxygen/Zhou+Beacom_2tau_incoh_O.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) * A_Ar/16  # Zhou+Beacom values are already normalized by Enu but are for oxygen. Use A_Ar/A_O for comparison since incoherent should scale by N.

        energy_2tau_incoh_Ar_Beacom.append(energy)
        xsec_2tau_incoh_Ar_Beacom.append(xsec_perE)

### vmu -> vmu mu+ mu- ; coherent ; Argon ; Magill and Plestid ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/argon/Magill_Plestid_coh_Ar.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) * (Z_Ar)**2 * 1e-45 # Ballett et al. values are normalized to Z^2*E_\nu and are in units of 1e-45 cm^2.

        energy_2mu_coh_Ar_Magill.append(energy)
        xsec_2mu_coh_Ar_Magill.append(xsec_perE)


####################
##### Plotting #####
####################

fig1, ax1 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Validation cross sections
fig2, ax2 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Tau cross sections + vmuCC
fig3, ax3 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Nucleon cross sections
fig4, ax4 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Tungsten tau cross sections

color_1tau = 'firebrick'
color_2tau = 'orange'
color_2mu  = 'navy'

### vmu -> vtau tau+ mu- ###
# Argon #
ax2.errorbar(energy_1tau_coh_Ar, xsec_1tau_coh_Ar, color = color_1tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.errorbar(energy_1tau_incoh_Ar, xsec_1tau_incoh_Ar, color = color_1tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.text(5,1e-52,r'$\tau^+ \mu^-$ {\bf Coh.}',color=color_1tau,rotation=62,fontsize=30)
ax2.text(18,8e-44,r'$\tau^+ \mu^-$ {\bf Incoh.}',color=color_1tau,rotation=10,fontsize=30)
#ax2.text(5.5,1e-46,r'$\tau^+ \mu^-$ {\bf Incoh.}',color=color_1tau,rotation=50,fontsize=30)

upper_curve_1tau_coh_Ar = [xsec + sigma for xsec, sigma in zip(xsec_1tau_coh_Ar, delta_1tau_coh_Ar)]
lower_curve_1tau_coh_Ar = [xsec - sigma for xsec, sigma in zip(xsec_1tau_coh_Ar, delta_1tau_coh_Ar)]

upper_curve_1tau_incoh_Ar = [xsec + sigma for xsec, sigma in zip(xsec_1tau_incoh_Ar, delta_1tau_incoh_Ar)]
lower_curve_1tau_incoh_Ar = [xsec - sigma for xsec, sigma in zip(xsec_1tau_incoh_Ar, delta_1tau_incoh_Ar)]

ax2.fill_between(energy_1tau_coh_Ar, lower_curve_1tau_coh_Ar, upper_curve_1tau_coh_Ar, color = color_1tau, alpha=0.50)
ax2.fill_between(energy_1tau_incoh_Ar, lower_curve_1tau_incoh_Ar, upper_curve_1tau_incoh_Ar, color = color_1tau, alpha=0.50) # due to large y-range, error seems smaller but it's not

# Tungsten #
ax4.errorbar(energy_1tau_coh_W, xsec_1tau_coh_W, color = color_1tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax4.errorbar(energy_1tau_incoh_W, xsec_1tau_incoh_W, color = color_1tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax4.text(10,5e-7,r'$\tau^+ \mu^-$ {\bf Coh.}',color=color_1tau,rotation=37,fontsize=30)
#ax4.text(10,8e-4,r'$\tau^+ \mu^-$ {\bf Incoh.}',color=color_1tau,rotation=20,fontsize=30) 

upper_curve_1tau_coh_W = [xsec + sigma for xsec, sigma in zip(xsec_1tau_coh_W, delta_1tau_coh_W)]
lower_curve_1tau_coh_W = [xsec - sigma for xsec, sigma in zip(xsec_1tau_coh_W, delta_1tau_coh_W)]

upper_curve_1tau_incoh_W = [xsec + sigma for xsec, sigma in zip(xsec_1tau_incoh_W, delta_1tau_incoh_W)]
lower_curve_1tau_incoh_W = [xsec - sigma for xsec, sigma in zip(xsec_1tau_incoh_W, delta_1tau_incoh_W)]

ax4.fill_between(energy_1tau_coh_W, lower_curve_1tau_coh_W, upper_curve_1tau_coh_W, color = color_1tau, alpha=0.50)
ax4.fill_between(energy_1tau_incoh_W, lower_curve_1tau_incoh_W, upper_curve_1tau_incoh_W, color = color_1tau, alpha=0.50) # due to large y-range, error seems smaller but it's not

# Nucleon #
ax3.errorbar(energy_1tau_p, xsec_1tau_p, color = color_1tau, linestyle = 'solid', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.errorbar(energy_1tau_n, xsec_1tau_n, color = color_1tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.text(290,3.5e-44,r'$\tau^+ \mu^-$ {\bf p}',color=color_1tau,rotation=10,fontsize=30)
ax3.text(18,1e-46,r'$\tau^+ \mu^-$ {\bf n}',color=color_1tau,rotation=33,fontsize=30)

upper_curve_1tau_p = [xsec + sigma for xsec, sigma in zip(xsec_1tau_p, delta_1tau_p)]
lower_curve_1tau_p = [xsec - sigma for xsec, sigma in zip(xsec_1tau_p, delta_1tau_p)]

upper_curve_1tau_n = [xsec + sigma for xsec, sigma in zip(xsec_1tau_n, delta_1tau_n)]
lower_curve_1tau_n = [xsec - sigma for xsec, sigma in zip(xsec_1tau_n, delta_1tau_n)]

ax3.fill_between(energy_1tau_p, lower_curve_1tau_p, upper_curve_1tau_p, color = color_1tau, alpha=0.50)
ax3.fill_between(energy_1tau_n, lower_curve_1tau_n, upper_curve_1tau_n, color = color_1tau, alpha=0.50)

### vmu -> vmu tau+ tau- ###
# Argon #
ax2.errorbar(energy_2tau_coh_Ar, xsec_2tau_coh_Ar, color = color_2tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.errorbar(energy_2tau_incoh_Ar, xsec_2tau_incoh_Ar, color = color_2tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.text(22,1e-52,r'$\tau^+ \tau^-$ {\bf Coh.}',color=color_2tau,rotation=53,fontsize=30)
ax2.text(35,1e-46,r'$\tau^+ \tau^-$ {\bf Incoh.}',color=color_2tau,rotation=19,fontsize=30)

upper_curve_2tau_coh_Ar = [xsec + sigma for xsec, sigma in zip(xsec_2tau_coh_Ar, delta_2tau_coh_Ar)]
lower_curve_2tau_coh_Ar = [xsec - sigma for xsec, sigma in zip(xsec_2tau_coh_Ar, delta_2tau_coh_Ar)]

upper_curve_2tau_incoh_Ar = [xsec + sigma for xsec, sigma in zip(xsec_2tau_incoh_Ar, delta_2tau_incoh_Ar)]
lower_curve_2tau_incoh_Ar = [xsec - sigma for xsec, sigma in zip(xsec_2tau_incoh_Ar, delta_2tau_incoh_Ar)]

ax2.fill_between(energy_2tau_coh_Ar, lower_curve_2tau_coh_Ar, upper_curve_2tau_coh_Ar, color = color_2tau, alpha=0.50)
ax2.fill_between(energy_2tau_incoh_Ar, lower_curve_2tau_incoh_Ar, upper_curve_2tau_incoh_Ar, color = color_2tau, alpha=0.50)

# Tungsten #
ax4.errorbar(energy_2tau_coh_W, xsec_2tau_coh_W, color = color_2tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax4.errorbar(energy_2tau_incoh_W, xsec_2tau_incoh_W, color = color_2tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax4.text(50,1e-7,r'$\tau^+ \tau^-$ {\bf Coh.}',color=color_2tau,rotation=34,fontsize=30)
#ax4.text(50,1e-4,r'$\tau^+ \tau^-$ {\bf Incoh.}',color=color_2tau,rotation=18,fontsize=30)

upper_curve_2tau_coh_W = [xsec + sigma for xsec, sigma in zip(xsec_2tau_coh_W, delta_2tau_coh_W)]
lower_curve_2tau_coh_W = [xsec - sigma for xsec, sigma in zip(xsec_2tau_coh_W, delta_2tau_coh_W)]

upper_curve_2tau_incoh_W = [xsec + sigma for xsec, sigma in zip(xsec_2tau_incoh_W, delta_2tau_incoh_W)]
lower_curve_2tau_incoh_W = [xsec - sigma for xsec, sigma in zip(xsec_2tau_incoh_W, delta_2tau_incoh_W)]

ax4.fill_between(energy_2tau_coh_W, lower_curve_2tau_coh_W, upper_curve_2tau_coh_W, color = color_2tau, alpha=0.50)
ax4.fill_between(energy_2tau_incoh_W, lower_curve_2tau_incoh_W, upper_curve_2tau_incoh_W, color = color_2tau, alpha=0.50)

# Nucleon #
ax3.errorbar(energy_2tau_p, xsec_2tau_p, color = color_2tau, linestyle = 'solid', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.errorbar(energy_2tau_n, xsec_2tau_n, color = color_2tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.text(290,6e-46,r'$\tau^+ \tau^-$ {\bf p}',color=color_2tau,rotation=13,fontsize=30)
ax3.text(65,3e-47,r'$\tau^+ \tau^-$ {\bf n}',color=color_2tau,rotation=28,fontsize=30)

upper_curve_2tau_p = [xsec + sigma for xsec, sigma in zip(xsec_2tau_p, delta_2tau_p)]
lower_curve_2tau_p = [xsec - sigma for xsec, sigma in zip(xsec_2tau_p, delta_2tau_p)]

upper_curve_2tau_n = [xsec + sigma for xsec, sigma in zip(xsec_2tau_n, delta_2tau_n)]
lower_curve_2tau_n = [xsec - sigma for xsec, sigma in zip(xsec_2tau_n, delta_2tau_n)]

ax3.fill_between(energy_2tau_p, lower_curve_2tau_p, upper_curve_2tau_p, color = color_2tau, alpha=0.50)
ax3.fill_between(energy_2tau_n, lower_curve_2tau_n, upper_curve_2tau_n, color = color_2tau, alpha=0.50)

### vmu -> vmu mu+ mu- ###
# Argon #
ax1.errorbar(energy_2mu_coh_Ar, xsec_2mu_coh_Ar, color = color_2mu, label='_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax1.errorbar(energy_2mu_n_Ar, xsec_2mu_incoh_Ar, color = color_2mu, label='_hidden', linestyle = 'dashed', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax1.text(50,2.1e-42,r'$\mu^+ \mu^-$ {\bf Coh.}',color=color_2mu,rotation=6,fontsize=30)
ax1.text(50,6.5e-43,r'$\mu^+ \mu^-$ {\bf Incoh.}',color=color_2mu,rotation=7,fontsize=30)

upper_curve_2mu_coh_Ar = [xsec + sigma for xsec, sigma in zip(xsec_2mu_coh_Ar, delta_2mu_coh_Ar)]
lower_curve_2mu_coh_Ar = [xsec - sigma for xsec, sigma in zip(xsec_2mu_coh_Ar, delta_2mu_coh_Ar)]

upper_curve_2mu_incoh_Ar = [xsec + sigma for xsec, sigma in zip(xsec_2mu_incoh_Ar, delta_2mu_incoh_Ar)]
lower_curve_2mu_incoh_Ar = [xsec - sigma for xsec, sigma in zip(xsec_2mu_incoh_Ar, delta_2mu_incoh_Ar)]

ax1.fill_between(energy_2mu_coh_Ar, lower_curve_2mu_coh_Ar, upper_curve_2mu_coh_Ar, color = color_2mu, alpha=0.50)
ax1.fill_between(energy_2mu_incoh_Ar, lower_curve_2mu_incoh_Ar, upper_curve_2mu_incoh_Ar, color = color_2mu, alpha=0.50)

#ax2.errorbar(energy_2mu_coh_Ar, xsec_2mu_coh_Ar, color = color_2mu, label='_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2.errorbar(energy_2mu_n_Ar, xsec_2mu_incoh_Ar, color = color_2mu, label='_hidden', linestyle = 'dashed', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2.text(50,2.1e-42,r'$\mu^+ \mu^-$ {\bf Coh.}',color=color_2mu,rotation=6,fontsize=30)
#ax2.text(50,6.5e-43,r'$\mu^+ \mu^-$ {\bf Incoh.}',color=color_2mu,rotation=5,fontsize=30)

#ax2.fill_between(energy_2mu_coh_Ar, lower_curve_2mu_coh_Ar, upper_curve_2mu_coh_Ar, color = color_2mu, alpha=0.50)
#ax2.fill_between(energy_2mu_incoh_Ar, lower_curve_2mu_incoh_Ar, upper_curve_2mu_incoh_Ar, color = color_2mu, alpha=0.50)

# Tungsten #
ax4.errorbar(energy_2mu_coh_W, xsec_2mu_coh_W, color = color_2mu, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax4.errorbar(energy_2mu_incoh_W, xsec_2mu_incoh_W, color = color_2mu, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax4.text(50,1e-7,r'$\tau^+ \tau^-$ {\bf Coh.}',color=color_2mu,rotation=34,fontsize=30)
#ax4.text(50,1e-4,r'$\tau^+ \tau^-$ {\bf Incoh.}',color=color_2mu,rotation=18,fontsize=30)

upper_curve_2mu_coh_W = [xsec + sigma for xsec, sigma in zip(xsec_2mu_coh_W, delta_2mu_coh_W)]
lower_curve_2mu_coh_W = [xsec - sigma for xsec, sigma in zip(xsec_2mu_coh_W, delta_2mu_coh_W)]

upper_curve_2mu_incoh_W = [xsec + sigma for xsec, sigma in zip(xsec_2mu_incoh_W, delta_2mu_incoh_W)]
lower_curve_2mu_incoh_W = [xsec - sigma for xsec, sigma in zip(xsec_2mu_incoh_W, delta_2mu_incoh_W)]

ax4.fill_between(energy_2mu_coh_W, lower_curve_2mu_coh_W, upper_curve_2mu_coh_W, color = color_2mu, alpha=0.50)
ax4.fill_between(energy_2mu_incoh_W, lower_curve_2mu_incoh_W, upper_curve_2mu_incoh_W, color = color_2mu, alpha=0.50)


# Nucleon #
ax3.errorbar(energy_2mu_p, xsec_2mu_p, color = color_2mu, linestyle = 'solid', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.errorbar(energy_2mu_n, xsec_2mu_n, color = color_2mu, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.text(0.8,7e-46,r'$\mu^+ \mu^-$ {\bf p}',color=color_2mu,rotation=36,fontsize=30)
ax3.text(1.7,1e-46,r'$\mu^+ \mu^-$ {\bf n}',color=color_2mu,rotation=28,fontsize=30)

upper_curve_2mu_p = [xsec + sigma for xsec, sigma in zip(xsec_2mu_p, delta_2mu_p)]
lower_curve_2mu_p = [xsec - sigma for xsec, sigma in zip(xsec_2mu_p, delta_2mu_p)]

upper_curve_2mu_n = [xsec + sigma for xsec, sigma in zip(xsec_2mu_n, delta_2mu_n)]
lower_curve_2mu_n = [xsec - sigma for xsec, sigma in zip(xsec_2mu_n, delta_2mu_n)]

ax3.fill_between(energy_2mu_p, lower_curve_2mu_p, upper_curve_2mu_p, color = color_2mu, alpha=0.50)
ax3.fill_between(energy_2mu_n, lower_curve_2mu_n, upper_curve_2mu_n, color = color_2mu, alpha=0.50)

### Digitized ###
ax1.scatter(energy_2mu_coh_Ar_Alt, xsec_2mu_coh_Ar_Alt, marker='o', color = 'gray', s=40, edgecolors='k', linewidths=1, alpha=0.75, zorder=3, label='Altmannshofer+ 2019') 
ax1.scatter(energy_2mu_coh_Ar_Ballett, xsec_2mu_coh_Ar_Ballett, marker='^', color = 'black', s=40, edgecolors='k', linewidths=1, alpha=0.75, zorder=4, label='Ballett+ 2018')
ax1.scatter(energy_2mu_coh_Ar_Magill, xsec_2mu_coh_Ar_Magill, marker='v', color = 'silver', s=40, edgecolors='k', linewidths=1, alpha=0.75, zorder=5, label=r'Magill $\&$ Plestid 2017')

#ax2.scatter(energy_2mu_coh_Ar_Alt, xsec_2mu_coh_Ar_Alt, marker='o', color = 'grey', s=40, edgecolors='k', linewidths=1, alpha=0.75, zorder=3, label='Altmannshofer+ 2019') 
#ax2.scatter(energy_2mu_coh_Ar_Ballett, xsec_2mu_coh_Ar_Ballett, marker='v', color = 'black', s=40, edgecolors='k', linewidths=1, alpha=0.75, zorder=4, label='Ballett+ 2018')

ax2.scatter(energy_2tau_incoh_Ar_Beacom, xsec_2tau_incoh_Ar_Beacom, marker='^', color = 'black', s=40, edgecolors='k', linewidths=1, alpha=0.75, zorder=4, label=r'Zhou $\&$ Beacom 2020')
ax2.plot(energy_numuCC_FZ, xsec_numuCC_FZ ,'-', color = 'grey', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.text(100,4e-38,r'$\nu_\mu$ {\bf CC}',color='black',rotation=0,fontsize=30)

### Plotting options ###
ax1.set_xlabel(r'{\bf Neutrino Energy }$E_\nu$ (GeV)') 
ax1.set_ylabel(r'{\bf Cross Section / }$E_\nu$ (cm$^{2}$ GeV$^{-1}$)') 
ax1.set_ylim(1e-47, 1e-40)
ax1.set_xlim(0.1,220)
ax1.legend(loc='lower right', fontsize=25)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_title(r'{\bf Argon $^{40}$Ar}', fontsize=40)

ax2.set_xlabel(r'{\bf Neutrino Energy }$E_\nu$ (GeV)') 
ax2.set_ylabel(r'{\bf Cross Section / }$E_\nu$ (cm$^{2}$ GeV$^{-1}$)') 
ax2.set_ylim(1e-56, 1e-36)
#ax2.set_ylim(1e-51, 1e-42)
ax2.set_xlim(0.1,11000)
ax2.legend(loc='lower right', fontsize=25)
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_title(r'{\bf Argon $^{40}$Ar}', fontsize=40)

ax3.set_xlabel(r'{\bf Neutrino Energy }$E_\nu$ (GeV)') 
ax3.set_ylabel(r'{\bf Cross Section / }$E_\nu$ (cm$^{2}$ GeV$^{-1}$)') 
ax3.set_ylim(1e-51, 1e-42)
ax3.set_xlim(0.1,11000)
#ax3.legend(loc='lower right') 
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.set_title(r'{\bf Fermi gas model}', fontsize=40)

ax4.set_xlabel(r'{\bf Neutrino Energy }$E_\nu$ (GeV)') 
ax4.set_ylabel(r'{\bf Cross Section / }$E_\nu$ (cm$^{2}$ GeV$^{-1}$)') 
ax4.set_ylim(1e-56, 1e-33)
ax4.set_xlim(0.1,11000)
#ax4.legend(loc='lower right') 
ax4.set_yscale('log')
ax4.set_xscale('log')
ax4.set_title(r'{\bf Tungsten $^{184}$W}', fontsize=40)

locmaj = mticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
locmin = mticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)
ax1.xaxis.set_major_locator(locmaj)
ax1.xaxis.set_minor_locator(locmin)
ax2.xaxis.set_major_locator(locmaj)
ax2.xaxis.set_minor_locator(locmin)
ax3.xaxis.set_major_locator(locmaj)
ax3.xaxis.set_minor_locator(locmin)
ax4.xaxis.set_major_locator(locmaj)
ax4.xaxis.set_minor_locator(locmin)


xmajor1 = [0.1, 1, 5, 10, 50, 100]
xmajor2 = [0.1, 1, 5, 10, 50, 100, 1000, 10000]
xmajor3 = [0.1, 1, 5, 10, 50, 100, 1000, 10000]
xmajor4 = [0.1, 1, 5, 10, 50, 100, 1000, 10000]

ax1.set_xticks(xmajor1, labels=['0.1', '1', '5', '10', '50', '100'])
ax1.grid(which='major', axis='both')
ax2.set_xticks(xmajor2, labels=['0.1', '1', '5', '10', '50', '100', '1000', '10000'])
ax2.grid(which='major', axis='both')
ax3.set_xticks(xmajor3, labels=['0.1', '1', '5', '10', '50', '100', '1000', '10000'])
ax3.grid(which='major', axis='both')
ax4.set_xticks(xmajor4, labels=['0.1', '1', '5', '10', '50', '100', '1000', '10000'])
ax4.grid(which='major', axis='both')

fig1.savefig("../plots/xsec_perE_mu_argon.png", dpi=400, bbox_inches='tight')
fig2.savefig("../plots/xsec_perE_tau_argon.png", dpi=400, bbox_inches='tight')
fig3.savefig("../plots/xsec_perE_nucleons.png", dpi=400, bbox_inches='tight')
fig4.savefig("../plots/xsec_perE_tau_tungsten.png", dpi=400, bbox_inches='tight')
