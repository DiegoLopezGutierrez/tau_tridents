import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, DrawingArea
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

Z_Fe = 26
A_Fe = 56

## Flux Ranges ##
DUNE_std_low = 0.125
DUNE_std_high = 60.0

DUNE_tau_opt_low = 1

## Cross Section Uncertainties ##
# Components # (see Altmannshofer et al. for details)
aEM = 1/137  # low q^2 usually, so using zero momentum value of fine structure
sigma_highQED_Ar = Z_Ar*aEM/(4*np.pi) + 0.02  # higher order QED corrections should give roughly Z*aEM/(4*pi) which is about 1% for Ar. Add 2% to be conservative.
sigma_highQED_W  = Z_W*aEM/(4*np.pi) + 0.02 # higher order QED corrections should give roughly Z*aEM/(4*pi) which is about 4% for W. Add 2% to be conservative.
sigma_highQED_Fe  = Z_Fe*aEM/(4*np.pi) + 0.02 # higher order QED corrections should give roughly Z*aEM/(4*pi) which is about 4% for W. Add 2% to be conservative.
sigma_highQED_p = 1*aEM/(4*np.pi) + 0.02
sigma_highQED_n = 0

sigma_form_factors_coh = 0.01 # 1% uncertainty on coherent cross sections due to form factors
sigma_form_factors_incoh = 0.03 # 3% uncertainty on incoherent cross sections due to form factors

sigma_highEW = 0.05 # 5% uncertainty due to higher order weak corrections

sigma_nuclear_modeling = 0.30 # for incoherent scattering, largest uncertainty due to other nuclear effects besides the included Pauli blocking. 30% to be conservative.

# Total #
sigma_total_coh_Ar = np.sqrt(sigma_highQED_Ar**2 + sigma_form_factors_coh**2 + sigma_highEW**2)
sigma_total_coh_W  = np.sqrt(sigma_highQED_W**2 + sigma_form_factors_coh**2 + sigma_highEW**2)
sigma_total_coh_Fe  = np.sqrt(sigma_highQED_Fe**2 + sigma_form_factors_coh**2 + sigma_highEW**2)

sigma_total_incoh_Ar = np.sqrt(sigma_highQED_Ar**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)
sigma_total_incoh_W  = np.sqrt(sigma_highQED_W**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)
sigma_total_incoh_Fe  = np.sqrt(sigma_highQED_Fe**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)

sigma_total_p = np.sqrt(sigma_highQED_p**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)
sigma_total_n = np.sqrt(sigma_highQED_n**2 + sigma_form_factors_incoh**2 + sigma_highEW**2 + sigma_nuclear_modeling**2)

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

energy_ve1tau_coh_Ar = []
xsec_ve1tau_coh_Ar = []
delta_ve1tau_coh_Ar = []

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

energy_ve1tau_coh_W = []
xsec_ve1tau_coh_W = []
delta_ve1tau_coh_W = []

### Iron ###
energy_1tau_coh_Fe = [] 
xsec_1tau_coh_Fe = []
delta_1tau_coh_Fe = []

energy_2tau_coh_Fe = []
xsec_2tau_coh_Fe = []
delta_2tau_coh_Fe = []

energy_2mu_coh_Fe = []
xsec_2mu_coh_Fe = []
delta_2mu_coh_Fe = []

energy_ve1tau_coh_Fe = []
xsec_ve1tau_coh_Fe = []
delta_ve1tau_coh_Fe = []

############################
##### Incoherent xsec ######
############################

##### Argon #####
### Proton ###
energy_1tau_p_Ar = []
xsec_1tau_p_Ar = []
delta_1tau_p_Ar = []

energy_2tau_p_Ar = []
xsec_2tau_p_Ar = []
delta_2tau_p_Ar = []

energy_2mu_p_Ar = []
xsec_2mu_p_Ar = []
delta_2mu_p_Ar = []

energy_ve1tau_p_Ar = []
xsec_ve1tau_p_Ar = []
delta_ve1tau_p_Ar = []

### Neutron ###
energy_1tau_n_Ar = []
xsec_1tau_n_Ar = []
delta_1tau_n_Ar = []

energy_2tau_n_Ar = []
xsec_2tau_n_Ar = []
delta_2tau_n_Ar = []

energy_2mu_n_Ar = []
xsec_2mu_n_Ar = []
delta_2mu_n_Ar = []

energy_ve1tau_n_Ar = []
xsec_ve1tau_n_Ar = []
delta_ve1tau_n_Ar = []

### Proton + Neutron ###
energy_1tau_incoh_Ar = []
xsec_1tau_incoh_Ar = []
delta_1tau_incoh_Ar = []

energy_2tau_incoh_Ar = []
xsec_2tau_incoh_Ar = []
delta_2tau_incoh_Ar = []

energy_2mu_incoh_Ar = []
xsec_2mu_incoh_Ar = []
delta_2mu_incoh_Ar = []

energy_ve1tau_incoh_Ar = []
xsec_ve1tau_incoh_Ar = []
delta_ve1tau_incoh_Ar = []

##### Tungsten #####
### Proton ###
energy_1tau_p_W = []
xsec_1tau_p_W = []
delta_1tau_p_W = []

energy_2tau_p_W = []
xsec_2tau_p_W = []
delta_2tau_p_W = []

energy_2mu_p_W = []
xsec_2mu_p_W = []
delta_2mu_p_W = []

energy_ve1tau_p_W = []
xsec_ve1tau_p_W = []
delta_ve1tau_p_W = []

### Neutron ###
energy_1tau_n_W = []
xsec_1tau_n_W = []
delta_1tau_n_W = []

energy_2tau_n_W = []
xsec_2tau_n_W = []
delta_2tau_n_W = []

energy_2mu_n_W = []
xsec_2mu_n_W = []
delta_2mu_n_W = []

energy_ve1tau_n_W = []
xsec_ve1tau_n_W = []
delta_ve1tau_n_W = []

### Proton + Neutron ###
energy_1tau_incoh_W = []
xsec_1tau_incoh_W = []
delta_1tau_incoh_W = []

energy_2tau_incoh_W = []
xsec_2tau_incoh_W = []
delta_2tau_incoh_W = []

energy_2mu_incoh_W = []
xsec_2mu_incoh_W = []
delta_2mu_incoh_W = []

energy_ve1tau_incoh_W = []
xsec_ve1tau_incoh_W = []
delta_ve1tau_incoh_W = []

##### Iron #####
### Proton ###
energy_1tau_p_Fe = []
xsec_1tau_p_Fe = []
delta_1tau_p_Fe = []

energy_2tau_p_Fe = []
xsec_2tau_p_Fe = []
delta_2tau_p_Fe = []

energy_2mu_p_Fe = []
xsec_2mu_p_Fe = []
delta_2mu_p_Fe = []

energy_ve1tau_p_Fe = []
xsec_ve1tau_p_Fe = []
delta_ve1tau_p_Fe = []

### Neutron ###
energy_1tau_n_Fe = []
xsec_1tau_n_Fe = []
delta_1tau_n_Fe = []

energy_2tau_n_Fe = []
xsec_2tau_n_Fe = []
delta_2tau_n_Fe = []

energy_2mu_n_Fe = []
xsec_2mu_n_Fe = []
delta_2mu_n_Fe = []

energy_ve1tau_n_Fe = []
xsec_ve1tau_n_Fe = []
delta_ve1tau_n_Fe = []

### Proton + Neutron ###
energy_1tau_incoh_Fe = []
xsec_1tau_incoh_Fe = []
delta_1tau_incoh_Fe = []

energy_2tau_incoh_Fe = []
xsec_2tau_incoh_Fe = []
delta_2tau_incoh_Fe = []

energy_2mu_incoh_Fe = []
xsec_2mu_incoh_Fe = []
delta_2mu_incoh_Fe = []

energy_ve1tau_incoh_Fe = []
xsec_ve1tau_incoh_Fe = []
delta_ve1tau_incoh_Fe = []

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

energy_ve1tau_p = []
xsec_ve1tau_p = []
delta_ve1tau_p = []

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

energy_ve1tau_n = []
xsec_ve1tau_n = []
delta_ve1tau_n = []

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

### ve -> vtau tau+ e- ; coherent ; Argon ###
with open(CROSS_SECTION_DIR + '/ve_to_vtau_tau+_e-_xsec/coherent/argon/ve_to_vtau_tau+_e-_coh_Ar_xsec.csv','r') as csvfile: 
    data = csv.reader(csvfile, delimiter = ',') 
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_ve1tau_coh_Ar.append(energy) 
        xsec_ve1tau_coh_Ar.append(xsec_perE)
        delta_ve1tau_coh_Ar.append(delta)

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

### ve -> vtau tau+ e- ; coherent ; Tungsten ###
with open(CROSS_SECTION_DIR + '/ve_to_vtau_tau+_e-_xsec/coherent/tungsten/ve_to_vtau_tau+_e-_coh_W_xsec.csv','r') as csvfile: 
    data = csv.reader(csvfile, delimiter = ',') 
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_W * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_ve1tau_coh_W.append(energy) 
        xsec_ve1tau_coh_W.append(xsec_perE)
        delta_ve1tau_coh_W.append(delta)

### vmu -> vtau tau+ mu- ; coherent ; Iron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/coherent/iron/vmu_to_vtau_tau+_mu-_coh_Fe_xsec.csv','r') as csvfile: 
    data = csv.reader(csvfile, delimiter = ',') 
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_1tau_coh_Fe.append(energy)
        xsec_1tau_coh_Fe.append(xsec_perE)
        delta_1tau_coh_Fe.append(delta)

### vmu -> vmu tau+ tau- ; coherent ; Iron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/coherent/iron/vmu_to_vmu_tau+_tau-_coh_Fe_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_2tau_coh_Fe.append(energy)
        xsec_2tau_coh_Fe.append(xsec_perE)
        delta_2tau_coh_Fe.append(delta)

### vmu -> vmu mu+ mu- ; coherent ; Iron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/coherent/iron/vmu_to_vmu_mu+_mu-_coh_Fe_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_2mu_coh_Fe.append(energy)
        xsec_2mu_coh_Fe.append(xsec_perE)
        delta_2mu_coh_Fe.append(delta)

### ve -> vtau tau+ e- ; coherent ; Iron ###
with open(CROSS_SECTION_DIR + '/ve_to_vtau_tau+_e-_xsec/coherent/iron/ve_to_vtau_tau+_e-_coh_Fe_xsec.csv','r') as csvfile: 
    data = csv.reader(csvfile, delimiter = ',') 
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_coh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_ve1tau_coh_Fe.append(energy) 
        xsec_ve1tau_coh_Fe.append(xsec_perE)
        delta_ve1tau_coh_Fe.append(delta)

##### Incoherent #####
### vmu -> vmu tau+ tau- ; incoherent ; proton ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv','r') as csvfile:
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
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv','r') as csvfile:
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
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/proton/vmu_to_vmu_mu+_mu-_nucleon_p_xsec.csv','r') as csvfile:
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

### ve -> vtau tau+ e- ; incoherent ; proton ; Argon ###
with open(CROSS_SECTION_DIR + '/ve_to_vtau_tau+_e-_xsec/nucleon/proton/ve_to_vtau_tau+_e-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * Z_Ar * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * Z_Ar * fb_to_cm2

        energy_ve1tau_p_Ar.append(energy)
        xsec_ve1tau_p_Ar.append(xsec_perE)
        delta_ve1tau_p_Ar.append(delta)

### vmu -> vmu mu+ mu- ; incoherent ; neutron ; Argon ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/neutron/vmu_to_vmu_mu+_mu-_nucleon_n_xsec.csv','r') as csvfile:
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
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/neutron/vmu_to_vmu_tau+_tau-_nucleon_n_xsec.csv','r') as csvfile:
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
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/neutron/vmu_to_vtau_tau+_mu-_nucleon_n_xsec.csv','r') as csvfile:
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

### ve -> vtau tau+ e- ; incoherent ; neutron ; Argon ###
with open(CROSS_SECTION_DIR + '/ve_to_vtau_tau+_e-_xsec/nucleon/neutron/ve_to_vtau_tau+_e-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * (A_Ar-Z_Ar) * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Ar * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * (A_Ar-Z_Ar) * fb_to_cm2

        energy_ve1tau_n_Ar.append(energy)
        xsec_ve1tau_n_Ar.append(xsec_perE)
        delta_ve1tau_n_Ar.append(delta)

### vmu -> vmu tau+ tau- ; incoherent ; proton ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv','r') as csvfile:
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
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv','r') as csvfile:
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

### vmu -> vmu mu+ mu- ; incoherent ; proton ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/proton/vmu_to_vmu_mu+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * Z_W * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_W * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * Z_W * fb_to_cm2

        energy_2mu_p_W.append(energy)
        xsec_2mu_p_W.append(xsec_perE)
        delta_2mu_p_W.append(delta)

### ve -> vtau tau+ e- ; incoherent ; proton ; Tungsten ###
with open(CROSS_SECTION_DIR + '/ve_to_vtau_tau+_e-_xsec/nucleon/proton/ve_to_vtau_tau+_e-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * Z_W * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_W * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * Z_W * fb_to_cm2

        energy_ve1tau_p_W.append(energy)
        xsec_ve1tau_p_W.append(xsec_perE)
        delta_ve1tau_p_W.append(delta)

### vmu -> vmu tau+ tau- ; incoherent ; neutron ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/neutron/vmu_to_vmu_tau+_tau-_nucleon_n_xsec.csv','r') as csvfile:
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
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/neutron/vmu_to_vtau_tau+_mu-_nucleon_n_xsec.csv','r') as csvfile:
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

### vmu -> vmu mu+ mu- ; incoherent ; neutron ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/neutron/vmu_to_vmu_mu+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * (A_W-Z_W) * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_W * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * (A_W-Z_W) * fb_to_cm2

        energy_2mu_n_W.append(energy)
        xsec_2mu_n_W.append(xsec_perE)
        delta_2mu_n_W.append(delta)

### ve -> vtau tau+ e- ; incoherent ; neutron ; Tungsten ###
with open(CROSS_SECTION_DIR + '/ve_to_vtau_tau+_e-_xsec/nucleon/neutron/ve_to_vtau_tau+_e-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * (A_W-Z_W) * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_W * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * (A_W-Z_W) * fb_to_cm2

        energy_ve1tau_n_W.append(energy)
        xsec_ve1tau_n_W.append(xsec_perE)
        delta_ve1tau_n_W.append(delta)

### vmu -> vtau tau+ mu- ; incoherent ; proton ; Iron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * Z_Fe * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * Z_Fe * fb_to_cm2

        energy_1tau_p_Fe.append(energy)
        xsec_1tau_p_Fe.append(xsec_perE)
        delta_1tau_p_Fe.append(delta)

### vmu -> vmu tau+ tau- ; incoherent ; proton ; Tungsten ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * Z_Fe * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * Z_Fe * fb_to_cm2

        energy_2tau_p_Fe.append(energy)
        xsec_2tau_p_Fe.append(xsec_perE)
        delta_2tau_p_Fe.append(delta)

### vmu -> vmu mu+ mu- ; incoherent ; proton ; Iron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/proton/vmu_to_vmu_mu+_mu-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * Z_Fe * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * Z_Fe * fb_to_cm2

        energy_2mu_p_Fe.append(energy)
        xsec_2mu_p_Fe.append(xsec_perE)
        delta_2mu_p_Fe.append(delta)

### ve -> vtau tau+ e- ; incoherent ; proton ; Iron ###
with open(CROSS_SECTION_DIR + '/ve_to_vtau_tau+_e-_xsec/nucleon/proton/ve_to_vtau_tau+_e-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * Z_Fe * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * Z_Fe * fb_to_cm2

        energy_ve1tau_p_Fe.append(energy)
        xsec_ve1tau_p_Fe.append(xsec_perE)
        delta_ve1tau_p_Fe.append(delta)

### vmu -> vtau tau+ mu- ; incoherent ; neutron ; Iron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/neutron/vmu_to_vtau_tau+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * (A_Fe-Z_Fe) * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * (A_Fe-Z_Fe) * fb_to_cm2

        energy_1tau_n_Fe.append(energy)
        xsec_1tau_n_Fe.append(xsec_perE)
        delta_1tau_n_Fe.append(delta)

### vmu -> vmu tau+ tau- ; incoherent ; neutron ; Iron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/neutron/vmu_to_vmu_tau+_tau-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * (A_Fe-Z_Fe) * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * (A_Fe-Z_Fe) * fb_to_cm2

        energy_2tau_n_Fe.append(energy)
        xsec_2tau_n_Fe.append(xsec_perE)
        delta_2tau_n_Fe.append(delta)

### vmu -> vmu mu+ mu- ; incoherent ; neutron ; Iron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/neutron/vmu_to_vmu_mu+_mu-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * (A_Fe-Z_Fe) * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * (A_Fe-Z_Fe) * fb_to_cm2

        energy_2mu_n_Fe.append(energy)
        xsec_2mu_n_Fe.append(xsec_perE)
        delta_2mu_n_Fe.append(delta)

### ve -> vtau tau+ e- ; incoherent ; neutron ; Iron ###
with open(CROSS_SECTION_DIR + '/ve_to_vtau_tau+_e-_xsec/nucleon/neutron/ve_to_vtau_tau+_e-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * (A_Fe-Z_Fe) * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_incoh_Fe * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * (A_Fe-Z_Fe) * fb_to_cm2

        energy_ve1tau_n_Fe.append(energy)
        xsec_ve1tau_n_Fe.append(xsec_perE)
        delta_ve1tau_n_Fe.append(delta)


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

### ve -> vtau tau+ e- ; incoherent ; total ; Argon ###
energy_ve1tau_incoh_Ar = energy_ve1tau_n_Ar
xsec_ve1tau_incoh_Ar   = [sum(i) for i in zip(xsec_ve1tau_n_Ar,xsec_ve1tau_p_Ar)]
delta_ve1tau_incoh_Ar  = [np.sqrt(i**2 + j**2) for i,j in zip(delta_ve1tau_n_Ar, delta_ve1tau_p_Ar)]

### vmu -> vtau tau+ mu- ; incoherent ; total ; Tungsten ###
energy_1tau_incoh_W = energy_1tau_n_W
xsec_1tau_incoh_W   = [sum(i) for i in zip(xsec_1tau_n_W,xsec_1tau_p_W)]
delta_1tau_incoh_W  = [np.sqrt(i**2 + j**2) for i,j in zip(delta_1tau_n_W, delta_1tau_p_W)]

### vmu -> vmu tau+ tau- ; incoherent ; total ; Tungsten ###
energy_2tau_incoh_W = energy_2tau_n_W
xsec_2tau_incoh_W   = [sum(i) for i in zip(xsec_2tau_n_W,xsec_2tau_p_W)]
delta_2tau_incoh_W  = [np.sqrt(i**2 + j**2) for i,j in zip(delta_2tau_n_W, delta_2tau_p_W)]

### vmu -> vmu mu+ mu- ; incoherent ; total ; Tungsten ###
energy_2mu_incoh_W = energy_2mu_p_W
xsec_2mu_incoh_W   = [sum(i) for i in zip(xsec_2mu_n_W,xsec_2mu_p_W)]
delta_2mu_incoh_W  = [np.sqrt(i**2 + j**2) for i,j in zip(delta_2mu_n_W, delta_2mu_p_W)]

### ve -> vtau tau+ e- ; incoherent ; total ; Tungsten ###
energy_ve1tau_incoh_W = energy_ve1tau_n_W
xsec_ve1tau_incoh_W   = [sum(i) for i in zip(xsec_ve1tau_n_W,xsec_ve1tau_p_W)]
delta_ve1tau_incoh_W  = [np.sqrt(i**2 + j**2) for i,j in zip(delta_ve1tau_n_W, delta_ve1tau_p_W)]

### vmu -> vtau tau+ mu- ; incoherent ; total ; Iron ###
energy_1tau_incoh_Fe = energy_1tau_n_Fe
xsec_1tau_incoh_Fe   = [sum(i) for i in zip(xsec_1tau_n_Fe,xsec_1tau_p_Fe)]
delta_1tau_incoh_Fe  = [np.sqrt(i**2 + j**2) for i,j in zip(delta_1tau_n_Fe, delta_1tau_p_Fe)]

### vmu -> vmu tau+ tau- ; incoherent ; total ; Iron ###
energy_2tau_incoh_Fe = energy_2tau_n_Fe
xsec_2tau_incoh_Fe   = [sum(i) for i in zip(xsec_2tau_n_Fe,xsec_2tau_p_Fe)]
delta_2tau_incoh_Fe  = [np.sqrt(i**2 + j**2) for i,j in zip(delta_2tau_n_Fe, delta_2tau_p_Fe)]

### vmu -> vmu mu+ mu- ; incoherent ; total ; Iron ###
energy_2mu_incoh_Fe = energy_2mu_p_Fe
xsec_2mu_incoh_Fe   = [sum(i) for i in zip(xsec_2mu_n_Fe,xsec_2mu_p_Fe)]
delta_2mu_incoh_Fe  = [np.sqrt(i**2 + j**2) for i,j in zip(delta_2mu_n_Fe, delta_2mu_p_Fe)]

### ve -> vtau tau+ e- ; incoherent ; total ; Iron ###
energy_ve1tau_incoh_Fe = energy_ve1tau_n_Fe
xsec_ve1tau_incoh_Fe   = [sum(i) for i in zip(xsec_ve1tau_n_Fe,xsec_ve1tau_p_Fe)]
delta_ve1tau_incoh_Fe  = [np.sqrt(i**2 + j**2) for i,j in zip(delta_ve1tau_n_Fe, delta_ve1tau_p_Fe)]


##### Nucleon #####
### vmu -> vmu tau+ tau- ; nucleon ; proton ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/proton/vmu_to_vmu_tau+_tau-_nucleon_p_xsec.csv','r') as csvfile:
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
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/proton/vmu_to_vtau_tau+_mu-_nucleon_p_xsec.csv','r') as csvfile:
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
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/proton/vmu_to_vmu_mu+_mu-_nucleon_p_xsec.csv','r') as csvfile:
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

### ve -> vtau tau+ e- ; nucleon ; proton ###
with open(CROSS_SECTION_DIR + '/ve_to_vtau_tau+_e-_xsec/nucleon/proton/ve_to_vtau_tau+_e-_nucleon_p_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_p * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_ve1tau_p.append(energy)
        xsec_ve1tau_p.append(xsec_perE)
        delta_ve1tau_p.append(delta)

### vmu -> vmu mu+ mu- ; nucleon ; neutron ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_mu+_mu-_xsec/nucleon/neutron/vmu_to_vmu_mu+_mu-_nucleon_n_xsec.csv','r') as csvfile:
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
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/nucleon/neutron/vmu_to_vmu_tau+_tau-_nucleon_n_xsec.csv','r') as csvfile:
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
with open(CROSS_SECTION_DIR + '/vmu_to_vtau_tau+_mu-_xsec/nucleon/neutron/vmu_to_vtau_tau+_mu-_nucleon_n_xsec.csv','r') as csvfile:
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

### ve -> vtau tau+ e- ; nucleon ; neutron ###
with open(CROSS_SECTION_DIR + '/ve_to_vtau_tau+_e-_xsec/nucleon/neutron/ve_to_vtau_tau+_e-_nucleon_n_xsec.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) / energy * fb_to_cm2
        numerical_delta = float(row[2])
        physical_delta = sigma_total_n * float(row[1])
        delta = (numerical_delta + physical_delta) / energy * fb_to_cm2

        energy_ve1tau_n.append(energy)
        xsec_ve1tau_n.append(xsec_perE)
        delta_ve1tau_n.append(delta)

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

### vmu -> vmu tau+ tau- ; coherent ; Argon ; Zhou and Beacom ###
with open(CROSS_SECTION_DIR + '/vmu_to_vmu_tau+_tau-_xsec/digitized/oxygen/Zhou+Beacom_2tau_coh_O.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec_perE = float(row[1]) * (8/Z_Ar)**2  # Zhou+Beacom values are already normalized by Enu but are for oxygen. Use (Z_O/Z_Ar)^2 for rough comparison.

        energy_2tau_coh_Ar_Beacom.append(energy)
        xsec_2tau_coh_Ar_Beacom.append(xsec_perE)

### vmu -> vmu tau+ tau- ; incoherent ; Argon ; Zhou and Beacom ###
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
fig2, ax2 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Argon cross sections + vmuCC
fig3, ax3 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Nucleon cross sections
fig4, ax4 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Tungsten cross sections + vmuCC
fig5, ax5 = plt.subplots(1, 1, figsize=(15,12), tight_layout=True)  ## Iron cross sections + vmuCC
fig6, (ax61, ax62, ax63) = plt.subplots(1, 3, figsize=(50,25), tight_layout=True, sharey=True, gridspec_kw={'wspace':0})  ## Combine multi-plot 

#color_1tau   = 'firebrick'
#color_2tau   = 'orange'
#color_2mu    = 'navy'
#color_ve1tau = 'darkgreen'

# Color-blind friendly palette: https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40
color_1tau   = '#DC3220'
#color_2tau   = '#FFB50A'
color_2tau   = '#FFA500'
color_2mu    = '#005AB5'
color_ve1tau = '#752A64'

def line_style_annotation(ax, text, x, y, linestyle, color='black', xshift=0.1, fontsize=20):
    line = Line2D([0,25], [0,0], linestyle=linestyle, color=color)
    line_box = DrawingArea(30,10,0,0)
    line_box.add_artist(line)
    annotation = AnnotationBbox(line_box, (x,y), frameon=False, box_alignment=(0, 0))
    ax.add_artist(annotation)
    ax.text(x + xshift, y, text, verticalalignment='center', fontsize=fontsize, color=color)

# Linestyle Annotation - Argon #
line_style_annotation(ax2, r'\textbf{Coherent}', 0.2, 1e-38, '-')
line_style_annotation(ax2, r'\textbf{Incoherent}', 0.2, 1e-39, '--')

# Linestyle Annotation - Tungsten #
line_style_annotation(ax4, r'\textbf{Coherent}', 0.2, 1e-40, '-', xshift=0.2)
line_style_annotation(ax4, r'\textbf{Incoherent}', 0.2, 1e-41, '--',xshift=0.2)

# Linestyle Annotation - Iron #
line_style_annotation(ax5, r'\textbf{Coherent}', 0.2, 1e-41, '-')
line_style_annotation(ax5, r'\textbf{Incoherent}', 0.2, 1e-42, '--')

# Linestyle Annotation - All #
line_style_annotation(ax61, r'\textbf{Coherent}', 0.2, 1e-38, '-', fontsize=40)
line_style_annotation(ax61, r'\textbf{Incoherent}', 0.2, 1e-39, '--', fontsize=40)

### vmu -> vtau tau+ mu- ###
# Argon #
ax2.errorbar(energy_1tau_coh_Ar, xsec_1tau_coh_Ar, color = color_1tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.errorbar(energy_1tau_incoh_Ar, xsec_1tau_incoh_Ar, color = color_1tau, linestyle = 'dashed', label ='_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2.text(5,1e-52,r'$\tau^+ \mu^-$ {\bf Coh.}',color=color_1tau,rotation=62,fontsize=20)
#ax2.text(18,8e-44,r'$\tau^+ \mu^-$ {\bf Incoh.}',color=color_1tau,rotation=10,fontsize=20)
#ax2.text(5.5,1e-46,r'$\tau^+ \mu^-$ {\bf Incoh.}',color=color_1tau,rotation=50,fontsize=30)
txt2_1tau = ax2.text(2,1e-39,r'\textbf{$\nu_\mu$ Ar $\to \nu_\tau \tau^+ \mu^-$ Ar}',color=color_1tau,fontsize=25)
txt2_1tau.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

upper_curve_1tau_coh_Ar = [xsec + sigma for xsec, sigma in zip(xsec_1tau_coh_Ar, delta_1tau_coh_Ar)]
lower_curve_1tau_coh_Ar = [xsec - sigma for xsec, sigma in zip(xsec_1tau_coh_Ar, delta_1tau_coh_Ar)]

upper_curve_1tau_incoh_Ar = [xsec + sigma for xsec, sigma in zip(xsec_1tau_incoh_Ar, delta_1tau_incoh_Ar)]
lower_curve_1tau_incoh_Ar = [xsec - sigma for xsec, sigma in zip(xsec_1tau_incoh_Ar, delta_1tau_incoh_Ar)]

ax2.fill_between(energy_1tau_coh_Ar, lower_curve_1tau_coh_Ar, upper_curve_1tau_coh_Ar, color = color_1tau, alpha=0.50)
ax2.fill_between(energy_1tau_incoh_Ar, lower_curve_1tau_incoh_Ar, upper_curve_1tau_incoh_Ar, color = color_1tau, alpha=0.50) # due to large y-range, error seems smaller but it's not

# Tungsten #
ax4.errorbar(energy_1tau_coh_W, xsec_1tau_coh_W, color = color_1tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax4.errorbar(energy_1tau_incoh_W, xsec_1tau_incoh_W, color = color_1tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax4.text(10,1e-45,r'$\tau^+ \mu^-$ {\bf Coh.}',color=color_1tau,rotation=25,fontsize=30)
#ax4.text(5,2e-52,r'$\tau^+ \mu^-$ {\bf Incoh.}',color=color_1tau,rotation=60,fontsize=30)
txt4_1tau = ax4.text(300,1e-51,r'\textbf{$\nu_\mu$ W $\to \nu_\tau \tau^+ \mu^-$ W}',color=color_1tau,fontsize=25)
txt4_1tau.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

upper_curve_1tau_coh_W = [xsec + sigma for xsec, sigma in zip(xsec_1tau_coh_W, delta_1tau_coh_W)]
lower_curve_1tau_coh_W = [xsec - sigma for xsec, sigma in zip(xsec_1tau_coh_W, delta_1tau_coh_W)]

upper_curve_1tau_incoh_W = [xsec + sigma for xsec, sigma in zip(xsec_1tau_incoh_W, delta_1tau_incoh_W)]
lower_curve_1tau_incoh_W = [xsec - sigma for xsec, sigma in zip(xsec_1tau_incoh_W, delta_1tau_incoh_W)]

ax4.fill_between(energy_1tau_coh_W, lower_curve_1tau_coh_W, upper_curve_1tau_coh_W, color = color_1tau, alpha=0.50)
ax4.fill_between(energy_1tau_incoh_W, lower_curve_1tau_incoh_W, upper_curve_1tau_incoh_W, color = color_1tau, alpha=0.50) # due to large y-range, error seems smaller but it's not

# Iron #
ax5.errorbar(energy_1tau_coh_Fe, xsec_1tau_coh_Fe, color = color_1tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax5.errorbar(energy_1tau_incoh_Fe, xsec_1tau_incoh_Fe, color = color_1tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax5.text(10,1e-45,r'$\tau^+ \mu^-$ {\bf Coh.}',color=color_1tau,rotation=25,fontsize=30)
#ax5.text(5,2e-52,r'$\tau^+ \mu^-$ {\bf Incoh.}',color=color_1tau,rotation=60,fontsize=30)
txt5_1tau = ax5.text(20,1e-51,r'\textbf{$\nu_\mu$ Fe $\to \nu_\tau \tau^+ \mu^-$ Fe}',color=color_1tau,fontsize=25)
txt5_1tau.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

upper_curve_1tau_coh_Fe = [xsec + sigma for xsec, sigma in zip(xsec_1tau_coh_Fe, delta_1tau_coh_Fe)]
lower_curve_1tau_coh_Fe = [xsec - sigma for xsec, sigma in zip(xsec_1tau_coh_Fe, delta_1tau_coh_Fe)]

upper_curve_1tau_incoh_Fe = [xsec + sigma for xsec, sigma in zip(xsec_1tau_incoh_Fe, delta_1tau_incoh_Fe)]
lower_curve_1tau_incoh_Fe = [xsec - sigma for xsec, sigma in zip(xsec_1tau_incoh_Fe, delta_1tau_incoh_Fe)]

ax5.fill_between(energy_1tau_coh_Fe, lower_curve_1tau_coh_Fe, upper_curve_1tau_coh_Fe, color = color_1tau, alpha=0.50)
ax5.fill_between(energy_1tau_incoh_Fe, lower_curve_1tau_incoh_Fe, upper_curve_1tau_incoh_Fe, color = color_1tau, alpha=0.50) # due to large y-range, error seems smaller but it's not

# Nucleon #
ax3.errorbar(energy_1tau_p, xsec_1tau_p, color = color_1tau, linestyle = 'solid', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.errorbar(energy_1tau_n, xsec_1tau_n, color = color_1tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.text(290,1.2e-44,r'$\tau^+ \mu^-$ {\bf p}',color=color_1tau,rotation=10,fontsize=30)
ax3.text(18,1e-46,r'$\tau^+ \mu^-$ {\bf n}',color=color_1tau,rotation=33,fontsize=30)

upper_curve_1tau_p = [xsec + sigma for xsec, sigma in zip(xsec_1tau_p, delta_1tau_p)]
lower_curve_1tau_p = [xsec - sigma for xsec, sigma in zip(xsec_1tau_p, delta_1tau_p)]

upper_curve_1tau_n = [xsec + sigma for xsec, sigma in zip(xsec_1tau_n, delta_1tau_n)]
lower_curve_1tau_n = [xsec - sigma for xsec, sigma in zip(xsec_1tau_n, delta_1tau_n)]

ax3.fill_between(energy_1tau_p, lower_curve_1tau_p, upper_curve_1tau_p, color = color_1tau, alpha=0.50)
ax3.fill_between(energy_1tau_n, lower_curve_1tau_n, upper_curve_1tau_n, color = color_1tau, alpha=0.50)

## All ##
# Argon #
ax61.errorbar(energy_1tau_coh_Ar, xsec_1tau_coh_Ar, color = color_1tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax61.errorbar(energy_1tau_incoh_Ar, xsec_1tau_incoh_Ar, color = color_1tau, linestyle = 'dashed', label ='_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
txt61_1tau = ax61.text(3,1e-39,r'\textbf{$\nu_\mu$ $\to \nu_\tau \tau^+ \mu^-$}',color=color_1tau,fontsize=50)
txt61_1tau.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

ax61.fill_between(energy_1tau_coh_Ar, lower_curve_1tau_coh_Ar, upper_curve_1tau_coh_Ar, color = color_1tau, alpha=0.50)
ax61.fill_between(energy_1tau_incoh_Ar, lower_curve_1tau_incoh_Ar, upper_curve_1tau_incoh_Ar, color = color_1tau, alpha=0.50) # due to large y-range, error seems smaller but it's not

# Tungsten #
ax62.errorbar(energy_1tau_coh_W, xsec_1tau_coh_W, color = color_1tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax62.errorbar(energy_1tau_incoh_W, xsec_1tau_incoh_W, color = color_1tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#txt62_1tau = ax62.text(300,1e-51,r'\textbf{$\nu_\mu$ W $\to \nu_\tau \tau^+ \mu^-$ W}',color=color_1tau,fontsize=25)
#txt62_1tau.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

ax62.fill_between(energy_1tau_coh_W, lower_curve_1tau_coh_W, upper_curve_1tau_coh_W, color = color_1tau, alpha=0.50)
ax62.fill_between(energy_1tau_incoh_W, lower_curve_1tau_incoh_W, upper_curve_1tau_incoh_W, color = color_1tau, alpha=0.50) # due to large y-range, error seems smaller but it's not

# Iron #
ax63.errorbar(energy_1tau_coh_Fe, xsec_1tau_coh_Fe, color = color_1tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax63.errorbar(energy_1tau_incoh_Fe, xsec_1tau_incoh_Fe, color = color_1tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#txt63_1tau = ax63.text(20,1e-51,r'\textbf{$\nu_\mu$ Fe $\to \nu_\tau \tau^+ \mu^-$ Fe}',color=color_1tau,fontsize=25)
#txt63_1tau.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

ax63.fill_between(energy_1tau_coh_Fe, lower_curve_1tau_coh_Fe, upper_curve_1tau_coh_Fe, color = color_1tau, alpha=0.50)
ax63.fill_between(energy_1tau_incoh_Fe, lower_curve_1tau_incoh_Fe, upper_curve_1tau_incoh_Fe, color = color_1tau, alpha=0.50) # due to large y-range, error seems smaller but it's not

### vmu -> vmu tau+ tau- ###
# Argon #
ax2.errorbar(energy_2tau_coh_Ar, xsec_2tau_coh_Ar, color = color_2tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.errorbar(energy_2tau_incoh_Ar, xsec_2tau_incoh_Ar, color = color_2tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2.text(22,1e-52,r'$\tau^+ \tau^-$ {\bf Coh.}',color=color_2tau,rotation=53,fontsize=30)
#ax2.text(35,1e-46,r'$\tau^+ \tau^-$ {\bf Incoh.}',color=color_2tau,rotation=19,fontsize=30)
txt2_2tau = ax2.text(2,1e-41,r'\textbf{$\nu_\mu$ Ar $\to \nu_\mu \tau^+ \tau^-$ Ar}',color=color_2tau,fontsize=25)
txt2_2tau.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

upper_curve_2tau_coh_Ar = [xsec + sigma for xsec, sigma in zip(xsec_2tau_coh_Ar, delta_2tau_coh_Ar)]
lower_curve_2tau_coh_Ar = [xsec - sigma for xsec, sigma in zip(xsec_2tau_coh_Ar, delta_2tau_coh_Ar)]

upper_curve_2tau_incoh_Ar = [xsec + sigma for xsec, sigma in zip(xsec_2tau_incoh_Ar, delta_2tau_incoh_Ar)]
lower_curve_2tau_incoh_Ar = [xsec - sigma for xsec, sigma in zip(xsec_2tau_incoh_Ar, delta_2tau_incoh_Ar)]

ax2.fill_between(energy_2tau_coh_Ar, lower_curve_2tau_coh_Ar, upper_curve_2tau_coh_Ar, color = color_2tau, alpha=0.50)
ax2.fill_between(energy_2tau_incoh_Ar, lower_curve_2tau_incoh_Ar, upper_curve_2tau_incoh_Ar, color = color_2tau, alpha=0.50)

# Tungsten #
ax4.errorbar(energy_2tau_coh_W, xsec_2tau_coh_W, color = color_2tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax4.errorbar(energy_2tau_incoh_W, xsec_2tau_incoh_W, color = color_2tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax4.text(40,1e-49,r'$\tau^+ \tau^-$ {\bf Coh.}',color=color_2tau,rotation=38,fontsize=30)
#ax4.text(30,2e-46,r'$\tau^+ \tau^-$ {\bf Incoh.}',color=color_2tau,rotation=18,fontsize=30)
txt4_2tau = ax4.text(300,1e-52,r'\textbf{$\nu_\mu$ W $\to \nu_\mu \tau^+ \tau^-$ W}',color=color_2tau,fontsize=25)
txt4_2tau.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

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

## All ##
# Argon #
ax61.errorbar(energy_2tau_coh_Ar, xsec_2tau_coh_Ar, color = color_2tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax61.errorbar(energy_2tau_incoh_Ar, xsec_2tau_incoh_Ar, color = color_2tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
txt61_2tau = ax61.text(3,1e-41,r'\textbf{$\nu_\mu$ $\to \nu_\mu \tau^+ \tau^-$}',color=color_2tau,fontsize=50)
txt61_2tau.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

ax61.fill_between(energy_2tau_coh_Ar, lower_curve_2tau_coh_Ar, upper_curve_2tau_coh_Ar, color = color_2tau, alpha=0.50)
ax61.fill_between(energy_2tau_incoh_Ar, lower_curve_2tau_incoh_Ar, upper_curve_2tau_incoh_Ar, color = color_2tau, alpha=0.50)

# Tungsten #
ax62.errorbar(energy_2tau_coh_W, xsec_2tau_coh_W, color = color_2tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax62.errorbar(energy_2tau_incoh_W, xsec_2tau_incoh_W, color = color_2tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#txt4_2tau = ax4.text(300,1e-52,r'\textbf{$\nu_\mu$ W $\to \nu_\mu \tau^+ \tau^-$ W}',color=color_2tau,fontsize=25)
#txt4_2tau.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

ax62.fill_between(energy_2tau_coh_W, lower_curve_2tau_coh_W, upper_curve_2tau_coh_W, color = color_2tau, alpha=0.50)
ax62.fill_between(energy_2tau_incoh_W, lower_curve_2tau_incoh_W, upper_curve_2tau_incoh_W, color = color_2tau, alpha=0.50)

# Iron #
ax63.errorbar(energy_2tau_coh_Fe, xsec_2tau_coh_Fe, color = color_2tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax63.errorbar(energy_2tau_incoh_Fe, xsec_2tau_incoh_Fe, color = color_2tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#txt63_2tau = ax63.text(20,1e-51,r'\textbf{$\nu_\mu$ Fe $\to \nu_\tau \tau^+ \mu^-$ Fe}',color=color_2tau,fontsize=25)
#txt63_2tau.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

upper_curve_2tau_coh_Fe = [xsec + sigma for xsec, sigma in zip(xsec_2tau_coh_Fe, delta_2tau_coh_Fe)]
lower_curve_2tau_coh_Fe = [xsec - sigma for xsec, sigma in zip(xsec_2tau_coh_Fe, delta_2tau_coh_Fe)]

upper_curve_2tau_incoh_Fe = [xsec + sigma for xsec, sigma in zip(xsec_2tau_incoh_Fe, delta_2tau_incoh_Fe)]
lower_curve_2tau_incoh_Fe = [xsec - sigma for xsec, sigma in zip(xsec_2tau_incoh_Fe, delta_2tau_incoh_Fe)]

ax63.fill_between(energy_2tau_coh_Fe, lower_curve_2tau_coh_Fe, upper_curve_2tau_coh_Fe, color = color_2tau, alpha=0.50)
ax63.fill_between(energy_2tau_incoh_Fe, lower_curve_2tau_incoh_Fe, upper_curve_2tau_incoh_Fe, color = color_2tau, alpha=0.50) # due to large y-range, error seems smaller but it's not

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

ax2.errorbar(energy_2mu_coh_Ar, xsec_2mu_coh_Ar, color = color_2mu, label='_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.errorbar(energy_2mu_n_Ar, xsec_2mu_incoh_Ar, color = color_2mu, label='_hidden', linestyle = 'dashed', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2.text(50,2.1e-42,r'$\mu^+ \mu^-$ {\bf Coh.}',color=color_2mu,rotation=6,fontsize=30)
#ax2.text(50,6.5e-43,r'$\mu^+ \mu^-$ {\bf Incoh.}',color=color_2mu,rotation=5,fontsize=30)
txt2_2mu = ax2.text(2,1e-38,r'\textbf{$\nu_\mu$ Ar $\to \nu_\mu \mu^+ \mu^-$ Ar}',color=color_2mu,fontsize=25)
txt2_2mu.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

ax2.fill_between(energy_2mu_coh_Ar, lower_curve_2mu_coh_Ar, upper_curve_2mu_coh_Ar, color = color_2mu, alpha=0.50)
ax2.fill_between(energy_2mu_incoh_Ar, lower_curve_2mu_incoh_Ar, upper_curve_2mu_incoh_Ar, color = color_2mu, alpha=0.50)

# Tungsten #
ax4.errorbar(energy_2mu_coh_W, xsec_2mu_coh_W, color = color_2mu, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax4.errorbar(energy_2mu_incoh_W, xsec_2mu_incoh_W, color = color_2mu, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax4.text(5,7e-42,r'$\mu^+ \mu^-$ {\bf Coh.}',color=color_2mu,rotation=10,fontsize=30)
#ax4.text(2,1.5e-44,r'$\mu^+ \mu^-$ {\bf Incoh.}',color=color_2mu,rotation=10,fontsize=30)
txt4_2mu = ax4.text(300,1e-50,r'\textbf{$\nu_\mu$ W $\to \nu_\mu \mu^+ \mu^-$ W}',color=color_2mu,fontsize=25)
txt4_2mu.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

upper_curve_2mu_coh_W = [xsec + sigma for xsec, sigma in zip(xsec_2mu_coh_W, delta_2mu_coh_W)]
lower_curve_2mu_coh_W = [xsec - sigma for xsec, sigma in zip(xsec_2mu_coh_W, delta_2mu_coh_W)]

upper_curve_2mu_incoh_W = [xsec + sigma for xsec, sigma in zip(xsec_2mu_incoh_W, delta_2mu_incoh_W)]
lower_curve_2mu_incoh_W = [xsec - sigma for xsec, sigma in zip(xsec_2mu_incoh_W, delta_2mu_incoh_W)]

ax4.fill_between(energy_2mu_coh_W, lower_curve_2mu_coh_W, upper_curve_2mu_coh_W, color = color_2mu, alpha=0.50)
ax4.fill_between(energy_2mu_incoh_W, lower_curve_2mu_incoh_W, upper_curve_2mu_incoh_W, color = color_2mu, alpha=0.50)

# Iron #
ax5.errorbar(energy_2mu_coh_Fe, xsec_2mu_coh_Fe, color = color_2mu, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax5.errorbar(energy_2mu_incoh_Fe, xsec_2mu_incoh_Fe, color = color_2mu, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax5.text(5,7e-42,r'$\mu^+ \mu^-$ {\bf Coh.}',color=color_2mu,rotation=10,fontsize=30)
#ax5.text(2,1.5e-44,r'$\mu^+ \mu^-$ {\bf Incoh.}',color=color_2mu,rotation=10,fontsize=30)
txt5_2mu = ax5.text(20,1e-50,r'\textbf{$\nu_\mu$ Fe $\to \nu_\mu \mu^+ \mu^-$ Fe}',color=color_2mu,fontsize=25)
txt5_2mu.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

upper_curve_2mu_coh_Fe = [xsec + sigma for xsec, sigma in zip(xsec_2mu_coh_Fe, delta_2mu_coh_Fe)]
lower_curve_2mu_coh_Fe = [xsec - sigma for xsec, sigma in zip(xsec_2mu_coh_Fe, delta_2mu_coh_Fe)]

upper_curve_2mu_incoh_Fe = [xsec + sigma for xsec, sigma in zip(xsec_2mu_incoh_Fe, delta_2mu_incoh_Fe)]
lower_curve_2mu_incoh_Fe = [xsec - sigma for xsec, sigma in zip(xsec_2mu_incoh_Fe, delta_2mu_incoh_Fe)]

ax5.fill_between(energy_2mu_coh_Fe, lower_curve_2mu_coh_Fe, upper_curve_2mu_coh_Fe, color = color_2mu, alpha=0.50)
ax5.fill_between(energy_2mu_incoh_Fe, lower_curve_2mu_incoh_Fe, upper_curve_2mu_incoh_Fe, color = color_2mu, alpha=0.50)

# Nucleon #
ax3.errorbar(energy_2mu_p, xsec_2mu_p, color = color_2mu, linestyle = 'solid', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.errorbar(energy_2mu_n, xsec_2mu_n, color = color_2mu, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.text(0.8,9e-46,r'$\mu^+ \mu^-$ {\bf p}',color=color_2mu,rotation=36,fontsize=30)
ax3.text(1.7,1e-46,r'$\mu^+ \mu^-$ {\bf n}',color=color_2mu,rotation=28,fontsize=30)

upper_curve_2mu_p = [xsec + sigma for xsec, sigma in zip(xsec_2mu_p, delta_2mu_p)]
lower_curve_2mu_p = [xsec - sigma for xsec, sigma in zip(xsec_2mu_p, delta_2mu_p)]

upper_curve_2mu_n = [xsec + sigma for xsec, sigma in zip(xsec_2mu_n, delta_2mu_n)]
lower_curve_2mu_n = [xsec - sigma for xsec, sigma in zip(xsec_2mu_n, delta_2mu_n)]

ax3.fill_between(energy_2mu_p, lower_curve_2mu_p, upper_curve_2mu_p, color = color_2mu, alpha=0.50)
ax3.fill_between(energy_2mu_n, lower_curve_2mu_n, upper_curve_2mu_n, color = color_2mu, alpha=0.50)

## All ##
# Argon #
ax61.errorbar(energy_2mu_coh_Ar, xsec_2mu_coh_Ar, color = color_2mu, label='_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax61.errorbar(energy_2mu_n_Ar, xsec_2mu_incoh_Ar, color = color_2mu, label='_hidden', linestyle = 'dashed', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
txt61_2mu = ax61.text(3,1e-38,r'\textbf{$\nu_\mu$ $\to \nu_\mu \mu^+ \mu^-$}',color=color_2mu,fontsize=50)
txt61_2mu.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

ax61.fill_between(energy_2mu_coh_Ar, lower_curve_2mu_coh_Ar, upper_curve_2mu_coh_Ar, color = color_2mu, alpha=0.50)
ax61.fill_between(energy_2mu_incoh_Ar, lower_curve_2mu_incoh_Ar, upper_curve_2mu_incoh_Ar, color = color_2mu, alpha=0.50)

# Tungsten #
ax62.errorbar(energy_2mu_coh_W, xsec_2mu_coh_W, color = color_2mu, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax62.errorbar(energy_2mu_incoh_W, xsec_2mu_incoh_W, color = color_2mu, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#txt62_2mu = ax62.text(300,1e-50,r'\textbf{$\nu_\mu$ W $\to \nu_\mu \mu^+ \mu^-$ W}',color=color_2mu,fontsize=25)
#txt62_2mu.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

ax62.fill_between(energy_2mu_coh_W, lower_curve_2mu_coh_W, upper_curve_2mu_coh_W, color = color_2mu, alpha=0.50)
ax62.fill_between(energy_2mu_incoh_W, lower_curve_2mu_incoh_W, upper_curve_2mu_incoh_W, color = color_2mu, alpha=0.50)

# Iron #
ax63.errorbar(energy_2mu_coh_Fe, xsec_2mu_coh_Fe, color = color_2mu, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax63.errorbar(energy_2mu_incoh_Fe, xsec_2mu_incoh_Fe, color = color_2mu, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#txt5_2mu = ax5.text(20,1e-50,r'\textbf{$\nu_\mu$ Fe $\to \nu_\mu \mu^+ \mu^-$ Fe}',color=color_2mu,fontsize=25)
#txt5_2mu.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

ax63.fill_between(energy_2mu_coh_Fe, lower_curve_2mu_coh_Fe, upper_curve_2mu_coh_Fe, color = color_2mu, alpha=0.50)
ax63.fill_between(energy_2mu_incoh_Fe, lower_curve_2mu_incoh_Fe, upper_curve_2mu_incoh_Fe, color = color_2mu, alpha=0.50)

### ve -> vtau tau+ e- ###
# Argon #
ax2.errorbar(energy_ve1tau_coh_Ar, xsec_ve1tau_coh_Ar, color = color_ve1tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.errorbar(energy_ve1tau_incoh_Ar, xsec_ve1tau_incoh_Ar, color = color_ve1tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#ax2.text(5,1e-52,r'$\tau^+ e^-$ {\bf Coh.}',color=color_ve1tau,rotation=62,fontsize=30)
#ax2.text(18,8e-44,r'$\tau^+ e^-$ {\bf Incoh.}',color=color_ve1tau,rotation=10,fontsize=30)
#ax2.text(5.5,1e-46,r'$\tau^+ \mu^-$ {\bf Incoh.}',color=color_ve1tau,rotation=50,fontsize=30)
txt2_ve1tau = ax2.text(2,1e-40,r'\textbf{$\nu_e$ Ar $\to \nu_\tau \tau^+ e^-$ Ar}',color=color_ve1tau,fontsize=25)
txt2_ve1tau.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

upper_curve_ve1tau_coh_Ar = [xsec + sigma for xsec, sigma in zip(xsec_ve1tau_coh_Ar, delta_ve1tau_coh_Ar)]
lower_curve_ve1tau_coh_Ar = [xsec - sigma for xsec, sigma in zip(xsec_ve1tau_coh_Ar, delta_ve1tau_coh_Ar)]

upper_curve_ve1tau_incoh_Ar = [xsec + sigma for xsec, sigma in zip(xsec_ve1tau_incoh_Ar, delta_ve1tau_incoh_Ar)]
lower_curve_ve1tau_incoh_Ar = [xsec - sigma for xsec, sigma in zip(xsec_ve1tau_incoh_Ar, delta_ve1tau_incoh_Ar)]

ax2.fill_between(energy_ve1tau_coh_Ar, lower_curve_ve1tau_coh_Ar, upper_curve_ve1tau_coh_Ar, color = color_ve1tau, alpha=0.50)
ax2.fill_between(energy_ve1tau_incoh_Ar, lower_curve_ve1tau_incoh_Ar, upper_curve_ve1tau_incoh_Ar, color = color_ve1tau, alpha=0.50) # due to large y-range, error seems smaller but it's not

# Nucleon #
ax3.errorbar(energy_ve1tau_p, xsec_ve1tau_p, color = color_ve1tau, linestyle = 'solid', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.errorbar(energy_ve1tau_n, xsec_ve1tau_n, color = color_ve1tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax3.text(290,1.2e-44,r'$\tau^+ e^-$ {\bf p}',color=color_ve1tau,rotation=10,fontsize=30)
ax3.text(18,1e-46,r'$\tau^+ e^-$ {\bf n}',color=color_ve1tau,rotation=33,fontsize=30)

upper_curve_ve1tau_p = [xsec + sigma for xsec, sigma in zip(xsec_ve1tau_p, delta_ve1tau_p)]
lower_curve_ve1tau_p = [xsec - sigma for xsec, sigma in zip(xsec_ve1tau_p, delta_ve1tau_p)]

upper_curve_ve1tau_n = [xsec + sigma for xsec, sigma in zip(xsec_ve1tau_n, delta_ve1tau_n)]
lower_curve_ve1tau_n = [xsec - sigma for xsec, sigma in zip(xsec_ve1tau_n, delta_ve1tau_n)]

ax3.fill_between(energy_ve1tau_p, lower_curve_ve1tau_p, upper_curve_ve1tau_p, color = color_ve1tau, alpha=0.50)
ax3.fill_between(energy_ve1tau_n, lower_curve_ve1tau_n, upper_curve_ve1tau_n, color = color_ve1tau, alpha=0.50)

## All ##
# Argon #
ax61.errorbar(energy_ve1tau_coh_Ar, xsec_ve1tau_coh_Ar, color = color_ve1tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax61.errorbar(energy_ve1tau_incoh_Ar, xsec_ve1tau_incoh_Ar, color = color_ve1tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
txt61_ve1tau = ax61.text(3,1e-40,r'\textbf{$\nu_e$ $\to \nu_\tau \tau^+ e^-$}',color=color_ve1tau,fontsize=50)
txt61_ve1tau.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

ax61.fill_between(energy_ve1tau_coh_Ar, lower_curve_ve1tau_coh_Ar, upper_curve_ve1tau_coh_Ar, color = color_ve1tau, alpha=0.50)
ax61.fill_between(energy_ve1tau_incoh_Ar, lower_curve_ve1tau_incoh_Ar, upper_curve_ve1tau_incoh_Ar, color = color_ve1tau, alpha=0.50) # due to large y-range, error seems smaller but it's not

# Tungsten #
ax62.errorbar(energy_ve1tau_coh_W, xsec_ve1tau_coh_W, color = color_ve1tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax62.errorbar(energy_ve1tau_incoh_W, xsec_ve1tau_incoh_W, color = color_ve1tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#txt62_ve1tau = ax62.text(2,1e-40,r'\textbf{$\nu_e$ $\to \nu_\tau \tau^+ e^-$}',color=color_ve1tau,fontsize=50)
#txt62_ve1tau.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

upper_curve_ve1tau_coh_W = [xsec + sigma for xsec, sigma in zip(xsec_ve1tau_coh_W, delta_ve1tau_coh_W)]
lower_curve_ve1tau_coh_W = [xsec - sigma for xsec, sigma in zip(xsec_ve1tau_coh_W, delta_ve1tau_coh_W)]

upper_curve_ve1tau_incoh_W = [xsec + sigma for xsec, sigma in zip(xsec_ve1tau_incoh_W, delta_ve1tau_incoh_W)]
lower_curve_ve1tau_incoh_W = [xsec - sigma for xsec, sigma in zip(xsec_ve1tau_incoh_W, delta_ve1tau_incoh_W)]

ax62.fill_between(energy_ve1tau_coh_W, lower_curve_ve1tau_coh_W, upper_curve_ve1tau_coh_W, color = color_ve1tau, alpha=0.50)
ax62.fill_between(energy_ve1tau_incoh_W, lower_curve_ve1tau_incoh_W, upper_curve_ve1tau_incoh_W, color = color_ve1tau, alpha=0.50) 

# Iron #
ax63.errorbar(energy_ve1tau_coh_Fe, xsec_ve1tau_coh_Fe, color = color_ve1tau, label = "_hidden", path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax63.errorbar(energy_ve1tau_incoh_Fe, xsec_ve1tau_incoh_Fe, color = color_ve1tau, linestyle = 'dashed', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#txt62_ve1tau = ax62.text(2,1e-40,r'\textbf{$\nu_e$ $\to \nu_\tau \tau^+ e^-$}',color=color_ve1tau,fontsize=50)
#txt62_ve1tau.set_path_effects([pe.withStroke(linewidth=1, foreground='k')])

upper_curve_ve1tau_coh_Fe = [xsec + sigma for xsec, sigma in zip(xsec_ve1tau_coh_Fe, delta_ve1tau_coh_Fe)]
lower_curve_ve1tau_coh_Fe = [xsec - sigma for xsec, sigma in zip(xsec_ve1tau_coh_Fe, delta_ve1tau_coh_Fe)]

upper_curve_ve1tau_incoh_Fe = [xsec + sigma for xsec, sigma in zip(xsec_ve1tau_incoh_Fe, delta_ve1tau_incoh_Fe)]
lower_curve_ve1tau_incoh_Fe = [xsec - sigma for xsec, sigma in zip(xsec_ve1tau_incoh_Fe, delta_ve1tau_incoh_Fe)]

ax63.fill_between(energy_ve1tau_coh_Fe, lower_curve_ve1tau_coh_Fe, upper_curve_ve1tau_coh_Fe, color = color_ve1tau, alpha=0.50)
ax63.fill_between(energy_ve1tau_incoh_Fe, lower_curve_ve1tau_incoh_Fe, upper_curve_ve1tau_incoh_Fe, color = color_ve1tau, alpha=0.50) 

### Digitized ###
ax1.scatter(energy_2mu_coh_Ar_Alt, xsec_2mu_coh_Ar_Alt, marker='o', color = 'gray', s=40, edgecolors='k', linewidths=1, alpha=0.75, zorder=3, label='Altmannshofer+ 2019') 
ax1.scatter(energy_2mu_coh_Ar_Ballett, xsec_2mu_coh_Ar_Ballett, marker='^', color = 'black', s=40, edgecolors='k', linewidths=1, alpha=0.75, zorder=4, label='Ballett+ 2018')
ax1.scatter(energy_2mu_coh_Ar_Magill, xsec_2mu_coh_Ar_Magill, marker='v', color = 'silver', s=40, edgecolors='k', linewidths=1, alpha=0.75, zorder=5, label=r'Magill $\&$ Plestid 2017')

ax3.scatter(energy_2mu_incoh_p_Ar_Alt, np.divide(xsec_2mu_incoh_p_Ar_Alt,Z_Ar), marker='^', color = 'black', s=40, edgecolors='k', linewidths=1, alpha=0.75, zorder=4, label='Altmannshofer+ 2019')
ax3.scatter(energy_2mu_incoh_n_Ar_Alt, np.divide(xsec_2mu_incoh_n_Ar_Alt,(A_Ar-Z_Ar)), marker='v', color = 'silver', s=40, edgecolors='k', linewidths=1, alpha=0.75, zorder=4, label='Altmannshofer+ 2019')

#ax2.scatter(energy_2mu_coh_Ar_Alt, xsec_2mu_coh_Ar_Alt, marker='o', color = 'grey', s=40, edgecolors='k', linewidths=1, alpha=0.75, zorder=3, label='Altmannshofer+ 2019') 
#ax2.scatter(energy_2mu_coh_Ar_Ballett, xsec_2mu_coh_Ar_Ballett, marker='v', color = 'black', s=40, edgecolors='k', linewidths=1, alpha=0.75, zorder=4, label='Ballett+ 2018')

#ax2.scatter(energy_2tau_incoh_Ar_Beacom, xsec_2tau_incoh_Ar_Beacom, marker='^', color = 'black', s=40, edgecolors='k', linewidths=1, alpha=0.75, zorder=4, label=r'Zhou $\&$ Beacom 2020')
ax2.plot(energy_numuCC_FZ, xsec_numuCC_FZ ,'-', color = 'grey', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.text(100,4e-38,r'$\nu_\mu$ {\bf CC}',color='black',rotation=0,fontsize=30)

ax61.plot(energy_numuCC_FZ, xsec_numuCC_FZ ,'-', color = 'grey', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax61.text(88,4e-38,r'$\nu_\mu$ {\bf CC}',color='black',rotation=0,fontsize=45)
ax62.plot(energy_numuCC_FZ, xsec_numuCC_FZ ,'-', color = 'grey', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax62.text(88,4e-38,r'$\nu_\mu$ {\bf CC}',color='black',rotation=0,fontsize=45)
ax63.plot(energy_numuCC_FZ, xsec_numuCC_FZ ,'-', color = 'grey', label = '_hidden', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax63.text(88,4e-38,r'$\nu_\mu$ {\bf CC}',color='black',rotation=0,fontsize=45)

### Plotting options ###
ax1.set_xlabel(r'{\bf Neutrino Energy }$E_\nu$ (GeV)') 
ax1.set_ylabel(r'{\bf Cross Section / }$E_\nu$ (cm$^{2}$ GeV$^{-1}$)') 
ax1.set_ylim(1e-47, 1e-40)
ax1.set_xlim(0.1,220)
#ax1.legend(loc='lower right', fontsize=25)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_title(r'{\bf Argon $^{40}$Ar}', fontsize=40)

ax2.set_xlabel(r'{\bf Neutrino Energy }$E_\nu$ (GeV)') 
ax2.set_ylabel(r'{\bf Cross Section / }$E_\nu$ (cm$^{2}$ GeV$^{-1}$)') 
ax2.set_ylim(1e-56, 1e-36)
#ax2.set_ylim(1e-51, 1e-42)
ax2.set_xlim(0.1,220)
#ax2.legend(loc='lower right', fontsize=25)
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
ax4.set_ylim(1e-56, 1e-39)
ax4.set_xlim(0.1,11000)
#ax4.legend(loc='lower right') 
ax4.set_yscale('log')
ax4.set_xscale('log')
ax4.set_title(r'{\bf Tungsten $^{184}$W}', fontsize=40)

ax5.set_xlabel(r'{\bf Neutrino Energy }$E_\nu$ (GeV)') 
ax5.set_ylabel(r'{\bf Cross Section / }$E_\nu$ (cm$^{2}$ GeV$^{-1}$)') 
ax5.set_ylim(1e-54, 1e-40)
ax5.set_xlim(0.1,220)
#ax4.legend(loc='lower right') 
ax5.set_yscale('log')
ax5.set_xscale('log')
ax5.set_title(r'{\bf Iron $^{56}$Fe}', fontsize=40)

#ax61.set_xlabel(r'{\bf Neutrino Energy }$E_\nu$ (GeV)', fontsize=50)
#ax62.set_xlabel(r'{\bf Neutrino Energy }$E_\nu$ (GeV)', fontsize=50) 
#ax63.set_xlabel(r'{\bf Neutrino Energy }$E_\nu$ (GeV)', fontsize=50) 
fig6.supxlabel(r'\textbf{Neutrino Energy }$E_\nu$ (GeV)', fontsize=75)
ax61.set_ylabel(r'{\bf Cross Section / }$E_\nu$ (cm$^{2}$ GeV$^{-1}$)', fontsize=75) 
ax61.set_ylim(1e-56, 1e-36)
ax61.set_xlim(0.1,220)
ax62.set_xlim(0.1,10000)
ax63.set_xlim(0.1,220)
ax61.set_yscale('log')
ax61.set_xscale('log')
ax62.set_xscale('log')
ax63.set_xscale('log')
ax61.set_title(r'{\bf Argon $^{40}$Ar}', fontsize=60)
ax62.set_title(r'{\bf Tungsten $^{184}$W}', fontsize=60)
ax63.set_title(r'{\bf Iron $^{56}$Fe}', fontsize=60)

locmaj11 = mticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
locmin11 = mticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)
locmaj21 = mticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
locmin21 = mticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)
locmaj31 = mticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
locmin31 = mticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)
locmaj41 = mticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
locmin41 = mticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)
locmaj51 = mticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
locmin51 = mticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)

#locmaj61 = mticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
#locmin61 = mticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)

ax1.xaxis.set_major_locator(locmaj11)
ax1.xaxis.set_minor_locator(locmin11)
ax2.xaxis.set_major_locator(locmaj21)
ax2.xaxis.set_minor_locator(locmin21)
ax3.xaxis.set_major_locator(locmaj31)
ax3.xaxis.set_minor_locator(locmin31)
ax4.xaxis.set_major_locator(locmaj41)
ax4.xaxis.set_minor_locator(locmin41)
ax5.xaxis.set_major_locator(locmaj51)
ax5.xaxis.set_minor_locator(locmin51)

locmaj12 = mticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
locmin12 = mticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)
locmaj22 = mticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
locmin22 = mticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)
locmaj32 = mticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
locmin32 = mticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)
locmaj42 = mticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
locmin42 = mticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)
locmaj52 = mticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
locmin52 = mticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)

ax1.yaxis.set_major_locator(locmaj12)
ax1.yaxis.set_minor_locator(locmin12)
ax2.yaxis.set_major_locator(locmaj22)
ax2.yaxis.set_minor_locator(locmin22)
ax3.yaxis.set_major_locator(locmaj32)
ax3.yaxis.set_minor_locator(locmin32)
ax4.yaxis.set_major_locator(locmaj42)
ax4.yaxis.set_minor_locator(locmin42)
ax5.yaxis.set_major_locator(locmaj52)
ax5.yaxis.set_minor_locator(locmin52)

ax1.tick_params(which='both',right=False, bottom=True)
ax2.tick_params(which='both',right=False, bottom=True)
ax3.tick_params(which='both',right=False, bottom=True)
ax4.tick_params(which='both',right=False, bottom=True)
ax5.tick_params(which='both',right=False, bottom=True)
ax61.tick_params(which='both',right=False, bottom=True)
ax62.tick_params(which='both',right=False, bottom=True)
ax63.tick_params(which='both',right=False, bottom=True)

xmajor1 = [0.1, 1, 5, 10, 50, 100]
xmajor2 = [0.1, 1, 5, 10, 50, 100]
xmajor3 = [0.1, 1, 5, 10, 50, 100, 1000, 10000]
xmajor4 = [0.1, 1, 5, 10, 50, 100, 1000, 10000]
xmajor5 = [0.1, 1, 5, 10, 50, 100]
xmajor62 = [0.1, 1, 5, 10, 50, 100, 1000]

ax1.set_xticks(xmajor1, labels=['0.1', '1', '5', '10', '50', '100'])
ax1.grid(which='major', axis='both')
ax2.set_xticks(xmajor2, labels=['0.1', '1', '5', '10', '50', '100'])
ax2.grid(which='major', axis='both')
ax3.set_xticks(xmajor3, labels=['0.1', '1', '5', '10', '50', '100', '1000', '10000'])
ax3.grid(which='major', axis='both')
ax4.set_xticks(xmajor4, labels=['0.1', '1', '5', '10', '50', '100', '1000', '10000'])
ax4.grid(which='major', axis='both')
ax5.set_xticks(xmajor5, labels=['0.1', '1', '5', '10', '50', '100'])
ax5.grid(which='major', axis='both')

ax61.set_xticks(xmajor2, labels=['0.1', '1', '5', '10', '50', '100'])
ax61.grid(which='major', axis='both')
ax62.set_xticks(xmajor62, labels=['0.1', '1', '5', '10', '50', '100', '1000'])
ax62.grid(which='major', axis='both')
ax63.set_xticks(xmajor5, labels=['0.1', '1', '5', '10', '50', '100'])
ax63.grid(which='major', axis='both')
ax61.tick_params(axis='y',labelsize=50)
ax61.tick_params(axis='x',labelsize=50)
ax62.tick_params(axis='x',labelsize=50)
ax63.tick_params(axis='x',labelsize=50)

fig1.savefig("../plots/xsec_perE_mu_argon.png", dpi=100, bbox_inches='tight')
fig2.savefig("../plots/xsec_perE_argon.png", dpi=100, bbox_inches='tight')
fig3.savefig("../plots/xsec_perE_nucleons.png", dpi=100, bbox_inches='tight')
fig4.savefig("../plots/xsec_perE_tungsten.png", dpi=100, bbox_inches='tight')
fig5.savefig("../plots/xsec_perE_iron.png", dpi=100, bbox_inches='tight')
fig6.savefig("../plots/xsec_perE.png", dpi=100, bbox_inches='tight')
