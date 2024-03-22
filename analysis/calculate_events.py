import numpy as np
import csv
from scipy.integrate import simpson

atomic_mass_unit = 1.6605e-27 # kg

N_POT = 1.1e21               # Number of POTs
MD = 1000                    # Mass of the DUNE ND detector in kg (1 tonne)
MAr = 39.95*atomic_mass_unit # Mass of argon in atomic mass units

N_tonne = 147 # Should be 50 now for liquid argon.
N_year = 3

Phi_Alt = 1.04e-3

DUNE_filename = 'csv/DUNE_diff_flux.csv'
Alt_flux = 'csv/vmu_normalized_flux_Altmannshofer_digitized.csv'
#DUNE_FD_tau-opt_flux = 'csv/DUNE_tau-opt_FD_flux.csv'

xsec_1tau_filename = 'csv/vmu_to_vtau_tau+_mu-_xsec/vmu_to_vtau_tau+_mu-_xsec.csv'
xsec_2tau_filename = 'csv/vmu_to_vmu_tau+_tau-_xsec/coherent/argon/vmu_to_vmu_tau+_tau-_xsec.csv'
xsec_2tau_filename_p = 'csv/vmu_to_vmu_tau+_tau-_xsec/incoherent/proton/vmu_p_to_vmu_tau+_tau-_xsec.csv'
xsec_2mu_filename  = 'csv/vmu_to_vmu_mu+_mu-_xsec/vmu_to_vmu_mu+_mu-_xsec.csv'
xsec_cc_filename = 'csv/vmuCC_xsec_perE.csv'
#xsec_2mu_filename  = 'csv/vmu_to_vmu_mu+_mu-_xsec_digitized.csv'

energy_DUNE = []
diff_flux_DUNE = []

energy_Alt = []
norm_flux_Alt = []

#energy_DUNE_tau-opt = []
#diff_flux_DUNE_tau-opt = []

energy_1tau = []
xsec_1tau = []

energy_2tau = []
xsec_2tau = []

energy_2mu = []
xsec_2mu = []

xsec_1tau_matched = []
xsec_2tau_matched = []
xsec_2mu_matched = []

xsec_1tau_matched_Alt = []
xsec_2tau_matched_Alt = []
xsec_2mu_matched_Alt = []

with open(DUNE_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_DUNE.append(float(row[0]))
        diff_flux_DUNE.append(float(row[1]) * 1e4) # DUNE differential flux is in [cm^-2 GeV^-1 POT^-1]; convert to [m^-2 GeV^-1 POT^-1]

with open(Alt_flux,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_Alt.append(float(row[0]))
        norm_flux_Alt.append(float(row[1]))

#with open(DUNE_FD_tau-opt_flux,'r') as csvfile:
#    data = csv.reader(csvfile, delimiter = ',')
#    for row in data:
#        energy_DUNE_tau-opt.append(float(row[0]))
#        diff_flux_DUNE_tau-opt.append(float(row[1]) * 1e-21)

with open(xsec_1tau_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_1tau.append(float(row[0]))
        xsec_1tau.append(float(row[1]) * 1e-43) # Convert to m^2.

with open(xsec_2tau_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2tau.append(float(row[0]))
        xsec_2tau.append(float(row[1]) * 1e-43) # Convert to m^2.

with open(xsec_2mu_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_2mu.append(float(row[0]))
        xsec_2mu.append(float(row[1]) * 1e-43) # Convert to m^2.

with open(xsec_cc_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_cc.append(float(row[0]))

# Match xsec to energy given by DUNE flux.
for i in range(len(energy_DUNE)):
    for j in range(len(energy_1tau)):
        if (energy_DUNE[i] < energy_1tau[0]): # Check that DUNE energy val is within the energy range of 1tau trident energies
            xsec_1tau_matched.append(0.0) # and set xsec to 0; otherwise, these values will be skipped altogether
            break
        elif ((energy_DUNE[i] >= energy_1tau[j]) and (energy_DUNE[i] <= energy_1tau[j+1])): # If the DUNE energy val is between two energy vals in 1tau trident
            xsec_1tau_matched.append(xsec_1tau[j]) # Add the 1tau trident cross section of the lower of the two energy vals to the xsec array matched at the DUNE energy val.

    for j in range(len(energy_2tau)):
        if (energy_DUNE[i] < energy_2tau[0]):
            xsec_2tau_matched.append(0.0)
            break
        elif ((energy_DUNE[i] >= energy_2tau[j]) and (energy_DUNE[i] <= energy_2tau[j+1])):
            xsec_2tau_matched.append(xsec_2tau[j])

    for j in range(len(energy_2mu)):
        if (energy_DUNE[i] < energy_2mu[0]):
            xsec_2mu_matched.append(0.0)
            break
        elif ((energy_DUNE[i] >= energy_2mu[j]) and (energy_DUNE[i] <= energy_2mu[j+1])):
            xsec_2mu_matched.append(xsec_2mu[j])

for i in range(len(energy_Alt)):
    for j in range(len(energy_1tau)):
        if (energy_Alt[i] < energy_1tau[0]):
            xsec_1tau_matched_Alt.append(0.0)
            break
        elif ((energy_Alt[i] >= energy_1tau[j]) and (energy_Alt[i] <= energy_1tau[j+1])):
            xsec_1tau_matched_Alt.append(xsec_1tau[j])

    for j in range(len(energy_2tau)):
        if (energy_Alt[i] < energy_2tau[0]):
            xsec_2tau_matched_Alt.append(0.0)
            break
        elif ((energy_Alt[i] >= energy_2tau[j]) and (energy_Alt[i] <= energy_2tau[j+1])):
            xsec_2tau_matched_Alt.append(xsec_2tau[j])

    for j in range(len(energy_2mu)):
        if (energy_Alt[i] < energy_2mu[0]):
            xsec_2mu_matched_Alt.append(0.0)
            break
        elif ((energy_Alt[i] >= energy_2mu[j]) and (energy_Alt[i] <= energy_2mu[j+1])):
            xsec_2mu_matched_Alt.append(xsec_2mu[j])

#print("Size of DUNE energy array: ", len(energy_DUNE))
#print("Size of Altmannshofer energy array: ", len(energy_Alt))

#print("Size of cross section array for 1 tau: ", len(xsec_1tau_matched))
#print("Size of cross section array for 2 tau: ", len(xsec_2tau_matched))
#print("Size of cross section array for 2 mu: ", len(xsec_2mu_matched))

#print("Size of cross section array for 1 tau (Alt): ", len(xsec_1tau_matched_Alt))
#print("Size of cross section array for 2 tau (Alt): ", len(xsec_2tau_matched_Alt))
#print("Size of cross section array for 2 mu (Alt): ", len(xsec_2mu_matched_Alt))

# Calculate the expected number of events
Phi = simpson(diff_flux_DUNE, x=energy_DUNE)
print("Integrated neutrino flux DUNE: ", Phi)
print("Integrated neutrino flux Altmannshofer: ", Phi_Alt)

sigma_1tau = np.multiply(diff_flux_DUNE, xsec_1tau_matched) / Phi
sigma_1tau_Alt = np.multiply(norm_flux_Alt, xsec_1tau_matched_Alt)
#print("Size of convoluted sigma integrand for 1 tau: ", len(sigma_1tau))
#print("Size of convoluted sigma integrand for 1 tau (Alt): ", len(sigma_1tau_Alt))

Sigma_1tau = simpson(sigma_1tau, x=energy_DUNE)
Sigma_1tau_Alt = simpson(sigma_1tau_Alt, x=energy_Alt)

#print(sigma_1tau)
#print(Sigma_1tau)

N_1tau = Phi * Sigma_1tau * MD * N_POT / MAr
N_1tau_Alt = Phi_Alt * Sigma_1tau_Alt * MD * N_POT / MAr

print("Events expected at DUNE ND for 1tau trident:")
print("\t1 tonne Ar, 1 year, DUNE flux: ", N_1tau)
print("\t147 tonne Ar, 3 year, DUNE flux: ", N_1tau * N_tonne * N_year)
print("\t147 tonne Ar, 10 year, DUNE flux: ", N_1tau * N_tonne * 10)
print("\t1 tonne Ar, 1 year, Altmannshofer flux: ", N_1tau_Alt)
print("\t147 tonne Ar, 3 year, Altmannshofer flux: ", N_1tau_Alt * N_tonne * N_year)
print("\t147 tonne Ar, 10 year, Altmannshofer flux: ", N_1tau_Alt * N_tonne * 10)

sigma_2tau = np.multiply(diff_flux_DUNE, xsec_2tau_matched) / Phi
sigma_2tau_Alt = np.multiply(norm_flux_Alt, xsec_2tau_matched_Alt)
#print("Size of convoluted sigma integrand for 2 tau: ", len(sigma_2tau))
#print("Size of convoluted sigma integrand for 2 tau (Alt): ", len(sigma_2tau_Alt))

Sigma_2tau = simpson(sigma_2tau, x=energy_DUNE)
Sigma_2tau_Alt = simpson(sigma_2tau_Alt, x=energy_Alt)

#print(sigma_2tau)
#print(Sigma_2tau)

N_2tau = Phi * Sigma_2tau * MD * N_POT / MAr
N_2tau_Alt = Phi_Alt * Sigma_2tau_Alt * MD * N_POT / MAr

print("Events expected at DUNE ND for 2tau trident:")
print("\t1 tonne Ar, 1 year, DUNE flux: ", N_2tau)
print("\t147 tonne Ar, 3 year, DUNE flux: ", N_2tau * N_tonne * N_year)
print("\t147 tonne Ar, 10 year, DUNE flux: ", N_2tau * N_tonne * 10)
print("\t1 tonne Ar, 1 year, Altmannshofer flux: ", N_2tau_Alt)
print("\t147 tonne Ar, 3 year, Altmannshofer flux: ", N_2tau_Alt * N_tonne * N_year)
print("\t147 tonne Ar, 10 year, Altmannshofer flux: ", N_2tau_Alt * N_tonne * 10)


sigma_2mu = np.multiply(diff_flux_DUNE, xsec_2mu_matched) / Phi
sigma_2mu_Alt = np.multiply(norm_flux_Alt, xsec_2mu_matched_Alt)
#print("Size of convoluted sigma integrand for 2 mu: ", len(sigma_2mu))
#print("Size of convoluted sigma integrand for 2 mu (Alt): ", len(sigma_2mu_Alt))

Sigma_2mu = simpson(sigma_2mu, x=energy_DUNE)
Sigma_2mu_Alt = simpson(sigma_2mu_Alt, x=energy_Alt)

#print(sigma_2mu)
#print(Sigma_2mu)

N_2mu = Phi * Sigma_2mu * MD * N_POT / MAr
N_2mu_Alt = Phi_Alt * Sigma_2mu_Alt * MD * N_POT / MAr

print("Events expected at DUNE ND for 2mu trident:")
print("\t1 tonne Ar, 1 year, DUNE flux: ", N_2mu)
print("\t147 tonne Ar, 3 year, DUNE flux: ", N_2mu * N_tonne * N_year)
print("\t147 tonne Ar, 10 year, DUNE flux: ", N_2mu * N_tonne * 10)
print("\t1 tonne Ar, 1 year, Altmannshofer flux: ", N_2mu_Alt)
print("\t147 tonne Ar, 3 year, Altmannshofer flux: ", N_2mu_Alt * N_tonne * N_year)
print("\t147 tonne Ar, 10 year, Altmannshofer flux: ", N_2mu_Alt * N_tonne * 10)

