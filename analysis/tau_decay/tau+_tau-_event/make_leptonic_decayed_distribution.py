import numpy as np
import math
import csv
import pyhepmc

from boosted import lorentz_factor, boost_matrix, apply_boost

EVENTS_DIR_1TAU = '../csv/events/vmu_to_vtau_tau+_mu-'

############################
##### Helper Functions #####
############################

### tau+ decay from Pythia8 ###
def is_final_state(particle):
    """ Check if the particle is in the final state """
    return particle.status == 1

def contains_leptonic_decay(products):
    """ Check if the decay products contain any leptonic decay channels """
    leptonic_pairs = [{-11, 12}, {11, -12}, {-13, 14}, {13, -14}]
    product_pids = {p.pid for p in products}
    return any(pair.issubset(product_pids) for pair in leptonic_pairs)

def find_final_state_particles(event, tau_plus):
    """ Find final state particles originating from the tau+ decay """
    final_state_particles = []
    to_process = [tau_plus]

    while to_process:
        particle = to_process.pop()
        if is_final_state(particle):
            final_state_particles.append(particle)
        else:
            # Add decay products to the list to process
            for vertex in particle.end_vertex.particles_out:
                to_process.append(vertex)
    
    return final_state_particles

def process_hepmc_file(filename):
    all_final_state_particles = []

    with pyhepmc.open(filename) as f:
        for event in f:
            # Find the tau+ particle
            tau_plus = None
            for particle in event.particles:
                if particle.pid == -15:  # PDG ID for tau+
                    tau_plus = particle
                    break

            if tau_plus is None:
                continue
            
            if tau_plus.end_vertex:
                decay_products = tau_plus.end_vertex.particles_out
                if not contains_leptonic_decay(decay_products):
                    continue
           
            # Find final state particles from tau+ decay
            final_state_particles = find_final_state_particles(event, tau_plus)
            all_final_state_particles.append(final_state_particles)
        
    return all_final_state_particles

def read_tau_decays(all_final_state_particles, pdg_neutrinos, pdg_leptons):
    """ 
    Reads all the leptonic decay products of the tau+ and returns arrays for the coherent and incoherent regimes.
    For coherent scattering:
        Only stores vtau_bar and ve/vmu information. Returns array of vectorial sum of neutrino four vectors for each event.
    For incoherent scattering:
        Stores the four vector sum of the particle PDG codes contained in the pdg_to_store list.
        Returns array of the summed four vector for each event.
    Outputs a tuple of coherent and incoherent arrays.
    """
    coherent_array = []
    incoherent_array = []
    for i, final_state_particles in enumerate(all_final_state_particles):
        sum_vector_coh = np.array([0.0, 0.0, 0.0, 0.0])
        sum_vector_incoh = np.array([0.0, 0.0, 0.0, 0.0])
        for particle in final_state_particles:
            if particle.pid in pdg_neutrinos:
                four_vector = np.array([particle.momentum[3], particle.momentum[0], particle.momentum[1], particle.momentum[2]])
                sum_vector_coh = four_vector + sum_vector_coh
            if particle.pid in pdg_leptons:
                four_vector = np.array([particle.momentum[3], particle.momentum[0], particle.momentum[1], particle.momentum[2]])
                sum_vector_incoh = four_vector + sum_vector_incoh
        coherent_array.append(sum_vector_coh)
        incoherent_array.append(sum_vector_incoh)
    return coherent_array, incoherent_array

### events from TEG generator ###
def read_four_vectors(file_path, skip_lines=37, coherent=True):
    with open(file_path, 'r') as txtfile:
        for _ in range(skip_lines):
            next(txtfile)
        if coherent: # There are 100,000 events. For hadronic decays, we used the first 65,150. Let's skip those and use the remaining 34,850.
            for _ in range(6*65150): # Each coherent event takes 4+2 lines.
                next(txtfile)
        else:
            for _ in range(7*65150): # Each incoherent event takes 5+2 lines.
                next(txtfile)
        incoming_nu_mu_array = []
        primary_nu_tau_array = []
        primary_muon_array = []
        boosted_tau_array = []
        outgoing_nucleon_array = [] # Only incoherent proton and neutron will have an outgoing nucleon
        for line in txtfile:
            if line.startswith('14'):  # Read and store incoming nu_mu
                data = line.split()
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                E  = float(data[5])
                incoming_nu_mu = np.array([E, px, py, pz])
                incoming_nu_mu_array.append(incoming_nu_mu)
            if line.startswith('16'):  # Read and store primary tau neutrino
                data = line.split()
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                E  = float(data[5])
                primary_nu_tau = np.array([E, px, py, pz])
                primary_nu_tau_array.append(primary_nu_tau)
            if line.startswith('13'):  # Read and store primary muon
                data = line.split()
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                E  = float(data[5])
                primary_muon = np.array([E, px, py, pz])
                primary_muon_array.append(primary_muon)
            if line.startswith('-15'):  # Read and store boosted tau
                data = line.split()
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                E  = float(data[5])
                boosted_tau = np.array([E, px, py, pz])
                boosted_tau_array.append(boosted_tau)
            if line.startswith('2112') or line.startswith('2212'):  # Read and store outgoing nucleon
                data = line.split()
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                E  = float(data[5])
                outgoing_nucleon = np.array([E, px, py, pz])
                outgoing_nucleon_array.append(outgoing_nucleon)
        return incoming_nu_mu_array, primary_nu_tau_array, primary_muon_array, boosted_tau_array, outgoing_nucleon_array

### Boosting of tau+ decays to events frame ###
def decay_tau(boosted_tau, fv):
    """ Boost four vector fv to the frame of the boosted_tau and return boosted_fv"""
    L = boost_matrix(boosted_tau)
    boosted_fv = apply_boost(L, fv)
    return boosted_fv

def make_tau_decayed_four_vectors(boosted_tau_array, fv_array):
    """ Produce array of four vectors boosted_fv for each boosted_tau and four vector fv in boosted_tau_array and fv_array """
    boosted_fv_array = []
    for boosted_tau, fv in zip(boosted_tau_array, fv_array):
        boosted_fv = decay_tau(boosted_tau, fv)
        boosted_fv_array.append(boosted_fv)
    return boosted_fv_array

### Exporting of information ###
def printout_incoh(filename, incoming_nu_mu_array, primary_nu_tau_array, primary_muon_array, boosted_nubar_tau_array, outgoing_nucleon_array, nucleon):
    with open(filename,'w',newline='') as outfile:
        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        if nucleon == 'proton':
            nucleon_PID = '2212'
            nucleon_mass = 0.93827200
        if nucleon == 'neutron':
            nucleon_PID = '2112'
            nucleon_mass = 0.93956500
        for incoming_nu_mu, primary_nu_tau, primary_muon, boosted_nubar_tau, outgoing_nucleon in zip(incoming_nu_mu_array,
                                                                                                                     primary_nu_tau_array,
                                                                                                                     primary_muon_array,
                                                                                                                     boosted_nubar_tau_array,
                                                                                                                     outgoing_nucleon_array):
            writer.writerow(['<event>'])
            writer.writerow(['14', '-1', incoming_nu_mu[1], incoming_nu_mu[2], incoming_nu_mu[3], incoming_nu_mu[0], '0.0'])
            writer.writerow(['16', '1', primary_nu_tau[1], primary_nu_tau[2], primary_nu_tau[3], primary_nu_tau[0], '0.0'])
            writer.writerow(['13', '1', primary_muon[1], primary_muon[2], primary_muon[3], primary_muon[0], '0.10565800'])
            writer.writerow(['999', '2', boosted_nubar_tau[1], boosted_nubar_tau[2], boosted_nubar_tau[3], boosted_nubar_tau[0], '0.0'])
            writer.writerow([nucleon_PID, '1', outgoing_nucleon[1], outgoing_nucleon[2], outgoing_nucleon[3], outgoing_nucleon[0], nucleon_mass])
            writer.writerow(['</event>'])

def printout_coh(filename, incoming_nu_mu_array, primary_nu_tau_array, primary_muon_array, boosted_nubar_tau_array):
    with open(filename,'w',newline='') as outfile:
        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for incoming_nu_mu, primary_nu_tau, primary_muon, boosted_nubar_tau in zip(incoming_nu_mu_array,
                                                                                                   primary_nu_tau_array,
                                                                                                   primary_muon_array,
                                                                                                   boosted_nubar_tau_array):
            writer.writerow(['<event>'])
            writer.writerow(['14', '-1', incoming_nu_mu[1], incoming_nu_mu[2], incoming_nu_mu[3], incoming_nu_mu[0], '0.0'])
            writer.writerow(['16', '1', primary_nu_tau[1], primary_nu_tau[2], primary_nu_tau[3], primary_nu_tau[0], '0.0'])
            writer.writerow(['13', '1', primary_muon[1], primary_muon[2], primary_muon[3], primary_muon[0], '0.10565800'])
            writer.writerow(['999', '2', boosted_nubar_tau[1], boosted_nubar_tau[2], boosted_nubar_tau[3], boosted_nubar_tau[0], '0.0'])
            writer.writerow(['</event>'])

############################
###### Distributions #######
############################

part_dict = {211:'pi+', 
             -211:'pi-',
             311:'K0',
             321:'K+',
             -321:'K-',
             22:'gamma',
             310:'K0S',
             111:'pi0',
             -311:'K0BAR',
             16:'vtau',
             -16:'vtaubar',
             14:'vmu',
             -14:'vmubar',
             13:'mu-',
             -13:'mu+',
             12:'ve',
             -12:'vebar',
             11:'e-',
             -11:'e+',
             130:'K0L'}

TAU_DECAY_FILENAME = EVENTS_DIR_1TAU+'/../tau_decay/hadronic/tag_1_pythia8_events.hepmc'

final_state_particles_list = process_hepmc_file(TAU_DECAY_FILENAME)
pdg_leptons = [-11,   # e+
               -13,   # mu+
               ]

pdg_neutrinos = [12,   # ve
                 14,   # vmu
                 -16,  # vtau_bar
                ]

coherent_array, incoherent_array = read_tau_decays(final_state_particles_list, pdg_neutrinos, pdg_leptons)

print(len(coherent_array))
print(len(incoherent_array))

incoming_nu_mu_p_33GeV, primary_nu_tau_p_33GeV, primary_muon_p_33GeV, boosted_tau_p_33GeV, outgoing_nucleon_p_33GeV = read_four_vectors(EVENTS_DIR_1TAU+'/nucleon/proton/1tau_1e5Events_p_33GeV.txt', coherent=False)
incoming_nu_mu_n_33GeV, primary_nu_tau_n_33GeV, primary_muon_n_33GeV, boosted_tau_n_33GeV, outgoing_nucleon_n_33GeV = read_four_vectors(EVENTS_DIR_1TAU+'/nucleon/neutron/1tau_1e5Events_n_33GeV.txt', coherent=False)
incoming_nu_mu_Ar_47GeV, primary_nu_tau_Ar_47GeV, primary_muon_Ar_47GeV, boosted_tau_Ar_47GeV, outgoing_nucleon_Ar_47GeV = read_four_vectors(EVENTS_DIR_1TAU+'/coherent/argon/1tau_1e5Events_Ar_47GeV.txt')

print(len(incoming_nu_mu_p_33GeV))
print(len(incoming_nu_mu_n_33GeV))
print(len(incoming_nu_mu_Ar_47GeV))

boosted_nubar_tau_p_33GeV  = make_tau_decayed_four_vectors(boosted_tau_p_33GeV, incoherent_array)
boosted_nubar_tau_n_33GeV  = make_tau_decayed_four_vectors(boosted_tau_n_33GeV, incoherent_array)
boosted_nubar_tau_Ar_47GeV  = make_tau_decayed_four_vectors(boosted_tau_Ar_47GeV, coherent_array)

printout_incoh('../csv/distributions/vmu_to_vtau_tau+_mu-/nucleon/proton/tau_leptonic_decayed_distribution_p_33GeV.txt',
               incoming_nu_mu_p_33GeV, primary_nu_tau_p_33GeV, primary_muon_p_33GeV, boosted_nubar_tau_p_33GeV, outgoing_nucleon_p_33GeV,
               'proton')
printout_incoh('../csv/distributions/vmu_to_vtau_tau+_mu-/nucleon/neutron/tau_leptonic_decayed_distribution_n_33GeV.txt',
               incoming_nu_mu_n_33GeV, primary_nu_tau_n_33GeV, primary_muon_n_33GeV, boosted_nubar_tau_n_33GeV, outgoing_nucleon_n_33GeV,
               'neutron')
printout_coh('../csv/distributions/vmu_to_vtau_tau+_mu-/coherent/argon/tau_leptonic_decayed_distribution_Ar_47GeV.txt',
               incoming_nu_mu_Ar_47GeV, primary_nu_tau_Ar_47GeV, primary_muon_Ar_47GeV, boosted_nubar_tau_Ar_47GeV)
