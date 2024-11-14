import numpy as np
import math
import csv
import pyhepmc

from boosted import lorentz_factor, boost_matrix, apply_boost, boost_to_cm_matrix

EVENTS_DIR_1TAU = '../../../csv/events/vmu_to_vtau_tau+_mu-'

Ev = input("Neutrino energy (GeV): ")

write_dist = True

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

def contains_tau(products):
    """ Check if the decay of the tau+ contains a tau+ """
    tau_pid = -15
    product_pids = {p.pid for p in products}
    return tau_pid in product_pids

def pop_tau(products):
    """ When the tau+ decays into a tau+, pop the new tau+ """
    for p in products:
        if p.pid == -15:
            return p

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
    all_taus = []
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
                decay_products = tau_plus.end_vertex.particles_out  # This will just check the immediate decay of the tau+.
                while contains_tau(decay_products):  # Some tau+ will emit a photon before they decay. This will cause the decay product to be tau+ -> tau+, gamma.
                    new_tau = pop_tau(decay_products)  # We wish to know the decay type of the "new" tau+.
                    decay_products = new_tau.end_vertex.particles_out
                if contains_leptonic_decay(decay_products):
                    continue

            # Find final state particles from tau+ decay
            final_state_particles = find_final_state_particles(event, tau_plus)
            all_final_state_particles.append(final_state_particles)

            # Prepare tau+'s
            four_vector = np.array([tau_plus.momentum[3], tau_plus.momentum[0], tau_plus.momentum[1], tau_plus.momentum[2]])
            all_taus.append(four_vector)
            
    return all_final_state_particles, all_taus

def read_tau_decays(all_final_state_particles, pdg_to_store):
    """ 
    Reads all the decay products of the tau+ and returns arrays for the coherent and incoherent regimes.
    For coherent scattering:
        Only stores vtau_bar information. Returns array of vtau_bar four vectors for each event.
    For incoherent scattering:
        Stores the four vector sum of the particle PDG codes contained in the pdg_to_store list.
        Returns array of the summed four vector for each event.
    Outputs a tuple of coherent and incoherent arrays.
    """
    coherent_array = []
    incoherent_array = []
    for i, final_state_particles in enumerate(all_final_state_particles):
        sum_vector = np.array([0.0, 0.0, 0.0, 0.0])
        for particle in final_state_particles:
            if particle.pid == -16:
                four_vector = np.array([particle.momentum[3], particle.momentum[0], particle.momentum[1], particle.momentum[2]])
                coherent_array.append(four_vector)
            if particle.pid in pdg_to_store:
                four_vector = np.array([particle.momentum[3], particle.momentum[0], particle.momentum[1], particle.momentum[2]])
                sum_vector = four_vector + sum_vector
        incoherent_array.append(sum_vector)
    return coherent_array, incoherent_array

### events from TEG generator ###
def read_four_vectors(file_path, skip_lines=37):
    with open(file_path, 'r') as txtfile:
        for _ in range(skip_lines):
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
def make_boost_matrix(tau_at_decay_frame, tau_at_trident_frame):
    """ Prepare the boost matrix that will take a four vector from the tau+ decay frame to the tau+ trident frame"""
    L_decay_to_rest = boost_to_cm_matrix(tau_at_decay_frame)
    L_rest_to_trident = boost_matrix(tau_at_trident_frame)
    L_decay_to_rest_to_trident = np.matmul(L_rest_to_trident, L_decay_to_rest)
    return L_decay_to_rest_to_trident

def boost_fvs_to_trident(boosted_tau_array, decaying_tau_array, fv_array):
    """ Boost the four vectors, which are in the decaying tau frame, to the boosted tau (trident) frame"""
    boosted_fv_array = []
    nflags = 0
    for boosted_tau, decaying_tau, fv in zip(boosted_tau_array, decaying_tau_array, fv_array):
        L_decay_to_rest_to_trident = make_boost_matrix(decaying_tau, boosted_tau)
        boosted_fv = apply_boost(L_decay_to_rest_to_trident, fv)
        if boosted_fv[0] >= boosted_tau[0]:  # Flag events that don't conserve energy
            nflags += 1
            continue
        boosted_fv_array.append(boosted_fv)
    print(f"Flagged {nflags} events due to violation of energy conservation.")
    return boosted_fv_array

#def decay_tau(boosted_tau, fv):
#    """ Boost four vector fv to the frame of the boosted_tau and return boosted_fv"""
#    L = boost_matrix(boosted_tau)
#    boosted_fv = apply_boost(L, fv)
#    return boosted_fv
#
#def make_tau_decayed_four_vectors(boosted_tau_array, fv_array):
#    """ Produce array of four vectors boosted_fv for each boosted_tau and four vector fv in boosted_tau_array and fv_array """
#    boosted_fv_array = []
#    for boosted_tau, fv in zip(boosted_tau_array, fv_array):
#        boosted_fv = decay_tau(boosted_tau, fv)
#        boosted_fv_array.append(boosted_fv)
#    return boosted_fv_array

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
            writer.writerow(['-16', '2', boosted_nubar_tau[1], boosted_nubar_tau[2], boosted_nubar_tau[3], boosted_nubar_tau[0], '0.0'])
            writer.writerow(['</event>'])

def printout_decays(filename, all_final_state_particles, all_taus, part_dict, at_rest=True):
    with open(filename,'w') as textfile:
        if at_rest:
            print("********** Tau+ decays in the CM frame of the tau+ **********", file=textfile)
            i = 0
            for final_state_particles, tau in zip(all_final_state_particles, all_taus):
                print(f"Event {i+1}:", file=textfile)
                L_decay_to_rest = boost_to_cm_matrix(tau)
                tau_at_rest = apply_boost(L_decay_to_rest,tau)
                print(f"  Particle ID: -15, Particle: tau+, Momentum: E = {tau_at_rest[0]}, p = [{tau_at_rest[1]}, {tau_at_rest[2]}, {tau_at_rest[3]}]", file=textfile)
                for particle in final_state_particles:
                    p = np.array([particle.momentum[3], particle.momentum[0], particle.momentum[1], particle.momentum[2]])
                    particle_in_rest_frame = apply_boost(L_decay_to_rest,p)
                    print(f"  Particle ID: {particle.pid}, Particle: {part_dict[particle.pid]}, Momentum: E = {particle_in_rest_frame[0]}, p = [{particle_in_rest_frame[1]}, {particle_in_rest_frame[2]}, {particle_in_rest_frame[3]}]", file=textfile)
    #                print(f"  Particle ID: {particle.pid}, Particle: {part_dict[particle.pid]}, Momentum: E = {particle.momentum[3]}, p = [{particle.momentum[0]}, {particle.momentum[1]}, {particle.momentum[2]}]", file=textfile)
                i += 1
        else:
            print("********** Tau+ decays in the decay frame of the tau+ **********", file=textfile)
            i = 0
            for final_state_particles, tau in zip(all_final_state_particles, all_taus):
                print(f"Event {i+1}:", file=textfile)
                print(f"  Particle ID: -15, Particle: tau+, Momentum: E = {tau[0]}, p = [{tau[1]}, {tau[2]}, {tau[3]}]", file=textfile)
                for particle in final_state_particles:
                    print(f"  Particle ID: {particle.pid}, Particle: {part_dict[particle.pid]}, Momentum: E = {particle.momentum[3]}, p = [{particle.momentum[0]}, {particle.momentum[1]}, {particle.momentum[2]}]", file=textfile)
                i += 1


def printout_boosted_decays(filename, boosted_fv_array):
    with open(filename,'w') as textfile:
        for i, boosted_fv in enumerate(boosted_fv_array):
            print(f"Event {i+1}: E = {boosted_fv[0]}, p = [{boosted_fv[1]}, {boosted_fv[2]}, {boosted_fv[3]}]", file=textfile)

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
             130:'K0L',
             15:'tau-',
             -15:'tau+'}

TAU_DECAY_FILENAME = EVENTS_DIR_1TAU+'/../tau_decay/ve_e+_to_vtau_tau+/tag_3_pythia8_events.hepmc'

final_state_particles_list, tau_list = process_hepmc_file(TAU_DECAY_FILENAME)
pdg_to_store = [211,   # pi+
                -211,  # pi-
                321,   # K+
                -321,  # K-
                130,   # K0L
                22,    # gamma
                11,    # e-
                -11,   # e+
                ]

coherent_array, incoherent_array = read_tau_decays(final_state_particles_list, pdg_to_store)

#if write_dist:
#    printout_decays('tau+_hadronic_decays_REST_FRAME.txt', final_state_particles_list, tau_list, part_dict, at_rest=True)
#    printout_decays('tau+_hadronic_decays_DECAY_FRAME.txt', final_state_particles_list, tau_list, part_dict, at_rest=False)

print(len(coherent_array))
print(len(incoherent_array))

incoming_nu_mu_p, primary_nu_tau_p, primary_muon_p, boosted_tau_p, outgoing_nucleon_p = read_four_vectors(EVENTS_DIR_1TAU+f'/nucleon/proton/1tau_1e5Events_p_{Ev}GeV.txt')
incoming_nu_mu_n, primary_nu_tau_n, primary_muon_n, boosted_tau_n, outgoing_nucleon_n = read_four_vectors(EVENTS_DIR_1TAU+f'/nucleon/neutron/1tau_1e5Events_n_{Ev}GeV.txt')
incoming_nu_mu_Ar, primary_nu_tau_Ar, primary_muon_Ar, boosted_tau_Ar, outgoing_nucleon_Ar = read_four_vectors(EVENTS_DIR_1TAU+f'/coherent/argon/1tau_1e5Events_Ar_{Ev}GeV.txt')

boosted_nubar_tau_p  = boost_fvs_to_trident(boosted_tau_p, tau_list, incoherent_array)
boosted_nubar_tau_n  = boost_fvs_to_trident(boosted_tau_n, tau_list, incoherent_array)
boosted_nubar_tau_Ar  = boost_fvs_to_trident(boosted_tau_Ar, tau_list, coherent_array)

if write_dist:
    printout_boosted_decays('tau+_boosted_hadronic_decays_p_33GeV.txt', boosted_nubar_tau_p_33GeV)
    printout_boosted_decays('tau+_boosted_hadronic_decays_n_33GeV.txt', boosted_nubar_tau_n_33GeV)
    printout_boosted_decays('tau+_boosted_hadronic_decays_Ar_47GeV.txt', boosted_nubar_tau_Ar_47GeV)

printout_incoh(f'../../../csv/distributions/vmu_to_vtau_tau+_mu-/nucleon/proton/tau+_vtau_events/tau_hadronic_decayed_distribution_p_{Ev}GeV.txt',
               incoming_nu_mu_p, primary_nu_tau_p, primary_muon_p, boosted_nubar_tau_p, outgoing_nucleon_p,
               'proton')
printout_incoh(f'../../../csv/distributions/vmu_to_vtau_tau+_mu-/nucleon/neutron/tau+_vtau_events/tau_hadronic_decayed_distribution_n_{Ev}GeV.txt',
               incoming_nu_mu_n, primary_nu_tau_n, primary_muon_n, boosted_nubar_tau_n, outgoing_nucleon_n,
               'neutron')
printout_coh(f'../../../csv/distributions/vmu_to_vtau_tau+_mu-/coherent/argon/tau+_vtau_events/tau_hadronic_decayed_distribution_Ar_{Ev}GeV.txt',
               incoming_nu_mu_Ar, primary_nu_tau_Ar, primary_muon_Ar, boosted_nubar_tau_Ar)
