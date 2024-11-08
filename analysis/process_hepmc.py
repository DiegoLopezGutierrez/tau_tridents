import numpy as np
import csv
import pyhepmc

TAU_DECAY_FILENAME = '../csv/events/tau_decay/hadronic/tag_1_pythia8_events.hepmc'

def is_final_state(particle):
    """ Check if the particle is in the final state """
    return particle.status == 1

def is_lepton(pid):
    lepton = {-11, 12, -12, 14, -14, 13, -13}
    return pid in lepton

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
            if particle.end_vertex:
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
                if contains_leptonic_decay(decay_products):
                    continue

            # Find final state particles from tau+ decay
            final_state_particles = find_final_state_particles(event, tau_plus)
            all_final_state_particles.append(final_state_particles)
        
    return all_final_state_particles

# Usage example
final_state_particles_list = process_hepmc_file(TAU_DECAY_FILENAME)
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

# Display information
with open('tau+_hadronic_decays.txt','w') as txtfile:
    for i, final_state_particles in enumerate(final_state_particles_list):
        print(f"Event {i+1}:", file=txtfile)
        for particle in final_state_particles:
            print(f"  Particle ID: {particle.pid}, Particle: {part_dict[particle.pid]}, Momentum: {particle.momentum}", file=txtfile)

    print("Total events processed:", len(final_state_particles_list), file=txtfile)
