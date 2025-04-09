import numpy as np
from scipy.spatial.transform import Rotation as R
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import csv
from apply_detector_effects import smear_mom_id, rotate_to_z, angle_smearing_chunk, angle_smearing, apply_einsum_chunk, smear_angle

n_chunks = 8 # this M1 2021 MacBook Pro laptop has 8 cores

############################################# Format of events ###################################################
#
#  0PID, 1status, 2px [GeV], 3py [GeV], 4pz [GeV], 5energy [GeV], 6mass [GeV](, 7p [GeV]; added during processing)
#
##################################################################################################################

COHERENT_HADRONIC_DECAY_FOLDER = "../csv/distributions/vmu_to_vtau_tau+_mu-/coherent/argon/tau+_vtau_events/hadronic"
PROTON_HADRONIC_DECAY_FOLDER = "../csv/distributions/vmu_to_vtau_tau+_mu-/nucleon/proton/tau+_vtau_events/hadronic"
NEUTRON_HADRONIC_DECAY_FOLDER = "../csv/distributions/vmu_to_vtau_tau+_mu-/nucleon/neutron/tau+_vtau_events/hadronic"

key_detector_effects = True
target = input("Target [Ar,p,n]: ")
energies = np.array([5,8,11,14,17,20,24,27,30,33,36,39,42,45,49,52,55,58,61,64,67,71,74,77,80,83,86,89,92,96,99,102,105,108,111,114,118,121,124]) # GeV
#energies = [47]

if target == "Ar":
    for Enu in energies:
        events = np.loadtxt(COHERENT_HADRONIC_DECAY_FOLDER+f"/tau_hadronic_decayed_distribution_Ar_{Enu}GeV.txt", comments="<", delimiter=',')

        # Add entry for particle momentum P = sqrt(Px^2 + Py^2 + Pz^2)
        events = np.hstack((events, np.sqrt(events[:,2]**2 + events[:,3]**2 + events[:,4]**2).reshape(-1, 1)))

        ######################################################
        ###             Processing the events              ###
        ######################################################

        particles_1 = [13] # mu-
        particles_2 = [16,-16] # nutau, nutaubar

        # 1. Smear momentum (in %)
        mom_smearing = [0.05, 0.10]
        events = smear_mom_id(events, particles_1, mom_smearing[0])
        events = smear_mom_id(events, particles_2, mom_smearing[1])

        # 2. Smear angles
        angle_smearing_values = [2, 10] # degrees
        events = smear_angle(events, particles_1, angle_smearing_values[0])
        events = smear_angle(events, particles_2, angle_smearing_values[1])

        # save smeared coherent hadronic distribution
        with open(COHERENT_HADRONIC_DECAY_FOLDER+f"/../smeared/hadronic/smeared_tau_hadronic_decayed_distribution_Ar_{Enu}GeV.txt",'w',newline='') as outfile:
            writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            count_event = 0
            for event in events:
                # in coherent hadronic, all events contain just 14,16,13,-16 enclosed by <event> and </event>
                if count_event == 0:
                    writer.writerow(['<event>'])
                writer.writerow([f"{event[0]:.0f}",
                                 f"{event[1]:.0f}",
                                 f"{event[2]:.8f}",
                                 f"{event[3]:.8f}",
                                 f"{event[4]:.8f}",
                                 f"{event[5]:.8f}",
                                 f"{event[6]:.8f}",
                                 f"{event[7]:.8f}"])
                count_event += 1

                if count_event == 4:
                    writer.writerow(['</event>'])
                    count_event = 0            

if target == "p" or target == "n":
    for Enu in energies:
        if target == "p":
            NUCLEON_HADRONIC_DECAY_FOLDER = PROTON_HADRONIC_DECAY_FOLDER
            nucleon = "p"
        if target == "n":
            NUCLEON_HADRONIC_DECAY_FOLDER = NEUTRON_HADRONIC_DECAY_FOLDER
            nucleon = "n"

        events = np.loadtxt(NUCLEON_HADRONIC_DECAY_FOLDER+f"/tau_hadronic_decayed_distribution_{nucleon}_{Enu}GeV.txt", comments="<", delimiter=',')

        # Add entry for particle momentum P = sqrt(Px^2 + Py^2 + Pz^2)
        events = np.hstack((events, np.sqrt(events[:,2]**2 + events[:,3]**2 + events[:,4]**2).reshape(-1, 1)))

        ######################################################
        ###             Processing the events              ###
        ######################################################

        particles_1 = [13] # mu-,mu+
        particles_2 = [16,999,2212] # nutau, visible particles, proton
#        particles_3 = [2112,2212] # neutron, proton

        # 1. Smear momentum (in %)
#        mom_smearing = [0.05, 0.10, 1]
        mom_smearing = [0.05, 0.10]
        events = smear_mom_id(events, particles_1, mom_smearing[0])
        events = smear_mom_id(events, particles_2, mom_smearing[1])
 #       events = smear_mom_id(events, particles_3, mom_smearing[2])

        # 2. Smear angles
        angle_smearing_values = [2, 10] # degrees
        events = smear_angle(events, particles_1, angle_smearing_values[0])
        events = smear_angle(events, particles_2, angle_smearing_values[1])
        with open(NUCLEON_HADRONIC_DECAY_FOLDER+f"/../smeared/hadronic/smeared_tau_hadronic_decayed_distribution_{nucleon}_{Enu}GeV.txt",'w',newline='') as outfile:
            writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            count_event = 0
            for event in events:
                # in incoherent hadronic, all events contain just 14,16,13,999,nucleon_pid enclosed by <event> and </event>
                if count_event == 0:
                    writer.writerow(['<event>'])
                writer.writerow([f"{event[0]:.0f}",
                                 f"{event[1]:.0f}",
                                 f"{event[2]:.8f}",
                                 f"{event[3]:.8f}",
                                 f"{event[4]:.8f}",
                                 f"{event[5]:.8f}",
                                 f"{event[6]:.8f}",
                                 f"{event[7]:.8f}"])
                count_event += 1

                if count_event == 5:
                    writer.writerow(['</event>'])
                    count_event = 0            
