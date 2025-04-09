import numpy as np
from scipy.spatial.transform import Rotation as R
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import csv

n_chunks = 8 # this M1 2021 MacBook Pro laptop has 8 cores

def smear_mom_id(events, pid, sigma):
    """
    Smear the momentum of particles with a given PID by a Gaussian with width sigma.
    """
    # if pid is not a list, make it a list
    if not isinstance(pid, list):
        pid = [pid]
    if len(pid) == 0:
        raise ValueError("PID list is empty")
    # find all indices of particles with the given ID
    idx = np.isin(events[:,0], pid)
    # generate random numbers for smearing
    smear = np.random.normal(1, sigma, len(events))
    # for each event in our code, we have summed the relevant tau decay products four-momenta (hadronic or leptonic).
    # we need to get the "mass" of this particle.
    squared_masses = events[idx,5]**2 - events[idx,2]**2 - events[idx,3]**2 - events[idx,4]**2
    # apply smearing on Px, Py, Pz and P (index 2, 3, 4, 7)
    events[idx,2] *= smear[idx]
    events[idx,3] *= smear[idx]
    events[idx,4] *= smear[idx]
    events[idx,7] *= smear[idx]
    # update energy E = sqrt(P^2 + M^2)
    events[idx,5] = np.sqrt(events[idx,7]**2 + squared_masses)
    return events

def rotate_to_z(momentum):
    '''Returns the matrix that rotates a 3-vector to the z axis. Returns identity if already aligned with the z-axis.'''
    # normalize momentum
    mom = momentum/np.linalg.norm(momentum)
    z_axis = np.array([0,0,1])
    # check if momentum is already aligned with the z-axis
    if mom.all() == z_axis.all():
        return np.array([[1,0,0],[0,1,0],[0,0,1]])
    cos = np.dot(mom, z_axis)
    u_cross_v = np.cross(mom, z_axis)
    sin = np.linalg.norm(u_cross_v)
    rot = R.from_rotvec(np.arccos(cos)*u_cross_v/sin).as_matrix()
    return rot

def angle_smearing_chunk(chunk, sigma):
    """Generate random rotation matrices for a chunk. sigma is in degrees. """
    smear_theta = np.random.normal(0, sigma, len(chunk))
    smear_phi = np.random.uniform(0, 360, len(chunk))
    # rotate about y axis by smear_theta and then about z axis by smear_phi
    random_rotations = R.from_euler('yz', np.array([smear_theta, smear_phi]).T, degrees=True).as_matrix()
    return random_rotations

def angle_smearing(events, ids, sigma):
    '''Parallel version of angle smearing. Returns matrices that smear angles of particles with a given ID by a Gaussian width sigma.'''
    # Filter events based on PID
    idx = np.isin(events[:,0], ids)
    evt_smear = events[idx]

    # Parallel generation of random rotation matrices
    # n_chunks = 8  # Adjust based on your CPU cores
    chunks = np.array_split(evt_smear, n_chunks)
    
    random_rotations = np.vstack(
        Parallel(n_jobs=-1)(
            delayed(angle_smearing_chunk)(chunk, sigma) for chunk in chunks
        )
    )
    return random_rotations

# Define a function to perform the nested `einsum` operations on each chunk
def apply_einsum_chunk(rot_chunk, random_rot_chunk, momentum_chunk):
    """
    Perform the einsum operations on each chunk. rot_chunk rotates momentum_chunk to the z-axis then this z-axis chunk is rotated by the random smear rotations
    random_rot_chunk and finally we reverse the z-axis rotation done by rot_chunk using rot_chunk transposed on the randomzied z-axis chunk.
    We end up with our original momentum_chunk but rotated in place by the random smear theta and phi.
    """
    return np.einsum(
        'ijk,ik->ij',
        np.transpose(rot_chunk, (0, 2, 1)),
        np.einsum('ijk,ik->ij', random_rot_chunk, np.einsum('ijk,ik->ij', rot_chunk, momentum_chunk))
    )

def smear_angle(events, ids, sigmas):
    '''Smears the angle of particles with a given ID by a Gaussian with width sigma. Essentially combines the two functions above.'''
    idx = np.isin(events[:,0], ids)
    evt_smear = events[idx]
    momenta = evt_smear[:,2:5]
    # rotations to z - paralellized
    rots = np.vstack(
        Parallel(n_jobs=-1)(
            delayed(lambda chunk: np.apply_along_axis(rotate_to_z, axis=1, arr=chunk))(chunk)
            for chunk in np.array_split(momenta, n_chunks)
        )
    )
    print('Found rotations to z-axis')
    # random rotations
    random_rots = angle_smearing(evt_smear, ids, sigmas)
    print('Generated random rotations')
        
    # Split rots, random_rots, and momenta into chunks for parallel processing
    rots_chunks = np.array_split(rots, n_chunks)
    random_rots_chunks = np.array_split(random_rots, n_chunks)
    momenta_chunks = np.array_split(momenta, n_chunks)
    
    # Perform the einsum operations in parallel
    smeared_momenta_chunks = Parallel(n_jobs=-1)(
        delayed(apply_einsum_chunk)(rot_chunk, random_rot_chunk, momentum_chunk)
        for rot_chunk, random_rot_chunk, momentum_chunk in zip(rots_chunks, random_rots_chunks, momenta_chunks)
    )
    # Combine the results back into a single array
    smeared_momenta = np.vstack(smeared_momenta_chunks)
    evt_smear[:,2:5] = smeared_momenta

    # smear the momenta: rots.T . random_rots . rots . momenta for each 3-vector
    # evt_smear[:,2:5] = np.einsum('ijk,ik->ij', np.transpose(rots, (0, 2, 1)), np.einsum('ijk,ik->ij', random_rots, np.einsum('ijk,ik->ij', rots, momenta)))
    print('Multiplied all rotations')
    events[idx,2:5] = evt_smear[:,2:5]
    return events

############################################# Format of events ###################################################
#
#  0PID, 1status, 2px [GeV], 3py [GeV], 4pz [GeV], 5energy [GeV], 6mass [GeV](, 7p [GeV]; added during processing)
#
##################################################################################################################

#key_detector_effects = True
#key_background_processing = False
#
#if key_detector_effects:
#    COHERENT_HADRONIC_DECAY_FOLDER = "../csv/distributions/vmu_to_vtau_tau+_mu-/coherent/argon/tau+_vtau_events/hadronic"
#    PROTON_HADRONIC_DECAY_FOLDER = "../csv/distributions/vmu_to_vtau_tau+_mu-/nucleon/proton/tau+_vtau_events/hadronic"
#    NEUTRON_HADRONIC_DECAY_FOLDER = "../csv/distributions/vmu_to_vtau_tau+_mu-/nucleon/neutron/tau+_vtau_events/hadronic"
#
#    events = np.loadtxt(COHERENT_HADRONIC_DECAY_FOLDER+"/tau_hadronic_decayed_distribution_Ar_17GeV.txt", comments="<", delimiter=',')
#
#    # Add entry for particle momentum P = sqrt(Px^2 + Py^2 + Pz^2)
#    events = np.hstack((events, np.sqrt(events[:,2]**2 + events[:,3]**2 + events[:,4]**2).reshape(-1, 1)))
#
#    ######################################################
#    ###             Processing the events              ###
#    ######################################################
#
#    particles_1 = [13,-13] # mu-,mu+
#    particles_2 = [16,-16] # nutau, nutaubar
#    particles_3 = [2112,2212] # neutron, proton
#
#    # 1. Smear momentum (in %)
#    mom_smearing = [0.02, 0.10, 1]
#    events = smear_mom_id(events, particles_1, mom_smearing[0])
#    events = smear_mom_id(events, particles_2, mom_smearing[1])
#    events = smear_mom_id(events, particles_3, mom_smearing[2])
#
#    # 2. Smear angles
#    angle_smearing_values = [2, 10] # degrees
#    events = smear_angle(events, particles_1, angle_smearing_values[0])
#    events = smear_angle(events, particles_2, angle_smearing_values[1])
#
#    # save smeared coherent hadronic distribution
#    with open("smeared_tau_hadronic_decayed_distribution_Ar_17GeV.txt",'w',newline='') as outfile:
#        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
#        count_event = 0
#        for event in events:
#            #if count_event % 6 == 1:
#            if count_event == 0:
#                writer.writerow(['<event>'])
#                #count_event += 1
#            #if count_event % 6 == 0:
#            writer.writerow([f"{event[0]:.0f}", f"{event[1]:.0f}", f"{event[2]:.8f}", f"{event[3]:.8f}", f"{event[4]:.8f}", f"{event[5]:.8f}", f"{event[6]:.8f}", f"{event[7]:.8f}"])
#            count_event += 1
#
#            if count_event == 4:
#                writer.writerow(['</event>'])
#                #count_event += 1
#                count_event = 0            
#
