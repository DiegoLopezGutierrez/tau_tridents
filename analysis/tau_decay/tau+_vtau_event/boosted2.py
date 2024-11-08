import numpy as np

def lorentz_boost_matrix(boost_vector):
    """ Calculate the Lorentz boost matrix for a given boost vector. """
    beta = np.array(boost_vector, dtype=np.float64)
    print("Beta = ", beta)
    beta2 = np.dot(beta, beta)
    print("Beta^2 = ", beta2)
    if beta2 >= 1.0:
        raise ValueError("Beta squared is greater than or equal to 1. Check the input values.")
    gamma = 1.0 / np.sqrt(1.0 - beta2)
    print("Relativistic factor = ", gamma)

    L = np.zeros((4, 4), dtype=np.float64)
    L[0, 0] = gamma
    L[1:, 0] = L[0, 1:] = -gamma * beta
    if beta2 > 0:
        L[1:, 1:] = np.eye(3, dtype=np.float64) + (gamma - 1) * np.outer(beta, beta) / beta2
    else:
        L[1:, 1:] = np.eye(3, dtype=np.float64)  # Identity matrix for zero beta
    
    print("Boost matrix = ", L)
   
    return L

def boost_from_rest(rest_four_vector, target_four_vector):
    """
    Boost the four-vector at rest to the frame defined by the target four-vector.
    
    Parameters:
    original_four_vector (array-like): The four-vector to be boosted [E, px, py, pz]
    target_four_vector (array-like): The four-vector defining the target frame [E, px, py, pz]

    Returns:
    np.ndarray: The boosted four-vector in the frame of the target four-vector
    """
    # Ensure high precision
    rest_four_vector = np.array(rest_four_vector, dtype=np.float64)
    target_four_vector = np.array(target_four_vector, dtype=np.float64)

    # Calculate the boost vector (velocity) from the target four-vector
    E_target, px_target, py_target, pz_target = target_four_vector
#    print("E_target = {}; px_target = {}, py_target = {}, pz_target = {}".format(E_target, px_target, py_target, pz_target))
    beta = np.array([px_target, py_target, pz_target], dtype=np.float64) / E_target
#    print("Beta = ", beta)
#    print("Rest four vector = ", rest_four_vector)
#    print("Target four vector = ", target_four_vector)
    
    # Calculate the Lorentz boost matrix
    L = lorentz_boost_matrix(beta)
    
    # Apply the boost to the original four-vector
    boosted_four_vector = np.dot(L, rest_four_vector)
#    print("Boosted rest four vector to target frame = ", boosted_four_vector)

    return boosted_four_vector

def boost_to_rest(original_four_vector, rest_four_vector):
    """
    Boost the four-vector at rest to the frame defined by the target four-vector.
    
    Parameters:
    original_four_vector (array-like): The four-vector to be boosted [E, px, py, pz]
    target_four_vector (array-like): The four-vector defining the target frame [E, px, py, pz]

    Returns:
    np.ndarray: The boosted four-vector in the frame of the target four-vector
    """
    # Ensure high precision
    original_four_vector = np.array(original_four_vector, dtype=np.float64)
    rest_four_vector = np.array(rest_four_vector, dtype=np.float64)

    # Calculate the boost vector (velocity) from the target four-vector
    E_original, px_original, py_original, pz_original = original_four_vector
#    print("E_original = {}; px_original = {}, py_original = {}, pz_original = {}".format(E_original, px_original, py_original, pz_original))
    beta = np.array([-px_original, -py_original, -pz_original], dtype=np.float64) / E_original
#    print("Beta = ", beta)
#    print("Original four vector = ", original_four_vector)
#    print("Rest four vector = ", rest_four_vector)
   
    # Calculate the Lorentz boost matrix
    L = lorentz_boost_matrix(beta)
    
    # Apply the boost to the original four-vector
    boosted_four_vector = np.dot(L, original_four_vector)
#    print("Boosted original four vector to rest frame = ", boosted_four_vector)

    return boosted_four_vector

def boost_to_frame(original_four_vector, target_four_vector):
    """
    Boost the original four-vector to the frame defined by the target four-vector.
    It is assumed that neither four vector is at rest.
    
    Parameters:
    original_four_vector (array-like): The four-vector to be boosted [E, px, py, pz]
    target_four_vector (array-like): The four-vector defining the target frame [E, px, py, pz]

    Returns:
    np.ndarray: The boosted four-vector in the frame of the target four-vector
    """
    # Ensure high precision
    original_four_vector = np.array(original_four_vector, dtype=np.float64)
    target_four_vector = np.array(target_four_vector, dtype=np.float64)

    # Calculate the boost vectors
    E_original, px_original, py_original, pz_original = original_four_vector
#    print("E_original = {}; px_original = {}, py_original = {}, pz_original = {}".format(E_original, px_original, py_original, pz_original))
    beta_to_rest = np.array([-px_original, -py_original, -pz_original], dtype=np.float64) / E_original

    E_target, px_target, py_target, pz_target = target_four_vector
#    print("E_target = {}; px_target = {}, py_target = {}, pz_target = {}".format(E_target, px_target, py_target, pz_target))
    beta_from_rest = np.array([px_target, py_target, pz_target], dtype=np.float64) / E_target

    L_to_rest = lorentz_boost_matrix(beta_to_rest)
    L_from_rest = lorentz_boost_matrix(beta_from_rest)

    # Boost original four vector to target frame
    original_at_target = np.dot(L_from_rest, np.dot(L_to_rest, original_four_vector))
#    print("Boosted original four vector to target frame = ", original_at_target)

    return original_at_target

def boost_to_frame_v2(original_four_vector, rest_four_vector, target_four_vector):
    """
    Boost the original four-vector to the frame defined by the target four-vector.
    It is assumed that neither four vector is at rest.
    
    Parameters:
    original_four_vector (array-like): The four-vector to be boosted [E, px, py, pz]
    target_four_vector (array-like): The four-vector defining the target frame [E, px, py, pz]

    Returns:
    np.ndarray: The boosted four-vector in the frame of the target four-vector
    """
    # Ensure high precision
    original_four_vector = np.array(original_four_vector, dtype=np.float64)
    target_four_vector = np.array(target_four_vector, dtype=np.float64)

    # Boost original four vector to the rest frame
    original_at_rest = boost_to_rest(original_four_vector, rest_four_vector)

    # Boost original four vector in the rest frame to the target frame
    original_at_target = boost_from_rest(original_at_rest, target_four_vector)

    return original_at_target

## Example usage:
tau_at_decay_frame = np.array([7.0830790682178829e+01,-2.0508444847000000e+01,2.1324713789000000e+01,6.4331201227999983e+01])
tau_at_trident_frame = np.array([22.75822814,-0.60381348,0.65876622,22.67106306])
tau_at_rest_frame = np.array([1.77800000,0.0,0.0,0.0])

# Boost tau in the decay frame to the rest frame
tau_decay_to_rest = boost_to_rest(tau_at_decay_frame, tau_at_rest_frame)
tau_decay_to_rest_to_trident = boost_from_rest(tau_decay_to_rest, tau_at_trident_frame)
tau_rest_to_trident = boost_from_rest(tau_at_rest_frame, tau_at_trident_frame)

print("Tau+ in the decay frame: ", tau_at_decay_frame)
print("Tau+ in the rest frame: ", tau_at_rest_frame)
print("Tau+ in the trident frame: ", tau_at_trident_frame)

print("Boosting tau+ in the decay frame to the rest frame:")
print("\tResult: ", tau_decay_to_rest)
print("\tExpected: ", tau_at_rest_frame)
print("\tDifference: ", tau_decay_to_rest - tau_at_rest_frame)

print("Boosting tau+ in the decay frame to the rest frame to the trident frame:")
print("\tResult: ", tau_decay_to_rest_to_trident)
print("\tExpected: ", tau_at_trident_frame)
print("\tDifference: ", tau_decay_to_rest_to_trident - tau_at_trident_frame)

print("Boosting tau+ in the rest frame to the trident frame:")
print("\tResult: ", tau_rest_to_trident)
print("\tExpected: ", tau_at_trident_frame)
print("\tDifference: ", tau_rest_to_trident - tau_at_trident_frame)

## Example usage:
## Define the original four-vector [E, px, py, pz] in its original frame
#decayed_four_momentum = np.array([5.0, 1.0, 1.0, 1.0])
#
## Define the target four-vector [E, px, py, pz] which represents the new frame
#trident_four_momentum = np.array([10.0, 2.0, 2.0, 2.0])
#
## Perform the boost
#boosted_four_momentum = boost_to_frame(decayed_four_momentum, trident_four_momentum)
#
#print("Original four-vector in its frame:", decayed_four_momentum)
#print("Target four-vector (frame to boost to):", trident_four_momentum)
#print("Boosted four-vector in the target frame:", boosted_four_momentum)
