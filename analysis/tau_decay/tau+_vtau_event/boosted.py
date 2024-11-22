import numpy as np

def lorentz_factor(v):
    """
    Calculate the Lorentz factor (gamma) for a given velocity vector.
    """
    c = 1  # speed of light in natural units
    v2 = np.dot(v, v)
    return 1 / np.sqrt(1 - v2 / c**2)

def boost_matrix(four_vector):
    """
    Calculate the boost matrix for a given four-vector.
    """
    c = 1  # speed of light in natural units
    E, px, py, pz = four_vector
    p = -np.array([px, py, pz], dtype=np.float64)
    p_magnitude = np.linalg.norm(p)
    
    if p_magnitude == 0:
        return np.eye(4)
    
    beta = p / E
    gamma = lorentz_factor(beta)
    beta_magnitude = np.linalg.norm(beta)
    
    L = np.zeros((4, 4), dtype=np.float64)
    L[0, 0] = gamma
    L[1:, 0] = -gamma * beta
    L[0, 1:] = -gamma * beta
    L[1:, 1:] = np.eye(3, dtype=np.float64) + (gamma - 1) * np.outer(beta, beta) / (beta_magnitude ** 2)
    
    return L

def boost_to_cm_matrix(four_vector):
    """
    Calculate the boost matrix for a given four-vector.
    """
    c = 1  # speed of light in natural units
    E, px, py, pz = four_vector
    p = np.array([px, py, pz], dtype=np.float64) # the sign here is the only difference with the previous function
    p_magnitude = np.linalg.norm(p)
    
    if p_magnitude == 0:
        return np.eye(4)
    
    beta = p / E
    gamma = lorentz_factor(beta)
    beta_magnitude = np.linalg.norm(beta)
    
    L = np.zeros((4, 4), dtype=np.float64)
    L[0, 0] = gamma
    L[1:, 0] = -gamma * beta
    L[0, 1:] = -gamma * beta
    L[1:, 1:] = np.eye(3, dtype=np.float64) + (gamma - 1) * np.outer(beta, beta) / (beta_magnitude ** 2)
    
    return L


def apply_boost(boost_matrix, four_vector):
    """
    Apply the boost matrix to a given four-vector.
    """
    return np.dot(boost_matrix, four_vector)


## When printing, only show 3 decimal places and suppress scientific notation
#np.set_printoptions(precision=3, suppress=True)
#
## Example usage
#angle = 0.8
#cos = np.cos(angle)
#sin = np.sin(angle)
#p = 0.142829
#mother_fv = np.array([0.2, 0, p*sin, p*cos], dtype=np.float64)  # Replace with your four-vector
#
## Calculate boost matrix
#L = boost_matrix(mother_fv)
## print("Boost Matrix:")
## print(L)
#
## cross check
#mother_at_rest = np.array([0.140, 0, 0, 0], dtype=np.float64)
#boosted_mother = apply_boost(L, mother_at_rest)
#print("************* First check: take mother at rest and boost it by matrix obtained with mother fv, compare with mother fv")
#print("Boosted mother Four-Vector:")
#print(boosted_mother)
#print("Difference:")
#print(boosted_mother - mother_fv)
#
#print("************* Second check: take mother fv and boost it back to CM frame")
#L_cm = boost_to_cm_matrix(mother_fv)
#mother_in_cm = apply_boost(L_cm, boosted_mother)
#print("Mother in CM frame:")
#print(mother_in_cm)
#
#
## Apply boost to another four-vector
#daughter_fv = np.array([0.030625,0.,0.0119259,0.0282075], dtype=np.float64)  # Replace with your four-vector
#boosted_daughter = apply_boost(L, daughter_fv)
#print("************* Boost daughter by matrix obtained with mother fv")
#print("Boosted daughter Four-Vector:")
#print(boosted_daughter)
