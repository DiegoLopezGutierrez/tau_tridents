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
   p = -np.array([px, py, pz])
   p_magnitude = np.linalg.norm(p)

   if p_magnitude == 0:
       return np.eye(4)

   beta = p / E
   gamma = lorentz_factor(beta)
   beta_magnitude = np.linalg.norm(beta)

   L = np.zeros((4, 4))
   L[0, 0] = gamma
   L[1:, 0] = -gamma * beta
   L[0, 1:] = -gamma * beta
   L[1:, 1:] = np.eye(3) + (gamma - 1) * np.outer(beta, beta) / (beta_magnitude ** 2)

   return L

def apply_boost(boost_matrix, four_vector):
   """
   Apply the boost matrix to a given four-vector.
   """
   return np.dot(boost_matrix, four_vector)


# Example usage
#angle = 0.8
#cos = np.cos(angle)
#sin = np.sin(angle)
#p = 0.142829
#mother_fv = np.array([0.2, 0, p*sin, p*cos])  # Replace with your four-vector
#
## Calculate boost matrix
#L = boost_matrix(mother_fv)
#
#print("Boost Matrix:")
#
#print(L)
#
## cross check
#mother_at_rest = np.array([0.140, 0, 0, 0])
#boosted_mother = apply_boost(L, mother_at_rest)
#print("Boosted mother Four-Vector:")
#print(boosted_mother)
#print("Difference:")
#print(boosted_mother - mother_fv)
#
#
## Apply boost to another four-vector
#daughter_fv = np.array([0.030625,0.,0.0119259,0.0282075])  # Replace with your four-vector
#boosted_daughter = apply_boost(L, daughter_fv)
#print("Boosted daughter Four-Vector:")
#
#print(boosted_daughter)

# how to apply that to tau decays

# L = boost_matrix(tau_fv) # ==> boosted tau (from Diego)

# nu_tau_boosted = apply_boost(L, nu_tau_fv) # (from Bhupal)

# nu_mu_boosted = apply_boost(L, nu_mu_fv) # (from Bhupal)

# muon_boosted = apply_boost(L, muon_fv) # (from Bhupal)

#boosted_tau_filename = ''
#
#for boosted_tau in boosted_tau_array:
#    L = boost_matrix(boosted_tau)
#    for nu_tau, nu_mu, muon in zip(nu_tau_array, nu_mu_array, muon_array):
#        boosted_nu_tau = apply_boost(L, nu_tau)
#        boosted_nu_mu  = apply_boost(L, nu_mu)
#        boosted_nu_tau_array.append(boosted_nu_tau)
#    for  
