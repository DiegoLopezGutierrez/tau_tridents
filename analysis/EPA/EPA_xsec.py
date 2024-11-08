import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline, Akima1DInterpolator
import csv
import matplotlib.pyplot as plt

#### Constants ####
GF = 0.0000116637 # Fermi constant [1/GeV^2]
aem = 1/137.036 # EM coupling
me = 5.109e-4
me2 = me*me
mmu = 0.105
mmu2 = mmu*mmu
mtau = 1.778 # Tau lepton mass [GeV]
mtau2 = mtau*mtau

MArgon = 39.95*0.9315 # Mass of argon nucleus
Mproton = 0.938272 # Mass of proton
Mneutron = 0.939565 # Mass of neutron

Z_Ar = 18
A_Ar = 40

LambdaQCD = 0.217 # Assuming 200 MeV for now
QMax = LambdaQCD / A_Ar**(1./3.)

invGeV2_to_cm2 = (1/5.06*1e-13)**2 # 1 GeV^-1 = 1/5.06 * 10^-13 cm

FORM_FACTOR_ARGON_FILENAME = "../../csv/form_factors/FF_Ar_3Fp-red_Alt.csv"
TRANSVERSE_XSEC_1TAU1MU_FILENAME = "./xsec_splines/vmu_to_vtau_tau+_mu-_transverse_xsec.csv"
#TRANSVERSE_XSEC_2MU_FILENAME = "./xsec_splines/vmu_to_vmu_mu+_mu-_transverse_xsec.csv"
TRANSVERSE_XSEC_2MU_FILENAME = "../../csv/cross_sections/vmu_gamma_to_vmu_mu+_mu-_xsec/BeacomZhou_vgamma_2mu_xsec_per_s.csv"

#### Initialize Arrays ####
q2_data_array = []
FF_integrand_data_array = []
s_data_array = []
transverse_xsec_data_array = []

#### Woods-Saxon Form Factor ####
### Form factor function from Woods-Saxon distribution ###
fm_to_invGeV = 5.068 # 1 fm = 5.068 GeV^-1
a = 0.523 * fm_to_invGeV # fm -> GeV^-1

def WS_form_factor(A, Q):
    if Q == 0:   # First value of Q2s is 0.0 which will lead to division by zero when integrating
        Q = 0.00001
    r0 = 1.126 * A**(1/3) * fm_to_invGeV   # Form factor should be dimensionless. Convert fm to GeV^-1.
    Q2 = Q*Q

    # WS form factor can be written analytically; see (A.3) in Ballett et al.
    WS = (3 * np.pi * a) / (r0**2 + np.pi**2 * a**2)

    WS *= np.pi * a * np.tanh(np.pi * Q * a)**(-1) * np.sin(Q * r0) - r0 * np.cos(Q * r0)

    WS /= Q * r0 * np.sinh(np.pi * Q * a)

    return WS

# Normalize to F(0) = 1; use 0.001 instead.
Ar_norm  = WS_form_factor(A_Ar, 0.00001)

#### Retrieve Data ####
with open(FORM_FACTOR_ARGON_FILENAME,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        q = float(row[0])
        FF = float(row[1])
        q2 = q*q
        integrand = FF*FF/q2
        if (len(q2_data_array) > 1 and q2 == q2_data_array[-1]):
            FF_integrand_data_array[-1] = integrand
            continue
        q2_data_array.append(q2)
        FF_integrand_data_array.append(integrand)

with open(TRANSVERSE_XSEC_2MU_FILENAME,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        s = float(row[0]) # GeV^2
        transverse_xsec_per_s = float(row[1]) # cm^2 / GeV^2
        transverse_xsec = transverse_xsec_per_s * s # cm^2
        s_data_array.append(s)
        transverse_xsec_data_array.append(transverse_xsec)

transverse_xsec_interpolator = CubicSpline(s_data_array, transverse_xsec_data_array)
form_factor_integrand_interpolator = Akima1DInterpolator(q2_data_array, FF_integrand_data_array)

fig, ax = plt.subplots(1,1)
ax.plot(q2_data_array, FF_integrand_data_array, label="data")
ax.plot(np.linspace(0.01, 0.5, num=100), form_factor_integrand_interpolator(np.linspace(0.01, 0.5, num=100)), label="Akima1DInterpolator")
ax.set_xlabel("Q2")
ax.set_ylabel("F^2/Q2")
ax.set_yscale('log')
ax.legend()
fig.savefig("form_factor_interpolator.png")

def form_factor_integral(s, Ev, ml):
    lower_limit = (s/(2*Ev))**2
    upper_limit = s - (ml)**2
    
    if (lower_limit > upper_limit):
        return 0.0
    q2_array = np.linspace(lower_limit, upper_limit, num=200)
    FF_integrand_array = [(WS_form_factor(A_Ar, q) / Ar_norm)**2 / (q*q) for q in q2_array]
    integral = simpson(FF_integrand_array, x=q2_array)
    return integral


#def form_factor_integral(s, Ev, ml):
#    lower_limit = (s/(2*Ev))**2
#    upper_limit = s - (ml)**2
#    min_q2 = min(q2_data_array)
#    max_q2 = max(q2_data_array)
#    
#    if (upper_limit > max_q2):
#        upper_limit = max_q2
#    if (lower_limit < min_q2):
#        lower_limit = min_q2
#    if (lower_limit > upper_limit):
#        return 0.0
##    print("Form factor integral...")
##    print("\tLower limit = ", lower_limit)
##    print("\tUpper limit = ", upper_limit)
#    q2_array = np.linspace(lower_limit, upper_limit, num=100)
#    FF_integrand_array = form_factor_integrand_interpolator(q2_array)
#    integral = simpson(FF_integrand_array, x=q2_array)
##    print("\tIntegral = ", integral)
#    return integral

#def form_factor_integral(s, Ev):
#    lower_limit = (s/(2*Ev))**2
#    upper_limit = QMax*QMax
#    if (lower_limit > upper_limit):
#        return 0.0
##    print("Form factor integral...")
##    print("\tLower limit = ", lower_limit)
##    print("\tUpper limit = ", upper_limit)
#    q2_array = np.linspace(lower_limit, upper_limit, num=100)
#    FF_integrand_array = [1/q2 for q2 in q2_array]
#    integral = simpson(FF_integrand_array, x=q2_array)
##    print("\tIntegral = ", integral)
#    return integral

#### Neutrino-nucleus Cross Section ####
def slimits(Ev, MTarget, ml):
    smax = Ev / (2*Ev + MTarget) *(ml*ml + 2*Ev*MTarget + np.sqrt((2*Ev*MTarget - ml*ml)**2 - 4*MTarget**2*ml*ml))
    smin = Ev / (2*Ev + MTarget) *(ml*ml + 2*Ev*MTarget - np.sqrt((2*Ev*MTarget - ml*ml)**2 - 4*MTarget**2*ml*ml))
    return smin, smax

def neutrino_nucleus_integrand(Ev, s, ml):
    integrand = Z_Ar**2*aem/s/np.pi
    integrand *= transverse_xsec_interpolator(s)
    integrand *= form_factor_integral(s, Ev, ml)
    return integrand

def neutrino_nucleus_xsec(Ev, MTarget, ml):
    lower_limit, upper_limit = slimits(Ev, MTarget, ml)
    upper_limit = MTarget*(MTarget + 2*Ev)  # This should be the center-of-mass of s = (p1 + P)^2 = 0 + MTarget^2 + 2p1*P = MTarget^2 + 2*Ev*MTarget
#    upper_limit = 2*Ev*QMax
#    print("Neutrino nucleus xsec...")
#    print("\tLower limit = ", lower_limit, " GeV^2")
#    print("\tUpper limit = ", upper_limit, " GeV^2")
    s_array = np.linspace(lower_limit, upper_limit, num=300)
    neutrino_nucleus_integrand_array = [neutrino_nucleus_integrand(Ev, s, ml) for s in s_array]
    xsec = simpson(neutrino_nucleus_integrand_array, x=s_array)
    return xsec
#    return xsec*invGeV2_to_cm2

Ev = float(input("Enter neutrino energy in GeV: "))
target_material = input("Enter target material [Ar/Fe/W]: ")
outfile = input("Enter name of output file: ")
if target_material == "Ar":
    MTarget = MArgon
    ml = mmu + mmu
    E_start = float(input("Enter start energy in GeV: "))
    E_end = float(input("Enter end energy in GeV: "))
    E_array = np.linspace(E_start, E_end, 500)
    xsec_array = [neutrino_nucleus_xsec(E, MTarget, ml) / E for E in E_array]
    #xsec = neutrino_nucleus_xsec(Ev, MTarget, ml)
    #xsec_perE = xsec/Ev
    #print(f"Calculated cross section over energy for Ev = {Ev} GeV off {target_material} is ", xsec_perE)

with open(outfile,'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for E, xsec in zip(E_array, xsec_array):
        writer.writerow([E,xsec])
