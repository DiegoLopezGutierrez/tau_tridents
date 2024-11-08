import matplotlib.pyplot as plt
import csv
import numpy as np

inverseGeV2_to_cm2 = (0.197e-13)**2 # 1 GeV^-1 = 0.197e-15 m

E_2mu_array = []
xsec_2mu_array = []

E_1tau1mu_array = []
xsec_1tau1mu_array = []

with open("./vmu_to_vmu_mu+_mu-_transverse_xsec.csv",'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        E_2mu_array.append(np.sqrt(float(row[0])))
        xsec_2mu_array.append(float(row[1])*inverseGeV2_to_cm2)

with open("./vmu_to_vtau_tau+_mu-_transverse_xsec.csv",'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        E_1tau1mu_array.append(np.sqrt(float(row[0])))
        xsec_1tau1mu_array.append(float(row[1])*inverseGeV2_to_cm2)

fig, ax = plt.subplots()

ax.plot(E_2mu_array, xsec_2mu_array, label='2mu')
ax.plot(E_1tau1mu_array, xsec_1tau1mu_array, label='1tau1mu')

ax.set_xlabel(r'$\sqrt{s}$ (GeV)')
ax.set_ylabel(r'$\sigma_{\gamma\nu}$ (cm^2)')
ax.legend()

fig.savefig('transverse_xsec.png')
