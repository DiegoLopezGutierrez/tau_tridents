import matplotlib.pyplot as plt
import csv
import numpy as np

inverseGeV2_to_cm2 = (0.197e-13)**2 # 1 GeV^-1 = 0.197e-15 m

E_2mu_array = []
xsec_2mu_array = []

E_2mu_WS_array = []
xsec_2mu_WS_array = []

with open("./EPA_xsec_2mu.csv",'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        E = float(row[0])
        xsec = float(row[1])
        E_2mu_array.append(E)
        xsec_2mu_array.append(xsec / E)

with open("./EPA_xsec_2mu_Woods-Saxon.csv",'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        E = float(row[0])
        xsec = float(row[1])
        E_2mu_WS_array.append(E)
        xsec_2mu_WS_array.append(xsec / E)

fig, ax = plt.subplots()

ax.plot(E_2mu_array, xsec_2mu_array, label='2mu')
ax.plot(E_2mu_WS_array, xsec_2mu_WS_array, label='2mu - Woods-Saxon')

ax.set_xlabel(r'$E_\nu$ (GeV)')
ax.set_ylabel(r'$\sigma / E_\nu$ (cm$^2$ GeV$^-1$)')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(1e-56,1e-36)
ax.set_xlim(1e-1,1e4)
ax.legend()

fig.savefig('EPA_xsec.png')
