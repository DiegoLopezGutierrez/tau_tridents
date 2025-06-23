import csv
import numpy as np

CROSS_SECTION_DIR = '../csv/cross_sections'
A_Ar = 40
fb_to_cm2 = 1e-39 # 1 fb = 1e-39 cm^2

energy_vmuCC = []
xsec_vmuCC = []

with open(CROSS_SECTION_DIR+'/vmuCC/vmuCC_xsec_perE_Formaggio.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy = float(row[0])
        xsec = float(row[1]) * 10 * A_Ar * fb_to_cm2 * energy * 1e-4 # [1e-38 cm^2 / GeV] -> [m^2]

        energy_vmuCC.append(energy)
        xsec_vmuCC.append(xsec)

for i,energy in enumerate(np.geomspace(energy_vmuCC[-1],10000.,num=100)):
    if i == 0:
        # avoid double counting the first value
        continue
    # this value is roughly the value of xsec/E at the last point in the original file; the vmuCC xsec / E should be roughly linear at high energies.
    xsec = 0.6400880713923851*10.* A_Ar * fb_to_cm2 * energy * 1e-4
    energy_vmuCC.append(energy)
    xsec_vmuCC.append(xsec)

with open('vmuCC_xsec.csv','w') as csvout:
    writer = csv.writer(csvout, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for energy, xsec in zip(energy_vmuCC, xsec_vmuCC):
        writer.writerow([energy,xsec]) # file exported as [GeV, m^2]
