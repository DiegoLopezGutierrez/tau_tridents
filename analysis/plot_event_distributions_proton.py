import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker
import csv

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

EVENTS_DIR_1TAU = '../csv/events/vmu_to_vtau_tau+_mu-/nucleon/proton'

RMiss_background_filename = '../csv/events/AlexSousa_RMiss_Background.csv'
energy_RMiss = []
background_RMiss = []

vertical_flag = True
horizontal_flag = False
write_dist = False


############################
#### Read Distributions ####
############################

def read_after_skip(file_path, E_v, skip_lines=37):
    with open(file_path, 'r') as txtfile:
        for _ in range(skip_lines):
            next(txtfile)
        pt_miss_array = []
        pt_muon_array = []
        pt_tau_array = []
        q2_array = []
        theta_array = []
        for line in txtfile:
            if line.startswith('16'):  # Read and store vtau transverse momentum
                data = line.split()
                assert int(data[0]) ==  16
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                pt_miss = np.sqrt(px**2 + py**2)
                theta = math.acos(pz / (np.sqrt(pz**2 + pt_miss**2)))
                pt_miss_array.append(pt_miss)
                theta_array.append(math.degrees(theta))
            if line.startswith('13'):  # Read and store muon transverse momentum
                data = line.split()
                assert int(data[0]) == 13
                px = float(data[2])
                py = float(data[3])
                pz = float(data[4])
                E  = float(data[5])
                M  = float(data[6])
                pt_muon = np.sqrt(px**2 + py**2)
                q2 = M**2 - E_v*(E - pz)  # Incoming v has E_v = pz_v and px_v = py_v = 0.
                pt_muon_array.append(pt_muon)
                q2_array.append(-q2)
            if line.startswith('-15'):  # Read and store tau transverse momentum
                data = line.split()
                assert int(data[0]) == -15
                px = float(data[2])
                py = float(data[3])
                pt_tau = np.sqrt(px**2 + py**2)
                pt_tau_array.append(pt_tau)
        r_miss_array = [pt_miss / (pt_miss + pt_muon) for pt_miss, pt_muon in zip(pt_miss_array, pt_muon_array)]  # Calculate missing momentum ratio
        return pt_miss_array, pt_muon_array, pt_tau_array, r_miss_array, q2_array, theta_array

############################
###### Distributions #######
############################

ptMiss_1tau_nucleon_p_5GeV, ptMuon_1tau_nucleon_p_5GeV, ptTau_1tau_nucleon_p_5GeV, RMiss_1tau_nucleon_p_5GeV, Fermi4Q2_1tau_nucleon_p_5GeV, ThetaMiss_1tau_nucleon_p_5GeV = read_after_skip(EVENTS_DIR_1TAU+'/1tau_1e5Events_p_5GeV.txt', 5)
ptMiss_1tau_nucleon_p_10GeV, ptMuon_1tau_nucleon_p_10GeV, ptTau_1tau_nucleon_p_10GeV, RMiss_1tau_nucleon_p_10GeV, Fermi4Q2_1tau_nucleon_p_10GeV, ThetaMiss_1tau_nucleon_p_10GeV = read_after_skip(EVENTS_DIR_1TAU+'/1tau_1e5Events_p_10GeV.txt', 10)
ptMiss_1tau_nucleon_p_20GeV, ptMuon_1tau_nucleon_p_20GeV, ptTau_1tau_nucleon_p_20GeV, RMiss_1tau_nucleon_p_20GeV, Fermi4Q2_1tau_nucleon_p_20GeV, ThetaMiss_1tau_nucleon_p_20GeV = read_after_skip(EVENTS_DIR_1TAU+'/1tau_1e5Events_p_20GeV.txt', 20)
ptMiss_1tau_nucleon_p_50GeV, ptMuon_1tau_nucleon_p_50GeV, ptTau_1tau_nucleon_p_50GeV, RMiss_1tau_nucleon_p_50GeV, Fermi4Q2_1tau_nucleon_p_50GeV, ThetaMiss_1tau_nucleon_p_50GeV = read_after_skip(EVENTS_DIR_1TAU+'/1tau_1e5Events_p_50GeV.txt', 50)
ptMiss_1tau_nucleon_p_200GeV, ptMuon_1tau_nucleon_p_200GeV, ptTau_1tau_nucleon_p_200GeV, RMiss_1tau_nucleon_p_200GeV, Fermi4Q2_1tau_nucleon_p_200GeV, ThetaMiss_1tau_nucleon_p_200GeV = read_after_skip(EVENTS_DIR_1TAU+'/1tau_1e5Events_p_200GeV.txt', 200)

if write_dist:
    with open('ptDist_1tau_nucleon_p_5GeV.txt','w',newline='') as outfile:
        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for pt_miss, pt_muon, pt_tau, r_miss, theta_miss in zip(ptMiss_1tau_nucleon_p_5GeV, ptMuon_1tau_nucleon_p_5GeV, ptTau_1tau_nucleon_p_5GeV, RMiss_1tau_nucleon_p_5GeV, ThetaMiss_1tau_nucleon_p_5GeV):
            writer.writerow([pt_miss, pt_muon, pt_tau, r_miss, theta_miss])

    with open('ptDist_1tau_nucleon_p_10GeV.txt','w',newline='') as outfile:
        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for pt_miss, pt_muon, pt_tau, r_miss, theta_miss in zip(ptMiss_1tau_nucleon_p_10GeV, ptMuon_1tau_nucleon_p_10GeV, ptTau_1tau_nucleon_p_10GeV, RMiss_1tau_nucleon_p_10GeV, ThetaMiss_1tau_nucleon_p_10GeV):
            writer.writerow([pt_miss, pt_muon, pt_tau, r_miss, theta_miss])

    with open('ptDist_1tau_nucleon_p_20GeV.txt','w',newline='') as outfile:
        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for pt_miss, pt_muon, pt_tau, r_miss, theta_miss in zip(ptMiss_1tau_nucleon_p_20GeV, ptMuon_1tau_nucleon_p_20GeV, ptTau_1tau_nucleon_p_20GeV, RMiss_1tau_nucleon_p_20GeV, ThetaMiss_1tau_nucleon_p_20GeV):
            writer.writerow([pt_miss, pt_muon, pt_tau, r_miss, theta_miss])

    with open('ptDist_1tau_nucleon_p_50GeV.txt','w',newline='') as outfile:
        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for pt_miss, pt_muon, pt_tau, r_miss, theta_miss in zip(ptMiss_1tau_nucleon_p_50GeV, ptMuon_1tau_nucleon_p_50GeV, ptTau_1tau_nucleon_p_50GeV, RMiss_1tau_nucleon_p_50GeV, ThetaMiss_1tau_nucleon_p_50GeV):
            writer.writerow([pt_miss, pt_muon, pt_tau, r_miss, theta_miss])

    with open('ptDist_1tau_nucleon_p_200GeV.txt','w',newline='') as outfile:
        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for pt_miss, pt_muon, pt_tau, r_miss, theta_miss in zip(ptMiss_1tau_nucleon_p_200GeV, ptMuon_1tau_nucleon_p_200GeV, ptTau_1tau_nucleon_p_200GeV, RMiss_1tau_nucleon_p_200GeV, ThetaMiss_1tau_nucleon_p_200GeV):
            writer.writerow([pt_miss, pt_muon, pt_tau, r_miss, theta_miss])

with open(RMiss_background_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        energy_RMiss.append(float(row[0]))
        background_RMiss.append(float(row[1]))

##################
#### Binning ####
##################
NBINS = 40

# vmu -> vtau tau+ mu- #
print("vmu -> vtau tau+ mu-")
print("---------- 5 GeV ----------")
ptMiss_1tau_nucleon_p_5GeV_max = max(ptMiss_1tau_nucleon_p_5GeV)
ptMiss_1tau_nucleon_p_5GeV_min = min(ptMiss_1tau_nucleon_p_5GeV)
ptMiss_1tau_nucleon_p_5GeV_bins, ptMiss_1tau_nucleon_p_5GeV_step = np.linspace(ptMiss_1tau_nucleon_p_5GeV_min, ptMiss_1tau_nucleon_p_5GeV_max, num=NBINS, retstep=True)
print("ptMiss:")
print("\tmin: ", ptMiss_1tau_nucleon_p_5GeV_min)
print("\tmax: ", ptMiss_1tau_nucleon_p_5GeV_max)
ptMuon_1tau_nucleon_p_5GeV_max = max(ptMuon_1tau_nucleon_p_5GeV)
ptMuon_1tau_nucleon_p_5GeV_min = min(ptMuon_1tau_nucleon_p_5GeV)
ptMuon_1tau_nucleon_p_5GeV_bins, ptMuon_1tau_nucleon_p_5GeV_step = np.linspace(ptMuon_1tau_nucleon_p_5GeV_min, ptMuon_1tau_nucleon_p_5GeV_max, num=NBINS, retstep=True)
print("ptMuon:")
print("\tmin: ", ptMuon_1tau_nucleon_p_5GeV_min)
print("\tmax: ", ptMuon_1tau_nucleon_p_5GeV_max)
ptTau_1tau_nucleon_p_5GeV_max = max(ptTau_1tau_nucleon_p_5GeV)
ptTau_1tau_nucleon_p_5GeV_min = min(ptTau_1tau_nucleon_p_5GeV)
ptTau_1tau_nucleon_p_5GeV_bins, ptTau_1tau_nucleon_p_5GeV_step = np.linspace(ptTau_1tau_nucleon_p_5GeV_min, ptTau_1tau_nucleon_p_5GeV_max, num=NBINS, retstep=True)
print("ptTau:")
print("\tmin: ", ptTau_1tau_nucleon_p_5GeV_min)
print("\tmax: ", ptTau_1tau_nucleon_p_5GeV_max)
RMiss_1tau_nucleon_p_5GeV_max = max(RMiss_1tau_nucleon_p_5GeV)
RMiss_1tau_nucleon_p_5GeV_min = min(RMiss_1tau_nucleon_p_5GeV)
RMiss_1tau_nucleon_p_5GeV_bins, RMiss_1tau_nucleon_p_5GeV_step = np.linspace(0, 1, num=NBINS, retstep=True)
print("bin width: ", RMiss_1tau_nucleon_p_5GeV_step)
print("RMiss:")
print("\tmin: ", RMiss_1tau_nucleon_p_5GeV_min)
print("\tmax: ", RMiss_1tau_nucleon_p_5GeV_max)
Fermi4Q2_1tau_nucleon_p_5GeV_max = max(Fermi4Q2_1tau_nucleon_p_5GeV)
Fermi4Q2_1tau_nucleon_p_5GeV_min = min(Fermi4Q2_1tau_nucleon_p_5GeV)
Fermi4Q2_1tau_nucleon_p_5GeV_bins, Fermi4Q2_1tau_nucleon_p_5GeV_step = np.linspace(Fermi4Q2_1tau_nucleon_p_5GeV_min, Fermi4Q2_1tau_nucleon_p_5GeV_max, num=NBINS, retstep=True)
print("Fermi4Q2:")
print("\tmin: ", Fermi4Q2_1tau_nucleon_p_5GeV_min)
print("\tmax: ", Fermi4Q2_1tau_nucleon_p_5GeV_max)
ThetaMiss_1tau_nucleon_p_5GeV_max = max(ThetaMiss_1tau_nucleon_p_5GeV)
ThetaMiss_1tau_nucleon_p_5GeV_min = min(ThetaMiss_1tau_nucleon_p_5GeV)
ThetaMiss_1tau_nucleon_p_5GeV_bins, ThetaMiss_1tau_nucleon_p_5GeV_step = np.linspace(0,180, num=NBINS, retstep=True)
print("ThetaMiss:")
print("\tmin: ", ThetaMiss_1tau_nucleon_p_5GeV_min)
print("\tmax: ", ThetaMiss_1tau_nucleon_p_5GeV_max)

print("---------- 10 GeV ----------")
ptMiss_1tau_nucleon_p_10GeV_max = max(ptMiss_1tau_nucleon_p_10GeV)
ptMiss_1tau_nucleon_p_10GeV_min = min(ptMiss_1tau_nucleon_p_10GeV)
ptMiss_1tau_nucleon_p_10GeV_bins, ptMiss_1tau_nucleon_p_10GeV_step = np.linspace(ptMiss_1tau_nucleon_p_10GeV_min, ptMiss_1tau_nucleon_p_10GeV_max, num=NBINS, retstep=True)
print("ptMiss:")
print("\tmin: ", ptMiss_1tau_nucleon_p_10GeV_min)
print("\tmax: ", ptMiss_1tau_nucleon_p_10GeV_max)
ptMuon_1tau_nucleon_p_10GeV_max = max(ptMuon_1tau_nucleon_p_10GeV)
ptMuon_1tau_nucleon_p_10GeV_min = min(ptMuon_1tau_nucleon_p_10GeV)
ptMuon_1tau_nucleon_p_10GeV_bins, ptMuon_1tau_nucleon_p_10GeV_step = np.linspace(ptMuon_1tau_nucleon_p_10GeV_min, ptMuon_1tau_nucleon_p_10GeV_max, num=NBINS, retstep=True)
print("ptMuon:")
print("\tmin: ", ptMuon_1tau_nucleon_p_10GeV_min)
print("\tmax: ", ptMuon_1tau_nucleon_p_10GeV_max)
ptTau_1tau_nucleon_p_10GeV_max = max(ptTau_1tau_nucleon_p_10GeV)
ptTau_1tau_nucleon_p_10GeV_min = min(ptTau_1tau_nucleon_p_10GeV)
ptTau_1tau_nucleon_p_10GeV_bins, ptTau_1tau_nucleon_p_10GeV_step = np.linspace(ptTau_1tau_nucleon_p_10GeV_min, ptTau_1tau_nucleon_p_10GeV_max, num=NBINS, retstep=True)
print("ptTau:")
print("\tmin: ", ptTau_1tau_nucleon_p_10GeV_min)
print("\tmax: ", ptTau_1tau_nucleon_p_10GeV_max)
RMiss_1tau_nucleon_p_10GeV_max = max(RMiss_1tau_nucleon_p_10GeV)
RMiss_1tau_nucleon_p_10GeV_min = min(RMiss_1tau_nucleon_p_10GeV)
RMiss_1tau_nucleon_p_10GeV_bins, RMiss_1tau_nucleon_p_10GeV_step = np.linspace(0, 1, num=NBINS, retstep=True)
print("RMiss:")
print("\tmin: ", RMiss_1tau_nucleon_p_10GeV_min)
print("\tmax: ", RMiss_1tau_nucleon_p_10GeV_max)
Fermi4Q2_1tau_nucleon_p_10GeV_max = max(Fermi4Q2_1tau_nucleon_p_10GeV)
Fermi4Q2_1tau_nucleon_p_10GeV_min = min(Fermi4Q2_1tau_nucleon_p_10GeV)
Fermi4Q2_1tau_nucleon_p_10GeV_bins, Fermi4Q2_1tau_nucleon_p_10GeV_step = np.linspace(Fermi4Q2_1tau_nucleon_p_10GeV_min, Fermi4Q2_1tau_nucleon_p_10GeV_max, num=NBINS, retstep=True)
print("Fermi4Q2:")
print("\tmin: ", Fermi4Q2_1tau_nucleon_p_10GeV_min)
print("\tmax: ", Fermi4Q2_1tau_nucleon_p_10GeV_max)
ThetaMiss_1tau_nucleon_p_10GeV_max = max(ThetaMiss_1tau_nucleon_p_10GeV)
ThetaMiss_1tau_nucleon_p_10GeV_min = min(ThetaMiss_1tau_nucleon_p_10GeV)
ThetaMiss_1tau_nucleon_p_10GeV_bins, ThetaMiss_1tau_nucleon_p_10GeV_step = np.linspace(0, 180, num=NBINS, retstep=True)
print("ThetaMiss:")
print("\tmin: ", ThetaMiss_1tau_nucleon_p_10GeV_min)
print("\tmax: ", ThetaMiss_1tau_nucleon_p_10GeV_max)

print("---------- 20 GeV ----------")
ptMiss_1tau_nucleon_p_20GeV_max = max(ptMiss_1tau_nucleon_p_20GeV)
ptMiss_1tau_nucleon_p_20GeV_min = min(ptMiss_1tau_nucleon_p_20GeV)
ptMiss_1tau_nucleon_p_20GeV_bins, ptMiss_1tau_nucleon_p_20GeV_step = np.linspace(ptMiss_1tau_nucleon_p_20GeV_min, ptMiss_1tau_nucleon_p_20GeV_max, num=NBINS, retstep=True)
print("ptMiss:")
print("\tmin: ", ptMiss_1tau_nucleon_p_20GeV_min)
print("\tmax: ", ptMiss_1tau_nucleon_p_20GeV_max)
ptMuon_1tau_nucleon_p_20GeV_max = max(ptMuon_1tau_nucleon_p_20GeV)
ptMuon_1tau_nucleon_p_20GeV_min = min(ptMuon_1tau_nucleon_p_20GeV)
ptMuon_1tau_nucleon_p_20GeV_bins, ptMuon_1tau_nucleon_p_20GeV_step = np.linspace(ptMuon_1tau_nucleon_p_20GeV_min, ptMuon_1tau_nucleon_p_20GeV_max, num=NBINS, retstep=True)
print("ptMuon:")
print("\tmin: ", ptMuon_1tau_nucleon_p_20GeV_min)
print("\tmax: ", ptMuon_1tau_nucleon_p_20GeV_max)
ptTau_1tau_nucleon_p_20GeV_max = max(ptTau_1tau_nucleon_p_20GeV)
ptTau_1tau_nucleon_p_20GeV_min = min(ptTau_1tau_nucleon_p_20GeV)
ptTau_1tau_nucleon_p_20GeV_bins, ptTau_1tau_nucleon_p_20GeV_step = np.linspace(ptTau_1tau_nucleon_p_20GeV_min, ptTau_1tau_nucleon_p_20GeV_max, num=NBINS, retstep=True)
print("ptTau:")
print("\tmin: ", ptTau_1tau_nucleon_p_20GeV_min)
print("\tmax: ", ptTau_1tau_nucleon_p_20GeV_max)
RMiss_1tau_nucleon_p_20GeV_max = max(RMiss_1tau_nucleon_p_20GeV)
RMiss_1tau_nucleon_p_20GeV_min = min(RMiss_1tau_nucleon_p_20GeV)
RMiss_1tau_nucleon_p_20GeV_bins, RMiss_1tau_nucleon_p_20GeV_step = np.linspace(0, 1, num=NBINS, retstep=True)
print("RMiss:")
print("\tmin: ", RMiss_1tau_nucleon_p_20GeV_min)
print("\tmax: ", RMiss_1tau_nucleon_p_20GeV_max)
Fermi4Q2_1tau_nucleon_p_20GeV_max = max(Fermi4Q2_1tau_nucleon_p_20GeV)
Fermi4Q2_1tau_nucleon_p_20GeV_min = min(Fermi4Q2_1tau_nucleon_p_20GeV)
Fermi4Q2_1tau_nucleon_p_20GeV_bins, Fermi4Q2_1tau_nucleon_p_20GeV_step = np.linspace(Fermi4Q2_1tau_nucleon_p_20GeV_min, Fermi4Q2_1tau_nucleon_p_20GeV_max, num=NBINS, retstep=True)
print("Fermi4Q2:")
print("\tmin: ", Fermi4Q2_1tau_nucleon_p_20GeV_min)
print("\tmax: ", Fermi4Q2_1tau_nucleon_p_20GeV_max)
ThetaMiss_1tau_nucleon_p_20GeV_max = max(ThetaMiss_1tau_nucleon_p_20GeV)
ThetaMiss_1tau_nucleon_p_20GeV_min = min(ThetaMiss_1tau_nucleon_p_20GeV)
ThetaMiss_1tau_nucleon_p_20GeV_bins, ThetaMiss_1tau_nucleon_p_20GeV_step = np.linspace(0, 180, num=NBINS, retstep=True)
print("ThetaMiss:")
print("\tmin: ", ThetaMiss_1tau_nucleon_p_20GeV_min)
print("\tmax: ", ThetaMiss_1tau_nucleon_p_20GeV_max)

print("---------- 50 GeV ----------")
ptMiss_1tau_nucleon_p_50GeV_max = max(ptMiss_1tau_nucleon_p_50GeV)
ptMiss_1tau_nucleon_p_50GeV_min = min(ptMiss_1tau_nucleon_p_50GeV)
ptMiss_1tau_nucleon_p_50GeV_bins, ptMiss_1tau_nucleon_p_50GeV_step = np.linspace(ptMiss_1tau_nucleon_p_50GeV_min, ptMiss_1tau_nucleon_p_50GeV_max, num=NBINS, retstep=True)
print("ptMiss:")
print("\tmin: ", ptMiss_1tau_nucleon_p_50GeV_min)
print("\tmax: ", ptMiss_1tau_nucleon_p_50GeV_max)
ptMuon_1tau_nucleon_p_50GeV_max = max(ptMuon_1tau_nucleon_p_50GeV)
ptMuon_1tau_nucleon_p_50GeV_min = min(ptMuon_1tau_nucleon_p_50GeV)
ptMuon_1tau_nucleon_p_50GeV_bins, ptMuon_1tau_nucleon_p_50GeV_step = np.linspace(ptMuon_1tau_nucleon_p_50GeV_min, ptMuon_1tau_nucleon_p_50GeV_max, num=NBINS, retstep=True)
print("ptMuon:")
print("\tmin: ", ptMuon_1tau_nucleon_p_50GeV_min)
print("\tmax: ", ptMuon_1tau_nucleon_p_50GeV_max)
ptTau_1tau_nucleon_p_50GeV_max = max(ptTau_1tau_nucleon_p_50GeV)
ptTau_1tau_nucleon_p_50GeV_min = min(ptTau_1tau_nucleon_p_50GeV)
ptTau_1tau_nucleon_p_50GeV_bins, ptTau_1tau_nucleon_p_50GeV_step = np.linspace(ptTau_1tau_nucleon_p_50GeV_min, ptTau_1tau_nucleon_p_50GeV_max, num=NBINS, retstep=True)
print("ptTau:")
print("\tmin: ", ptTau_1tau_nucleon_p_50GeV_min)
print("\tmax: ", ptTau_1tau_nucleon_p_50GeV_max)
RMiss_1tau_nucleon_p_50GeV_max = max(RMiss_1tau_nucleon_p_50GeV)
RMiss_1tau_nucleon_p_50GeV_min = min(RMiss_1tau_nucleon_p_50GeV)
RMiss_1tau_nucleon_p_50GeV_bins, RMiss_1tau_nucleon_p_50GeV_step = np.linspace(0, 1, num=NBINS, retstep=True)
print("RMiss:")
print("\tmin: ", RMiss_1tau_nucleon_p_50GeV_min)
print("\tmax: ", RMiss_1tau_nucleon_p_50GeV_max)
Fermi4Q2_1tau_nucleon_p_50GeV_max = max(Fermi4Q2_1tau_nucleon_p_50GeV)
Fermi4Q2_1tau_nucleon_p_50GeV_min = min(Fermi4Q2_1tau_nucleon_p_50GeV)
Fermi4Q2_1tau_nucleon_p_50GeV_bins, Fermi4Q2_1tau_nucleon_p_50GeV_step = np.linspace(Fermi4Q2_1tau_nucleon_p_50GeV_min, Fermi4Q2_1tau_nucleon_p_50GeV_max, num=NBINS, retstep=True)
print("Fermi4Q2:")
print("\tmin: ", Fermi4Q2_1tau_nucleon_p_50GeV_min)
print("\tmax: ", Fermi4Q2_1tau_nucleon_p_50GeV_max)
ThetaMiss_1tau_nucleon_p_50GeV_max = max(ThetaMiss_1tau_nucleon_p_50GeV)
ThetaMiss_1tau_nucleon_p_50GeV_min = min(ThetaMiss_1tau_nucleon_p_50GeV)
ThetaMiss_1tau_nucleon_p_50GeV_bins, ThetaMiss_1tau_nucleon_p_50GeV_step = np.linspace(0, 180, num=NBINS, retstep=True)
print("ThetaMiss:")
print("\tmin: ", ThetaMiss_1tau_nucleon_p_50GeV_min)
print("\tmax: ", ThetaMiss_1tau_nucleon_p_50GeV_max)

print("---------- 200 GeV ----------")
ptMiss_1tau_nucleon_p_200GeV_max = max(ptMiss_1tau_nucleon_p_200GeV)
ptMiss_1tau_nucleon_p_200GeV_min = min(ptMiss_1tau_nucleon_p_200GeV)
ptMiss_1tau_nucleon_p_200GeV_bins, ptMiss_1tau_nucleon_p_200GeV_step = np.linspace(ptMiss_1tau_nucleon_p_200GeV_min, ptMiss_1tau_nucleon_p_200GeV_max, num=NBINS, retstep=True)
print("ptMiss:")
print("\tmin: ", ptMiss_1tau_nucleon_p_200GeV_min)
print("\tmax: ", ptMiss_1tau_nucleon_p_200GeV_max)
ptMuon_1tau_nucleon_p_200GeV_max = max(ptMuon_1tau_nucleon_p_200GeV)
ptMuon_1tau_nucleon_p_200GeV_min = min(ptMuon_1tau_nucleon_p_200GeV)
ptMuon_1tau_nucleon_p_200GeV_bins, ptMuon_1tau_nucleon_p_200GeV_step = np.linspace(ptMuon_1tau_nucleon_p_200GeV_min, ptMuon_1tau_nucleon_p_200GeV_max, num=NBINS, retstep=True)
print("ptMuon:")
print("\tmin: ", ptMuon_1tau_nucleon_p_200GeV_min)
print("\tmax: ", ptMuon_1tau_nucleon_p_200GeV_max)
ptTau_1tau_nucleon_p_200GeV_max = max(ptTau_1tau_nucleon_p_200GeV)
ptTau_1tau_nucleon_p_200GeV_min = min(ptTau_1tau_nucleon_p_200GeV)
ptTau_1tau_nucleon_p_200GeV_bins, ptTau_1tau_nucleon_p_200GeV_step = np.linspace(ptTau_1tau_nucleon_p_200GeV_min, ptTau_1tau_nucleon_p_200GeV_max, num=NBINS, retstep=True)
print("ptTau:")
print("\tmin: ", ptTau_1tau_nucleon_p_200GeV_min)
print("\tmax: ", ptTau_1tau_nucleon_p_200GeV_max)
RMiss_1tau_nucleon_p_200GeV_max = max(RMiss_1tau_nucleon_p_200GeV)
RMiss_1tau_nucleon_p_200GeV_min = min(RMiss_1tau_nucleon_p_200GeV)
RMiss_1tau_nucleon_p_200GeV_bins, RMiss_1tau_nucleon_p_200GeV_step = np.linspace(0, 1, num=NBINS, retstep=True)
print("RMiss:")
print("\tmin: ", RMiss_1tau_nucleon_p_200GeV_min)
print("\tmax: ", RMiss_1tau_nucleon_p_200GeV_max)
Fermi4Q2_1tau_nucleon_p_200GeV_max = max(Fermi4Q2_1tau_nucleon_p_200GeV)
Fermi4Q2_1tau_nucleon_p_200GeV_min = min(Fermi4Q2_1tau_nucleon_p_200GeV)
Fermi4Q2_1tau_nucleon_p_200GeV_bins, Fermi4Q2_1tau_nucleon_p_200GeV_step = np.linspace(Fermi4Q2_1tau_nucleon_p_200GeV_min, Fermi4Q2_1tau_nucleon_p_200GeV_max, num=NBINS, retstep=True)
print("Fermi4Q2:")
print("\tmin: ", Fermi4Q2_1tau_nucleon_p_200GeV_min)
print("\tmax: ", Fermi4Q2_1tau_nucleon_p_200GeV_max)
ThetaMiss_1tau_nucleon_p_200GeV_max = max(ThetaMiss_1tau_nucleon_p_200GeV)
ThetaMiss_1tau_nucleon_p_200GeV_min = min(ThetaMiss_1tau_nucleon_p_200GeV)
ThetaMiss_1tau_nucleon_p_200GeV_bins, ThetaMiss_1tau_nucleon_p_200GeV_step = np.linspace(0, 180, num=NBINS, retstep=True)
print("ThetaMiss:")
print("\tmin: ", ThetaMiss_1tau_nucleon_p_200GeV_min)
print("\tmax: ", ThetaMiss_1tau_nucleon_p_200GeV_max)

### Weights ###
NEVENTS = 1e5

# 5 GeV #
wts_ptMiss_1tau_nucleon_p_5GeV = np.ones_like(ptMiss_1tau_nucleon_p_5GeV) / (NEVENTS * ptMiss_1tau_nucleon_p_5GeV_step)
wts_ptMuon_1tau_nucleon_p_5GeV = np.ones_like(ptMuon_1tau_nucleon_p_5GeV) / (NEVENTS * ptMuon_1tau_nucleon_p_5GeV_step)
wts_ptTau_1tau_nucleon_p_5GeV = np.ones_like(ptTau_1tau_nucleon_p_5GeV) / (NEVENTS * ptTau_1tau_nucleon_p_5GeV_step)
wts_RMiss_1tau_nucleon_p_5GeV = np.ones_like(RMiss_1tau_nucleon_p_5GeV) / (NEVENTS * RMiss_1tau_nucleon_p_5GeV_step)
wts_Fermi4Q2_1tau_nucleon_p_5GeV = np.ones_like(Fermi4Q2_1tau_nucleon_p_5GeV) / (NEVENTS * Fermi4Q2_1tau_nucleon_p_5GeV_step)
wts_ThetaMiss_1tau_nucleon_p_5GeV = np.ones_like(ThetaMiss_1tau_nucleon_p_5GeV) / (NEVENTS * ThetaMiss_1tau_nucleon_p_5GeV_step)

# 10 GeV #
wts_ptMiss_1tau_nucleon_p_10GeV = np.ones_like(ptMiss_1tau_nucleon_p_10GeV) / (NEVENTS * ptMiss_1tau_nucleon_p_10GeV_step)
wts_ptMuon_1tau_nucleon_p_10GeV = np.ones_like(ptMuon_1tau_nucleon_p_10GeV) / (NEVENTS * ptMuon_1tau_nucleon_p_10GeV_step)
wts_ptTau_1tau_nucleon_p_10GeV = np.ones_like(ptTau_1tau_nucleon_p_10GeV) / (NEVENTS * ptTau_1tau_nucleon_p_10GeV_step)
wts_RMiss_1tau_nucleon_p_10GeV = np.ones_like(RMiss_1tau_nucleon_p_10GeV) / (NEVENTS * RMiss_1tau_nucleon_p_10GeV_step)
wts_Fermi4Q2_1tau_nucleon_p_10GeV = np.ones_like(Fermi4Q2_1tau_nucleon_p_10GeV) / (NEVENTS * Fermi4Q2_1tau_nucleon_p_10GeV_step)
wts_ThetaMiss_1tau_nucleon_p_10GeV = np.ones_like(ThetaMiss_1tau_nucleon_p_10GeV) / (NEVENTS * ThetaMiss_1tau_nucleon_p_10GeV_step)

# 20 GeV #
wts_ptMiss_1tau_nucleon_p_20GeV = np.ones_like(ptMiss_1tau_nucleon_p_20GeV) / (NEVENTS * ptMiss_1tau_nucleon_p_20GeV_step)
wts_ptMuon_1tau_nucleon_p_20GeV = np.ones_like(ptMuon_1tau_nucleon_p_20GeV) / (NEVENTS * ptMuon_1tau_nucleon_p_20GeV_step)
wts_ptTau_1tau_nucleon_p_20GeV = np.ones_like(ptTau_1tau_nucleon_p_20GeV) / (NEVENTS * ptTau_1tau_nucleon_p_20GeV_step)
wts_RMiss_1tau_nucleon_p_20GeV = np.ones_like(RMiss_1tau_nucleon_p_20GeV) / (NEVENTS * RMiss_1tau_nucleon_p_20GeV_step)
wts_Fermi4Q2_1tau_nucleon_p_20GeV = np.ones_like(Fermi4Q2_1tau_nucleon_p_20GeV) / (NEVENTS * Fermi4Q2_1tau_nucleon_p_20GeV_step)
wts_ThetaMiss_1tau_nucleon_p_20GeV = np.ones_like(ThetaMiss_1tau_nucleon_p_20GeV) / (NEVENTS * ThetaMiss_1tau_nucleon_p_20GeV_step)

# 50 GeV #
wts_ptMiss_1tau_nucleon_p_50GeV = np.ones_like(ptMiss_1tau_nucleon_p_50GeV) / (NEVENTS * ptMiss_1tau_nucleon_p_50GeV_step)
wts_ptMuon_1tau_nucleon_p_50GeV = np.ones_like(ptMuon_1tau_nucleon_p_50GeV) / (NEVENTS * ptMuon_1tau_nucleon_p_50GeV_step)
wts_ptTau_1tau_nucleon_p_50GeV = np.ones_like(ptTau_1tau_nucleon_p_50GeV) / (NEVENTS * ptTau_1tau_nucleon_p_50GeV_step)
wts_RMiss_1tau_nucleon_p_50GeV = np.ones_like(RMiss_1tau_nucleon_p_50GeV) / (NEVENTS * RMiss_1tau_nucleon_p_50GeV_step)
wts_Fermi4Q2_1tau_nucleon_p_50GeV = np.ones_like(Fermi4Q2_1tau_nucleon_p_50GeV) / (NEVENTS * Fermi4Q2_1tau_nucleon_p_50GeV_step)
wts_ThetaMiss_1tau_nucleon_p_50GeV = np.ones_like(ThetaMiss_1tau_nucleon_p_50GeV) / (NEVENTS * ThetaMiss_1tau_nucleon_p_50GeV_step)

# 200 GeV #
wts_ptMiss_1tau_nucleon_p_200GeV = np.ones_like(ptMiss_1tau_nucleon_p_200GeV) / (NEVENTS * ptMiss_1tau_nucleon_p_200GeV_step)
wts_ptMuon_1tau_nucleon_p_200GeV = np.ones_like(ptMuon_1tau_nucleon_p_200GeV) / (NEVENTS * ptMuon_1tau_nucleon_p_200GeV_step)
wts_ptTau_1tau_nucleon_p_200GeV = np.ones_like(ptTau_1tau_nucleon_p_200GeV) / (NEVENTS * ptTau_1tau_nucleon_p_200GeV_step)
wts_RMiss_1tau_nucleon_p_200GeV = np.ones_like(RMiss_1tau_nucleon_p_200GeV) / (NEVENTS * RMiss_1tau_nucleon_p_200GeV_step)
wts_Fermi4Q2_1tau_nucleon_p_200GeV = np.ones_like(Fermi4Q2_1tau_nucleon_p_200GeV) / (NEVENTS * Fermi4Q2_1tau_nucleon_p_200GeV_step)
wts_ThetaMiss_1tau_nucleon_p_200GeV = np.ones_like(ThetaMiss_1tau_nucleon_p_200GeV) / (NEVENTS * ThetaMiss_1tau_nucleon_p_200GeV_step)

##################
#### Plotting ####
##################

# vertical #
if vertical_flag:
    #fig, ax = plt.subplots(4, 2, figsize=(28, 36), tight_layout=True)
    fig1, (ax11, ax12, ax13, ax14, ax15) = plt.subplots(1, 5, sharey=True, figsize=(60,15), tight_layout=True, gridspec_kw={'wspace':0}) # ptMiss
    fig2, (ax21, ax22, ax23, ax24, ax25) = plt.subplots(1, 5, sharey=True, figsize=(60,15), tight_layout=True, gridspec_kw={'wspace':0}) # ptMuon
    fig3, (ax31, ax32, ax33, ax34, ax35) = plt.subplots(1, 5, sharey=True, figsize=(60,15), tight_layout=True, gridspec_kw={'wspace':0}) # ptTau
    fig4, (ax41, ax42, ax43, ax44, ax45) = plt.subplots(1, 5, sharey=True, figsize=(60,15), tight_layout=True, gridspec_kw={'wspace':0}) # RMiss
    fig5, (ax51, ax52, ax53, ax54, ax55) = plt.subplots(1, 5, sharey=True, figsize=(60,15), tight_layout=True, gridspec_kw={'wspace':0}) # Fermi4Q2
    fig6, (ax61, ax62, ax63, ax64, ax65) = plt.subplots(1, 5, sharey=True, figsize=(60,15), tight_layout=True, gridspec_kw={'wspace':0}) # ThetaMiss

    ### ptMiss ###
    ax11.hist(ptMiss_1tau_nucleon_p_5GeV, density=False, bins=ptMiss_1tau_nucleon_p_5GeV_bins, weights=wts_ptMiss_1tau_nucleon_p_5GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax11.hist(ptMiss_1tau_nucleon_p_5GeV, density=False, bins=ptMiss_1tau_nucleon_p_5GeV_bins, weights=wts_ptMiss_1tau_nucleon_p_5GeV, histtype='step', color='black', alpha=1, lw=2)
    ax12.hist(ptMiss_1tau_nucleon_p_10GeV, density=False, bins=ptMiss_1tau_nucleon_p_10GeV_bins, weights=wts_ptMiss_1tau_nucleon_p_10GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax12.hist(ptMiss_1tau_nucleon_p_10GeV, density=False, bins=ptMiss_1tau_nucleon_p_10GeV_bins, weights=wts_ptMiss_1tau_nucleon_p_10GeV, histtype='step', color='black', alpha=1, lw=2)
    ax13.hist(ptMiss_1tau_nucleon_p_20GeV, density=False, bins=ptMiss_1tau_nucleon_p_20GeV_bins, weights=wts_ptMiss_1tau_nucleon_p_20GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax13.hist(ptMiss_1tau_nucleon_p_20GeV, density=False, bins=ptMiss_1tau_nucleon_p_20GeV_bins, weights=wts_ptMiss_1tau_nucleon_p_20GeV, histtype='step', color='black', alpha=1, lw=2)
    ax14.hist(ptMiss_1tau_nucleon_p_50GeV, density=False, bins=ptMiss_1tau_nucleon_p_50GeV_bins, weights=wts_ptMiss_1tau_nucleon_p_50GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax14.hist(ptMiss_1tau_nucleon_p_50GeV, density=False, bins=ptMiss_1tau_nucleon_p_50GeV_bins, weights=wts_ptMiss_1tau_nucleon_p_50GeV, histtype='step', color='black', alpha=1, lw=2)
    ax15.hist(ptMiss_1tau_nucleon_p_200GeV, density=False, bins=ptMiss_1tau_nucleon_p_200GeV_bins, weights=wts_ptMiss_1tau_nucleon_p_200GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax15.hist(ptMiss_1tau_nucleon_p_200GeV, density=False, bins=ptMiss_1tau_nucleon_p_200GeV_bins, weights=wts_ptMiss_1tau_nucleon_p_200GeV, histtype='step', color='black', alpha=1, lw=2)
   
    ### ptMuon ###
    ax21.hist(ptMuon_1tau_nucleon_p_5GeV, density=False, bins=ptMuon_1tau_nucleon_p_5GeV_bins, weights=wts_ptMuon_1tau_nucleon_p_5GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax21.hist(ptMuon_1tau_nucleon_p_5GeV, density=False, bins=ptMuon_1tau_nucleon_p_5GeV_bins, weights=wts_ptMuon_1tau_nucleon_p_5GeV, histtype='step', color='black', alpha=1, lw=2)
    ax22.hist(ptMuon_1tau_nucleon_p_10GeV, density=False, bins=ptMuon_1tau_nucleon_p_10GeV_bins, weights=wts_ptMuon_1tau_nucleon_p_10GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax22.hist(ptMuon_1tau_nucleon_p_10GeV, density=False, bins=ptMuon_1tau_nucleon_p_10GeV_bins, weights=wts_ptMuon_1tau_nucleon_p_10GeV, histtype='step', color='black', alpha=1, lw=2)
    ax23.hist(ptMuon_1tau_nucleon_p_20GeV, density=False, bins=ptMuon_1tau_nucleon_p_20GeV_bins, weights=wts_ptMuon_1tau_nucleon_p_20GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax23.hist(ptMuon_1tau_nucleon_p_20GeV, density=False, bins=ptMuon_1tau_nucleon_p_20GeV_bins, weights=wts_ptMuon_1tau_nucleon_p_20GeV, histtype='step', color='black', alpha=1, lw=2)
    ax24.hist(ptMuon_1tau_nucleon_p_50GeV, density=False, bins=ptMuon_1tau_nucleon_p_50GeV_bins, weights=wts_ptMuon_1tau_nucleon_p_50GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax24.hist(ptMuon_1tau_nucleon_p_50GeV, density=False, bins=ptMuon_1tau_nucleon_p_50GeV_bins, weights=wts_ptMuon_1tau_nucleon_p_50GeV, histtype='step', color='black', alpha=1, lw=2)
    ax25.hist(ptMuon_1tau_nucleon_p_200GeV, density=False, bins=ptMuon_1tau_nucleon_p_200GeV_bins, weights=wts_ptMuon_1tau_nucleon_p_200GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax25.hist(ptMuon_1tau_nucleon_p_200GeV, density=False, bins=ptMuon_1tau_nucleon_p_200GeV_bins, weights=wts_ptMuon_1tau_nucleon_p_200GeV, histtype='step', color='black', alpha=1, lw=2)
  
    ### ptTau ###
    ax31.hist(ptTau_1tau_nucleon_p_5GeV, density=False, bins=ptTau_1tau_nucleon_p_5GeV_bins, weights=wts_ptTau_1tau_nucleon_p_5GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax31.hist(ptTau_1tau_nucleon_p_5GeV, density=False, bins=ptTau_1tau_nucleon_p_5GeV_bins, weights=wts_ptTau_1tau_nucleon_p_5GeV, histtype='step', color='black', alpha=1, lw=2)
    ax32.hist(ptTau_1tau_nucleon_p_10GeV, density=False, bins=ptTau_1tau_nucleon_p_10GeV_bins, weights=wts_ptTau_1tau_nucleon_p_10GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax32.hist(ptTau_1tau_nucleon_p_10GeV, density=False, bins=ptTau_1tau_nucleon_p_10GeV_bins, weights=wts_ptTau_1tau_nucleon_p_10GeV, histtype='step', color='black', alpha=1, lw=2)
    ax33.hist(ptTau_1tau_nucleon_p_20GeV, density=False, bins=ptTau_1tau_nucleon_p_20GeV_bins, weights=wts_ptTau_1tau_nucleon_p_20GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax33.hist(ptTau_1tau_nucleon_p_20GeV, density=False, bins=ptTau_1tau_nucleon_p_20GeV_bins, weights=wts_ptTau_1tau_nucleon_p_20GeV, histtype='step', color='black', alpha=1, lw=2)
    ax34.hist(ptTau_1tau_nucleon_p_50GeV, density=False, bins=ptTau_1tau_nucleon_p_50GeV_bins, weights=wts_ptTau_1tau_nucleon_p_50GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax34.hist(ptTau_1tau_nucleon_p_50GeV, density=False, bins=ptTau_1tau_nucleon_p_50GeV_bins, weights=wts_ptTau_1tau_nucleon_p_50GeV, histtype='step', color='black', alpha=1, lw=2)
    ax35.hist(ptTau_1tau_nucleon_p_200GeV, density=False, bins=ptTau_1tau_nucleon_p_200GeV_bins, weights=wts_ptTau_1tau_nucleon_p_200GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax35.hist(ptTau_1tau_nucleon_p_200GeV, density=False, bins=ptTau_1tau_nucleon_p_200GeV_bins, weights=wts_ptTau_1tau_nucleon_p_200GeV, histtype='step', color='black', alpha=1, lw=2)

    ### RMiss ###
    ax41.hist(RMiss_1tau_nucleon_p_5GeV, density=False, bins=RMiss_1tau_nucleon_p_5GeV_bins, weights=wts_RMiss_1tau_nucleon_p_5GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax41.hist(RMiss_1tau_nucleon_p_5GeV, density=False, bins=RMiss_1tau_nucleon_p_5GeV_bins, weights=wts_RMiss_1tau_nucleon_p_5GeV, histtype='step', color='black', alpha=1, lw=2)
    ax41.hist(energy_RMiss, density=False, bins=RMiss_1tau_nucleon_p_5GeV_bins, weights=background_RMiss, histtype='step',color='red',alpha=1,lw=2)
    ax42.hist(RMiss_1tau_nucleon_p_10GeV, density=False, bins=RMiss_1tau_nucleon_p_10GeV_bins, weights=wts_RMiss_1tau_nucleon_p_10GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax42.hist(RMiss_1tau_nucleon_p_10GeV, density=False, bins=RMiss_1tau_nucleon_p_10GeV_bins, weights=wts_RMiss_1tau_nucleon_p_10GeV, histtype='step', color='black', alpha=1, lw=2)
    ax43.hist(RMiss_1tau_nucleon_p_20GeV, density=False, bins=RMiss_1tau_nucleon_p_20GeV_bins, weights=wts_RMiss_1tau_nucleon_p_20GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax43.hist(RMiss_1tau_nucleon_p_20GeV, density=False, bins=RMiss_1tau_nucleon_p_20GeV_bins, weights=wts_RMiss_1tau_nucleon_p_20GeV, histtype='step', color='black', alpha=1, lw=2)
    ax44.hist(RMiss_1tau_nucleon_p_50GeV, density=False, bins=RMiss_1tau_nucleon_p_50GeV_bins, weights=wts_RMiss_1tau_nucleon_p_50GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax44.hist(RMiss_1tau_nucleon_p_50GeV, density=False, bins=RMiss_1tau_nucleon_p_50GeV_bins, weights=wts_RMiss_1tau_nucleon_p_50GeV, histtype='step', color='black', alpha=1, lw=2)
    ax45.hist(RMiss_1tau_nucleon_p_200GeV, density=False, bins=RMiss_1tau_nucleon_p_200GeV_bins, weights=wts_RMiss_1tau_nucleon_p_200GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax45.hist(RMiss_1tau_nucleon_p_200GeV, density=False, bins=RMiss_1tau_nucleon_p_200GeV_bins, weights=wts_RMiss_1tau_nucleon_p_200GeV, histtype='step', color='black', alpha=1, lw=2)
    
    ### Fermi4Q2 ###
    ax51.hist(Fermi4Q2_1tau_nucleon_p_5GeV, density=False, bins=Fermi4Q2_1tau_nucleon_p_5GeV_bins, weights=wts_Fermi4Q2_1tau_nucleon_p_5GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax51.hist(Fermi4Q2_1tau_nucleon_p_5GeV, density=False, bins=Fermi4Q2_1tau_nucleon_p_5GeV_bins, weights=wts_Fermi4Q2_1tau_nucleon_p_5GeV, histtype='step', color='black', alpha=1, lw=2)
    ax52.hist(Fermi4Q2_1tau_nucleon_p_10GeV, density=False, bins=Fermi4Q2_1tau_nucleon_p_10GeV_bins, weights=wts_Fermi4Q2_1tau_nucleon_p_10GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax52.hist(Fermi4Q2_1tau_nucleon_p_10GeV, density=False, bins=Fermi4Q2_1tau_nucleon_p_10GeV_bins, weights=wts_Fermi4Q2_1tau_nucleon_p_10GeV, histtype='step', color='black', alpha=1, lw=2)
    ax53.hist(Fermi4Q2_1tau_nucleon_p_20GeV, density=False, bins=Fermi4Q2_1tau_nucleon_p_20GeV_bins, weights=wts_Fermi4Q2_1tau_nucleon_p_20GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax53.hist(Fermi4Q2_1tau_nucleon_p_20GeV, density=False, bins=Fermi4Q2_1tau_nucleon_p_20GeV_bins, weights=wts_Fermi4Q2_1tau_nucleon_p_20GeV, histtype='step', color='black', alpha=1, lw=2)
    ax54.hist(Fermi4Q2_1tau_nucleon_p_50GeV, density=False, bins=Fermi4Q2_1tau_nucleon_p_50GeV_bins, weights=wts_Fermi4Q2_1tau_nucleon_p_50GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax54.hist(Fermi4Q2_1tau_nucleon_p_50GeV, density=False, bins=Fermi4Q2_1tau_nucleon_p_50GeV_bins, weights=wts_Fermi4Q2_1tau_nucleon_p_50GeV, histtype='step', color='black', alpha=1, lw=2)
    ax55.hist(Fermi4Q2_1tau_nucleon_p_200GeV, density=False, bins=Fermi4Q2_1tau_nucleon_p_200GeV_bins, weights=wts_Fermi4Q2_1tau_nucleon_p_200GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax55.hist(Fermi4Q2_1tau_nucleon_p_200GeV, density=False, bins=Fermi4Q2_1tau_nucleon_p_200GeV_bins, weights=wts_Fermi4Q2_1tau_nucleon_p_200GeV, histtype='step', color='black', alpha=1, lw=2)

    ### ThetaMiss ###
    ax61.hist(ThetaMiss_1tau_nucleon_p_5GeV, density=False, bins=ThetaMiss_1tau_nucleon_p_5GeV_bins, weights=wts_ThetaMiss_1tau_nucleon_p_5GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax61.hist(ThetaMiss_1tau_nucleon_p_5GeV, density=False, bins=ThetaMiss_1tau_nucleon_p_5GeV_bins, weights=wts_ThetaMiss_1tau_nucleon_p_5GeV, histtype='step', color='black', alpha=1, lw=2)
    ax62.hist(ThetaMiss_1tau_nucleon_p_10GeV, density=False, bins=ThetaMiss_1tau_nucleon_p_10GeV_bins, weights=wts_ThetaMiss_1tau_nucleon_p_10GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax62.hist(ThetaMiss_1tau_nucleon_p_10GeV, density=False, bins=ThetaMiss_1tau_nucleon_p_10GeV_bins, weights=wts_ThetaMiss_1tau_nucleon_p_10GeV, histtype='step', color='black', alpha=1, lw=2)
    ax63.hist(ThetaMiss_1tau_nucleon_p_20GeV, density=False, bins=ThetaMiss_1tau_nucleon_p_20GeV_bins, weights=wts_ThetaMiss_1tau_nucleon_p_20GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax63.hist(ThetaMiss_1tau_nucleon_p_20GeV, density=False, bins=ThetaMiss_1tau_nucleon_p_20GeV_bins, weights=wts_ThetaMiss_1tau_nucleon_p_20GeV, histtype='step', color='black', alpha=1, lw=2)
    ax64.hist(ThetaMiss_1tau_nucleon_p_50GeV, density=False, bins=ThetaMiss_1tau_nucleon_p_50GeV_bins, weights=wts_ThetaMiss_1tau_nucleon_p_50GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax64.hist(ThetaMiss_1tau_nucleon_p_50GeV, density=False, bins=ThetaMiss_1tau_nucleon_p_50GeV_bins, weights=wts_ThetaMiss_1tau_nucleon_p_50GeV, histtype='step', color='black', alpha=1, lw=2)
    ax65.hist(ThetaMiss_1tau_nucleon_p_200GeV, density=False, bins=ThetaMiss_1tau_nucleon_p_200GeV_bins, weights=wts_ThetaMiss_1tau_nucleon_p_200GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax65.hist(ThetaMiss_1tau_nucleon_p_200GeV, density=False, bins=ThetaMiss_1tau_nucleon_p_200GeV_bins, weights=wts_ThetaMiss_1tau_nucleon_p_200GeV, histtype='step', color='black', alpha=1, lw=2)

    ### Labels ###
    ax11.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / bin / N)', fontsize=40)
    ax21.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / bin / N)', fontsize=40)
    ax31.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / bin / N)', fontsize=40)
    ax41.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / bin / N)', fontsize=40)
    ax51.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / bin / N)', fontsize=40)
    ax61.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / bin / N)', fontsize=40)
  
    for ax in [ax11, ax12, ax13, ax14, ax15]:
        ax.set_xlabel(r'\textbf{Missing Transverse Momentum} $p^T_{\mathrm{Miss}}$ (GeV) ', fontsize=40)

    for ax in [ax21, ax22, ax23, ax24, ax25]:
        ax.set_xlabel(r'\textbf{Muon Transverse Momentum} $p^T_\mu$ (GeV)', fontsize=40)

    for ax in [ax31, ax32, ax33, ax34, ax35]:
        ax.set_xlabel(r'\textbf{Tau Transverse Momentum} $p^T_\tau$ (GeV)', fontsize=40)
    
    for ax in [ax41, ax42, ax43, ax44, ax45]:
        ax.set_xlabel(r'\textbf{Missing Momentum Ratio} $R^T_{\mathrm{Miss}}$', fontsize=40)
    
    for ax in [ax51, ax52, ax53, ax54, ax55]:
        ax.set_xlabel(r'\textbf{4-Fermi Momentum Transfer} $Q^2_{\mathrm{4-Fermi}}$', fontsize=40)
        
    for ax in [ax61, ax62, ax63, ax64, ax65]:
        ax.set_xlabel(r'\textbf{Missing Transverse Polar Angle} $\theta^T_{\mathrm{Miss}}$', fontsize=40)

    for ax in [ax11, ax21, ax31, ax61, ax51]:
        ax.text(0.90,0.90,r'\textbf{5 GeV}',ha='right',transform=ax.transAxes, fontsize=40)

    for ax in [ax12, ax22, ax32, ax62, ax52]:
        ax.text(0.90,0.90,r'\textbf{10 GeV}',ha='right',transform=ax.transAxes, fontsize=40)

    for ax in [ax13, ax23, ax33, ax63, ax53]:
        ax.text(0.90,0.90,r'\textbf{20 GeV}',ha='right',transform=ax.transAxes, fontsize=40)

    for ax in [ax14, ax24, ax34, ax64, ax54]:
        ax.text(0.90,0.90,r'\textbf{50 GeV}',ha='right',transform=ax.transAxes, fontsize=40)

    for ax in [ax15, ax25, ax35, ax65, ax55]:
        ax.text(0.90,0.90,r'\textbf{200 GeV}',ha='right',transform=ax.transAxes, fontsize=40)

    ax41.text(0.10,0.90,r'\textbf{5 GeV}',ha='left',transform=ax41.transAxes,fontsize=40)
    ax42.text(0.10,0.90,r'\textbf{10 GeV}',ha='left',transform=ax42.transAxes,fontsize=40)
    ax43.text(0.10,0.90,r'\textbf{20 GeV}',ha='left',transform=ax43.transAxes,fontsize=40)
    ax44.text(0.10,0.90,r'\textbf{50 GeV}',ha='left',transform=ax44.transAxes,fontsize=40)
    ax45.text(0.10,0.90,r'\textbf{200 GeV}',ha='left',transform=ax45.transAxes,fontsize=40)

    ### Limits ###
    # ptMiss #
    ax11.set_xlim(ptMiss_1tau_nucleon_p_5GeV_min, ptMiss_1tau_nucleon_p_5GeV_max)
    ax12.set_xlim(ptMiss_1tau_nucleon_p_10GeV_min, ptMiss_1tau_nucleon_p_10GeV_max)
    ax13.set_xlim(ptMiss_1tau_nucleon_p_20GeV_min, ptMiss_1tau_nucleon_p_20GeV_max)
    ax14.set_xlim(ptMiss_1tau_nucleon_p_50GeV_min, ptMiss_1tau_nucleon_p_50GeV_max)
    ax15.set_xlim(ptMiss_1tau_nucleon_p_200GeV_min, ptMiss_1tau_nucleon_p_200GeV_max)

    # ptMuon #
    ax21.set_xlim(ptMuon_1tau_nucleon_p_5GeV_min, ptMuon_1tau_nucleon_p_5GeV_max)
    ax22.set_xlim(ptMuon_1tau_nucleon_p_10GeV_min, ptMuon_1tau_nucleon_p_10GeV_max)
    ax23.set_xlim(ptMuon_1tau_nucleon_p_20GeV_min, ptMuon_1tau_nucleon_p_20GeV_max)
    ax24.set_xlim(ptMuon_1tau_nucleon_p_50GeV_min, ptMuon_1tau_nucleon_p_50GeV_max)
    ax25.set_xlim(ptMuon_1tau_nucleon_p_200GeV_min, ptMuon_1tau_nucleon_p_200GeV_max)
    
    # ptTau #
    ax31.set_xlim(ptTau_1tau_nucleon_p_5GeV_min, ptTau_1tau_nucleon_p_5GeV_max)
    ax32.set_xlim(ptTau_1tau_nucleon_p_10GeV_min, ptTau_1tau_nucleon_p_10GeV_max)
    ax33.set_xlim(ptTau_1tau_nucleon_p_20GeV_min, ptTau_1tau_nucleon_p_20GeV_max)
    ax34.set_xlim(ptTau_1tau_nucleon_p_50GeV_min, ptTau_1tau_nucleon_p_50GeV_max)
    ax35.set_xlim(ptTau_1tau_nucleon_p_200GeV_min, ptTau_1tau_nucleon_p_200GeV_max)
    
    ax41.set_xlim(0, 1)
    ax42.set_xlim(0, 1)
    ax43.set_xlim(0, 1)
    ax44.set_xlim(0, 1)
    ax45.set_xlim(0, 1)
   
    # Fermi4Q2 #
    ax51.set_xlim(Fermi4Q2_1tau_nucleon_p_5GeV_min, Fermi4Q2_1tau_nucleon_p_5GeV_max)
    ax52.set_xlim(Fermi4Q2_1tau_nucleon_p_10GeV_min, Fermi4Q2_1tau_nucleon_p_10GeV_max)
    ax53.set_xlim(Fermi4Q2_1tau_nucleon_p_20GeV_min, Fermi4Q2_1tau_nucleon_p_20GeV_max)
    ax54.set_xlim(Fermi4Q2_1tau_nucleon_p_50GeV_min, Fermi4Q2_1tau_nucleon_p_50GeV_max)
    ax55.set_xlim(Fermi4Q2_1tau_nucleon_p_200GeV_min, Fermi4Q2_1tau_nucleon_p_200GeV_max)
    
    # ThetaMiss #
    ax61.set_xlim(0, 180)
    ax62.set_xlim(0, 180)
    ax63.set_xlim(0, 180)
    ax64.set_xlim(0, 180)
    ax65.set_xlim(0, 180)

    # Tick labels #
    xmajor_RMiss = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    xmajor_ThetaMiss = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    
    ax41.set_xticks(xmajor_RMiss, labels=['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'])
    ax41.tick_params(which='both',right=False, bottom=True)
    ax42.set_xticks(xmajor_RMiss, labels=['', '0.2', '0.4', '0.6', '0.8', '1.0'])
    ax42.tick_params(which='both',right=False, bottom=True)
    ax43.set_xticks(xmajor_RMiss, labels=['', '0.2', '0.4', '0.6', '0.8', '1.0'])
    ax43.tick_params(which='both',right=False, bottom=True)
    ax44.set_xticks(xmajor_RMiss, labels=['', '0.2', '0.4', '0.6', '0.8', '1.0'])
    ax44.tick_params(which='both',right=False, bottom=True)
    ax45.set_xticks(xmajor_RMiss, labels=['', '0.2', '0.4', '0.6', '0.8', '1.0'])
    ax45.tick_params(which='both',right=False, bottom=True)

    ax61.set_xticks(xmajor_ThetaMiss, labels=['0', '20', '40', '60', '80', '100', '120', '140', '160', '180'])
    ax61.tick_params(which='both',right=False, bottom=True)
    ax62.set_xticks(xmajor_ThetaMiss, labels=['', '20', '40', '60', '80', '100', '120', '140', '160', '180'])
    ax62.tick_params(which='both',right=False, bottom=True)
    ax63.set_xticks(xmajor_ThetaMiss, labels=['', '20', '40', '60', '80', '100', '120', '140', '160', '180'])
    ax63.tick_params(which='both',right=False, bottom=True)
    ax64.set_xticks(xmajor_ThetaMiss, labels=['', '20', '40', '60', '80', '100', '120', '140', '160', '180'])
    ax64.tick_params(which='both',right=False, bottom=True)
    ax65.set_xticks(xmajor_ThetaMiss, labels=['', '20', '40', '60', '80', '100', '120', '140', '160', '180'])
    ax65.tick_params(which='both',right=False, bottom=True)

# horizontal #
if horizontal_flag:
    fig, ax = plt.subplots(2, 4, figsize=(50, 26), tight_layout=True)

    ### vmu -> vtau tau+ mu- ###
    # Plot Ar for 10, 20, 50 and 200 GeV #
    ax[0,0].hist(Q_1tau_nucleon_p_20GeV, density=False, bins=Q_1tau_nucleon_p_20GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[0,0].hist(Q_1tau_nucleon_p_20GeV, density=False, bins=Q_1tau_nucleon_p_20GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[0,0].text(0.935,0.90,r'$E_\nu=$\textbf{ 20 GeV}',ha='right',fontsize=30, transform=ax[0,0].transAxes)

    ax[0,1].hist(Q_1tau_nucleon_p_50GeV, density=False, bins=Q_1tau_nucleon_p_50GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[0,1].hist(Q_1tau_nucleon_p_50GeV, density=False, bins=Q_1tau_nucleon_p_50GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[0,1].text(0.935,0.90,r'$E_\nu=$\textbf{ 50 GeV}',ha='right',fontsize=30, transform=ax[0,1].transAxes)

    ax[0,2].hist(Q_1tau_nucleon_p_200GeV, density=False, bins=Q_1tau_nucleon_p_200GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[0,2].hist(Q_1tau_nucleon_p_200GeV, density=False, bins=Q_1tau_nucleon_p_200GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[0,2].text(0.935,0.90,r'$E_\nu=$\textbf{ 200 GeV}',ha='right',fontsize=30, transform=ax[0,2].transAxes)

    ax[0,3].hist(Q_1tau_nucleon_p_1000GeV, density=False, bins=Q_1tau_nucleon_p_1000GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[0,3].hist(Q_1tau_nucleon_p_1000GeV, density=False, bins=Q_1tau_nucleon_p_1000GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[0,3].text(0.935,0.90,r'$E_\nu=$\textbf{ 1000 GeV}',ha='right',fontsize=30, transform=ax[0,3].transAxes)

    ### vmu -> vmu tau+ tau- ###
    # Plot Ar for 10, 20, 50 and 200 GeV #
    ax[1,0].hist(Q_2tau_nucleon_p_20GeV, density=False, bins=Q_2tau_nucleon_p_20GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[1,0].hist(Q_2tau_nucleon_p_20GeV, density=False, bins=Q_2tau_nucleon_p_20GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[1,0].text(0.935,0.90,r'$E_\nu=$\textbf{ 20 GeV}',ha='right',fontsize=30, transform=ax[1,0].transAxes)

    ax[1,1].hist(Q_2tau_nucleon_p_50GeV, density=False, bins=Q_2tau_nucleon_p_50GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[1,1].hist(Q_2tau_nucleon_p_50GeV, density=False, bins=Q_2tau_nucleon_p_50GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[1,1].text(0.935,0.90,r'$E_\nu=$\textbf{ 50 GeV}',ha='right',fontsize=30, transform=ax[1,1].transAxes)

    ax[1,2].hist(Q_2tau_nucleon_p_200GeV, density=False, bins=Q_2tau_nucleon_p_200GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[1,2].hist(Q_2tau_nucleon_p_200GeV, density=False, bins=Q_2tau_nucleon_p_200GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[1,2].text(0.935,0.90,r'$E_\nu=$\textbf{ 200 GeV}',ha='right',fontsize=30, transform=ax[1,2].transAxes)

    ax[1,3].hist(Q_2tau_nucleon_p_1000GeV, density=False, bins=Q_2tau_nucleon_p_1000GeV_bins, color='green',alpha=0.3, edgecolor='black', lw=0.5)
    ax[1,3].hist(Q_2tau_nucleon_p_1000GeV, density=False, bins=Q_2tau_nucleon_p_1000GeV_bins, histtype='step', color='black',lw=2, alpha=1)
    ax[1,3].text(0.935,0.90,r'$E_\nu=$\textbf{ 1000 GeV}',ha='right',fontsize=30, transform=ax[1,3].transAxes)

 
    ### Styling and save ###
    # Axis labels #
    ax[0,0].set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$', fontsize=60)

    ax[1,0].set_ylabel(r'$\nu_\mu \to \nu_\mu \tau^+ \tau^-$', fontsize=60)
    
    locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
    locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)

    xmin_values = [Q_1tau_nucleon_p_20GeV_min, Q_1tau_nucleon_p_50GeV_min, Q_1tau_nucleon_p_200GeV_min, Q_1tau_nucleon_p_1000GeV_min,
                   Q_2tau_nucleon_p_20GeV_min, Q_2tau_nucleon_p_50GeV_min, Q_2tau_nucleon_p_200GeV_min, Q_2tau_nucleon_p_1000GeV_min]

    k = 0
    for i in range(2):
        for j in range(4):        
            xmin = xmin_values[k]
            print(i, j, k, xmin)
            ax[i,j].xaxis.set_major_locator(locmaj)
            ax[i,j].xaxis.set_minor_locator(locmin)
            ax[i,j].xaxis.set_major_locator(locmaj)
            ax[i,j].xaxis.set_minor_locator(locmin)
            ax[i,j].set_ylim(1,1e7)
            ax[i,j].set_xlim(xmin,20)
            
            ax[i,j].set_xscale('log')
            ax[i,j].set_yscale('log')
            k += 1
# Save #
fig1.savefig("../plots/ptMiss_1tau_incoh_p.png", dpi=100, bbox_inches='tight')
fig2.savefig("../plots/ptMuon_1tau_incoh_p.png", dpi=100, bbox_inches='tight')
fig3.savefig("../plots/ptTau_1tau_incoh_p.png", dpi=100, bbox_inches='tight')
fig4.savefig("../plots/RMiss_1tau_incoh_p.png", dpi=100, bbox_inches='tight')
fig5.savefig("../plots/Fermi4Q2_1tau_incoh_p.png", dpi=100, bbox_inches='tight')
fig6.savefig("../plots/ThetaMiss_1tau_incoh_p.png", dpi=100, bbox_inches='tight')
