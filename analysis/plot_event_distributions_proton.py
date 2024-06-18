import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import csv

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

EVENTS_DIR_1TAU = '../csv/events/vmu_to_vtau_tau+_mu-/nucleon/proton'

vertical_flag = True
horizontal_flag = False
write_dist = False

##########################
#### Q2 Distributions ####
##########################

###### Nucleon Proton ######
### 1tau ; proton ; 5 GeV ###
ptMiss_1tau_nucleon_p_5GeV = []
ptMuon_1tau_nucleon_p_5GeV = []
ptTau_1tau_nucleon_p_5GeV = []
RMiss_1tau_nucleon_p_5GeV = []

### 1tau ; proton ; 10 GeV ###
ptMiss_1tau_nucleon_p_10GeV = []
ptMuon_1tau_nucleon_p_10GeV = []
ptTau_1tau_nucleon_p_10GeV = []
RMiss_1tau_nucleon_p_10GeV = []

### 1tau ; proton ; 20 GeV ###
ptMiss_1tau_nucleon_p_20GeV = []
ptMuon_1tau_nucleon_p_20GeV = []
ptTau_1tau_nucleon_p_20GeV = []
RMiss_1tau_nucleon_p_20GeV = []

### 1tau ; proton ; 50 GeV ###
ptMiss_1tau_nucleon_p_50GeV = []
ptMuon_1tau_nucleon_p_50GeV = []
ptTau_1tau_nucleon_p_50GeV = []
RMiss_1tau_nucleon_p_50GeV = []

### 1tau ; proton ; 200 GeV ###
ptMiss_1tau_nucleon_p_200GeV = []
ptMuon_1tau_nucleon_p_200GeV = []
ptTau_1tau_nucleon_p_200GeV = []
RMiss_1tau_nucleon_p_200GeV = []


############################
#### Read Distributions ####
############################

def read_after_skip(file_path, skip_lines=37):
    with open(file_path, 'r') as txtfile:
        for _ in range(skip_lines):
            next(txtfile)
        pt_miss_array = []
        pt_muon_array = []
        pt_tau_array = []
        for line in txtfile:
            if line.startswith('16'):  # Read and store vtau transverse momentum
                data = line.split()
                assert int(data[0]) ==  16
                px = float(data[2])
                py = float(data[3])
                pt_miss = np.sqrt(px**2 + py**2)
                pt_miss_array.append(pt_miss)
            if line.startswith('13'):  # Read and store muon transverse momentum
                data = line.split()
                assert int(data[0]) == 13
                px = float(data[2])
                py = float(data[3])
                pt_muon = np.sqrt(px**2 + py**2)
                pt_muon_array.append(pt_muon)
            if line.startswith('-15'):  # Read and store tau transverse momentum
                data = line.split()
                assert int(data[0]) == -15
                px = float(data[2])
                py = float(data[3])
                pt_tau = np.sqrt(px**2 + py**2)
                pt_tau_array.append(pt_tau)
        r_miss_array = [pt_miss / (pt_miss + pt_muon + pt_tau) for pt_miss, pt_muon, pt_tau in zip(pt_miss_array, pt_muon_array, pt_tau_array)]  # Calculate missing momentum ratio
        return pt_miss_array, pt_muon_array, pt_tau_array, r_miss_array

ptMiss_1tau_nucleon_p_5GeV, ptMuon_1tau_nucleon_p_5GeV, ptTau_1tau_nucleon_p_5GeV, RMiss_1tau_nucleon_p_5GeV = read_after_skip(EVENTS_DIR_1TAU+'/1tau_1e5Events_p_5GeV.txt')
ptMiss_1tau_nucleon_p_10GeV, ptMuon_1tau_nucleon_p_10GeV, ptTau_1tau_nucleon_p_10GeV, RMiss_1tau_nucleon_p_10GeV = read_after_skip(EVENTS_DIR_1TAU+'/1tau_1e5Events_p_10GeV.txt')
ptMiss_1tau_nucleon_p_20GeV, ptMuon_1tau_nucleon_p_20GeV, ptTau_1tau_nucleon_p_20GeV, RMiss_1tau_nucleon_p_20GeV = read_after_skip(EVENTS_DIR_1TAU+'/1tau_1e5Events_p_20GeV.txt')
ptMiss_1tau_nucleon_p_50GeV, ptMuon_1tau_nucleon_p_50GeV, ptTau_1tau_nucleon_p_50GeV, RMiss_1tau_nucleon_p_50GeV = read_after_skip(EVENTS_DIR_1TAU+'/1tau_1e5Events_p_50GeV.txt')
ptMiss_1tau_nucleon_p_200GeV, ptMuon_1tau_nucleon_p_200GeV, ptTau_1tau_nucleon_p_200GeV, RMiss_1tau_nucleon_p_200GeV = read_after_skip(EVENTS_DIR_1TAU+'/1tau_1e5Events_p_200GeV.txt')

if write_dist:
    with open('ptDist_1tau_nucleon_p_5GeV.txt','w',newline='') as outfile:
        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for pt_miss, pt_muon, pt_tau, r_miss in zip(ptMiss_1tau_nucleon_p_5GeV, ptMuon_1tau_nucleon_p_5GeV, ptTau_1tau_nucleon_p_5GeV, RMiss_1tau_nucleon_p_5GeV):
            writer.writerow([pt_miss, pt_muon, pt_tau, r_miss])

    with open('ptDist_1tau_nucleon_p_10GeV.txt','w',newline='') as outfile:
        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for pt_miss, pt_muon, pt_tau, r_miss in zip(ptMiss_1tau_nucleon_p_10GeV, ptMuon_1tau_nucleon_p_10GeV, ptTau_1tau_nucleon_p_10GeV, RMiss_1tau_nucleon_p_10GeV):
            writer.writerow([pt_miss, pt_muon, pt_tau, r_miss])

    with open('ptDist_1tau_nucleon_p_20GeV.txt','w',newline='') as outfile:
        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for pt_miss, pt_muon, pt_tau, r_miss in zip(ptMiss_1tau_nucleon_p_20GeV, ptMuon_1tau_nucleon_p_20GeV, ptTau_1tau_nucleon_p_20GeV, RMiss_1tau_nucleon_p_20GeV):
            writer.writerow([pt_miss, pt_muon, pt_tau, r_miss])

    with open('ptDist_1tau_nucleon_p_50GeV.txt','w',newline='') as outfile:
        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for pt_miss, pt_muon, pt_tau, r_miss in zip(ptMiss_1tau_nucleon_p_50GeV, ptMuon_1tau_nucleon_p_50GeV, ptTau_1tau_nucleon_p_50GeV, RMiss_1tau_nucleon_p_50GeV):
            writer.writerow([pt_miss, pt_muon, pt_tau, r_miss])

    with open('ptDist_1tau_nucleon_p_200GeV.txt','w',newline='') as outfile:
        writer = csv.writer(outfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for pt_miss, pt_muon, pt_tau, r_miss in zip(ptMiss_1tau_nucleon_p_200GeV, ptMuon_1tau_nucleon_p_200GeV, ptTau_1tau_nucleon_p_200GeV, RMiss_1tau_nucleon_p_200GeV):
            writer.writerow([pt_miss, pt_muon, pt_tau, r_miss])


##################
#### Binning ####
##################
NBINS = 100

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
RMiss_1tau_nucleon_p_5GeV_bins, RMiss_1tau_nucleon_p_5GeV_step = np.linspace(RMiss_1tau_nucleon_p_5GeV_min, RMiss_1tau_nucleon_p_5GeV_max, num=NBINS, retstep=True)
print("RMiss:")
print("\tmin: ", RMiss_1tau_nucleon_p_5GeV_min)
print("\tmax: ", RMiss_1tau_nucleon_p_5GeV_max)

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
RMiss_1tau_nucleon_p_10GeV_bins, RMiss_1tau_nucleon_p_10GeV_step = np.linspace(RMiss_1tau_nucleon_p_10GeV_min, RMiss_1tau_nucleon_p_10GeV_max, num=NBINS, retstep=True)
print("RMiss:")
print("\tmin: ", RMiss_1tau_nucleon_p_10GeV_min)
print("\tmax: ", RMiss_1tau_nucleon_p_10GeV_max)

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
RMiss_1tau_nucleon_p_20GeV_bins, RMiss_1tau_nucleon_p_20GeV_step = np.linspace(RMiss_1tau_nucleon_p_20GeV_min, RMiss_1tau_nucleon_p_20GeV_max, num=NBINS, retstep=True)
print("RMiss:")
print("\tmin: ", RMiss_1tau_nucleon_p_20GeV_min)
print("\tmax: ", RMiss_1tau_nucleon_p_20GeV_max)

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
RMiss_1tau_nucleon_p_50GeV_bins, RMiss_1tau_nucleon_p_50GeV_step = np.linspace(RMiss_1tau_nucleon_p_50GeV_min, RMiss_1tau_nucleon_p_50GeV_max, num=NBINS, retstep=True)
print("RMiss:")
print("\tmin: ", RMiss_1tau_nucleon_p_50GeV_min)
print("\tmax: ", RMiss_1tau_nucleon_p_50GeV_max)

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
RMiss_1tau_nucleon_p_200GeV_bins, RMiss_1tau_nucleon_p_200GeV_step = np.linspace(RMiss_1tau_nucleon_p_200GeV_min, RMiss_1tau_nucleon_p_200GeV_max, num=NBINS, retstep=True)
print("RMiss:")
print("\tmin: ", RMiss_1tau_nucleon_p_200GeV_min)
print("\tmax: ", RMiss_1tau_nucleon_p_200GeV_max)

### Weights ###
NEVENTS = 1e5

# 5 GeV #
wts_ptMiss_1tau_nucleon_p_5GeV = np.ones_like(ptMiss_1tau_nucleon_p_5GeV) / (NEVENTS * ptMiss_1tau_nucleon_p_5GeV_step)
wts_ptMuon_1tau_nucleon_p_5GeV = np.ones_like(ptMuon_1tau_nucleon_p_5GeV) / (NEVENTS * ptMuon_1tau_nucleon_p_5GeV_step)
wts_ptTau_1tau_nucleon_p_5GeV = np.ones_like(ptTau_1tau_nucleon_p_5GeV) / (NEVENTS * ptTau_1tau_nucleon_p_5GeV_step)
wts_RMiss_1tau_nucleon_p_5GeV = np.ones_like(RMiss_1tau_nucleon_p_5GeV) / (NEVENTS * RMiss_1tau_nucleon_p_5GeV_step)

# 10 GeV #
wts_ptMiss_1tau_nucleon_p_10GeV = np.ones_like(ptMiss_1tau_nucleon_p_10GeV) / (NEVENTS * ptMiss_1tau_nucleon_p_10GeV_step)
wts_ptMuon_1tau_nucleon_p_10GeV = np.ones_like(ptMuon_1tau_nucleon_p_10GeV) / (NEVENTS * ptMuon_1tau_nucleon_p_10GeV_step)
wts_ptTau_1tau_nucleon_p_10GeV = np.ones_like(ptTau_1tau_nucleon_p_10GeV) / (NEVENTS * ptTau_1tau_nucleon_p_10GeV_step)
wts_RMiss_1tau_nucleon_p_10GeV = np.ones_like(RMiss_1tau_nucleon_p_10GeV) / (NEVENTS * RMiss_1tau_nucleon_p_10GeV_step)

# 20 GeV #
wts_ptMiss_1tau_nucleon_p_20GeV = np.ones_like(ptMiss_1tau_nucleon_p_20GeV) / (NEVENTS * ptMiss_1tau_nucleon_p_20GeV_step)
wts_ptMuon_1tau_nucleon_p_20GeV = np.ones_like(ptMuon_1tau_nucleon_p_20GeV) / (NEVENTS * ptMuon_1tau_nucleon_p_20GeV_step)
wts_ptTau_1tau_nucleon_p_20GeV = np.ones_like(ptTau_1tau_nucleon_p_20GeV) / (NEVENTS * ptTau_1tau_nucleon_p_20GeV_step)
wts_RMiss_1tau_nucleon_p_20GeV = np.ones_like(RMiss_1tau_nucleon_p_20GeV) / (NEVENTS * RMiss_1tau_nucleon_p_20GeV_step)

# 50 GeV #
wts_ptMiss_1tau_nucleon_p_50GeV = np.ones_like(ptMiss_1tau_nucleon_p_50GeV) / (NEVENTS * ptMiss_1tau_nucleon_p_50GeV_step)
wts_ptMuon_1tau_nucleon_p_50GeV = np.ones_like(ptMuon_1tau_nucleon_p_50GeV) / (NEVENTS * ptMuon_1tau_nucleon_p_50GeV_step)
wts_ptTau_1tau_nucleon_p_50GeV = np.ones_like(ptTau_1tau_nucleon_p_50GeV) / (NEVENTS * ptTau_1tau_nucleon_p_50GeV_step)
wts_RMiss_1tau_nucleon_p_50GeV = np.ones_like(RMiss_1tau_nucleon_p_50GeV) / (NEVENTS * RMiss_1tau_nucleon_p_50GeV_step)

# 200 GeV #
wts_ptMiss_1tau_nucleon_p_200GeV = np.ones_like(ptMiss_1tau_nucleon_p_200GeV) / (NEVENTS * ptMiss_1tau_nucleon_p_200GeV_step)
wts_ptMuon_1tau_nucleon_p_200GeV = np.ones_like(ptMuon_1tau_nucleon_p_200GeV) / (NEVENTS * ptMuon_1tau_nucleon_p_200GeV_step)
wts_ptTau_1tau_nucleon_p_200GeV = np.ones_like(ptTau_1tau_nucleon_p_200GeV) / (NEVENTS * ptTau_1tau_nucleon_p_200GeV_step)
wts_RMiss_1tau_nucleon_p_200GeV = np.ones_like(RMiss_1tau_nucleon_p_200GeV) / (NEVENTS * RMiss_1tau_nucleon_p_200GeV_step)

##################
#### Plotting ####
##################

# vertical #
if vertical_flag:
    #fig, ax = plt.subplots(4, 2, figsize=(28, 36), tight_layout=True)
    fig1, (ax11, ax12, ax13, ax14, ax15) = plt.subplots(1, 5, sharey=True, figsize=(60,15), tight_layout=True) # ptMiss
    fig2, (ax21, ax22, ax23, ax24, ax25) = plt.subplots(1, 5, sharey=True, figsize=(60,15), tight_layout=True) # ptMuon
    fig3, (ax31, ax32, ax33, ax34, ax35) = plt.subplots(1, 5, sharey=True, figsize=(60,15), tight_layout=True) # ptTau
    fig4, (ax41, ax42, ax43, ax44, ax45) = plt.subplots(1, 5, sharey=True, figsize=(60,15), tight_layout=True) # RMiss

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
    ax41.hist(RMiss_1tau_nucleon_p_5GeV, density=False, bins=RMiss_1tau_nucleon_p_10GeV_bins, weights=wts_RMiss_1tau_nucleon_p_10GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax41.hist(RMiss_1tau_nucleon_p_5GeV, density=False, bins=RMiss_1tau_nucleon_p_10GeV_bins, weights=wts_RMiss_1tau_nucleon_p_10GeV, histtype='step', color='black', alpha=1, lw=2)
    ax42.hist(RMiss_1tau_nucleon_p_10GeV, density=False, bins=RMiss_1tau_nucleon_p_10GeV_bins, weights=wts_RMiss_1tau_nucleon_p_10GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax42.hist(RMiss_1tau_nucleon_p_10GeV, density=False, bins=RMiss_1tau_nucleon_p_10GeV_bins, weights=wts_RMiss_1tau_nucleon_p_10GeV, histtype='step', color='black', alpha=1, lw=2)
    ax43.hist(RMiss_1tau_nucleon_p_20GeV, density=False, bins=RMiss_1tau_nucleon_p_20GeV_bins, weights=wts_RMiss_1tau_nucleon_p_20GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax43.hist(RMiss_1tau_nucleon_p_20GeV, density=False, bins=RMiss_1tau_nucleon_p_20GeV_bins, weights=wts_RMiss_1tau_nucleon_p_20GeV, histtype='step', color='black', alpha=1, lw=2)
    ax44.hist(RMiss_1tau_nucleon_p_50GeV, density=False, bins=RMiss_1tau_nucleon_p_50GeV_bins, weights=wts_RMiss_1tau_nucleon_p_50GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax44.hist(RMiss_1tau_nucleon_p_50GeV, density=False, bins=RMiss_1tau_nucleon_p_50GeV_bins, weights=wts_RMiss_1tau_nucleon_p_50GeV, histtype='step', color='black', alpha=1, lw=2)
    ax45.hist(RMiss_1tau_nucleon_p_200GeV, density=False, bins=RMiss_1tau_nucleon_p_200GeV_bins, weights=wts_RMiss_1tau_nucleon_p_200GeV, color='royalblue', alpha=0.3, edgecolor='black', lw=0.5)
    ax45.hist(RMiss_1tau_nucleon_p_200GeV, density=False, bins=RMiss_1tau_nucleon_p_200GeV_bins, weights=wts_RMiss_1tau_nucleon_p_200GeV, histtype='step', color='black', alpha=1, lw=2)

    ### Labels ###
    ax11.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / bin / N)', fontsize=40)
    ax21.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / bin / N)', fontsize=40)
    ax31.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / bin / N)', fontsize=40)
    ax41.set_ylabel(r'$\nu_\mu \to \nu_\tau \tau^+ \mu^-$ (dN / bin / N)', fontsize=40)
    
    for ax in [ax11, ax12, ax13, ax14, ax15]:
        ax.set_xlabel(r'\textbf{Missing Transverse Momentum} $p^T_{\mathrm{Miss}}$ (GeV) ', fontsize=40)

    for ax in [ax21, ax22, ax23, ax24, ax25]:
        ax.set_xlabel(r'\textbf{Muon Transverse Momentum} $p^T_\mu$ (GeV)', fontsize=40)

    for ax in [ax31, ax32, ax33, ax34, ax35]:
        ax.set_xlabel(r'\textbf{Tau Transverse Momentum} $p^T_\tau$ (GeV)', fontsize=40)
    
    for ax in [ax41, ax42, ax43, ax44, ax45]:
        ax.set_xlabel(r'\textbf{Missing Momentum Ratio} $R^T_{\mathrm{Miss}}$', fontsize=40)

    for ax in [ax11, ax21, ax31, ax41]:
        ax.text(0.90,0.90,r'\textbf{5 GeV}',ha='right',transform=ax.transAxes)

    for ax in [ax12, ax22, ax32, ax42]:
        ax.text(0.90,0.90,r'\textbf{10 GeV}',ha='right',transform=ax.transAxes)

    for ax in [ax13, ax23, ax33, ax43]:
        ax.text(0.90,0.90,r'\textbf{20 GeV}',ha='right',transform=ax.transAxes)

    for ax in [ax14, ax24, ax34, ax44]:
        ax.text(0.90,0.90,r'\textbf{50 GeV}',ha='right',transform=ax.transAxes)

    for ax in [ax15, ax25, ax35, ax45]:
        ax.text(0.90,0.90,r'\textbf{200 GeV}',ha='right',transform=ax.transAxes)

    ### Ticks ###
    locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
    locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)
    locmaj2 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
    locmin2 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)
    locmaj3 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
    locmin3 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)
    locmaj4 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
    locmin4 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)

    ax11.yaxis.set_major_locator(locmaj)
    ax11.yaxis.set_minor_locator(locmin)
    ax21.yaxis.set_major_locator(locmaj2)
    ax21.yaxis.set_minor_locator(locmin2)
    ax31.yaxis.set_major_locator(locmaj3)
    ax31.yaxis.set_minor_locator(locmin3)
    ax41.yaxis.set_major_locator(locmaj4)
    ax41.yaxis.set_minor_locator(locmin4)

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
    
    # RMiss #
#    ax41.set_xlim(RMiss_1tau_nucleon_p_5GeV_min, RMiss_1tau_nucleon_p_5GeV_max)
#    ax42.set_xlim(RMiss_1tau_nucleon_p_10GeV_min, RMiss_1tau_nucleon_p_10GeV_max)
#    ax43.set_xlim(RMiss_1tau_nucleon_p_20GeV_min, RMiss_1tau_nucleon_p_20GeV_max)
#    ax44.set_xlim(RMiss_1tau_nucleon_p_50GeV_min, RMiss_1tau_nucleon_p_50GeV_max)
#    ax45.set_xlim(RMiss_1tau_nucleon_p_200GeV_min, RMiss_1tau_nucleon_p_200GeV_max)
    
    ax41.set_xlim(0, 1)
    ax42.set_xlim(0, 1)
    ax43.set_xlim(0, 1)
    ax44.set_xlim(0, 1)
    ax45.set_xlim(0, 1)

    ax11.set_yscale('log')
    ax21.set_yscale('log')
    ax31.set_yscale('log')
    ax41.set_yscale('log')

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
#            ax[i,j].yaxis.set_major_locator(locmaj)
#            ax[i,j].yaxis.set_minor_locator(locmin)
#            ax[i,j].yaxis.set_major_locator(locmaj)
#            ax[i,j].yaxis.set_minor_locator(locmin)


#fig.suptitle(r'\textbf{Distribution for Incoherent Scattering (proton) off $^{40}$Ar}', fontsize=50)
#fig.supxlabel(r'\textbf{Momentum Transfer} $Q$ (GeV)', fontsize=60)
# Save #
fig1.savefig("../plots/ptMiss_1tau_incoh_p.png", dpi=100, bbox_inches='tight')
fig2.savefig("../plots/ptMuon_1tau_incoh_p.png", dpi=100, bbox_inches='tight')
fig3.savefig("../plots/ptTau_1tau_incoh_p.png", dpi=100, bbox_inches='tight')
fig4.savefig("../plots/RMiss_1tau_incoh_p.png", dpi=100, bbox_inches='tight')
