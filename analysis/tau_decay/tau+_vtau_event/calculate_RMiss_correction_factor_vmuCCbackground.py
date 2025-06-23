import numpy as np
import csv

DIST_DIR_1TAU = '../../../csv/distributions/vmu_to_vtau_tau+_mu-'

RMiss_background_filename = '../../../csv/events/AlexSousa_RMiss_Background.csv'
background_value_RMiss = []
background_weight_RMiss = []

RMiss_NuWro_background_filename = DIST_DIR_1TAU+'/RT-thresholds-no-single.txt'
background_value_RMiss_NuWro = []

R_cut = float(input("Cut value for RMiss: "))

RMiss_cut_area = 0
RMiss_total_area = 0

RMiss_NuWro_pass = 0
RMiss_NuWro_total = 0

## TO-DO: Consider calculating the correction factor for background events as well.
with open(RMiss_background_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for row in data:
        rmiss = float(row[0])
        wgt = float(row[1])
#        if rmiss >= R_cut:
#            RMiss_cut_area += rmiss * wgt
#        RMiss_total_area += rmiss * wgt

        background_value_RMiss.append(rmiss)     # This is the actual RMiss value (i.e. x-value)
        background_weight_RMiss.append(wgt)    # Since the histogram from Alex Sousa is digitized, these are the weights (i.e. y-value) of the histogram.
#    correction_factor = RMiss_cut_area / RMiss_total_area
#    print(f"For Alex Sousa's vmuCC background, for an R_cut value of {R_cut}, only {correction_factor} of the total background events pass the cut.")

for i in range(39):
    rmiss_low = background_value_RMiss[i]
    rmiss_high = background_value_RMiss[i+1]
    delta_rmiss = rmiss_high - rmiss_low
    wgt = background_weight_RMiss[i]
    if rmiss_low >= R_cut:
        RMiss_cut_area += delta_rmiss * wgt
    RMiss_total_area += delta_rmiss * wgt
correction_factor = RMiss_cut_area / RMiss_total_area
print(f"For Alex Sousa's vmuCC background, for an R_cut value of {R_cut}, only {correction_factor} of the total background events pass the cut.")

with open(RMiss_NuWro_background_filename,'r') as txtfile:
    for line in txtfile:
        rmiss = float(line)
        if rmiss >= R_cut:
            RMiss_NuWro_pass += 1
        RMiss_NuWro_total += 1

        background_value_RMiss_NuWro.append(rmiss)  # Pedro's NuWro simulation is just the actual RMiss values (i.e. x-values) that need to be binned and histogramed8
    correction_factor_NuWro = RMiss_NuWro_pass / RMiss_NuWro_total
    print(f"For NuWro's vmuCC background, for an R_cut value of {R_cut}, only {correction_factor_NuWro} of the total background events pass the cut.")
