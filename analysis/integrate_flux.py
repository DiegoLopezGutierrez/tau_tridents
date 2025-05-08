import numpy as np
from scipy.integrate import simpson
import csv

FLUX_DIR = '../csv/fluxes'

SBND_filename = FLUX_DIR + '/SBND/BNB_SBND_flux.csv'

# SBND #
energy_low_SBND = []
energy_high_SBND = []
bin_widths_SBND = []
flux_SBND_neutrino_vmu = []
flux_SBND_neutrino_vmubar = []
flux_SBND_neutrino_ve = []
flux_SBND_neutrino_vebar = []

with open(SBND_filename,'r') as txtfile:
    data = csv.reader(txtfile, delimiter = ',')
    for row in data:
        elow = float(row[0])
        ehigh = float(row[1])
        bin_width = ehigh - elow
        bin_widths_SBND.append(bin_width)
        energy_low_SBND.append(elow)
        energy_high_SBND.append(ehigh)
        flux_SBND_neutrino_vmu.append(float(row[2]) / 50 / 1e6 * 1e3) # Histogram was in units of [v / m^2 (50 MeV) (1e6 POT)] and is now [v / m^2 GeV POT]
        flux_SBND_neutrino_ve.append(float(row[4]) / 50 / 1e6 * 1e3)
        flux_SBND_neutrino_vmubar.append(float(row[3]) / 50 / 1e6 * 1e3)
        flux_SBND_neutrino_vebar.append(float(row[5]) / 50 / 1e6 * 1e3)

def integrate_hist(widths, heights):
    assert len(widths)==len(heights)
    return sum([widths[i] * heights[i] for i in range(len(widths))])

SBND_neutrino_vmu_integrated_flux = integrate_hist(bin_widths_SBND, flux_SBND_neutrino_vmu)  # Integrated flux [v / m^2 POT]
SBND_neutrino_ve_integrated_flux = integrate_hist(bin_widths_SBND, flux_SBND_neutrino_ve)  # Integrated flux [v / m^2 POT]
SBND_neutrino_vmubar_integrated_flux = integrate_hist(bin_widths_SBND, flux_SBND_neutrino_vmubar)  # Integrated flux [v / m^2 POT]
SBND_neutrino_vebar_integrated_flux = integrate_hist(bin_widths_SBND, flux_SBND_neutrino_vebar)  # Integrated flux [v / m^2 POT]

print("--SBND--")
print("\tvmu: {} [neutrinos / m^2 POT]".format(SBND_neutrino_vmu_integrated_flux))
print("\tve: {} [neutrinos / m^2 POT]".format(SBND_neutrino_ve_integrated_flux))
print("\tvmubar: {} [neutrinos / m^2 POT]".format(SBND_neutrino_vmubar_integrated_flux))
print("\tvebar: {} [neutrinos / m^2 POT]".format(SBND_neutrino_vebar_integrated_flux))

#def probability_distribution(diff_flux, widths, total_flux):
#    assert len(diff_flux) == len(widths)
#    p_array = []
#    for i in range(len(diff_flux)):
#        p = diff_flux[i] * widths[i] / total_flux
#        p_array.append(p)
#
#    return p_array
#
#SBND_probability_distribution_vmu = probability_distribution(flux_SBND_neutrino_vmu, bin_widths_SBND, SBND_neutrino_vmu_integrated_flux)
#SBND_probability_distribution_vmubar = probability_distribution(flux_SBND_neutrino_vmubar, bin_widths_SBND, SBND_neutrino_vmubar_integrated_flux)
#SBND_probability_distribution_ve = probability_distribution(flux_SBND_neutrino_ve, bin_widths_SBND, SBND_neutrino_ve_integrated_flux)
#SBND_probability_distribution_vebar = probability_distribution(flux_SBND_neutrino_vebar, bin_widths_SBND, SBND_neutrino_vebar_integrated_flux)
#
#with open("probability_distribution_SBND_vmu.csv",'w') as fout:
#    writer = csv.writer(fout)
#    assert len(SBND_probability_distribution_vmu) == len(energy_low_SBND)
#    for i in range(len(energy_low_SBND)):
#        writer.writerow([energy_low_SBND[i],energy_high_SBND[i],SBND_probability_distribution_vmu[i]])
