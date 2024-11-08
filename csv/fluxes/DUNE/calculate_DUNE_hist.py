import csv

DUNE_filename = 'DUNE_diff_flux.csv'  # dPhi/dE

energy_DUNE = []
flux_DUNE = []

with open(DUNE_filename,'r') as csvfile:
    data = csv.reader(csvfile, delimiter = ',')
    for row in data:
        energy_DUNE.append(float(row[0]))
        flux_DUNE.append(float(row[1]) * 1e4) # DUNE differential flux is in [cm^-2 GeV^-1 POT^-1]; convert to [m^-2 GeV^-1 POT^-1]

bin_edges = []
with open('bin_edges_120.txt') as txt_file:
    for row in txt_file:
        bin_edges.append(float(row))

with open('DUNE_hist.csv','w',newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    for i, b in enumerate(bin_edges):
        if i == (len(bin_edges)-1):
            continue
        low_bin = bin_edges[i]
        high_bin = bin_edges[i+1]
        flux_bin = 0
        for energy, flux in zip(energy_DUNE, flux_DUNE):
            if energy >= low_bin and energy <= high_bin:
                flux_bin = flux
            else:
                continue
        writer.writerow([low_bin, "\t"+str(high_bin), "\t"+str(flux_bin)])
