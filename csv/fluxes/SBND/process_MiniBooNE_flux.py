import csv

with open('MiniBooNE_Fluxes_NuMode.dat.txt','r') as fin, open('BNB_MiniBooNE_flux.csv','w') as fout:
    print("Reading MiniBooNE_Fluxes_NuMode.dat.txt")
    writer = csv.writer(fout)
    i = 1
    for line in fin:
        print(f"Line {i}")
        if i <= 9:
            i += 1
            continue
        writer.writerow(line.split())
        print(f"Writing line {i}")
        i += 1
