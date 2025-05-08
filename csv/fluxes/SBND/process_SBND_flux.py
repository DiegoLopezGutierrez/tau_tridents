import csv

with open('SBND_Fluxes_NeutrinoMode.dat.txt','r') as fin, open('BNB_SBND_flux.csv','w') as fout:
    print("Reading SBND_Fluxes_NeutrinoMode.dat.txt")
    writer = csv.writer(fout)
    i = 1
    for line in fin:
        print(f"Line {i}")
        if i <= 6:
            i += 1
            continue
        writer.writerow(line.split())
        print(f"Writing line {i}")
        i += 1
