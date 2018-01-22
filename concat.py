import fileinput
outfilename = "WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6.txt"
filenames = ['WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6_ext1-v1.txt','WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6_ext2-v1.txt','WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6_ext3-v1.txt']
with open(outfilename, 'w') as fout:
    fin = fileinput.input(filenames)
    for line in fin:
        fout.write(line)
    fin.close()

file1 = open("WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6_ext1-v1_weight.log").readlines()
for lines in file1:
    k=lines.split()
    print k[1]

file2 = open("WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6_ext2-v1_weight.log").readlines()
for lines in file2:
    l=lines.split()
    print l[1]

file3 = open("WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6_ext3-v1_weight.log").readlines()
for lines in file3:
    n=lines.split()
    print n[1]

m = str(format(float(k[1])+float(l[1])+float(n[1]), '.6f'))#+float(n[1])

outfile = "WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6_weight.log"
with open(outfile, 'w') as fout:
    fout.write(k[0] + ' ' + m)
