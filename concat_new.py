import fileinput
outfilename = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6.txt"
filenames = ['DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext1-v2.txt','DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext2-v1.txt']#,'EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_v6_ext2-v1.txt']
with open(outfilename, 'w') as fout:
    fin = fileinput.input(filenames)
    for line in fin:
        fout.write(line)
    fin.close()

file1 = open("DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext2-v1_weight.log").readlines()
for lines in file1:
    k=lines.split()
    print k[1]

file2 = open("DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext1-v2_weight.log").readlines()
for lines in file2:
    l=lines.split()
    print l[1]

#file3 = open("EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_v6_ext2-v1_weight.log").readlines()
#for lines in file3:
#    n=lines.split()
#    print n[1]

m = str(format(float(k[1])+float(l[1]), '.6f'))#+float(n[1])

outfile = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_weight.log"
with open(outfile, 'w') as fout:
    fout.write(k[0] + ' ' + m)
