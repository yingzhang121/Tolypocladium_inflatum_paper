#!/usr/bin/env python

import sys, glob

#Assembly                    Ti_714_70.pilon.filter  Ti_31671_revcomp.pilon.filter  Ti_31975_revcomp.pilon.filter  824_70_revcomp.pilon.filter  Ti_CBS567_84_revcomp.pilon.filter  NRRL_8044_revcomp.pilon.filter
#Total length                29861265                30053927                       30421992                       31582666                     31764937                           31230432

gsize = {"714":29861265, "824":31582666, "8044":31230432, "567":31764937, "31671":30053927, "31975":30421992}
genomes = ["714", "31671", "31975", "824", "567", "8044"]
results = {}

for ref in genomes:
    files = glob.glob(ref+"*/*.inversions")
    result = {}
    result[ref] = 0
    for f in files:
        qry = f.split("/")[1].split(".")[0].split("_")[1]
        ref_inv = 0
        qry_inv = 0
        for line in open(f):
            fields = line.split()
            ref_inv += int(fields[2])-int(fields[1])
            qry_inv += int(fields[5])-int(fields[4])
        result[qry] = (ref_inv/gsize[ref] + qry_inv/gsize[qry])/2
    results[ref] = result

print("\t".join(genomes))
for ref in genomes:
    printline = []
    for qry in genomes:  
        printline.append("%.4f" % results[ref][qry])
    print("\t".join(printline))
