#!/bin/sh

a="8044 31975 824 567 714"
ref=31671

module load mummer python3

for b in $a; do
    nucmer -prefix=$ref\_$b ../$ref*.fasta ../$b*.fasta
    delta-filter -1 -r $ref\_$b\.delta | show-coords -r -T -l -d -c /dev/stdin | awk 'NR>4' | sed 's/quiver//g' | sed 's/pilon//g' | sed 's/|//g' | sort -k14,14 -k15,15 -k1n,1 > raw.$ref\_$b\.txt
    python3 crude_search_inversion.py raw.$ref\_$b\.txt $ref\_$b\.inversions 
done
