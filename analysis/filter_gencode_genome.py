#!/usr/bin/env python
import sys
from Bio import SeqIO

def filter_gencode_genome(ifa, ofa):
    ifile = open(ifa, "rU")
    ofile = open(ofa,'w')
    for rec in SeqIO.parse(ifile, "fasta"):
        if rec.id.startswith("chr"): 
            SeqIO.write(rec,ofile, "fasta")
    ifile.close()
    ofile.close()

def main():
    if len(sys.argv)!=3:
        print "Usage: python filter_gencode_genome.py input_fa output_fa"
        exit(1)
    ifa = sys.argv[1]
    ofa = sys.argv[2]
    filter_gencode_genome(ifa, ofa)

if __name__ == "__main__": main()
