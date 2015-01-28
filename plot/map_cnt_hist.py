#!/usr/bin/env python
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 20
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5

def read_numbers_from_file(fn):
    tf = open(fn)
    counts = map(int, tf.readline().strip().split())
    tf.close()
    return counts

def main():
    if len(sys.argv)!=3:
        print "Usage: python map_cnt_hist.py map_cnt.txt output_hist_name"
        exit(1)
    fn = sys.argv[1]
    ofn= sys.argv[2]
    counts = read_numbers_from_file(fn)
    plt.figure()
    ns,bins,patches=plt.hist(counts, range(1,21), normed=True, alpha=0.5)
    plt.bar(1,ns[0], 1, color='r')
    plt.xlabel("number of mapped location")
    plt.ylabel("percentage of mapped reads")
    plt.savefig(ofn,bbox_inches="tight")

if __name__=="__main__": main()
