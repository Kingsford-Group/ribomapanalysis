#!/usr/bin/env python
import numpy as np
import scipy.stats
from ribomap_result_parser import *

def pearsonr(vec_cmp, vec_true):
    r = scipy.stats.pearsonr(np.array(vec_true), np.array(vec_cmp))[0]
    return 0 if np.isnan(r) else r

def get_tot_rcnt(tid_order, profile_list):
    return np.array([sum(profile_list[tid]['rprofile']) if tid in profile_list else 0 for tid in tid_order])

if __name__ == "__main__":
    # result_dir="../data/hela/"
    # rmp_fname = result_dir+"GSM546920_filtered_sequence_ribomap.codon"
    # spp_fname = result_dir+"GSM546920_filtered_sequence_starprime.codon"
    # rms_fname = result_dir+"GSM546920_filtered_sequence_ribomap.stats"

    result_dir="../data/mouse/"
    file_core="ribo_mesc_yeslif"
    rmp_fname = result_dir+file_core+"_ribomap.codon"
    spp_fname = result_dir+file_core+"_starprime.codon"
    rms_fname = result_dir+file_core+"_ribomap.stats"

    print "getting ribomap profiles..."
    rmp = parse_estimated_profile(rmp_fname)
    print "getting star prime profiles..."
    spp = parse_estimated_profile(spp_fname)
    print "getting transcript abundance..."
    stats = parse_stats(rms_fname)

    tabd_list = np.array([stats[rid]["tabd"] for rid in stats if stats[rid]["tabd"]!=0 ])
    rid_list = np.array([rid for rid in stats if stats[rid]["tabd"]!=0 ])
    rabd_rm = get_tot_rcnt(rid_list, rmp)
    rabd_sp = get_tot_rcnt(rid_list, spp)
    print "pearson correlation between ribosome loads and transcript abundance"
    print "ribomap: {0}".format(pearsonr(tabd_list, rabd_rm))
    print "star prime: {0}".format(pearsonr(tabd_list, rabd_sp))
