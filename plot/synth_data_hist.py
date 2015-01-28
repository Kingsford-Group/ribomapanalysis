#!/usr/bin/env python
from ribomap_result_parser import *
import numpy as np
import scipy.stats
import os
from operator import itemgetter
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 20
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5
colormap = plt.cm.cool

#=======================================
# routine functions
#=======================================
def normalize_vector(vec):
    s = float(sum(vec))
    return [ e/s for e in vec] if s!=0 else vec

def norm(vec):
    return np.sqrt(np.mean(np.square(np.array(vec))))

def RMSE(vec_cmp, vec_true):
    return np.sqrt(np.mean(np.square(np.array(vec_true) - np.array(vec_cmp))))

def pearsonr(vec_cmp, vec_true):
    r = scipy.stats.pearsonr(np.array(vec_true), np.array(vec_cmp))[0]
    return 0 if np.isnan(r) else r

def valid_profile(profile):
    return np.all((profile>0)&(profile<1))

#=======================================
# analysis and plot
#=======================================
def get_tot_rcnt(tid_order, profile_list):
    return np.array([sum(profile_list[tid]['rprofile']) if tid in profile_list else 0 for tid in tid_order])

def ptruth_stats(tid_order, ptruth):
    rcnt_true = np.array([sum(ptruth[tid]['rprofile'])  if tid in ptruth else 0 for tid in tid_order])
    irate_true = np.array([ptruth[tid]['irate'] if tid in ptruth else 0 for tid in tid_order])
    trate_true = np.array([ptruth[tid]['trate'] if tid in ptruth else 0 for tid in tid_order])
    lens = np.array([ptruth[tid]['len'] if tid in ptruth else 0 for tid in tid_order])
    tcnt = np.array([ptruth[tid]['tcnt'] if tid in ptruth else 0 for tid in tid_order])
    return rcnt_true, irate_true, trate_true, lens, tcnt 

def compare_assignment_profile(tid_order, ptruth):
    print "comparing models with counts..."
    true_cm_list = []
    for tid in tid_order:
        if tid not in ptruth:
            true_cm = 0
        else:
            true_cnt = np.array(ptruth[tid]["rprofile"])
            true_mod = np.array(ptruth[tid]["mprofile"])
            true_cm = pearsonr(true_cnt, true_mod)
        true_cm_list.append(true_cm)
    true_cm_list = np.array(true_cm_list)
    return true_cm_list

def compare_bb_with_pmodel(tid_order, ptruth, pmodel, pbb):
    print "comparing ribofish with bowtie best..."
    corr_mod = []
    corr_bb = []
    rmse_mod = []
    rmse_bb = []
    i = 0
    for tid in tid_order:
        # synthetic data index inconsistent with assignment data index
        # quick fix off-by-one
        if tid not in ptruth:
            cr_mod = 0
            cr_bb = 0
            if tid not in pmodel:
                emod = 0
            else:
                emod = norm(pmodel[tid]["rprofile"][1:-1])
            if tid not in pbb:
                ebb = 0
            else:
                ebb = norm(pbb[tid]["rprofile"][1:-1])
        else:
            true_cnt = np.array(ptruth[tid]["rprofile"])
            true_mod = np.array(ptruth[tid]["mprofile"])
            if tid not in pmodel:
                cr_mod = 0
                emod = norm(true_cnt)
            else:
                tasep_cnt = np.array(pmodel[tid]["rprofile"][1:-1])
                cr_mod = pearsonr(tasep_cnt, true_cnt)
                emod = RMSE(tasep_cnt, true_cnt)
            if tid not in pbb:
                cr_bb = 0
                ebb = norm(true_cnt)
            else:
                bb_cnt = np.array(pbb[tid]['rprofile'][1:-1])
                cr_bb = pearsonr(bb_cnt, true_cnt)
                ebb = RMSE(bb_cnt, true_cnt)
        corr_mod.append(cr_mod)
        corr_bb.append(cr_bb)
        rmse_mod.append(emod)
        rmse_bb.append(ebb)
    corr_mod = np.array(corr_mod)
    corr_bb = np.array(corr_bb)
    rmse_mod = np.array(rmse_mod)
    rmse_bb = np.array(rmse_bb)
    return corr_mod, corr_bb, rmse_mod, rmse_bb


def plot_hist_compare_bb(corr_mod, corr_bb, rmse_mod, rmse_bb, true_cm, prefix):
    print "pearsonr mean: ribomap: {0} star prime: {1}".format(np.mean(corr_mod), np.mean(corr_bb))
    print "pearsonr median: ribomap: {0} star prime: {1}".format(np.median(corr_mod), np.median(corr_bb))
    print "tasep beat bb on pearsonr (all):", np.mean(corr_mod > corr_bb)
    print "mwu:", scipy.stats.mannwhitneyu(corr_mod, corr_bb)
    print "tasep beat bb on rmse (all):", np.mean(rmse_mod < rmse_bb)
    print "mwu:", scipy.stats.mannwhitneyu(rmse_mod, rmse_bb)
    #cond = (objs>0) & (true_cm>0.8)
    cond = true_cm> -1
    print "tasep beat bb on pearsonr (valid):", np.mean(corr_mod[cond] > corr_bb[cond])
    print "tasep beat bb on rmse (valid):", np.mean(rmse_mod[cond] < rmse_bb[cond])
    plt.figure()
    ns, bins, patches = plt.hist(corr_mod[cond], bins=np.arange(-0.1,1.1,0.05), histtype='stepfilled', color = 'c', edgecolor='c', alpha = 0.5, label="Ribomap")
    ns, bins, patches = plt.hist(corr_bb[cond], bins=np.arange(-0.1,1.1,0.05), histtype='stepfilled', color = 'b', hatch='/', edgecolor='b', alpha = 0.5, label = "Star prime")
    plt.xlabel("pearson correlation")
    plt.ylabel("number of transcripts")
    plt.legend(loc=9, frameon=False, prop={'size':20})
    plt.xlim(0,1.05)
    plt.savefig(prefix+"pearsonr_hist.png", bbox_inches="tight")
    plt.figure()
    ns, bins, patches = plt.hist(rmse_mod[cond], np.arange(0,2,0.05), histtype='stepfilled', color = 'c', edgecolor='c', alpha = 0.5, label="Ribomap")
    ns, bins, patches = plt.hist(rmse_bb[cond], np.arange(0,2,0.05), histtype='stepfilled', color = 'b', hatch='/', edgecolor='b', alpha = 0.5, label = "Star prime")
    plt.xlabel("RMSE")
    plt.ylabel("number of transcripts")
    plt.legend(loc=1, frameon=False, prop={'size':20})
    plt.savefig(prefix+"rmse_hist.png", bbox_inches="tight")

def plot_synthetic_profile(ptruth):
    i = 0
    for k,v in ptruth.iteritems():
        ptarget = v['mprofile']
        ptarget = normalize_vector(ptarget)
        psimulate = v['rprofile']
        psimulate = normalize_vector(psimulate)
        width = int(len(psimulate)/10.0*1.5)
        plt.figure(figsize=(width,6))
        plt.xlabel('codon position')
        plt.ylabel('probability')
        plt.plot(ptarget,'--o', alpha=0.5)
        plt.plot(psimulate,'--o', alpha=0.5)
        plt.xlim(-1,len(psimulate)+1)
        plt.title("len: {0} abundance {1:.3e}".format(len(psimulate), v['tabd']))
        plt.legend(['target', 'simulation'])
        plt.show()
        i += 1
        if i>10: break

if __name__ == "__main__":
    result_dir="../data/synth_results/"
    truth_fname = result_dir+"synth_riboseq.profile"
    ribomap_fname = lambda e: result_dir+"synth_riboseq_{0}_ribomap.codon".format(e)
    sp_fname = lambda e: result_dir+"synth_riboseq_{0}_starprime.codon".format(e)
    for e in [ '005', '01', '02' ]:
        print "error rate: 0.{0}".format(e)
        ptruth = parse_true_profile(truth_fname)
        #plot_synthetic_profile(ptruth)
        print "getting ribomap results..."
        pribomap = parse_estimated_profile(ribomap_fname(e))
        print "getting star prime results..."
        psp = parse_estimated_profile(sp_fname(e))
        print "making result list..."
        tid_list = np.array([ tid for tid in ptruth if sum(ptruth[tid]['rprofile'])!=0 ])
        rcnt_true, irate_true, trate_true, lens, tcnt = ptruth_stats(tid_list, ptruth)
        order = np.argsort(rcnt_true)[::-1]
        tid_order = tid_list[order]
        rcnt_guess = get_tot_rcnt(tid_order, pribomap)
        print "pearson correlation between true and esimated ribosome loads:", scipy.stats.pearsonr(rcnt_true[order], rcnt_guess)


        # # eye balling whether the estimated profiles are fishy
        # i = 0
        # for tid in tid_order:
        #     if tid not in ptruth or tid not in pribomap or tid not in psp: continue
        #     print ptruth[tid]['tid']
        #     true_cnt = np.array(ptruth[tid]["rprofile"])
        #     true_mod = np.array(ptruth[tid]["mprofile"])
        #     tasep_cnt = np.array(pribomap[tid]["rprofile"][1:-1])
        #     bb_cnt = np.array(psp[tid]['rprofile'][1:-1])
        #     width = int(len(true_cnt)/10.0*1.5)
        #     plt.figure(figsize=(width,6))
        #     plt.xlabel('codon position')
        #     plt.ylabel('footprint count')
        #     plt.plot(true_cnt,'g--o', alpha=0.5)
        #     plt.plot(tasep_cnt,'b--o', alpha=0.5)
        #     plt.plot(bb_cnt, 'r--o', alpha=0.5)
        #     plt.xlim(-1,len(true_cnt)+1)
        #     plt.legend(['truth', 'ribomap', 'starprime'], fontsize=10, frameon=False)
        #     plt.show()
        #     i += 1
        #     if i>10: break

        true_cm_list = compare_assignment_profile(tid_order, ptruth)
        corr_mod, corr_bb, rmse_mod, rmse_bb = compare_bb_with_pmodel(tid_order, ptruth, pribomap, psp)
        prefix = "synth_{0}_".format(e)
        plot_hist_compare_bb(corr_mod, corr_bb, rmse_mod, rmse_bb, true_cm_list, prefix)
