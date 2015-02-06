#!/usr/bin/env python
import numpy as np
from ribomap_result_parser import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 20
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5

turquoise = (64,224,208)
LightSlateBlue = (132,112,255)


def get_sailfish_info(fn,cutoff=0):
    tid2abd = {}
    gid2tids = {}
    tid2gid = {}
    sf_file = open(fn)
    for line in sf_file:
        if line.startswith("#"): continue
        words = line.rstrip().split()
        header = words[0]
        tpm = float(words[2])
        if tpm > cutoff: 
            twords = header.split("|")
            tid = twords[0]
            gid = twords[1]
            tid2abd[tid] = tpm
            gid2tids.setdefault(gid, []).append(tid)
            tid2gid[tid] = gid
    return tid2abd, gid2tids, tid2gid

def build_tid_rid_map(profile_list):
    tid2rid = {}
    rid2tid = {}
    for rid,t in profile_list.iteritems():
        tid = t["tid"].split("|")[0]
        tid2rid[tid] = rid
        rid2tid[rid] = tid
    return tid2rid, rid2tid

def get_nz_cvg(rid_order, profile_list):
    return np.array([ np.mean(np.array(profile_list[rid]['rprofile'])!=0) for rid in rid_order])

def get_genome_corrd(gtfname):
    """
    gtf fields:
    0: seqname
    1: source
    2: feature
    3: start
    4: end
    5: score
    6: strand
    7: frame
    8: attribute
      gene_id ENSGXXXXXXXXXXX *
      transcript_id ENSTXXXXXXXXXXX *
      gene_type list of biotypes
      gene_status {KNOWN, NOVEL, PUTATIVE}
      gene_name string
      transcript_type list of biotypes
      transcript_status {KNOWN, NOVEL, PUTATIVE}
      transcript_name string
      level 1 (verified loci), 2 (manually annotated loci), 3 (automatically annotated loci)
    """
    print "getting exon list in CDS regions"
    gtfile = open(gtfname)
    tid2exon = {}
    tid2cds = {}
    gid2range = {}
    for line in gtfile:
        if line.startswith("#"): continue
        words = line.rstrip('\n').split('\t')
        if words[2] == "CDS":
            start = int(words[3])
            end = int(words[4])
            strand = words[6]
            frame = int(words[7])
            assert words[8].find("exon_number") != -1
            attributes = words[8].split("; ")
            for a in attributes:
                aname, aval = a.lstrip().split(" ")
                aval = aval.lstrip('"').rstrip('"')
                if aname == "transcript_id":
                    tid = aval
                if aname == "exon_number":
                    enum = int(aval)
            tid2cds.setdefault(tid,{0:strand})
            tid2cds[tid][enum] = (start,end,frame)
        elif words[2] == "gene":
            chrm = words[0]
            start = int(words[3])
            end = int(words[4])
            attributes = words[8].split("; ")
            for a in attributes:
                aname, aval = a.lstrip().split(" ")
                aval = aval.lstrip('"').rstrip('"')
                if aname == "gene_id":
                    gid = aval
            gid2range[gid] = (chrm, start, end)
        elif words[2] == "exon":
            start = int(words[3])
            end = int(words[4])
            strand = words[6]
            assert words[8].find("exon_number") != -1
            attributes = words[8].split("; ")
            for a in attributes:
                aname, aval = a.lstrip().split(" ")
                aval = aval.lstrip('"').rstrip('"')
                if aname == "transcript_id":
                    tid = aval
                if aname == "exon_number":
                    enum = int(aval)
            tid2exon.setdefault(tid,{0:strand})
            tid2exon[tid][enum] = (start,end)
            
    gtfile.close()
    return tid2exon, tid2cds, gid2range

def tcoord_to_gcoord(plen, elist):
    exon_list = sorted(elist.keys())
    strand = elist[0]
    coord = []
    ei = 1
    e = exon_list[ei]
    start, end, frame = elist[e]
    if strand == '+':
        cur = start+frame
        coord.append(cur)
        for i in xrange(plen-1):
            cur += 3
            if cur>end:
                f = cur-end
                ei += 1
                if ei == len(exon_list): 
                    coord.append(cur)
                    break
                e = exon_list[ei]
                start, end, frame = elist[e]
                if (f != frame+1):
                    print "frame{0} error! {1}".format(frame, f)
                cur = start+frame
            coord.append(cur)
    elif strand == '-':
        cur = end-frame
        coord.append(cur)
        for i in xrange(plen-1):
            cur -= 3
            if cur<start:
                f = start - cur
                ei += 1
                if ei == len(exon_list):
                    coord.append(cur)
                    break
                e = exon_list[ei]
                start, end, frame = elist[e]
                if (f != frame+1):
                    print "frame{0} error! {1}".format(frame, f)
                cur = end - frame
            coord.append(cur)
    else:
        print "strand type not recognized! {0}".format(strand)
        return []
    return coord

def tcoord_to_gcoord(elist):
    exon_list = sorted(elist.keys())
    strand = elist[0]
    coord = []
    for e in exon_list[1:]:
        start, end = elist[e]
        if strand == '+':
            coord.extend(range(start, end+1))
        elif strand == '-':
            coord.extend(range(end,start-1,-1))
        else:
            print "strand type not recognized! {0}".format(strand)
            return []
    return coord

def check_gene_isoform(gid, gid2topt, tid2exon):
    exon_list = [ set(tid2exon[tid].values())-set(tid2exon[tid][0]) for tid in gid2topt[gid] ]
    exon_include = [ exon_list[0] ]
    for i in xrange(1,len(exon_list)):
        if exon_list[i] in exon_include:
            return False
        exon_include.append(exon_list[i])
    return True

def write_transcript_wiggle(tf, tid, ptype, chrm, tid2rid, tid2exon, profile, color, priority):
    rid = tid2rid[tid]
    if rid not in profile:
        return
    p = profile[rid]['rprofile']
    if sum(np.array(p)>1) == 0: 
        return
    track_line = 'track type=wiggle_0 name="{0}_{1}" visibility=2 color={2},{3},{4} db=hg19 priority={5} graphType=bar viewLimits=1:{6}\n'.format(tid, ptype, color[0], color[1], color[2], priority, max(p))
    tf.write(track_line)
    elist = tid2exon[tid]
    gcoord=tcoord_to_gcoord(elist)
    if gcoord[0] > gcoord[-1]:
        gcoord = gcoord[::-1]
        p = p[::-1]
    step_line = "variableStep chrom={0}\n".format(chrm)
    tf.write(step_line)
    value_lines = [ "{0}\t{1}\n".format(gcoord[i], int(p[i])) for i in xrange(len(p)) if p[i]>1 ]
    tf.writelines(value_lines)

def write_gene_wiggle(fname, chrm, start, end, tid_list, tid2rid, tid2exon, rmp, spp):
    tf = open(fname, 'w')
    i = 1    
    browser_line = "browser position: {0}:{1}-{2}\n".format(chrm,start,end)
    tf.write(browser_line)
    for tid in tid_list:
        write_transcript_wiggle(tf, tid, "rm", chrm, tid2rid, tid2exon, rmp, turquoise, i)
        i += 1
        write_transcript_wiggle(tf, tid, "sp", chrm, tid2rid, tid2exon, spp, LightSlateBlue, i)
        i += 1
    tf.close()

if __name__ == "__main__":
    result_dir = "../data/hela/"
    file_core = "GSM546920_filtered_sequence"
    rmp_cds_fname = result_dir + file_core + "_ribomap.codon"
    spp_cds_fname = result_dir + file_core + "_starprime.codon"
    rmp_fname = result_dir + file_core + "_ribomap.base"
    spp_fname = result_dir + file_core + "_starprime.base"
    sf_fname = result_dir + "quant_bias_corrected.sf"
    gtf_fname = "/Users/hao/Study/ckgroup/riboseq/GenCode/gencode.v18.annotation.gtf"

    tid2exon, tid2cds, gid2range = get_genome_corrd(gtf_fname)
    print "getting sailfish result..."
    tid2abd, gid2tids, tid2gid = get_sailfish_info(sf_fname)
    tid_list = tid2abd.keys()
    print "getting ribomap cds profiles..."
    rmp_cds = parse_estimated_profile(rmp_cds_fname)
    # print "getting star prime cds profiles..."
    # spp_cds = parse_estimated_profile(spp_cds_fname)
    print "getting ribomap profiles..."
    rmp = parse_estimated_profile(rmp_fname)
    print "getting star prime profiles..."
    spp = parse_estimated_profile(spp_fname)
    print "building tid2rid"
    tid2rid, rid2tid=build_tid_rid_map(rmp)
    print "getting rid list"
    rid_list = np.array([ tid2rid[tid] for tid in tid_list if tid in tid2rid ])
    print "getting coverage..."
    nz_cvg=get_nz_cvg(rid_list, rmp_cds)
    print "collecting gene count"
    gid2topt = {}
    rid_cand = rid_list[nz_cvg>0.5]
    tid_cand = [ rid2tid[rid] for rid in rid_cand ]
    for tid in tid_cand:
        gid = tid2gid[tid]
        gid2topt.setdefault(gid, []).append(tid)
    gid_list = [ gid for gid,tlist in gid2topt.iteritems() if len(tlist)>1 ]
    print len(gid_list)
    gid_list = [ gid for gid in gid_list if check_gene_isoform(gid, gid2topt, tid2cds) ]
    print len(gid_list)

    for gid in gid_list:
        chrm,start,end = gid2range[gid]
        tid_list = gid2tids[gid]
        fname = "wig_folder/{0}.wig".format(gid)
        write_gene_wiggle(fname, chrm, start, end, tid_list, tid2rid, tid2exon, rmp, spp)

