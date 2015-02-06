#!/usr/bin/env python
def parse_true_profile(fname):
    print "start parsing ground truth profile..."
    ptruth = {}
    transcript = {}
    tf = open(fname)
    line = tf.readline()
    while line:
        if line.startswith("refID: "):
            if transcript:
                rid = transcript["id"]
                ptruth[rid] = transcript.copy()
                transcript.clear()
            rid = int(line.lstrip("refID: ").rstrip("\n"))
            transcript["id"] = rid
        elif line.startswith("tid: "):
            tid = line.lstrip("tid: ").rstrip("\n")
            transcript["tid"] = tid
        elif line.startswith("irate: "):
            irate = float(line.lstrip("irate: ").rstrip("\n"))
            transcript["irate"] = irate
        elif line.startswith("trate: "):
            trate = float(line.lstrip("trate: ").rstrip("\n"))
            transcript["trate"] = trate
        elif line.startswith("plen: "):
            length = int(line.lstrip("plen: ").rstrip("\n"))
            transcript["len"] = length
        elif line.startswith("tcnt: "):
            # tcnt: total number of transcripts
            tcnt = int(line.lstrip("tcnt: ").rstrip("\n"))
            transcript["tcnt"] = tcnt
        elif line.startswith("tabd: "):
            # tabd: tabd per base X tlen
            tabd = float(line.lstrip("tabd: ").rstrip("\n"))
            transcript["tabd"] = tabd
        elif line.startswith("rprofile: "):
            rprofile = map(int, line.lstrip("rprofile: ").rstrip("\n").split())
            transcript["rprofile"] = rprofile
        elif line.startswith("mprofile: "):
            mprofile = map(float, line.lstrip("mprofile: ").rstrip("\n").split())
            transcript["mprofile"] = mprofile
        else:
            print "cannot add property {0}!".format(line.split()[0][:-1])
            pass
        line = tf.readline()
    if transcript:
        rid = transcript["id"]
        ptruth[rid] = transcript
    tf.close()
    return ptruth

def parse_estimated_profile(fname):
    pest = {}
    transcript = {}
    tf = open(fname)
    line = tf.readline()
    while line:
        if line.startswith("refID: "):
            if transcript:
                rid = transcript["id"]
                pest[rid] = transcript.copy()
                transcript.clear()
            rid = int(line.lstrip("refID: ").rstrip("\n"))
            transcript["id"] = rid
        elif line.startswith("tid: "):
            tid = line.lstrip("tid: ").rstrip("\n")
            transcript["tid"] = tid
        elif line.startswith("ribo profile: "):
            rprofile = map(float, line.lstrip("ribo profile: ").rstrip("\n").split())
            transcript["rprofile"] = rprofile
        else:
            pass
        line = tf.readline()
    if transcript:
        rid = transcript["id"]
        pest[rid] = transcript
    tf.close()
    return pest

def parse_stats(fname):
    stats = {}
    transcript = {}
    tf = open(fname)
    line = tf.readline()
    while line:
        if line.startswith("refID: "):
            if transcript:
                rid = transcript["id"]
                stats[rid] = transcript.copy()
                transcript.clear()
            rid = int(line.lstrip("refID: ").rstrip("\n"))
            transcript["id"] = rid
        elif line.startswith("tid: "):
            tid = line.lstrip("tid: ").rstrip("\n")
            transcript["tid"] = tid
        elif line.startswith("rabd: "):
            rabd = float(line.lstrip("rabd: ").rstrip("\n"))
            transcript["rabd"] = rabd
        elif line.startswith("tabd: "):
            tabd = float(line.lstrip("tabd: ").rstrip("\n"))
            transcript["tabd"] = tabd
        elif line.startswith("te: "):
            te = float(line.lstrip("te: ").rstrip("\n"))
            transcript["te"] = te
        else:
            pass
        line = tf.readline()
    if transcript:
        rid = transcript["id"]
        stats[rid] = transcript
    tf.close()
    return stats
