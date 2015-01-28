def get_sailfish_info(fn,cutoff=1):
    tid2abd = {}
    sf_file = open(fn)
    for line in sf_file:
        if line.startswith("#"): continue
        words = line.rstrip().split()
        header = words[0]
        tpm = float(words[2])
        if tpm >= cutoff: 
            twords = header.split("|")
            tid = twords[0]
            gid = twords[1]
            tid2abd[tid] = tpm
    return tid2abd


tid2abd_truth = get_sailfish_info("quant_truth.sf")
tid2abd_synth = get_sailfish_info("quant_synth.sf")


tpm_truth = []
tpm_synth = []

for tid in tid2abd_truth:
    tpm_truth.append(tid2abd_truth[tid])
    if tid in tid2abd_synth:
        tpm_synth.append(tid2abd_synth[tid])
    else:
        tpm_synth.append(0)

import scipy.stats

print scipy.stats.pearsonr(tpm_truth, tpm_synth)
