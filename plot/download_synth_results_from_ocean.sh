#!/bin/bash
idir="ocean:~/scratch/scratch2/ribomap-playground/hela/synthetic_data/"
odir=../data/synth_results/
for e in 005 01 02; do
    scp ${idir}ribomap/synth_riboseq_$e.codon ${odir}synth_riboseq_${e}_ribomap.codon
    scp ${idir}star_prime/synth_riboseq_$e.codon ${odir}synth_riboseq_${e}_starprime.codon
done

scp ${idir}synth_riboseq.profile ${odir}
scp ${idir}quant_bias_corrected.sf ${odir}quant_truth.sf
scp ${idir}sm_align_quant/quant_bias_corrected.sf ${odir}quant_synth.sf


