
#!/bin/bash
# idir="ocean:~/scratch/scratch2/ribomap-playground/hela/"
# file_core=GSM546920_filtered_sequence
# odir=../data/hela/
# mkdir -p $odir
# scp ${idir}ribomap/${file_core}.base ${odir}${file_core}_ribomap.base
# scp ${idir}ribomap/${file_core}.codon ${odir}${file_core}_ribomap.codon
# scp ${idir}ribomap/${file_core}.stats ${odir}${file_core}_ribomap.stats
# scp ${idir}star_prime/${file_core}.base ${odir}${file_core}_starprime.base
# scp ${idir}star_prime/${file_core}.codon ${odir}${file_core}_starprime.codon
#scp ${idir}alignment/map_cnt.txt ${odir}

idir="ocean:~/scratch/scratch2/ribomap-playground/mouse/"
file_core=ribo_mesc_yeslif
odir=../data/mouse/
mkdir -p $odir
scp ${idir}ribomap/${file_core}.codon ${odir}${file_core}_ribomap.codon
scp ${idir}star_prime/${file_core}.codon ${odir}${file_core}_starprime.codon
scp ${idir}ribomap/${file_core}.stats ${odir}${file_core}_ribomap.stats
# scp ${idir}alignment/map_cnt.txt ${odir}






