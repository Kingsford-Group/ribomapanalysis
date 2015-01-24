#!/bin/bash
cd ~/scratch/scratch2/ribomapanalysis/
./hela_ribo_analysis.sh
cd ~/scratch/scratch2/software_testing/ribomap/scripts/footprint_generator/
./footprint_generation_pipe.sh
cd ~/scratch/scratch2/ribomapanalysis/
./synth_ribo_analysis.sh
./mouse_ribo_analysis.sh
