CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wno-deprecated-declarations -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 
INC = -I/home/hw1/.local/include/ -I/opt/local/include
CXXFLAGS += $(INC)
COMPILE.c = $(CXX) $(CXXFLAGS)
LDFLAGS = -lz -lbz2
OUTPUT_OPTION = -o $@

bam_info: bam_info_extraction.o
	$(COMPILE.c) $(OUTPUT_OPTION) $^ $(LDFLAGS)
	chmod u+x $@

soft_clip_stats: bam_soft_clip.o
	$(COMPILE.c) $(OUTPUT_OPTION) $^ $(LDFLAGS)
	chmod u+x $@

multimapper_info: genome_common_transcript_bam.o
	$(COMPILE.c) $(OUTPUT_OPTION) $^ $(LDFLAGS)
	chmod u+x $@

