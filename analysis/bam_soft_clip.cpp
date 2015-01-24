#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <numeric>

#include <seqan/bam_io.h>

using namespace std;

bool get_soft_clips_from_bam(const char *fn)
{
  long sbegin(0), send(0), total(0);
  double sb_len(0), se_len(0);
  // open bam file
  seqan::BamStream bamIn(fn);
  if (!isGood(bamIn)){
    cerr << "ERROR: Could not open "<<fn<<"!"<<endl;
    exit(1);
  }
  seqan::BamAlignmentRecord bam_rec;
  while(!atEnd(bamIn)){
    if(readRecord(bam_rec,bamIn)!=0){
      cerr << "ERROR: Could not read bam record!\n";
      return true;
    }
    // get mapped reads
     if (!hasFlagUnmapped(bam_rec)){
       if (hasFlagSecondary(bam_rec)) continue;
       // parse cigar string to figure out read length
       bool skip(false);
       for (int i=0; i<length(bam_rec.cigar); ++i) {
	 switch (bam_rec.cigar[i].operation) {
	 case 'S':
	   if (i==0) {
	     ++sbegin;
	     sb_len += bam_rec.cigar[i].count;
	   }
	   else if (i==length(bam_rec.cigar)-1) {
	     ++send;
	     se_len += bam_rec.cigar[i].count;
	   }
	   break;
	 case 'N':
	   // skipped regions indicate novel splices
	   // since reads are mapped to the transcriptome
	   // skip these alignments for now
	   skip = true;
	   break;
	 case 'M':
	 case 'D':
	 case 'I':
	 default:
	   break;
	 }
	 if (skip) break;
       }
       ++total;
    }//if(!hasFlagUnmapped(bam_rec))
  }//while(!atEnd(bamIn))
  if (sbegin!=0) sb_len /= sbegin;
  if (send!=0) se_len /= send;
  cout<<"soft clipping at the beginning: "<<sbegin<<"("<<sb_len<<")"<<endl;
  cout<<"soft clipping at the end: "<<send<<"("<<se_len<<")"<<endl;
  cout<<"total: "<<total;
  return false;
}

int main(int argc, char** argv)
{
  if (argc!=2) {
    cout<<"Usage ./soft_clip_stats input_bam"<<endl;
  }
  char* bam_fn = argv[1];
  get_soft_clips_from_bam(bam_fn);
  return 1;
}
