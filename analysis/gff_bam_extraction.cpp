#include <string>
#include <unordered_set>
#include <iostream>
#include <algorithm>

#include <seqan/bam_io.h>

using namespace std;

bool get_multimapper_from_bam(const char *fn, unordered_set<string>& name_list)
{
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
       seqan::BamTagsDict tag(bam_rec.tags);
       unsigned i(0);
       int val(0);
       if (findTagKey(i, tag, seqan::CharString("NH")) and 
	   extractTagValue(val, tag, i)) {
	 if (val==1) continue;
	 else
	   name_list.emplace(string(toCString(bam_rec.qName)));
       }// if NH
     }//if(!hasFlagUnmapped(bam_rec))
  }//while(!atEnd(bamIn))
  return false;
}

bool get_multimapper_from_gffbam(const char *fn, unordered_set<string>& name_list)
{
  int tot_cnt(0), uniq_cnt(0);
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
       seqan::BamTagsDict tag(bam_rec.tags);
       unsigned i(0);
       int val(0);
       if (findTagKey(i, tag, seqan::CharString("HI")) and 
	   extractTagValue(val, tag, i)) {
	 if (val!=1) continue;
	 if (findTagKey(i, tag, seqan::CharString("NH")) and
	     extractTagValue(val, tag, i)) {
	   if (val==1)
	     ++uniq_cnt;
	   else
	     name_list.emplace(string(toCString(bam_rec.qName)));
	   ++tot_cnt;
	 }//if NH
       }// if HI
    }//if(!hasFlagUnmapped(bam_rec))
  }//while(!atEnd(bamIn))
  cout<<"transcripts bam: "<<endl;
  cout<<"total reads: "<<tot_cnt<<" unique mappers: "<<uniq_cnt<<endl;
  return false;
}

int main(int argc, char** argv)
{
  if (argc!=3) {
    cout<<"Usage ./multimapper_info genome_bam transcript_bam"<<endl;
  }
  char* genome_bam = argv[1];
  char* transcript_bam = argv[2];
  unordered_set<string> genome_list, transcript_list, common_list;
  cout<<"getting multimappers from genome bam...";
  get_multimapper_from_bam(genome_bam, genome_list);
  cout<<"multimapped: "<<genome_list.size()<<endl;
  cout<<"getting multimappers from transcript bam...";
  get_multimapper_from_gffbam(transcript_bam, transcript_list);
  cout<<"multimapped: "<<transcript_list.size()<<endl;
  int common(0);
  for (auto r: genome_list)
    if (transcript_list.find(r)!=transcript_list.end())
      ++common;
  cout<<"multimapper caused by repetitive sequence: "<<common<<endl;
  return 1;
}
