#include <string>
#include <unordered_set>
#include <iostream>
#include <algorithm>

#include <seqan/bam_io.h>

using namespace std;

int get_multimapper_from_bam(const char *fn, unordered_set<string>& name_list)
{
  int tot_cnt(0);
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
	 ++tot_cnt;
	 if (val==1) continue;
	 else
	   name_list.emplace(string(toCString(bam_rec.qName)));
       }// if NH
     }//if(!hasFlagUnmapped(bam_rec))
  }//while(!atEnd(bamIn))
  return tot_cnt;
}

int get_uniqmapper_from_bam(const char *fn, unordered_set<string>& name_list)
{
  int tot_cnt(0);
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
	 ++tot_cnt;
	 if (val==1) 
	   name_list.emplace(string(toCString(bam_rec.qName)));
       }// if NH
     }//if(!hasFlagUnmapped(bam_rec))
  }//while(!atEnd(bamIn))
  return tot_cnt;
}

int main(int argc, char** argv)
{
  if (argc!=3) {
    cout<<"Usage ./multimapper_info genome_bam transcript_bam"<<endl;
  }
  char* genome_bam = argv[1];
  char* transcript_bam = argv[2];
  unordered_set<string> genome_list, transcript_list, common_list;
  cout<<"getting multimappers from genome bam..."<<endl;
  get_multimapper_from_bam(genome_bam, genome_list);
  cout<<"multimapped: "<<genome_list.size()<<endl;
  cout<<"getting multimappers from transcript bam..."<<endl;
  int tot_cnt = get_multimapper_from_bam(transcript_bam, transcript_list);
  cout<<"multimapped: "<<transcript_list.size()<<endl;
  int common(0);
  for (auto r: genome_list)
    if (transcript_list.find(r)!=transcript_list.end())
      ++common;
  cout<<"multimapper caused by alternative splicing: "<<common<<endl;
  cout<<"percentage:"<<endl;
  cout<<"unique: "<<(tot_cnt-transcript_list.size())/double(tot_cnt)<<endl;
  cout<<"rep seq: "<<common/double(tot_cnt)<<endl;
  cout<<"isoform: "<<(transcript_list.size()-common)/double(tot_cnt)<<endl;
  // cout<<"getting uniqmappers from genome bam..."<<endl;
  // get_uniqmapper_from_bam(genome_bam, genome_list);
  // cout<<"uniqmapped: "<<genome_list.size()<<endl;
  // cout<<"getting multimappers from transcript bam..."<<endl;
  // int tot_cnt = get_multimapper_from_bam(transcript_bam, transcript_list);
  // cout<<"multimapped: "<<transcript_list.size()<<endl;
  // int common(0);
  // for (auto r: transcript_list)
  //   if (genome_list.find(r)!=genome_list.end())
  //     ++common;
  // cout<<"multimapper caused by alternative splicing: "<<common<<endl;
  // cout<<"percentage:"<<endl;
  // cout<<"unique: "<<(tot_cnt-transcript_list.size())/double(tot_cnt)<<endl;
  // cout<<"rep seq: "<<(transcript_list.size()-common)/double(tot_cnt)<<endl;
  // cout<<"isoform: "<<common/double(tot_cnt)<<endl;
  return 1;
}
