#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <numeric>

#include <seqan/bam_io.h>

using namespace std;

bool get_map_cnt_from_bam(const char *fn, vector<int>& vec)
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
       findTagKey(i, tag, seqan::CharString("NH"));
       int val(1);
       extractTagValue(val, tag, i);
       vec.push_back(val);
    }//if(!hasFlagUnmapped(bam_rec))
  }//while(!atEnd(bamIn))
  return false;
}

bool get_mismatch_cnt_from_bam(const char *fn, vector<int>& vec)
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
       unsigned i(1);
       findTagKey(i, tag, seqan::CharString("NM"));
       int val(1);
       extractTagValue(val, tag, i);
       vec.push_back(val);
    }//if(!hasFlagUnmapped(bam_rec))
  }//while(!atEnd(bamIn))
  return false;
}

bool get_read_len_from_bam(const char *fn, vector<int>& vec)
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
       // read length is the length on the transcript covered by the read
       int read_len(0);
       // parse cigar string to figure out read length
       bool skip(false);
       for (int i=0; i<length(bam_rec.cigar); ++i) {
	 switch (bam_rec.cigar[i].operation) {
	 case 'M':
	   // add match length to read length
	   read_len += bam_rec.cigar[i].count;
	   break;
	 case 'D':
	   // Deletion is a gap to the read
	   // read covered one base longer on the transcript
	   read_len += 1;
	   break;
	 case 'N':
	   // skipped regions indicate novel splices
	   // since reads are mapped to the transcriptome
	   // skip these alignments for now
	   skip = true;
	   break;
	 case 'I':
	   // Insertion is a gap to the reference
	   // no change in read length
	 case 'S':
	   // softclipping won't change read length
	 default:
	   break;
	 }
	 if (skip) break;
       }
       vec.push_back(read_len);
    }//if(!hasFlagUnmapped(bam_rec))
  }//while(!atEnd(bamIn))
  return false;
}

int main(int argc, char** argv)
{
  if (argc!=3) {
    cout<<"Usage ./bam_info input_bam output_dir"<<endl;
  }
  char* bam_fn = argv[1];
  char* odir = argv[2];
  vector<int> vec;
  cout<<"getting map_cnt statistics..."<<endl;
  get_map_cnt_from_bam(bam_fn, vec);
  cout<<"writing map counts to file..."<<endl;
  cout<<string(odir)+"map_cnt.txt"<<" "<<vec.size()<<endl;
  ofstream ofile(string(odir)+"map_cnt.txt");
  for (auto v: vec)
    ofile<<v<<"\t";
  ofile<<endl;
  ofile.close();
  cout<<"getting mismatch statistics..."<<endl;
  vec.clear();
  get_mismatch_cnt_from_bam(bam_fn, vec);  
  long nmismatch=accumulate(vec.begin(), vec.end(), long(0));
  cout<<"nmismatch: "<<nmismatch<<" "<<vec.size()<<endl;
  cout<<"getting read length statistics..."<<endl;
  vec.clear();
  get_read_len_from_bam(bam_fn, vec);
  long read_len=accumulate(vec.begin(), vec.end(), long(0));
  cout<<"error rate for mapped reads: "<<double(nmismatch)/read_len<<endl;
  cout<<"writing read length to file.."<<endl;
  cout<<string(odir)+"read_len.txt"<<" "<<vec.size()<<endl;
  ofile.open(string(odir)+"read_len.txt");
  for (auto v: vec)
    ofile<<v<<"\t";
  ofile<<endl;
  ofile.close();
  return 1;
}
