/*************************************************************************
BEGIN OF LICENSE
Copyright (c) Yinghua Wu and Hongzhe Li (rsicnv project).

This program is free software; you can redistribute and/or modify
the codes written by the authors under the terms of the GNU General 
Public License as published by the Free Software Foundation 
(www.fsf.org); either version 2 of the License, or (at your option) 
any later version. 

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; See the GNU General Public License for 
more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
END OF LICENSE

CONTACT: wu_yinghua@hotmail.com; hongzhe@upenn.edu
*************************************************************************/

/**** system headers ****/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <unistd.h>
using namespace std;

/* samtools header */
#include "samfunctions.h"
#include "readref.h"

/**** user headers ****/
#include "wu2.h"
#include "rsi.h"
#include "wufunctions.h"
#include "gccontent.h"
//#include "linasminterface.h"
#include "loaddata.h"
/**** user headers ****/

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/* remove the N regions givn in rsi::noncodelist vector */
void concatenate_data(Array<int>& RD)
{
  if ( rsi::noncodelist.size()==0 ) return;
  if ( RD.size() != RD.capacity() ) {
    cerr << "RD array not packed full\n";
    exit(0);
  }
  
  int dx=0;
  for(int k=0; k<(int)rsi::noncodelist.size(); ++k) {
    dx+=rsi::noncodelist[k].end-rsi::noncodelist[k].start+1;
    for(int i=rsi::noncodelist[k].start; i<=rsi::noncodelist[k].end; ++i) 
      RD[i]=-9;
  }
  
  int idx=0;
  for(int i=0; i<RD.size(); ++i) {
    if ( RD[i] == -9 )  continue;
    RD[idx]=RD[i];
    ++idx;
  }
  RD.reuse(idx);
  
  if ( RD.size()+dx != RD.capacity() ) {
    cerr << "Something wrong concatenating RD array\n"
	 << idx << "\t" << RD.size() << "\t" << dx << "\t" <<  RD.capacity() << endl;
    exit(0);
  } 
  
  if ( rsi::debug ) 
    cerr << "done concatenating RD " <<  RD.size() << "\t" << dx << "\t" 
	 << RD.capacity() << endl;
  
  rsi::start=1;
  rsi::end=RD.size();
  
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/* remove the N regions givn in rsi::noncodelist vector */
void concatenate_data_safe(Array<int>& RD)
{
  if ( rsi::noncodelist.size()==0 ) return;
  
  int dx=0;
  for(int k=0; k<(int)rsi::noncodelist.size(); ++k) {
    dx+=rsi::noncodelist[k].end-rsi::noncodelist[k].start+1;
    for(int i=rsi::noncodelist[k].start; i<=rsi::noncodelist[k].end; ++i) 
      RD[i]=-9;
  }
  
  Array<int> RD2(RD.size()-dx);
  int k=0;
  for(int i=0; i<RD.size(); ++i) {
    if ( RD[i] != -9 ) { RD2[k]=RD[i]; ++k; }
  }
  if ( k!=RD2.size() ) {
    cerr << "something wrong in concatenating array" << endl;
    exit(0);
  }
  
  RD=RD2;
  RD2.resize(0);
  
  rsi::start=1;
  rsi::end=RD.size();
  
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/* patch the N regions givn in rsi::noncodelist vector with mean */
void patch_N_regions(Array<int>& RD)
{
  if ( rsi::noncodelist.size()==0 ) return;
  
  double RDmean=0;
  for(int i=0;i<RD.size();++i) RDmean+=RD[i];
  RDmean/=((double)RD.size()+0.00001);
  
  for(int k=0; k<(int)rsi::noncodelist.size(); ++k) 
    for(int i=rsi::noncodelist[k].start; i<=rsi::noncodelist[k].end; ++i) 
      RD[i]=RDmean;
  
  rsi::start=1;
  rsi::end=RD.size();
  
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void expand_data(Array<int>& RD)
{
  /* add the N regions given in rsi::noncodelist vector */
  if ( rsi::noncodelist.size()==0 ) return;
  
  int dx=0;
  for(int k=0; k<(int)rsi::noncodelist.size(); ++k) 
    dx+=rsi::noncodelist[k].end-rsi::noncodelist[k].start+1;
  
  if ( RD.size()+dx != RD.capacity() ) {
    cerr << "cannot expand RD array\n"
	 << RD.size() << "+" << dx << "!=" << RD.capacity() << endl;
    exit(0);
  }
  int size0=RD.size();

  int idx=RD.size()-1;
  RD.reuse( RD.capacity() );
  int k=rsi::noncodelist.size()-1;
  for(int i=RD.capacity()-1; i>=0; --i) {
    if ( i>rsi::noncodelist.back().end || i<rsi::noncodelist[0].start ) {
      RD[i]=RD[idx];
      --idx;
      continue;
    }
    if ( i>=rsi::noncodelist[k].start && i<=rsi::noncodelist[k].end ) 
      RD[i]=0;
    if ( i<rsi::noncodelist[k].start ) --k;
    if ( i>rsi::noncodelist[k].end && i<rsi::noncodelist[k+1].start ) {
      RD[i]=RD[idx];
      --idx;
      continue;
    }
  }
  
  if ( idx != -1 ) {
    cerr << "Not expanded correctly\n"
	 << idx << "\t" << k << endl;
    exit(0);
  }
  
  if ( rsi::debug ) 
    cerr << "done expanding RD " <<  size0 << "\t" << dx 
	 << "\t" << RD.size() << "\t" << RD.capacity() << endl;
  
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void expand_data_safe(Array<int>& RD)
{
  /* add the N regions given in rsi::noncodelist vector */
  if ( rsi::noncodelist.size()==0 ) return;
  
  
  int dx=0;
  for(int k=0; k<(int)rsi::noncodelist.size(); ++k) 
    dx+=rsi::noncodelist[k].end-rsi::noncodelist[k].start+1;
  
  Array<int> RD2(RD.size()+dx);
  
  dx=0;
  int idx=0, ibeg=0;
  for(int k=0; k<(int)rsi::noncodelist.size(); ++k) {
    for(int i=ibeg+dx; i<rsi::noncodelist[k].start; ++i) {
      RD2[idx]=RD[ibeg];
      ++ibeg;
      ++idx;
    }
    for(int i=rsi::noncodelist[k].start; i<=rsi::noncodelist[k].end; ++i) {
      RD2[idx]=0;
      ++idx;
    }
    dx+=rsi::noncodelist[k].end-rsi::noncodelist[k].start+1;
  }
  for(int i=ibeg; i<RD.size(); ++i, ++idx) RD2[idx]=RD[i];
  
  if ( RD2.size() != idx ) {
    cerr << "something wrong in inserting NONSEQ regions" << endl;
    exit(0);
  }
  
  RD=RD2;
  RD2.resize(0);
  
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void apply_cap(Array<int>& RD)
{
  /* applying a cap for read depth, default is 4 times the mean */
  if ( rsi::cap<=1 ) return;
  rsi::RDmedian = _median(&RD[0], RD.size());
  rsi::dout << "applying cap " << rsi::cap 
	    << " times of mean " << rsi::RDmedian << endl
	    << "cap = " << rsi::cap*rsi::RDmedian << endl;
  for(int i=0;i<RD.size();++i) 
    if ( RD[i]>rsi::RDmedian*rsi::cap ) RD[i]=rsi::RDmedian*rsi::cap;
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void get_noseq_regions(string& FASTA)
{
  vector<int> N_beg(0);
  vector<int> N_end(0);
  get_N_regions(FASTA, N_beg, N_end);    
  
  int dx=max(50, rsi::m/4);    // a save room for N regions
  for(int k=0;k<(int)N_beg.size();++k) {
    N_beg[k]-=dx;
    N_end[k]+=dx;
    if ( N_beg[k]<0 ) N_beg[k]=0;
    if ( N_end[k]>=(int)FASTA.size() ) N_end[k]=FASTA.size()-1;
    for(int i=N_beg[k]; i<=N_end[k]; ++i) FASTA[i]='N';
  }
  
  get_N_regions(FASTA, N_beg, N_end);    
  
  rsi::dout << "#Noseq regions excluded\n";
  rsi::noncodelist.clear();
  for(int k=0;k<(int)N_beg.size();++k) {
    rsi::dout << rsi::target_name[rsi::tid] << "\t" 
	      << N_beg[k] << "\t" << N_end[k] << endl;
    cnv_st icnv;
    icnv.tid=rsi::tid;
    icnv.start=N_beg[k];
    icnv.end=N_end[k];
    rsi::noncodelist.push_back(icnv);
  }
  
  return;
}

// load all chromosome
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void load_data_from_bam(samfile_t *fp_in, bam_index_t *bamidx, int ref,
			Array<int>& RD) 
{
  //* load reference
  string FASTA;
  rsi::chr=fp_in->header->target_name[ref];
  read_fasta(rsi::reffile, rsi::chr, FASTA);
  if ( FASTA.size() != (size_t)fp_in->header->target_len[ref] ) 
    rsi::dout << "reference and target not same size " 
	 << FASTA.size() << "\t" <<  fp_in->header->target_len[ref]
	 << endl;
  
  //* check GC from reference
  Array<bool> GC(FASTA.size());
  GC.assign(false);
  for(int k=0;k<(int)FASTA.size();++k) GC[k]=( FASTA[k]=='G' || FASTA[k]=='C' );
  
  //* rsi::noncodelist updated
  get_noseq_regions(FASTA);
  
  FASTA.clear();
  FASTA.reserve(0);
  
  bam1_t *b=NULL;   
  bam_iter_t iter=0;
  b = bam_init1();
  int beg=0, end=0x7fffffff;
  
  iter = bam_iter_query(bamidx, ref, beg, end);
  RD.resize(fp_in->header->target_len[ref]);
  RD.assign(0);
  rsi::start=1;
  rsi::end=RD.size();
  
  size_t count=0;
  while( bam_iter_read(fp_in->x.bam, iter, b)>0 ) {
    count++;
    if ( count%1000000==0 ) 
      rsi::dout << "#processed " << count/1000000 << "M reads" << '\xd';
    if ( (int)b->core.pos == 0 ) continue;
    if ( (int)b->core.tid < 0 ) continue;
    if ( (int)b->core.qual < rsi::minq ) continue;
    if ( b->core.flag & BAM_FSECONDARY ) continue;
    if ( b->core.flag & BAM_FDUP ) continue;
    //if ( b->core.flag & BAM_FQCFAIL ) continue;
    //if ( b->core.mtid != b->core.tid && b->core.mtid>=0 ) continue;
    
    POSCIGAR_st bm;
    resolve_cigar_pos(b, bm);  
    for(int k=0;k<(int)bm.op.size();++k) {
      if ( bm.op[k] != BAM_CMATCH && bm.op[k] != BAM_CEQUAL ) continue;
      //      for(size_t p1=bm.cop[k]-1; p1<bm.cop[k]-1+bm.nop[k] && (int)p1<RD.size(); ++p1) 
      //	++RD[p1];
      int p1=bm.cop[k]-1;
      int q1=bm.qop[k];
      for(size_t i=0; i<bm.nop[k] && p1<RD.size(); ++i, ++p1, ++q1)      
	if ( bam1_qual(b)[q1]>=rsi::min_baseQ ) ++RD[p1];
    }
  }
  bam_destroy1(b);
  bam_iter_destroy(iter);
  rsi::dout << endl;
  
  if ( rsi::saverd ) {
    string tmps=rsi::outfile+"."+rsi::chr+"_rd";
    rsi::dout << "RD of " << rsi::chr << " is saved to " << tmps << endl;
    write_rd_to_file(RD, tmps);
  }
  
  double rdsum=0.0;
  for(int k=0;k<(int)RD.size();++k) rdsum+=RD[k];
  rdsum/=(double)RD.size();
  rsi::dout << RD.size() << "\t" << rdsum << endl;
  
  if ( rsi::gcadjust ) checkgccontent(RD, GC);
  
  rdsum=0.0;
  for(int k=0;k<(int)RD.size();++k) rdsum+=RD[k];
  rdsum/=(double)RD.size();
  rsi::dout << RD.size() << "\t" << rdsum << endl;
  
  if ( rsi::cap>1 ) apply_cap(RD);
  
  return;
}

// load a region
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void load_data_from_bam(samfile_t *fp_in, bam_index_t *bamidx, 
			int ref, int beg, int end,
			Array<int>& RD) 
{
  
  bam1_t *b=NULL;   
  bam_iter_t iter=0;
  b = bam_init1();
  
  if ( RD.size() != (int)fp_in->header->target_len[ref] ) {
    RD.resize(fp_in->header->target_len[ref]);
    RD.assign(0);
  }
  rsi::start=1;
  rsi::end=RD.size();
  
  if ( beg<0 ) beg=0;
  if ( end>RD.size() ) end=RD.size();
  
  size_t count=0;
  iter = bam_iter_query(bamidx, ref, beg, end);
  while( bam_iter_read(fp_in->x.bam, iter, b)>0 ) {
    count++;
    if ( count%1000000==0 ) 
      rsi::dout << "#processed " << count/1000000 << "M reads" << '\xd';
    if ( (int)b->core.pos == 0 ) continue;
    if ( (int)b->core.tid < 0 ) continue;
    if ( (int)b->core.qual < rsi::minq ) continue;
    if ( b->core.flag & BAM_FSECONDARY ) continue;
    if ( b->core.flag & BAM_FDUP ) continue;
    //if ( b->core.flag & BAM_FQCFAIL ) continue;
    //if ( b->core.mtid != b->core.tid && b->core.mtid>=0 ) continue;
    
    POSCIGAR_st bm;
    resolve_cigar_pos(b, bm);  
    for(int k=0;k<(int)bm.op.size();++k) {
      if ( bm.op[k] != BAM_CMATCH && bm.op[k] != BAM_CEQUAL ) continue;
      int p1=bm.cop[k]-1;
      int q1=bm.qop[k];
      for(size_t i=0; i<bm.nop[k] && p1<RD.size(); ++i, ++p1, ++q1)      
	if ( bam1_qual(b)[q1]>=rsi::min_baseQ ) ++RD[p1];
    }
  }
  bam_destroy1(b);
  bam_iter_destroy(iter);
  
  double rdmean=0.0;
  for(int k=beg; k<end; ++k) rdmean+=RD[k];
  rdmean/=(double)(end-beg+0.000000001);
  rsi::dout << "Loaded length\t" << end-beg << "\t" << rdmean << endl;
  
  // check GC content only when loading whole chromosome
  if ( rsi::gcadjust && end-beg>=(int)fp_in->header->target_len[ref] ) { 
    string FASTA;
    rsi::chr=fp_in->header->target_name[ref];
    read_fasta(rsi::reffile, rsi::chr, FASTA);
    if ( FASTA.size() != (size_t)fp_in->header->target_len[ref] ) 
      rsi::dout << "reference and target not same size " 
		<< FASTA.size() << "\t" <<  fp_in->header->target_len[ref]
		<< endl;
    
    //* check GC from reference
    Array<bool> GC(FASTA.size());
    GC.assign(false);
    for(int k=0;k<(int)FASTA.size();++k) GC[k]=( FASTA[k]=='G' || FASTA[k]=='C' );
    
    //* rsi::noncodelist updated
    get_noseq_regions(FASTA);
    
    FASTA.clear();
    FASTA.reserve(0);
    
    rsi::dout << endl;
    
    double rdmean=0.0;
    int len=0;
    for(int k=beg;k<(int)RD.size() && k<end; ++k,++len) rdmean+=RD[k];
    rdmean/=(double)len;
    rsi::dout << "length\t" << len << "\t" << rdmean << endl;
    
    if ( rsi::gcadjust ) checkgccontent(RD, GC);
    
    rdmean=0.0;
    len=0;
    for(int k=beg;k<(int)RD.size() && k<end; ++k,++len) rdmean+=RD[k];
    rdmean/=(double)len;
    rsi::dout << "length\t" << len << "\t" << rdmean << endl;
  }
  else {
    rsi::dout << "Note: no GC adjustment" << endl;
  }
  
  if ( rsi::cap>1 && end-beg>=(int)fp_in->header->target_len[ref] ) 
    apply_cap(RD);
  
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void write_rd_to_file(Array<int>& RD, string inpfile) 
{
  ofstream FOUT(inpfile.c_str());
  for(int i=0; i<RD.size(); ++i) FOUT << i+1 << "\t" << RD[i] << "\n";
  FOUT.close();
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void load_data_from_text(string inpfile, Array<int>& RD) 
{
  if ( rsi::chr == "" ) {
    rsi::dout << "chromosome name not set\nexit" << endl;
    exit(0);
  }
  string FASTA;
  read_fasta(rsi::reffile, rsi::chr, FASTA);
  Array<bool> GC(FASTA.size());
  GC.assign(false);
  for(int k=0;k<(int)FASTA.size();++k) GC[k]=( FASTA[k]=='G' || FASTA[k]=='C' );
  
  //* rsi::noncodelist updated
  get_noseq_regions(FASTA);
  
  FASTA.clear();
  FASTA.reserve(0);
  
  RD.resize(GC.size());
  RD.assign(0);
  rsi::start=1;
  rsi::end=RD.size();
  
  ifstream inp(inpfile.c_str());
  if ( !inp ) {
    cerr << "Cannot open file " << inpfile << endl;
    exit(0);
  }
  inp.clear();
  inp.seekg(0);
  while ( !inp.eof() ) {
    string tmps;
    getline(inp,tmps);
    if (tmps.length()< 1) continue;
    if (tmps[0] == '#' ) continue;
    
    istringstream iss(tmps);
    int pos, d;
    iss >> pos >> d;
    
    if ( pos<1 ) continue;
    if ( pos>=RD.size() ) break;
    RD[pos-1]=d;
  }
  inp.close();
  
  double rdsum=0.0;
  for(int k=0;k<(int)RD.size();++k) rdsum+=RD[k];
  rdsum/=(double)RD.size();
  rsi::dout << RD.size() << "\t" << rdsum << endl;
  
  if ( rsi::gcadjust ) checkgccontent(RD, GC);
  
  rdsum=0.0;
  for(int k=0;k<(int)RD.size();++k) rdsum+=RD[k];
  rdsum/=(double)RD.size();
  rsi::dout << RD.size() << "\t" << rdsum << endl;
  
  if ( rsi::cap>1 ) apply_cap(RD);
  
  rdsum=0.0;
  for(int k=0;k<(int)RD.size();++k) rdsum+=RD[k];
  rdsum/=(double)RD.size();
  rsi::dout << RD.size() << "\t" << rdsum << endl;
  
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void read_cnvlist(string inpfile, vector<cnv_st>& cnvlist)
{
  cnvlist.clear();
  cnv_st icnv;

  int i1,i2;
  string tmps,chr,tmps1;
  
  ifstream INP(inpfile.c_str());
  if ( !INP ) {
    rsi::dout << "File " << inpfile << " not exist" << endl;
    exit(0);
  }
  
  i1=i2=0;
  while ( !INP.eof() ) {
    vector<string> cell;

    string tmps;
    getline(INP,tmps);
    istringstream iss(tmps);
    string tmps1;
    
    if (tmps.length() < 5) continue;
    if (tmps[0] == '#' ) continue;
    
    cell.clear();
    while ( iss >> tmps1 ) if ( tmps1.length()>0 ) cell.push_back(tmps1);
    
    if ( cell.size() < 4 ) {
      cerr << "skipping lines with less than 4 fields\n"
	   << "expecting RNAME cnv_begin_pos cnv_end_pos cnv_type in each line" << endl;
      continue;
    }

    std::vector<string>::iterator it;    
    it=find(rsi::target_name.begin(), rsi::target_name.end(), cell[0]);
    if ( it!=rsi::target_name.end() ) 
      icnv.tid=std::distance(rsi::target_name.begin(), it);
    else icnv.tid=-1;
    icnv.start=atoi(cell[1].c_str());
    icnv.end=atoi(cell[2].c_str());
    icnv.type=TYPE_UNKNOWN;
    if ( ci_find(cell[3],"del")!=string::npos ) icnv.type=TYPE_DEL;
    if ( ci_find(cell[3],"loss")!=string::npos ) icnv.type=TYPE_DEL;
    if ( ci_find(cell[3],"dup")!=string::npos ) icnv.type=TYPE_DUP;
    if ( ci_find(cell[3],"gain")!=string::npos ) icnv.type=TYPE_DUP;
    if ( ci_find(cell[3],"add")!=string::npos ) icnv.type=TYPE_DUP;
    cnvlist.push_back(icnv);
  }
  INP.close();
  
  if ( rsi::debug )
    cerr << "#Found " << cnvlist.size() << " CNVs in " << inpfile << endl;
  
  if ( cnvlist.size() == 0 ) {
    cerr << "I can't find any CNV from " << inpfile << endl;
    exit(0);
  }
  
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void read_cnvlist(string inpfile, vector<cnv_st>& cnvlist, 
		  string & text, vector<size_t>& text_idx)
{
  cnvlist.clear();
  text="";
  text_idx.clear();
  cnv_st icnv;
  
  int i1,i2;
  string tmps,chr,tmps1;
  
  ifstream INP(inpfile.c_str());
  if ( !INP ) {
    rsi::dout << "File " << inpfile << " not exist" << endl;
    exit(0);
  }
  
  i1=i2=0;
  while ( !INP.eof() ) {
    vector<string> cell;

    string tmps;
    getline(INP,tmps);
    istringstream iss(tmps);
    string tmps1;
    
    if (tmps.length() < 5) continue;
    if (tmps[0] == '#' ) continue;
    
    cell.clear();
    while ( iss >> tmps1 ) if ( tmps1.length()>0 ) cell.push_back(tmps1);
    
    if ( cell.size() < 4 ) {
      cerr << "skipping lines with less than 4 fields\n"
	   << "expecting RNAME cnv_begin_pos cnv_end_pos cnv_type in each line" << endl;
      continue;
    }
    
    text_idx.push_back( text.size() );
    text+=tmps+"\n";

    std::vector<string>::iterator it;    
    it=find(rsi::target_name.begin(), rsi::target_name.end(), cell[0]);
    if ( it!=rsi::target_name.end() ) 
      icnv.tid=std::distance(rsi::target_name.begin(), it);
    else icnv.tid=-1;
    icnv.start=atoi(cell[1].c_str());
    icnv.end=atoi(cell[2].c_str());
    icnv.type=TYPE_UNKNOWN;
    if ( ci_find(cell[3],"del")!=string::npos ) icnv.type=TYPE_DEL;
    if ( ci_find(cell[3],"loss")!=string::npos ) icnv.type=TYPE_DEL;
    if ( ci_find(cell[3],"dup")!=string::npos ) icnv.type=TYPE_DUP;
    if ( ci_find(cell[3],"gain")!=string::npos ) icnv.type=TYPE_DUP;
    if ( ci_find(cell[3],"add")!=string::npos ) icnv.type=TYPE_DUP;
    cnvlist.push_back(icnv);
  }
  INP.close();
  
  if ( rsi::debug )
    cerr << "#Found " << cnvlist.size() << " CNVs in " << inpfile << endl;
  
  if ( cnvlist.size() == 0 ) {
    cerr << "I can't find any CNV from " << inpfile << endl;
    exit(0);
  }
  
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

