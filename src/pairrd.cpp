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
#include <complex>
#include <algorithm>
#include <string>
#include <vector>
#include <unistd.h>
using namespace std;

/* samtools header */
#include "samfunctions.h"
#define _cop(c) ((c)&BAM_CIGAR_MASK)
#define _cln(c) ((c)>>BAM_CIGAR_SHIFT)
/**** rsi headers ****/
#include "rsi.h"
/**** user headers ****/
#include "wufunctions.h"
#include "pairrd.h"

struct bamstat_st {   
  int ref;        // ref number
  int l_qseq;        // read length
  double rd;         //read depth		      
  double rd_sd;         //read depth SD	      
  double drd;        //read depth change+-len/4    
  double drd_sd;     //read depth change standard deviation	      
  bool is_paired;
  int isize;
  int isize_sd;
  bamstat_st(): ref(-1),
		l_qseq(-1),
		rd(0),
		drd(0),
		drd_sd(0),
		is_paired(false),
		isize(-1),
		isize_sd(-1) {};
};

struct depths_st {   
  int pos;
  int rd;           //read depth		      
  int drd;          //read depth change+-len/4    
  int m_beg;	    //reads start with M	      
  int s_beg;	    //reads start with S	      
  int m_end;	    //reads end with M	      
  int s_end;	    //reads end with S           
  depths_st():pos(-1),
	      rd(0),
	      drd(0),
	      m_beg(0),
	      s_beg(0), 
	      m_end(0),
	      s_end(0) {};
};


void print_RD_counts(int& pos_printed, int chr, int cur, vector<depths_st>& RDC)
{
  int r_len=100;
  while( pos_printed < cur - r_len/4 ) {
    ++pos_printed;
    
    int dsum0=0;
    int dsum1=0;
    for (int i=max(pos_printed-r_len/4, 0); i<pos_printed; ++i) 
      dsum0+=RDC[i].rd;
    for (int i=pos_printed; i<pos_printed+r_len/4 && i< cur; ++i) 
      dsum1+=RDC[i].rd;
    RDC[pos_printed].drd=(dsum1-dsum0)*4/r_len;
    
    if ( RDC[pos_printed].m_beg==0 && RDC[pos_printed].m_end==0 &&
	 RDC[pos_printed].s_beg==0 && RDC[pos_printed].s_end==0 && 
	 RDC[pos_printed].rd==0 ) continue;
    // cout << b->core.tid+1 << "\t" 
    cout << chr << "\t" 
	 << pos_printed+1 << "\t"
	 << RDC[pos_printed].rd << "\t"
	 << RDC[pos_printed].drd << "\t"
	 << RDC[pos_printed].m_beg << "\t"
	 << RDC[pos_printed].s_beg << "\t"
	 << RDC[pos_printed].m_end << "\t"
	 << RDC[pos_printed].s_end << endl;
  }
  return;
}

void bam_rd_pr_stats(samfile_t *fp_in, bam_index_t *bamidx, 
		     int ref, int beg, int end, bamstat_st& baminfo)
{
  bam1_t *b=NULL; b = bam_init1();
  bam_iter_t iter;
  
  int sample_len=1000000;
  size_t count=0;
  int pos_start=0,pos_end=-10000;
  int offset=0;
  int tid_len=-1;
  int l_qseq=100; double l_qseq_sum=0.0;
  double isize=0.0, isize2=0.0, isize_c=0, isize_sd=0;

  vector<depths_st> RDC(0);
  
  
  // cerr << ref << "\t" << beg << "\t" << end << endl;
  iter = bam_iter_query(bamidx, ref, beg, end);
  tid_len=fp_in->header->target_len[ ref ];
  
  RDC.reserve( sample_len +2000 );
  RDC.resize(  sample_len +1000 );
  
  while(  bam_iter_read(fp_in->x.bam, iter, b)  > 0 ) {
    
    if ( (int)b->core.tid < 0 ) continue;
    if ( b->core.mtid != b->core.tid && b->core.mtid>0 ) continue;
    if ( (b->core.flag & BAM_FSECONDARY)>0 ) continue;
    if ( (b->core.flag & BAM_FDUP)>0 ) continue;
    
    //*! BAM_FPROPER_PAIR may flag long/negative inserts as not properly paired
    if ( (b->core.flag & BAM_FPROPER_PAIR) && (b->core.mtid == b->core.tid) ) {
      //int ISIZE=abs(b->core.isize);
      isize+=abs(b->core.isize);
      isize2+=b->core.isize*b->core.isize;
      isize_c+=1;
    }
    
    int read_pos_end=bam_calend(&b->core, bam1_cigar(b));
    if ( b->core.pos >= tid_len ) break;
    if ( read_pos_end >= tid_len ) break;
    
    // reset and avoid breakage
    if ( b->core.pos > pos_end+1000 ) {
      count=0;
      pos_start=b->core.pos;
      offset=pos_start;
      pos_end=b->core.pos;
      l_qseq_sum=0;
    }
    
    if ( read_pos_end-offset > (int)RDC.size()-100 ) {
      // RDC.reserve( 2*RDC.size() );
      cerr << "expanding" << RDC.size() << endl;
      cerr << offset << "\t" << read_pos_end << endl;
      RDC.insert(RDC.end(), RDC.size(), RDC.back());
    }
    
    pos_end=b->core.pos;
    l_qseq_sum+=b->core.l_qseq;
    
    count++;
    if ( count>1000000 || (pos_end-pos_start)>sample_len ) break;
    
    uint32_t *cigar = bam1_cigar(b);
    int ncigar=b->core.n_cigar;
    if ( ncigar<=1 ) {
      if ( _cop(cigar[0]) == BAM_CMATCH ) {
	++RDC[b->core.pos-offset].m_beg;
	++RDC[read_pos_end-offset].m_end;
      }
    }
    else {
      int op = _cop(cigar[0]);
      if ( op == BAM_CMATCH || 
	   op == BAM_CDEL || 
	   op == BAM_CEQUAL || 
	   op == BAM_CDIFF ) ++RDC[b->core.pos-offset].m_beg;
      else  ++RDC[b->core.pos-offset].s_beg;
      
      op = _cop(cigar[ncigar-1]);
      if ( op == BAM_CMATCH || 
	   op == BAM_CDEL || 
	   op == BAM_CEQUAL || 
	   op == BAM_CDIFF ) ++RDC[read_pos_end-offset].m_end;
      else  ++RDC[read_pos_end-offset].s_end;
    }
    
    POSCIGAR_st bm;
    resolve_cigar_pos(b, bm);  
    for(int k=0;k<(int)bm.op.size();++k) {
      if ( bm.op[k] != BAM_CMATCH && bm.op[k] != BAM_CEQUAL ) continue;
      int p1=bm.cop[k]-1;
      int q1=bm.qop[k];
      for(size_t i=0; i<bm.nop[k] && p1<tid_len; ++i, ++p1, ++q1)      
	if ( bam1_qual(b)[q1]>=0 ) ++RDC[p1-offset].rd;
    }
    
  }
  
  l_qseq=l_qseq_sum/count+0.5;
  double rdmean=0.0;
  double drdmean=0.0;
  for(int k=pos_start+l_qseq/4; k<pos_end-l_qseq/4; ++k) {
    int dsum0=0;
    int dsum1=0;
    for (int i=max(k-l_qseq/4, 0); i<k; ++i) dsum0+=RDC[i-offset].rd;
    for (int i=k; i<k+l_qseq/4; ++i) dsum1+=RDC[i-offset].rd;
    RDC[k-offset].drd=(dsum1-dsum0)*4/l_qseq;
    rdmean+=RDC[k-offset].rd;
    drdmean+=RDC[k-offset].drd;
  }
  rdmean/=(pos_end-pos_start-l_qseq/2);
  drdmean/=(pos_end-pos_start-l_qseq/2);
  
  double drdsd=0.0, rdsd=0.0;
  for(int k=pos_start+l_qseq/4; k<pos_end-l_qseq/4; ++k) {
    drdsd+=(drdmean-RDC[k-offset].drd)*(drdmean-RDC[k-offset].drd);
    rdsd+=(rdmean-RDC[k-offset].rd)*(rdmean-RDC[k-offset].rd);
  }
  drdsd=sqrt(drdsd/(pos_end-pos_start-l_qseq/2));
  rdsd=sqrt(rdsd/(pos_end-pos_start-l_qseq/2));
  
  if ( isize_c>2 ) {
    isize/=isize_c;
    isize_sd =sqrt( (isize2-isize_c*isize*isize)/isize_c );
    baminfo.isize=isize;
    baminfo.isize_sd=isize_sd;
  }
  
  //update main control
  rsi::bam_ref=ref;
  rsi::bam_l_qseq=l_qseq;
  rsi::bam_rd=rdmean;
  rsi::bam_rd_sd=rdsd;
  rsi::bam_drd=drdmean;
  rsi::bam_drd_sd=drdsd;
  
  //return baminfo
  baminfo.ref=ref;
  baminfo.l_qseq=l_qseq;
  baminfo.rd=rdmean;
  baminfo.rd_sd=rdsd;
  baminfo.drd=drdmean;
  baminfo.drd_sd=drdsd;

  return;
}


void get_rds_stats(samfile_t *fp_in, bam_index_t *bamidx, cnv_st icnv, int pend,
		   vector<depths_st>& RDC)
{
  bam_iter_t iter;
  bam1_t *b=NULL; b = bam_init1();   
  int ref=-1;
  int beg=0;
  int end=0;
  // bam_parse_region(fp_in->header, icnv.chr.c_str(), &ref, &beg, &end); 
  ref=icnv.tid;
  if ( ref<0 ) {
    cerr << "chromosome cnv.tid < 0 " << endl;
    exit(0);
  }
  
  size_t count=0;
  int tid_len=-1;
  int r_len=100;
  r_len=max(rsi::m,rsi::bam_l_qseq);
  r_len=max(50,r_len);
  
  if ( pend==0 ) {
    beg=icnv.start-r_len;
    end=icnv.start+r_len;
  }
  else if ( pend==1 ) {
    beg=icnv.end-r_len;
    end=icnv.end+r_len;
  }
  else {
    cerr << "pend==0; first break point\n"
	 << "pend==1; second break point\n"
	 << "pend=" << end << " illegal" << endl;
    exit(0);
  }
  RDC.clear();
  RDC.reserve(2*r_len+2);  // break point +/- r_len
  RDC.resize(2*r_len+1);   // break point +/- r_len
  
  b = bam_init1();
  iter = bam_iter_query(bamidx, ref, beg, end);
  tid_len=fp_in->header->target_len[ ref ];
  // cerr << ref << "\t" << beg << "\t" << end << endl;

  while(  bam_iter_read(fp_in->x.bam, iter, b)  > 0 ) {
    
    if ( (int)b->core.tid < 0 ) continue;
    if ( b->core.mtid != b->core.tid && b->core.mtid>0 ) continue;
    if ( (b->core.flag & BAM_FSECONDARY)>0 ) continue;
    if ( (b->core.flag & BAM_FDUP)>0 ) continue;
    
    int read_pos_end=bam_calend(&b->core, bam1_cigar(b));
    if ( b->core.pos >= tid_len ) continue;
    if ( read_pos_end >= tid_len ) continue;
    
    count++;
    uint32_t *cigar = bam1_cigar(b);
    int ncigar=b->core.n_cigar;
    int ipos=-1;
    if ( ncigar<=1 ) {
      if ( _cop(cigar[0]) == BAM_CMATCH ) {
	ipos=b->core.pos+1-beg;
	if ( ipos>=0 && ipos<=2*r_len ) ++RDC[ipos].m_beg;
	ipos=read_pos_end+1-beg;
	if ( ipos>=0 && ipos<=2*r_len ) ++RDC[ipos].m_end;
      }
    }
    else {
      int op = _cop(cigar[0]);
      if ( op == BAM_CMATCH || 
	   op == BAM_CDEL || 
	   op == BAM_CEQUAL || 
	   op == BAM_CDIFF ) 
	{
	  ipos=b->core.pos+1-beg;
	  if ( ipos>=0 && ipos<=2*r_len ) ++RDC[ipos].m_beg;
	}
      else {
	ipos=b->core.pos+1-beg;
	if ( ipos>=0 && ipos<=2*r_len ) ++RDC[ipos].s_beg;
      }
      
      op = _cop(cigar[ncigar-1]);
      if ( op == BAM_CMATCH || 
	   op == BAM_CDEL || 
	   op == BAM_CEQUAL || 
	   op == BAM_CDIFF ) 
	{
	  ipos=read_pos_end+1-beg;
	  if ( ipos>=0 && ipos<=2*r_len ) ++RDC[ipos].m_end;
	}
      else {
	ipos=read_pos_end+1-beg;
	if ( ipos>=0 && ipos<=2*r_len ) ++RDC[ipos].s_end;
      }
    }
    
    POSCIGAR_st bm;
    resolve_cigar_pos(b, bm);  
    for(int k=0;k<(int)bm.op.size();++k) {
      if ( bm.op[k] != BAM_CMATCH && bm.op[k] != BAM_CEQUAL ) continue;
      int p1=bm.cop[k]-1;
      int q1=bm.qop[k];
      for(size_t i=0; i<bm.nop[k]; ++i, ++p1, ++q1)      
	if ( bam1_qual(b)[q1]>=0 ) {
	  int ipos=p1+1-beg;
	  if ( ipos>=0 && ipos<=2*r_len ) ++RDC[ipos].rd;
	}
    }
    
  }
  
  int rlen=rsi::bam_l_qseq;
  if ( rlen<1 ) {
    cerr << "Error reading rsi::bam_l_qseq " << rlen << endl;
    exit(0);
  }
  for(int k=0+rlen/4; k<(int)RDC.size()-rlen/4;++k) {
    double dsum0=0;
    double dsum1=0;
    for (int i=max(k-rlen/4, 0); i<k; ++i) dsum0+=RDC[i].rd;
    for (int i=k; i<min(k+rlen/4, (int)RDC.size()); ++i) dsum1+=RDC[i].rd;
    RDC[k].drd=(dsum1-dsum0)*4/rlen;
  }
  
  return;
}

// note: cnvs in cnvlist are 1 based
struct rdstats_st {
  int pos;
  int rd;
  int drd;
  int ss;
  int es;
  double score;
  rdstats_st(): pos(-1),
		rd(0),
		drd(0),
		ss(0),
		es(0),
		score(0) {};
} ;

void optimize_resolution(samfile_t* fp_in, bam_index_t *bamidx, 
			 vector<cnv_st>& cnvlist)
{
  if ( cnvlist.size()==0 ) return;

  bamstat_st baminfo,baminfo0;
  double drd_sd=3.5;
  double s_weight=3.0;

  cnv_st icnv;
  rdstats_st bp1, bp2, rdstat0;
  vector<depths_st>* RDC;
  vector<depths_st> RDC0(0);
  vector<depths_st> RDC1(0);
  vector<double> score5(0);
  vector<double> score3(0);
  
  vector<bool> is_keep(cnvlist.size(),false);
  for(int k=0; k<(int)cnvlist.size(); ++k) {
    
    icnv=cnvlist[k];
    // remove short cnvs
    if ( abs(icnv.end-icnv.start) < 800 ) continue;
    
    // sample bam for rd infomation
    int ref=icnv.tid;
    if ( baminfo.ref != ref ) {
      bam_rd_pr_stats(fp_in, bamidx, ref, 10000000, 349250621, baminfo);
      cerr << "sample chr : " << rsi::target_name[ref] << endl;
      cerr << "read length: " << baminfo.l_qseq << endl;
      cerr << "rd mean    : " << baminfo.rd << endl;
      cerr << "rd sd      : " << baminfo.rd_sd << endl;
      cerr << "drd mean   : " << baminfo.drd << endl;
      cerr << "drd sd     : " << baminfo.drd_sd << endl;
      drd_sd=baminfo.drd_sd;
    }
    
    if ( icnv.type==TYPE_UNKNOWN ) {
      cerr << rsi::target_name[icnv.tid] << "\t"
	   << icnv.start << "\t"
	   << icnv.end << "\t"
	   << rsi::cnv_type[icnv.type] << endl;
      cerr << "unknown cnv type : " << icnv.type << "\n"
	   << "keywords are del, loss, dup, gain, add" << endl;
      continue;
    }
    
    // gather info at the both ends
    get_rds_stats(fp_in, bamidx, icnv, 0, RDC0 );
    get_rds_stats(fp_in, bamidx, icnv, 1, RDC1 );
    
    double score_threshold=s_weight*2;
    score5.clear();
    score3.clear();
    score5.resize( (*RDC).size(), 0 );
    score3.resize( (*RDC).size(), 0 );
    
#define _get_score_from_drd_ss_es {					\
      score_threshold=s_weight*2;					\
      score5[i]=0;							\
      score5[i]=-(*RDC)[i].drd/drd_sd+s_weight*(*RDC)[i].s_end;		\
      score3[i]=0;							\
      score3[i]=(*RDC)[i].drd/drd_sd+s_weight*(*RDC)[i].s_beg;		\
    }
    
    // process first break point
    RDC=&RDC0;
    double score_max=0.0;
    int i_score_max=-1;
    bp1=bp2=rdstat0;
    for(int i=0; i<(int)(*RDC).size(); ++i) {
      score5[i]=score3[i]=0;			\
      _get_score_from_drd_ss_es;
      double score=0;
      if ( icnv.type==TYPE_DEL ) score=score5[i];
      if ( icnv.type==TYPE_DUP ) score=score3[i];
      if ( score > score_max ) {
	i_score_max=i;
	score_max=score;
      }
      
      /*
	int ipos=icnv.start-(*RDC).size()/2+i ;
	cout << ipos << "\t" 
	   << (*RDC)[i].rd << "\t"
	   << (*RDC)[i].drd << "\t"
	   << (*RDC)[i].m_beg << "\t"
	   << (*RDC)[i].s_beg << "\t"
	   << (*RDC)[i].m_end << "\t"
	   << (*RDC)[i].s_end << "\t" 
	   << score5[i] << "\t" 
	   << score3[i] 
	   << endl;
      */
    }
    if ( i_score_max >=0 ) {
      bp1.pos=icnv.start-(*RDC).size()/2+i_score_max ;
      bp1.rd=(*RDC)[i_score_max].rd;
      bp1.drd=(*RDC)[i_score_max].drd;
      bp1.ss=(*RDC)[i_score_max].s_beg;
      bp1.es=(*RDC)[i_score_max].s_end;
      bp1.score=score_max;
    }
    if ( icnv.type==TYPE_DEL ) cnvlist[k].sc1=bp1.es;
    if ( icnv.type==TYPE_DUP ) cnvlist[k].sc1=bp1.ss;
    
    // process second break point
    RDC=&RDC1;
    score_max=0;
    i_score_max=-1;
    for(int i=0; i<(int)(*RDC).size(); ++i) {
      _get_score_from_drd_ss_es;
      double score=0;
      if ( icnv.type==TYPE_DEL ) score=score3[i];
      if ( icnv.type==TYPE_DUP ) score=score5[i];
      if ( score > score_max ) {
	i_score_max=i;
	score_max=score;
      }
      
      /*
	int ipos=icnv.end-(*RDC).size()/2+i ;
	cout << ipos << "\t" 
	<< (*RDC)[i].rd << "\t"
	<< (*RDC)[i].drd << "\t"
	   << (*RDC)[i].m_beg << "\t"
	   << (*RDC)[i].s_beg << "\t"
	   << (*RDC)[i].m_end << "\t"
	   << (*RDC)[i].s_end << "\t"
	   << score5[i] << "\t" 
	   << score3[i] 
	   << endl;
      */
    }
    if ( i_score_max >=0 ) {
      bp2.pos=icnv.end-(*RDC).size()/2+i_score_max ;
      bp2.rd=(*RDC)[i_score_max].rd;
      bp2.drd=(*RDC)[i_score_max].drd;
      bp2.ss=(*RDC)[i_score_max].s_beg;
      bp2.es=(*RDC)[i_score_max].s_end;
      bp2.score=score_max;
    }
    if ( icnv.type==TYPE_DEL ) cnvlist[k].sc2=bp2.ss;
    if ( icnv.type==TYPE_DUP ) cnvlist[k].sc2=bp2.es;
    
    if ( cnvlist[k].sc1 >= 2 ) cnvlist[k].start=bp1.pos;
    if ( cnvlist[k].sc2 >= 2 ) cnvlist[k].end=bp2.pos;

    /*
    cerr << icnv.chr << "\t" 
	 << icnv.start << "\t"
	 << icnv.end << "\t" 
	 << rsi::cnv_type[icnv.type] << endl;
    cerr << icnv.chr << "\t" 
	 << bp1.pos << "\t"
	 << bp2.pos << "\t" 
	 << rsi::cnv_type[icnv.type] << endl;
    cerr << bp1.pos << "\t"
	 << bp1.rd << "\t"
	 << bp1.drd << "\t"
	 << bp1.ss << "\t"
	 << bp1.es << "\t" 
	 << bp1.score << endl;
    cerr << bp2.pos << "\t"
	 << bp2.rd << "\t"
	 << bp2.drd << "\t"
	 << bp2.ss << "\t"
	 << bp2.es << "\t" 
	 << bp2.score << endl;
    cerr << "------------------\t" << is_keep[k] << endl;
    */

  }
  
  // The following needs to be done
  
  return;
}


void q0check(samfile_t* fp_in, bam_index_t *bamidx, vector<cnv_st>& cnvlist)
{
  if ( cnvlist.size()==0 ) return;
  
  bam1_t *b=NULL; b = bam_init1();   
  bam_iter_t iter=NULL;
 
  for(int k=0; k<(int)cnvlist.size(); ++k) {
    
    int ref=cnvlist[k].tid;
    int beg=cnvlist[k].start;
    int end=cnvlist[k].end;
    if ( beg>end ) swap(beg, end);
    
    iter = bam_iter_query(bamidx, ref, beg, end);
    
    b = bam_init1();
    double count=0;
    double q0_count=0;
    while(  bam_iter_read(fp_in->x.bam, iter, b)  > 0 ) {
      count+=1;
      if ( b->core.qual==0 ) q0_count+=1;
    }
    
    if ( count>0 ) cnvlist[k].q0=q0_count/count;
    else cnvlist[k].q0=-1;
    
  }
  
  bam_destroy1(b);
  bam_iter_destroy(iter);
  return;
}


void cnv_stat(samfile_t* fp_in, bam_index_t *bamidx, vector<cnv_st>& cnvlist)
{
  if ( cnvlist.size()==0 ) return;
  
  bam1_t *b=NULL; b = bam_init1();
  bam_iter_t iter=NULL;
  bamstat_st baminfo;
  int DIS=1000;
  
  for(int k=0; k<(int)cnvlist.size(); ++k) {
    
    // sample bam for rd infomation
    int ref=cnvlist[k].tid;
    if ( baminfo.ref != ref ) {
      bam_rd_pr_stats(fp_in, bamidx, ref, 10000000, 349250621, baminfo);
      cerr << "sample chr : " << rsi::target_name[ref] << endl;
      cerr << "read length: " << baminfo.l_qseq << endl;
      cerr << "rd mean    : " << baminfo.rd << endl;
      cerr << "rd sd      : " << baminfo.rd_sd << endl;
      cerr << "drd mean   : " << baminfo.drd << endl;
      cerr << "drd sd     : " << baminfo.drd_sd << endl;
      cerr << "ISIZE mean : " << baminfo.isize << endl;
      cerr << "ISIZE sd   : " << baminfo.isize_sd << endl;
    }
    
    int ISIZE_mean=baminfo.isize;
    int ISIZE_std=baminfo.isize_sd;
    
    int beg=cnvlist[k].start;
    int end=cnvlist[k].end;
    int TYPE=cnvlist[k].type;
    if ( beg>end ) swap(beg, end);
    int LEN=end-beg+1;
    DIS=max(DIS,LEN);
    DIS=min(DIS,5000);
    int p1e=beg-DIS;
    int p2e=end+DIS;
    if ( p1e<1 ) p1e=1;
    
    
    iter = bam_iter_query(bamidx, ref, p1e, p2e);
    
    double q0_all_count=0;
    double q0_q0_count=0;
    size_t rp_count=0;
    while(  bam_iter_read(fp_in->x.bam, iter, b)  > 0 ) {
      
      if ( (int)b->core.n_cigar <=1 ) continue;
      if ( (int)b->core.tid < 0 ) continue;
      //if ( NOSECONDARY && (b->core.flag & 0x100)>0 ) continue;
      
      int rbeg = b->core.pos;
      int rend = b->core.n_cigar ? 
	bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + 1;
      
      // check q0 within cnv region
      if ( rend > beg && rbeg < end ) {
	q0_all_count+=1;
	if ( b->core.qual==0 ) q0_q0_count+=1;
      }
      
      // not on same chr
      if ( b->core.mtid != b->core.tid && b->core.mtid>0 ) continue;
      
      int FLAG=b->core.flag;
      // both reads forward
      if ( (FLAG & BAM_FREVERSE)==0 && (FLAG & BAM_FMREVERSE)==0 ) {
	//cerr << "both reads forward\n"
	//   << read << "\n"
	//   << (FLAG & BAM_FREVERSE) << "\t" << (FLAG & BAM_FMREVERSE) << endl;
	continue;
      }
      // both reads reverse
      if ( (FLAG & BAM_FREVERSE)>0 && (FLAG & BAM_FMREVERSE)>0 ) {
	//cerr << "both reads reverse\n"
	//   << read << "\n"
	//   << (FLAG & BAM_FREVERSE) << "\t" << (FLAG & BAM_FMREVERSE) << endl;
	continue;
      }
      
      // RF orientation
      if ( (FLAG & BAM_FREVERSE)>0 && (FLAG & BAM_FMREVERSE)==0 ) {
	//cerr << "RF orientation\n"
	//   << read << "\n"
	//   << (FLAG & BAM_FREVERSE) << "\t" << (FLAG & BAM_FMREVERSE) << endl;
	// continue;
      }
      
      // only FR orientation left
      // inner positions of pair
      int r1= rend;
      int r2= b->core.mpos;
      
      double overlapratio=0.5;
      if ( TYPE==TYPE_DEL ) {
	if ( r2-r1 < ISIZE_mean+ISIZE_std*3 ) continue;
	int overlap=min(r2, end)-max(r1, beg);
	if ( overlap < 0 ) continue; // not overlap
	if ( abs(r1-beg)+abs(r2-end) < ISIZE_mean+ISIZE_std*3 ) {
	  ++rp_count;
	  continue;
	}
	if ( overlap < LEN*overlapratio ) continue;
	if ( overlap < (r2-r1)*overlapratio ) continue;
	++rp_count;
	continue;
      }
      
      if ( TYPE==TYPE_DUP ) {
	if ( r2-r1 > ISIZE_mean-ISIZE_std*3 ) continue;
	if ( abs(r1-beg)+abs(r2-end) < ISIZE_mean+ISIZE_std*3 ) {
	  ++rp_count;
	  continue;
	}
	if ( r1>r2 ) swap(r1,r2);
	int overlap=min(r2, end)-max(r1, beg);
	if ( overlap < LEN*overlapratio ) continue;
	if ( overlap < (r2-r1)*overlapratio ) continue;
	++rp_count;
	continue;
      }
    }
    
    cnvlist[k].q0=q0_q0_count/(q0_all_count+0.00001);
    cnvlist[k].rp=rp_count;
  }
}
