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
//#include <sys/types.h>
//#include <cstdlib>
using namespace std;

/* samtools header */
#include "samfunctions.h"
#include "readref.h"

/* mathlib header */
#include "alglibinterface.h"
//#include "linasminterface.h"

/**** user headers ****/
#include "pstream.h"
#include "wu2.h"
#include "rsi.h"
#include "wufunctions.h"
#include "gccontent.h"
#include "loaddata.h"
#include "plotcnv.h"
#include "pairrd.h"
/**** user headers ****/

bool rsi::saverd=false;
bool rsi::debug=false;       // print more output
bool rsi::histonly=false;    // check histogram
bool rsi::overlaponly=false; // check overlap between two set of files
bool rsi::combine=false;     // combine cnvs from a set of files
bool rsi::nocodeonly=false;
bool rsi::regiononly=false;  
bool rsi::plot=true;    // plot cnvs given a cnv file
string rsi::plotfolder="cnv_plots";  // plot folder
bool rsi::genotype=false;    // genotype or not
bool rsi::merge=true;        // is mering needed?
bool rsi::gcadjust=true;    // adjust rd according to gc?
int rsi::pid=0;              // process id
long rsi::seed=-19104;       // seed for random number
int rsi::start=0;            // pileup start position
int rsi::end=0;              // pileup end position
int rsi::m=101;              // median block size
string rsi::trans="NBN";     // which transformation to use
double rsi::minmlen=3.01;     // minimal number of blocks considered cnv
double rsi::chklen=2.5;      // when testing cnv, check before and after 
int rsi::maxchkbp=100000;    // maximum points used to test cnv  
int rsi::L=-1;               // if given and greater than 0, maximum length checked
int rsi::Lmax=-1;            // if greater than 0, maximum length checked
int rsi::minq=0;            // if greater than 0, maximum length checked
int rsi::min_baseQ=13;            // if greater than 0, maximum length checked
double rsi::cap=4.0;         // if readdepth > 4.0*mean(RD) RD=4.0*mean(RD)
double rsi::p=0.05;
double rsi::ci=0.95;
double rsi::buffer=0.05;
double rsi::RDmedian=0;
double rsi::RDsd=0;
double rsi::RDsigma=0.0;
double rsi::median=0;
double rsi::medmedian=0;
double rsi::nbnmedian=0;
double rsi::sigma=0.0;
double rsi::lamda=0.0;
double rsi::medlamda=0.0;
double rsi::nbnlamda=0.0;
double rsi::threshold=-1.0;  // if given and >0, lamda*=threshold
double rsi::epsilon=1.5;     // if given, lamda*=sqrt(1.0+e)
double rsi::factor=6.6;      // lamda=factor*sigma
int rsi::hetnorm=1000;
ofstream rsi::fout;
string rsi::inpfile="";    // input read depths file
string rsi::reffile="";    // list of known cnv
string rsi::rdfile="";    // list of known cnv
string rsi::bamfile="";    // list of known cnv
string rsi::outfile="rsiout.txt";    // output file
string rsi::cnvfile="";    // list of known cnv
string rsi::logfile="";
string rsi::function="rsi";  // rsi or plot
dual_stream rsi::dout(cerr,rsi::fout);
string rsi::pidfile="";    // process id file for *nix system
int rsi::bam_ref=-1;
int rsi::bam_l_qseq=0;
double rsi::bam_rd=0;
double rsi::bam_rd_sd=0;
double rsi::bam_drd=0;
double rsi::bam_drd_sd=0;
string rsi::cnv_type[]={"DEL","DUP","UNKNOWN"};
string rsi::chr="1-22XY";          // chromosome name
int rsi::tid=-1;
vector<string> rsi::target_name(0);
vector<cnv_st> rsi::noncodelist(0);   // noncoding regions, 0 based; RD is rsi::start based

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void isitcnv(Array<int>& RDref, Array<int>& RDcnv, cnv_st& icnv)
{
  int i;
  double ref_runmean_med,ref_runmean_sd;
  
  int d=RDcnv.size();
  Array<float> RDrefrunmean(RDref.size()-d);
  
  double sum=0,sum2=0,mn=0,var=0;
  for(i=0;i<d;++i) {
    sum+=RDref[i];
    sum2+=RDref[i]*RDref[i];
  }
  mn=sum/double(d);
  var=sum2/double(d)-mn*mn;
  
  RDrefrunmean[0]=mn;
  for(i=1;i<RDref.size()-d;++i) {
    sum=sum-RDref[i-1]+RDref[i-1+d];
    sum2=sum2-RDref[i-1]*RDref[i-1]+RDref[i-1+d]*RDref[i-1+d];
    mn=sum/double(d);
    RDrefrunmean[i]=mn;
  }
  
  ref_runmean_med=_median(&RDrefrunmean[0],RDrefrunmean.size());
  ref_runmean_sd=sqrt(variance(RDrefrunmean,0,RDrefrunmean.size()-1,0.0,-1));
  // ref_runmean_sd=_interquartilerange(&RDrefrunmean[0],RDrefrunmean.size())/1.349;
  if ( ref_runmean_sd<1E-3 ) ref_runmean_sd=ref_runmean_med/40.0+1E-3;
  
  icnv.length=icnv.end-icnv.start+1;
  icnv.cnvmed=_median( &RDcnv[0],RDcnv.size() );
  icnv.cnvsd=sqrt(variance(RDcnv,0,RDcnv.size()-1,0.0,-1));
  // icnv.cnvsd=_interquartilerange(&RDcnv[0], RDcnv.size())/1.349;
  icnv.cnviqr=_interquartilerange(&RDcnv[0], RDcnv.size());
  icnv.refmed=ref_runmean_med;
  // icnv.refsd=ref_runmean_sd;
  icnv.refsd=_interquartilerange(&RDrefrunmean[0], RDrefrunmean.size())/1.349;
  icnv.refiqr=_interquartilerange(&RDrefrunmean[0], RDrefrunmean.size());
  icnv.geno=1; 
  icnv.status=1;
  
  int flag = icnv.cnvmed>rsi::RDmedian ? TYPE_DUP : TYPE_DEL;
  if ( icnv.type==TYPE_UNKNOWN ) icnv.type=flag;
  if ( icnv.type!=flag ) {
    rsi::dout << "basic assignment error\n"
	      << "RD=" << icnv.cnvmed << "\t" 
	      << "REF=" << icnv.refmed << "\t" 
	      << "CHRMED=" << rsi::RDmedian << "\t" 
	      << rsi::cnv_type[icnv.type] << "\n"
	      << icnv.start << "\t" << icnv.end << "\t" << icnv.end - icnv.start << endl;
    icnv.status=-9 ; return;
  }
  
  double reference=rsi::RDmedian;
  if ( icnv.type==TYPE_DEL )  {
    reference=min(ref_runmean_med, rsi::RDmedian);
    reference=max(reference, 0.8*rsi::RDmedian);
    double nu=(double)(3.0*icnv.cnvmed-2.0*reference)/ref_runmean_sd;
    icnv.p1=alglib::pnorm( nu );
    if ( nu>0 ) { icnv.status=-9 ; icnv.geno=0; }
  }
  if ( icnv.type==TYPE_DUP )  {
    reference=max(ref_runmean_med, rsi::RDmedian);
    double nu=(double)(2.5*icnv.cnvmed-3.0*reference)/ref_runmean_sd/1.5;
    icnv.p1=1.0-alglib::pnorm( nu );
    if ( nu<0 ) { icnv.status=-9 ; icnv.geno=0; }
  }
  if ( icnv.type==TYPE_UNKNOWN )  
    rsi::dout << TYPE_UNKNOWN  << "\t" << icnv.start << "\t" << icnv.end << endl;
  
  return;
} 

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void isitcnvwrap(Array<int>& RD, vector<cnv_st>& cnvlist, int cnvidx)
{
  int flag=cnvlist[cnvidx].type;
  if ( flag != TYPE_DUP && flag !=TYPE_DEL ) {
    rsi::dout << "flag = TYPE_DUP | TYPE_DEL " << endl
	      << "isitcnvwrap()" << endl;
    exit(0);
  }
  
  /* if test cnv is too long, interpolate it */
  int cnvlen=cnvlist[cnvidx].end-cnvlist[cnvidx].start+1;
  int pts=rsi::maxchkbp*10;
  
  int i,j,k;
  int d;
  int idx;
  
  d=cnvlist[cnvidx].end-cnvlist[cnvidx].start+1;
  if ( RD.size() == rsi::end-rsi::start+1 ) {
    if ( d < rsi::m*rsi::minmlen ) d=rsi::m*rsi::minmlen ;
  }
  if ( RD.size() < (rsi::end-rsi::start+1)/2 ) {
    if ( d < rsi::minmlen ) d=(int)rsi::minmlen+1 ;
  }
  
  Array<int> RDref(rsi::chklen*d*2,0);  
  
  int buffer = int(cnvlen*rsi::buffer+1);
  i=cnvlist[cnvidx].start - buffer;
  idx=cnvidx-1;
  while ( i>0 && idx>0 && i<cnvlist[idx].start ) --idx ;
  while( idx>0 && cnvlist[idx].status==-9 ) --idx;
  j=0;
  k=rsi::chklen*d-1;
  if ( RD.size()-cnvlist[cnvidx].end < rsi::chklen*d )
    k=RDref.size()-1-RD.size()+cnvlist[cnvidx].end;
  int stopper=k;
  double upper=3.0;
  double lower=0.15;
  while( i>2 && k>=0 ) {
    --i;
    if ( flag ==TYPE_DEL && RD[i] > rsi::RDmedian*upper ) continue;
    if ( flag ==TYPE_DUP && RD[i] < rsi::RDmedian*lower ) continue;
    if ( idx >= 0 ) {
      if ( i>=cnvlist[idx].start && i<=cnvlist[idx].end ) {
	i=cnvlist[idx].start-1;
	--idx;
	while( idx>0 && cnvlist[idx].status==-9 ) --idx;
	continue;
      }
    }
    RDref[k]=RD[i];
    --k;
  }
  if ( k>=0 ) {
    int k1,i1;
    for(k1=k+1,i1=0;k1<=stopper;++k1,++i1) RDref[i1]=RDref[k1];
    k=i1;
  }
  else k=stopper+1;
  
  int nghidx=k-1; 
  
  i=cnvlist[cnvidx].end + buffer;
  idx=cnvidx+1;
  while ( i<RD.size()-2 && idx<(int)cnvlist.size() && i>cnvlist[idx].end ) ++idx ;
  while ( idx<(int)cnvlist.size()-1 && cnvlist[idx].status==-9 ) ++idx;
  while( i<RD.size()-2 && k<2*rsi::chklen*d ) {
    ++i;
    if ( flag ==TYPE_DEL && RD[i] > rsi::RDmedian*upper ) continue;
    if ( flag ==TYPE_DUP && RD[i] < rsi::RDmedian*lower ) continue;
    if ( idx < (int)cnvlist.size() ) {
      if ( i>=cnvlist[idx].start && i<=cnvlist[idx].end ) {
	i=cnvlist[idx].end+1;
	++idx;
	while ( idx<(int)cnvlist.size()-1 && cnvlist[idx].status==-9 ) ++idx;
	continue;
      }
    }
    RDref[k]=RD[i];
    ++k;
  }
  
  if ( k<RDref.size() ) RDref.reuse(k);
  
  Array<int> RDcnv(cnvlist[cnvidx].end-cnvlist[cnvidx].start+1);
  for(i=cnvlist[cnvidx].start,k=0;i<=cnvlist[cnvidx].end;++i,++k) 
    RDcnv[k]=RD[i];
  
  d=RDref.size()+RDcnv.size();
  if ( d>pts ) {
    int dref=0,dcnv=0;
    dref=(double)RDref.size()/(double)d*(double)pts;
    dcnv=(double)RDcnv.size()/(double)d*(double)pts;
    Array<int> RD1(dref);
    for(i=0;i<dref;++i) {
      j=double(i)/double(dref)*double(RDref.size());
      RD1[i]=RDref[j];
    }
    nghidx=double(nghidx)/double(RDref.size())*double(dref);
    RDref=RD1;
    RD1.resize(dcnv);
    for(i=0;i<dcnv;++i) {
      j=double(i)/double(dcnv)*double(RDcnv.size());
      RD1[i]=RDcnv[j];
    }
    RDcnv=RD1;
  }
  
  isitcnv(RDref, RDcnv, cnvlist[cnvidx] );
  
  return;
} 


/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void get_continuous_segments(Array<int>& RDstatus, int d, 
			     vector<cnv_st>& seglist)  
{
  seglist.clear();
  
  cnv_st iseg;
  int i;
  int istart=0;
  int iend=0;
  int icount=0;
  
  for(i=0;i<RDstatus.size();++i) {
    if ( RDstatus[i]!=0 ) {
      if ( icount == 0 ) {
	istart=i ;
	iend=i;
        icount++;
	continue;
      }

      if ( (double)RDstatus[i]*(double)RDstatus[iend] >0 
	   && (i-iend)<=d ) {
	iend=i; 
	continue; 
      }
      
      iseg.start=istart;
      iseg.end=iend;
      seglist.push_back(iseg);
      
      istart=i;
      iend=i;
      icount++;
    }
  }
  
}
void get_continuous_segments(Array<char>& RDstatus, int d, 
			     vector<cnv_st>& seglist)  
{
  seglist.clear();
  
  cnv_st iseg;
  int i;
  int istart=0;
  int iend=0;
  int icount=0;
  
  for(i=0;i<RDstatus.size();++i) {
    if ( RDstatus[i]!=0 ) {
      if ( icount == 0 ) {
	istart=i ;
	iend=i;
        icount++;
	continue;
      }

      if ( (double)RDstatus[i]*(double)RDstatus[iend] >0 
	   && (i-iend)<=d ) {
	iend=i; 
	continue; 
      }
      
      iseg.start=istart;
      iseg.end=iend;
      seglist.push_back(iseg);
      
      istart=i;
      iend=i;
      icount++;
    }
  }
  
}


/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void multisegments(cnv_st& iseg, Array<int>& RDmedmean_status, 
		   vector<cnv_st>& mseg)
{
  mseg.clear();
  
  cnv_st icnv;
  
  int i,k;
  
  Array<int> status2(iseg.end-iseg.start+1,0);
  for(i=iseg.start,k=0;i<=iseg.end;++i,++k) status2[k]=RDmedmean_status[i];
  Array<int> status(iseg.end-iseg.start+1,0);
  
  int minlevel=min(status2, 0, status2.size()-1);
  int maxlevel=max(status2, 0, status2.size()-1);
  
  for(int ilevel=minlevel; ilevel<maxlevel ; ++ilevel) {
    if ( ilevel == 0 ) continue;
    
    status.assign(0);
    int levelcount=0;
    for(i=0;i<status2.size();++i) {
      if ( status2[i] == 0 ) { status[i]=0; continue; }
      if ( status2[i] == ilevel ) ++levelcount;
      if ( ilevel<0 && status2[i]<0 && status2[i] >= ilevel) status[i]=1;
      if ( ilevel>0 && status2[i]>0 && status2[i] <= ilevel) status[i]=1;
    }
    if ( levelcount==0 ) continue;
    
    vector<cnv_st> subseg;
    get_continuous_segments(status, 1, subseg);
    for(i=0;i<(int)subseg.size();++i) {
      subseg[i].start+=iseg.start;
      subseg[i].end+=iseg.start;;
      if ( subseg[i].start < iseg.start ) subseg[i].start=iseg.start;
      if ( subseg[i].end > iseg.end ) subseg[i].end=iseg.end;
    }
    
    mseg.insert(mseg.end(),subseg.begin(),subseg.end());
  }
  
  return;
}


/* test if cnv with localmedian transfermation data */ 
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void areblockscnv(Array<int>& RDmedint,Array<int>& RDtrans_status,
		  vector<cnv_st>& rsiseglist)
{
  int flag;
  int i,j;
  vector<cnv_st> testseglist;
  testseglist=rsiseglist;
  
  for(i=0;i<(int)testseglist.size();++i) {
    cnv_st iseg;
    iseg=testseglist[i];
    flag = RDtrans_status[iseg.start]>0 ? TYPE_DUP : TYPE_DEL; 
    if ( flag != iseg.type ) {
      rsi::dout << "error : some thing wrong with type" << endl
		<< iseg.start << "\t"
		<< RDmedint[iseg.start] << "\t" 
		<< (iseg.start+iseg.end)/2 << "\t"
		<< RDmedint[(iseg.start+iseg.end)/2] << "\t" 
		<< iseg.end << "\t"
		<< RDmedint[iseg.end] << "\t" 
		<< rsi::cnv_type[iseg.type] << endl 
		<< "areblockscnv()" << endl;
    }
    
    isitcnvwrap(RDmedint, testseglist, i);
    if ( rsi::debug) 
      rsi::dout << i+1 << "\t" 
		<< testseglist[i].start << "\t" << testseglist[i].end << "\t"
		<< testseglist[i].length << "\t" 
		<< rsi::cnv_type[testseglist[i].type] << "\t"
		<< testseglist[i].geno << "\t"
		<< testseglist[i].p1 << "\t"
		<< testseglist[i].p2 << "\t"
		<< testseglist[i].cnvmed << "\t"
		<< testseglist[i].refmed << "\t"
		<< "rsi"
		<< endl;
  }
  
  /* for not cnvs test subset with localmedian transfermation data */ 
  for(i=0;i<(int)testseglist.size();++i) {
    if ( testseglist[i].status != -9 ) continue;
    
    if ( testseglist[i].type == TYPE_DEL &&
	 testseglist[i].cnvmed < 0.7*testseglist[i].refmed ) {
      testseglist[i].geno = 1;
      testseglist[i].p1 = rsi::p;
      continue;
    }
    if ( testseglist[i].type == TYPE_DUP &&
	 testseglist[i].cnvmed > 1.3*testseglist[i].refmed ) {
      testseglist[i].geno = 1;
      testseglist[i].p1 = rsi::p;
      continue;
    }
    
    if ( rsi::debug ) {
      rsi::dout << i+1 << "\t" 
		<< testseglist[i].start << "\t" 
		<< testseglist[i].end << "\t" 
		<< testseglist[i].length << "\t" 
                << rsi::cnv_type[testseglist[i].type] << "\t"
		<< testseglist[i].geno << "\t"
                << testseglist[i].cnvmed << "\t"
		<< testseglist[i].refmed << "\t"
		<< endl;
    }
    cnv_st iseg,oseg;
    oseg=testseglist[i];
    iseg=testseglist[i];
    flag=iseg.type;
    
    vector<cnv_st> mseg;
    mseg.clear();
    multisegments(iseg,RDtrans_status,mseg);
    for(j=mseg.size()-1;j>=0;--j) {
      mseg[j].type=iseg.type;
      testseglist[i]=mseg[j];
      isitcnvwrap(RDmedint, testseglist, i);
      mseg[j]=testseglist[i];
      if ( rsi::debug) 
	rsi::dout << "\t" 
		  << mseg[j].start << "\t" 
		  << mseg[j].end << "\t" 
		  << mseg[j].end-mseg[j].start+1 << "\t" 
		  << RDtrans_status[mseg[j].start] << "\t" 
		  << RDtrans_status[mseg[j].end] << "\t" 
	  //<< mean(RDmedint,mseg[j].start,mseg[j].end) << "\t"
		  << testseglist[i].cnvmed << "\t" 
		  << testseglist[i].refmed << "\t" 
		  << testseglist[i].refsd << "\t" 
		  << testseglist[i].geno << "\t"
		  << testseglist[i].p1 << endl;
    }
    
    for(j=mseg.size()-1;j>=0;--j) {
      if ( mseg[j].geno==0 ) continue;
      if ( iseg.geno==0 ) iseg=mseg[j];
      if ( mseg[j].length>iseg.length ) iseg=mseg[j];
    }
    
    if ( iseg.geno==0 ) iseg=oseg;
    testseglist[i]=iseg;
    
    if ( testseglist[i].p1<rsi::p && rsi::debug ) {
      rsi::dout << "*" << "\t"
		<< testseglist[i].start << "\t" 
		<< testseglist[i].end << "\t"
		<< testseglist[i].length << "\t" 
		<< rsi::cnv_type[testseglist[i].type] << "\t"
		<< testseglist[i].geno << "\t"
		<< testseglist[i].cnvmed << "\t"
		<< testseglist[i].refmed << "\t"
		<< testseglist[i].p1 << "\t"
		<< "-seg-opt-" << endl;
    }
    if ( rsi::debug ) {
      rsi::dout << i+1 << "\t" 
		<< testseglist[i].start << "\t" 
		<< testseglist[i].end << "\t" 
		<< testseglist[i].length << "\t" 
                << rsi::cnv_type[testseglist[i].type] << "\t"
		<< testseglist[i].geno << "\t"
                << testseglist[i].cnvmed << "\t"
		<< testseglist[i].refmed << "\t"
		<< "-final-" << endl;
    }
  }

  rsiseglist=testseglist;
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void sortcnvstartposition(vector<cnv_st>& cnvlist)
{
  int i;
  
  Array<int> starts(cnvlist.size());
  for(i=0;i<(int)cnvlist.size();++i) {
    if ( cnvlist[i].start>cnvlist[i].end ) SWAP(cnvlist[i].start,cnvlist[i].end);
    starts[i]=cnvlist[i].start;
  }
  
  Array<int> starts_rk(cnvlist.size(),-1);
  arrayindex(starts,starts_rk,1); 
  if ( rsi::debug ) cout << "indexing done" << endl;
  vector<cnv_st> testlist;
  testlist.clear();
  for(i=0;i<(int)cnvlist.size();++i) testlist.push_back(cnvlist[starts_rk[i]]);
  cnvlist=testlist;
  for(i=0;i<(int)cnvlist.size()-1;++i) 
    if ( cnvlist[i].start > cnvlist[i+1].start ) {
      rsi::dout << "unsorted ? << " << endl
		<< i << "\t" 
		<< cnvlist[i].start << "\t" << cnvlist[i].end << endl
		<< i+1 << "\t" 
		<< cnvlist[i+1].start << "\t" << cnvlist[i+1].end
		<< endl; 
      exit(0);
    }

}


/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
string cnv_format1(cnv_st& icnv)
{
  double q1,q2;
  if ( icnv.p1<1.0E-10 ) q1=99;
  else q1=-10.0*log(icnv.p1)/log(10.0);
  if ( icnv.p2<1.0E-10 ) q2=99;
  else q2=-10.0*log(icnv.p2)/log(10.0);
  
  std::stringstream ss;
  
  if ( icnv.start<0 && icnv.end<0 ) {
    ss << "#CHROM" << "\t"
       << "START" << "\t"
       << "END" << "\t"
       << "TYPE" << "\t"
       << "SCORE" << "\t"
       << "LENGTH" << "\t"
       << "CNV_MED" << "("
       << "CNV_SD" << ");"
       << "NEIGHBOR_MED" << "("
       << "NEIGHBOR_RUNMEANSD" << ");"
       << "CHR_MED" << "("
       << "CHR_SD" << ")\t"
       << "RP=#support_read_pairs" << ";"
       << "Q0=#fraction_of_Q0_reads" << "\t"
       << "METHOD";
  }
  else {
    string RNAME= rsi::tid >= 0 && rsi::tid < (int) rsi::target_name.size() ? 
      rsi::target_name[rsi::tid] : "chr";
    ss << RNAME << "\t"
       << icnv.start << "\t"
       << icnv.end << "\t"
       << rsi::cnv_type[icnv.type] << "\t"
       << (int)q1 << "\t"
       << icnv.end-icnv.start+1 << "\t"
       << icnv.cnvmed << "("
       << icnv.cnviqr/1.349 << ");"
       << icnv.refmed << "("
       << icnv.refiqr/1.349 << ");"
       << rsi::RDmedian << "("
       << rsi::RDsd << ")\t"
      // << icnv.p1 << "\t"
       << "RP=" << icnv.rp << ";"
       << "Q0=" << icnv.q0 << "\t"
      // << icnv.score << "\t"
       << "rsi";
  }
  
  return(ss.str());
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
string cnv_format_all(cnv_st& icnv)
{
  double q1,q2;
  if ( icnv.p1<1.0E-10 ) q1=99;
  else q1=-10.0*log(icnv.p1)/log(10.0);
  if ( icnv.p2<1.0E-10 ) q2=99;
  else q2=-10.0*log(icnv.p2)/log(10.0);
  
  std::stringstream ss;
  
  if ( icnv.start<0 && icnv.end<0 ) {
    ss << "#CHROM" << "\t"
       << "START" << "\t"
       << "END" << "\t"
       << "TYPE" << "\t"
       << "LENGTH" << "\t"
       << "CNV_MEAN" << "("
       << "CNV_SD" << ");"
       << "NEIGHBOR_MEAN" << "("
       << "NEIGHBOR_RUNMEANSD" << ");"
       << "CHR_MEDIAN" << "("
       << "CHR_SD" << ")\t"
       << "LOSS/GAIN" << "\t"
       << "QUAL+-1" << "\t"
       << "QUAL+-2" << "\t"
       << "q0" << "\t"
       << "RSIscore" << "\t"
       << "METHOD";
  }
  else {
    string RNAME= rsi::tid >= 0 && rsi::tid < (int)rsi::target_name.size() ? 
      rsi::target_name[rsi::tid] : "chr";
    ss << RNAME << "\t"
       << icnv.start << "\t"
       << icnv.end << "\t"
       << rsi::cnv_type[icnv.type] << "\t"
       << icnv.end-icnv.start+1 << "\t"
       << icnv.cnvmed << "("
      //<< icnv.cnvsd << ","
       << icnv.cnviqr/1.349 << ");"
       << icnv.refmed << "("
      //<< icnv.refsd << ","
       << icnv.refiqr/1.349 << ");"
       << rsi::RDmedian << "("
       << rsi::RDsd << ")\t"
      //<< icnv.sc1 << ";"
      //<< icnv.sc2 << "\t"
      //<< icnv.geno << "\t"
      //<< (int)q1 << "\t"
      //<< (int)q2 << "\t"
       << "RP=" << icnv.rp << ";"
       << "Q0=" << icnv.q0 << "\t"
       << icnv.score << "\t"
       << "rsi";
  }
  
  return(ss.str());
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void mergesegments(Array<int>& RD, vector<cnv_st>& cnvlist)
{
  int i;
  cnv_st icnv;
  vector<cnv_st> testcnvlist;
  testcnvlist.clear();
  
  if (rsi::debug) 
    rsi::dout << "merging " << cnvlist.size() << " segments" << endl;
  
  // first test for overlap
  //
  for(i=0;i<(int)cnvlist.size()-1;++i) {
    if ( cnvlist[i].type != cnvlist[i+1].type ) continue;
    if ( max(cnvlist[i].start, cnvlist[i+1].start) <
	 min(cnvlist[i].end, cnvlist[i+1].end) ) {
      
      if (rsi::debug)
	rsi::dout << "cnv overlap detected" << endl
		  << cnvlist[i].start << "\t" 
		  << cnvlist[i].end << "\t" 
		  << cnvlist[i].cnvmed << "\t"
		  << cnvlist[i].refmed << "\t"
		  << cnvlist[i].p1 << "\t"
		  << cnvlist[i].geno << endl
		  << cnvlist[i+1].start << "\t" 
		  << cnvlist[i+1].end << "\t" 
		  << cnvlist[i+1].cnvmed << "\t"
		  << cnvlist[i+1].refmed << "\t"
		  << cnvlist[i+1].p1 << "\t"
		  << cnvlist[i+1].geno << endl;
      
      icnv=cnvlist[i];
      icnv.start=cnvlist[i+1].start>cnvlist[i].start ? 
	cnvlist[i].start:cnvlist[i+1].start;
      icnv.end=cnvlist[i+1].end>cnvlist[i].end ? 
	cnvlist[i+1].end:cnvlist[i].end;
      
      icnv.start=min(cnvlist[i].start, cnvlist[i+1].start);
      icnv.end=max(cnvlist[i].end, cnvlist[i+1].end);
      
      testcnvlist.clear();
      testcnvlist=cnvlist;
      testcnvlist[i]=icnv;
      testcnvlist[i+1]=icnv;
      testcnvlist[i+1].status=-9;
      isitcnvwrap(RD, testcnvlist, i);
      
      if ( rsi::debug) 
	rsi::dout << "merged cnv test result" << endl
		  << testcnvlist[i].start << "\t" 
		  << testcnvlist[i].end << "\t" 
		  << testcnvlist[i].cnvmed << "\t"
		  << testcnvlist[i].refmed << "\t"
		  << testcnvlist[i].p1 << "\t"
		  << testcnvlist[i].geno << endl;
      
      if ( testcnvlist[i].geno == 0 ) {
	
	testcnvlist=cnvlist;
	testcnvlist[i+1].status=-9;
	isitcnvwrap(RD, testcnvlist, i);
	
	testcnvlist[i].status=-9;
	testcnvlist[i+1].status=0;
	isitcnvwrap(RD, testcnvlist, i+1);
	
	if ( testcnvlist[i+1].p1 < testcnvlist[i].p1 ) 
	  testcnvlist[i]=testcnvlist[i+1];
	
	if ( testcnvlist[i].p1 > rsi::p ) {
	  cnvlist[i].status=-9;
	  cnvlist[i+1].status=-9;
	  if ( rsi::debug) 
	    rsi::dout << "neither segment is a good cnv, both are removed" << endl;
	      
	}
	
	if ( rsi::debug) 
	  rsi::dout << "merged cnv test result" << endl
		    << testcnvlist[i].start << "\t" 
		    << testcnvlist[i].end << "\t" 
		    << testcnvlist[i].cnvmed << "\t"
		    << testcnvlist[i].refmed << "\t"
		    << testcnvlist[i].p1 << "\t"
		    << testcnvlist[i].geno << endl;
      }
	
      if ( testcnvlist[i].geno == 0 ) {
	if ( rsi::debug) 
	  rsi::dout 
	    << "merged segment is not a cnv, segments are not merged" << endl;
	continue;
      }
      
      if (rsi::debug) rsi::dout << "overlapped cnvs are merged" << endl;
      cnvlist[i]=testcnvlist[i];
      cnvlist[i].status=-9;
      cnvlist[i+1]=testcnvlist[i];
      cnvlist[i+1].status=0;
      
    }
  }
  
  vector<cnv_st> cnvlist1;
  cnvlist1.clear();
  for(i=0;i<(int)cnvlist.size();++i) 
    if ( cnvlist[i].status != -9 ) cnvlist1.push_back(cnvlist[i]);
  cnvlist=cnvlist1; // tested for overlap
  
  //  return;
  if ( !rsi::merge ) return;
  // test for connecting two separated segments
  //  
  for(i=0;i<(int)cnvlist.size()-1;++i) {
    if ( cnvlist[i].type!=cnvlist[i+1].type ) continue;
    if ( cnvlist[i].geno==0 || cnvlist[i+1].geno==0 ) continue; 
    int gap=cnvlist[i+1].start-cnvlist[i].end;
    if ( gap > (cnvlist[i].end-cnvlist[i].start)*rsi::chklen*0.7 &&
         gap > (cnvlist[i+1].end-cnvlist[i+1].start)*rsi::chklen*0.7 ) continue;
    
    double gap_mean,merged_mean,cnv1_mean,cnv2_mean,cnv_mean;
    cnv1_mean = mean(RD, cnvlist[i].start, cnvlist[i].end);
    cnv2_mean = mean(RD, cnvlist[i+1].start, cnvlist[i+1].end);
    cnv_mean=(cnv1_mean*(cnvlist[i].end-cnvlist[i].start)+
	      cnv2_mean*(cnvlist[i+1].end-cnvlist[i+1].start))/
      ((cnvlist[i].end-cnvlist[i].start)+(cnvlist[i+1].end-cnvlist[i+1].start));
    gap_mean = mean(RD, cnvlist[i].end, cnvlist[i+1].start);
    merged_mean = mean(RD, cnvlist[i].start, cnvlist[i+1].end);
    
    if ( cnvlist[i].type==TYPE_DEL && 
	 merged_mean > cnv_mean+1.5*cnvlist[i+1].refsd+1.5*cnvlist[i].refsd )
      continue;
    if ( cnvlist[i].type==TYPE_DUP && 
	 merged_mean < cnv_mean-1.5*cnvlist[i+1].refsd-1.5*cnvlist[i].refsd ) 
      continue;
    
    if ( rsi::debug) 
      rsi::dout << "considering merging two separated segments" << endl
		<< cnvlist[i].start << "\t" 
		<< cnvlist[i].end << "\t" 
		<< cnvlist[i].cnvmed << "\t"
		<< cnvlist[i].refmed << "\t"
		<< cnvlist[i].geno << endl
		<< cnvlist[i+1].start << "\t" 
		<< cnvlist[i+1].end << "\t" 
		<< cnvlist[i+1].cnvmed << "\t"
		<< cnvlist[i+1].refmed << "\t"
		<< cnvlist[i+1].geno << "\t"
		<< merged_mean << endl;
    
    icnv=cnvlist[i];
    icnv.start=cnvlist[i].start;
    icnv.end=cnvlist[i+1].end;
    testcnvlist.clear();
    testcnvlist=cnvlist;
    testcnvlist[i]=icnv;
    testcnvlist[i+1]=icnv;
    testcnvlist[i+1].status=-9;
    isitcnvwrap(RD, testcnvlist, i);
    
    if (rsi::debug) rsi::dout << "merged cnv test result" << endl
			      << testcnvlist[i].start << "\t" 
			      << testcnvlist[i].end << "\t" 
			      << testcnvlist[i].cnvmed << "\t"
			      << testcnvlist[i].refmed << "\t"
			      << testcnvlist[i].p1 << "\t"
			      << testcnvlist[i].geno << endl;
    
    if ( testcnvlist[i].geno == 0 ) {
      if (rsi::debug)
	rsi::dout << "merged segment is not a cnv, segments not merged" << endl;
      continue;
    }
    if (rsi::debug) rsi::dout << "merged" << endl;
    cnvlist[i]=testcnvlist[i];
    cnvlist[i+1]=testcnvlist[i];
    cnvlist[i].status=-9;
    
  }
  
  cnvlist1.clear();
  for(i=0;i<(int)cnvlist.size();++i) 
    if ( cnvlist[i].status != -9 ) cnvlist1.push_back(cnvlist[i]);
  cnvlist=cnvlist1;  // return vector
  
  if (rsi::debug) 
    rsi::dout << "Done merging " 
	      << "got " << cnvlist.size() << " segments" 
	      << endl << endl;
  
} //void mergesegments


/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void optimize_with_derivative(Array<int>& RD, cnv_st& icnv)
{
  int len=icnv.end-icnv.start+1;
  //  int displacement=max(500,len/2);
  int displacement=max(250,len/4);
  
  int nstart=icnv.start-displacement;
  int nend=icnv.end+displacement;
  
  // too close to ends
  if ( nstart<2*len ) return;
  if ( nend>RD.size()-2*len ) return;
  
  vector<double> dd(0);
  double diff=0.0;
  for(int k=nstart-len; k<nstart; ++k) diff+=RD[k];
  for(int k=nstart; k<nstart+len; ++k) diff-=RD[k];
  dd.push_back(diff);
  for(int i=nstart+1; i<nend; ++i ) {
    diff=diff-RD[i-1-len]+RD[i-1]+RD[i-1]-RD[i-1+len];
    dd.push_back(diff);
  }
  
  int imax=-1;
  double ddmax=0;
  for(int i=0; i<2*displacement; ++i) {
    if ( icnv.type==TYPE_DEL ) {
      if ( dd[i] > ddmax ) { ddmax=dd[i]; imax=i; }
    }
    if ( icnv.type==TYPE_DUP ) {
      if ( dd[i] < ddmax ) { ddmax=dd[i]; imax=i; }
    }
  }
  if ( imax>0 ) icnv.start=nstart+imax;
  
  imax=-1;
  ddmax=0;
  for(int i=dd.size()-2*displacement; i<(int)dd.size(); ++i) {
    if ( icnv.type==TYPE_DEL ) {
      if ( dd[i] < ddmax ) { ddmax=dd[i]; imax=i; }
    }
    if ( icnv.type==TYPE_DUP ) {
      if ( dd[i] > ddmax ) { ddmax=dd[i]; imax=i; }
    }
  }
  if ( imax>0 ) icnv.end=nend-dd.size()+imax;

  return;
}

void optimize_with_derivative(Array<int>& RD, vector<cnv_st>& cnvlist)
{
  for(size_t i=0; i<cnvlist.size(); ++i) 
    optimize_with_derivative(RD, cnvlist[i]);
  return;
}


/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
template<class Etype>
void filterstatus_tp(Array<Etype>& RDtrans, double dev, Array<int>& RDtrans_status) 
{
  cnv_st icnv;
  
  int i;
  
  int minlevel=min(RDtrans_status, 0, RDtrans_status.size()-1);
  int maxlevel=max(RDtrans_status, 0, RDtrans_status.size()-1);
  
  int mlevels=maxlevel-minlevel+1;
  Array<int> level(mlevels,0);
  
  for(i=0;i<RDtrans_status.size();++i) level[ RDtrans_status[i]-minlevel ]=1;
  
  vector<int> LL;
  LL.clear();
  for(i=0;i<mlevels;++i) if ( level[i] == 1 ) LL.push_back(i);
  
  Array<float> RDsum(mlevels,0.0);
  Array<int> RDsumcount(mlevels,0);
  int ilevel;
  for(i=0;i<RDtrans.size();++i) {
    ilevel = RDtrans_status[i]-minlevel;
    RDsum[ ilevel ] += RDtrans[i];
    ++RDsumcount[ ilevel ];
  }
  for(ilevel=0;ilevel<RDsum.size();++ilevel) 
    if ( RDsumcount[ilevel] !=0 ) RDsum[ilevel]/=(double)RDsumcount[ilevel];
  
  int leveldel=minlevel;
  for(ilevel=0;ilevel<RDsum.size();++ilevel) 
    if ( RDsum[ilevel] < RDsum[-minlevel]-dev ) {
      leveldel=ilevel+minlevel;
      break;
    }
  int leveladd=maxlevel;
  for(ilevel=RDsum.size()-1;ilevel>=0;ilevel--) 
    if ( RDsum[ilevel] > RDsum[-minlevel]+dev ) {
      leveladd=ilevel+minlevel;
      break;
    }
  
  for(ilevel=0;ilevel<RDsum.size();++ilevel) 
    if ( RDsumcount[ilevel] !=0 ) 
      rsi::dout << ilevel + minlevel << "\t"
		<< RDsumcount[ilevel] << "\t"
		<< RDsum[ilevel] << endl;
  rsi::dout << leveldel << "\t" << RDsum[leveldel-minlevel] << endl
	    << leveladd << "\t" << RDsum[leveladd-minlevel] << endl;
  
  if ( leveldel>0 || leveladd<0 || leveldel>leveladd ) {
    rsi::dout << "warning level error, status not filtered" << endl
	      << leveldel << "\t" << leveladd << endl
	      << "filterstatus()" << endl;
    return;
  } 
  
  // remove bins close to median
  //
  double delthreshold=RDsum[-minlevel]-dev;
  double addthreshold=RDsum[-minlevel]+dev;
  
  /*
  int beforecount=0,aftercount=0;
  for(i=0;i<RDtrans_status.size();++i) {
    if ( RDtrans_status[i] != 0 ) ++beforecount;
    // if ( RDtrans_status[i] < leveldel || RDtrans_status[i] > leveladd )
    //  RDtrans_status[i]=0;
    
    if ( RDtrans_status[i]<=leveldel && RDtrans[i]>delthreshold ) RDtrans_status[i]=0;
    if ( RDtrans_status[i]>=leveladd && RDtrans[i]<addthreshold ) RDtrans_status[i]=0;
      
    if ( RDtrans_status[i] != 0 ) ++aftercount;
  }
  rsi::dout << beforecount << "\t-->\t" << aftercount << endl;
  */
  
  // trim edge
  vector<cnv_st> seglist;
  get_continuous_segments(RDtrans_status, 1, seglist);
  for(i=0; i<(int)seglist.size(); ++i) {
    int i1=seglist[i].start;
    int i2=seglist[i].end;
    while ( ( RDtrans[i1]>delthreshold && RDtrans_status[i1]<0 ) ||
	    ( RDtrans[i1]<addthreshold && RDtrans_status[i1]>0 ) ) {
      RDtrans_status[i1] = 0;
      ++i1;
      if ( i1>=i2 ) break;
    }
    while ( ( RDtrans[i2]>delthreshold && RDtrans_status[i2]<0 ) ||
	    ( RDtrans[i2]<addthreshold && RDtrans_status[i2]>0 ) ) {
      RDtrans_status[i2] = 0;
      --i2;
      if ( i2<=i1 ) break;
    }
  }
  
  return;
}
void filterstatus(Array<int>& RDtrans, double dev, Array<int>& RDtrans_status) 
{
  filterstatus_tp(RDtrans, dev, RDtrans_status);
  return;
}
void filterstatus(Array<float>& RDtrans, double dev, Array<int>& RDtrans_status) 
{
  filterstatus_tp(RDtrans, dev, RDtrans_status);
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void get_rsi_segments(Array<float>& RDmed, 
		      Array<int>& RDstatus, 
		      double tmedian,
		      vector <cnv_st>& rsiseglist)
{
  cnv_st icnv;
  vector<cnv_st> seglist;
  seglist.clear();
  get_continuous_segments(RDstatus, 1, seglist);
  
  int i, j, L, istart, iend, nstart, nend;
  int optistart,optlength;
  
  for(i=0;i<(int)seglist.size();++i) {
    istart=seglist[i].start;
    iend=seglist[i].end;
    
    nstart=istart;
    nend=iend;
    optistart=istart;
    optlength=iend-istart+1;
    double score=0,optscore=0, sum=0;
    for(L=1;L<=iend-istart+1;++L) {
      sum=0.0;
      for(j=istart;j<=iend && j<istart+L;++j) sum+=RDmed[j];
      score= abs( sum/(double)L - tmedian ) * sqrt(double(L));
      if ( score > optscore ) { 
	nstart=istart; 
	nend=nstart+L-1; 
	optscore=score; 
      }
      for(j=istart+1;j+L-1<=iend;++j) {
	sum=sum-RDmed[j-1]+RDmed[j+L-1];
	score= abs( sum/(double)L - tmedian ) * sqrt(double(L));
	if ( score > optscore ) { 
	  nstart=j; 
	  nend=nstart+L-1; 
	  optscore=score; 
	}
      }
    }
    
    icnv.start=nstart;
    icnv.end=nend;
    // if ( mean(RDstatus,icnv.start,icnv.end) > 0 ) {
    if ( _median(&RDstatus[icnv.start],icnv.end-icnv.start+1) > 0 ) {
      icnv.type = TYPE_DUP;
      icnv.score = optscore;
    }
    else {
      icnv.type = TYPE_DEL;
      icnv.score = -optscore;
    }
    
    rsiseglist.push_back(icnv);
  }
  
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void negative_binomial_transfer(Array<int>& RD, int m, Array<float>& RDt)
{
  int i,j,k;
  int i1,i2;
  double RDmedian,RDmad,r;
  double del_nbt,med_nbt,dup_nbt;
  
  RDmedian=_median(&RD[0],RD.size());  
  
  int ns=31;
  Array<int> RDtmp(RD.size()/ns);
  Array<double> RDMAD(ns,0.0);
  for(j=0;j<ns;j++) {
    for(i=0+j,k=0; i<RD.size() && k<RDtmp.size(); i+=ns,++k) 
      RDtmp[k]=abs((float)RD[i]-RDmedian);
    RDmad=_median(&RDtmp[0],RDtmp.size());  
    RDMAD[j]=RDmad;
  }
  RDmad=_median(&RDMAD[0], RDMAD.size());
  
  rsi::dout << "RD median : " << RDmedian << endl;  
  rsi::dout << "RD median absolute deviation : " << RDmad << endl;  
  
  r=RDmedian/RDmad;
  
  RDt.resize(RD.size()/m);
  double sum=0.0;
  for(i=0;i<RDt.size();++i) {
    i1=i*m;
    i2=i*m+m-1;
    if ( i2 > (RD.size()-1) ) i2=RD.size()-1;
    double m2;
    m2=double(i2-i1+1);
    sum=0.0;
    for(j=i1;j<=i2;++j) sum+=RD[j];
    RDt[i]=2.0*sqrt(r)*log( sqrt( (sum+0.25)/(m2*r-0.5) ) 
			    + sqrt( 1.0+(sum+0.25)/(m2*r-0.5) ) );
  }
  
  sum=RDmedian*m;
  med_nbt=2.0*sqrt(r)*log( sqrt( (sum+0.25)/((double)m*r-0.5) ) 
			   + sqrt( 1.0+(sum+0.25)/((double)m*r-0.5) ) );
  sum=RDmedian/2.0*(double)m;
  del_nbt=2.0*sqrt(r)*log( sqrt( (sum+0.25)/((double)m*r-0.5) ) 
			   + sqrt( 1.0+(sum+0.25)/((double)m*r-0.5) ) );
  sum=RDmedian*1.5*(double)m;
  dup_nbt=2.0*sqrt(r)*log( sqrt( (sum+0.25)/((double)m*r-0.5) ) 
			   + sqrt( 1.0+(sum+0.25)/((double)m*r-0.5) ) );
  
  /* scale back to original RD scale */
  // double tmin=min(RDt,0,RDt.size()-1);
  double tmin=RDt[0];
  for(i=0;i<RDt.size();++i) if ( RDt[i]<tmin ) tmin=RDt[i];
  for(i=0;i<RDt.size();++i) RDt[i]-=tmin;
  med_nbt-=tmin;
  for(i=0;i<RDt.size();++i) RDt[i]=RDt[i]/med_nbt*RDmedian;
  
  del_nbt-=tmin;
  dup_nbt-=tmin;
  del_nbt=del_nbt/med_nbt*RDmedian;
  dup_nbt=dup_nbt/med_nbt*RDmedian;
  med_nbt=med_nbt/med_nbt*RDmedian;
  
  RDt[0]=del_nbt;
  RDt[1]=dup_nbt;
  RDt[2]=med_nbt;  
  /* --- */
  
} 

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void rsistatus(Array<float>& RDtrans, Array<int>& RDmedint, 
	       double tmedian, double tlamda, int Lmax,
	       Array<int>& RDtrans_status )
{
  int i,j;
  Array<float> RDtransmean(RDtrans.size(),0.0);
  RDtrans_status.assign(0);
  
  /* first mark deletions */
  for(int L=1;L<=Lmax;L+=1) { 
    RDtransmean.assign(0.0);
    runmean(RDtrans, RDtransmean, RDtrans.size(), L, 1);
    
    for(i=L/2+1;i<RDtransmean.size()-L/2-1;++i) {
      double score=(RDtransmean[i]-tmedian)*sqrt(double(L));
      if ( score > -tlamda ) continue;
      int i1=i-L/2;
      int i2=i1+L-1;
      if ( alglib::median(&RDmedint[i1],L) >  rsi::RDmedian*0.75 ) continue;
      // if ( alglib::percentile(&RDmedint[i1],L,0.50) >  rsi::RDmedian*0.75 ) continue;
      while ( RDtrans[i1] > tmedian ) i1++;
      while ( RDmedint[i1] > rsi::RDmedian*0.75 ) i1++;
      while ( RDtrans[i2] > tmedian ) i2--;
      while ( RDmedint[i2] > rsi::RDmedian*0.75 ) i2--;
      for(j=i1;j<=i2;++j) if ( RDtrans_status[j]==0 ) RDtrans_status[j]=-L;
    }
    
    int icount=0;
    for(i=0;i<RDtrans.size();++i) icount+=(int)(RDtrans_status[i]<0); 
    double portion=double(icount)/double(RDtrans.size());
    rsi::dout << "DEL-\t" << L << "\t" 
	      << icount << "\t" 
	      << RDtrans.size() << "\t" 
	      << portion << endl;
    if ( portion > 0.2 ) break;
  }
  
  /* second mark duplications */
  for(int L=1;L<=Lmax;L+=1) { 
    RDtransmean.assign(0.0);
    runmean(RDtrans, RDtransmean, RDtrans.size(), L, 1);
    
    for(i=L/2+1;i<RDtransmean.size()-L/2-1;++i) {
      double score=(RDtransmean[i]-tmedian)*sqrt(double(L));
      if ( score < tlamda ) continue;
      int i1=i-L/2;
      int i2=i1+L-1;
      if ( alglib::median(&RDmedint[i1],L) <  rsi::RDmedian*1.25 ) continue;
      //if ( alglib::percentile(&RDmedint[i1],L,0.50) <  rsi::RDmedian*1.3 ) continue;

      while ( RDtrans[i1] < tmedian ) i1++;
      while ( RDmedint[i1] < rsi::RDmedian*1.25 ) i1++;
      while ( RDtrans[i2] < tmedian ) i2--;
      while ( RDmedint[i2] < rsi::RDmedian*1.25 ) i2--;
      for(j=i1;j<=i2;++j) if ( RDtrans_status[j]==0 ) RDtrans_status[j]=L;
    }
    
    int icount=0;
    for(i=0;i<RDtrans.size();++i) icount+=(int)(RDtrans_status[i]>0); 
    double portion=double(icount)/double(RDtrans.size());
    rsi::dout << "DUP+\t" << L << "\t" 
	      << icount << "\t" 
	      << RDtrans.size() << "\t" 
	      << portion << endl;
    if ( portion > 0.2 ) break;
  }

  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void rsicnvnbn(Array<float>& RDtrans, Array<int>& RDmedint, 
	       Array<int>& RDtrans_status, 
	       vector<cnv_st>& cnvlist )
{
  
  cnvlist.clear();
  RDtrans_status.assign(0);
  
  int i,k;
  double tmedian,tlamda,tsigma;
  
  tmedian = _median(&RDtrans[0], RDtrans.size()) ;
  rsi::nbnmedian=tmedian;
  Array<float> RDtmp((RDtrans.size()));
  for(i=0;i<RDtmp.size();++i) RDtmp[i]= abs(RDtrans[i]-tmedian);
  tsigma=_median(&RDtmp[0], RDtmp.size())/0.6745;
  tlamda=rsi::factor*tsigma;
  // single bin is not accepted
  double target_tlamda = (RDtrans[2]-RDtrans[0])*sqrt(2.5);
  tlamda=max(tlamda,target_tlamda);
  
  rsi::nbnlamda=tlamda;
  
  /* set max length */
  int Lmax=rsi::Lmax;
  double dnb=abs(RDtrans[2]-RDtrans[0])+0.0001;
  int cal_max=pow(tlamda*2/dnb, 2);
  if ( Lmax < cal_max )  Lmax=cal_max;
  
  rsi::dout << "Negative binormial transformation : " << endl
	    << "\tmedian of transformations : " <<  tmedian << endl
	    << "\tsigma : " << tsigma << endl
	    << "\tmedian/sigma : " << tmedian/tsigma << endl
	    << "\trsifactor: " << rsi::factor << endl
	    << "\tlamda : " << tlamda << endl
	    << "\ttarget_tlamda : " << target_tlamda << endl
	    << "\tMax L needed  : " << cal_max << endl;
  
  
  /* calculate rsi status */
  rsistatus(RDtrans, RDmedint, tmedian, tlamda, Lmax, RDtrans_status);
  /* filter out points too close to mean, +- dev */
  double dev=tsigma*3.0;
  filterstatus(RDtrans, dev, RDtrans_status);
  
  k=0;
  for(i=0;i<RDtrans.size();++i) 
    if ( RDtrans_status[i]==0 ) {
      RDtmp[k]=RDtrans[i];
      k++;
    }
  if ( k>RDtrans.size()/2 ) {
    tmedian=_median(&RDtmp[0], k);
    for(i=0;i<k;++i) RDtmp[i]= abs(RDtmp[i]-tmedian);
    tsigma=_median(&RDtmp[0],k)/0.6745;
    tlamda=rsi::factor*tsigma;
    tlamda=max(tlamda, target_tlamda);
  }
  rsi::dout << "second pass" << endl
	    << "Negative binormial transformation : " << endl
	    << "\tmedian of transformations : " <<  tmedian << endl
	    << "\tsigma : " << tsigma << endl
	    << "\tmedian/sigma : " << tmedian/tsigma << endl
	    << "\tlamda : " << tlamda << endl
	    << "\target_tlamda : " << target_tlamda << endl;
  rsi::nbnlamda=tlamda;
  rsi::nbnmedian=tmedian;
  rsistatus(RDtrans, RDmedint, tmedian, tlamda, Lmax, RDtrans_status);
  // filterstatus(RDtrans, tsigma, RDtrans_status);
  // filterstatus(RDmedint, rsi::RDmedian*0.25, RDtrans_status);
  
  vector<cnv_st> seglist;
  get_continuous_segments(RDtrans_status, 1, seglist);   
  
  /* generate continuous segments with max score */ 
  vector<cnv_st> rsiseglist;
  rsiseglist.clear();
  get_rsi_segments(RDtrans, RDtrans_status, tmedian, rsiseglist);
  
  /* make sure selected segments comply with rsi algorithm */ 
  cnvlist.clear();
  for(i=0;i<(int)rsiseglist.size();++i) {
    if ( abs(rsiseglist[i].score) < tlamda*0.5 ) rsiseglist[i].status=-9;
    if ( rsiseglist[i].status!=-9 ) cnvlist.push_back(rsiseglist[i]);
  }
  if (rsi::debug){
    for(i=0;i<(int)cnvlist.size();++i) {
      rsi::dout << i << "\t" 
		<< cnvlist[i].start << "\t" 
		<< cnvlist[i].end << "\t" 
		<< cnvlist[i].end-cnvlist[i].start+1 << "\t"
		<< mean(RDtrans_status,cnvlist[i].start,cnvlist[i].end) << "\t"
		<< cnvlist[i].score << "\t"
		<< cnvlist[i].status << endl;
    }
  }
  
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void median_transfer(Array<int>& RD, int m, Array<float>& RDt)
{
  
  if ( m<1 ) { 
    rsi::dout << "m<1"  << endl
	      << "median_transfer()" << endl;
    exit(0);
  }
  RDt.resize(RD.size()/m);
  for(int i=0;i<RDt.size();++i) {
    Array<int> RDtmp(m);
    for(int j=i*m;j<i*m+m;++j) RDtmp[j-i*m]=RD[j];
    // RDt[i]=_median(&RDtmp[0],m);
    RDt[i]=alglib::median(&RDtmp[0],m);
  }
  
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void median_transfer(Array<int>& RD, int m, Array<double>& RDt)
{
  
  if ( m<1 ) { 
    rsi::dout << "m<1"  << endl
	      << "median_transfer()" << endl;
    exit(0);
  }
  RDt.resize(RD.size()/m);
  for(int i=0;i<RDt.size();++i) {
    Array<int> RDtmp(m);
    for(int j=i*m;j<i*m+m;++j) RDtmp[j-i*m]=RD[j];
    // RDt[i]=_median(&RDtmp[0],m);
    RDt[i]=alglib::median(&RDtmp[0],m);
  }
  

}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void rsicnvmed(Array<float>& RDtrans, 
	       Array<int>& RDmedint, 
	       Array<int>& RDtrans_status, 
	       vector<cnv_st>& cnvlist )
{
  
  cnvlist.clear();
  int i,k;
  double tmedian,tlamda,tsigma;
  
  //tmedian = samplemedian(RDtrans, 0, RDtrans.size()-1) ;
  tmedian = rsi::RDmedian;
  rsi::median=tmedian;
  rsi::medmedian=tmedian;
  Array<float> RDtmp((RDtrans.size()));
  for(i=0;i<RDtmp.size();++i) RDtmp[i]= abs(RDtrans[i]-tmedian);
  tsigma=_median(&RDtmp[0],RDtmp.size())/0.6745;
  //tsigma=alglib::median(&RDtmp[0],RDtmp.size())/0.6745;
  tlamda=rsi::factor*tsigma;
  double target_tlamda = max( 0.5*tmedian*sqrt(800.0/(double)rsi::m),
			      tmedian*sqrt(2.0) );
  
  // two zero bins not accepted
  target_tlamda = tmedian*sqrt(2.0);
  tlamda=max(tlamda,target_tlamda);
  if ( rsi::threshold > 0 ) tlamda=tmedian*rsi::threshold;
  rsi::medlamda=tlamda;
  
  /* set max length */
  int Lmax=rsi::Lmax;
  int cal_max=pow(tlamda*4/(tmedian+0.001),2);
  if ( Lmax < cal_max )  Lmax=cal_max;
  
  
  rsi::dout << "local median transformation : " << endl
	    << "\tmedian of transformations : " <<  tmedian << endl
	    << "\tsigma : " << tsigma << endl
	    << "\tmedian/sigma : " << tmedian/tsigma << endl
	    << "\trsifactor: " << rsi::factor << endl
	    << "\tlamda : " << tlamda << endl
	    << "\ttarget_tlamda : " << target_tlamda << endl
	    << "\tMax L needed  : " << cal_max << endl;
  
  if ( tlamda < tmedian*sqrt(2.5) ) 
    rsi::dout << "Warning : lamda might be low" << endl;
  if ( tlamda > tmedian*sqrt(3.0) ) 
    rsi::dout << "Warning : lamda might be too large" << endl;
  
  /* calculate rsi status */
  /* filter out points too close to mean, +- 0.2 mean */
  rsistatus(RDtrans, RDmedint, tmedian, tlamda, Lmax, RDtrans_status);
  double dev=tsigma;
  dev = tmedian*0.6;
  filterstatus(RDtrans, dev, RDtrans_status);
  
  k=0;
  for(i=0;i<RDtrans.size();++i) 
    if ( RDtrans_status[i]==0 ) {
      RDtmp[k]=RDtrans[i];
      k++;
    }
  if ( k>RDtrans.size()/2 ) {
    tmedian=_median(&RDtmp[0],k);
    //tmedian=alglib::median(&RDtmp[0],k);
    for(i=0;i<k;++i) RDtmp[i]= abs(RDtmp[i]-tmedian);
    tsigma=_median(&RDtmp[0],k)/0.6745;
    //tsigma=alglib::median(&RDtmp[0],k)/0.6745;
    tlamda=rsi::factor*tsigma;
    tlamda=max(tlamda,target_tlamda);
  }
  rsi::dout << "second pass" << endl
	    << "local median transformation : " << endl
	    << "\tmedian of transformations : " <<  tmedian << endl
	    << "\tsigma : " << tsigma << endl
	    << "\tmedian/sigma : " << tmedian/tsigma << endl
	    << "\tlamda : " << tlamda << endl
	    << "\target_tlamda : " << target_tlamda << endl;
  rsi::medlamda=tlamda;
  rsi::median=tmedian;
  rsi::medmedian=tmedian;
  rsistatus(RDtrans, RDmedint, tmedian, tlamda, Lmax, RDtrans_status);
  // filterstatus(RDtrans, dev, RDtrans_status);
  // filterstatus(RDmedint, rsi::RDmedian*0.25, RDtrans_status);
  
  vector<cnv_st> seglist;
  get_continuous_segments(RDtrans_status, 1, seglist);   
  
  if ( rsi::debug ) rsi::dout << procpidstatus( rsi::pidfile, "vm");
  
  /* generate continuous segments with max score */ 
  vector<cnv_st> rsiseglist;
  rsiseglist.clear();
  get_rsi_segments(RDtrans, RDtrans_status, tmedian, rsiseglist);

  /* make sure selected segments comply with rsi algorithm */ 
  cnvlist.clear();
  for(i=0;i<(int)rsiseglist.size();++i) {
    if ( abs(rsiseglist[i].score) < tlamda*0.5 ) rsiseglist[i].status=-9;
    if ( rsiseglist[i].status!=-9 ) cnvlist.push_back(rsiseglist[i]);
  }
  if (rsi::debug){
    for(i=0;i<(int)cnvlist.size();++i) {
      rsi::dout << i << "\t" 
		<< cnvlist[i].start << "\t" 
		<< cnvlist[i].end << "\t" 
		<< cnvlist[i].end-cnvlist[i].start+1 << "\t"
		<< mean(RDtrans_status,cnvlist[i].start,cnvlist[i].end) << "\t"
		<< cnvlist[i].score << "\t"
		<< cnvlist[i].status << endl;
    }
  }

  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*! 
  @abstract convert index of the concatenated RD array to original
  @param  p1    index for the concatenated array
  @return p2    index on the reference
  @discussion   both p1 and p2 are 0-based; needs rsi::noncodelist
*/
int expand_coordinate(int p1)
{
  if ( rsi::noncodelist.size()==0 ) return( p1 );
  
  vector<int> breaks(0);
  vector<int> incrments(0);
  int dx=0;
  for(int i=0;i<(int)rsi::noncodelist.size();++i) {
    if ( rsi::noncodelist[i].end<=rsi::noncodelist[i].start ) {
      rsi::dout << "N region error " << rsi::noncodelist[i].start 
	   << "\t" << rsi::noncodelist[i].end << endl;
      exit(0);
    }
    dx+=rsi::noncodelist[i].end-rsi::noncodelist[i].start+1;
    breaks.push_back(rsi::noncodelist[i].end+1-dx);
    incrments.push_back(dx);
  }
  
  if ( p1<breaks[0] ) return( p1 );
  if ( p1>=breaks.back() ) return( p1+incrments.back() );
  
  for(int i=0; i<(int)breaks.size()-1;++i) {
    if ( p1>=breaks[i] && p1<breaks[i+1] ) return( p1+incrments[i] );
  }
  
  rsi::dout << "[xpand_coordinate] do't know what to return" << endl;
  exit(0);
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*! 
  @abstract convert pos coordinate to the index of concatenated array
  @param  p1    index for the concatenated array
  @return p2    index on the reference
  @discussion   both p1 and p2 are 0-based; needs rsi::noncodelist
*/
int shrink_coordinate(int p1)
{
  return(p1);
  if ( rsi::noncodelist.size()==0 ) return( p1 );
  
  vector<int> breaks(0);
  vector<int> incrments(0);
  int dx=0;
  for(int i=0;i<(int)rsi::noncodelist.size();++i) {
    if ( rsi::noncodelist[i].end<=rsi::noncodelist[i].start ) {
      rsi::dout << "N region error " << rsi::noncodelist[i].start 
		<< "\t" << rsi::noncodelist[i].end << endl;
      exit(0);
    }
    dx+=rsi::noncodelist[i].end-rsi::noncodelist[i].start+1;
    breaks.push_back(rsi::noncodelist[i].end+1);
    incrments.push_back(dx);
    if ( p1+1 >= rsi::noncodelist[i].start && p1+1 <= rsi::noncodelist[i].end )
      return(-1);
  }
  
  if ( p1<breaks[0] ) return( p1 );
  if ( p1>=breaks.back() ) return( p1-incrments.back() );
  
  for(int i=breaks.size()-1; i>=0; --i)
    if ( p1>=breaks[i] ) return( p1-incrments[i] );
  
  rsi::dout << "[xpand_coordinate] do't know what to return" << endl;
  exit(0);
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void write_cnv_to_file(vector<cnv_st>& cnvlist, string outfile) 
{
  static bool printed=false;
  ofstream FOUT;
  if ( !printed ) FOUT.open(outfile.c_str());
  else FOUT.open(outfile.c_str(), std::ofstream::app);
  
  if ( ! printed ) {
    if ( rsi::rdfile != "" )
      FOUT << "#input " << rsi::rdfile << " " << rsi::chr << endl;
    if ( rsi::bamfile != "" )
      FOUT << "#input " << rsi::bamfile << endl;
    if ( rsi::gcadjust ) FOUT << "#GC adjusted\n";
    cnv_st header;
    header.start=header.end=-1;;
    FOUT << cnv_format1( header ) << endl;
  }
  
  for(int i=0; i<(int)cnvlist.size(); ++i) 
    FOUT << cnv_format1(cnvlist[i]) << endl;
  FOUT.close();
  
  printed=true;
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*! 
  @abstract filter cnvs
  @param  cnvlist    index for the concatenated array
  @return cnvlist    index on the reference
  @discussion appearance filter, if neighbor_runmean_sd is on the same
  scale as the global sd, the region is highly ossilitory and should be
  removed, note refsd is sd of running mean and should be much less 
  than sd of RD
*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void rd_sc_filters(vector<cnv_st>& cnvlist)  
{
  vector<bool> is_kept(cnvlist.size(), true);
  
  int minlen=max(rsi::m*2, 500);
  int minsc=(rsi::RDmedian-2*rsi::RDsd)/2/2/4;
  
  // filter_1 : get rid of highly osscilitary regions
  // note refsd is sd of runmean, should be small
  for(size_t k=0; k<cnvlist.size(); ++k) {
    if ( cnvlist[k].type==TYPE_DEL ) {
      if ( cnvlist[k].refsd > 0.5*rsi::RDsd ) is_kept[k]=false;
      if ( cnvlist[k].cnvsd > 1.5*rsi::RDsd ) is_kept[k]=false;
    }
    
    if ( cnvlist[k].type==TYPE_DUP ) {
      if ( cnvlist[k].refsd > 0.5*rsi::RDsd ) is_kept[k]=false;
      if ( cnvlist[k].cnvsd > 1.5*rsi::RDsd ) is_kept[k]=false;
    }
  }
  for(size_t k=0; k<cnvlist.size(); ++k) {
    if ( cnvlist[k].refsd > 0.5*rsi::RDsd ) is_kept[k]=false;
    if ( cnvlist[k].type==TYPE_DEL &&
	 cnvlist[k].cnvsd  > 1.5*rsi::RDsd ) is_kept[k]=false;
    if ( cnvlist[k].type==TYPE_DUP &&
	 cnvlist[k].cnvsd*rsi::RDmedian  > 2.0*cnvlist[k].cnvmed*rsi::RDsd ) is_kept[k]=false;
    if ( cnvlist[k].type==TYPE_DEL &&
	 cnvlist[k].cnvsd*rsi::RDmedian  > 2.0*cnvlist[k].cnvmed*rsi::RDsd ) is_kept[k]=false;
    if ( abs(cnvlist[k].cnvmed-cnvlist[k].refmed) <
	 3.0*cnvlist[k].refiqr ) is_kept[k]=false;
    if ( abs(cnvlist[k].cnvmed-cnvlist[k].refmed) <
	 1.25*(cnvlist[k].refiqr+cnvlist[k].cnviqr) ) is_kept[k]=false;
  }
  
  for(size_t k=0; k<cnvlist.size(); ++k) {
    if ( cnvlist[k].end-cnvlist[k].start<500 ) is_kept[k]=false;
    if ( cnvlist[k].end-cnvlist[k].start<1000 &&
	 cnvlist[k].cnvsd > rsi::RDsd  && 
	 cnvlist[k].type==TYPE_DEL ) is_kept[k]=false;
    if ( cnvlist[k].end-cnvlist[k].start<1000 &&
	 cnvlist[k].cnvsd > 2.0*rsi::RDsd  && 
	 cnvlist[k].type==TYPE_DUP ) is_kept[k]=false;
  }
  goto SUMMARIZE_FILTERS;
  
  // filter_2 : check distance between cnv and neighbor
  // if refsd not small, cnv should be well below neighbors
  for(size_t k=0; k<cnvlist.size(); ++k) {
    if ( cnvlist[k].refsd > 0.3*cnvlist[k].cnvsd &&
	 abs(cnvlist[k].refmed-cnvlist[k].cnvmed)<2.5*(cnvlist[k].cnvsd+1E-6) ) 
      is_kept[k]=false;;
  }
  
  // filter_3 : softclip reads 
  // roughly requires a minimum of 2
  // for highcov, expect one copy of chromosome(/2), half mapped correctly(/2),
  // and a 1:3 ratio between two break points(/4)
  minsc=max(2,minsc);
  for(size_t k=0; k<cnvlist.size(); ++k) {
    if ( cnvlist[k].sc1>=minsc && cnvlist[k].sc2>=minsc ) continue;
    if ( cnvlist[k].sc1<minsc && cnvlist[k].sc2<minsc ) {
      is_kept[k]=false;
      continue;
    }
  }
  
  // filter_4: length
  for(size_t k=0; k<cnvlist.size(); ++k) {
    if ( abs(cnvlist[k].end-cnvlist[k].start)<minlen &&
	 cnvlist[k].sc1 < 2*minsc && 
	 cnvlist[k].sc2 < 2*minsc ) is_kept[k]=false;
  }
  
  // salvage: perfect looking CNVs
  for(size_t k=0; k<cnvlist.size(); ++k) {
    if ( abs(cnvlist[k].end-cnvlist[k].start)<2000 ) continue; 
    if ( is_kept[k] ) continue;
    if ( cnvlist[k].type==TYPE_DEL &&
	 cnvlist[k].cnvsd < rsi::RDsd &&
	 (cnvlist[k].refmed-cnvlist[k].cnvmed)>3*cnvlist[k].cnvsd &&
	 (cnvlist[k].refmed-cnvlist[k].cnvmed)>10*cnvlist[k].refsd &&
	 cnvlist[k].cnvsd > 5*cnvlist[k].refsd ) is_kept[k]=true;
    if ( cnvlist[k].type==TYPE_DUP &&
	 cnvlist[k].cnvsd > rsi::RDsd &&
	 cnvlist[k].cnvsd < rsi::RDsd*2 &&
	 (cnvlist[k].cnvmed-cnvlist[k].refmed)>3*cnvlist[k].cnvsd &&
	 (cnvlist[k].cnvmed-cnvlist[k].refmed)>10*cnvlist[k].refsd &&
	 cnvlist[k].cnvsd > 5*cnvlist[k].refsd ) is_kept[k]=true;
  }
  
 SUMMARIZE_FILTERS:
  vector<cnv_st> tmplist(0);
  for(size_t k=0; k<cnvlist.size(); ++k) {
    //    if ( ! is_kept[k] ) tmplist.push_back( cnvlist[k] );
    if ( is_kept[k] ) tmplist.push_back( cnvlist[k] );
  }

  cnvlist=tmplist;
  return;
}

/* perl code
    $chrsd=$chrsd/1.2;
    my $is_keep=1;
    
    $is_keep=0 if $cell[2]-$cell[1]<1000;
    
    if ( $type eq "DEL" ) {
	$is_keep=0 if $p>0.2;
	$is_keep=0 if $refsd>0.6*$chrsd;
	$is_keep=0 if $cnvsd>1.3*$chrsd ;
	$is_keep=0 if ($cnvsd/$chrsd)>2.5*($cnv/$chr) ;
	$is_keep=1 if $cnv<0.66*($ref,$chr)[$ref>$chr] && $cnvsd<$chrsd;
    }
    
    if ( $type eq "DUP" ) {
	$is_keep=0 if $p>0.05;
	$is_keep=0 if $refsd>0.6*$chrsd;
	$is_keep=0 if ($cnvsd/$chrsd)>2.0*($cnv/$chr) ;
    }
    
    $is_keep=0 if $cell[2]-$cell[1]<500;
*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void sd_filters(vector<cnv_st>& cnvlist)  
{
  
  int minlen=max(rsi::m*2, 500);
  double target_sd=rsi::RDsd/1.2;

  vector<bool> is_kept(cnvlist.size(), true);
  
  for(size_t k=0; k<cnvlist.size(); ++k) {
    if ( abs(cnvlist[k].end-cnvlist[k].start)<1000 ) is_kept[k]=false;
    
    if ( cnvlist[k].type==TYPE_DEL ) {
      if ( cnvlist[k].p1 > 0.2 ) is_kept[k]=false;
      if ( cnvlist[k].refsd > 0.6*target_sd ) is_kept[k]=false;
      if ( cnvlist[k].cnvsd > 1.3*target_sd ) is_kept[k]=false;
      if ( cnvlist[k].cnvsd*rsi::RDmedian > 2.5*cnvlist[k].cnvmed*target_sd ) 
	is_kept[k]=false;
      if ( cnvlist[k].cnvmed < 0.66*min(rsi::RDmedian, cnvlist[k].refmed) &&
	   cnvlist[k].cnvsd < target_sd && 
	   abs(cnvlist[k].end-cnvlist[k].start)>800 ) is_kept[k]=true ;
    }
    
    if ( cnvlist[k].type==TYPE_DUP ) {
      if ( cnvlist[k].p1 > 0.05 ) is_kept[k]=false;
      if ( cnvlist[k].refsd > 0.6*target_sd ) is_kept[k]=false;
      if ( cnvlist[k].cnvsd*rsi::RDmedian > 2.0*cnvlist[k].cnvmed*target_sd ) 
	is_kept[k]=false;
    }
    
    if ( abs(cnvlist[k].end-cnvlist[k].start)<minlen ) is_kept[k]=false;
  }
  
  vector<cnv_st> tmplist(0);
  for(size_t k=0; k<cnvlist.size(); ++k) {
    if ( is_kept[k] ) tmplist.push_back( cnvlist[k] );
  }
  
  cnvlist=tmplist;
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void detectcnv(Array<int>& RD,  vector<cnv_st>& cnvlist )
{
  
  if ( rsi::debug ) rsi::dout << procpidstatus( rsi::pidfile, "vm");
  
  cnvlist.clear();
  cnv_st icnv;
  int i;
  
  rsi::dout << "region  : " << rsi::target_name[rsi::tid] << ":" << rsi::start << "-" << rsi::end << endl
	    << "median  : " <<  rsi::RDmedian << endl
	    << "rs::m   : " << rsi::m << endl
	    << "rs::cap : " << rsi::cap << endl;
  
  if ( rsi::RDmedian<5 ) {
    rsi::dout << "Read depths too low, cannot call" << endl;
    return;
  }
  if ( rsi::RDmedian<10 ) rsi::dout << "Read depths low, not reliable" << endl;
  
  /* median transformation */
  Array<float> RDmed(RD.size()/rsi::m, -1.0);
  Array<int> RDmedint(RDmed.size());  // int version of RDmed
  median_transfer(RD, rsi::m, RDmed);
  for(i=0;i<RDmed.size();++i) RDmedint[i]=(int)(RDmed[i]+0.5);
  
  rsi::RDmedian = _median( &RD[0], RD.size() );
  // rsi::RDsd=sqrt(variance(RDmed,0,RDmed.size()-1,0.0,-1));
  
  /* negative binomial transformation */
  Array<float> RDnbn(RD.size()/rsi::m,-1.0);
  negative_binomial_transfer(RD, rsi::m, RDnbn);
  
  /* search length */
  rsi::factor=sqrt(2.0*(1.0+rsi::epsilon)*log(3.1E9));
  rsi::Lmax = 10000/rsi::m;
  if ( rsi::Lmax<20 ) rsi::Lmax=20;
  
  /* rsi::trans=="MED" */
  Array<int> RDtrans_status(RD.size()/rsi::m,0);
  vector<cnv_st> rsiseglist;
  rsiseglist.clear();
  if ( rsi::trans!="NBN" ) {
    rsicnvmed(RDmed,RDmedint,RDtrans_status,rsiseglist);
    areblockscnv(RDmedint,RDtrans_status,rsiseglist);
  }
  
  /* rsi::trans=="NBN" */
  Array<int> RDtrans_status_nbn(RD.size()/rsi::m,0);
  vector<cnv_st> rsiseglistnbn;
  if ( rsi::trans=="NBN" ) {
    rsicnvnbn(RDnbn,RDmedint,RDtrans_status_nbn,rsiseglistnbn);
    areblockscnv(RDmedint,RDtrans_status_nbn,rsiseglistnbn);
    rsiseglist=rsiseglistnbn;
  }
  
  /* rsi::trans=="ALL" */
  if ( rsi::trans=="ALL" ) {
    vector<cnv_st> rsiseglistnbn;
    rsicnvnbn(RDnbn,RDmedint,RDtrans_status_nbn,rsiseglistnbn);
    areblockscnv(RDmedint,RDtrans_status_nbn,rsiseglistnbn);
    rsiseglist.insert(rsiseglist.end(),
		      rsiseglistnbn.begin(),rsiseglistnbn.end());
  }
  
  sortcnvstartposition(rsiseglist);
  
  /* convert coordinates from bins to bases */
  cnvlist.clear();
  for(i=0;i<(int)rsiseglist.size();++i) {
    if ( rsiseglist[i].geno == 0 ) continue;
    if ( rsiseglist[i].start==rsiseglist[i].end) continue;
    
    rsiseglist[i].start=(rsiseglist[i].start)*rsi::m+rsi::m/2;
    rsiseglist[i].end=(rsiseglist[i].end)*rsi::m+rsi::m/2;
    if ( rsiseglist[i].start < 0 ) rsiseglist[i].start=0;
    if ( rsiseglist[i].end > RD.size()-1 ) rsiseglist[i].end=RD.size()-1;
    rsiseglist[i].length=rsiseglist[i].end-rsiseglist[i].start+1;
    
    cnvlist.push_back(rsiseglist[i]);
  }
  for(i=0;i<(int)cnvlist.size();++i) cnvlist[i].tid=rsi::tid;
  
  /* optimize the start and end positions of the possible cnvs based on RD
   * by maximizing sqrt(length)*mean/sqrt(rsi::chklen*cnvvar+refvar)
   */
  optimize_with_derivative(RD,cnvlist);
  optimize_with_derivative(RD,cnvlist);
  
  rsi::dout << "Selected " << cnvlist.size() << " segments for testing" << endl;
  
  /* merge overlapped cnv segments or near by segments */
  /*
  vector<cnv_st> cnvdel,cnvadd;
  cnvdel.clear();
  cnvadd.clear();
  for(i=0;i<(int)cnvlist.size();++i) {
    if ( cnvlist[i].status==-9 ) continue;
    if ( cnvlist[i].type==TYPE_DEL ) cnvdel.push_back(cnvlist[i]);
    if ( cnvlist[i].type==TYPE_DUP ) cnvadd.push_back(cnvlist[i]);
  }
  sortcnvstartposition(cnvdel);
  mergesegments(RD, cnvdel);
  sortcnvstartposition(cnvadd);
  mergesegments(RD, cnvadd);
  cnvlist=cnvadd;
  cnvlist.insert(cnvlist.end(),cnvdel.begin(),cnvdel.end());
  sortcnvstartposition(cnvlist);
  */
  sortcnvstartposition(cnvlist);
  mergesegments(RD, cnvlist);
  sortcnvstartposition(cnvlist);
  
  /* final filtering with rsi::p and  rsi score */
  vector<cnv_st> finallist(0); 
  for(i=0;i<(int)cnvlist.size();++i) {
    double len=double(cnvlist[i].end-cnvlist[i].start+1)/double(rsi::m);
    isitcnvwrap(RD, cnvlist, i);
    cnvlist[i].score=(cnvlist[i].cnvmed-rsi::RDmedian)*sqrt(len);
    // remove cnv overlaps with N region
    int p1=expand_coordinate(cnvlist[i].start);
    int p2=expand_coordinate(cnvlist[i].end);
    for(int k=0; k<(int)rsi::noncodelist.size(); ++k) 
      if ( max(p1,rsi::noncodelist[k].start) <= min(p2,rsi::noncodelist[k].end) )
	cnvlist[i].status=-9;
    
    if ( cnvlist[i].status!=-9 )finallist.push_back( cnvlist[i] );
  }
  cnvlist=finallist;
  
  // convert coordinate back onto reference
  for(i=0;i<(int)cnvlist.size();++i) {
    int p1=expand_coordinate(cnvlist[i].start);
    int p2=expand_coordinate(cnvlist[i].end);
    cnvlist[i].start=p1;
    cnvlist[i].end=p2;
  }
  
  int nadd=0,ndel=0;
  for(i=0;i<(int)cnvlist.size();++i) {
    if ( cnvlist[i].type==TYPE_DUP ) nadd++;
    if ( cnvlist[i].type==TYPE_DEL ) ndel++;
  }
  
  rsi::dout << "Done checking overlaps, after merging, " 
	    << cnvlist.size() << " segments left" << endl;
  rsi::dout << "Found " << ndel+nadd << " CNVs : " 
	    << ndel << " DEL + " << nadd << " DUP" << endl;
  
  return;
} // void detectcnv( vector<cnv_st>& cnvlist )


/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
int rsi_usage() 
{
  cerr << "Usage:\n"
       << "\n1. detect CNV\n\n"
       << "   rsicnv rsi <options> [-b BAMFILE | -d RDFILE -c RNAME ] -f REFFILE \n"
       << "\nOptions:\n"
       << "   -m   INT  bin size, default=101\n"
       << "   -q   INT  minimum mapping quality, default=0\n"
       << "   -Q   INT  minimum base quality, default=10\n"
       << "   -cap INT  cap read depth at INT*median, dafault=4\n"
       << "             if INT<0, do not cap read depth\n"
       << "   -NOGC     do not adjust GC content, default=adjust\n"
       << "   -MED      only use median transformation \n"
       << "   -NB       only use negative binomial transformation (default)\n"
       << "   -s        save read depths before capping and GC adjustment to output and exit \n"
       << "   -o   STR  output file, default=rsiout.txt \n"
       << "   -p   STR  output plot folder, default=cnv_plots\n"
       << "   -np       do not plot CNV\n"
       << "\nNote:\n"
       << "   The reference file that was used for mapping is required.\n"
       << "   Input can be given either as a BAM file or a read depth file. In the latter case, the chromosome name must be given with -c RNAME. The read depth file can be generated with samtools mpileup BAM | cut -f2,4 , keeping only the position and read depth fields. If -c RNAME is given with -b BAMFILE, only chromosome RNAME will be processed. By default, samtools mpileup does not count duplicated reads, bases with base quality less than 13, or paired reads far apart. Rsicnv adopts the same rules except the last one.\n"
       << "\nExample:\n"
       <<"    rsicnv rsi -f b37.fasta -b $BAM -o all.rsi\n"
       <<"    rsicnv rsi -f b37.fasta -b $BAM -c $CHR -o $CHR.rsi\n"
       <<"    rsicnv rsi -f b37.fasta -d $RD  -c $CHR -o $CHR.rsi\n"
       << "\n2. plot CNV\n\n"
       << "   rsicnv plot <options> [-b BAMFILE | -d RDFILE -c RNAME ] -f REFFILE -v CNVFILE \n"
       << "\nOptions:\n"
       << "   -p   STR  output folder, default=cnv_plots\n"
       << "\nNote:\n"
       << "Rsicnv requires gnuplot to plot the figures. The figures are saved in ./cnv_plots directory in postscript format. If ImageMagic is available, the ps files will be converted to png format. GC contents are not adjusted.\n"
       << "\nExample:\n"
       <<"    rsicnv plot -f b37.fasta -b $BAM -v $CNV -p $PLOTFOLDER\n"
       << endl;
  
  return 0;
}
int get_parameters(int argc, char* argv[])
{
#define _next2 ARGV[i]=""; ARGV[i+1]=""; continue;
#define _next1 ARGV[i]=""; continue;
  
  vector<string> ARGV(0);
  for(int i=0;i<argc;++i) ARGV.push_back(string(argv[i]));
  
  if ( ARGV.size()<2 ) exit( rsi_usage() );
  if ( ARGV[1][0] != '-' ) {
    rsi::function=ARGV[1];
    ARGV[1]="";
    if ( rsi::function != "rsi" &&
	 rsi::function != "plot" && 
	 rsi::function != "stat" && 
	 rsi::function != "pin" ) {
      cerr << "no such function " << rsi::function << endl;
      exit( rsi_usage() );
    }
  }
  
  for(int i=1;i<(int)ARGV.size();++i) {
    if ( ARGV[i]=="-d" ) { rsi::rdfile=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-b" ) { rsi::bamfile=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-f" ) { rsi::reffile=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-v" ) { rsi::cnvfile=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-o" ) { rsi::outfile=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-c" ) { rsi::chr=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-s" ) { rsi::saverd=true; _next1; }
    if ( ARGV[i]=="-m" ) { rsi::m=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-q" ) { rsi::minq=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-Q" ) { rsi::min_baseQ=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-L" ) { rsi::L=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-p" ) { rsi::plotfolder=ARGV[i+1]; _next2; }
    if ( ARGV[i]=="-np" ) { rsi::plot=false; _next1; }
    if ( ARGV[i]=="-threshold" ) { rsi::threshold=atof(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-e" ) { rsi::epsilon=atof(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-cap" ) { rsi::cap=atof(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-reflen" ) { rsi::chklen=atof(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-maxchkbp" ) { rsi::maxchkbp=atoi(ARGV[i+1].c_str()); _next2; }
    if ( ARGV[i]=="-debug" ) { rsi::debug=true; _next1; }
    if ( ARGV[i]=="-MED" ) { rsi::trans="MED"; _next1; }
    if ( ARGV[i]=="-NB" ) { rsi::trans="NBN"; _next1; }
    if ( ARGV[i]=="-ALL" ) { rsi::trans="ALL"; _next1; }
    if ( ARGV[i]=="-nomerge" ) { rsi::merge=false; _next1; }
    if ( ARGV[i]=="-hist" ) { rsi::histonly=true; _next1; }
    if ( ARGV[i]=="-overlap" ) { rsi::overlaponly=true; _next1; }
    if ( ARGV[i]=="-combine" ) { rsi::combine=true; _next1; }
    if ( ARGV[i]=="-nocode" ) { rsi::nocodeonly=true; _next1; }
    if ( ARGV[i]=="-NOGC" ) { rsi::gcadjust=false; _next1; }
  }
  
  for(int i=1;i<(int)ARGV.size();++i) 
    if ( ARGV[i]!="" ) rsi::dout << "unknown option " << ARGV[i] << endl;
  for(int i=1;i<(int)ARGV.size();++i) 
    if ( ARGV[i]!="" ) exit( rsi_usage() );
  
  if ( rsi::rdfile=="" && rsi::bamfile=="" ) {
    rsi::dout << "need input file " << endl;
    exit( rsi_usage() );
  }
  if ( rsi::reffile=="" && rsi::function=="rsi" ) {
    rsi::dout << "need reference file " << endl;
    exit( rsi_usage() );
  }
  if ( rsi::outfile==rsi::bamfile || rsi::outfile==rsi::rdfile ) {
    rsi::dout << "output file is same as input file " << endl;
    exit( rsi_usage() );
  }
  if ( rsi::rdfile != "" && rsi::chr=="" ) {
    rsi::dout << "readdepth file and chromosome must be specified together" << endl;
    exit( rsi_usage() );
  }
  
  
  if ( (rsi::m % 2)!=1 ) {
    rsi::m+=1;
    rsi::dout << "m is changed to " << rsi::m << endl; 
  }
  if ( rsi::function=="plot" ) { rsi::gcadjust=false; rsi::cap=-1; }
  
  return 0;
}
int main(int argc, char* argv[])
{
  int i;
  string mycommand="",RNAME;
  
  rsi::pid = getpid();
  rsi::pidfile="/proc/"+to_string(rsi::pid)+"/status";
  
  get_parameters(argc, argv);
  
  rsi::logfile=rsi::outfile+".log";
  rsi::fout.open(rsi::logfile.c_str() ); 
  
  for(i=0;i<argc;++i) mycommand+=string(argv[i])+" ";
  
  ostringstream oss;
  oss << "#command:   " << mycommand << "\n"
      << "#bamfile:   " << rsi::bamfile << "\n"
      << "#rdfile:    " << rsi::rdfile << "\n"
      << "#reffile:   " << rsi::reffile << "\n"
      << "#cnvfile:   " << rsi::cnvfile << "\n"
      << "#chrom:     " << rsi::chr << "\n"
      << "#min_mapq:  " << rsi::minq << "\n"
      << "#min_baseQ: " << rsi::min_baseQ << "\n"
      << "#binsize:   " << rsi::m << "\n"
      << "#adjustGC:  " << rsi::gcadjust << "\n"
      << "#plots:     " << plot::folder << "\n" 
      << "#output:    " << rsi::outfile << "\n" ;
  
  rsi::dout << oss.str();
  
  vector<int> populated(0);
  Array<int> RD(1);
  vector <cnv_st> cnvlist(0);   
  
  //! open bam file and get a list of chromosomes
  samfile_t *fp_in = NULL;
  bam1_t *b=NULL;   
  bam_index_t *bamidx=NULL;
  bam_iter_t iter=0;
  b = bam_init1();
  int ref=0, beg=0, end=0x7fffffff;
  if ( rsi::bamfile != "" ) {
    fp_in = samopen(rsi::bamfile.c_str(), "rb", 0);
    bamidx = bam_index_load(rsi::bamfile.c_str()); // load BAM index
    rsi::dout << "#Check bam header for 1-22XY \n";
    rsi::target_name.clear();
    for(int i=0;i<fp_in->header->n_targets;++i) {
      RNAME=string(fp_in->header->target_name[i]);
      rsi::target_name.push_back(RNAME);
      if ( RNAME.find("MT") != string::npos ) continue;
      if ( RNAME.find(".") != string::npos ) continue;
      //if ( RNAME.find("X") != string::npos ) continue;
      //if ( RNAME.find("Y") != string::npos ) continue;
      bam_parse_region(fp_in->header, RNAME.c_str(), &ref, &beg, &end); 
      iter = bam_iter_query(bamidx, ref, beg, end);
      if ( bam_iter_read(fp_in->x.bam, iter, b)>0 ) populated.push_back(ref);
      rsi::dout << RNAME << "\t" << fp_in->header->target_len[i] << "\t" << ref << "\t" << beg << "\t" << end << endl;
    }
    rsi::dout << "#BAM has reads on :";
    for(int i=0;i<(int)populated.size();++i) 
      rsi::dout << " " << rsi::target_name[populated[i]];
    rsi::dout << endl;
  }
  else {
    rsi::target_name.clear();
    rsi::target_name.push_back(rsi::chr);
  }
  if ( rsi::chr != "1-22XY" ) {
    populated.clear();
    for(int i=0;i<(int)rsi::target_name.size();++i) 
      if ( rsi::target_name[i]==rsi::chr ) {
	populated.push_back(i);
	break;
      }
    if ( populated.empty() ) {
      cerr << "BAM file doesn't have " << rsi::chr << endl;
      exit(0);
    }
  }
  
  //! plot cnv
  if ( rsi::function=="plot" ) {
    rsi::tid=-1;
    if ( gnuplot_version() < 0 ) {
      cerr << "gnuplot not found" << endl;
      goto CLOSEBAM;
    }
    read_cnvlist(rsi::cnvfile, cnvlist);
    if ( cnvlist.size()==0 ) {
      cerr << "no cnv found" << endl;
      goto CLOSEBAM;
    }
    if ( fp_in ) {
      plot_cnv_from_bam(cnvlist, fp_in, bamidx);
      goto CLOSEBAM;
    }
    
    for(int k=0; k<(int)cnvlist.size(); ++k) {
      if ( cnvlist[k].tid<0 ) continue;
      if ( cnvlist[k].tid!=rsi::tid ) {
	if ( rsi::rdfile!="" ) {
	  rsi::tid=0;
	  load_data_from_text(rsi::rdfile,RD);
	}
      }
      
      if ( cnvlist[k].tid==rsi::tid ) {
	vector<cnv_st> tmplist(0);
	tmplist.clear();
	tmplist.push_back( cnvlist[k] );
	plot_cnv(tmplist, RD);
      }
    }
    
    goto CLOSEBAM;
  }
  
  //! detect CNV
  if ( rsi::function=="rsi" ) {
    for(int i=0;i<(int)populated.size();++i) {
      rsi::dout << "#processing " << rsi::target_name[populated[i]] << endl;
      if ( rsi::rdfile!="" ) {
	rsi::tid=0;
	load_data_from_text(rsi::rdfile,RD);
      }
      else {
	rsi::tid=populated[i];
	load_data_from_bam(fp_in, bamidx, rsi::tid, RD);
      }
      
      concatenate_data(RD);
      
      rsi::RDmedian = _median( &RD[0], RD.size() );
      rsi::RDsd=sqrt(variance(RD,0,RD.size()-1,0.0,-1));
      
      cnvlist.clear();
      detectcnv(RD,cnvlist);
      
      sd_filters(cnvlist);
      if ( fp_in ) cnv_stat(fp_in, bamidx, cnvlist);
      
      write_cnv_to_file(cnvlist, rsi::outfile);
      rsi::dout << "output written to " << rsi::outfile << endl;
      if ( gnuplot_version()>0 && rsi::plot ) {
	expand_data(RD); // add back N regions
	plot_cnv(cnvlist, RD);
      }
    }
    goto CLOSEBAM;
  }
  
  //! improve precision based on other info
  if ( rsi::function=="pin" ) {
    cnvlist.clear();
    read_cnvlist(rsi::cnvfile, cnvlist);
    optimize_resolution(fp_in, bamidx, cnvlist);
    write_cnv_to_file(cnvlist, rsi::outfile);
    goto CLOSEBAM;
  }
  
  //! get basic cnv region info
  if ( rsi::function=="stat" ) {
    cnvlist.clear();
    string text; vector<size_t> text_idx;
    read_cnvlist(rsi::cnvfile, cnvlist, text, text_idx);
    cnv_stat(fp_in, bamidx, cnvlist);
    ofstream FOUT(rsi::outfile.c_str());
    for(size_t i=0; i<cnvlist.size(); ++i) {
      size_t pos=text.find("\n",text_idx[i]);
      FOUT << text.substr(text_idx[i], pos-text_idx[i]) << "\t"
	   << "RP=" << cnvlist[i].rp << ";"
	   << "Q0=" << cnvlist[i].q0 
	   << endl;
    }
    FOUT.close();
    goto CLOSEBAM;
  }
  
 CLOSEBAM:
  if ( fp_in) samclose(fp_in);
  if ( b ) bam_destroy1(b);
  if ( iter ) bam_iter_destroy(iter);
  if ( bamidx) bam_index_destroy(bamidx);
  
  rsi::dout << "\n" << oss.str() << endl;
  rsi::dout << procpidstatus( rsi::pidfile, "vm");  
  rsi::dout << "exit" << endl;  
  
  exit(0);
} 
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
