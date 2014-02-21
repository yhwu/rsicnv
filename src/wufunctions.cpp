/*************************************************************************
BEGIN OF LICENSE
Copyright (c) Yinghua Wu and Hongzhe Li (rsicnv project).

This program is free software; you can redistribute and/or modify
the codes written by the authors under the terms of the GNU General 
Public License as published by the Free Software Foundation 
(www.fsf.org); either version 2 of the License, or (at your option) 
any later version. 

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
END OF LICENSE

CONTACT: wu_yinghua@hotmail.com; hongzhe@upenn.edu
*************************************************************************/
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <complex>
#include <string>
#include <limits>
#include <algorithm>
using namespace std;

#include "wu2.h"  //general array matrix headers
#include "alglibinterface.h"

#include "wufunctions.h"  // functions defined here

string procpidstatus(string file, string fields)
{
  string value="";
  ifstream inp(file.c_str());
  if ( !inp ) {
    cout << "file " << file << " does not exist" << endl;
    return value;
  }
  
  while ( !inp.eof() ) {
    string tmps;
    getline(inp,tmps);
    if ( ci_find(tmps, fields ) != string::npos ) value+=tmps+"\n";
  }
  inp.close();
  return(value);
}
string procpidstatus(int pid, string fields)
{
  string value="";
  
  string file="/proc/"+to_string(pid)+"/status";
  ifstream inp(file.c_str());
  if ( !inp ) {
    cout << "file " << file << " does not exist" << endl;
    return value;
  }
  
  while ( !inp.eof() ) {
    string tmps;
    getline(inp,tmps);
    if ( ci_find(tmps, fields ) != string::npos ) value+=tmps+"\n";
  }
  inp.close();
  return(value);
}

string to_lower(string word) {
  std::transform(word.begin(), word.end(), word.begin(), ::tolower);
  return word;
}
string to_upper(string word) {
  std::transform(word.begin(), word.end(), word.begin(), ::toupper);
  return word;
}

/* case insensitive find 
 * find the position in first string that match the second
 */
bool ci_equal_char(char ch1, char ch2)
{ return toupper((unsigned char)ch1) == toupper((unsigned char)ch2); }
size_t ci_find(const string& str1, const string& str2)
{
  string::const_iterator pos = 
    std::search(str1.begin(), str1.end(), 
		str2.begin(), str2.end(), 
		ci_equal_char);
  if (pos == str1.end() )  return string::npos;
  else return pos-str1.begin();
}
/* case insensitive find */

/* case insensitive equal
 */
bool ci_equal(const string& str1, const string& str2)
{
  if ( str1.length() != str2.length() ) return(false); 
  for(int i=0;i<(int)str1.length();++i)
    if ( ! ci_equal_char(str1[i], str2[i]) ) return(false); 
  return(true);
}
/* case insensitive equal */


/* begin runmed 
 *
 **/
template<class Etype>
void Srunmed(Array<Etype>& y, Array<Etype>& smo, int n, int band,
	     int end_rule, int debug)
{
  Etype rmed, rmin, temp, rnew, yout, yi;
  Etype rbe, rtb, rse, yin, rts;
  int imin, ismo, i, j, first, last, band2, kminus, kplus;
  int bw;

  if ( band % 2 != 1 ) {
    cout << "band must be odd" << endl;
    cout << "band =" << band << endl;
    cout << "runmed" << endl;
    exit(0);
  }
  if( band > n) {
    cout << "band = " << band << endl
	 << "n = " << n 
	 << "bandwidth/span of running medians is larger than n" << endl
	 << endl;
    exit(0);
  }
  if ( n<= 0 ) n=y.size();

  bw=band;

  Array<Etype> scrat(bw);
  

  /* 1. Compute  'rmed' := Median of the first 'band' values
     ======================================================== */
  
  for (i = 0; i < bw; ++i ) scrat[i] = y[i];
  
  /* find minimal value  rmin = scrat[imin] <= scrat[j] */
  rmin = scrat[0]; imin = 0;
  for (i = 1; i < bw; ++i )
    if (scrat[i] < rmin) { rmin = scrat[i]; imin = i; }

  /* swap scrat[0] <-> scrat[imin] */
  temp = scrat[0]; scrat[0] = rmin; scrat[imin] = temp;
  
  /* sort the rest of of scrat[] by bubble (?) sort -- */
  for (i = 2; i < bw; ++i) {
    if (scrat[i] < scrat[i - 1]) {/* find the proper place for scrat[i] */
      temp = scrat[i];
      j = i;
      do {
	scrat[j] = scrat[j - 1];
	--j;
      } while (scrat[j - 1] > temp); /* now:  scrat[j-1] <= temp */
      scrat[j] = temp;
    }
  }
  
  band2 = bw / 2;
  rmed = scrat[band2]; /* == Median( y[(1:band2)-1] ) */
  /* "malloc" had  free( (char*) scrat);*/ /*-- release scratch memory --*/
  
  if( end_rule == 0) { /*-- keep DATA at end values */
    for (i = 0; i < band2; ++i) smo[i] = y[i];
  }
  else { /* if(*end_rule == 1)  copy median to CONSTANT end values */
    for (i = 0; i < band2; ++i) smo[i] = rmed;
  }
  smo[band2] = rmed;
  band2++; /* = bw / 2 + 1*/;
  
  if( debug ) cout << bw << "\t" << band2 << endl;
  
  /* Big	FOR Loop: RUNNING median, update the median 'rmed'
     ======================================================= */
  for (first = 1, last = bw, ismo = band2;
       last < n;
       ++first, ++last, ++ismo) {
    
    yin = y[last];
    yout = y[first - 1];
    
    if( debug ) 
      cout << "is=" << ismo << " y(in/out)= " << yin << "\t" << yout << endl;
    
    rnew = rmed; /* New median = old one   in all the simple cases --*/
	
    if (yin < rmed) {
      if (yout >= rmed) {
	kminus = 0;
	if (yout > rmed) {/*	--- yin < rmed < yout --- */
	  if( debug) cout << ": yin < rmed < yout " << endl;
	  rnew = yin;/* was -rinf */
	  for (i = first; i <= last; ++i) {
	    yi = y[i];
	    if (yi < rmed) {
	      ++kminus;
	      if (yi > rnew) rnew = yi;
	    }
	  }
	  if (kminus < band2)		rnew = rmed;
	}
	else {/*		--- yin < rmed = yout --- */
	  if( debug) cout << ": yin < rmed == yout " << endl;
	  rse = rts = yin;   /* was -rinf */
	  for (i = first; i <= last; ++i) {
	    yi = y[i];
	    if (yi <= rmed) {
	      if (yi < rmed) {
		++kminus;
		if (yi > rts)	rts = yi;
		if (yi > rse)	rse = yi;
	      } else		rse = yi;
	      
	    }
	  }
	  rnew = (kminus == band2) ? rts : rse ;
	  if(debug) cout << "k- : " << kminus << endl;
	}
      } /* else: both  yin, yout < rmed -- nothing to do .... */
    }
    else if (yin != rmed) { /* yin > rmed */
      if (yout <= rmed) {
	kplus = 0;
	if (yout < rmed) {/* -- yout < rmed < yin --- */
	  if(debug) cout << ": yout < rmed < yin " << endl;
	  rnew = yin; /* was rinf */
	  for (i = first; i <= last; ++i) {
	    yi = y[i];
	    if (yi > rmed) {
	      ++kplus;
	      if (yi < rnew)	rnew = yi;
	    }
	  }
	  if (kplus < band2)	rnew = rmed;
	  
	} else { /* -- yout = rmed < yin --- */
	  if(debug) cout << ": yout == rmed < yin " << endl;
	  rbe = rtb = yin; /* was rinf */
	  for (i = first; i <= last; ++i) {
	    yi = y[i];
	    if (yi >= rmed) {
	      if (yi > rmed) {
		++kplus;
		if (yi < rtb)	rtb = yi;
		if (yi < rbe)	rbe = yi;
	      } else		rbe = yi;
	    }
	  }
	  rnew = (kplus == band2) ? rtb : rbe;
	  if(debug) cout << "k+ : " << kplus << endl;
	}
      } /* else: both  yin, yout > rmed --> nothing to do */
    } /* else: yin == rmed -- nothing to do .... */
    if(debug) cout << "=> " << rmed << "\t" <<  rnew << endl;
    rmed = rnew;
    smo[ismo] = rmed;
  } /*     END FOR ------------ big Loop -------------------- */
  
  if(end_rule == 0) { /*-- keep DATA at end values */
    for (i = ismo; i < n; ++i)  smo[i] = y[i];
  }
  else { /* if(*end_rule == 1)  copy median to CONSTANT end values */
    for (i = ismo; i < n; ++i) smo[i] = rmed;
  }
} /* Srunmed */

void runmed(Array<double>& y, Array<double>& smo, int n, int band,
	    int end_rule, int debug)
{ Srunmed(y, smo, n, band, end_rule, debug); }
void runmed(Array<float>& y, Array<float>& smo, int n, int band,
	    int end_rule, int debug)
{ Srunmed(y, smo, n, band, end_rule, debug); }
void runmed(Array<int>& y, Array<int>& smo, int n, int band,
	    int end_rule, int debug)
{ Srunmed(y, smo, n, band, end_rule, debug); }

/* end runmed 
 *
 **/


template<class Etype>
double partition_percentile_tp(Etype *x, size_t n, double dy, double  p)
{
  if ( dy<=0 ) {
    std::cerr << "dy cannot be zero or negative for partion_median\n";
    exit(0);
  }
  if ( p<0 || p>1 ) {
    std::cerr << "0<p<1   p=" << p << "\n"
	      << "partition_percentile_tp()" << endl;
    exit(0);
  }
  
  double ymin,ymax,mean=0;
  ymin=ymax=x[0];
  for(size_t i=0;i<n;++i) {
    mean+=x[i];
    if ( x[i]<ymin ) ymin=x[i];
    if ( x[i]>ymax ) ymax=x[i];
  }
  mean/=(double)n;
  if ( (ymax-ymin)<dy ) return(mean);
  
  size_t np=(ymax-ymin)/dy+2;
  if ( np>10000000 ) std::cerr << "extremely big read depths:\t" 
			       << ymin << "\t" << ymax << "\n"
			       << "partition_stat_tp()\n" ;
  
  size_t *dat = new (std::nothrow)  size_t[np];
  if ( dat==0 ) { 
    std::cerr << "not enough memory, request size_t: " << n << "\n"
	      << "partition_percentile_tp()\n" ;
    exit(0);
  }
  for(size_t i=0;i<np;++i) dat[i]=0;
  
  for(size_t i=0;i<n;++i) {
    double idx=(x[i]-ymin)/dy+0.5;
    if ( idx<0 || idx>np ) {
      std::cerr << "idx negative" << x[i] << "\t" << ymax << "\n"
		<< "partition_percentile_tp()\n" ;
      exit(0);
    }
    dat[(size_t)idx]+=1;
  }
  
  size_t count=0;
  size_t imed=(double)n*p;
  double v_p=ymin;
  for(size_t i=0;i<np;++i) {
    if ( count < imed && count+dat[i] >= imed ) {
      v_p=ymin+i*dy;
      break;
    }
    count+=dat[i];
  }
  delete[] dat;
  
  return(v_p);
}
double partition_percentile(double *x, size_t n, double  p) 
{  return partition_percentile_tp(x, n, 0.01, p); }
double partition_percentile(float *x, size_t n, double  p) 
{  return partition_percentile_tp(x, n, 0.01, p); }
double partition_percentile(int *x, size_t n, double  p) 
{  return partition_percentile_tp(x, n, 1, p); }
double partition_percentile(size_t *x, size_t n, double  p) 
{  return partition_percentile_tp(x, n, 1, p); }



template<class Etype>
void partition_stat_tp(Etype *x, size_t n, double dy,
		       double& lqt, double& med, double& uqt)
{
  if ( dy<=0 ) {
    std::cerr << "dy cannot be zero or negative for partion_median\n";
    exit(0);
  }
  
  double ymin,ymax,mean=0;
  ymin=ymax=x[0];
  for(size_t i=0;i<n;++i) {
    mean+=x[i];
    if ( x[i]<ymin ) ymin=x[i];
    if ( x[i]>ymax ) ymax=x[i];
  }
  mean/=(double)n;
  lqt=ymin;
  med=mean;
  uqt=ymax;
  if ( (ymax-ymin)<dy ) return;
  
  size_t np=(ymax-ymin)/dy+2;
  if ( np>10000000 ) std::cerr << "extremely big read depths:\t" 
			       << ymin << "\t" << ymax << "\n"
			       << "partition_stat_tp()\n" ;
  
  size_t *dat = new (std::nothrow)  size_t[np];
  if ( dat==0 ) { 
    std::cerr << "not enough memory, request size_t: " << n << "\n"
	      << "partition_stat_tp()\n" ;
    exit(0);
  }
  for(size_t i=0;i<np;++i) dat[i]=0;
  
  for(size_t i=0;i<n;++i) {
    double idx=(x[i]-ymin)/dy+0.5;
    if ( idx<0 || idx>np ) {
      std::cerr << "idx negative" << x[i] << "\t" << ymax << "\n"
		<< "partition_stat_tp()\n" ;
      exit(0);
    }
    dat[(size_t)idx]+=1;
  }
  
  size_t count=0;
  size_t ilqt=n/4, imed=n/2, iuqt=n*3/4;
  for(size_t i=0;i<np;++i) {
    if ( count < ilqt && count+dat[i] >= ilqt ) lqt=ymin+i*dy;
    if ( count < imed && count+dat[i] >= imed ) med=ymin+i*dy;
    if ( count < iuqt && count+dat[i] >= iuqt ) uqt=ymin+i*dy;
    count+=dat[i];
  }
  if ( count!=n ) {
    cerr << "counting error\t" << count << "\t" << n << "\n" 
	 << "partition_median_tp()\n" ;
    exit(0);
  }
  
  delete[] dat;
  return;
}
double partition_interquartilerange(double* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 0.01, lqt, med, uqt);  
  return(uqt-lqt); 
}
double partition_interquartilerange(float* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 0.01, lqt, med, uqt);  
  return(uqt-lqt); 
}
double partition_interquartilerange(int* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 1, lqt, med, uqt);  
  return(uqt-lqt); 
}
double partition_interquartilerange(size_t* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 1, lqt, med, uqt);  
  return(uqt-lqt); 
}

double partition_lowerquartile(double* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 0.01, lqt, med, uqt);  
  return(lqt); 
}
double partition_lowerquartile(float* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 0.01, lqt, med, uqt);  
  return(lqt); 
}
double partition_lowerquartile(int* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 1, lqt, med, uqt);  
  return(lqt); 
}
double partition_lowerquartile(size_t* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 1, lqt, med, uqt);  
  return(lqt); 
}

double partition_upperquartile(double* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 0.01, lqt, med, uqt);  
  return(uqt); 
}
double partition_upperquartile(float* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 0.01, lqt, med, uqt);  
  return(uqt); 
}
double partition_upperquartile(int* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 1, lqt, med, uqt);  
  return(uqt); 
}
double partition_upperquartile(size_t* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 1, lqt, med, uqt);  
  return(uqt); 
}

double partition_median(double* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 0.01, lqt, med, uqt);  
  return(med); 
}
double partition_median(float* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 0.01, lqt, med, uqt);  
  return(med); 
}
double partition_median(int* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 1, lqt, med, uqt);  
  return(med); 
}
double partition_median(size_t* x, size_t n) { 
  double lqt,med,uqt;
  partition_stat_tp(x, n, 1, lqt, med, uqt);  
  return(med); 
}

/*
template<class Etype>
double partition_median_tp(Etype *x, size_t n, double dy)
{
  if ( dy<=0 ) {
    std::cerr << "dy cannot be negative for partion_median\n";
    exit(0);
  }
  
  double med,ymin,ymax;
  ymin=ymax=x[0];
  for(size_t i=0;i<n;++i) {
    if ( x[i]<ymin ) ymin=x[i];
    if ( x[i]>ymax ) ymax=x[i];
  }
  if ( ymax-ymin < dy ) return ( (ymax+ymin)*0.5 );
  
  size_t np=(ymax-ymin)/dy+2;
  size_t *dat = new (std::nothrow)  size_t[np];
  if ( dat==0 ) { 
    std::cerr << "not enough memory, request size_t: " << n << "\n"
	      << "partition_median_tp()\n" ;
    exit(0);
  }
  for(size_t i=0;i<np;++i) dat[i]=0;
  
  for(size_t i=0;i<n;++i) {
    double idx=(x[i]-ymin)/dy+0.5;
    if ( idx<0 ) {
      std::cerr << "idx negative" << x[i] << "\t" << ymax << "\n"
		<< "partition_median_tp()\n" ;
      exit(0);
    }
    dat[(size_t)idx]+=1;
  }
  
  med=ymin-100;
  size_t count=0;
  for(size_t i=0;i<np;++i) {
    if ( count < n/2 && count+dat[i] > n/2 ) {
      med=ymin+i*dy;
      // break;
    }
    count+=dat[i];
  }
  if ( count!=n ) {
    cerr << "counting error\t" << count << "\t" << n << "\n" 
	 << "partition_median_tp()\n" ;
    exit(0);
  }
  
  delete[] dat;
  return(med);
}
double partition_median(double *x, size_t n) { return partition_median_tp(x, n, 0.01);  }
double partition_median(float *x, size_t n) { return partition_median_tp(x, n, 0.01);  }
double partition_median(size_t *x, size_t n) { return partition_median_tp(x, n, 1);  }
double partition_median(int *x, size_t n) { return partition_median_tp(x, n, 1);  }
*/
/* begin runmean 
 *
 **/

template<class Etype>
void runmeantp(Array<Etype>& y, Array<double>& smo, int n, int band,int end_rule)
{
  int ismo, i, first, last, band2;
  double sum,mean;
  int bw;

  //if ( band % 2 != 1 ) {
  //  cout << "band must be odd" << endl;
  //  cout << "runmeantp()" << endl;
  //  exit(0);
  //}
  if(n > y.size() ) {
    cout << "n out of subscript of y" << endl;
    cout << "runmeantp()" << endl;
    exit(0);
  }
  if( band > n) {
    cout << "bandwidth/span of running medians is larger than n" << endl
	 << "band = " << band << endl
	 << "n = " << n 
	 << endl;
    cout << "runmeantp()" << endl;
    exit(0);
  }
  if ( y.size() > smo.size() ) {
    cout << "Array1's length is longer than Array2's length()" << endl
	 << "y.size()=" << y.size() << "t"
	 << "smo.size()=" << smo.size() << endl;
    cout << "runmeantp()" << endl;
    exit(0);
  }
  if ( n <= 0 ) {
    cout << "n is negative" << endl;
    cout << "runmeantp()" << endl;
    exit(0);
  }
  bw=band;

  sum=0;
  for (i = 0; i < bw; ++i ) sum+=(double)y[i];
  mean=sum/double(bw);

  band2 = bw / 2;
  if( end_rule == 0) { /*-- keep DATA at end values */
    for (i = 0; i < band2; ++i) smo[i] = (double)y[i];
  }
  else { /* if(*end_rule == 1)  copy median to CONSTANT end values */
    for (i = 0; i < band2; ++i) smo[i] = mean;
  }
  smo[band2] = mean;
  band2++; /* = bw / 2 + 1*/;
  
  for (first = 1, last = bw, ismo = band2;
       last < n;
       ++first, ++last, ++ismo) {
    sum = sum - (double)y[first-1] + (double)y[last];
    mean = sum / double(bw);
    smo[ismo] = mean;
  } 
  
  if(end_rule == 0) { /*-- keep DATA at end values */
    for (i = ismo; i < n; ++i)  smo[i] = (double)y[i];
  }
  else { /* if(*end_rule == 1)  copy median to CONSTANT end values */
    for (i = ismo; i < n; ++i) smo[i] = mean;
  }
} /* runmeantp */
template<class Etype>
void runmeantp(Array<Etype>& y, Array<float>& smo, int n, int band,int end_rule)
{
  Array<double> smo1(n);
  for(int i=0;i<n;++i) smo1[i]=smo[i];
  runmeantp(y, smo1, n, band, end_rule);
  for(int i=0;i<n;++i) smo[i]=smo1[i];
} /* runmeantp */
void runmean(Array<double>& y, Array<double>& smo, int n, int band,
	     int end_rule)
{ runmeantp(y, smo, n, band, end_rule); }
void runmean(Array<float>& y, Array<float>& smo, int n, int band,
	     int end_rule)
{ runmeantp(y, smo, n, band, end_rule); }
void runmean(Array<int>& y, Array<float>& smo, int n, int band,
	     int end_rule)
{ runmeantp(y, smo, n, band, end_rule); }

/* end runmean 
 *
 **/

/* begin mean 
 *
 **/
template<class Etype>
double mean_tp( Array<Etype>& y, int low, int high )
{
  int i,bw;
  double sum;
  
  if ( low < 0 ) {
    cout << "low < 0, low = " << low << endl
	 << "meantp()" << endl;
    exit(0);
  }
  if(high > y.size()-1 ) {
    cout << "high out of subscript of y" << endl;
    cout << "high=" << high << "\t" << "y.size=" << y.size() << endl;
    cout << "meantp()" << endl;
    exit(0);
  }

  bw=high-low+1;

  sum=0;
  for (i = low; i <= high; ++i ) sum+=(double)y[i];
  return( sum/double(bw) );

} /* mean */
double mean(Array<double>& y, int low, int high) {return(mean_tp(y,low,high));}
double mean(Array<float>& y, int low, int high) {return(mean_tp(y,low,high));}
double mean(Array<int>& y, int low, int high) {return(mean_tp(y,low,high));}
/* end mean 
 *
 **/

/* begin min 
 *
 **/
template<class Etype>
Etype mintp( Array<Etype>& y, int low, int high)
{
  int i;
  Etype min;

  if ( low < 0 ) {
    cout << "low < 0, low = " << low << endl
	 << "mintp()" << endl;
    exit(0);
  }
  if(high > y.size()-1 ) {
    cout << "high out of subscript of y" << endl;
    cout << "high=" << high << "\t" << "y.size=" << y.size() << endl;
    cout << "mintp()" << endl;
    exit(0);
  }

  min=y[low];
  for (i = low; i <= high; ++i ) if ( y[i] < min ) min=y[i];
  return( min );
} /* min */
double min(Array<double>& y, int low, int high) {return(mintp(y,low,high));}
float min(Array<float>& y, int low, int high) {return(mintp(y,low,high));}
int min(Array<int>& y, int low, int high) {return(mintp(y,low,high));}
/* end min 
 *
 **/

/* begin max 
 *
 **/
template<class Etype>
Etype maxtp( Array<Etype>& y, int low, int high)
{
  int i;
  Etype max;

  if ( low < 0 ) {
    cout << "low < 0, low = " << low << endl
	 << "maxtp()" << endl;
    exit(0);
  }
  if(high > y.size()-1 ) {
    cout << "high out of subscript of y" << endl;
    cout << "high=" << high << "\t" << "y.size=" << y.size() << endl;
    cout << "maxtp()" << endl;
    exit(0);
  }

  max=y[low];
  for (i = low; i <= high; ++i ) if ( y[i] > max ) max=y[i];
  return( max );
} /* max */
double max(Array<double>& y, int low, int high) {return(maxtp(y,low,high));}
float max(Array<float>& y, int low, int high) {return(maxtp(y,low,high));}
int max(Array<int>& y, int low, int high) {return(maxtp(y,low,high));}
/* end max 
 *
 **/


/* begin variance 
 *
 **/
template<class Etype>
double variancetp( Array<Etype>& y, int low, int high, double d, int end_rule )
{
  /*
   * d, data outside d * standard dev is filtered out
   * end_rule == -1, no filtering
   * end_rule == 0, both sides
   * end_rule == 1, lower sides
   * end_rule == 2, higher sides
   *
   */

  int i,bw;
  double sum=0.0,sum2=0.0,v=0.0,std=0.0,mean=0.0;
  
  if ( low < 0 ) {
    cout << "low is negative" << endl;
    cout << "low=" << low << endl;
    cout << "variance2tp()" << endl;
    exit(0);
  }
  if(high > y.size()-1 ) {
    cout << "high out of subscript of y" << endl;
    cout << "high=" << high << "\t" << "y.size=" << y.size() << endl;
    cout << "variance2tp()" << endl;
    exit(0);
  }
  if( d < 0 ) {
    cout << "d is supposed to be positive" << endl;
    cout << "using d=" << abs(d) << endl;
    cout << "variance2tp()" << endl;
    d=abs(d);
    exit(0);
  }
  bw=high-low+1;
  
  sum=0;
  for (i = low; i <= high; ++i ) {
    sum+=(double)y[i];
    sum2+=(double)y[i]*(double)y[i];
  }
  mean=sum/double(bw);
  v=sum2/double(bw)-mean*mean;
  std=sqrt(v);
  if ( end_rule == -1 ) return(v);

  int del_end1=0,del_end2=0;
  double offset_end1,offset2_end1;
  double offset_end2,offset2_end2;
  offset_end1=offset2_end1=0.0;
  offset_end2=offset2_end2=0.0;
  for (i = low; i <= high; ++i ) {
    if ( y[i] < mean-d*std ) {
      del_end1++;
      offset_end1+=(double)y[i];
      offset2_end1+=(double)y[i]*(double)y[i];
    }
    if ( y[i] > mean+d*std ) {
      del_end2++;
      offset_end2+=(double)y[i];
      offset2_end2+=(double)y[i]*(double)y[i];
    }
  }

  if ( del_end1+del_end2 >= bw ) {
    cout << "filtering too many points" << endl;
    cout << "variance2tp()" << endl;
    exit(0);
  }
 
  if ( end_rule == 0 ) {
    sum=sum-offset_end1-offset_end2;
    sum2=sum2-offset2_end1-offset2_end2;
    mean=sum/double(bw-del_end1-del_end2);
    v=sum2/double(bw-del_end1-del_end2)-mean*mean;
    return(v);
  }

  if ( end_rule == 1 ) {
    sum=sum-offset_end1;
    sum2=sum2-offset2_end1;
    mean=sum/double(bw-del_end1);
    v=sum2/double(bw-del_end1)-mean*mean;
    return(v);
  }

  if ( end_rule == 2 ) {
    sum=sum-offset_end2;
    sum2=sum2-offset2_end2;
    mean=sum/double(bw-del_end2);
    v=sum2/double(bw-del_end2)-mean*mean;
    return(v);
  }
  cout << "variancetp(): it shouldn't get here" << endl;
  return(-1);
} /* variancetp */
double variance(Array<double>& y, int low, int high, double d, int end_rule) 
{return(variancetp(y,low,high,d,end_rule));}
double variance(Array<float>& y, int low, int high, double d, int end_rule) 
{return(variancetp(y,low,high,d,end_rule));}
double variance(Array<int>& y, int low, int high, double d, int end_rule) 
{return(variancetp(y,low,high,d,end_rule));}
/* end variance 
 *
 **/

/* begin remove_outlier 
 *
 **/
template<class Etype>
void remove_outliertp(Array<Etype>& y, Array<Etype>& smo, int low, int high, double d, int end_rule )
{
  /*
   * d, data outside d * standard dev is filtered out
   * end_rule == -1, no filtering
   * end_rule == 0, both sides
   * end_rule == 1, lower sides
   * end_rule == 2, higher sides
   *
   */

  int i,bw;
  double sum=0.0,sum2=0.0,v=0.0,std=0.0,mean=0.0;
  
  if ( low < 0 ) {
    cout << "low is negative" << endl;
    cout << "low=" << low << endl;
    cout << "remove_outliertp()" << endl;
    exit(0);
  }
  if(high > y.size() ) {
    cout << "high out of subscript of y" << endl;
    cout << "high=" << high << "\t" << "y.size=" << y.size() << endl;
    cout << "remove_outliertp()" << endl;
    exit(0);
  }
  if( d < 0 ) {
    cout << "d is supposed to be positive" << endl;
    cout << "using d=" << abs(d) << endl;
    cout << "remove_outliertp()" << endl;
    d=abs(d);
    exit(0);
  }
  bw=high-low+1;

  Array<Etype> tmp(bw);

  if ( end_rule == -1 ) {
    smo.resize(bw);
    for(i=low;i<=high;++i) smo[i-low]=y[i];
    return;
  }
  
  sum=0;
  for (i = low; i <= high; ++i ) {
    sum+=(double)y[i];
    sum2+=(double)y[i]*(double)y[i];
  }
  mean=sum/double(bw);
  v=sum2/double(bw)-mean*mean;
  std=sqrt(v);

  if ( end_rule==0 ) {
    int j=0;
    for(i=low;i<=high;++i) 
      if ( y[i] >= mean-d*std  && y[i] <= mean+d*std ) { tmp[j]=y[i]; ++j; }
    smo.resize(j);
    for(i=0;i<j;++i) smo[i]=tmp[i];
    return;
  }
  if ( end_rule==1 ) {
    int j=0;
    for(i=low;i<=high;++i) 
      if ( y[i] >= mean-d*std  ) { tmp[j]=y[i]; ++j; }
    smo.resize(j);
    for(i=0;i<j;++i) smo[i]=tmp[i];
    return;
  }
  if ( end_rule==2 ) {
    int j=0;
    for(i=low;i<=high;++i) 
      if ( y[i] <= mean+d*std  ) { tmp[j]=y[i]; ++j; }
    smo.resize(j);
    for(i=0;i<j;++i) smo[i]=tmp[i];
    return;
  }
  cout << "remove_outliertp(): it shouldn't get here" << endl;
  return;
} /* remove_outliertp */
void remove_outlier(Array<double>& y, Array<double>& smo, int low, int high, double d, int end_rule) 
{remove_outliertp(y,smo,low,high,d,end_rule);}
void remove_outlier(Array<float>& y, Array<float>& smo, int low, int high, double d, int end_rule) 
{remove_outliertp(y,smo,low,high,d,end_rule);}
void remove_outlier(Array<int>& y, Array<int>& smo, int low, int high, double d, int end_rule) 
{remove_outliertp(y,smo,low,high,d,end_rule);}
/* end remove_outlier 
 *
 **/

/* begin median 
 *
 */
template < class Etype >
void insertionSort(Etype arr[], int length) 
{
  int i, j, tmp;
  for (i = 1; i < length; i++) {
    j = i;
    while (j > 0 && arr[j - 1] > arr[j]) {
      tmp = arr[j];
      arr[j] = arr[j - 1];
      arr[j - 1] = tmp;
      j--;
    }
  }
}

template < class Etype >
void QuickSelect (Etype A[], int Low, int High, int k)
{
  int Cutoff=1;
  if (Low + Cutoff > High)
    insertionSort (&A[Low], High - Low + 1);
  else {
    
    // Sort Low, Middle, High
    int Middle = (Low + High) / 2;
    
    if (A[Middle] < A[Low]) Swap (A[Low], A[Middle]);
    if (A[High] < A[Low]) Swap (A[Low], A[High]);
    if (A[High] < A[Middle]) Swap (A[Middle], A[High]);

    // Place pivot at Position High-1
    Etype Pivot = A[Middle];
    Swap (A[Middle], A[High - 1]);
    
    // Begin partitioning
    int i, j;
    for (i = Low, j = High - 1;;)
      {
	while (A[++i] < Pivot);
	while (Pivot < A[--j]);
	if (i < j)
	  Swap (A[i], A[j]);
	else
	  break;
      }
    
    // Restore pivot
    Swap (A[i], A[High - 1]);
    
    // Recurse: only this part changes
    if (k < i) QuickSelect (A, Low, i - 1, k);
    else if (k > i) QuickSelect (A, i + 1, High, k);

  }
}

template < class Etype >
void QuickSelect (Etype A[], int N, int k)
{
  QuickSelect (A, 0, N - 1, k - 1);
}

template < class Etype >
Etype quickmedian (Etype A[], int Low, int High)
{
  int N;
  N=High-Low+1;
  Etype * B = new Etype[N];

  Etype m;
  for( int i=Low; i<=High; i++) B[i-Low]=A[i]; 
  //  for( int i=0; i<N; i++) cout << B[i] << "  "; cout << endl;
  int k;
  k=N/2;
  QuickSelect (B, 0, N - 1, k);
  m=B[k];
  //  cout << k << endl;
  delete[] B;
  return m;
}

template < class Etype >
Etype quickmedian ( Array<Etype>& A, int Low, int High)
{
  if ( High >= A.size() ) {
    cout << "High > A.size()" << endl;
    exit(0);
  }

  int N;
  N=High-Low+1;
  Etype * B = new Etype[N];

  Etype m;
  for( int i=Low; i<=High; i++) B[i-Low]=A[i]; 
  //  for( int i=0; i<N; i++) cout << B[i] << "  "; cout << endl;
  int k;
  k=N/2;
  QuickSelect (B, 0, N - 1, k);
  m=B[k];
  //  cout << k << endl;
  delete[] B;
  return m;
}

int median(Array<int>& A, int Low, int High) 
{return quickmedian(A,Low,High);}
float median(Array<float>& A, int Low, int High) 
{return quickmedian(A,Low,High);}
double median(Array<double>& A, int Low, int High) 
{return quickmedian(A,Low,High);}
/* begin median 
 *
 */


/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
template < class Etype >
double binmedian_tp ( Array<Etype>& A, int Low, int High)
{
  int i;
  int nbin=200;
  double m=-1.0;
  double xmin,xmax,bin;
  Etype Amin,Amax;

  if ( High >= A.size() ) {
    cout << "High > A.size()" << endl
	 << "binmedian_tp()" << endl;
    exit(0);
  }
  if ( Low < 0  ) {
    cout << "Low < 0 " << endl
	 << "binmedian_tp()" << endl;
    exit(0);
  }
  
  int N;
  N=High-Low+1;

  if ( N<1000000 ) {
    m=alglib::median(&A[Low],High-Low+1);
    return(m);
  }
  
  double Amean,Avar,Asd;
  Amean=mean(A,Low,High);
  Avar=variance(A,Low,High,0.0,-1);
  Asd=sqrt(Avar);
  if ( Asd/Amean < 1.0E-10 ) return(Amean);
  Amin=min(A, Low, High);
  Amax=max(A, Low, High);
  
  //  cout << Amean << "\t" << Asd << "\t" << xmin << "\t" << xmax << endl;
  
  double r=3.0;
  int C=0;
  do {
    xmin=Amean-r*Asd;
    if ( xmin < Amin ) xmin=Amin;
    xmax=Amean+r*Asd;
    if ( xmax > Amax ) xmax=Amax;
    C=0;
    for(i=Low;i<=High;++i) if (A[i]>=xmin && A[i]<=xmax) ++C;
    //cout << Amean << "\t" << Asd << "\t" << xmin << "\t" << xmax << "\t" << C << "\t" << N << "\t" << (double)C/(double)N << endl;
    if ( C < 0.95*(double)N ) r*=1.1;
  } while ( C<0.95*(double)N );
  
  do {
    xmin=Amean-r*Asd;
    if ( xmin < Amin ) xmin=Amin;
    xmax=Amean+r*Asd;
    if ( xmax > Amax ) xmax=Amax;
    C=0;
    for(i=Low;i<=High;++i) if (A[i]>=xmin && A[i]<=xmax) ++C;
    // cout << Amean << "\t" << Asd << "\t" << xmin << "\t" << xmax << "\t" << C << "\t" << N << "\t" << (double)C/(double)N << endl;
    if ( C > 0.99*(double)N ) r*=0.90;
  } while ( C>0.99*(double)N );
  
  Array<int> bincount(nbin+2);
  bincount.assign(0);
  bin=(xmax-xmin)/(double)nbin;
  for(i=Low;i<=High;++i) {
    if ( A[i]<xmin ) { bincount[0]++; continue;}
    if ( A[i]>xmax ) { bincount[nbin+1]++; continue;}
    int ibin=(A[i]-xmin)/(double)bin+1;
    bincount[ibin]++;
  }
  
  C=0;
  int ipri,ipos,npri,npos;
  ipri=ipos=nbin/2;
  npri=npos=N/2;
  for(i=0;i<=nbin+1;++i) {
    if ( bincount[i]==0 ) continue;
    C+=bincount[i];
    if ( C<N/2 ) {ipri=i; npri=C;}
    if ( C>N/2 ) {ipos=i; npos=C; break;}
  }
  
  if ( bincount[ipos] < 4 ) return(xmin+double(ipos-0.5)*bin);
  
  double sum,sum2,std;
  sum=sum2=0;C=0;
  for (i = Low; i <= High; ++i ) {
    if ( A[i]>=xmin+(double)(ipos-1)*bin && A[i]<=xmin+(double)(ipos)*bin ) {
      sum+=(double)A[i];
      sum2+=(double)A[i]*(double)A[i];
      ++C;
    }
  }
  sum/=double(C);
  sum2=sum2/double(C)-sum*sum;
  std=sqrt(sum2);
  if ( std/sum < 1.0E-10 ) return(sum);

  m=((xmin+double(ipos-1)*bin)*double(N/2-npri)
     +(xmin+double(ipos)*bin)*double(npos-N/2))/double(npos-npri);
  
  //cout << ipri << "\t" << bincount[ipri] << "\t" 
  //     << ipos << "\t" << bincount[ipos] << "\t" << m << endl;
  // exit(0);

  return m;
}
double binmedian(Array<int>& A, int Low, int High) 
{return binmedian_tp(A,Low,High);}
double binmedian(Array<float>& A, int Low, int High) 
{return binmedian_tp(A,Low,High);}
double binmedian(Array<double>& A, int Low, int High) 
{return binmedian_tp(A,Low,High);}
/* end binmedian 
 *
 */

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
template < class Etype >
double samplemedianestimator_tp ( Array<Etype>& A)
{
  double m=0.0;
  int i,j,k;

  if ( A.size()<1000000 ) {
    m=alglib::median(&A[0],A.size());
    return(m);
  }
  
  int ns=31;
  Array<Etype> RDtmp(A.size()/ns);
  Array<double> RDMAD(ns,0.0);
  for(j=0;j<ns;j++) {
    for(i=0+j,k=0; i<A.size() && k<RDtmp.size(); i+=ns,++k) RDtmp[k]=A[i];
    RDMAD[j]=alglib::median(&RDtmp[0],RDtmp.size());  
  }
  m=alglib::median(&RDMAD[0],RDMAD.size());

  return m;
}
template < class Etype >
double samplemedianestimator_tp ( Array<Etype>& A, int low, int high)
{
  if ( low<0 || high<0 || low>high || high>=A.size() ) {
    cout << "check low high parameters" << endl
	 << low << "\t" << high << "\t" << A.size() << endl
	 << "samplemedianestimator_tp()" << endl;
    exit(0);
  }
  
  double m=0.0;
  int i,j,k;
  
  if ( A.size()<1000000 ) {
    m=alglib::median(&A[0],A.size());
    return(m);
  }
  if ( high-low < 1000000 ) {
    m=alglib::median(&A[low],high-low+1);
    return(m);
  }
  
  int ns=31;
  Array<Etype> RDtmp((high-low+1)/ns);
  Array<double> RDMAD(ns,0.0);
  for(j=0;j<ns;j++) {
    for(i=low+j,k=0; i<=high && k<RDtmp.size(); i+=ns,++k) RDtmp[k]=A[i];
    RDMAD[j]=alglib::median(&RDtmp[0],RDtmp.size());  
  }
  m=alglib::median(&RDMAD[0],RDMAD.size());

  return m;
}
double samplemedianestimator(Array<int>& A) 
{ return samplemedianestimator_tp(A);}
double samplemedianestimator(Array<float>& A) 
{ return samplemedianestimator_tp(A);}
double samplemedianestimator(Array<double>& A) 
{ return samplemedianestimator_tp(A);}
double samplemedianestimator(Array<int>& A, int low, int high) 
{ return samplemedianestimator_tp(A, low, high);}
double samplemedianestimator(Array<float>& A, int low, int high) 
{ return samplemedianestimator_tp(A, low, high);}
double samplemedianestimator(Array<double>& A, int low, int high) 
{ return samplemedianestimator_tp(A, low, high);}


/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
string lastline(string file) {
  ifstream inp(file.c_str());
  if ( !inp ) {
    cout << "File " << file << " not exist" << endl; 
    exit(0);
  }
  //cout << file << endl;

  inp.seekg(0, ios::end);
  long int current;
  current=inp.tellg();
  if ( current < 0 ) {
    cout << "file end position " << current << endl;
    cout << "change type of current to longer bytes to accomodate size of the file" << endl;
    cout << "change type of indices in other places too!" << endl;
    exit(0);
  }

  /*
    cout << current << endl;
    cout << std::numeric_limits<int>::min() << endl;
    cout << std::numeric_limits<int>::max() << endl;
    cout << std::numeric_limits<long int>::min() << endl;
    cout << std::numeric_limits<long int>::max() << endl;
    exit(0);
  */

  string lastLine;
  long int i=100;
  i= current>i? i:current-1;
  inp.seekg(-i, ios::end);
  char b;
  inp.read(&b,1);
  int foundnr=0;
  while ( b != '\n' || foundnr<2 ){ 
    // inp.seekg(length-i,ios::beg);
    i++;
    if (b != '\n') foundnr++;
    inp.seekg(-i,ios::end);
    inp.read(&b,1);
    //cout << b ;
  }
  while( !inp.eof() ) {
    string tmps1;
    getline(inp,tmps1);
    // cout << "---[" << tmps1 << "]---" << endl;
    if ( tmps1.length()>2 ) lastLine=tmps1;
  }
  inp.close();

  return(lastLine);
}  

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
string firstline(string file) {
  string firstLine;
  
  ifstream inp(file.c_str());
  if ( !inp ) {
    cout << "File " << file << " not exist" << endl; 
    exit(0);
  }
  
  while ( !inp.eof() ) {
    string tmps;
    getline(inp,tmps);
    if (tmps.length() < 2) continue;
    if (tmps[0] == '#' ) continue;
    firstLine=tmps;
    break; 
  };
  inp.close();
  
  return(firstLine);
}  

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/* GNUPLOT - interpol.c */
/**/
/*[
 * Copyright 1986 - 1993, 1998, 2004   Thomas Williams, Colin Kelley
 *
 * Permission to use, copy, and distribute this software and its
 * documentation for any purpose with or without fee is hereby granted,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 *
 * Permission to modify the software is granted, but not the right to
 * distribute the complete modified source code.  Modifications are to
 * be distributed as patches to the released version.  Permission to
 * distribute binaries produced by compiling modified sources is granted,
 * provided you
 *   1. distribute the corresponding source modifications from the
 *    released version in the form of a patch file along with the binaries,
 *   2. add special version identification to distinguish your version
 *    in addition to the base release version number,
 *   3. provide your name and address as the primary contact for the
 *    support of your modified version, and
 *   4. retain our contact information in regard to use of the base
 *    software.
 * Permission to distribute the released version of the source code along
 * with corresponding source modifications in the form of a patch file is
 * granted with same provisions 2 through 4 for binary distributions.
 *
 * This software is provided "as is" without express or implied warranty
 * to the extent permitted by applicable law.
 ]*/
/**/
/**/
/*
 * C-Source file identification Header
 *
 * This file belongs to a project which is:
 *
 * done 1993 by MGR-Software, Asgard  (Lars Hanke)
 * written by Lars Hanke
 *
 * Contact me via:
 *
 *  InterNet: mgr@asgard.bo.open.de
 *      FIDO: Lars Hanke @ 2:243/4802.22   (as long as they keep addresses)
 *
 **************************************************************************
 *
 *   Project: gnuplot
 *    Module:
 *      File: interpol.c
 *
 *   Revisor: Lars Hanke
 *   Revised: 26/09/93
 *  Revision: 1.0
 *
 **************************************************************************
 *
 * LEGAL
 *  This module is part of gnuplot and distributed under whatever terms
 *  gnuplot is or will be published, unless exclusive rights are claimed.
 *
 * DESCRIPTION
 *  Supplies 2-D data interpolation and approximation routines
 *
 * IMPORTS
 *  plot.h
 *    - cp_extend()
 *    - structs: curve_points, coordval, coordinate
 *
 *  setshow.h
 *    - samples, axis array[] variables
 *    - plottypes
 *
 *  proto.h
 *    - solve_tri_diag()
 *    - typedef tri_diag
 *
 * EXPORTS
 *  gen_interp()
 *  sort_points()
 *  cp_implode()
 *
 * BUGS and TODO
 *  I would really have liked to use Gershon Elbers contouring code for
 *  all the stuff done here, but I failed. So I used my own code.
 *  If somebody is able to consolidate Gershon's code for this purpose
 *  a lot of gnuplot users would be very happy - due to memory problems.
 *
 **************************************************************************
 *
 * HISTORY
 * Changes:
 *  Nov 24, 1995  Markus Schuh (M.Schuh@meteo.uni-koeln.de):
 *      changed the algorithm for csplines
 *      added algorithm for approximation csplines
 *      copied point storage and range fix from plot2d.c
 *
 *  Dec 12, 1995 David Denholm
 *      oops - at the time this is called, stored co-ords are
 *      internal (ie maybe log of data) but min/max are in
 *      user co-ordinates.
 *      Work with min and max of internal co-ords, and
 *      check at the end whether external min and max need to
 *      be increased. (since samples_1 is typically 100 ; we
 *      dont want to take more logs than necessary)
 *      Also, need to take into account which axes are active
 *
 *  Jun 30, 1996 Jens Emmerich
 *      implemented handling of UNDEFINED points
 */

/* HBB 990205: rewrote the 'bezier' interpolation routine,
 * to prevent numerical overflow and other undesirable things happening
 * for large data files (num_data about 1000 or so), where binomial
 * coefficients would explode, and powers of 'sr' (0 < sr < 1) become
 * extremely small. Method used: compute logarithms of these
 * extremely large and small numbers, and only go back to the
 * real numbers once they've cancelled out each other, leaving
 * a reasonable-sized one. */
/**/
/*
 * cp_binomial() computes the binomial coefficients needed for BEZIER stuff
 *   and stores them into an array which is hooked to sdat.
 * (MGR 1992)
 */
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void bezier_smooth(Array<double>& RD, Array<int>& idx, Array<double>& RDt)
{
  
  if ( idx.size() > RD.size() ) {
    cout << "too many reference points" << endl;
    cout << "bezier_transfer()" << endl;
    exit(0);
  }
  
  Array<double> coeff(idx.size());
  
  int points=idx.size();
  int n = points - 1;
  int e = n / 2;
  
  coeff[0] = 0.0;
  for (int k = 0; k < e; ++k) 
    coeff[k + 1] = coeff[k] + log(((double) (n - k)) / ((double) (k + 1)));
  for (int k = n; k >= e; k--)  coeff[k] = coeff[n - k];
  
  RDt[0]=RD[idx[0]];
  RDt[RDt.size()-1]=RD[idx[idx.size()-1]];

  for (int i=1; i<RDt.size()-1; ++i) {
    double sr;
    sr=(double)i / (double)(RDt.size()-1);
    
    double ly = 0.0;
    double log_dsr_to_the_n = (double)n * log(1.0 - sr);
    double log_sr_over_dsr = log(sr) - log(1.0 - sr);
    
    for (int j=0; j < idx.size(); ++j) {
      double u = exp(coeff[j] + log_dsr_to_the_n + (double)j * log_sr_over_dsr);
      ly += RD[idx[j]] * u;
    }
    
    RDt[i]=ly;    
  }
  
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
template<class Etype>
void arrayindex_tp(Array<Etype>& val, Array<int>& rk, int order)
{
  
  if ( val.size() != rk.size() ) {
    cout << "array sizes not same" << endl
	 << "val.size()=" << val.size() << endl
	 << "rk.size()=" << rk.size() << endl
	 << "arrayindex_tp()" << endl;
    exit(0);
  }
  if ( order==0 ) {
    cout << "int order>0 : increasing order from idx 0~n" << endl
	 << "int order<0 : decreasing order from idx 0~n" << endl
	 << "int order=0 : no such option" << endl
	 << "order=" << order << endl;
    exit(0);
  }

  order = order>0 ? 1:-1;  

  int i,k;
  int itmp;
  for(i=0;i<rk.size();++i) rk[i]=i;
  
  bool swapped = false;
  k=0;
  do {
    swapped = false;
    if ( k<0 ) k=0;
    for(i=k;i>=0 && i<(int)val.size()-1;++i) {
      if ( order>0 && val[ rk[i] ] > val[ rk[i+1] ] ) swapped = true;
      if ( order<0 && val[ rk[i] ] < val[ rk[i+1] ] ) swapped = true;
      if ( swapped ) { itmp=rk[i]; rk[i]=rk[i+1]; rk[i+1]=itmp; k=i-order; break;}
    }
  } while ( swapped );
  
} //void arrayindex_tp()
void arrayindex(Array<double>& val, Array<int>& rk, int order)
{  arrayindex_tp(val,rk,order); }
void arrayindex(Array<float>& val, Array<int>& rk, int order)
{  arrayindex_tp(val,rk,order); }
void arrayindex(Array<int>& val, Array<int>& rk, int order)
{  arrayindex_tp(val,rk,order); }

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
template<class Etype>
void arrayrank_tp(Array<Etype>& val, Array<int>& rk, int order)
{
  cout << "code not checked yet!!!!!!!!!!!!!!!" << endl;  
  cout << "code not checked yet!!!!!!!!!!!!!!!" << endl;
  cout << "code not checked yet!!!!!!!!!!!!!!!" << endl; exit(0);  
  if ( val.size() != rk.size() ) {
    cout << "array sizes not same" << endl
	 << "val.size()=" << val.size() << endl
	 << "rk.size()=" << rk.size() << endl
	 << "arrayrank_tp()" << endl;
    exit(0);
  }
  if ( order==0 ) {
    cout << "int order>0 : increasing order from idx 0~n" << endl
	 << "int order<0 : decreasing order from idx 0~n" << endl
	 << "int order=0 : no such option" << endl
	 << "order=" << order << endl;
    exit(0);
  }
  
  int i;
  int itmp;
  for(i=0;i<rk.size();++i) rk[i]=i;
  
  bool swapped = false;
  do {
    swapped = false;
    for(i=0;i<(int)val.size()-1;++i) {
      if ( order>0 && val[ rk[i] ] > val[ rk[i+1] ] ) swapped = true;
      if ( order<0 && val[ rk[i] ] < val[ rk[i+1] ] ) swapped = true;
      if ( swapped ) {itmp=rk[i]; rk[i]=rk[i+1]; rk[i+1]=itmp; }
    }
  } while ( swapped );
  
  Array<int> idx(rk.size());
  for(i=0;i<rk.size();++i) idx[ rk[i] ] = i ;
  rk=idx;  
} //void arrayrank_tp()
void arrayrank(Array<double>& val, Array<int>& rk, int order)
{  arrayrank_tp(val,rk,order); }
void arrayrank(Array<float>& val, Array<int>& rk, int order)
{  arrayrank_tp(val,rk,order); }
void arrayrank(Array<int>& val, Array<int>& rk, int order)
{  arrayrank_tp(val,rk,order); }

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
