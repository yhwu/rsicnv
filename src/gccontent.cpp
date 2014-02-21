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
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <complex>
#include <algorithm>
#include <string>
#include <vector>
using namespace std;

#include "wu2.h"
#include "rsi.h"
#include "wufunctions.h"
#include "gccontent.h"

/*
  Adjust read depths of RD in [low,low+RDA.size()-1] according to GC 
  content and return them in RDA
*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void adjustgccontent(Array<int>& RD,  Array<bool>& GC, Array<double>& GCRD,
		     int bin, double RDmean, 
		     int low, Array<int>& RDA)
{
  
  if ( RD.size() != GC.size() ) {
    rsi::dout << "size of RD not equal to size of GC array" << endl
	 << RD.size() << "\t" << GC.size() << endl
	 << "adjustgccontent()" << endl;
    exit(0);
  }
  if ( low<0 || low>=RD.size() ) {
    rsi::dout << "check parameter low " << endl
	 << low << endl
	 << "adjustgccontent()" << endl;
    exit(0);
  }
  if ( RDA.size() > RD.size() ) {
    rsi::dout << "check the size of RDA, it should be much smaller" << endl
	 << RDA.size() << endl
	 << "adjustgccontent()" << endl;
    exit(0);
  }
  if ( RDA.size() <= bin ) {
    rsi::dout << "size of RDA should be much larger than bin size" << endl
	 << RDA.size() << "\t" << bin << endl
	 << "adjustgccontent()" << endl;
    exit(0);
  }
  if ( low+RDA.size() > RD.size() ) RDA.resize(RD.size()-low);
  
  int i=0,j=0,k=0;
  
  int i1,i2;
  int nGC=0;
  for(i=low,k=0; i<low+RDA.size() && i<RD.size() && k<RDA.size(); ++i,++k) {
    i1=i-bin/2;
    if (i1<0) i1=0;
    i2=i1+bin-1;
    if (i2>=RD.size()) {i2=RD.size()-1; i1=i2-bin+1;}
    if ( i1<0 ) { cout << "i1<0\nadjustgccontent()" << endl; exit(0); }
    
    if (i==low) for(j=i1;j<=i2;++j) nGC+=GC[j];
    else if ( i1>0 && i2<RD.size()-1 ) nGC=nGC-GC[i1-1]+GC[i2];
    
    RDA[k]=RD[i]*RDmean/GCRD[ nGC ]+0.5;
  }
  
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void checkgccontent(Array<int>& RD,  Array<bool>& GC )
{
  
  if ( RD.size() != GC.size() ) {
    rsi::dout << "size of RD not equal to size of GC array" << endl
	 << RD.size() << "\t" << GC.size() << endl
	 << "checkgccontent()" << endl;
    exit(0);
  }
  
  int i;
  int bin=201;
  string flag;
  
  double RDmean=0.0;
  int C=0;
  for(i=0; i<RD.size(); ++i) if (RD[i]>0) { RDmean+=RD[i]; ++C;}
  if ( C>0 ) RDmean/=(double)C;  
  if ( rsi::debug ) cout << "RDmean=" << RDmean << endl;

  
  Array<double> GCRD(bin+1);
  Array<int> GCRDn(bin+1);
  GCRD.assign(0.0);
  GCRDn.assign(0);
  
  int i1,i2;
  int istart=0;
  int nGC=0;
  for(istart=0; istart<RD.size(); ++istart) {
    i1=istart-bin/2;
    if (i1<0) i1=0;
    i2=i1+bin-1;
    if (i2>=RD.size()) {i2=RD.size()-1; i1=i2-bin+1;}
    if ( i1<0 ) { cout << "n1<0\nadjustgccontent()" << endl; exit(0); }
    
    if ( istart==0 ) for(i=i1;i<=i2;++i) nGC+=GC[i];
    else if ( i1>0 && i2<RD.size()-1 ) nGC=nGC-GC[i1-1]+GC[i2];
    
    if ( nGC<0 || nGC>bin ) 
      cout << "error: GC count < 0 " << endl
	   << i1 << "\t" << istart << "\t" << i2 << "\t" << nGC 
	   << endl;
    GCRD[ nGC ]+=RD[istart];
    GCRDn[ nGC ]++;
  }
  for(i=0; i<=bin; ++i) {
    if (GCRDn[i]>0) GCRD[i]/=double(GCRDn[i]); 
    else GCRD[i]=RDmean;
    if ( GCRD[i] < 1 ) GCRD[i]=RDmean;
  }
  
  if ( rsi::debug) {
    string outfile=rsi::inpfile+".GCsum.txt";
    ofstream FOUT(outfile.c_str());
    for(i=0;i<=bin;++i) FOUT << i << "\t" << GCRD[i] << endl;
    FOUT.close();
    cout << "GC content summary written to " << outfile << endl;
  }
  rsi::dout << "RD mean before GC adjust = " << RDmean << endl;
  
  Array<int> RDslice(RD.size()/20);
  Array<int> RDslicepre;
  
  int ipre=0;
  istart=0;
  adjustgccontent(RD, GC, GCRD, bin, RDmean, istart, RDslice);
  RDslicepre=RDslice;
  
  while( istart < RD.size()-1 ) {
    //cout << istart << endl;
    int k;
    for(i=istart,k=0; i<istart+RDslice.size()-bin && i<RD.size(); ++i,++k) 
      RD[i]=RDslice[k];
    ipre=istart+RDslice.size()-bin;
    istart+=RDslice.size();
    if ( istart<RD.size() ) 
      adjustgccontent(RD, GC, GCRD, bin, RDmean, istart, RDslice);
    for(i=ipre; k<RDslicepre.size() && i<RD.size(); ++i,++k) RD[i]=RDslicepre[k];
    RDslicepre=RDslice;
  }
  
  RDmean=0.0;
  C=0;
  for(i=0; i<RD.size(); ++i) if (RD[i]>0) { RDmean+=RD[i]; ++C;}
  RDmean/=(double)C;  
  rsi::dout << "RD mean after GC adjust = " << RDmean << endl;
  
  //  exit(0);
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void read_ref_gc(string inpfile, Array<bool>& GC)
{
  
  ifstream INP(inpfile.c_str());
  if ( !INP ) {
    cout << "File " << inpfile << " not exist" << endl;
    exit(0);
  }
  
  if ( rsi::debug ) 
    rsi::dout << "read GC content from human reference file" << endl;
  
  long size=0;
  bool found=false;
  while ( !INP.eof() ) {
    string tmps;
    getline(INP,tmps);
    if (tmps[0] != '>' ) continue;
    
    for (unsigned i=0;i<tmps.length();++i) if ( tmps[i] == ':' ) tmps[i]=' ';
    istringstream iss(tmps);
    string tmps1;
    
    vector<string> cell;
    cell.clear();
    while ( iss >> tmps1 ) if ( tmps1.length()>0 ) cell.push_back(tmps1);
    if ( rsi::chr!=cell[0].substr(1) ) continue;
    found=true;
    size=atoi(cell[7].c_str());
    if ( size < 0 ) {
      cout << "possible integer overflow of long type size" << endl
	   << "read_ref_GC()" << endl;
      exit(0);
    }
    if ( size < 100000 ) {
      cout << "does not look like a human reference, size = " << size << endl
	   << "read_ref_GC()" << endl;
      exit(0);
    }
    if ( GC.size()+rsi::start-1 > size ) {
      cout << "size of RD larger than reference  " << size << endl
	   << GC.size() << " + " << rsi::start << " = " << GC.size()+rsi::start
	   << "\t" << size << endl
	   << "read_ref_GC() warning" << endl;
    }
    if ( GC.size()+rsi::start-1 > size + 10000 ) {
      cout << "size of RD too much larger than reference  " << size << endl
	   << GC.size() << " + " << rsi::start << " = " << GC.size()+rsi::start
	   << "\t" << size << endl
	   << "read_ref_GC() warning" << endl;
      exit(0);
    }
    
    cout << tmps << endl;
    cout << rsi::chr << " : " << size << endl;
    break;
  }
  if ( !found ) {
    cout << "chr" << rsi::chr << " not found" << endl;
    exit(0);
  }
  
  long int pos=0;
  int k=-1;
  while ( !INP.eof() ) {
    string tmps;
    getline(INP,tmps);
    if (tmps[0] == '>' ) break;
    
    istringstream iss(tmps);
    string tmps1;
    iss >> tmps1;
    for (int i=0;i<(int)tmps1.size();++i) {
      ++pos;
      if ( pos==rsi::start ) k=0;
      if ( k>=0 && k<GC.size() ) { 
	GC[k] = (tmps1[i] == 'G' || tmps1[i] == 'C' );
	++k;
      }
    }
  }
  
  if ( pos != size ) {
    cout << "something is wrong reading human reference sequence" << endl
	 << pos << "\t" << size << endl
	 << "read_human_reference()" << endl;
    exit(0);
  }

  INP.close();

  if ( rsi::debug ) 
    rsi::dout << "done; read GC content from human reference file" << endl;

} //void read_human_reference

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
