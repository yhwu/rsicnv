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

/**** user headers ****/
#include "pstream.h"
#include "wu2.h"
#include "rsi.h"
#include "wufunctions.h"
#include "alglibinterface.h"
//#include "linasminterface.h"
#include "gccontent.h"
/**** user headers ****/

#include "sam.h"
#include "loaddata.h"
#include "plotcnv.h"

int plot::pts=30000;
Array<int>* plot::RD=NULL;
int plot::count=0;
double plot::RDmed=0;
string plot::format="ps";
string plot::term="png xffffff x222222";
string plot::convert="convert -limit thread 1 -limit area 256MB -limit disk 512MB -density 72 -rotate 90 -background white -render -antialias -flatten ";
string plot::folder="cnv_plots";

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
std::string exec(const char* cmd) 
{
  FILE* pipe = popen(cmd, "r");
  if (!pipe) return "ERROR";
  char buffer[1024];
  std::string result = "";
  while( !feof(pipe) ) {
    if(fgets(buffer, 128, pipe) != NULL) result+=buffer;
  }
  pclose(pipe);
  return result;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
float gnuplot_version(){
  string tmps="gnuplot -V 2>/dev/null | cut -d' ' -f2 ";
  float v=-1.0;
  string tmp1=exec( tmps.c_str() );
  istringstream iss(tmp1);
  
  if ( iss >> v ) return v;
  else {
    cerr << "gnuplot not found" << endl;
    return(-1.0);
  }
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void plot_void(cnv_st& icnv, string title, string datfile, string gpfile, string imgfile)
{
  double RDmed=plot::RDmed;
  
  if ( icnv.tid != rsi::tid ) {
    cerr << "Warning, chr numbers don't match, "
	 << rsi::target_name[icnv.tid] << " | " 
	 << rsi::target_name[rsi::tid] << endl;
  }
  
  if ( plot::format=="ps" ) plot::term="postscript color enhanced solid";
  if ( plot::format=="eps" ) plot::term="postscript eps enhanced solid";
  if ( plot::format=="png" ) plot::term="png xffffff x222222";
  
  ofstream FTMP(datfile.c_str());
  if ( !FTMP ) {
    cerr << "file " << datfile << " can't be created" << "\n"
	 << "plot_void()" << endl;
    return;
  }
  
  FTMP << "#" << icnv.start << " ~ " << icnv.end << "  "
       << icnv.length << "  "
       << rsi::cnv_type[icnv.type] << "  " << icnv.p1 << endl;
  
  // CNV and reference : pos RD refrunmean //
  FTMP << icnv.start << "\t0\tNaN" << endl;
  FTMP << icnv.end << "\t0\tNaN" << endl;
  
  FTMP << endl;
  FTMP << endl;
  
  FTMP << icnv.start << "\t0\tNaN" << endl;
  FTMP << icnv.end << "\t0\tNaN" << endl;
  
  FTMP << endl;
  FTMP << endl;
  
  FTMP << icnv.start << "\t0\tNaN" << endl;
  FTMP << icnv.end << "\t0\tNaN" << endl;

  FTMP << endl;
  FTMP << endl;

  FTMP << icnv.start << "\t0\tNaN" << endl;
  FTMP << icnv.end << "\t0\tNaN" << endl;

  FTMP << endl;
  FTMP << endl;

  FTMP << icnv.start << "\t0\tNaN" << endl;
  FTMP << icnv.end << "\t0\tNaN" << endl;

  FTMP << endl;
  FTMP << endl;

  // nothing //
  FTMP << "NaN\tNaN\n"
       << "NaN\tNaN\n";
  FTMP << endl;
  FTMP << endl;
  
  FTMP.close(); // test save
  
  replace(title.begin(), title.end(), '~', '-');
  title+=" Cross border";
  
  int d=icnv.end-icnv.start+1;
  int i1=icnv.start-rsi::chklen*d;
  int i2=icnv.end+rsi::chklen*d;
  int ymax=10;
  double y2max=1.0;

  d=(i2-i1+1)/6;
  i1=i1+d/2;
  i2=i2-d/2;
  
  ofstream FGP(gpfile.c_str());
  float version=gnuplot_version();
  //version=4.0;
  string offset="offset";
  if ( version < 4.19 ) offset="";
  
  if ( version > 4.19 ) {
    FGP << "f=\"" << datfile << "\"" << endl
	<< "set datafile missing \'NaN\'" << endl
	<< "info=\"" <<  title << " \"" << endl
	<< "set terminal " << plot::term << endl
	<< "set output \"" << imgfile << "\"" << endl
	<< "#set nokey" << endl
	<< "##set label 1 info at graph  0.25, graph  0.9" << endl
	<< "set title info offset 0,-0.5" << endl
	<< "set xrange [" << i1 << ":" << i2 << "]" << endl
	<< "set xtics " << i1 << "," << d << "," << i2 << endl
	<< "set yrange [0:" << ymax << "]" << endl;
    if ( plot::format==".ps" ) {
      FGP << "set xtics font \"helvetica bold,18\"" << endl
	  << "set ytics font \"helvetica bold,18\"" << endl
	  << "set title info font \"helvetica bold,24\"" << endl;
    }
    FGP << "set y2range[0:" << y2max << "]" << endl;
    FGP << "set ytics nomirror" << endl;
    FGP << "plot \\" << endl
	<< "f in 0 u 1:2 w p pt 7 ps 0.5 lt rgb \"blue\" not, \\" << endl
	<< "f in 1 u 1:2 w p pt 7 ps 0.5 lt rgb \"red\" not, \\" << endl
	<< "f in 2 u 1:2 w p pt 7 ps 0.5 lt rgb \"blue\" not, \\" << endl
	<< "f in 0 u 1:3 w l lt 1 lw 8 lc rgb \"green\" t \"Runmean\", \\" << endl
	<< "f in 2 u 1:3 w l lt 1 lw 8 lc rgb \"green\" not, \\" << endl
	<< "f in 3 u 1:2 w l lt 0 lw 8 lc rgb \"green\" not, \\" << endl
	<< "f in 4 u 1:2 w l lt 1 lw 8 lc rgb \"cyan\" t \"CNV mean\", \\" << endl;
    if ( RDmed>0 ) 
      FGP << RDmed << " w l lt 4 lw 4 t \"CHROM med\", \\" << endl;
    FGP << "f in 5 u 1:2 not" << endl;
    FGP << "set output" << endl
	<< "quit"
	<< endl;
  } 
  else {
    FGP << "#f=\"" << datfile << "\"" << endl
	<< "#info=\"" <<  title << " \"" << endl
	<< "set terminal " << plot::term << endl
	<< "set output \"" << imgfile << "\"" << endl
	<< "#set nokey" << endl
	<< "set title \"" << title << "\"" << offset << " 0,-0.5" << endl
	<< "set xrange [" << i1 << ":" << i2 << "]" << endl
	<< "set xtics " << i1 << "," << d << "," << i2 << endl
	<< "set yrange [0:" << ymax << "]" << endl;
    if ( plot::format==".ps" ) {
      FGP << "set xtics font \"helvetica bold,18\"" << endl
	  << "set ytics font \"helvetica bold,18\"" << endl
	  << "set title \"" << title << "\" 0,-0.5 font \"helvetica bold,24\""  << endl;
    }
    FGP << "set ytics nomirror" << endl;
    FGP << "set y2range[0:" << y2max << "]" << endl;
    FGP << "plot \\" << endl
	<< "\"" << datfile << "\" in 0 u 1:2 w p pt 7 ps 0.5 lt 3 not, \\" << endl
	<< "\"" << datfile <<  "\" in 1 u 1:2 w p pt 7 ps 0.5 lt 1 not, \\" << endl
	<< "\"" << datfile << "\" in 2 u 1:2 w p pt 7 ps 0.5 lt 3 not, \\" << endl
	<< "\"" << datfile << "\" in 0 u 1:3 w l lt 2 lw 8 t \"runmean\", \\" << endl
	<< "\"" << datfile << "\" in 2 u 1:3 w l lt 2 lw 8 not, \\" << endl
	<< "\"" << datfile << "\" in 3 u 1:2 w l lt 0 lw 8 not, \\" << endl
	<< "\"" << datfile << "\" in 4 u 1:2 w l lt 5 lw 8 t \"cnvmed\", \\" << endl;
    if ( RDmed>0 ) 
      FGP << RDmed << " w l lt 4 lw 4 t \"WG mean\", \\" << endl;
    FGP << "\"" << datfile << "\" in 5 u 1:2 not" << endl;
    FGP << "set output" << endl
	<< "quit"
	<< endl;
  }      
  FGP.close();
  
  string plotcommand="gnuplot < " + gpfile;
  redi::ipstream proc(plotcommand);
  string tmps;
  while (getline(proc, tmps))  std::cerr << tmps << endl;
  proc.close();
  //  system(plotcommand.c_str());
  
  /*
  if ( imgfile.find(".ps") == imgfile.size()-3 ) {
    string convertcom="convert -version | grep Image";
    if ( exec( convertcom.c_str() ).find("ImageMagick") != string::npos ) {
      string png=imgfile.substr(0,imgfile.size()-3)+".png";
      convertcom=plot::convert + imgfile + " " + png; 
      proc.open(convertcom);
      while ( getline(proc,tmps) ) cerr << tmps << endl;
      proc.close();
      //system( convertcom.c_str() );
    }
  }
  */
  
  remove(datfile.c_str());
  remove(gpfile.c_str());
  
  return;  
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void plot_icnv(cnv_st icnv, string title, 
	       string datfile, string gpfile, string imgfile)
{
  int pts=plot::pts;
  double RDmed=plot::RDmed;
  Array<int>* RD=plot::RD;
  
  if ( icnv.start>icnv.end ) swap(icnv.start, icnv.end) ;
  
  if ( icnv.start > (*RD).size() || icnv.end > (*RD).size() ) {
    cerr << "CNV exceeds reference length" << endl;
    plot_void(icnv, title, datfile, gpfile, imgfile);
    return;
  }
  
  if ( icnv.tid != rsi::tid ) {
    cerr << "Warning, chr numbers don't match, "
	 << rsi::target_name[icnv.tid] << " | " 
	 << rsi::target_name[rsi::tid] << endl;
  }
  
  if ( plot::format=="ps" ) plot::term="postscript color enhanced solid";
  if ( plot::format=="eps" ) plot::term="postscript eps enhanced solid";
  if ( plot::format=="png" ) plot::term="png xffffff x222222";
  
  int d=icnv.end-icnv.start+1;
  if ( d < rsi::m*rsi::minmlen ) d=rsi::m*rsi::minmlen;
  int i1=icnv.start-rsi::chklen*d;
  int i2=icnv.end+rsi::chklen*d;
  if ( i1<1 ) i1=1;
  if ( i1>(*RD).size()-1 ) i1=(*RD).size()-1;
  if ( i2>(*RD).size()-1 ) i2=(*RD).size()-1;
  
  int i1c=i1-1;
  int i2c=i2-1;
  if ( i1c<0 ) i1c=0;
  if ( i1c>=(*RD).size() ) i1c=(*RD).size()-1;
  if ( i2c<0 ) i2c=0;
  if ( i2c>=(*RD).size() ) i2c=(*RD).size()-1;
  int c1=icnv.start-1;
  int c2=icnv.end-1;
  if ( c1<0 ) c1=0;
  if ( c1>=(*RD).size() ) c1=(*RD).size()-1;
  if ( c2>=(*RD).size() ) c2=(*RD).size()-1;
  
  if ( i2c-i1c != i2-i1 ) {
    cerr << "coordinate converting error\n"
	 << i1 << "\t" << i1c << "\n"
	 << i2 << "\t" << i2c << endl;
  }
  if ( c2-c1 != icnv.end-icnv.start ) {
    cerr << "coordinate converting error\n"
	 << icnv.start << "\t" << c1 << "\n"
	 << icnv.end << "\t" << c2 << endl;
  }
  
  Array<int> ref(abs(icnv.start-i1+i2-icnv.end));
  int i,k;
  double sum=0;
  for(i=i1, k=0; i<icnv.start; ++i,++k) {
    int ic=i-1;
    if ( ic>=(*RD).size() ) break;
    if ( ic<0 || ic>=(*RD).size() ) ref[k]=0;
    else ref[k]=(*RD)[ic];
    sum+=ref[k];
  }
  for(i=icnv.end+1; i<=i2; ++i,++k) {
    int ic=i-1;
    if ( ic>=(*RD).size() ) break;
    if ( ic<0 || ic>=(*RD).size() ) ref[k]=0;
    else ref[k]=(*RD)[ic];
    sum+=ref[k];
  }
  icnv.refmed=sum/double(k);
  if ( k!=ref.size() ) { 
    cerr << "size error\n"
	 << "plot_icnv()" << endl;
    return;
  }
  
  sum=0.0;
  for (i=c1; i<=c2 && i<(*RD).size(); ++i) sum+=(*RD)[i];
  icnv.cnvmed=sum/double(c2-c1+1);
  
  double icnv_ref_med,icnv_ref_lqt,icnv_ref_uqt,icnv_cnv_med,icnv_cnv_lqt,icnv_cnv_uqt;
  int n=c2-c1+1;
  if (n<0) {
    cerr << c1 << "\t" << c2 << "\n"
	 << "Cannot plot it" << endl
      	 << "plot_icnv()" << endl;
    return;
  }
  icnv_cnv_med=_median(&(*RD)[c1],n);
  icnv_cnv_lqt=_lowerquartile(&(*RD)[c1],n);
  icnv_cnv_uqt=_upperquartile(&(*RD)[c1],n);
  
  Array<float> refrunmean(ref.size());
  Array<int> refrunmed(ref.size());
  d=abs(icnv.end-icnv.start)+1;
  int band=d+((d+1)%2)*1;
  if ( band > ref.size() ) band=ref.size()/2+((ref.size()/2+1)%2)*1;
  while ( band > ref.size() ) band-=2;
  if ( band<0 ) {
    cerr << "negative band size " << band << "\t" << ref.size() << "\n"
	 << "plot_icnv()" << endl;
    return;
  }
  // runmed(ref, refrunmed, ref.size(), band, 1, 0);
  // for(i=0;i<ref.size();++i) refrunmean[i]=refrunmed[i];
  runmean(ref, refrunmean, ref.size(), band, 1);
  icnv_ref_med=_median(&(ref[0]),ref.size());
  icnv_ref_lqt=_lowerquartile(&(ref[0]),ref.size());
  icnv_ref_uqt=_upperquartile(&(ref[0]),ref.size());
  
  int ymax=icnv_ref_med*2.25;
  double y2max=(double)ymax / (double)icnv.refmed / 2.0;
  
  if ( icnv_cnv_med < icnv_ref_med ) ymax=icnv_ref_med*2;
  else { ymax= icnv_cnv_uqt + 1.5*(icnv_ref_uqt-icnv_ref_lqt); }
  
  
  ofstream FTMP(datfile.c_str());
  if ( !FTMP ) {
    cerr << "file " << datfile << " can't be created" << "\n"
	 << "plot_cnv()" << endl;
    exit(0);
  }
  
  FTMP << "#" << icnv.start << " ~ " << icnv.end << "  "
       << icnv.length << "  "
       << rsi::cnv_type[icnv.type] << "  " << icnv.p1 << endl;
  
  // dat file for gnuplot
  // data block 0
  // before CNV  : pos RD refrunmean //
  int RDmax=0;
  double step=(double)(i2-i1+1)/(double)pts;
  if ( step < 1.0 ) step=1.0;
  // before cnv
  for (double ir=(double)i1+0.00001; ir<(double)icnv.start+0.000011; ir+=step) {
    i=(int) ir;
    int ic=i-1;
    if ( ic<0 || ic>=(*RD).size() ) {
      FTMP << i << "\t0\tNaN" << endl;
      continue;
    }
    FTMP << i << "\t" << (*RD)[ic] << "\t" << refrunmean[i-i1] << endl;
    if ( ic>0 ) if ( (*RD)[ic]>RDmax ) RDmax=(*RD)[ic];
  }
  FTMP << endl;
  FTMP << endl;
  
  // data block 1
  // cnv  : pos RD NaN //
  for (double ir=(double)icnv.start+0.00001; ir<=(double)icnv.end+0.000011; ir+=step) {
    int ic=ir-1;
    if ( ic<0 || ic>=(*RD).size() ) {
      FTMP << i << "\t0\tNaN" << endl;
      continue;
    }
    FTMP << (int)ir << "\t" << (*RD)[ic] << "\t" << "NaN" << endl;
    if ( ic>0 ) if ( (*RD)[ic]>RDmax ) RDmax=(*RD)[ic];
  }
  FTMP << endl;
  FTMP << endl;
  
  // data block 2
  // after cnv  : pos RD refrunmean //
  for (double ir=(double)icnv.end+1.00001; ir<=(double)i2+0.000011; ir+=step) {
    int ic=ir-1;
    if ( ic<0 || ic>=(*RD).size() ) {
      FTMP << i << "\t0\tNaN" << endl;
      continue;
    }
    FTMP << (int)ir << "\t" << (*RD)[ic] << "\t" << refrunmean[ir-i1-d] << endl;
    if ( ic>0 ) if ( (*RD)[ic]>RDmax ) RDmax=(*RD)[ic];
  }
  FTMP << endl;
  FTMP << endl;
  
  // data block 3
  // dotted line to connect runmean //
  int imid1=icnv.start-3-i1;
  int imid2=icnv.end+1-i1-d;
  if ( imid1<0 ) imid1=0;
  if ( imid2>=refrunmean.size() ) imid2=refrunmean.size()-1;
  FTMP << icnv.start << "\t" << refrunmean[imid1] << endl
       << icnv.end << "\t" << refrunmean[imid2] << endl;
  FTMP << endl;
  FTMP << endl;
  
  // data block 4
  // line for median of CNV //
  FTMP << icnv.start << "\t" << icnv_cnv_med << endl
       << icnv.end << "\t" << icnv_cnv_med << endl;
  FTMP << endl;
  FTMP << endl;
  
  // data block 5
  // lower quartile of CNV //
  FTMP << icnv.start << "\t" << icnv_cnv_lqt << endl
       << icnv.end << "\t" << icnv_cnv_lqt << endl;
  FTMP << endl;
  FTMP << endl;
  
  // data block 6
  // upper quartile of CNV //
  FTMP << icnv.start << "\t" << icnv_cnv_uqt << endl
       << icnv.end << "\t" << icnv_cnv_uqt << endl;
  FTMP << endl;
  FTMP << endl;
  
  // data block 7
  // line for median of neighbor //
  FTMP << i1 << "\t" << icnv_ref_med << endl
       << i2 << "\t" << icnv_ref_med << endl;
  FTMP << endl;
  FTMP << endl;
  
  // data block 8
  // line for lower quartile of neighbor //
  FTMP << i1 << "\t" << icnv_ref_lqt << endl
       << i2 << "\t" << icnv_ref_lqt << endl;
  FTMP << endl;
  FTMP << endl;
  
  // data block 9
  // line for upper of neighbor //
  FTMP << i1 << "\t" << icnv_ref_uqt << endl
       << i2 << "\t" << icnv_ref_uqt << endl;
  FTMP << endl;
  FTMP << endl;
  
  // data block 10
  // nothing //
  FTMP << "NaN\tNaN\n"
       << "NaN\tNaN\n";
  FTMP << endl;
  FTMP << endl;
  
  FTMP.close(); // test save
  
  if ( plot::format==".ps" || plot::format==".eps" ) {
    for(int it=0; it<(int)title.length();++it) 
      if ( title[it] == '~' ) title[it]='-';
  }
  replace(title.begin(), title.end(), '~', '-');
  
  if ( ymax>RDmax ) ymax=RDmax+RDmax/10;
  d=(i2-i1+1)/6;
  i1=i1+d/2;
  i2=i2-d/2;
  ofstream FGP(gpfile.c_str());
  float version=gnuplot_version();
  //version=4.0;
  string offset="offset";
  if ( version < 4.19 ) offset="";

  if ( version > 4.19 ) {
    FGP << "f=\"" << datfile << "\"" << endl
	<< "set datafile missing \'NaN\'" << endl
	<< "info=\"" <<  title << " \"" << endl
	<< "set terminal " << plot::term << endl
	<< "set output \"" << imgfile << "\"" << endl
	<< "#set nokey" << endl
	<< "##set label 1 info at graph  0.25, graph  0.9" << endl
	<< "set title info offset 0,-0.5" << endl
	<< "set xrange [" << i1 << ":" << i2 << "]" << endl
	<< "set xtics " << i1 << "," << d << "," << i2 << endl
	<< "set yrange [0:" << ymax << "]" << endl;
    if ( plot::format==".ps" ) {
      FGP << "set xtics font \"helvetica bold,18\"" << endl
	  << "set ytics font \"helvetica bold,18\"" << endl
	  << "set title info font \"helvetica bold,24\"" << endl;
    }
    FGP << "set y2range[0:" << y2max << "]" << endl;
    FGP << "set ytics nomirror" << endl;
    FGP << "plot \\" << endl
	<< "f in 0 u 1:2 w p pt 7 ps 0.5 lt rgb \"blue\" t \"Neighbor\", \\" << endl
	<< "f in 1 u 1:2 w p pt 7 ps 0.5 lt rgb \"red\" t \"CNV\", \\" << endl
	<< "f in 2 u 1:2 w p pt 7 ps 0.5 lt rgb \"blue\" not, \\" << endl
	<< "f in 0 u 1:3 w l lt 1 lw 8 lc rgb \"green\" t \"Runmean\", \\" << endl
	<< "f in 2 u 1:3 w l lt 1 lw 8 lc rgb \"green\" not, \\" << endl
	<< "f in 3 u 1:2 w l lt 0 lw 8 lc rgb \"green\" not, \\" << endl
	<< "f in 4 u 1:2 w l lt 1 lw 8 lc rgb \"cyan\" t \"CNV med\", \\" << endl
	<< "f in 5 u 1:2 w l lt 0 lw 8 lc rgb \"cyan\" t \"CNV 1st,3rd quart\", \\" << endl
	<< "f in 6 u 1:2 w l lt 0 lw 8 lc rgb \"cyan\" not, \\" << endl
      //<< icnv_ref_med << " w l lt 1 lw 8 lc rgb \"green\" t \"Neighbor med\", \\" << endl
	<< icnv_ref_lqt << " w l lt 0 lw 8 lc rgb \"green\" t \"Neighbor 1st,3rd quar\", \\" << endl
	<< icnv_ref_uqt << " w l lt 0 lw 8 lc rgb \"green\" not, \\" << endl;
    if ( RDmed>0 ) 
      FGP << RDmed << " w l lt 4 lw 4 t \"CHROM med\", \\" << endl;
    FGP << "f in 10 u 1:2 not" << endl;
    FGP << "set output" << endl
	<< "quit"
	<< endl;
  } 
  else {
    FGP << "#f=\"" << datfile << "\"" << endl
	<< "#info=\"" <<  title << " \"" << endl
	<< "set terminal " << plot::term << endl
	<< "set output \"" << imgfile << "\"" << endl
	<< "#set nokey" << endl
	<< "set title \"" << title << "\"" << offset << " 0,-0.5" << endl
	<< "set xrange [" << i1 << ":" << i2 << "]" << endl
	<< "set xtics " << i1 << "," << d << "," << i2 << endl
	<< "set yrange [0:" << ymax << "]" << endl;
    if ( plot::format==".ps" ) {
      replace(title.begin(), title.end(), '~', '-');
      FGP << "set xtics font \"helvetica bold,18\"" << endl
	  << "set ytics font \"helvetica bold,18\"" << endl
	  << "set title \"" << title << "\" 0,-0.5 font \"helvetica bold,24\""  << endl;
    }
    FGP << "set ytics nomirror" << endl;
    FGP << "set y2range[0:" << y2max << "]" << endl;
    FGP << "plot \\" << endl
	<< "\"" << datfile << "\" in 0 u 1:2 w p pt 7 ps 0.5 lt 3 not, \\" << endl
	<< "\"" << datfile <<  "\" in 1 u 1:2 w p pt 7 ps 0.5 lt 1 not, \\" << endl
	<< "\"" << datfile << "\" in 2 u 1:2 w p pt 7 ps 0.5 lt 3 not, \\" << endl
	<< "\"" << datfile << "\" in 0 u 1:3 w l lt 2 lw 8 t \"runmean\", \\" << endl
	<< "\"" << datfile << "\" in 2 u 1:3 w l lt 2 lw 8 not, \\" << endl
	<< "\"" << datfile << "\" in 3 u 1:2 w l lt 2 lw 2 not, \\" << endl
	<< "\"" << datfile << "\" in 4 u 1:2 w l lt 5 lw 8 t \"CNV med\", \\" << endl
	<< "\"" << datfile << "\" in 5 u 1:2 w l lt 5 lw 4 t \"CNV 1st,3rd quart\", \\" << endl
	<< "\"" << datfile << "\" in 6 u 1:2 w l lt 5 lw 4 not, \\" << endl
      //<< icnv_ref_med << " w l lt 5 lw 2 t \"Neighbor med\", \\" << endl
	<< icnv_ref_lqt << " w l lt 2 lw 2 t \"Neighbor 1st,3rd quart\", \\" << endl
	<< icnv_ref_uqt << " w l lt 2 lw 2 not , \\" << endl;
    if ( RDmed>0 ) 
      FGP << RDmed << " w l lt 4 lw 4 t \"WG mean\", \\" << endl;
    FGP << "\"" << datfile << "\" in 10 u 1:2 not" << endl;
    FGP << "set output" << endl
	<< "quit"
	<< endl;
  }      
  FGP.close();
  
  string plotcommand="gnuplot < " + gpfile ;
  
  //  redi::ipstream proc(plotcommand, redi::pstreams::pstderr);
  redi::ipstream proc(plotcommand);
  string tmps;
  while (getline(proc, tmps))  std::cerr << tmps << endl;
  proc.close();
  //  system(plotcommand.c_str());
  
  /*
    if ( imgfile.find(".ps") == imgfile.size()-3 ) {
    string convertcom="convert -version | grep Image";
    if ( exec( convertcom.c_str() ).find("ImageMagick") != string::npos ) {
    string png=imgfile.substr(0,imgfile.size()-3)+".png";
      convertcom=plot::convert + imgfile + " " + png ; 
      proc.open(convertcom);
      while (getline(proc, tmps))  std::cerr << tmps << endl;
      proc.close();
      // system( convertcom.c_str() );
    }
  }
  */

  remove(datfile.c_str());
  remove(gpfile.c_str());
  
  return;  
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void plot_cnv(vector<cnv_st> &cnvlist, Array<int>& RD)
{
  plot::folder=rsi::plotfolder;
  string command="mkdir -p "+plot::folder;
  redi::pstream proc(command);
  proc.close();
  // system( command.c_str() );
  
  plot::RD=reinterpret_cast<Array<int>*>(&RD);
  plot::RDmed=_median(&RD[0],RD.size());
  
  vector<string> psfiles(0);

  for(int i=0; i<(int)cnvlist.size(); ++i) {
    if ( cnvlist[i].tid<0 ) {
      cerr << "cannot plot cnv " 
	   << cnvlist[i].tid << "\t" 
	   << cnvlist[i].start << "\t" 
	   << cnvlist[i].end << "\t"
	   << rsi::cnv_type[cnvlist[i].type] << endl;
      continue;
    }

    plot::count+=1;
    
    ostringstream oss;
    oss.clear();
    oss.str("");
    oss << plot::folder << "/rsi_" << rsi::target_name[cnvlist[i].tid] << "_" 
      //<< setw(digits) << setfill('0') << plot::count 	<< "_" 
	<< cnvlist[i].start << "_" 
	<< cnvlist[i].end << "_" 
	<< rsi::cnv_type[cnvlist[i].type];
    
    string datfile=oss.str()+".dat";
    string gpfile=oss.str()+".gp";
    string imgfile=oss.str()+"."+plot::format;
    
    oss.clear();
    oss.str("");
    oss << rsi::target_name[cnvlist[i].tid] << ":" 
	<< cnvlist[i].start << "-" 
	<< cnvlist[i].end << " "
	<< cnvlist[i].end-cnvlist[i].start+1 << " "
	<< rsi::cnv_type[cnvlist[i].type];
    string title=oss.str();
    cerr << "plotting: " << title << endl;
    psfiles.push_back(imgfile);
    plot_icnv(cnvlist[i], title, datfile, gpfile, imgfile);
  }
  
  string convertcom="convert -version | grep Image";
  if ( exec( convertcom.c_str() ).find("ImageMagick") != string::npos ) {
    for(int i=0; i<(int)psfiles.size(); ++i) {
      string png=psfiles[i].substr(0,psfiles[i].size()-3)+".png";
      convertcom=plot::convert + psfiles[i] + " " + png ; 
      cerr << "converting: " << psfiles[i] << endl;
      redi::pstream proc(convertcom);
      string tmps;
      while (getline(proc, tmps))  std::cerr << tmps << endl;
      proc.close();
    }
  }
  
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
void plot_cnv_from_bam(vector<cnv_st> &cnvlist, samfile_t *fp_in, bam_index_t *bamidx)
{
  if ( rsi::outfile!="rsiout.txt" ) plot::folder=rsi::outfile;
  string command="mkdir -p "+plot::folder;
  redi::pstream proc(command);
  proc.close();
  //  system( command.c_str() );
  
  Array<int> RD(1);
  plot::RD=reinterpret_cast<Array<int>*>(&RD);
  plot::RDmed=0;
  
  rsi::tid=-1;

  vector<string> psfiles(0);

  for(int i=0; i<(int)cnvlist.size(); ++i) {
    
    if ( cnvlist[i].tid<0 ) {
      cerr << "cannot plot cnv " 
	   << cnvlist[i].tid << "\t" 
	   << cnvlist[i].start << "\t" 
	   << cnvlist[i].end << "\t"
	   << rsi::cnv_type[cnvlist[i].type] << endl;
      continue;
    }
    
    int d=cnvlist[i].end-cnvlist[i].start+1;
    if ( d < rsi::m*rsi::minmlen ) d=rsi::m*rsi::minmlen;
    int ref=cnvlist[i].tid;
    int beg=cnvlist[i].start-(rsi::chklen+0.5)*d;
    if ( beg<0 ) beg=0;
    int end=cnvlist[i].end+(rsi::chklen+0.5)*d;
    rsi::tid=ref;
    load_data_from_bam(fp_in, bamidx, ref, beg, end, RD);
    plot::RD=reinterpret_cast<Array<int>*>(&RD);
    
    ostringstream oss;
    oss.clear();
    oss.str("");
    oss << plot::folder << "/rsi_" << rsi::target_name[cnvlist[i].tid] << "_" 
      //<< setw(digits) << setfill('0') << plot::count 	<< "_" 
	<< cnvlist[i].start << "_" 
	<< cnvlist[i].end << "_" 
	<< rsi::cnv_type[cnvlist[i].type];
    
    string datfile=oss.str()+".dat";
    string gpfile=oss.str()+".gp";
    string imgfile=oss.str()+"."+plot::format;
    
    oss.clear();
    oss.str("");
    oss << rsi::target_name[cnvlist[i].tid] << ":" 
	<< cnvlist[i].start << "-" 
	<< cnvlist[i].end << " "
	<< cnvlist[i].end-cnvlist[i].start+1 << " "
	<< rsi::cnv_type[cnvlist[i].type];
    string title=oss.str();
    cerr << "plotting: " << title << endl;
    psfiles.push_back(imgfile);
    plot_icnv(cnvlist[i], title, datfile, gpfile, imgfile);
  }
  
  string convertcom="convert -version | grep Image";
  if ( exec( convertcom.c_str() ).find("ImageMagick") != string::npos ) {
    for(int i=0; i<(int)psfiles.size(); ++i) {
      string png=psfiles[i].substr(0,psfiles[i].size()-3)+".png";
      convertcom=plot::convert + psfiles[i] + " " + png ; 
      cerr << "converting: " << psfiles[i] << endl;
      redi::pstream proc(convertcom);
      string tmps;
      while (getline(proc, tmps))  std::cerr << tmps << endl;
      proc.close();
    }
  }
  
  return;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
