#ifndef _RSI_H
#define _RSI_H

#define TYPE_DEL 0
#define TYPE_DUP 1
#define TYPE_UNKNOWN 2

struct cnv_st {
  int tid;
  int type;
  int geno;
  int status;
  int start;
  int end;
  int length;
  int sc1;
  int sc2;
  int pair;
  double score;
  double p1;
  double p2;
  double cnvmed;
  double cnvsd;
  double cnviqr;
  double refmed;
  double refsd;
  double refiqr;
  double q0;
  int rp;
cnv_st():  tid(-1),
           type(TYPE_UNKNOWN),
           geno(0),
	   status(0),
	   start(0),
	   end(0),
	   length(0),
	   sc1(0),
	   sc2(0),
	   pair(0),
	   score(0.0),
	   p1(1.0),
	   p2(1.0),
	   cnvmed(0),
	   cnvsd(0),
	   cnviqr(0),
	   refmed(0),
           refsd(0),
           refiqr(0),
           q0(-1.0),
           rp(-1) {};
} ;


class rsi {
public: 
  static bool debug;
  static bool saverd;
  static bool histonly;
  static bool overlaponly;
  static bool combine;
  static bool nocodeonly;
  static bool regiononly;
  static bool plot;
  static bool genotype;
  static bool merge;
  static bool gcadjust;
  static int pid;
  static long seed;
  static int start;
  static int end;
  static int m;
  static int L;
  static int Lmax;
  static int minq;
  static int min_baseQ;
  static string trans;
  static double minmlen;
  static double chklen;
  static double buffer;
  static int maxchkbp;
  static double cap;
  static double p;
  static double ci;
  static double RDmedian;
  static double RDsd;
  static double RDsigma;
  static double median;
  static double medmedian;
  static double nbnmedian;
  static double sigma;
  static double lamda;
  static double medlamda;
  static double nbnlamda;
  static double threshold;
  static double epsilon;
  static double factor;
  static int hetnorm;
  static ofstream fout;
  static string inpfile;
  static string outfile;
  static string bamfile;
  static string rdfile;
  static string reffile;
  static string cnvfile;
  static string logfile;
  static string function;
  static string plotfolder;
  static dual_stream dout;
  static string pidfile;
  static int bam_ref;
  static int bam_l_qseq;
  static double bam_rd;
  static double bam_rd_sd;
  static double bam_drd;
  static double bam_drd_sd;
  static string cnv_type[];
  static string chr;
  static int tid;
  static vector<string> target_name;
  static vector<cnv_st> noncodelist;
  ~rsi(){};
};


#endif
