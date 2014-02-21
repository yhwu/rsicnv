#ifndef _PLOT_H
#define _PLOT_H

class plot {
public: 
  static int pts;
  static Array<int>* RD;
  static int count;
  static double RDmed;
  static string format;
  static string term;
  static string convert;
  static string folder;
  ~plot(){};
};

float gnuplot_version();
void plot_cnv(vector<cnv_st> &cnvlist, Array<int>& RD);
void plot_cnv_from_bam(vector<cnv_st> &cnvlist, samfile_t *fp_in, bam_index_t *bamidx);

#endif
