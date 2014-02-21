#ifndef _ALGLIBINTERFACE_H
#define _ALGLIBINTERFACE_H

namespace alglib {
  double pnorm(double x) ;
  double pnorm(float x) ;
  double pnorm(int x) ;
  
  // note this median function use double internally
  // original x array is not changed
  double median(double *x, size_t n);
  double median(float *x, size_t n);
  double median(int *x, size_t n);

  double percentile(double *x, size_t n, double p);
  double percentile(float  *x, size_t n, double p);
  double percentile(int    *x, size_t n, double p);
  double percentile(size_t *x, size_t n, double p);
  
}


#endif

