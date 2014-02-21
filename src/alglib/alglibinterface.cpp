#include "stdafx.h"
#include<stdlib.h>
#include<stdio.h>
#include "specialfunctions.h"
#include "statistics.h"

#include "alglibinterface.h"

namespace alglib {
  double pnorm(double x){ return(alglib::normaldistribution(x)); }
  double pnorm(float x){ return(alglib::normaldistribution(double(x)));}
  double pnorm(int x){ return(alglib::normaldistribution(double(x))); }
  
  // note this median function use double internally
  // original x array is not changed
  template<class Etype>
  double samplemedian_tp(Etype *x,  ae_int_t n)
  {
    double med;
    alglib::ae_int_t nn;
    nn=n;
    alglib::real_1d_array xx;
    xx.setlength(n);
    for(ae_int_t i=0;i<n;++i) xx[i]=x[i];
    alglib::samplemedian(xx,med);
    return(med);
  }
  double median(double *x, size_t n) { return samplemedian_tp(x, n);  }
  double median(float *x, size_t n) { return samplemedian_tp(x, n);  }
  double median(int *x, size_t n) { return samplemedian_tp(x, n);  }
  
  template<class Etype>
  double samplepercentile_tp(Etype *x,  ae_int_t n, double p)
  {
    double v;
    alglib::real_1d_array xx;
    xx.setlength(n);
    for(ae_int_t i=0;i<n;++i) xx[i]=x[i];
    alglib::samplepercentile(xx, p, v);
    return(v);
  }
  double percentile(double *x, size_t n, double p) 
  { return samplepercentile_tp(x, n, p);  }
  double percentile(float *x, size_t n, double p) 
  { return samplepercentile_tp(x, n, p);  }
  double percentile(int *x, size_t n, double p) 
  { return samplepercentile_tp(x, n, p);  }
  double percentile(size_t *x, size_t n, double p) 
  { return samplepercentile_tp(x, n, p);  }
  
  
}
