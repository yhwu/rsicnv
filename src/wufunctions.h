#ifndef _WUFUNCTIONS_H
#define _WUFUNCTIONS_H

/* check information from a pid file, or any file 
 * equivalent to grep $fields file */
string procpidstatus(string file, string fields);
/* check information from a pid file, /proc/$pid/status 
 * equivalent to grep $fields  /proc/$pid/status */
string procpidstatus(int pid, string fields);

string to_lower(string word);
string to_upper(string word);

/* case insensitive equal find */
size_t ci_find(const string& str1, const string& str2);
bool ci_equal(const string& str1, const string& str2);

/* bezier smoothing */
void bezier_smooth(Array<double>& RD, Array<int>& idx, Array<double>& RDt);

/* get last/first line of a file */
string lastline(string file);
string firstline(string file);

/* partition_median */
double partition_median(double *x, size_t n);
double partition_median(float *x, size_t n);
double partition_median(size_t *x, size_t n);
double partition_median(int *x, size_t n);

double partition_interquartilerange(double* x, size_t n);
double partition_interquartilerange(float* x, size_t n);
double partition_interquartilerange(int* x, size_t n);
double partition_interquartilerange(size_t* x, size_t n);

double partition_lowerquartile(double* x, size_t n); 
double partition_lowerquartile(float* x, size_t n); 
double partition_lowerquartile(int* x, size_t n); 
double partition_lowerquartile(size_t* x, size_t n); 

double partition_upperquartile(double* x, size_t n);
double partition_upperquartile(float* x, size_t n);
double partition_upperquartile(int* x, size_t n);
double partition_upperquartile(size_t* x, size_t n);

double partition_percentile(double *x, size_t n, double  p);
double partition_percentile(float *x, size_t n, double  p);
double partition_percentile(int *x, size_t n, double  p);
double partition_percentile(size_t *x, size_t n, double  p);


/* runmed */
void runmed(Array<double>& y, Array<double>& smo, int n, int band,
	    int end_rule, int debug);
void runmed(Array<float>& y, Array<float>& smo, int n, int band,
	    int end_rule, int debug);
void runmed(Array<int>& y, Array<int>& smo, int n, int band,
	    int end_rule, int debug);
/* runmed */

/* runmean */ 
void runmean(Array<double>& y, Array<double>& smo, int n, int band,
	     int end_rule);
void runmean(Array<float>& y, Array<float>& smo, int n, int band,
	     int end_rule);
void runmean(Array<int>& y, Array<float>& smo, int n, int band,
	     int end_rule);
/* runmean */ 

/* mean */
double mean(Array<double>& y, int low, int high);
double mean(Array<float>& y, int low, int high);
double mean(Array<int>& y, int low, int high);
/* mean */

/* min */
double min(Array<double>& y, int low, int high);
float min(Array<float>& y, int low, int high);
int min(Array<int>& y, int low, int high);
/* min */

/* max */
double max(Array<double>& y, int low, int high);
float max(Array<float>& y, int low, int high);
int max(Array<int>& y, int low, int high);
/* max */

/* variance */
double variance(Array<double>& y, int low, int high, double d, int end_rule);
double variance(Array<float>& y, int low, int high, double d, int end_rule);
double variance(Array<int>& y, int low, int high, double d, int end_rule);
/* variance */

/* remove_outlier */
void remove_outlier(Array<double>& y, Array<double>& smo, int low, int high, double d, int end_rule);
void remove_outlier(Array<float>& y, Array<float>& smo, int low, int high, double d, int end_rule);
void remove_outlier(Array<int>& y, Array<int>& smo, int low, int high, double d, int end_rule);
/* remove_outlier */

/* median */
int median(Array<int>& A, int Low, int High);
float median(Array<float>& A, int Low, int High);
double median(Array<double>& A, int Low, int High);
/* median */

double binmedian(Array<int>& A, int Low, int High);
double binmedian(Array<float>& A, int Low, int High); 
double binmedian(Array<double>& A, int Low, int High); 

double samplemedianestimator(Array<int>& A);
double samplemedianestimator(Array<float>& A);
double samplemedianestimator(Array<double>& A);
double samplemedianestimator(Array<int>& A, int low, int high);
double samplemedianestimator(Array<float>& A, int low, int high);
double samplemedianestimator(Array<double>& A, int low, int high);

void arrayindex(Array<double>& val, Array<int>& rk, int order);
void arrayindex(Array<float>& val, Array<int>& rk, int order);
void arrayindex(Array<int>& val, Array<int>& rk, int order);
void arrayrank(Array<double>& val, Array<int>& rk, int order);
void arrayrank(Array<float>& val, Array<int>& rk, int order);
void arrayrank(Array<int>& val, Array<int>& rk, int order);

#endif
