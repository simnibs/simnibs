/*
 * Christian Gaser
 * $Id: Amap.h 404 2011-04-11 10:03:40Z gaser $ 
 *
 */

#ifdef MATLAB_MEX
#include "mex.h"
#endif

#define SQRT2PI 2.506628
#define G 6

#define MAX_NC 6
#define TH_COLOR 1
#define TH_CHANGE 0.001
#ifndef TINY
#define TINY 1e-15 
#endif
#ifndef HUGE
#define HUGE 1e15 
#endif
#ifndef NULL
#define NULL ((void *) 0)
#endif

#define NOPVE 0
#define KMEANS 1

#define BKGCSFLABEL 0
#define CSFLABEL    1
#define GMCSFLABEL  2
#define GMLABEL     3
#define WMGMLABEL   4
#define WMLABEL     5

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

#ifndef ROUND
#define ROUND( x ) ((long) ((x) + ( ((x) >= 0) ? 0.5 : (-0.5) ) ))
#endif

#ifndef MIN3
#define MIN3(a,b,c) (MIN(a,MIN(b,c)))
#endif

extern double Kmeans(double *src, unsigned char *label, unsigned char *mask, int NI, int n_clusters, double *voxelsize, int *dims, int thresh_mask, int thresh_kmeans, int iters_nu, int pve, double bias_fwhm);
extern void Amap(double *src, unsigned char *label, unsigned char *prob, double *mean, int nc, int niters, int sub, int *dims, int pve, double weight_MRF, double *voxelsize, int niters_ICM, double offset);
extern void Pve5(double *src, unsigned char *prob, unsigned char *label, double *mean, int *dims);
extern void Pve6(double *src, unsigned char *prob, unsigned char *label, double *mean, int *dims);
extern void MrfPrior(unsigned char *label, int nc, double *alpha, double *beta, int init, int *dims);
#ifdef SPLINESMOOTH
  extern int splineSmooth( double *src, double lambda, double distance, int subsample, double *separations, int *dims);
#endif
void smooth_subsample_double(double *vol, int dims[3], double separations[3], double s[3], int use_mask, int samp);

struct point {
  double mean;
  double var;
};

struct ipoint {
  int n;
  double s;
  double ss;
};
