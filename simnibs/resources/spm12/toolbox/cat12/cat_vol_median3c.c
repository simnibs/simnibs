/* Median Filter
 * ________________________________________________________________________
 * Median Filter for 3d single image D. Bi is used to mask voxels for the 
 * filter process, whereas Bn is used to mask voxels that are used as 
 * neighbors in the filter process. 
 *
 *  This is a special subversion to filter label maps!
 *
 *  M = median3c(D[,Bi,Bn])
 *
 *  D      (single)   3d matrix for filter process 
 *  Bi     (logical)  3d matrix that mark voxel that should be filtered
 *  Bn     (logical)  3d matrix that mark voxel that are used as neighbors 
 *
 * ________________________________________________________________________
 * Robert Dahnke 2011_01
 * Center of Neuroimaging 
 * University Jena
 *
 * $Id: cat_vol_median3c.c 764 2015-11-17 13:11:53Z gaser $ 
 */

#include "mex.h"   
#include "matrix.h"
#include "math.h"
#include "float.h"

#ifndef isnan
#define isnan(a) ((a)!=(a)) 
#endif

#define index(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))

/* qicksort */
void swap(float *a, float *b)
{
  float t=*a; *a=*b; *b=t;
}

void sort(float arr[], int beg, int end)
{
  if (end > beg + 1)
  {
    float piv = arr[beg];
    int l = beg + 1, r = end;
    while (l < r)
    {
      if (arr[l] <= piv)
        l++;
      else
        swap(&arr[l], &arr[--r]);
    }
    swap(&arr[--l], &arr[beg]);
    sort(arr, beg, l);
    sort(arr, r, end);
  }
}

float abs2(float n) {	if (n<0) return -n; else return n; }        

/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs<1) mexErrMsgTxt("ERROR:median3: not enough input elements\n");
  if (nrhs>8) mexErrMsgTxt("ERROR:median3: too many input elements\n");
  if (nlhs<1) mexErrMsgTxt("ERROR:median3: not enough output elements\n");
  if (nlhs>1) mexErrMsgTxt("ERROR:median3: too many output elements\n");

  /* main informations about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]);
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = (int) mxGetNumberOfElements(prhs[0]);

  if ( dL != 3 || mxIsSingle(prhs[0])==0)        mexErrMsgTxt("ERROR:median3: first input must be a single 3d matrix\n");
  if ( nrhs>1) {
    if ( mxGetNumberOfDimensions(prhs[1]) != 3 ) mexErrMsgTxt("ERROR:median3: second input must be 3d - to use a later parameter use ''true(size( input1 ))''\n");
    if (mxIsLogical(prhs[1])==0)                 mexErrMsgTxt("ERROR:median3: second input must be a logical 3d matrix\n");
  } else if ( nrhs>2) {
    if ( mxGetNumberOfDimensions(prhs[2]) != 3 ) mexErrMsgTxt("ERROR:median3: third input must be 3d - to use a later parameter use ''true(size( input1 ))'\n");
    if ( mxIsLogical(prhs[2])==0)                mexErrMsgTxt("ERROR:median3: third input must be a logical 3d matrix\n"); 
  }

  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  float sf, bil, bih, bnl, bnh;
  int NV[256],i,j,k,ind,ni,x,y,z,n;
  bool *Bi, *Bn;
        
  /* in- and output */
  float *D = (float *) mxGetPr(prhs[0]);
  if (nrhs>1)  Bi = (bool *) mxGetPr(prhs[1]); 
  if (nrhs>2)  Bn =  (bool *) mxGetPr(prhs[2]); 
  if (nrhs<4)  sf = 0; 
  else         sf = (float) *mxGetPr(prhs[3]);

  if (nrhs<5) bil = -FLT_MAX;   
  else        bil = (float) *mxGetPr(prhs[4]);
  if (nrhs<6) bih =  FLT_MAX;   
  else        bih = (float) *mxGetPr(prhs[5]);
  if (nrhs<7) bnl = -FLT_MAX;   
  else        bnl = (float) *mxGetPr(prhs[6]);  
  if (nrhs<8) bnh =  FLT_MAX;   
  else        bnh = (float) *mxGetPr(prhs[7]);

  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  float *M = (float *) mxGetPr(plhs[0]);
  
  /* filter process */
  for (z=0;z<sL[2];z++) for (y=0;y<sL[1];y++) for (x=0;x<sL[0];x++) {
    ind = index(x,y,z,sL);
    if ((nrhs==1 || (nrhs>=2 && Bi[ind])) && D[ind]>=0 && D[ind]>=bil && D[ind]<=bih) {
      n = 0;
      /* go through all elements in a 3x3x3 box */
      for (i=0;i<255;i++) NV[i]=0;  
      for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (k=-1;k<=1;k++) {
        /* check borders */ 
        if ( ((x+i)>=0) && ((x+i)<sL[0]) && ((y+j)>=0) && ((y+j)<sL[1]) && ((z+k)>=0) && ((z+k)<sL[2])) {
          ni = index(x+i,y+j,z+k,sL);
          /* check masks and NaN or Infinities */
          if ((nrhs>=3 && Bn[ni]==0) || D[ni]<0 || D[ni]<bnl || D[ni]>bnh || isnan(D[ni]) ||
              D[ni]==FLT_MAX || D[ind]==-FLT_MAX ) ni = ind;
          NV[(int) D[ni]]++;
        }
      }
      /* find maximum occurent value */
      M[ind]=0; n=0; for (i=0;i<255;i++) {if ((NV[i]>=0) && (NV[i]>n)) {n=NV[i]; M[ind]=i;}}
     }
    else {
      M[ind] = D[ind];
    }
  }
  
  /* selective filter settings - only big changes (only change extremly noisy data) */
  if (sf>0) {
    for (i=0;i<nL;i++) {
      if ( (nrhs>=2 && Bi[i]) && D[i]>bil && D[i]<bih && (abs2(D[i]-M[i])<sf) ) M[i]=D[i];
    }
  }
  /* selective filter settings - only small changes */
  if (sf<0) { 
    for (i=0;i<nL;i++) {
      if ( (nrhs>=2 && Bi[i]) && D[i]>bil && D[i]<bih && (abs2(D[i]-M[i])>-sf) ) M[i]=D[i];
    }
  }
 
}


