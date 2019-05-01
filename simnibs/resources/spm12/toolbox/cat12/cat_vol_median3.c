/* Median Filter
 * _____________________________________________________________________________
 * Median Filter for a 3d single image D. Bi is used to mask voxels for the 
 * filter process, whereas Bn is used to mask voxels that are used as 
 * neighbors in the filter process. Both mask can changed by intensity 
 * threshold (Bi_low,Bi_high,Bn_low,Bn_high) for D.
 *
 *  M = cat_vol_median3(D[,Bi,Bn,sf,Bi_low,Bi_high,Bn_low,Bn_high])
 *
 *  D      (single)   3d matrix for filter process 
 *  Bi     (logical)  3d matrix that marks voxels that should be filtered
 *  Bn     (logical)  3d matrix that marks voxels that are used as neighbors 
 *  sf     (double)   threshold that is used to filter the result
 *                      sf=0: no filter
 *                      sf<0: only smaller changes
 *                      sf>0: only bigger changes
 *  Bi_low  (double)  low  threshold in D for filtering (add to Bi)
 *  Bi_high (double)  high threshold in D for filtering (add to Bi)
 *  Bn_low  (double)  low  threshold in D for neighbors (add to Bn)
 *  Bn_high (double)  high threshold in D for neighbors (add to Bn)
 *
 * Used slower quicksort for median calculation, because the faster median 
 * of the median application implementation leads to wrong results. 
 *
 * TODO: check all input elements... 
 * _____________________________________________________________________________
 * Robert Dahnke 
 * Center of Neuroimaging 
 * University Jena
 *
 * $Id: cat_vol_median3.c 1185 2017-09-13 15:02:47Z gaser $ 
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

/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs<1) mexErrMsgTxt("ERROR:cat_vol_median3: not enough input elements\n");
  if (nrhs>8) mexErrMsgTxt("ERROR:cat_vol_median3: too many input elements\n");
  if (nlhs<1) mexErrMsgTxt("ERROR:cat_vol_median3: not enough output elements\n");
  if (nlhs>1) mexErrMsgTxt("ERROR:cat_vol_median3: too many output elements\n");

  /* main informations about input data (size, dimensions, ...) */
  const mwSize *sL  = mxGetDimensions(prhs[0]);
  const int     dL  = mxGetNumberOfDimensions(prhs[0]);
  const int     nL  = (int) mxGetNumberOfElements(prhs[0]);
  
  if ( dL  != 3 || mxIsSingle(prhs[0])==0)        mexErrMsgTxt("ERROR:cat_vol_median3: first input must be a single 3d matrix\n");
  if ( nrhs>1) {
    const int     nBi = (int) mxGetNumberOfElements(prhs[1]);
    
    if ( mxGetNumberOfDimensions(prhs[1]) != 3 ) mexErrMsgTxt("ERROR:cat_vol_median3: second input must be 3d - to use a later parameter use ''true(size( input1 ))''\n");
    if ( mxIsLogical(prhs[1])==0)                mexErrMsgTxt("ERROR:cat_vol_median3: second input must be a logical 3d matrix\n");
    if ( nL != nBi)                              mexErrMsgTxt("ERROR:cat_vol_median3: second input must be a logical 3d matrix with equal size than input 1\n");
  } 
  if ( nrhs>2) {
    const int     nBn = (int) mxGetNumberOfElements(prhs[2]); 
    
    if ( mxGetNumberOfDimensions(prhs[2]) != 3 ) mexErrMsgTxt("ERROR:cat_vol_median3: third input must be 3d - to use a later parameter use ''true(size( input1 ))'\n");
    if ( mxIsLogical(prhs[2])==0)                mexErrMsgTxt("ERROR:cat_vol_median3: third input must be a logical 3d matrix\n"); 
    if ( nL != nBn)                              mexErrMsgTxt("ERROR:cat_vol_median3: third input must be a logical 3d matrix with equal size than input 1\n");
  }
  
  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  float NV[27], sf, bil, bih, bnl, bnh; 
  int i,j,k,ind,ni,x,y,z,n;
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
    if ((nrhs==1 || (nrhs>=2 && Bi[ind])) && D[ind]>=bil && D[ind]<=bih) {
      n = 0;
      /* go through all elements in a 3x3x3 box */
      for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (k=-1;k<=1;k++) {
        /* check borders */ 
        if ( ((x+i)>=0) && ((x+i)<sL[0]) && ((y+j)>=0) && ((y+j)<sL[1]) && ((z+k)>=0) && ((z+k)<sL[2])) {
          ni = index(x+i,y+j,z+k,sL);
          /* check masks and NaN or Infinities */
          if ((nrhs>=3 && Bn[ni]==0) || D[ni]<bnl || D[ni]>bnh || isnan(D[ni]) ||
              D[ni]==FLT_MAX || D[ind]==-FLT_MAX ) ni = ind;
          NV[n] = D[ni];
          n++;
        }
      }
      /* get correct n */
      /* n--; */
      
      /* sort and get the median by finding the element in the middle of the sorting */
      if (n>1) { if (n==2) {
          M[ind] = (NV[0] + NV[1]) / 2;  
        }
        else {
          sort(NV,0,n); 
          M[ind] = NV[(int)(n/2)];
        }
      }
    }
    else {
      M[ind] = D[ind];
    }
  }
  
  /* selective filter settings - only big changes (only change extremly noisy data) */
  if (sf>0) {
    for (i=0;i<nL;i++) {
      if ( (nrhs>=2 && Bi[i]) && D[i]>bil && D[i]<bih && (fabs(D[i]-M[i])<sf) ) M[i]=D[i];
    }
  }
  /* selective filter settings - only small changes */
  if (sf<0) { 
    for (i=0;i<nL;i++) {
      if ( (nrhs>=2 && Bi[i]) && D[i]>bil && D[i]<bih && (fabs(D[i]-M[i])>-sf) ) M[i]=D[i];
    }
  }
 
}


