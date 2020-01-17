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


/* Modified by Guilherme Saturnino, 2019
 * Removed parts of the code refering to MATLAB
*/



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
void vol_median_3(float *M, int nrhs, float *D, int *sL, unsigned char *Bi, unsigned char *Bn, float sf, float bil, float bih, float bnl, float bnh)
{
  float NV[27]; 
  int i,j,k,ind,ni,x,y,z,n;
  int nL = sL[0]*sL[1]*sL[2];
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


