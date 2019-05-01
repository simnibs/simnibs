/*
 * Christian Gaser
 * $Id: cat_amap.c 1172 2017-08-31 13:19:22Z gaser $ 
 *
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "stdio.h"
#include "Amap.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  unsigned char *label, *prob, *mask;
  double *src, *mean, *voxelsize;
  double max_vol = -1e15, weight_MRF, bias_fwhm, offset;
  const mwSize *dims;
  mwSize dims3[4];
  int dims2[4];
  int i, n_classes, pve, nvox, iters_icm;
  int niters, iters_nu, sub, init, thresh, thresh_kmeans_int;
    
  if (nrhs!=11)
    mexErrMsgTxt("11 inputs required.");
  else if (nlhs>2)
    mexErrMsgTxt("Too many output arguments.");
  
  if (!mxIsDouble(prhs[0]))
    mexErrMsgTxt("First argument must be double.");
  if (!mxIsUint8(prhs[1]))
    mexErrMsgTxt("Second argument must be uint8.");
  if (!mxIsDouble(prhs[2]))
    mexErrMsgTxt("Third argument must be double.");
  if (!mxIsDouble(prhs[3]))
    mexErrMsgTxt("4th argument must be double.");
  if (!mxIsDouble(prhs[4]))
    mexErrMsgTxt("5th argument must be double.");
  if (!mxIsDouble(prhs[5]))
    mexErrMsgTxt("6th argument must be double.");
  if (!mxIsDouble(prhs[6]))
    mexErrMsgTxt("7th argument must be double.");
  if (!mxIsDouble(prhs[7]))
    mexErrMsgTxt("8th argument must be double.");
  if (!mxIsDouble(prhs[8]))
    mexErrMsgTxt("9th argument must be double.");
  if (!mxIsDouble(prhs[9]))
    mexErrMsgTxt("10th argument must be double.");
  if (!mxIsDouble(prhs[10]))
    mexErrMsgTxt("11th argument must be double.");


  
  src        = (double*)mxGetPr(prhs[0]);
  label      = (unsigned char*)mxGetPr(prhs[1]);
  n_classes  = (int)mxGetScalar(prhs[2]);
  niters     = (int)mxGetScalar(prhs[3]);
  sub        = (int)mxGetScalar(prhs[4]);
  pve        = (int)mxGetScalar(prhs[5]);
  init       = (int)mxGetScalar(prhs[6]);
  weight_MRF = (double)mxGetScalar(prhs[7]);
  voxelsize  = (double*)mxGetPr(prhs[8]);
  iters_icm  = (int)mxGetScalar(prhs[9]);
  bias_fwhm  = (double)mxGetScalar(prhs[10]);

  if ( mxGetM(prhs[8])*mxGetN(prhs[8]) != 3) 
    mexErrMsgTxt("Voxelsize should have 3 values.");

  dims = mxGetDimensions(prhs[0]);
  dims2[0] = (int)dims[0]; dims2[1] = (int)dims[1]; dims2[2] = (int)dims[2];
  dims2[3] = n_classes;
  
  /* for PVE we need more classes */
  if(pve == 6) dims2[3] += 3;
  if(pve == 5) dims2[3] += 2;

  /* mxCreateNumericArray expects mwSize data type */
  for(i = 0; i < 4; i++) dims3[i] = (mwSize)dims2[i]; 

  plhs[0] = mxCreateNumericArray(4, dims3, mxUINT8_CLASS, mxREAL);
  plhs[1] = mxCreateNumericMatrix(1, n_classes+3, mxDOUBLE_CLASS, mxREAL);
  prob  = (unsigned char *)mxGetPr(plhs[0]);
  mean  = (double *)mxGetPr(plhs[1]);

  nvox = dims[0]*dims[1]*dims[2];
  
  for(i = 0; i < nvox; i++) {
    max_vol = MAX(src[i], max_vol);
  }

  offset = 0.2*max_vol;
  
  /* add offset to ensure that CSF values are much larger than background noise */
  for (i=0; i<nvox; i++) {
    if (label[i] > 0) src[i] += offset;
  }

  /* initial labeling using Kmeans */
  if (init) {
    mask = (unsigned char *)mxMalloc(sizeof(unsigned char)*nvox);
    if(mask == NULL) {
      mexErrMsgTxt("Memory allocation error\n");
      exit(EXIT_FAILURE);
    }
    for (i=0; i<nvox; i++)
      mask[i] = (src[i]>0) ? 255 : 0;

    thresh = 0;
    thresh_kmeans_int = 128;
    iters_nu = 40;

    /* initial Kmeans estimation with 6 classes */
    max_vol = Kmeans( src, label, mask, 25, n_classes, voxelsize, dims2, thresh, thresh_kmeans_int, iters_nu, KMEANS, bias_fwhm);
    /* final Kmeans estimation with 3 classes */
    max_vol = Kmeans( src, label, mask, 25, n_classes, voxelsize, dims2, thresh, thresh_kmeans_int, iters_nu, NOPVE, bias_fwhm);

    mxFree(mask);
  }
    
  Amap(src, label, prob, mean, n_classes, niters, sub, dims2, pve, weight_MRF, voxelsize, iters_icm, offset);
  if(pve==6) Pve6(src, prob, label, mean, dims2);
  if(pve==5) Pve5(src, prob, label, mean, dims2);

}

