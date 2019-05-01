/*
 * Christian Gaser
 * $Id: cat_sanlm.c 1172 2017-08-31 13:19:22Z gaser $ 
 *
 */

#include "math.h"
#include "mex.h"
#include <stdlib.h>
#include "matrix.h"

extern void anlm(float* ima, int v, int f, int rician, const int* dims);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
/* Declarations */
float *ima;
int   i, v,f, ndim, rician, dims2[3];
const mwSize *dims;

/* check inputs */
if (nrhs<3)
  mexErrMsgTxt("At least 3 inputs required.");
else if (nlhs>0)
  mexErrMsgTxt("No output arguments allowed.");
  
if (!mxIsSingle(prhs[0]))
	mexErrMsgTxt("First argument must be float.");

/* get input image */
ima = (float*)mxGetPr(prhs[0]);

ndim = mxGetNumberOfDimensions(prhs[0]);
if (ndim!=3)
  mexErrMsgTxt("Images does not have 3 dimensions.");
  
dims = mxGetDimensions(prhs[0]);

/* get parameters */
v = (int)(mxGetScalar(prhs[1]));
f = (int)(mxGetScalar(prhs[2]));

if (nrhs==4)
  rician = (int)(mxGetScalar(prhs[3]));
else  rician = 0;

/* we need to convert dims to int */
for(i = 0; i < 3; i++) dims2[i] = (int)dims[i]; 

anlm(ima, v, f, rician, dims2); 

return;

}

