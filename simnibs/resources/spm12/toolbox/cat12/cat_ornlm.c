/*
 * Christian Gaser
 * $Id: cat_ornlm.c 1172 2017-08-31 13:19:22Z gaser $ 
 *
 */

#include "math.h"
#include "mex.h"
#include <stdlib.h>
#include "matrix.h"

extern void ornlm(float* ima, float* fima, int v, int f, float h, const int* dims);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

/* Declarations */
float *ima, *fima;
float h;
int   i, v, f, ndim, dims2[3];
const mwSize *dims;

/* check inputs */
if (nrhs!=4)
  mexErrMsgTxt("4 inputs required.");
else if (nlhs>2)
  mexErrMsgTxt("Too many output arguments.");
  
if (!mxIsSingle(prhs[0]))
	mexErrMsgTxt("First argument must be single.");

/* get input image */
ima = (float*)mxGetPr(prhs[0]);

ndim = mxGetNumberOfDimensions(prhs[0]);
if (ndim!=3)
  mexErrMsgTxt("Images does not have 3 dimensions.");
  
dims = mxGetDimensions(prhs[0]);

/* get parameters */
v = (int)(mxGetScalar(prhs[1]));
f = (int)(mxGetScalar(prhs[2]));
h = (float)(mxGetScalar(prhs[3]));

/*Allocate memory and assign output pointer*/
plhs[0] = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxREAL);

/*Get a pointer to the data space in our newly allocated memory*/
fima = (float*)mxGetPr(plhs[0]);

/* we need to convert dims to int */
for(i = 0; i < 3; i++) dims2[i] = (int)dims[i]; 

ornlm(ima, fima, v, f, h, dims2); 

return;

}

