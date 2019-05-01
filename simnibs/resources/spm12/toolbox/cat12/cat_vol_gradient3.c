/* gradient calculation
 * _____________________________________________________________________
 * Estimation of gradient of a volume L (within a ROI M). 
 * 
 * [gi,gj,gk] = cat_vol_gradient3(L[,M])
 *
 * L          = 3d single input matrix
 * M          = 3d logical input matrix
 * [gi,gj,gk] = 3d single output matrix
 *
 * _____________________________________________________________________
 * Robert Dahnke 2014_06
 * Center of Neuroimaging 
 * University Jena
 *
 * $Id: cat_vol_gradient3.c 1185 2017-09-13 15:02:47Z gaser $ 
 */

#include "mex.h"   
#include "matrix.h"
#include "math.h"


/* estimate x,y,z position of index i in an array size sx,sxy=sx*sy... */
void ind2sub(int i,int *x,int *y, int *z, int sxy, int sy) {
  *z = (int)floor( i / (double)sxy ) +1; 
   i = i % (sxy);
  *y = (int)floor( i / (double)sy ) +1;        
  *x = i % sy + 1;
}


/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs<1)                  mexErrMsgTxt("ERROR:cat_vol_gradient3: not enought input elements\n");
  if (nrhs>2)                  mexErrMsgTxt("ERROR:cat_vol_gradient3: to many input elements\n");
  if (nlhs<3)                  mexErrMsgTxt("ERROR:cat_vol_gradient3: to less output elements\n");
  if (nlhs>3)                  mexErrMsgTxt("ERROR:cat_vol_gradient3: to many output elements\n");
  if (mxIsSingle(prhs[0])==0)  mexErrMsgTxt("ERROR:cat_vol_gradient3: input must be an 3d single matrix\n");
 
  
  /* main informations about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]);
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = (int) sL[0];
  const int     y  = (int) sL[1];
  const int     xy = x*y;

  if ( dL != 3 ) mexErrMsgTxt("ERROR:cat_vol_gradient3: input must be 3d\n");

  
  /* in- and output  */
  float*I = (float *)mxGetPr(prhs[0]);
  bool *M;
  if ( nrhs>1) {
    const int     nM = (int) mxGetNumberOfElements(prhs[1]);
    
    if ( mxGetNumberOfDimensions(prhs[1]) != 3 ) 
      mexErrMsgTxt("ERROR:cat_vol_gradient3: second input must be 3d - to use a later parameter use ''true(size( input1 ))''\n");
    if ( mxIsLogical(prhs[1])==0)                
      mexErrMsgTxt("ERROR:cat_vol_gradient3: second input must be a logical 3d matrix\n");
    if ( nL != nM)                               
      mexErrMsgTxt("ERROR:cat_vol_gradient3: second input must be a logical 3d matrix with equal size than input 1\n");
  
    M = (bool *)mxGetPr(prhs[1]);
  } 
  
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[2] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);

  float *G1 = (float *)mxGetPr(plhs[0]);
  float *G2 = (float *)mxGetPr(plhs[1]);
  float *G3 = (float *)mxGetPr(plhs[2]);
  
  int i,u,v,w,nu,nv,nw,n1i,n2i; 
  for (i=0;i<nL;i++) 
  {
    ind2sub(i,&u,&v,&w,xy,x);

    n1i=i-1; ind2sub(n1i,&nu,&nv,&nw,xy,x); if ( (n1i<0) || (n1i>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) n1i=i;
    n2i=i+1; ind2sub(n2i,&nu,&nv,&nw,xy,x); if ( (n2i<0) || (n2i>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) n2i=i;
    if ( nrhs==1 ) {
      G1[i] = ( I[n2i] - I[n1i] ) / 2;
    }
    else {
      if ( M[n1i] && M[n2i] ) 
        G1[i] = ( I[n2i] - I[n1i] ) / 2;
      else
        G1[i] = 0;  
    }
    
    n1i=i-x; ind2sub(n1i,&nu,&nv,&nw,xy,x); if ( (n1i<0) || (n1i>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) n1i=i;
    n2i=i+x; ind2sub(n2i,&nu,&nv,&nw,xy,x); if ( (n2i<0) || (n2i>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) n2i=i;
    if ( nrhs==1 )
      G2[i] = ( I[n2i] - I[n1i] ) / 2;
    else {
      if ( M[n1i] && M[n2i] )
        G2[i] = ( I[n2i] - I[n1i] ) / 2;
      else
        G2[i] = 0;  
    }
    
    n1i=i-xy; ind2sub(n1i,&nu,&nv,&nw,xy,x); if ( (n1i<0) || (n1i>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) n1i=i;
    n2i=i+xy; ind2sub(n2i,&nu,&nv,&nw,xy,x); if ( (n2i<0) || (n2i>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) n2i=i;
    if ( nrhs==1 )
      G3[i] = ( I[n2i] - I[n1i] ) / 2;
    else {
      if ( M[n1i] && M[n2i] )
        G3[i] = ( I[n2i] - I[n1i] ) / 2;
      else
        G3[i] = 0;
    }
  }

}


