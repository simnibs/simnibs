/* laplace calculation
 * _____________________________________________________________________
 * L = laplace3(SEG,low,high,TH)
 *
 * SEG  = 3d sinlge input matrix
 * low  = low boundary threshold
 * high = high boundary threshold
 * TH   = threshold to controll the number of interations
 *        maximum change of an element after interation
 *
 * TODO: change of L1 and L2 by pointer
 * _____________________________________________________________________
 * Robert Dahnke 
 * Structural Brain Mapping Group
 * University Jena
 *
 * $Id: cat_vol_laplace3R.c 1172 2017-08-31 13:19:22Z gaser $ 
 */

#include "mex.h"   
#include "matrix.h"
#include "math.h"
#include "float.h"


/* estimate x,y,z position of index i in an array size sx,sxy=sx*sy... */
void ind2sub(int i,int *x,int *y, int *z, int sxy, int sy) {
  *z = (int)floor( i / (double)sxy ) +1; 
   i = i % (sxy);
  *y = (int)floor( i / (double)sy ) +1;        
  *x = i % sy + 1;
}
float abs2(float n) {	if (n<0) return -n; else return n; }

/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs<2)                                       mexErrMsgTxt("ERROR:laplace3R: not enought input elements\n");
  if (nrhs>4)                                       mexErrMsgTxt("ERROR:laplace3R: to many input elements\n");
  if (nlhs>1)                                       mexErrMsgTxt("ERROR:laplace3R: to many output elements\n");
  if (mxIsSingle(prhs[0])==0)                       mexErrMsgTxt("ERROR:laplace3R: 1st input must be an 3d single matrix\n");
  if (mxIsLogical(prhs[1])==0)                      mexErrMsgTxt("ERROR:laplace3R: 2nd input must be an 3d logical matrix\n");
  if (nrhs==3 && mxIsDouble(prhs[2])==0)            mexErrMsgTxt("ERROR:laplace3R: 3rd input must be an double matrix\n");
  if (nrhs==4 && mxIsDouble(prhs[3])==0)            mexErrMsgTxt("ERROR:laplace3R: 4th input (voxelsize) must be a double matrix\n");
  if (nrhs==4 && mxGetNumberOfElements(prhs[3])!=3) mexErrMsgTxt("ERROR:laplace3R: 4th input (voxelsize) must have 3 Elements");
  
  /* main informations about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]);
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = (int)sL[0];
  const int     y  = (int)sL[1];
  const int     xy = x*y;
  
  /* input data */
  float*SEG = (float *)mxGetPr(prhs[0]);
  bool *M   = (bool  *)mxGetPr(prhs[1]); 
  float TH  = (float) mxGetScalar(prhs[2]); if ( TH>=0.5 || TH<0.000001 ) mexErrMsgTxt("ERROR:laplace3R: threshhold must be >0.000001 and smaller than 0.5\n");
  const mwSize sS[] = {1,3}; 
  mxArray *SS = mxCreateNumericArray(2,sS,mxDOUBLE_CLASS,mxREAL);
  double*S = mxGetPr(SS);
  if (nrhs<4) {S[0]=1; S[1]=1; S[2]=1;} else {S=mxGetPr(prhs[3]);}
  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  const int   NI[]  = { -1, 1, -x, x, -xy, xy};  
  const int   sN = sizeof(NI)/4;    
  int i, n;
  
  /* output data */
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[2] = mxCreateLogicalArray(dL,sL);
  
  float  *L1 = (float *)mxGetPr(plhs[0]);
  float  *L2 = (float *)mxGetPr(plhs[1]);
  bool   *LN = (bool  *)mxGetPr(plhs[2]);
  
  /* intitialisiation */
  for (i=0;i<nL;i++) {
    if ( mxIsNaN(SEG[i]) ) { L1[i] = FLT_MAX; } else { L1[i] = SEG[i]; }
    L2[i] = L1[i];
    LN[i] = M[i];
  }

  int u,v,w,nu,nv,nw,ni,iter=0,maxiter=2000;
  float Nn, diff, maxdiffi, maxdiff=1;
  while ( maxdiff > TH && iter < maxiter) {
    maxdiffi=0; iter++;
    for (i=0;i<nL;i++) {
      if ( M[i] && LN[i] ) {  
        ind2sub(i,&u,&v,&w,xy,x);

        /* read neighbor values */
        L2[i]=0; Nn=0;
        for (n=0;n<sN;n++) {
          ni = i + NI[n]; 
          ind2sub(ni,&nu,&nv,&nw,xy,x);
          if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (L1[ni]==-FLT_MAX) || (L1[ni]==FLT_MAX) )==0) {L2[i] = L2[i] + L1[ni]; Nn++;}
        }
        if (Nn>0) {L2[i]/=Nn;} else {L2[i]=L1[i];}
        
        diff  = abs2( L1[i] - L2[i] ); /*printf("%f %f %f\n",L1[i],L2[i],diff); */
        if ( diff>(TH/10) ) { 
          for (n=0;n<sN;n++) {
            ni = i + NI[n]; ind2sub(ni,&nu,&nv,&nw,xy,x);
            if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (L1[ni]==-FLT_MAX) || (L1[ni]==FLT_MAX) )==0) LN[ni] = 1; /* if i change his neigbors has to be recalculated */
          }
        }
      
        LN[i]=0;
        if ( maxdiffi<diff ) maxdiffi=diff; 
      }
    }
    maxdiff = maxdiffi;

    /* update of L1 */
    for (i=0;i<nL;i++) { L1[i]=L2[i]; }
  }
  /* printf("%d\n",iter); */
}


