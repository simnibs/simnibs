/* laplace calculation
 * ________________________________________________________________________
 * Filter SEG within the intensity range of low and high until the changes
 * are below TH. 
 *
 * L = laplace3(SEG,low,high,TH)
 *
 * SEG  = 3d sinlge input matrix
 * low  = low boundary threshold
 * high = high boundary threshold
 * TH   = threshold to controll the number of interations
 *        maximum change of an element after interation
 *
 * ________________________________________________________________________
 * Robert Dahnke
 * Center of Neuroimaging 
 * University Jena
 *
 * $Id: cat_vol_laplace3.c 1172 2017-08-31 13:19:22Z gaser $ 
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
float abs2(float n) {	if (n<0) return -n; else return n; }

/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs<4)                                       mexErrMsgTxt("ERROR:laplace3: not enought input elements\n");
  if (nrhs>5)                                       mexErrMsgTxt("ERROR:laplace3: to many input elements\n");
  if (nlhs>1)                                       mexErrMsgTxt("ERROR:laplace3: to many output elements\n");
  if (mxIsSingle(prhs[0])==0)                       mexErrMsgTxt("ERROR:laplace3: first input must be an 3d single matrix\n");
  if (nrhs==5 && mxIsDouble(prhs[4])==0)            mexErrMsgTxt("ERROR:laplace3: 5th input (voxelsize) must be a double matrix\n");
  if (nrhs==5 && mxGetNumberOfElements(prhs[4])!=3) mexErrMsgTxt("ERROR:laplace3: 5th input (voxelsize) must have 3 Elements");
  
  /* main informations about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]);
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = (int)sL[0];
  const int     y  = (int)sL[1];
  const int     xy = x*y;
  
  /* input data */
  float*SEG = (float *)mxGetPr(prhs[0]);
  float LB  = (float) mxGetScalar(prhs[1]);
  float HB  = (float) mxGetScalar(prhs[2]);
  float TH  = (float) mxGetScalar(prhs[3]); if ( TH>=0.5 || TH<0.0001 ) mexErrMsgTxt("ERROR:laplace3: threshhold must be >0.0001 and smaller than 0.5\n");
  const mwSize sS[] = {1,3}; 
  mxArray *SS = mxCreateNumericArray(2,sS,mxDOUBLE_CLASS,mxREAL);
  double*S = mxGetPr(SS);
  if (nrhs<3) {S[0]=1; S[1]=1; S[2]=1;} else {S=mxGetPr(prhs[2]);}
  
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
  for (i=0;i<nL;i++) 
  {
    if ( SEG[i]<LB ) L1[i]=0; else {if ( SEG[i]>HB ) L1[i]=1;  else L1[i] = SEG[i]; } /*(SEG[i]-LB) / BD;} */
    L2[i] = L1[i];
    if ( (SEG[i]>LB) && (SEG[i]<HB) ) LN[i]=1; else LN[i]=0;
  }

  int u,v,w,nu,nv,nw,ni;
  float Nn, diff, maxdiffi, maxdiff=1;
  while ( maxdiff > TH ) 
  {
    maxdiffi=0;
    for (i=0;i<nL;i++)
    {
      if ( L1[i]>LB && L1[i]<HB && LN[i]==1) 
      {
        ind2sub(i,&u,&v,&w,xy,x);

        /* read neighbor values */
        L2[i]=0; Nn=0;
        for (n=0;n<sN;n++) {
          ni = i + NI[n];
          ind2sub(ni,&nu,&nv,&nw,xy,x);
          if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) )==0 && SEG[ni]>=LB && SEG[ni]<=HB) {L2[i] = L2[i] + L1[ni]; Nn++;}
        }
        L2[i] /= Nn;
        diff  = abs2( L1[i] - L2[i] ); /*printf("%f %f %f\n",L1[i],L2[i],diff); */
        if ( diff>(TH/10) ) { 
          for (n=0;n<sN;n++) {
            ni = i + NI[n];
            ind2sub(ni,&nu,&nv,&nw,xy,x);
            if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) )==0 && SEG[ni]>LB && SEG[ni]<HB) LN[ni] = 1; /* if i change his neigbors has to be recalculated */
          }
        }
        LN[i]=0;
        if ( maxdiffi<diff ) maxdiffi=diff; 
      }
    }
    maxdiff = maxdiffi;
   /* printf("%f",maxdiff); */
    /* update of L1 */
    for (i=0;i<nL;i++) 
    {
      L1[i]=L2[i];
    }
    
    /* stop++;
       if (stop == 10) continue; 
       Lt=L1; L1=L2; L2=Lt;
       update L1 */
   
 }
}


