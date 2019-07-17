/* eikonal distance calculation
 * _____________________________________________________________________________
 * Estimates eikonal distance to an object (negative values) with the speed given
 * by the positive values
 *
 *  D = cat_vol_eidist(O,F) 
 *
 *  O      (single) object
 *  F      (single) speed map
 * 
 * _____________________________________________________________________________
 * Robert Dahnke 2010_01
 * Structural Brain Mapping Group
 * University Jena
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

int sub2ind(int x,int y, int z, int s[]) {
  int i=(z-1)*s[0]*s[1] + (y-1)*s[0] + (x-1);
  if (i<0 || i>s[0]*s[1]*s[2]) i=1;
  return i;
}

float abs2(float n) { if (n<0) return -n; else return n; }


/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs<1)                                       mexErrMsgTxt("ERROR:eikonal3: not enought input elements\n");
  if (nrhs>2)                                       mexErrMsgTxt("ERROR:eikonal3: to many input elements.\n");
  if (nlhs>2)                                       mexErrMsgTxt("ERROR:eikonal3: to many output elements.\n");
  if (mxIsSingle(prhs[0])==0)                       mexErrMsgTxt("ERROR:eikonal3: first  input must be an 3d single matrix\n");
  if (nrhs==2 && mxIsDouble(prhs[1])==0)            mexErrMsgTxt("ERROR:eikonal3: fourth  input must be an double matrix\n");
  if (nrhs==2 && mxGetNumberOfElements(prhs[1])!=3) mexErrMsgTxt("ERROR:eikonal3: fourth input must have 3 Elements");
  
  /* main informations about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]); int sSEG[] = {sL[0],sL[1],sL[2]}; 
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = sL[0];
  const int     y  = sL[1];
  const int     xy = x*y;

  const int sS[] = {1,3}; 
  mxArray *SS = mxCreateNumericArray(2,sS,mxDOUBLE_CLASS,mxREAL);
  double*S = mxGetPr(SS);
  if (nrhs<2) {S[0]=1; S[1]=1; S[2]=1;} else {S=mxGetPr(prhs[1]);}
  
  float s1 = abs2((float)S[0]),s2 = abs2((float)S[1]),s3 = abs2((float)S[2]);
  const float   s12  = sqrt( s1*s1  + s2*s2); /* xy - voxel size */
  const float   s13  = sqrt( s1*s1  + s3*s3); /* xz - voxel size */
  const float   s23  = sqrt( s2*s2  + s3*s3); /* yz - voxel size */
  const float   s123 = sqrt(s12*s12 + s3*s3); /* xyz - voxel size */
  const int     nr = nrhs;
  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  const int   NI[]  = {  1, -1,  x, -x, xy,-xy, -x-1,-x+1,x-1,x+1, -xy-1,-xy+1,xy-1,xy+1, -xy-x,-xy+x,xy-x,xy+x,  -xy-x-1,-xy-x+1,-xy+x-1,-xy+x+1, xy-x-1,xy-x+1,xy+x-1,xy+x+1};  
  const float ND[]  = { s1, s1, s2, s2, s3, s3,  s12, s12,s12,s12,   s13,  s13, s13, s13,   s23,  s23, s23, s23,     s123,   s123,   s123,   s123,   s123,  s123,  s123,  s123};
  
  const int   sN  = sizeof(NI)/4;  
  float       DN[sN],DI[sN],GMTN[sN],WMDN[sN],SEGN[sN],DNm;
    
  /* main volumes - actual without memory optimation ... */
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  float *D = (float *)mxGetPr(plhs[0]);  
  float *L = (float *)mxGetPr(plhs[1]);  
    
  /* input variables */
  float*SEG  = (float *)mxGetPr(prhs[0]);
  float*F    = (float *)mxGetPr(prhs[1]);

  float TH=0.01;
  if ( TH>=0.5 || TH<0.0001 ) mexErrMsgTxt("ERROR:eikonal3: threshold must be >0.0001 and smaller than 0.5\n");
  
  int i, n; 
  /* initialisiation */
  for (i=0;i<nL;i++) {
    if ( SEG[i]<0 ) {D[i]=0; L[i]=SEG[i];} 
    else            {D[i]=FLT_MAX; L[i]=0;} /*{ if ( SEG[i]<0 ) D[i]=SEG[i]-0.5;}else */
  }
  int u,v,w,nu,nv,nw,ni,kll=0;
  float diff, maxdiffi, maxdiff=1, Dio, DNi;
  
  while ( ( maxdiff > TH ) && kll<2000 ) {
    maxdiffi=0;
    kll++;
    for (i=0;i<nL;i++) {
      if ( SEG[i]>0 && SEG[i]<FLT_MAX) {
        ind2sub(i,&u,&v,&w,xy,x);
        /* read neighbor values */
        Dio = D[i];
        
        for (n=0;n<sN;n++) {
          if ( SEG[i]!=FLT_MAX ) {
            ni = i + NI[n]; ind2sub(ni,&nu,&nv,&nw,xy,x);
            if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) )==0 ) {
              DNi = D[ni] + ND[n]*SEG[i]; 
              if ( DNi<D[i] ) {D[i]=DNi; L[i]=L[ni]; } 
            } 
          }
        }
        
        diff  = abs2( Dio - D[i] );
        if ( maxdiffi<diff ) maxdiffi=diff;  
      }
    }
    for (i=nL-1;i<0;i--) {
      if ( SEG[i]>=0 && SEG[i]<FLT_MAX) {
        ind2sub(i,&u,&v,&w,xy,x);

        /* read neighbor values */
        Dio = D[i];
        
        for (n=0;n<sN;n++) {
          if ( SEG[i]!=FLT_MAX ) {
            ni = i + NI[n]; ind2sub(ni,&nu,&nv,&nw,xy,x);
            if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) )==0 ) {
              DNi = D[ni] + ND[n]*SEG[i]; 
              if ( DNi<D[i] ) {D[i]=DNi; L[i]=L[ni]; } 
            } 
          }
        }
        
        diff  = abs2( Dio - D[i] ); 
        if ( maxdiffi<diff ) maxdiffi=diff;  
      }
    }
    maxdiff = maxdiffi;   
  }
  /*printf("%d",kll); */
}