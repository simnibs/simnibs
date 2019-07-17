/* intensity limited region growing 
 * _____________________________________________________________________________
 * Region growing of integer objects in O depending on a distance map that is 
 * created on the path distance to the object and the intensity in L. The 
 * region growing for the intensity map is limited by the lim where a neighbor
 * value have fit for the following equation:
 *   L(neigbor(x)) + limit <= L(x)
 * But not only the intensity of L is important, also the object distance can be 
 * used. You can set up the relation of both with the dd value, where dd(1) is
 * used for the distance component and dd(2) for the intensity component.
 * This regions growing was orignialy used for skull-stripping and to alL used 
 *
 *
 * [D,I] = cat_vol_downcut(O,L,lim,vx,dd)
 * 
 * O      (3d single)   initial object (integer values for different objects)
 * L      (3d single)   intensity image 
 * lim    (1x1 double)  limit for the neighbor intensity in the regions growing
 *                        L(neigbor(x)) + limit <= L(x)
 * vx     (1x3 double)  voxel size
 * dd     (1x2 double)  distance weighting for the path length of the regions 
 *                      growing (dd(1)) and the intensity (dd(2))
 *                      (default [0.1 10])
 *
 *
 * Examples:
 *
 *
 *
 * _____________________________________________________________________________
 * Robert Dahnke 
 * Structural Brain Mapping Group
 * University Jena
 *
 * $Id: cat_vol_downcut.c 1185 2017-09-13 15:02:47Z gaser $ 
 */
 
#include "mex.h"   
#include "matrix.h"
#include "math.h"
#include "float.h"

/*
#ifndef isnan
#define isnan(a) ((a)!=(a)) 
#endif
*/

/* HELPFUNCTIONS */

/* estimate x,y,z position of index i in an array size sx,sxy=sx*sy... */
void ind2sub(int i,int *x,int *y, int *z, int sxy, int sy) {
  *z = (int)floor( i / (double)sxy ) +1; 
   i = i % (sxy);
  *y = (int)floor( i / (double)sy ) +1;        
  *x = i % sy + 1;
}

int sub2ind(int x,int y, int z, const int s[]) {
  int i=(z-1)*s[0]*s[1] + (y-1)*s[0] + (x-1);
  if (i<0 || i>s[0]*s[1]*s[2]) i=1;
  return i;
}

float abs2(float n) {	if (n<0) return -n; else return n; }
float sign(float n) {	if (n<0) return 1; else return 0; }
float max2(float a, float b) { if (a>b) return a; else return b; }

/* MAINFUNCTION */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs<1)                                                      mexErrMsgTxt("ERROR:cat_vol_downcut: not enought input elements\n");
  if (nrhs>5)                                                      mexErrMsgTxt("ERROR:cat_vol_downcut: to many input elements.\n");
  if (nlhs>2)                                                      mexErrMsgTxt("ERROR:cat_vol_downcut: to many output elements.\n");
  if (mxIsSingle(prhs[0])==0)                                      mexErrMsgTxt("ERROR:cat_vol_downcut: first  input must be an 3d single matrix\n");
  if (mxIsSingle(prhs[1])==0)                                      mexErrMsgTxt("ERROR:cat_vol_downcut: second input must be an 3d single matrix\n");
  if (mxIsDouble(prhs[2])==0 || mxGetNumberOfElements(prhs[2])!=1) mexErrMsgTxt("ERROR:cat_vol_downcut: third input must one double value\n");
  if (nrhs==4 && mxIsDouble(prhs[3])==0)                           mexErrMsgTxt("ERROR:cat_vol_downcut: fourth  input must be an double matrix\n");
  if (nrhs==4 && mxGetNumberOfElements(prhs[3])!=3)                mexErrMsgTxt("ERROR:cat_vol_downcut: fourth input must have 3 Elements");
  if (nrhs==5 && mxIsDouble(prhs[4])==0)                           mexErrMsgTxt("ERROR:cat_vol_downcut: fifth input must be an double matrix\n");
  if (nrhs==5 && mxGetNumberOfElements(prhs[4])!=2)                mexErrMsgTxt("ERROR:cat_vol_downcut: fifht input must have 2 Elements");
    
  /* main informations about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]); 
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = (int)sL[0];
  const int     y  = (int)sL[1];
  const int     xy = x*y;

  const mwSize sSS[] = {1,3}, sdsv[] = {1,2}; 
  mxArray *SS  = mxCreateNumericArray(2,sSS, mxDOUBLE_CLASS,mxREAL); double*S  = mxGetPr(SS);
  mxArray *dsv = mxCreateNumericArray(2,sdsv,mxDOUBLE_CLASS,mxREAL); double*dd = mxGetPr(dsv);
  float dI=0.0; double*SEGd; if (nrhs>=3) {SEGd=mxGetPr(prhs[2]); dI=(float) SEGd[0];}; 
  if (nrhs<4) {S[0]=1; S[1]=1; S[2]=1;} else {S=mxGetPr(prhs[3]);}
  if (nrhs<5) {dd[0]=0.1; dd[1]=10;}    else {dd=mxGetPr(prhs[4]);}
  
  float s1 = abs2((float)S[0]),s2 = abs2((float)S[1]),s3 = abs2((float)S[2]);
  const float   s12  = sqrt( s1*s1  + s2*s2); /* xy - voxel size */
  const float   s13  = sqrt( s1*s1  + s3*s3); /* xz - voxel size */
  const float   s23  = sqrt( s2*s2  + s3*s3); /* yz - voxel size */
  const float   s123 = sqrt(s12*s12 + s3*s3); /* xyz - voxel size */
  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  const int   NI[]  = {  1, -1,  x, -x, xy,-xy, -x-1,-x+1,x-1,x+1, -xy-1,-xy+1,xy-1,xy+1, -xy-x,-xy+x,xy-x,xy+x,  -xy-x-1,-xy-x+1,-xy+x-1,-xy+x+1, xy-x-1,xy-x+1,xy+x-1,xy+x+1};  
  const float ND[]  = { s1, s1, s2, s2, s3, s3,  s12, s12,s12,s12,   s13,  s13, s13, s13,   s23,  s23, s23, s23,     s123,   s123,   s123,   s123,   s123,  s123,  s123,  s123};
 
  int         i,n,ni,u,v,w,nu,nv,nw;
    
  /* main volumes - actual without memory optimation ... */
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); /* label map  */
  plhs[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); /* tissue map (speed) */
    
  /* input variables */
  float*ALAB = (float *)mxGetPr(prhs[0]);	/* label map */
  float*SEG  = (float *)mxGetPr(prhs[1]);	/* tissue map */
  
  /* output variables */
  float*SLAB = (float *)mxGetPr(plhs[0]);	/* label map */
  float*DIST = (float *)mxGetPr(plhs[1]);	/* distance map */
  
  int	  nCV = 0; 		/* # voxel of interest (negative voxel that have to processed)  */
	int 	kll;
	int 	kllv = 2000;
	float DISTN;
	
	/* initialisation of parameter volumes */
	for (i=0;i<nL;i++) { 
		SLAB[i] = ALAB[i]; 
		if (mxIsNaN(SLAB[i])) SLAB[i]=0;
		if (SLAB[i]==0) 					{DIST[i]= FLT_MAX; } 
		else {
			if (SLAB[i]==-FLT_MAX) {DIST[i]=-FLT_MAX; }
			else 										{DIST[i]=0; nCV++;} } 
	}
	/* diffusion */
	int   nC = nCV;
	kll=0;
	while ( nCV>0 && kll<kllv && nC>0 ) {
		kll++; nC=0;
	
		for (i=0;i<nL;i++) { 
			if ( (DIST[i]<=0) && (DIST[i]!=-FLT_MAX ) ) { 
				if (DIST[i]<0) DIST[i]=-DIST[i]; 
				nCV--; /* demark points - also the with zero distance */
				
				ind2sub(i,&u,&v,&w,xy,x); 
				for (n=0;n<26;n++) {
					ni = i+NI[n]; ind2sub(ni,&nu,&nv,&nw,xy,x);
					
					if ( ((ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) )==0 && (SEG[i]+dI>=(SEG[ni])) && ALAB[ni]==0 ) {
            if (nrhs==5) DISTN = DIST[i] + (float)dd[0]*ND[n] + (float)dd[1]*max2(0,SEG[ni]);  
            else         DISTN = DIST[i] + (float)dd[0]*ND[n] + (float)dd[1]*max2(0,4-SEG[ni]);
            
						/*if ( DISTN>0 && (abs2(DIST[ni])>abs2(DISTN)) &&  (( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) )==0)  )  { */
						if (  (DIST[ni]!=-FLT_MAX) && (abs2(DIST[ni])>abs2(DISTN)) )  {	
							if (DIST[ni]>0) nCV++; nC++;
							DIST[ni] = -DISTN;
							SLAB[ni] = SLAB[i];
						}
					}				
				}
				if (DIST[i]==0) DIST[i]=-FLT_MAX; /* demark start points */
				
			}
		}
		
	}

	/*for (int i=0;i<nL;i++) { if (DIST[i]==-FLT_MAX) {DIST[i] = 0; SLAB[i]++; } } */

	/*printf("(%d,%d)",nCV,kll); */
}
