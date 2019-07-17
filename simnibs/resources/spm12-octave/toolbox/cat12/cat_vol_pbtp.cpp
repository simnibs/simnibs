/* quasi-euclidean distance calculation
 * _____________________________________________________________________________
 * [GMT,RPM,WMD,CSFD,II] = eqdist(SEG,WMD,CSFD[,opt])
 *
 * SEG  = (single) segment image with low and high boundary bd
 * GMT  = (single) thickness image
 * RPM  = (single) radial position map
 * WMD  = (single) CSF distance map
 * CSFD = (single) CSF distance map
 * II  = (uint32) index of the inner (WM)  boundary voxel
 *
 * opt.bd   = (single) [low,high] boundary values (default 1.5 and 2.5)
 * opt.CSFD = calculate CSFD
 * opt.PVE  = use PVE informations (0=none,1=fast,2=exact)
 *
 * TODO:
 *  - eikonal distance for subsegmentation (region growing)
 *  - own labeling (
 * _____________________________________________________________________________
 * Robert Dahnke 2009_11
 * Center of Neuroimaging 
 * University Jena
 */

#include "mex.h"   
#include "matrix.h"
#include "math.h"
#include <stdlib.h>

struct opt_type {
	int   CSFD;													/* use CSFD */
	int   PVE;													/* 0, 1=fast, 2=exact */
	float LB, HB, LLB, HLB, LHB, HHB;  	/* boundary */
	int   sL[3];
	// ...
	} opt;

float min(float a, float b) {
  if (a<b) return a; else return b; 
}

float max(float a, float b) {
  if (a>b) return a; else return b; 
}

// get all values of the voxels which are in WMD-range (children of this voxel)  
void pmax(const float GMT[], const float RPM[], const float SEG[], const float ND[], const float WMD, const float SEGI, const int sA, float & maximum) {
  float T[27]; for (int i=0;i<27;i++) T[i]=-1; float n=0.0; maximum=WMD; 

  /* the pure maximum */
  /* (GMT[i]<1e15) && (maximum < GMT[i]) && ((RPM[i]-ND[i]*1.25)<=WMD) && ((RPM[i]-ND[i]*0.5)>WMD) && (SEGI)>=SEG[i] && SEG[i]>1 && SEGI>1.66) */
  for (int i=0;i<=sA;i++) {
    if (  ( GMT[i] < 1e15 ) && ( maximum < GMT[i] ) &&                  /* thickness/WMD of neighbors should be larger */
          ( SEG[i] >= 0.0 ) && ( SEGI>1.2 && SEGI<=2.75 ) &&           /* projection range */
          ( ( ( RPM[i] - ND[i] * 1.2 ) <= WMD ) ) &&                    /* upper boundary - maximum distance */
          ( ( ( RPM[i] - ND[i] * 0.5 ) >  WMD ) || ( SEG[i]<1.5 ) ) &&  /* lower boundary - minimum distance - corrected values outside */
          ( ( ( (SEGI * max(1.0,min(1.2,SEGI-1)) ) >= SEG[i] ) ) || ( SEG[i]<1.5 ) ) )  /* for high values will project data over sulcal gaps */
      { maximum = GMT[i]; }
  }

  /* the mean of the highest values*/
  float maximum2=maximum; float m2n=0; 
  for (int i=0;i<=sA;i++) {
    if ( ( GMT[i] < 1e15 ) && ( (maximum - 1) < GMT[i] ) && 
         ( SEG[i] >= 0.0 ) && ( SEGI>1.2 && SEGI<=2.75 ) && 
         ( ( (RPM[i] - ND[i] * 1.2 ) <= WMD ) ) && 
         ( ( (RPM[i] - ND[i] * 0.5 ) >  WMD ) || ( SEG[i]<1.5 ) ) &&
         ( ( ( (SEGI * max(1.0,min(1.2,SEGI-1)) ) >= SEG[i] ) ) || ( SEG[i]<1.5 ) ) ) 
      { maximum2 = maximum2 + GMT[i]; m2n++; } 
  }
  if ( m2n > 0 )  maximum = (maximum2 - maximum)/m2n;

}




// estimate x,y,z position of index i in an array size sx,sxy=sx*sy...
void ind2sub(int i, int *x, int *y, int *z, int snL, int sxy, int sy) {
  /* not here ... 
   *  if (i<0) i=0; 
   *  if (i>=snL) i=snL-1;
  */
  
  *z = (int)floor( (double)i / (double)sxy ) ; 
   i = i % (sxy);
  *y = (int)floor( (double)i / (double)sy ) ;        
  *x = i % sy ;
}



// main function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs<3) mexErrMsgTxt("ERROR: not enought input elements\n");
  if (nrhs>4) mexErrMsgTxt("ERROR: to many input elements.\n");
  if (nlhs>2) mexErrMsgTxt("ERROR: to many output elements.\n");
  if (mxIsSingle(prhs[0])==0) mexErrMsgTxt("ERROR: first  input must be an 3d single matrix\n");
 
  
  // main informations about input data (size, dimensions, ...)
  const mwSize *sL = mxGetDimensions(prhs[0]); 
  mwSize sSEG[] = {sL[0],sL[1],sL[2]}; 
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = sL[0];
  const int     y  = sL[1];
  const int     xy = x*y;
  const float   s2 = sqrt(2.0);
  const float   s3 = sqrt(3.0);
  const int     nr = nrhs;
  
  // indices of the neighbor Ni (index distance) and euclidean distance NW
  const int   NI[]  = {  0, -1,-x+1, -x,-x-1,  -xy+1,-xy,-xy-1,  -xy+x+1,-xy+x,-xy+x-1,  -xy-x+1,-xy-x,-xy-x-1};  
  const float ND[]  = {0.0,1.0,  s2,1.0,  s2,     s2,1.0,   s2,       s3,   s2,     s3,       s3,   s2,     s3};
  const int   sN  = sizeof(NI)/4;  
  float       DN[sN],DI[sN],GMTN[sN],WMDN[sN],SEGN[sN],DNm;
  
  float 	    du, dv, dw, dnu, dnv, dnw, d, dcf, WMu, WMv, WMw, GMu, GMv, GMw, SEGl, SEGu, tmpfloat;
  int         ni,u,v,w,nu,nv,nw, tmpint, WMC=0, CSFC=0;
    
  // main volumes - actual without memory optimation ...
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[2] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[3] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[4] = mxCreateNumericArray(dL,sL,mxUINT32_CLASS,mxREAL);  
  
  // input variables
  float*SEG  = (float *)mxGetPr(prhs[0]);
  float*WMD  = (float *)mxGetPr(prhs[1]);
  float*CSFD = (float *)mxGetPr(prhs[2]);
  
  /*if ( nrhs>1) {
		tmpint   = (int)mxGetScalar(mxGetField(prhs[1],1,"CSFD"));  printf("X=%d", tmpint);	if ( tmpint!=NULL && (tmpint>=0 && tmpint<=1) ) opt.CSFD = tmpint;   else opt.CSFD 	= 1;
		tmpint   = (int)mxGetScalar(mxGetField(prhs[1],1,"PVE"));   printf("X=%d", tmpint);	if ( tmpint!=NULL && (tmpint>=0 && tmpint<=2) ) opt.PVE  = tmpint; 	 else opt.PVE		= 2;
		tmpfloat = (float)mxGetScalar(mxGetField(prhs[1],1,"LB"));  printf("X=%d", tmpfloat);	if ( tmpfloat!=NULL ) 													opt.LB   = tmpfloat; else opt.LB  	= 1.5;
		tmpfloat = (float)mxGetScalar(mxGetField(prhs[1],1,"HB"));  printf("X=%d", tmpfloat);	if ( tmpfloat!=NULL ) 													opt.HB   = tmpfloat; else opt.HB 		= 2.5;
	} 
	else */{ opt.CSFD = 1;opt.PVE = 2;opt.LB = 1.5;opt.HB	= 2.5; }
	opt.LLB=floor(opt.LB), opt.HLB=ceil(opt.LB), opt.LHB=floor(opt.HB), opt.HHB=ceil(opt.HB);
  
  // output variables
  float        *GMT  = (float *)mxGetPr(plhs[0]);
  float        *RPM  = (float *)mxGetPr(plhs[1]);
  
  // intitialisiation
  for (int i=0;i<nL;i++) {
  	GMT[i] = WMD[i];
    RPM[i] = WMD[i];
		// proof distance input
    if ( SEG[i]>=opt.HB ) WMC++;
    if ( SEG[i]<=opt.LB ) CSFC++;
  }
	if (WMC==0)  mexErrMsgTxt("ERROR: no WM voxel\n");
	if (CSFC==0) opt.CSFD = 0;

  
  
// thickness calcuation
// =============================================================================
  for (int i=0;i<nL;i++) {
    if (SEG[i]>opt.LLB && SEG[i]<opt.HHB) {
      ind2sub(i,&u,&v,&w,nL,xy,x);
      
      // read neighbor values
      for (int n=0;n<sN;n++) {
        ni = i + NI[n];
        ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
        GMTN[n] = GMT[ni]; WMDN[n] = RPM[ni]; SEGN[n] = SEG[ni];
      }

      // find minimum distance within the neighborhood
      pmax(GMTN,WMDN,SEGN,ND,WMD[i],SEG[i],sN,DNm);
      GMT[i] = DNm;
    }
  }
  
  for (int i=nL-1;i>=0;i--) {
    if (SEG[i]>opt.LLB && SEG[i]<opt.HHB) {
      ind2sub(i,&u,&v,&w,nL,xy,x);
      
      // read neighbor values
      for (int n=0;n<sN;n++) {
        ni = i - NI[n];
        ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
        GMTN[n] = GMT[ni]; WMDN[n] = RPM[ni]; SEGN[n] = SEG[ni];
      }

      // find minimum distance within the neighborhood
      pmax(GMTN,WMDN,SEGN,ND,WMD[i],SEG[i],sN,DNm);
      if ( GMT[i] < DNm && DNm>0 ) GMT[i] = DNm;
    }
  }
  
  for (int i=0;i<nL;i++) if (SEG[i]<opt.LB || SEG[i]>opt.HB) GMT[i]=0; //WMD[i]

 
  


// final setings...
// =============================================================================
	float CSFDc = 0, GMTi, CSFDi; // 0.125
	for (int i=0;i<nL;i++) { 
		/* GMT[i] = min(CSFD[i] + WMD[i],GMT[i]); */
		if (SEG[i]>=opt.LB & SEG[i]<=opt.LB) {
			GMTi   = CSFD[i] + WMD[i];	
			CSFDi  = GMT[i]  - WMD[i];
		
			if ( CSFD[i]>CSFDi )	CSFD[i] = CSFDi; 					
			else               		GMT[i]  = GMTi;
		}
	}

 
// estimate RPM
// =============================================================================
	for (int i=0;i<nL;i++) {
		if ( SEG[i]>=opt.HB ) 	
      RPM[i]=1.0; 
		else {
			if ( SEG[i]<=opt.LB || GMT[i]==0.0 ) 
        RPM[i]=0.0;
			else {
				RPM[i] = (GMT[i] - WMD[i]) / GMT[i];
				if (RPM[i]>1.0) RPM[i]=1.0;
				if (RPM[i]<0.0) RPM[i]=0.0;	
			}
		} 
	}
  
}


