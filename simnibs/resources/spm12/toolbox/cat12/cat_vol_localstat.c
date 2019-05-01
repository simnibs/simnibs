/* Local Mean, Minimum, Maximum, SD, and Peak estimation
 * ________________________________________________________________________
 * Estimates specific functions volume V in a mask region M. For every 
 * voxel vEM the values of the neigbors of v within a distance nb that 
 * also belong to M were used. 
 * 
 *
 * S = cat_vol_localstat(V,M,nb,stat)
 * 
 * V    (single)    input volume
 * M    (logical)   mask volume
 * nb   (double)    neigbhour distance (1 .. 10)
 * stat (double)    1-mean, 2-min, 3-max, 4-std 
 *                  5-peak1, 6-peak2, 7-peak3
 *                  8-median
 *                  9-?
 *                  10-?
 *
 * ________________________________________________________________________
 * Robert Dahnke 2012_01
 * Center of Neuroimaging 
 * University Jena
 *
 * $Id: cat_vol_localstat.c 1185 2017-09-13 15:02:47Z gaser $ 
 */

#include "mex.h"   
#include "matrix.h"
#include "math.h"
#include "float.h"

#ifndef ROUND
#define ROUND( x ) ((long) ((x) + ( ((x) >= 0) ? 0.5 : (-0.5) ) ))
#endif


#define index(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))


/* qicksort for median */
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

float min(float a, float b) {	if (a<b) return a; else return b; }
float max(float a, float b) {	if (a>b) return a; else return b; }
float abs2(float n) {	if (n<0) return -n; else return n; }        
float pow2(float n) { return n*n;}

/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs<2) mexErrMsgTxt("ERROR:cat_vol_localstat: not enough input elements\n");
  if (nrhs>5) mexErrMsgTxt("ERROR:cat_vol_localstat: too many input elements\n");
  if (nlhs<1) mexErrMsgTxt("ERROR:cat_vol_localstat: not enough output elements\n");
  if (nlhs>2) mexErrMsgTxt("ERROR:cat_vol_localstat: too many output elements\n");

  /* main informations about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]);
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = (int) mxGetNumberOfElements(prhs[0]);
  const int     dB = mxGetNumberOfDimensions(prhs[0]);
  const int     nB = (int) mxGetNumberOfElements(prhs[0]);
  int nh, st;
  
  
  
  if ( dL != 3 || mxIsSingle(prhs[0])==0)                mexErrMsgTxt("ERROR:cat_vol_localstat: first input must be a single 3d matrix\n");
  if ( dB != 3 || mxIsLogical(prhs[1])==0 || nL != nB )  mexErrMsgTxt("ERROR:cat_vol_localstat: second input must be a logical 3d matrix with equal size than input 1\n");
  if (nrhs<3)  nh = 1; else  nh = (int) *mxGetPr(prhs[2]);
  if ( nh > 10 )                                         mexErrMsgTxt("ERROR:cat_vol_localstat: number of neighbors is limited to 10. (Use reduce resolution instead.) \n");
  if (nrhs<4)  st = 1; else st = (int) *mxGetPr(prhs[3]);
  if ( st<1 || st>8 )                                    mexErrMsgTxt("ERROR:cat_vol_localstat: fourth input has to be 1=mean, 2=min, 3=max, 4=std. \n");
  int verb=0;   double*dverb;   if (nrhs>=5) {dverb=mxGetPr(prhs[4]);   verb   = (int) dverb[0]>0.5;};    

  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */

  float nx; 
  float NV[9261],DN[9261], NVmd, NVmn, NVstd; /* nmax ==10 */
  float stdd[3],stdp[3],stdn[3];
  int   stddc[3],stdpc[3],stdnc[3];
  int   di,i,j,k,ind,ni,x,y,z,n,nn,md; /*,HIST1[1000],HIST2[1000]; */
  
  /* in- and output */
  float *D = (float *) mxGetPr(prhs[0]);
  bool  *B = (bool  *) mxGetPr(prhs[1]);
  
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); float *M   = (float *) mxGetPr(plhs[0]);
  plhs[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); float *M2  = (float *) mxGetPr(plhs[1]);
  plhs[2] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); float *M3  = (float *) mxGetPr(plhs[2]);
  for (i=0;i<nL;i++) { M2[i] = 0.0; } 
  for (i=0;i<nL;i++) { M3[i] = 0.0; }  
  for (i=0;i<nL;i++) { 
		if (D[i]==-FLT_MAX || mxIsNaN(D[i])) D[i]=0.0; 	/* correction of non-visited or other incorrect voxels */
    M[i]  = 0.0;
  } 
  
  
  int HISTmax=256; /* HISTmax1=40; HISTmax2=40; (int) (((float) (NVs))/20); if (HISTmax>1000) HISTmax=1000; */ 
  int HISTmin=0;

  
  /*
   * Display Initial Parameter
   */
  if ( verb ) printf("\ncat_vol_localstat.c debuging mode:\n  Initialize Parameter: \n");
  if ( verb ) printf("    size(B) = %d %d %d\n",sL[0],sL[1],sL[2]); 
  if ( verb ) printf("    nb      = %d\n",nh); 
  if ( verb ) printf("    stat    = %d\n",st); 
  
  
  if (st==7) {
    for (i=0;i<nL;i++) { if (HISTmin>D[i]) HISTmin=(int)D[i]; };  HISTmin--;
    for (i=0;i<nL;i++) { if (HISTmax<D[i]) HISTmax=(int)D[i]; };  HISTmax++;
    for (i=0;i<nL;i++) {
      D[i] = (float)ROUND(max(min(D[i],(float)HISTmax),(float)HISTmin));
    }
  }
 
  int HISTn  = HISTmax - HISTmin + 1;
  float *HIST;
  HIST = (float *) malloc(sizeof(float)*HISTn);
  for (nn=0;nn<HISTn;nn++) HIST[nn]=0.0;

  
  if ( verb ) printf(" .. filtering ");   
  /* filter process */
  for (z=0;z<sL[2];z++) for (y=0;y<sL[1];y++) for (x=0;x<sL[0];x++) {
    ind = index(x,y,z,sL);
    if ( B[ind] && mxIsNaN(D[ind])==0 && mxIsInf(D[ind])==0  ) {
      n  = 0;
      
      /* go through all elements in a nh x nh x nh box */
      for (i=-nh;i<=nh;i++) for (j=-nh;j<=nh;j++) for (k=-nh;k<=nh;k++) {
        /* check borders, masks, NaN or Infinities */
        if ( ((x+i)>=0) && ((x+i)<sL[0]) && ((y+j)>=0) && ((y+j)<sL[1]) && ((z+k)>=0) && ((z+k)<sL[2])) {
          ni = index(x+i,y+j,z+k,sL);
          if ( ( B[ni] || st==7 ) && mxIsNaN(D[ni])==0 && mxIsInf(D[ni])==0 ) {
            DN[n] = sqrtf( ( (float) ((i * i) + (j * j) + (k * k)) ) );
            if ( DN[n]<=(float) nh )  { /* && (DN[ni]<=(float)nh) ) { */
              NV[n] = D[ni];
              n++;
            }
          }
        }
      }
     
      NVmn = 0.0; nx = 0.0; 
      
      /* mean */
      if (st==1) { M[ind]=   0.0; for (nn=0;nn<n;nn++) { M[ind]+=NV[nn]; nx++;}; M[ind]/=nx;};
      
      /* mean, min, max */
      if (st==10) {
        M[ind]=0.0; 
        M2[ind]=D[ind]; 
        M3[ind]=D[ind]; 
        for (nn=0;nn<n;nn++) { 
          M[ind]+=NV[nn]; nx++;
          if (NV[nn]<M2[ind]) M2[ind]=NV[nn];
          if (NV[nn]>M3[ind]) M3[ind]=NV[nn];
        }
        M[ind]/=nx;
      }
      
      /* minimum */
      if (st==2) { M[ind]=D[ind]; for (nn=0;nn<n;nn++) { if (NV[nn]<M[ind]) M[ind]=NV[nn];};}; 
      
      /* maximum */
      if (st==3) { M[ind]=D[ind]; for (nn=0;nn<n;nn++) { if (NV[nn]>M[ind]) M[ind]=NV[nn];};}; 
      
      /* standard deviation */
      if (st==4) {
        M[ind]  = 0.0; for (nn=0;nn<n;nn++) { M[ind]+= NV[nn]; nx++;}; M[ind]/=nx; NVmn=M[ind];
        M2[ind] = M[ind];
        NVstd   = 0.0; for (nn=0;nn<n;nn++) { NVstd += (NV[nn]-NVmn)*(NV[nn]-NVmn);};
        M[ind]  = sqrtf((NVstd/(nx-1.0)));
      };
     
      /* ===============================================================
       * experimental histogram functions 
       * ===============================================================
       */ 
      /* max in histogram 1 */
      if (st==5) {
        for (nn=0;nn<HISTn;nn++) {HIST[nn]=0;}
        for (nn=0;nn<n;nn++) {HIST[(int) ROUND( NV[nn]* (float) HISTmax) ]++;}
        M[ind]=0.0; for (nn=0;nn<HISTn;nn++) { if (HIST[nn]>M[ind]) M[ind]=(float) HIST[nn]; }
        M[ind]=M[ind]/(float) HISTmax;
      };  
      
      /* max in histogram 2 */
      if (st==6) {
        NVmn = D[ni]; /*= 0.0; for (nn=0;nn<n;nn++) { NVmn  +=  NV[nn];}; NVmn/=n; */
        
        for (nn=0;nn<HISTn;nn++) {HIST[nn]=0;}
        for (nn=0;nn<n;nn++) {HIST[(int) ROUND( NV[nn]* (float) HISTmax) ]++;}
        /*for (nn=0;nn<n;nn++) if (NV[nn]<MVmn) {HIST1[(int) ROUND( NV[nn]* (float) HISTmax) ]++;} */
        M[ind]=0; for (nn=(int) ROUND( NVmn * (float) HISTmax);nn<HISTn;nn++) { if (HIST[nn]>M[ind]) M[ind]=(float) nn;}
        /*M[ind]=0; for (nn=0;nn<200;nn++) { if (HIST[nn]>M[ind]) M[ind]=(float) nn;} */
        M[ind]=M[ind]/(float) HISTmax;
        M2[ind]=0.0; for (nn=0;nn<(int) ROUND( NVmn * (float) HISTmax);nn++) { if (HIST[nn]>M2[ind]) M2[ind]=(float) nn;}
        M2[ind]=M2[ind]/(float) HISTmax;
      };  
      
      /* max in histogram 3 */
      if (st==7) {
        for (nn=0;nn<HISTn;nn++) {HIST[nn]=0.0;} /* init histogram */
        for (nn=0;nn<n;nn++) {HIST[(int) (NV[nn] - HISTmin)]++;}  /* estimate histogram */
        NVmn=0.0; M[ind]=0.0; for (nn=0;nn<HISTn;nn++) { if (HIST[nn]>NVmn) {NVmn=HIST[nn]; M[ind]=(float) nn;}}
        M[ind]=M[ind] + HISTmin;
      };   
     
      
     /* median */
     if (st==8) {
        if (n>nh*nh*nh) { 
          sort(NV,0,n); 
          md=(int)ROUND(n*0.75);
          NVmd = NV[md];

          
          nx=0;
          M[ind] = 0.0; for (nn=0;nn<md;nn++) { M[ind]+= NV[nn]; nx++;}; M[ind]/=nx; NVmn=M[ind]; /* NVmn=NVmd; */ 
          NVstd  = 0.0; for (nn=0;nn<md;nn++) { NVstd += (NV[nn]-NVmn)*(NV[nn]-NVmn);};
          if ( nx>nh*nh ) M2[ind] = sqrtf((NVstd/(nx-1.0))); else M[ind]=0.0;
          
          nx=0;
          M[ind] = 0.0; for (nn=md;nn<n;nn++) { M[ind]+= NV[nn]; nx++;}; M[ind]/=nx; NVmn=M[ind]; /* NVmn=NVmd; */
          NVstd  = 0.0; for (nn=md;nn<n;nn++) { NVstd += (NV[nn]-NVmn)*(NV[nn]-NVmn);};
          if ( nx>nh*nh ) M[ind] = sqrtf((NVstd/(nx-1.0))); else M[ind]=0.0;

          M[ind]  = M[ind]*2.0;
          M2[ind] = M2[ind]*2.0;
          if (M[ind]>M2[ind]) {NVmn=M2[ind]; M2[ind]=M[ind]; M[ind]=NVmn;}
                  
          /* M2[ind] = max(0,M2[ind] - M[ind]*5); */
          /* M2[ind] = M2[ind] - M[ind]; */
          /* if ( M2[ind]>0.3 )  M2[ind]=0; */
        }

      }
       /* ===============================================================
       * experimental noise/signal functions 
       * ===============================================================
       */ 
      if (st==9) {
        /* estimation of noise and strukture a subregions (normally the WM) */
        
        /* mean (or better median) intensity in the area */
        /*  M[ind] = 0.0; for (nn=0;nn<n;nn++) { M[ind]+= NV[nn]; nx++;}; M[ind]/=nx; NVmn=M[ind]; */
        if (n>1) { if (n==2) {
            NVmd = (NV[0] + NV[1]) / 2.0;  
          }
          else {
            sort(NV,0,n); 
           NVmd = NV[(int)(n/2.0)];
          }
        }
        M[ind] = 0.0; for (nn=0;nn<n;nn++) { M[ind]+= NV[nn]; nx++;}; M[ind]/=nx; NVmn=M[ind];
        NVstd  = 0.0; for (nn=0;nn<n;nn++) { NVstd += (NV[nn]-NVmn)*(NV[nn]-NVmn);};
        M[ind] = (float)sqrt((double)(NVstd/(nx-1.0)));
        
        /* estimation of noise for each dimention */
        stdd[0]=0.0; stdd[1]=0.0;stdd[2]=0.0;
        stdp[0]=0.0; stdp[1]=0.0;stdp[2]=0.0;
        stdn[0]=0.0; stdn[1]=0.0;stdn[2]=0.0;

        stddc[0]=0; stddc[1]=0;stddc[2]=0;
        stdpc[0]=0; stdpc[1]=0;stdpc[2]=0;
        stdnc[0]=0; stdnc[1]=0;stdnc[2]=0;
        
        for (di=0;di<=2;di++) { 
          for (i=-nh*2;i<=nh*2;i++) {
            if (di==0)
              ni = index(x+i,y,z,sL);
            else {
              if (di==1)
                ni = index(x,y+i,z,sL);
              else
                ni = index(x,y,z+1,sL);
            }

            if ( (x+i)>=0 && (x+i)<sL[0] ) {
              if (B[ni] && mxIsNaN(D[ni])==0 && mxIsInf(D[ni])==0 ) {
                stdd[di] += (D[ni]-NVmn)*(D[ni]-NVmn); stddc[di]++;
                if (D[ind]>NVmd) {
                  stdp[di] += (D[ni]-NVmn)*(D[ni]-NVmn); stdpc[di]++;
                }
                else {
                  stdn[di] += (D[ni]-NVmn)*(D[ni]-NVmn); stdnc[di]++;
                }
              }
            }
          }
          
          if ( stddc[di]>1 ) {stdd[di]=sqrtf((stdd[di]/(stddc[di]-1)));} else {stdd[di] = 0.0;}
          if ( stdpc[di]>1 ) {stdp[di]=sqrtf((stdp[di]/(stdpc[di]-1)));} else {stdp[di] = 0.0;}
          if ( stdnc[di]>1 ) {stdn[di]=sqrtf((stdn[di]/(stdnc[di]-1)));} else {stdn[di] = 0.0;}
        }
       /* sort(stdd,0,2); */
        
        M[ind]=stdd[0];
        M2[ind]=stdd[1];
        /* M[ind]  = (stdd[0] + stdd[1] + stdd[2])/3 * 2; */
        /* M[ind]  = (stdp[0] + stdp[1] + stdp[2])/3 * 2; */
        /* M2[ind] = (stdn[0] + stdn[1] + stdn[2])/3 * 2 -  M[ind]; */ 
         
      }
      
      
      if ((M[ind]==-FLT_MAX) || (D[ind]==FLT_MAX) || (mxIsNaN(D[ind])) ) M[ind]=0.0;
    

    }
    else {
      M[ind] = 0.0;
      /*D[ind];*/
     /* SD[ind] = 0; */      
    }
    
    
  } 
  for (i=0;i<nL;i++) { 
		if ( M[i]==-FLT_MAX || M[i]==FLT_MAX || mxIsNaN(M[i]) ) M[i]=0.0; 	/* correction of non-visited or other incorrect voxels */
	} 
  
  /* free(HIST); */
  if ( verb ) printf("done. \n");   
}


