/* 
 * This function estimates the Euclidean distance D to the closest  
 * boundary voxel I given by the 0.5 isolevel in B that should contain
 * values between 0 and 1 to define the boundary using PVE. 
 *
 * To align the closest voxel a modified Eikonal distances is estimated 
 * on the field L. 
 * For a correct side alignment a harder option "csf" is used that avoids
 * the diffusion to voxels with greater values of L (L[i]>L[ni] for
 * diffusion).
 *
 * Voxels with NaN and -INF are ignored and produce NaN in D, whereas 
 * in I the default index is used, because I has to be a integer matrix 
 * where NaN is not available. With setnan=0, the result of NAN voxel in 
 * D is changed to INF.
 * 
 *  [D,I] = cat_vol_eidist(B,L,[vx_vol,euclid,csf,setnan,verb])
 * 
 *  D         distance map   (3d-single-matrix)
 *  I         index map      (3d-uint32-matrix)
 *  B         boundary map   (3d-single-matrix)
 *  L         speed map      (3d-single-matrix)
 *  vx_vol    voxel size     (1x3-double-matrix): default=[1 1 1]
 *            ... not tested yet
 *  csf       option         (1x1-double-value): 0-no; 1-yes; default=1
 *  euclid    option         (1x1-double-value): 0-no; 1-yes; default=1
 *            output euclidean or speed map 
 *  setnan    option         (1x1-double-value): 0-no; 1-yes; default=1
 *  verb      option         (1x1-double-value): 0-no, 1-yes; default=0
 * 
 *
 * Small test examples:
 * 1) 
 *   A=zeros(50,50,50,'single'); A(20:30,5:15,20:30)=1; A(20:30,35:45,20:30)=1; 
 *   A=smooth3(A); A(1:5,1:25,:)=nan; A(1:5,26:50,:)=-inf; A(45:50,26:50,:)=inf;
 *   F=ones(size(A),'single'); F(10:40,20,:)=0.5; F(40,10:40,:)=0;
 *
 * 2) 1D-examples
 *   A=zeros(10,20,10,'single'); A(:,1:5,:)=1; A(:,15:end,:)=nan; F=ones(size(A),'single');
 *   A=zeros(10,20,10,'single'); A(:,1:5,:)=1; A(:,6,:)=0.2; A(:,15:end,:)=nan; F=ones(size(A),'single');
 *   A=zeros(10,20,10,'single'); A(:,1:5,:)=1; A(:,15:end,:)=nan; F=ones(size(A),'single');
 *
 *   [D,I]=cat_vol_eidist(A,F,[1 1 1],1,1,0,1);
 *   ds('x2','',1,A,D,D,I,25)
 *
 * see also compile.m
 * _____________________________________________________________________
 * Robert Dahnke 
 * Structural Brain Mapping Group
 * University Jena
 *
 * $Id: cat_vol_eidist.c 1185 2017-09-13 15:02:47Z gaser $ 
 */


/*
 * ToDo:
 * - check anisotropic distance estimation 
 */

/* Modifications, Guilherme Saturnino, 2019
 * Removed referces to MATLAB
 */
 

#include "math.h"   
#include "float.h"
#include "limits.h"

#ifdef _MSC_VER
  #define FINFINITY (FLT_MAX+FLT_MAX);
  static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
  #define FNAN (*(const float *) __nan)
#else
  #define FINFINITY 1.0f/0.0f;
  #define FNAN 0.0f/0.0f
#endif



/* 
 * Estimate x,y,z position of index i in an array size sx,sxy.
 * The index is given in c-notation for a matrix with n elements i is 
 * element of 0,1,...,n-1. Also the x,y,z axis are described in the 
 * c-notation and a matrix A with 20x30x40 elements will have coordinates 
 * with x=0,1,...,19; y=0,1,...,30; and z=0,1,...40. The index i=0 is 
 * identical with x=0,y=0,z=0; whereas the last index in the matrix A 
 * is i=n-1=(20*30*40)-1 and has the coordinates x=19,y=29, and z=39.
 *
 */

void ind2sub(int i, int *x, int *y, int *z, int snL, int sxy, int sy) {
  if (i<0) i=0; 
  if (i>=snL) i=snL-1;
  
  *z = (int)floor( (double)i / (double)sxy ) ; 
   i = i % (sxy);
  *y = (int)floor( (double)i / (double)sy ) ;        
  *x = i % sy ;
}

/* 
 * Estimate index i of a voxel x,y,z in an array size s.
 * See also for ind2sub.
 */
int sub2ind(int x, int y, int z, int s[]) {
  if (x<0) x=0; if (x>s[0]-1) x=s[0]-1; 
  if (y<0) y=0; if (y>s[1]-1) y=s[1]-1; 
  if (z<0) z=0; if (z>s[2]-1) z=s[2]-1; 
  
  return (z)*s[0]*s[1] + (y)*s[0] + (x);
}

float fpow(float x, float y) {
  return (float) pow((double) x,(double) y); 
}

float fsqr(float x) {
  return x*x; 
}

float ffloor(float x) {
  return (float) floor((double) x); 
}

/* 
 * Read out the linear interpolated value of a volume SEG with the size
 * s on the position x,y,z (c-notation). See also ind2sub for details of
 * the c-noation.
 */
float isoval(float SEG[], float x, float y, float z, int s[]){

  int i;
  float seg=0.0, n=0.0;
  float fx = ffloor(x),   fy = ffloor(y),   fz = ffloor(z);
  float cx = ffloor(x+1), cy = ffloor(y+1), cz = ffloor(z+1);
    
  float wfx = cx-x, wfy = cy-y, wfz = cz-z;
  float wcx = x-fx, wcy = y-fy, wcz = z-fz;

  /* value of the 8 neighbors and there distance weight */
  float N[8], W[8];  
  N[0]=SEG[sub2ind((int)fx,(int)fy,(int)fz,s)];  W[0]=wfx * wfy * wfz; 
  N[1]=SEG[sub2ind((int)cx,(int)fy,(int)fz,s)];  W[1]=wcx * wfy * wfz;
  N[2]=SEG[sub2ind((int)fx,(int)cy,(int)fz,s)];  W[2]=wfx * wcy * wfz;
  N[3]=SEG[sub2ind((int)cx,(int)cy,(int)fz,s)];  W[3]=wcx * wcy * wfz;
  N[4]=SEG[sub2ind((int)fx,(int)fy,(int)cz,s)];  W[4]=wfx * wfy * wcz;
  N[5]=SEG[sub2ind((int)cx,(int)fy,(int)cz,s)];  W[5]=wcx * wfy * wcz;
  N[6]=SEG[sub2ind((int)fx,(int)cy,(int)cz,s)];  W[6]=wfx * wcy * wcz; 
  N[7]=SEG[sub2ind((int)cx,(int)cy,(int)cz,s)];  W[7]=wcx * wcy * wcz;
    
  for (i=0; i<8; i++) {
    if ( isnan(N[i])==0 || isinf(N[i])==0 )
      seg = seg + N[i] * W[i]; n+= W[i];
  }
  if ( n>0.0 )
    return seg/n; 
  else 
    return FNAN;
}


/* 
 * MAINFUNCTION [D,I] = cat_vol_eidist(B,L,vx_vol,euclid,csf,setnan,verb])
 */
void vol_eidist(float *D, unsigned int *I, float *B, float *L, int *sL, float *S, int euclid, int csf, int setnan, int verb) {
  const int nL = sL[0]*sL[1]*sL[2];
  const int x  = (int)sL[0];
  const int y  = (int)sL[1];
  const int xy = x*y;
  int sizeL[] = {(int)sL[0],(int)sL[1],(int)sL[2]}; 
  float nanres;
  if (setnan>=1) nanres = FNAN; else nanres = FINFINITY; 
  float s1 = fabs((float)S[0]);
  float s2 = fabs((float)S[1]);
  float s3 = fabs((float)S[2]); /* x,y,z - voxel size */
  const float s12  = sqrt( s1*s1  + s2*s2); /* xy  - voxel size */
  const float s13  = sqrt( s1*s1  + s3*s3); /* xz  - voxel size */
  const float s23  = sqrt( s2*s2  + s3*s3); /* yz  - voxel size */
  const float s123 = sqrt(s12*s12 + s3*s3); /* nL - voxel size */
  
  const int   NI[]  = { 1, -1, x, -x, xy, -xy, -x-1, -x+1, x-1, x+1,      
                       -xy-1, -xy+1, xy-1, xy+1, -xy-x, -xy+x, xy-x, xy+x,
                       -xy-x-1, -xy-x+1, -xy+x-1, -xy+x+1, xy-x-1, xy-x+1,
                       xy+x-1,xy+x+1};  
  const float ND[]  = { s1, s1, s2, s2, s3, s3, s12, s12,s12,s12,
                        s13, s13, s13, s13, s23, s23, s23, s23,  
                        s123, s123, s123, s123, s123, s123, s123, s123};
  const int kllv = (sL[0]+sL[1]+sL[2]);

  
  /* 
   * other variables 
   */
  float dinu, dinv, dinw, dcf, WMu, WMv, WMw, WM, DIN, DINE;  /* distance and intensity variables */
  int   i,n,ni,u,v,w,nu,nv,nw,iu,iv,iw;                       /* nL and index-values of a voxel, one neighbor and the nearest boundary voxel */
  int   nC=sL[0]*sL[1]*sL[2], nCo=INT_MAX;                    /* runtime variables */
  int   kll=0;                                                /* runtime variables */
  int   fast=1;                                               /* stop if number of unvisited points stay constant */

  
  /*
   * Display Initial Parameter
   */
  if ( verb ) {
    printf("\ncat_vol_eidist.c debuging mode:\n  Initialize Parameter: \n");
    printf("    size(B) = %d %d %d\n",sL[0],sL[1],sL[2]); 
    printf("    vx_vol  = %0.0f %0.0f %0.0f\n",s1,s2,s3); 
    printf("    euclid  = %d\n",(int) euclid); 
    printf("    setnan  = %d\n",(int) setnan); 
  }
      
 
  /* 
   * Check input values.
   */
  int vx=0;
  for (i=0;i<nL;i++) { 
    if ( B[i]>=0.5 ) vx++;                                              /* count object voxel */
    
    if ( (L[i]<=0.0 || isnan(L[i])) && B[i]<0.5) B[i] = FNAN;         /* ignore voxel that canot be visited */
    
    if ( isinf(B[i]) && B[i]<0.0 ) B[i] = FNAN;                       /* voxel to ignore */
    if ( B[i]>1.0 ) B[i] = 1.0;                                         /* normalize object */
    if ( B[i]<0.0 ) B[i] = 0.0;                                         /* normalize object */
    I[i] = (unsigned int) i;                                            /* initialize index map */
    
    if ( L[i]<=0.0 || isnan(B[i]) || isnan(L[i]) ) 
      D[i] = FNAN;
    else 
      D[i] = FINFINITY;
  }

  /* 
   * Check if there is an object and return FINFINITY and i=i+1 (for matlab) 
   * for all points if there is no object.
   */
  if ( vx==0 ) { 
    for (i=0;i<nL;i++) { D[i]=nanres; I[i]=(unsigned int)i+1; } 
    printf("WARNING:cat_vol_eidist: Found no object for distance estimation!\n");
    return;
  }
    

  /* 
   * Eikonal distance initialization:
   * In 3D this can be really complex, especially for anisotropic volumes. 
   * So only a simple approximation by the PVE of a voxel is used. 
   */ 
  if ( verb ) printf("  Initialize Eikonal distance estimation and index alignment \n");
  for (i=0;i<nL;i++) { 
    if ( D[i]>0.0 && isnan(D[i])==0 ) { 
      if ( B[i]>=0.5 ) {
        D[i] = 0.0; 
        ind2sub(i,&u,&v,&w,nL,xy,x);
        for (n=0;n<26;n++) {
          ni = i + NI[n];
          ind2sub(ni,&nu,&nv,&nw,nL,xy,x); 

          if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || 
               (abs(nw-w)>1) || (ni==i) )==0 && B[ni]<0.5) {
              DIN = ND[n] * (B[ni] - 0.5) / ( B[ni] - B[i] ); 
              
              if ( fabs(D[ni])>DIN ) {
                D[ni] = -DIN; 
                I[ni] = I[i];
              }
           
          }
        }
      }
    }
  }
 
  
  /* 
   * iterative Eikonal distance estimation
   */ 
  if ( verb ) printf("  Eikonal distance estimation and index alignment \n");
  while ( nC>0 && kll<kllv && (nC!=nCo || fast==0)) {

    kll++; 
    nCo=nC;
    nC=0;

    /*
     * Forward direction:
     */
    for (i=0;i<nL;i++) { 
      if ( D[i]<=0.0 && isnan(D[i])==0) { 
        if ( D[i]<0.0 && isinf(D[i])==0 ) D[i]*=-1.0; /* mark voxel as visited */

        ind2sub(i,&u,&v,&w,nL,xy,x);

        for (n=0;n<26;n++) {
          ni = i + NI[n];
          ind2sub(ni,&nu,&nv,&nw,nL,xy,x); 


          /* Only process this part for real neighbor indices and
           * only if L is lower (to avoid region-growing over sulci)! 
           */
          if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || 
               (abs(nw-w)>1) || (ni==i) || isnan(D[ni]) || B[ni]>=0.5 )==0 &&
               ( csf==0 || L[ni]<=(L[i]+0.1) ) ) { 

            /* new distance */
            DIN = fabs(D[i]) + ND[n] / (FLT_MIN + L[ni]);

            /* use DIN, if the actual value is larger */
            if ( fabs(D[ni])>DIN ) {
              if (D[ni]>0.0) nC++;
              D[ni] = -DIN; 
              I[ni] = I[i];
            }
          }
        }
      }
    }


    /*
     * Backward direction:
     * Same as forward, with the small difference of demarking the start 
     * voxel at the end.
     */
    for (i=nL-1;i>=0;i--) { 
      if ( D[i]<=0.0 && isnan(D[i])==0) { 
        if ( D[i]<0.0 && isinf(D[i])==0 ) D[i]*=-1.0; /* mark voxel as visited */

        ind2sub(i,&u,&v,&w,nL,xy,x);

        for (n=0;n<26;n++) {
          ni = i + NI[n];
          ind2sub(ni,&nu,&nv,&nw,nL,xy,x); 

          /* Only process this part for real neighbor indices and
           * only if L is lower (to avoid region-growing over sulci)! 
           */
          if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || 
               (abs(nw-w)>1) || (ni==i) || isnan(D[ni]) || B[ni]>=0.5 )==0 && 
               ( csf==0 || L[ni]<=(L[i]+0.1) ) ) { 

            /* new distance */
            DIN = fabs(D[i]) + ND[n] / (FLT_MIN + L[ni]);

            /* use DIN, if the actual value is larger */
            if ( fabs(D[ni])>DIN ) {
              if (D[ni]>0.0) nC++;
              D[ni] = -DIN; 
              I[ni] = I[i];
            }
          }
        }

        /* 
         * Demarking the start voxels
         */
        if (D[i]==0.0) D[i]=-FINFINITY; 

      }
    }

    if ( verb && kll<30 ) 
      printf("    nC=%10d, kll=%4d, kllv=%4d\n",nC,kll,kllv);
    if ( nC==nCo ) { csf=0; nCo++; } /* further growing??? */
  }


  /* 
   * Correction of non-visited points due to miss-estimations of the 
   * exact boundary.
   */
  if ( verb ) printf("  Correction of unvisited points \n");  
  for (i=0;i<nL;i++) {
    if ( isinf(D[i]) && D[i]<0.0 ) D[i]=0.0;
    if ( D[i]<0.0 ) D[i]=-FINFINITY; 
  }


  /*
   * Euclidean distance estimation based on the nearest voxel in I.
   */
  if ( verb ) printf("  Euclidean distance estimation \n"); 
  if ( euclid ) {
    for (i=0;i<nL;i++) { 
      if ( isnan(B[i])==0 && isinf(B[i])==0 && D[i]>0.0 && I[i]!=(unsigned int)i ) { 

        ni = (int) I[i];

        ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
        ind2sub(i ,&iu,&iv,&iw,nL,xy,x);

        /* standard euclidean distance between i and closest object point I[i] */
        dinu = (float)iu - (float)nu; dinu *= s1;
        dinv = (float)iv - (float)nv; dinv *= s2;
        dinw = (float)iw - (float)nw; dinw *= s3;
        DIN  = sqrtf(fsqr(dinu) + fsqr(dinv) + fsqr(dinw)); 

        /* For voxels that are not too close to the object the exact 
         * Euclidean distance should be estimated. For closer points
         * the previous distance value is good enough.
         */
        if ( 1 ) {
          /* Estimation of the exact PVE boundary by estimation of a point
           * next to the boundary voxel.
           */
          dinu /= DIN; dinv /= DIN; dinw /= DIN; /* normal vector in normal space */
          dinu /= s1;  dinv /= s2;  dinw /= s3;  /* normal vector for anisotropic space */
          WM = isoval(B,(float)nu + dinu,(float)nv + dinv,(float)nw + dinw,sizeL);

          if ( B[ni]!=WM ) {
            /* estimate new point before border to get the gradient based on this the exact HB */
            dcf = (B[ni] - 0.5) / ( B[ni] - WM );
            WMu = (float)nu + dinu*dcf; 
            WMv = (float)nv + dinv*dcf; 
            WMw = (float)nw + dinw*dcf; 
            WM  = isoval(B,WMu,WMv,WMw,sizeL);

            /* new exact distance to interpolated boundary */ 
            dinu = (float)iu - WMu; dinu *= s1;
            dinv = (float)iv - WMv; dinv *= s2;
            dinw = (float)iw - WMw; dinw *= s3;
            DINE = sqrtf(fsqr(dinu) + fsqr(dinv) + fsqr(dinw)); 

            if ( WM<0.4 || WM>0.6 ) {
              WMu = (float)nu + 0.5*dinu*dcf; 
              WMv = (float)nv + 0.5*dinv*dcf; 
              WMw = (float)nw + 0.5*dinw*dcf; 
              WM  = isoval(B,WMu,WMv,WMw,sizeL);
            }
            
            /* use the new estimated distance only if it is a meanful */
            if ( DINE>0.0 && isnan(DINE)==0 && isinf(DINE)==0 ) {
              D[i] = DINE;
            }
            else {
              /* use the voxelboundary corrected euclidean distance DIN
               * in case of larger values, or otherwise use the initial
               * value used before ... for higher values the speed map
               * lead to non euclidean distances! 
               */
              if ( DIN>2 )
                D[i] = DIN; /* -0.5 */
            }
          }
        }
        else {
        /* simple euclidean distance without PVE for tests */
          D[i] = DIN; 
        }
      }
    }
  }

  
  /*
   * Final corrections
   */
  if ( verb ) printf("  Final corrections \n");  
  for (i=0;i<nL;i++) { 
    /* correct for matlab index */
    if ( I[i]>0 ) I[i]++; else I[i]=1;                        

    /* correction of non-visited or other incorrect voxels */ 
    /* if ( D[i]<0.0 || isnan(D[i]) || isinf(D[i]) ) D[i]=nanres; */

  } 
  if ( verb ) printf("done. \n");   
}
