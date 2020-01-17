/*
 * Robert Dahnke
 * $Id: cat_vbdist.c 1172 2017-08-31 13:19:22Z gaser $ 
 *
 */

/* voxelbased euclidean distance calculation
 * ________________________________________________________________________
 * Calculates the euclidean distance without PVE to an object in P with a 
 * boundary of 0.5.
 * 
 *  [D,I,L] = vbdist(P[,R])
 *  
 *  P (single)  input image with zero for non elements and uint8 values for labels
 *  R (logical) range for distance calculation
 *  D (single)  distance image
 *  L (uint8)   label map
 *  I (uint32)  index of nearest point
 * ________________________________________________________________________
 * Robert Dahnke 2010_01
 * Center of Neuroimaging 
 * University Jena
 */

/* Modifications, Guilherme Satutnino, 2019
 * Removed MATLAB references
 *
*/

#include "math.h"
#include "float.h"


/* estimate minimum of A and its index in A */
void pmin(float A[], int sA, float *minimum, int *index)
{
  int i; 
  *minimum = FLT_MAX; *index = 0; /* printf("%d ",sizeof(A)/8); */
  for(i=0; i<sA; i++) {
    if ((A[i]>0) && (*minimum>A[i]))
    { 
      *minimum = A[i]; 
      *index   = i;
    }
  }
}


/* estimate x,y,z position of index i in an array size sx,sxy=sx*sy... */
void ind2sub_vb(int i, int *x, int *y, int *z, int sxy, int sy) {
  *z = (int)floor( (double)i / (double)sxy ) +1; 
   i = i % (sxy);
  *y = (int)floor( (double)i / (double)sy ) +1;        
  *x = i % sy + 1;
}

/* main function */
void vbdist(float *D, unsigned int *I, unsigned char *L, float *V, unsigned char *R, int *sL, float *S){
  /* main informations about input data (size, dimensions, ...) */
  const int     nL = sL[0]*sL[1]*sL[2];
  const int     x  = (int)sL[0];
  const int     y  = (int)sL[1];
  const int     xy = x*y;

  float s1 = fabs((float)S[0]),s2 = fabs((float)S[1]),s3 = fabs((float)S[2]);
  const float   s12  = sqrt( s1*s1  + s2*s2); /* xy - voxel size */
  const float   s13  = sqrt( s1*s1  + s3*s3); /* xz - voxel size */
  const float   s23  = sqrt( s2*s2  + s3*s3); /* yz - voxel size */
  const float   s123 = sqrt(s12*s12 + s3*s3); /* xyz - voxel size */
  /*printf("%1.2f,%1.2f,%1.2f - %1.2f,%1.2f,%1.2f - %1.2f",s1,s2,s3,s12,s23,s13,s123); */
  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  const int   NI[] = {  0, -1,-x+1, -x,-x-1,  -xy+1,-xy,-xy-1,  -xy+x+1,-xy+x,-xy+x-1,  -xy-x+1,-xy-x,-xy-x-1};  
  const float ND[] = {0.0, s1, s12, s2, s12,    s13, s3,  s13,     s123,  s23,   s123,     s123,  s23,   s123};
  const int   sN = sizeof(NI)/4;    
  float       DN[sN];
  float       DNm = FLT_MAX;
  int i, n, ni, DNi = 0;

  unsigned char e255 = 0; 
  
  /* intitialisation */
  for (i=0;i<nL;i++) 
  {
    if (V[i]>=0.5) D[i]=0.0; else D[i]=FLT_MAX; 
    if (V[i]>255.0)  
    {
      if (e255==0) 
      {
        printf("Warning: First parameter of vbdist > 255!\n"); 
        e255 = 1;
      }
      V[i] = 255;
      
    }
    L[i]=(unsigned char) ceil(V[i]);
    I[i]=(unsigned int)i;
  }

  int u,v,w,nu,nv,nw; 
  for (i=0;i<nL;i++) 
  {
    if ( (D[i]>0) && R[i] )
    {
      ind2sub_vb(i,&u,&v,&w,xy,x);
      
      /* read neighbor values */
      for (n=0;n<sN;n++)
      {
        ni = i + NI[n];
        ind2sub_vb(ni,&nu,&nv,&nw,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) ni=i;
        DN[n] = D[ni] + ND[n];
      }

      /* find minimum distance within the neighborhood */
      pmin(DN,sN,&DNm,&DNi);

      /* update values */
      if (DNi>0) {
        L[i] = L[i+NI[DNi]];
        I[i] = (unsigned int) I[i+NI[DNi]];
        D[i] = DNm; 
        ind2sub_vb((int)I[i],&nu,&nv,&nw,xy,x); 
        D[i] = sqrt(pow((float)(u-nu)*s1,2) + pow((float)(v-nv)*s2,2) + pow((float)(w-nw)*s3,2));
        /* hier muss die genauerer Berechnung mit pve rein! */
      }
    }

  }
  for (i=nL-1;i>0;i--)
  {
    if ( (D[i]>0) && R[i] )
    {
      ind2sub_vb(i,&u,&v,&w,xy,x);

      /* read neighbor values */
      for (n=0;n<sN;n++)
      {
        ni = i - NI[n];
        ind2sub_vb(ni,&nu,&nv,&nw,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) ni=i;
        DN[n] = D[ni] + ND[n];
      }

      /* find minimum distance within the neighborhood */
      pmin(DN,sN,&DNm,&DNi);

      /* update values */
      if (DNi>0) {
        L[i] = L[i-NI[DNi]];
        I[i] = (unsigned int)  I[i-NI[DNi]];
        D[i] = DNm; 
        ind2sub_vb((int)I[i],&nu,&nv,&nw,xy,x); 
        D[i] = sqrt(pow((float)(u-nu)*s1,2) + pow((float)(v-nv)*s2,2) + pow((float)(w-nw)*s3,2));
      }
    }
  }

  
  /* euclidean calcuation + PVE information */
  for (i=0;i<nL;i++) 
  {
  /*  if ( (D[i]>0) && (nrhs==1 || (R[i]>=1) ) )  D[i] = D[i] + (1 - V[I[i]]);
    if (D[i]<=0.5) D[i]=0; else D[i]=D[i]-0.5; */
    I[i]++;
  }

}


