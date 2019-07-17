/*
 * Christian Gaser
 * $Id: vollib.c 510 2013-10-14 16:33:37Z gaser $ 
 *
 * This code is a substantially modified version of spm_conv_vol.c 
 * from J. Ashburner
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <float.h>

#define RINT(A) floor((A)+0.5)
#ifndef isfinite
#define isfinite(x) ((x) * (x) >= 0.) /* check for NaNs */
#endif

static void 
convxy(double out[], int xdim, int ydim, double filtx[], double filty[], int fxdim, int fydim, int xoff, int yoff, double buff[])
{
  int x,y,k;
  for(y=0; y<ydim; y++)
  {
    for(x=0; x<xdim; x++)
    {
      buff[x] = out[x+y*xdim];
      if (!isfinite(buff[x]))
        buff[x] = 0.0;
    }
    for(x=0; x<xdim; x++)
    {
      double sum1 = 0.0;
      int fstart, fend;
      fstart = ((x-xoff >= xdim) ? x-xdim-xoff+1 : 0);
      fend = ((x-(xoff+fxdim) < 0) ? x-xoff+1 : fxdim);

      for(k=fstart; k<fend; k++)
        sum1 += buff[x-xoff-k]*filtx[k];
      out[x+y*xdim] = sum1;
    }
  }
  for(x=0; x<xdim; x++)
  {
    for(y=0; y<ydim; y++)
      buff[y] = out[x+y*xdim];

    for(y=0; y<ydim; y++)
    {
      double sum1 = 0.0;
      int fstart, fend;
      fstart = ((y-yoff >= ydim) ? y-ydim-yoff+1 : 0);
      fend = ((y-(yoff+fydim) < 0) ? y-yoff+1 : fydim);

      for(k=fstart; k<fend; k++)
        sum1 += buff[y-yoff-k]*filty[k];
      out[y*xdim+x] = sum1;
    }
  }
}

static void 
convxy_float(float out[], int xdim, int ydim, float filtx[], float filty[], int fxdim, int fydim, int xoff, int yoff, float buff[])
{
  int x,y,k;
  for(y=0; y<ydim; y++)
  {
    for(x=0; x<xdim; x++)
    {
      buff[x] = out[x+y*xdim];
      if (!isfinite(buff[x]))
        buff[x] = 0.0;
    }
    for(x=0; x<xdim; x++)
    {
      float sum1 = 0.0;
      int fstart, fend;
      fstart = ((x-xoff >= xdim) ? x-xdim-xoff+1 : 0);
      fend = ((x-(xoff+fxdim) < 0) ? x-xoff+1 : fxdim);

      for(k=fstart; k<fend; k++)
        sum1 += buff[x-xoff-k]*filtx[k];
      out[x+y*xdim] = sum1;
    }
  }
  for(x=0; x<xdim; x++)
  {
    for(y=0; y<ydim; y++)
      buff[y] = out[x+y*xdim];

    for(y=0; y<ydim; y++)
    {
      float sum1 = 0.0;
      int fstart, fend;
      fstart = ((y-yoff >= ydim) ? y-ydim-yoff+1 : 0);
      fend = ((y-(yoff+fydim) < 0) ? y-yoff+1 : fydim);

      for(k=fstart; k<fend; k++)
        sum1 += buff[y-yoff-k]*filty[k];
      out[y*xdim+x] = sum1;
    }
  }
}

int 
convxyz_double(double *iVol, double filtx[], double filty[], double filtz[],
  int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
  double *oVol, int dims[3])
{
  double *tmp, *buff, **sortedv, *obuf;
  int xy, z, y, x, k, fstart, fend, startz, endz;
  int xdim, ydim, zdim;

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  tmp = (double *)malloc(sizeof(double)*xdim*ydim*fzdim);
  buff = (double *)malloc(sizeof(double)*((ydim>xdim) ? ydim : xdim));
  sortedv = (double **)malloc(sizeof(double *)*fzdim);

  if((tmp == NULL) || (buff == NULL) || (sortedv == NULL)) {
    printf("Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
  endz   = zdim+fzdim+zoff-1;

  for (z=startz; z<endz; z++)
  {
    double sum2 = 0.0;

    if (z >= 0 && z<zdim)
    {
      for (y=0;y<ydim;y++) for (x=0;x<xdim;x++)
        tmp[((z%fzdim)*xdim*ydim)+(y*xdim)+x] = iVol[(z*xdim*ydim)+(y*xdim)+x];   
      convxy(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
        filtx, filty, fxdim, fydim, xoff, yoff, buff);
    }
    if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim)
    {
      fstart = ((z >= zdim) ? z-zdim+1 : 0);
      fend = ((z-fzdim < 0) ? z+1 : fzdim);

      for(k=0; k<fzdim; k++)
      {
        int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
        sortedv[k] = &(tmp[z1*xdim*ydim]);
      }

      for(k=fstart, sum2=0.0; k<fend; k++)
        sum2 += filtz[k];

      obuf = &oVol[(z-fzdim-zoff+1)*ydim*xdim];
      if (sum2)
      {
        for(xy=0; xy<xdim*ydim; xy++)
        {
          double sum1=0.0;
          for(k=fstart; k<fend; k++)
            sum1 += filtz[k]*sortedv[k][xy];

          obuf[xy] = sum1/sum2;
        }
      }
      else
        for(xy=0; xy<xdim*ydim; xy++)
          obuf[xy] = 0.0;
    }
  }
  free(tmp);
  free(buff);
  free(sortedv);
  return(0);
}

int 
convxyz_float(float *iVol, float filtx[], float filty[], float filtz[],
  int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
  float *oVol, int dims[3])
{
  float *tmp, *buff, **sortedv, *obuf;
  int xy, z, y, x, k, fstart, fend, startz, endz;
  int xdim, ydim, zdim;

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  tmp = (float *)malloc(sizeof(float)*xdim*ydim*fzdim);
  buff = (float *)malloc(sizeof(float)*((ydim>xdim) ? ydim : xdim));
  sortedv = (float **)malloc(sizeof(float *)*fzdim);

  if((tmp == NULL) || (buff == NULL) || (sortedv == NULL)) {
    printf("Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
  endz   = zdim+fzdim+zoff-1;

  for (z=startz; z<endz; z++)
  {
    float sum2 = 0.0;

    if (z >= 0 && z<zdim)
    {
      for (y=0;y<ydim;y++) for (x=0;x<xdim;x++)
        tmp[((z%fzdim)*xdim*ydim)+(y*xdim)+x] = iVol[(z*xdim*ydim)+(y*xdim)+x];   
        convxy_float(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
            filtx, filty, fxdim, fydim, xoff, yoff, buff);
    }
    if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim)
    {
      fstart = ((z >= zdim) ? z-zdim+1 : 0);
      fend = ((z-fzdim < 0) ? z+1 : fzdim);

      for(k=0; k<fzdim; k++)
      {
        int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
        sortedv[k] = &(tmp[z1*xdim*ydim]);
      }

      for(k=fstart, sum2=0.0; k<fend; k++)
        sum2 += filtz[k];

      obuf = &oVol[(z-fzdim-zoff+1)*ydim*xdim];
      if (sum2)
      {
        for(xy=0; xy<xdim*ydim; xy++)
        {
          float sum1=0.0;
          for(k=fstart; k<fend; k++)
            sum1 += filtz[k]*sortedv[k][xy];

          obuf[xy] = sum1/sum2;
        }
      }
      else
        for(xy=0; xy<xdim*ydim; xy++)
          obuf[xy] = 0.0;
    }
  }
  free(tmp);
  free(buff);
  free(sortedv);
  return(0);
}

int
convxyz_uint8(unsigned char *iVol, double filtx[], double filty[], double filtz[],
  int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
  unsigned char *oVol, int dims[3])
{
  double *tmp, *buff, **sortedv;
  int xy, z, y, x, k, fstart, fend, startz, endz;
  int xdim, ydim, zdim;
  double tmp2;
  unsigned char *obuf;

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  tmp = (double *)malloc(sizeof(double)*xdim*ydim*fzdim);
  buff = (double *)malloc(sizeof(double)*((ydim>xdim) ? ydim : xdim));
  sortedv = (double **)malloc(sizeof(double *)*fzdim);

  if((tmp == NULL) || (buff == NULL) || (sortedv == NULL)) {
    printf("Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
  endz   = zdim+fzdim+zoff-1;

  for (z=startz; z<endz; z++)
  {
    double sum2 = 0.0;

    if (z >= 0 && z<zdim)
    {
      for (y=0;y<ydim;y++) for (x=0;x<xdim;x++)
        tmp[((z%fzdim)*xdim*ydim)+(y*xdim)+x] = (double)iVol[(z*xdim*ydim)+(y*xdim)+x];   
      convxy(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
        filtx, filty, fxdim, fydim, xoff, yoff, buff);
    }
    if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim)
    {
      fstart = ((z >= zdim) ? z-zdim+1 : 0);
      fend = ((z-fzdim < 0) ? z+1 : fzdim);

      for(k=0; k<fzdim; k++)
      {
        int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
        sortedv[k] = &(tmp[z1*xdim*ydim]);
      }

      for(k=fstart, sum2=0.0; k<fend; k++)
        sum2 += filtz[k];

      obuf = oVol;
      obuf = &obuf[(z-fzdim-zoff+1)*ydim*xdim];
      if (sum2)
      {
        for(xy=0; xy<xdim*ydim; xy++)
        {
          double sum1=0.0;
          for(k=fstart; k<fend; k++)
            sum1 += filtz[k]*sortedv[k][xy];
          tmp2 = sum1/sum2;
          if (tmp2<0.0) tmp2 = 0.0;
          else if (tmp2>255.0) tmp2 = 255.0;
          obuf[xy] = (unsigned char)RINT(tmp2);
        }
      }
      else
        for(xy=0; xy<xdim*ydim; xy++)
          obuf[xy] = 0;
    }
  }
  free(tmp);
  free(buff);
  free(sortedv);
  return(0);
}

int 
convxyz_int16(signed short *iVol, double filtx[], double filty[], double filtz[],
  int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
  signed short *oVol, int dims[3])
{
  double *tmp, *buff, **sortedv;
  int xy, z, y, x, k, fstart, fend, startz, endz;
  int xdim, ydim, zdim;
  double tmp2;
  signed short *obuf;

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  tmp = (double *)malloc(sizeof(double)*xdim*ydim*fzdim);
  buff = (double *)malloc(sizeof(double)*((ydim>xdim) ? ydim : xdim));
  sortedv = (double **)malloc(sizeof(double *)*fzdim);

  if((tmp == NULL) || (buff == NULL) || (sortedv == NULL)) {
    printf("Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
  endz   = zdim+fzdim+zoff-1;

  for (z=startz; z<endz; z++)
  {
    double sum2 = 0.0;

    if (z >= 0 && z<zdim)
    {
      for (y=0;y<ydim;y++) for (x=0;x<xdim;x++)
        tmp[((z%fzdim)*xdim*ydim)+(y*xdim)+x] = (double)iVol[(z*xdim*ydim)+(y*xdim)+x];   
      convxy(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
        filtx, filty, fxdim, fydim, xoff, yoff, buff);
    }
    if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim)
    {
      fstart = ((z >= zdim) ? z-zdim+1 : 0);
      fend = ((z-fzdim < 0) ? z+1 : fzdim);

      for(k=0; k<fzdim; k++)
      {
        int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
        sortedv[k] = &(tmp[z1*xdim*ydim]);
      }

      for(k=fstart, sum2=0.0; k<fend; k++)
        sum2 += filtz[k];

      obuf = oVol;
      obuf = &obuf[(z-fzdim-zoff+1)*ydim*xdim];
      if (sum2)
      {
        for(xy=0; xy<xdim*ydim; xy++)
        {
          double sum1=0.0;
          for(k=fstart; k<fend; k++)
            sum1 += filtz[k]*sortedv[k][xy];
          tmp2 = sum1/sum2;
          if (tmp2<-32768.0) tmp2 = -32768.0;
          else if (tmp2>32767.0) tmp2 = 32767.0;
          obuf[xy] = (signed short)RINT(tmp2);
        }
      }
      else
        for(xy=0; xy<xdim*ydim; xy++)
          obuf[xy] = 0;
    }
  }
  free(tmp);
  free(buff);
  free(sortedv);
  return(0);
}

int 
convxyz_int32(signed int *iVol, double filtx[], double filty[], double filtz[],
  int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
  signed int *oVol, int dims[3])
{
  double *tmp, *buff, **sortedv;
  int xy, z, y, x, k, fstart, fend, startz, endz;
  int xdim, ydim, zdim;
  double tmp2;
  signed int *obuf;

  xdim = dims[0];
  ydim = dims[1];
  zdim = dims[2];

  tmp = (double *)malloc(sizeof(double)*xdim*ydim*fzdim);
  buff = (double *)malloc(sizeof(double)*((ydim>xdim) ? ydim : xdim));
  sortedv = (double **)malloc(sizeof(double *)*fzdim);

  if((tmp == NULL) || (buff == NULL) || (sortedv == NULL)) {
    printf("Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
  endz   = zdim+fzdim+zoff-1;

  for (z=startz; z<endz; z++)
  {
    double sum2 = 0.0;

    if (z >= 0 && z<zdim)
    {
      for (y=0;y<ydim;y++) for (x=0;x<xdim;x++)
        tmp[((z%fzdim)*xdim*ydim)+(y*xdim)+x] = (double)iVol[(z*xdim*ydim)+(y*xdim)+x];   
      convxy(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
        filtx, filty, fxdim, fydim, xoff, yoff, buff);
    }
    if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim)
    {
      fstart = ((z >= zdim) ? z-zdim+1 : 0);
      fend = ((z-fzdim < 0) ? z+1 : fzdim);

      for(k=0; k<fzdim; k++)
      {
        int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
        sortedv[k] = &(tmp[z1*xdim*ydim]);
      }

      for(k=fstart, sum2=0.0; k<fend; k++)
        sum2 += filtz[k];

      obuf = oVol;
      obuf = &obuf[(z-fzdim-zoff+1)*ydim*xdim];
      if (sum2)
      {
        for(xy=0; xy<xdim*ydim; xy++)
        {
          double sum1=0.0;
          for(k=fstart; k<fend; k++)
            sum1 += filtz[k]*sortedv[k][xy];
          tmp2 = sum1/sum2;
          if (tmp2<-2147483648.0) tmp2 = -2147483648.0;
          else if (tmp2>2147483647.0) tmp2 = 2147483647.0;
          obuf[xy] = (signed int)RINT(tmp2);
        }
      }
      else
        for(xy=0; xy<xdim*ydim; xy++)
          obuf[xy] = 0;
    }
  }
  free(tmp);
  free(buff);
  free(sortedv);
  return(0);
}

void
morph_erode_uint8(unsigned char *vol, int dims[3], int niter, int th)
{
  double filt[3]={1,1,1};
  int i,j;
  
  /* threshold input */
  for (j=0;j<dims[2]*dims[1]*dims[0];j++)
    vol[j] = (vol[j]>th);

  for (i=0;i<niter;i++) {
    convxyz_uint8(vol,filt,filt,filt,3,3,3,-1,-1,-1,vol,dims);
    for (j=0;j<dims[2]*dims[1]*dims[0];j++)
      vol[j] = (vol[j]>=9);
  }
}


void
morph_dilate_uint8(unsigned char *vol, int dims[3], int niter, int th)
{
  double filt[3]={1,1,1};
  int i,x,y,z,j,band,dims2[3];
  unsigned char *buffer;

  /* add band with zeros to image to avoid clipping */  
  band = niter;
  band = 0;
  for (i=0;i<3;i++) dims2[i] = dims[i] + 2*band;
  
  buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);

  if(buffer == NULL) {
    printf("Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  memset(buffer,0,sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);
  
  /* threshold input */
  for (z=0;z<dims[2];z++) for (y=0;y<dims[1];y++) for (x=0;x<dims[0];x++) 
    buffer[((z+band)*dims2[0]*dims2[1])+((y+band)*dims2[0])+x+band] = (vol[(z*dims[0]*dims[1])+(y*dims[0])+x]>th);

  for (i=0;i<niter;i++) {
    convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims);
    for (j=0;j<dims2[2]*dims2[1]*dims2[0];j++)
      buffer[j] = (buffer[j]>0);
  }
  
  /* return image */
  for (z=0;z<dims[2];z++) for (y=0;y<dims[1];y++) for (x=0;x<dims[0];x++) 
    vol[(z*dims[0]*dims[1])+(y*dims[0])+x] = buffer[((z+band)*dims[0]*dims[1])+((y+band)*dims[0])+x+band];
  
  free(buffer);
}

void
morph_close_uint8(unsigned char *vol, int dims[3], int niter, int th)
{
  morph_dilate_uint8(vol, dims, niter, th);
  morph_erode_uint8(vol, dims, niter, 0);  
}

void
morph_open_uint8(unsigned char *vol, int dims[3], int niter, int th)
{
  morph_erode_uint8(vol, dims, niter, th);
  morph_dilate_uint8(vol, dims, niter, 0);
}

void
morph_close_double(double *vol, int dims[3], int niter, double th)
{
  unsigned char *buffer;
  int i;

  buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims[0]*dims[1]*dims[2]);

  if(buffer == NULL) {
    printf("Memory allocation error\n");
    exit(EXIT_FAILURE);
  }
  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    buffer[i] = (unsigned char) (vol[i] > th);
        
  morph_dilate_uint8(buffer, dims, niter, 0);
  morph_erode_uint8(buffer, dims, niter, 0);

  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    vol[i] = (double)buffer[i];
    
  free(buffer);
}

void
morph_open_double(double *vol, int dims[3], int niter, double th)
{
  unsigned char *buffer;
  int i;

  buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims[0]*dims[1]*dims[2]);

  if(buffer == NULL) {
    printf("Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    buffer[i] = (unsigned char) (vol[i] > th);
        
  morph_erode_uint8(buffer, dims, niter, 0);
  morph_dilate_uint8(buffer, dims, niter, 0);

  for (i=0;i<dims[2]*dims[1]*dims[0];i++)
    vol[i] = (double)buffer[i];
    
  free(buffer);
}

/* First order hold resampling - trilinear interpolation */
void 
subsample_double(double *in, double *out, int dim_in[3], int dim_out[3])
{
  int i, x, y, z;
  double k111,k112,k121,k122,k211,k212,k221,k222;
  double dx1, dx2, dy1, dy2, dz1, dz2, xi, yi, zi, samp[3];
  int off1, off2, xcoord, ycoord, zcoord;

  for (i=0; i<3; i++) {
    if(dim_out[i] > dim_in[i]) samp[i] = ceil((double)dim_out[i]/(double)dim_in[i]);
    else                       samp[i] = 1.0/(ceil((double)dim_in[i]/(double)dim_out[i]));
  }
  
  for (z=0; z<dim_out[2]; z++) {
    zi = 1.0+(double)z/samp[2];
    for (y=0; y<dim_out[1]; y++) {
      yi = 1.0+(double)y/samp[1];
      for (x=0; x<dim_out[0]; x++) {
        xi = 1.0+(double)x/samp[0];
        i = z*dim_out[0]*dim_out[1] + y*dim_out[0] + x;

        if (zi>=0 && zi<dim_in[2] && yi>=0 && yi<dim_in[1] && xi>=0 && xi<dim_in[0])  {
          xcoord = (int)floor(xi); dx1=xi-(double)xcoord; dx2=1.0-dx1;
          ycoord = (int)floor(yi); dy1=yi-(double)ycoord; dy2=1.0-dy1;
          zcoord = (int)floor(zi); dz1=zi-(double)zcoord; dz2=1.0-dz1;

          off1 = xcoord-1 + dim_in[0]*(ycoord-1 + dim_in[1]*(zcoord-1));
          k222 = (double)in[off1]; k122 = (double)in[off1+1]; off2 = off1+dim_in[0];
          k212 = (double)in[off2]; k112 = (double)in[off2+1]; off1+= dim_in[0]*dim_in[1];
          k221 = (double)in[off1]; k121 = (double)in[off1+1]; off2 = off1+dim_in[0];
          k211 = (double)in[off2]; k111 = (double)in[off2+1];

          out[i] = ((((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1))*dz2
                         + (((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1))*dz1);
                 
        } else out[i] = 0;
      }
    }
  }
}

/* First order hold resampling - trilinear interpolation */
void 
subsample_float(float *in, float *out, int dim_in[3], int dim_out[3])
{
  int i, x, y, z;
  float k111,k112,k121,k122,k211,k212,k221,k222;
  float dx1, dx2, dy1, dy2, dz1, dz2, xi, yi, zi, samp[3];
  int off1, off2, xcoord, ycoord, zcoord;

  for (i=0; i<3; i++) {
    if(dim_out[i] > dim_in[i]) samp[i] = ceil((float)dim_out[i]/(float)dim_in[i]);
    else                       samp[i] = 1.0/(ceil((float)dim_in[i]/(float)dim_out[i]));
  }
  
  for (z=0; z<dim_out[2]; z++) {
    zi = 1.0+(float)z/samp[2];
    for (y=0; y<dim_out[1]; y++) {
      yi = 1.0+(float)y/samp[1];
      for (x=0; x<dim_out[0]; x++) {
        xi = 1.0+(float)x/samp[0];
        i = z*dim_out[0]*dim_out[1] + y*dim_out[0] + x;

        if (zi>=0 && zi<dim_in[2] && yi>=0 && yi<dim_in[1] && xi>=0 && xi<dim_in[0])  {
          xcoord = (int)floor(xi); dx1=xi-(float)xcoord; dx2=1.0-dx1;
          ycoord = (int)floor(yi); dy1=yi-(float)ycoord; dy2=1.0-dy1;
          zcoord = (int)floor(zi); dz1=zi-(float)zcoord; dz2=1.0-dz1;

          off1 = xcoord-1 + dim_in[0]*(ycoord-1 + dim_in[1]*(zcoord-1));
          k222 = (float)in[off1]; k122 = (float)in[off1+1]; off2 = off1+dim_in[0];
          k212 = (float)in[off2]; k112 = (float)in[off2+1]; off1+= dim_in[0]*dim_in[1];
          k221 = (float)in[off1]; k121 = (float)in[off1+1]; off2 = off1+dim_in[0];
          k211 = (float)in[off2]; k111 = (float)in[off2+1];

          out[i] = ((((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1))*dz2
                         + (((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1))*dz1);
                 
        } else out[i] = 0;
      }
    }
  }
}

void
smooth_double(double *vol, int dims[3], double separations[3], double s[3], int use_mask)
{
  int i;
  double xsum, ysum, zsum;
  double *x, *y, *z;
  int xyz[3], nvol, sum_mask;
  double *mask;
  unsigned char *mask2;
  
  nvol = dims[0]*dims[1]*dims[2];

  for(i=0; i<3; i++) {
    s[i] /= separations[i];
    if(s[i] < 1.0) s[i] = 1.0;
    s[i] /= sqrt(8.0*log(2.0));
    xyz[i] = (int) RINT(6.0*s[i]);
  }
  
  x = (double *) malloc(sizeof(double)*((2*xyz[0])+1));
  y = (double *) malloc(sizeof(double)*((2*xyz[1])+1));
  z = (double *) malloc(sizeof(double)*((2*xyz[2])+1));
  
  /* build mask for masked smoothing */
  if(use_mask) {
    mask  = (double *) malloc(sizeof(double)*nvol);
    mask2 = (unsigned char *) malloc(sizeof(unsigned char)*nvol);
    sum_mask = 0;
    for(i=0; i<nvol; i++) {
      if(vol[i] > 0.0) {
        mask[i]  = 1.0;
        mask2[i] = 1;
        sum_mask++;
      } else {
        mask[i]  = 0.0;
        mask2[i] = 0;
      }
    }
  }
  
  for(i=-xyz[0]; i <= xyz[0]; i++) x[i+xyz[0]] = (double)i;
  for(i=-xyz[1]; i <= xyz[1]; i++) y[i+xyz[1]] = (double)i;
  for(i=-xyz[2]; i <= xyz[2]; i++) z[i+xyz[2]] = (double)i;
  
  xsum = 0.0; ysum = 0.0; zsum = 0.0;
  for(i=0; i < ((2*xyz[0])+1); i++) {
    x[i] = exp(-pow(x[i],2) / (2.0*pow(s[0],2)));
    xsum += x[i];
  }
  for(i=0; i < ((2*xyz[1])+1); i++) {
    y[i] = exp(-pow(y[i],2) / (2.0*pow(s[1],2)));
    ysum += y[i];
  }
  for(i=0; i < ((2*xyz[2])+1); i++) {
    z[i] = exp(-pow(z[i],2) / (2.0*pow(s[2],2)));
    zsum += z[i];
  }
  
  for(i=0; i < ((2*xyz[0])+1); i++) x[i] /= xsum;
  for(i=0; i < ((2*xyz[1])+1); i++) y[i] /= ysum;
  for(i=0; i < ((2*xyz[2])+1); i++) z[i] /= zsum;
  
  convxyz_double(vol,x,y,z,((2*xyz[0])+1),((2*xyz[1])+1),((2*xyz[2])+1),-xyz[0],-xyz[1],-xyz[2],vol,dims);
  if(use_mask) {
    /* only smooth mask if mask has values > 0 */
    if(sum_mask>0) {
      convxyz_double(mask,x,y,z,((2*xyz[0])+1),((2*xyz[1])+1),((2*xyz[2])+1),-xyz[0],-xyz[1],-xyz[2],mask,dims);
      for(i=0; i<nvol; i++) {
        if(mask2[i]>0) vol[i] /= (double)mask[i];   
        else vol[i] = 0.0; 
      }
    }
    free(mask);
    free(mask2);
  }
  
  free(x);
  free(y);
  free(z);

}

void
smooth_float(float *vol, int dims[3], float separations[3], float s[3], int use_mask)
{
  int i;
  float xsum, ysum, zsum;
  float *x, *y, *z;
  int xyz[3], nvol, sum_mask;
  float *mask;
  unsigned char *mask2;
  
  nvol = dims[0]*dims[1]*dims[2];

  for(i=0; i<3; i++) {
    s[i] /= separations[i];
    if(s[i] < 1.0) s[i] = 1.0;
    s[i] /= sqrt(8.0*log(2.0));
    xyz[i] = (int) RINT(6.0*s[i]);
  }
  
  x = (float *) malloc(sizeof(float)*((2*xyz[0])+1));
  y = (float *) malloc(sizeof(float)*((2*xyz[1])+1));
  z = (float *) malloc(sizeof(float)*((2*xyz[2])+1));
  
  /* build mask for masked smoothing */
  if(use_mask) {
    mask  = (float *) malloc(sizeof(float)*nvol);
    mask2 = (unsigned char *) malloc(sizeof(unsigned char)*nvol);
    sum_mask = 0;
    for(i=0; i<nvol; i++) {
      if(vol[i] > 0.0) {
        mask[i]  = 1.0;
        mask2[i] = 1;
        sum_mask++;
      } else {
        mask[i]  = 0.0;
        mask2[i] = 0;
      }
    }
  }
  
  for(i=-xyz[0]; i <= xyz[0]; i++) x[i+xyz[0]] = (float)i;
  for(i=-xyz[1]; i <= xyz[1]; i++) y[i+xyz[1]] = (float)i;
  for(i=-xyz[2]; i <= xyz[2]; i++) z[i+xyz[2]] = (float)i;
  
  xsum = 0.0; ysum = 0.0; zsum = 0.0;
  for(i=0; i < ((2*xyz[0])+1); i++) {
    x[i] = exp(-pow(x[i],2) / (2.0*pow(s[0],2)));
    xsum += x[i];
  }
  for(i=0; i < ((2*xyz[1])+1); i++) {
    y[i] = exp(-pow(y[i],2) / (2.0*pow(s[1],2)));
    ysum += y[i];
  }
  for(i=0; i < ((2*xyz[2])+1); i++) {
    z[i] = exp(-pow(z[i],2) / (2.0*pow(s[2],2)));
    zsum += z[i];
  }
  
  for(i=0; i < ((2*xyz[0])+1); i++) x[i] /= xsum;
  for(i=0; i < ((2*xyz[1])+1); i++) y[i] /= ysum;
  for(i=0; i < ((2*xyz[2])+1); i++) z[i] /= zsum;
  
  convxyz_float(vol,x,y,z,((2*xyz[0])+1),((2*xyz[1])+1),((2*xyz[2])+1),-xyz[0],-xyz[1],-xyz[2],vol,dims);
  if(use_mask) {
    /* only smooth mask if mask has values > 0 */
    if(sum_mask>0) {
      convxyz_float(mask,x,y,z,((2*xyz[0])+1),((2*xyz[1])+1),((2*xyz[2])+1),-xyz[0],-xyz[1],-xyz[2],mask,dims);
      for(i=0; i<nvol; i++) {
        if(mask2[i]>0) vol[i] /= (float)mask[i];   
        else vol[i] = 0.0; 
      }
    }
    free(mask);
    free(mask2);
  }
  
  free(x);
  free(y);
  free(z);

}

void
smooth_subsample_double(double *vol, int dims[3], double separations[3], double s[3], int use_mask, int samp)
{
  int i, nvol_samp, nvol;
  int dims_samp[3];
  double *vol_samp, separations_samp[3];
  
  /* define grid dimensions */
  for(i=0; i<3; i++) dims_samp[i] = (int) ceil((dims[i]-1)/((double) samp))+1;
  for(i=0; i<3; i++) separations_samp[i] = separations[i]*((double)dims[i]/(double)dims_samp[i]);

  nvol  = dims[0]*dims[1]*dims[2];
  nvol_samp  = dims_samp[0]*dims_samp[1]*dims_samp[2];
  vol_samp  = (double *)malloc(sizeof(double)*nvol_samp);

  subsample_double(vol, vol_samp, dims, dims_samp);    
  smooth_double(vol_samp, dims_samp, separations_samp, s, use_mask);
  subsample_double(vol_samp, vol, dims_samp, dims);    

  free(vol_samp);
}

void
smooth_subsample_float(float *vol, int dims[3], float separations[3], float s[3], int use_mask, int samp)
{
  int i, nvol_samp, nvol;
  int dims_samp[3];
  float *vol_samp;
  float separations_samp[3];
  
  /* define grid dimensions */
  for(i=0; i<3; i++) dims_samp[i] = (int) ceil((dims[i]-1)/((float) samp))+1;
  for(i=0; i<3; i++) separations_samp[i] = separations[i]*((float)dims[i]/(float)dims_samp[i]);

  nvol  = dims[0]*dims[1]*dims[2];
  nvol_samp  = dims_samp[0]*dims_samp[1]*dims_samp[2];
  vol_samp  = (float *)malloc(sizeof(float)*nvol_samp);

  subsample_float(vol, vol_samp, dims, dims_samp);    
  smooth_float(vol_samp, dims_samp, separations_samp, s, use_mask);
  subsample_float(vol_samp, vol, dims_samp, dims);    

  free(vol_samp);
}
