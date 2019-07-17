/*
 * Christian Gaser
 * $Id: Pve.c 603 2014-09-19 08:31:06Z gaser $ 
 *
 */

/* This PVE calculation is a modified version from
 * the PVE software bundle:
 * Copyright (C) Jussi Tohka, Institute of Signal Processing, Tampere University of
 * Technology, 2002 - 2004.
 * P.O. Box 553, FIN-33101, Finland
 * E-mail: jussi.tohka@tut.fi
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Amap.h"

void Pve5(double *src, unsigned char *prob, unsigned char *label, double *mean, int *dims)
{
  int x,y,z,z_area,y_dims,ind;
  double w;
  unsigned char new_val[MAX_NC];
  
  int area = dims[0]*dims[1];
  int vol = area*dims[2];
    
  for (z = 1; z < dims[2]-1; z++) {
    z_area = z*area;
    for (y = 1; y < dims[1]-1; y++) {
      y_dims = y*dims[0];
      for (x = 1; x < dims[0]-1; x++) {
        ind = z_area + y_dims + x;

        switch(label[ind]) {
        case 0: /* BG */
          new_val[CSFLABEL-1] = 0;
          new_val[GMLABEL-1]  = 0;
          new_val[WMLABEL-1]  = 0;
          break;
        case CSFLABEL: /* CSF */
          new_val[CSFLABEL-1] = 255;
          new_val[GMLABEL-1]  = 0;
          new_val[WMLABEL-1]  = 0;
          label[ind] = (unsigned char) ROUND(255.0/3.0);
          break;
        case GMLABEL: /* GM */
          new_val[CSFLABEL-1] = 0;
          new_val[GMLABEL-1]  = 255;
          new_val[WMLABEL-1]  = 0;
          label[ind] = (unsigned char) ROUND(2.0*255.0/3.0);
          break;
        case WMLABEL: /* WM */
          new_val[CSFLABEL-1] = 0;
          new_val[GMLABEL-1]  = 0;
          new_val[WMLABEL-1]  = 255;
          label[ind] = 255;
          break;
        case GMCSFLABEL: /* GMCSF */
          w = (src[ind] - mean[CSFLABEL-1])/(mean[GMLABEL-1]-mean[CSFLABEL-1]);
          if(w > 1.0) w = 1.0; if(w < 0.0) w = 0.0;
          new_val[CSFLABEL-1] = (unsigned char) ROUND(255.0*(1-w));
          new_val[GMLABEL-1]  = (unsigned char) ROUND(255.0*w);
          new_val[WMLABEL-1]  = 0;
          label[ind] = ROUND(255.0/3.0*(1.0 + w));
          break;
        case WMGMLABEL: /* WMGM */
          w = (src[ind] - mean[GMLABEL-1])/(mean[WMLABEL-1]-mean[GMLABEL-1]);
          if(w > 1.0) w = 1.0; if(w < 0.0) w = 0.0;
          new_val[CSFLABEL-1] = 0;
          new_val[GMLABEL-1]  = (unsigned char) ROUND(255.0*(1-w));
          new_val[WMLABEL-1]  = (unsigned char) ROUND(255.0*w);
          label[ind] = ROUND(255.0/3.0*(2.0 + w));
          break;
        }

        prob[          ind] = new_val[CSFLABEL-1];
        prob[vol +     ind] = new_val[GMLABEL-1];
        prob[(2*vol) + ind] = new_val[WMLABEL-1];
        
        /* set old probabilities for mixed classes to zero */
        prob[(3*vol) + ind] = 0;
        prob[(4*vol) + ind] = 0;
        
      }
    }
  }  
}

void Pve6(double *src, unsigned char *prob, unsigned char *label, double *mean, int *dims)
{
  int x,y,z,z_area,y_dims,ind;
  double w;
  unsigned char new_val[MAX_NC];
  
  int area = dims[0]*dims[1];
  int vol = area*dims[2];
    
  for (z = 1; z < dims[2]-1; z++) {
    z_area = z*area;
    for (y = 1; y < dims[1]-1; y++) {
      y_dims = y*dims[0];
      for (x = 1; x < dims[0]-1; x++) {
        ind = z_area + y_dims + x;

        switch(label[ind]) {
        case 0: /* BG */
          new_val[CSFLABEL] = 0;
          new_val[GMLABEL]  = 0;
          new_val[WMLABEL]  = 0;
          break;
        case CSFLABEL+1: /* CSF */
          new_val[CSFLABEL] = 255;
          new_val[GMLABEL]  = 0;
          new_val[WMLABEL]  = 0;
          label[ind] = (unsigned char) ROUND(255.0/3.0);
          break;
        case GMLABEL+1: /* GM */
          new_val[CSFLABEL] = 0;
          new_val[GMLABEL]  = 255;
          new_val[WMLABEL]  = 0;
          label[ind] = (unsigned char) ROUND(2.0*255.0/3.0);
          break;
        case WMLABEL+1: /* WM */
          new_val[CSFLABEL] = 0;
          new_val[GMLABEL]  = 0;
          new_val[WMLABEL]  = 255;
          label[ind] = 255;
          break;
        case BKGCSFLABEL+1: /* BKGCSF */
          w = src[ind]/mean[CSFLABEL];
          if(w > 1.0) w = 1.0; if(w < 0.0) w = 0.0;
          new_val[CSFLABEL] = (unsigned char) ROUND(255.0*w);
          new_val[GMLABEL]  = 0;
          new_val[WMLABEL]  = 0;
          label[ind] = ROUND(255.0/3.0*w);
          break;
        case GMCSFLABEL+1: /* GMCSF */
          w = (src[ind] - mean[CSFLABEL])/(mean[GMLABEL]-mean[CSFLABEL]);
          if(w > 1.0) w = 1.0; if(w < 0.0) w = 0.0;
          new_val[CSFLABEL] = (unsigned char) ROUND(255.0*(1-w));
          new_val[GMLABEL]  = (unsigned char) ROUND(255.0*w);
          new_val[WMLABEL]  = 0;
          label[ind] = ROUND(255.0/3.0*(1.0 + w));
          break;
        case WMGMLABEL+1: /* WMGM */
          w = (src[ind] - mean[GMLABEL])/(mean[WMLABEL]-mean[GMLABEL]);
          if(w > 1.0) w = 1.0; if(w < 0.0) w = 0.0;
          new_val[CSFLABEL] = 0;
          new_val[GMLABEL]  = (unsigned char) ROUND(255.0*(1-w));
          new_val[WMLABEL]  = (unsigned char) ROUND(255.0*w);
          label[ind] = ROUND(255.0/3.0*(2.0 + w));
          break;
        }

        prob[          ind] = new_val[CSFLABEL];
        prob[vol +     ind] = new_val[GMLABEL];
        prob[(2*vol) + ind] = new_val[WMLABEL];
        
        /* set old probabilities for mixed classes to zero */
        prob[(3*vol) + ind] = 0;
        prob[(4*vol) + ind] = 0;
        prob[(5*vol) + ind] = 0;
        
      }
    }
  }  
}
