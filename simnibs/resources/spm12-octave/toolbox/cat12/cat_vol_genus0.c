/* Genus0 topology correction
 * _____________________________________________________________________________
 *
 * $Id: cat_vol_genus0.c 1185 2017-09-13 15:02:47Z gaser $ 
 */

#include "mex.h"   
#include "matrix.h"
#include "math.h"
#include "float.h"
#include "genus0.h"

/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs<2) mexErrMsgTxt("ERROR:cat_vol_genus0: At least two input elements necessary.\n");
  if (nlhs<1) mexErrMsgTxt("ERROR:cat_vol_genus0: At least one output element necessary.\n");

  /* main informations about input data (size, dimensions, ...) */
  const mwSize *sL  = mxGetDimensions(prhs[0]);
  const int     dL  = mxGetNumberOfDimensions(prhs[0]);
    
  unsigned short *input;
  
  if ( dL  != 3 || mxIsSingle(prhs[0])==0)        mexErrMsgTxt("ERROR:cat_vol_genus0: first input must be a single 3d matrix\n");
  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  int i, j, sz = sL[0]*sL[1]*sL[2];
        
  /* in- and output */
  float *vol = (float *) mxGetPr(prhs[0]);
  float th   = (float) mxGetScalar(prhs[1]);

  plhs[0]  = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  float *M = (float *) mxGetPr(plhs[0]);
    
  input = (unsigned short *) mxMalloc(sizeof(unsigned short)*sz);

  /* use orientation from isosurface */
  float ijk2ras[16] = {0, 1, 0, 0,
                       1, 0, 0, 0,
                       0, 0, 1, 0,
                       0, 0, 0, 1};
  
  for (i=0; i<sz; i++) {
    if (vol[i] > th)
         input[i] = 1;
    else input[i] = 0;
  }
  
  genus0parameters g0[1];  /* need an instance of genus0parameters */

  genus0init(g0);  /* initialize the instance, set default parameters */

  /* set some parameters/options */
  for(j= 0; j <dL; j++ ) g0->dims[j] = sL[j];
  
  g0->return_adjusted_label_map = 1;
  g0->cut_loops = 0;
  g0->connectivity = 6;
  g0->connected_component = 1;
  g0->input = input;
  g0->value = 1;
  g0->alt_value = 1;
  g0->contour_value = 1;
  g0->alt_contour_value = 1;
  g0->biggest_component = 1;
  g0->pad[0] = g0->pad[1] = g0->pad[2] = 2;
  g0->ijk2ras = ijk2ras;
  g0->verbose = 1;
  g0->return_surface = 0;
  g0->extraijkscale[0] = g0->extraijkscale[1] = g0->extraijkscale[2] = 1;

  if (nrhs==3) g0->any_genus = (int) mxGetScalar(prhs[2]);
  else         g0->any_genus = 0;

  if (nlhs==3) g0->return_surface = 1;

  /* in case of error return empty variables */
  if (genus0(g0)) {
    mwSize dims[2];

    dims[0] = 0; dims[1] = 0;
    plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);

    dims[0] = 0; dims[1] = 0;
    plhs[2] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
        
    return;
  }

  g0->cut_loops = 1;
  g0->connectivity = 18;
  g0->value = 1;
  g0->alt_value = 0;
  g0->contour_value = 1;
  g0->alt_contour_value = 0;

  for (i = 0; i < sz; i++)
    input[i] = (unsigned int)g0->output[i];

  if (genus0(g0)) {
    mwSize dims[2];

    dims[0] = 0; dims[1] = 0;
    plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);

    dims[0] = 0; dims[1] = 0;
    plhs[2] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
        
    return;
  }

  for (i=0; i<sz; i++) 
    M[i] = (float)g0->output[i];

  if (nlhs==3) {
    mwSize dims[2];

    dims[0] = g0->tri_count; dims[1] = 3;
    plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);

    dims[0] = g0->vert_count; dims[1] = 3;
    plhs[2] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);

    double *Tris  = (double *) mxGetPr(plhs[1]);
    double *Verts = (double *) mxGetPr(plhs[2]);  
    
    /* return Tris and Verts and add 1 for matlab use */
    for (i=0; i<3*g0->tri_count; i++) 
      Tris[i] = (double)g0->triangles[i] + 1;
    for (i=0; i<3*g0->vert_count; i++) 
      Verts[i] = (double)g0->vertices[i] + 1.0;
  }

  genus0destruct(g0);
  mxFree(input);
 
}


