#include <stdio.h>
#include <stdlib.h>
#include <math.h>    
#include <float.h>

/*

Program to extract simply connected surfaces ("no holes, handles, loops, or self intersections") from segmented 3D images.
Steve Haker, Surgical Planning Lab, Brigham and Women's Hospital, Harvard University.
Distance transform by Andre Robatino, Surgical Planning Lab, Brigham and Women's Hospital, Harvard University.

Please send bug reports to haker@bwh.harvard.edu.

*/

typedef struct _genus0parameters
  {

  /* INPUT parameters */

  unsigned short * input;          /* Pointer to data from which to extract surface. */

  int return_adjusted_label_map ;  /* return an adjusted label map?  I.e. fill in output[]?  Default 1. */
  int return_surface ;             /* fill in triangles[] and vertices[]?  Default 1. */

  int dims[3];                     /* Dimensions of data set (cols,rows,depth). */

  unsigned short value;                     /* Value of label defining surface to extract. */

  unsigned short alt_value;                 /* Value to assign to adjusted voxels */
  unsigned short contour_value;             /* Value to set contour voxels to, if they are part of original labelmap.  Default alt_value */
  unsigned short alt_contour_value;         /* Value to set contour voxels to, if they are not part of original labelmap.  Default alt_value */

  int cut_loops;                   /* cut loops instead of patching holes? Sorry, it's either/or globally... Default 0.*/

  int pad[3];                      /* pad around volume after cropping.  Default is {2,2,2}. */

  int connectivity;                /* Connectivity model. 6 or 18.  Default is 6. */

  int any_genus;                   /* If not zero, just get the surface regardless of its genus. Essentially marching cubes. */
                                   /* If zero, we want genus zero surfaces.  Default is 0 */

                                   /* Only used if return_surface not zero... */
  int biggest_component;           /* If not zero, extract biggest connected triangulated component. */
                                   /* Default is 1 */
  
  int connected_component;         /* extract largest connected component from input[] data _before_ processing. Default 1 */

  float *ijk2ras;                  /* 4x4 matrix, ijk2ras[12..15]={0,0,0,1},... you know what I mean. */
                                   /* Multiply vertices by this, i.e. ijk2ras*vertices.  RAS means Right, Anterior, Superior.*/ 
                                   /* You can pass NULL to use the identity matrix. */
                                   /* i=[0..dim[0]], j=[0...dims[1]] ("down the side of each image"), k=[0...dims[2]] */
                                   /* Default is NULL (identity matrix) */


  float extraijkscale[3];          /* Voxel scaling during process.  Does not affect the coordinates of vertices created. */
                                   /* Just seems to help reduce the corrections needed sometimes if you scale in the */
                                   /* direction of the scan (normally k, i.e. set extraijkscale[2]>1 ).  Default is {1,1,1} */

  int verbose;                     /* If not zero, show progress of algorithm. Default is 0. */


  /* OUTPUT results */

                                   /* Only used if  return_adjusted_label_map not zero... */
  unsigned short * output;         /* Pointer to the place where adjusted label map will be put. */
                                   /* If NULL when genus0 is called, space will be allocated. */
                                   /* Else, it's assumed to be allocated already. */

                                   /* Only used if return_surface not zero... */
  int vert_count;                  /* Will hold number of vertices. */
  float * vertices;                /* Will hold output vertices. (vert_count cols x 3 rows).  Space will be allocated */
  int tri_count;                   /* Will hold number of triangles. */
  int * triangles;                 /* Will hold output triangles. (tri_count cols x 3 rows).  Space will be allocated */

  /* PRIVATE */
  int calloced_output;             /* A flag to remember if *output was calloced, so it can be freed upon destruction */

  } genus0parameters;


/* public stuff */
extern void genus0init(genus0parameters * g0);        /* Initialize fields in *g0 to their defaults.  Must be called before genus0(). */
extern int genus0(genus0parameters * g0);             /* Call the algorithm.  Do the work.  Returns 0 on success, 1 on failure. */
extern void genus0destruct(genus0parameters * g0);    /* Frees *vertices and *triangles, and frees *output if it was calloced */
