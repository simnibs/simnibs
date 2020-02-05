# cython: language_level=3
import sys
import numpy as np
cimport numpy as np
from libc.math cimport isnan
from libc.float cimport FLT_MAX


cdef extern from "sanlm_float.c":
    void anlm(float* ima, int v, int f, int use_rician, const int* dims)


cdef extern from "cat_vol_median3.c":
    void vol_median_3(float *M, int nrhs, float *D, int *sL,
                      unsigned char *Bi, unsigned char *Bn,
                      float sf, float bil, float bih, float bnl, float bnh)

cdef extern from "cat_vol_eidist.c":
    void vol_eidist(float *D, unsigned int *I, float *B, float *L, int *sL, float *S, int euclid, int csf, int setnan, int verb)

cdef extern from "cat_vol_localstat.c":
    void vol_localstat(float *M, float *M2, float*M3, float *D, unsigned char *B,  int *sL, int nh, int st, int verb)

cdef extern from "cat_vol_pbtp.c":
    void vol_pbtp(float *GMT, float *RMP, float * SEG, float *WMD, float *CSFD, int *sL)

cdef extern from "genus0.h":
    struct _genus0parameters:
        unsigned short * input
        int return_adjusted_label_map
        int return_surface
        int dims[3]
        unsigned short value
        unsigned short alt_value
        unsigned short contour_value
        unsigned short alt_contour_value
        int cut_loops
        int pad[3]
        int connectivity
        int any_genus
        int biggest_component
        int connected_component
        float *ijk2ras
        float extraijkscale[3]
        int verbose
        unsigned short * output
        int vert_count
        float * vertices
        int tri_count
        int * triangles
    ctypedef _genus0parameters genus0parameters
    void genus0init(genus0parameters * g0)
    int genus0(genus0parameters * g0)
    void genus0destruct(genus0parameters * g0)

cdef extern from "cat_vbdist.c":
    void vbdist(float *D, unsigned int *I, unsigned char *L, float *V, unsigned char *R, int *sL, float *S)


def sanlm(image, v, f, use_rician=False):
    ''' ANLM Denoising
    Unlike the CAT12 version, does not operate in-place. Instead, it creates a new image
    '''
    cdef np.ndarray[float, ndim=3] image_f = image.astype(np.float32)
    cdef np.ndarray[int, ndim=1] dims = np.array(image.shape, np.int32)
    anlm(&image_f[0, 0, 0], int(v), int(f), int(use_rician), &dims[0])
    return image_f



def cat_vol_median3(D, Bi=None, Bn=None, sf=0, Bi_low=-FLT_MAX, Bi_high=FLT_MAX,
                    Bn_low=-FLT_MAX, Bn_high=FLT_MAX):
    ''' Interface to the the cat_vol_median3 function

     Median Filter for a 3d single image D. Bi is used to mask voxels for the 
     filter process, whereas Bn is used to mask voxels that are used as 
     neighbors in the filter process. Both mask can changed by intensity 
     threshold (Bi_low,Bi_high,Bn_low,Bn_high) for D.
    
      M = cat_vol_median3(D[,Bi,Bn,sf,Bi_low,Bi_high,Bn_low,Bn_high])
    
      D      (single)   3d matrix for filter process 
      Bi     (logical)  3d matrix that marks voxels that should be filtered
      Bn     (logical)  3d matrix that marks voxels that are used as neighbors 
      sf     (double)   threshold that is used to filter the result
                          sf=0: no filter
                          sf<0: only smaller changes
                          sf>0: only bigger changes
      Bi_low  (double)  low  threshold in D for filtering (add to Bi)
      Bi_high (double)  high threshold in D for filtering (add to Bi)
      Bn_low  (double)  low  threshold in D for neighbors (add to Bn)
      Bn_high (double)  high threshold in D for neighbors (add to Bn)
    
     Used slower quicksort for median calculation, because the faster median 
     of the median application implementation leads to wrong results. 

 
    Original code by Robert Dahnke, Center of Neuroimaging University Jena
    '''
    # main informations about input data (size, dimensions, ...)
    cdef int nrhs = 1
    cdef np.ndarray[float, ndim=3] D_ = np.array(D, np.float32, order='F')
    cdef np.ndarray[int, ndim=1] sL_ = np.array(D.shape, np.int32)
    cdef np.ndarray[np.uint8_t, ndim=3] Bi_, Bn_
    cdef float sf_ = float(sf)
    cdef float Bi_low_ = float(Bi_low)
    cdef float Bi_high_ = float(Bi_high)
    cdef float Bn_low_ = float(Bn_low)
    cdef float Bn_high_ = float(Bn_high)

    cdef np.ndarray[float, ndim=3] M = np.zeros_like(D_)
    
    # Optional args
    if Bi is not None:
        assert Bi.shape == D.shape
        Bi_ = np.array(Bi, np.uint8, order='F')
        nrhs = 2
    else:
        Bi_ = np.zeros((1,1,1), np.uint8, order='F')

    if Bn is not None:
        assert Bn.shape == D.shape
        Bn_ = np.array(Bn, np.uint8, order='F')
        nrhs = 3
    else:
        Bn_ = np.zeros((1,1,1), np.uint8)

    vol_median_3(&M[0, 0, 0], nrhs, &D_[0, 0 ,0], &sL_[0],
                 &Bi_[0, 0 ,0], &Bn_[0, 0 ,0], sf_,
                 Bi_low_, Bi_high_, Bn_low_, Bn_high_)

    return np.ascontiguousarray(M)


 
def cat_vol_eidist(B, L, vx_vol=(1, 1, 1), csf=1, euclid=1, setnan=1, verb=0):
    ''' Interface to the the cat_vol_eidist function
    Original code by Robert Dahnke, Center of Neuroimaging University Jena

     [D,I] = cat_vol_eidist(B,L,[vx_vol,euclid,csf,setnan,verb])
     
      D         distance map   (3d-single-matrix)
      I         index map      (3d-uint32-matrix) ! 1-based, C notation
      B         boundary map   (3d-single-matrix)
      L         speed map      (3d-single-matrix)
      vx_vol    voxel size     (1x3-double-matrix): default=[1 1 1]
                ... not tested yet
      csf       option         (1x1-double-value): 0-no; 1-yes; default=1
      euclid    option         (1x1-double-value): 0-no; 1-yes; default=1
                output euclidean or speed map 
      setnan    option         (1x1-double-value): 0-no; 1-yes; default=1
      verb      option         (1x1-double-value): 0-no, 1-yes; default=0

    '''
    assert B.shape == L.shape, 'Distance and Speed matrix must have the same size'

    cdef np.ndarray[float, ndim=3] B_ = np.array(B, np.float32, order='F')
    cdef np.ndarray[float, ndim=3] L_ = np.array(L, np.float32, order='F')
    cdef np.ndarray[int, ndim=1] sL_ = np.array(B.shape, np.int32)
    cdef np.ndarray[float, ndim=1] S_ = np.array(vx_vol, np.float32)
    cdef int csf_ = int(csf)
    cdef int euclid_ = int(euclid)
    cdef int setnan_ = int(setnan)
    cdef int verb_ = int(verb)

    cdef np.ndarray[float, ndim=3] D = np.zeros_like(B_, order='F')
    cdef np.ndarray[np.uint32_t, ndim=3] I = np.zeros_like(B_, dtype=np.uint32,
                                                           order='F')

    vol_eidist(&D[0, 0, 0], &I[0, 0, 0],
               &B_[0, 0 ,0], &L_[0, 0 ,0],
               &sL_[0], &S_[0],
               csf_, euclid_, setnan_, verb_)
    
    return np.ascontiguousarray(D), np.ascontiguousarray(I)


def cat_vol_localstat(V, M, nb=1, stat=1, verb=False):
    ''' Interface to the the cat_vol_localstat function
    Original code by Robert Dahnke, Center of Neuroimaging University Jena

    Local Mean, Minimum, Maximum, SD, and Peak estimation
    ________________________________________________________________________
    Estimates specific functions volume V in a mask region M. For every 
    voxel vEM the values of the neigbors of v within a distance nb that 
    also belong to M were used. 
    
  
    S = cat_vol_localstat(V,M,nb,stat)
    
    V    (single)    input volume
    M    (logical)   mask volume
    nb   (double)    neigbhour distance (1 .. 10)
    stat (double)    1-mean, 2-min, 3-max, 4-std 
                     5-peak1, 6-peak2, 7-peak3
                     8-median
                     9-?
                     10-?
  
    ________________________________________________________________________
    Robert Dahnke 2012_01
    Center of Neuroimaging 
    University Jena
  
    $Id: cat_vol_localstat.c 1185 2017-09-13 15:02:47Z gaser $ 
    '''
    assert V.shape == M.shape, 'Volume and mask must have the same shape'
    if nb > 10:
        raise ValueError("number of neighbors is limited to 10. (Use reduce resolution instead.)")
    if stat < 1 or stat > 8:
        raise ValueError("fourth input has to be 1=mean, 2=min, 3=max, 4=std")

    cdef np.ndarray[float, ndim=3] V_ = np.array(V, np.float32, order='F')
    cdef np.ndarray[unsigned char, ndim=3] M_ = np.array(M, np.uint8, order='F')
    cdef np.ndarray[int, ndim=1] sL_ = np.array(V.shape, np.int32)
    cdef int nb_ = int(nb)
    cdef int stat_ = int(stat)
    cdef int verb_ = int(verb)

    cdef np.ndarray[float, ndim=3] Out1 = np.zeros_like(V_)
    cdef np.ndarray[float, ndim=3] Out2 = np.zeros_like(V_)
    cdef np.ndarray[float, ndim=3] Out3 = np.zeros_like(V_)

    vol_localstat(&Out1[0, 0, 0], &Out2[0, 0, 0], &Out3[0, 0, 0],
                  &V_[0, 0, 0], &M_[0, 0, 0], &sL_[0],
                  nb_, stat_, verb_)

    return np.ascontiguousarray(Out1), np.ascontiguousarray(Out2), np.ascontiguousarray(Out3)


def cat_vol_pbtp(SEG, WMD, CSFD):
    ''' Interface to the the cat_vol_pbtp function
    Original code by Robert Dahnke, Center of Neuroimaging University Jena


    quasi-euclidean distance calculation
    _____________________________________________________________________________
    [GMT,RPM,WMD,CSFD,II] = eqdist(SEG,WMD,CSFD[,opt])
    
    SEG  = (single) segment image with low and high boundary bd
    GMT  = (single) thickness image
    RPM  = (single) radial position map
    WMD  = (single) CSF distance map
    CSFD = (single) CSF distance map
    II  = (uint32) index of the inner (WM)  boundary voxel
    
    opt.bd   = (single) [low,high] boundary values (default 1.5 and 2.5)
    opt.CSFD = calculate CSFD
    opt.PVE  = use PVE informations (0=none,1=fast,2=exact)
    
    TODO:
     - eikonal distance for subsegmentation (region growing)
     - own labeling (
    _____________________________________________________________________________
    Robert Dahnke 2009_11
    Center of Neuroimaging 
    University Jena

    GBS: The "OPT" variables is hard-coded in the C code
    GBS: The "II" variables is NOT PART OF THE OUTPUT
    '''
    assert SEG.shape == WMD.shape
    assert SEG.shape == CSFD.shape


    cdef np.ndarray[float, ndim=3] SEG_ = np.array(SEG, np.float32, order='F')
    cdef np.ndarray[float, ndim=3] WMD_ = np.array(WMD, np.float32, order='F')
    cdef np.ndarray[float, ndim=3] CSFD_ = np.array(CSFD, np.float32, order='F')
    cdef np.ndarray[int, ndim=1] sL_ = np.array(SEG.shape, np.int32)

    cdef np.ndarray[float, ndim=3] GMT = np.zeros_like(SEG_, order='F')
    cdef np.ndarray[float, ndim=3] RPM = np.zeros_like(SEG_, order='F')


    vol_pbtp(&GMT[0, 0, 0], &RPM[0, 0, 0], &SEG_[0, 0, 0],
             &WMD_[0, 0, 0], &CSFD_[0, 0, 0], &sL_[0])
    
    return np.ascontiguousarray(GMT), np.ascontiguousarray(RPM)


def cat_vol_genus0(vol, th, any_genus=0):
    ''' Interface with the gnus0 code

    Based in the CAT12 code by Robert Dahnke, Center of Neuroimaging University Jena

    Some changes from the original:
    In genus0 fails, raises an exception
    argument


    Creates a closed surface from the isosurface in the "vol" with the threshold "th".
    '''
    assert vol.ndim == 3, 'Volume must be 3D'
    cdef np.ndarray[unsigned short, ndim=3] _input = np.array(vol > th, dtype=np.uint16, order='F')
    # use orientation from isosurface
    cdef float[16] ijk2ras = [0, 1, 0, 0,
                              1, 0, 0, 0,
                              0, 0, 1, 0,
                              0, 0, 0, 1]
    cdef genus0parameters[1] g0 # need an instance of genus0parameters
    cdef int i
    cdef int sz = int(_input.size)


    # Output
    cdef np.ndarray[unsigned short, ndim=1] M = np.zeros(sz, dtype=np.uint16)
    cdef np.ndarray[int, ndim=1] Tris
    cdef np.ndarray[double, ndim=1] Verts


    genus0init(g0) # initialize the instance, set default parameters */

    # set some parameters/options
    for i in range(3):
        g0.dims[i] = vol.shape[i]

    g0.return_adjusted_label_map = 1
    g0.cut_loops = 0
    g0.connectivity = 6
    g0.connected_component = 1
    g0.input = &_input[0, 0, 0]
    g0.value = 1
    g0.alt_value = 1
    g0.contour_value = 1
    g0.alt_contour_value = 1
    g0.biggest_component = 1
    for i in range(3):
        g0.pad[i] = 2
        g0.extraijkscale[i] = 1
    g0.ijk2ras = ijk2ras
    g0.verbose = 1
    g0.return_surface = 1
    g0.any_genus = int(any_genus)

    # in case of error raise RuntimeError
    if (genus0(g0)):
        genus0destruct(g0)
        raise RuntimeError('genus0 failed')

    g0.cut_loops = 1
    g0.connectivity = 18
    g0.value = 1
    g0.alt_value = 0
    g0.contour_value = 1
    g0.alt_contour_value = 0

    for i in range(sz):
      M[i] = g0.output[i]

    g0.input = &M[0]

    if (genus0(g0)):
        genus0destruct(g0)
        raise RuntimeError('genus0 failed')

    for i in range(sz):
      M[i] = g0.output[i];

    # return Tris and Verts
    Tris = np.zeros(3*g0.tri_count, dtype=np.int32)
    Verts = np.zeros(3*g0.vert_count, dtype=np.float)

    for i in range(3*g0.tri_count):
      Tris[i] = g0.triangles[i] + 1;

    for i in range(3*g0.vert_count):
      Verts[i] = g0.vertices[i] + 1.0;

    genus0destruct(g0)

    return (
        np.ascontiguousarray(M.reshape(vol.shape, order='F')),
        np.ascontiguousarray(Tris.reshape((-1, 3), order='F')),
        np.ascontiguousarray(Verts.reshape((-1, 3), order='F'))
    )


def cat_vbdist(P, R=None, S=(1, 1, 1)):
    ''' Interface to the the cat_vbdist function
    Original code by Robert Dahnke, Center of Neuroimaging University Jena

    voxelbased euclidean distance calculation
    ________________________________________________________________________
    Calculates the euclidean distance without PVE to an object in P with a 
    boundary of 0.5.

     [D,I,L] = vbdist(P[,R])
     
     P (single)  input image with zero for non elements and uint8 values for labels
     R (logical) range for distance calculation
     D (single)  distance image
     L (uint8)   label map
     I (uint32)  index of nearest point
    ________________________________________________________________________
    Robert Dahnke 2010_01
    Center of Neuroimaging 
    University Jena


    GBS: There is now an additional argument, S, for voxel sizes
    '''
    if R is None:
        R = np.ones_like(P, dtype=np.uint8)
    else:
        assert P.shape == R.shape

    cdef np.ndarray[float, ndim=3] P_ = np.array(P, np.float32, order='F')
    cdef np.ndarray[unsigned char, ndim=3] R_ = np.array(R, np.uint8, order='F')
    cdef np.ndarray[float, ndim=1] S_ = np.array(S, np.float32, order='F')
    cdef np.ndarray[int, ndim=1] sL_ = np.array(P.shape, np.int32)

    cdef np.ndarray[float, ndim=3] D = np.zeros_like(P_)
    cdef np.ndarray[unsigned int, ndim=3] I = np.zeros_like(P_, dtype=np.uint32)
    cdef np.ndarray[unsigned char, ndim=3] L = np.zeros_like(P_, dtype=np.uint8)

    vbdist(
        &D[0, 0, 0], &I[0,0, 0], &L[0, 0, 0],
        &P_[0, 0, 0], &R_[0, 0,0], &sL_[0], &S_[0]
    )

    return (
        np.ascontiguousarray(D),
        np.ascontiguousarray(I),
        np.ascontiguousarray(L)
    )
