#include "genus0.h"

#ifdef MATLAB_MEX_FILE
#include <mex.h> 
#endif

static int verbose,invconnectivity, connectivity, autocrop[3][2];
static int img_horiz,img_vert,img_depth,paddeddims[3];
static int nbrs[6],offs[27],nbrs18[18],nbrs26[26],pass[27];
static int elist18[27][19],elist[27][7],elist26[27][27];
static int *status,*cm,cm_size,que_size,*que,que_len,que_pos;
static int maxlevels,comp_count=10,cut_loops;
static int calloced=0,max_calloced=32;
static int *persist=NULL;

static void **calloc_list=NULL;
static unsigned char *zpic;
static float voxelsize[3],*fzpic,fzpicmax;

static mwSize *g_axis_len, *g_stride;
static mwSize paddeddims0[3];
static float *g_deltax, *g_tmp, *g_tmp_row, **g_j, *g_x, **g_recip, **g_square;
 
static void print_msg(char * msg)
  {
  printf("%s",msg); 
  }

static void *basic_calloc(mwSize nelem, mwSize elsize)
  {
  return(calloc(nelem,elsize));
  }

static void basic_free(void *ptr)
  {
  if (ptr!=NULL) free(ptr);
  }


static void Gfree(void *ptr)
  {
  int i,j;
  /* free something in the calloc_list */
  for(i=0;i<calloced;i++)
    if (calloc_list[i]==ptr)
      {
      basic_free(ptr);
      for(j=i+1;j<calloced;j++)
        {
        calloc_list[j-1]=calloc_list[j];
        persist[j-1]=persist[j];
        }
      calloced--;
      break;
      }
  }

static void Gfree_all(int keep_persist)
  {
  int i;
  if (calloc_list==NULL) return; /* already done freeing all */
  for(i=0;i<calloced;i++)
    if ((!keep_persist)||(persist[i]==0))
      {
      Gfree(calloc_list[i]);
      i--; /* don't advance i, since we compressed the list */
      }
  if (calloc_list!=NULL) basic_free(calloc_list);
  calloc_list=NULL;
  if (persist!=NULL) basic_free(persist);
  persist=NULL;
  }


static void error_msg(char * msg, int line)
  {
  char line_msg[100];
  if (calloc_list==NULL) return; /* must have been an error already */
  print_msg(msg);
  sprintf(line_msg,"Line: %d.\n",line);
  print_msg(line_msg);
  Gfree_all(0);
  }

static void *Gcalloc(mwSize nelem, mwSize elsize, int make_persist)
  {
  void **cl;
  int *p,j;

  if (calloced==max_calloced) /* we need more room */
    {
    cl=calloc_list;
    calloc_list=(void**)basic_calloc(max_calloced*2,sizeof(void*));
    if (calloc_list==NULL)
      {
      calloc_list=cl;
      error_msg("Memory error.\n",__LINE__);
      return(NULL);
      }    
    for(j=0;j<calloced;j++) calloc_list[j]=cl[j];
    basic_free(cl);

    p=persist;
    persist=(int*)basic_calloc(max_calloced*2,sizeof(int));
    if (persist==NULL) 
      {
      persist=p;
      error_msg("Memory error.\n",__LINE__);
      return(NULL);
      }
    for(j=0;j<calloced;j++) persist[j]=p[j];
    basic_free(p);
    max_calloced*=2;
    }

  calloc_list[calloced]=basic_calloc(nelem,elsize);
  if (calloc_list[calloced]==NULL) 
    {
    error_msg("Memory error.\n",__LINE__);
    return(NULL);
    }
  else persist[calloced]=make_persist;
  return(calloc_list[calloced++]);
  }
 
 
static void calc_elist(void)
  {
  int vv,i,j,k,h,i1,j1,k1,count,thecase;
  int cases[6][3]={{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
  int cases1826[26][3];

  /* calculate number of neighbors for each point in 3x3x3 cube */               

  h=0;
  for (k=0;k<3;k++)  /* for each spot in 3x3x3 cube */
    for(j=0;j<3;j++)
      for(i=0;i<3;i++)
        {
        count=0; /* count of nbrs */
        for (thecase=0;thecase<6;thecase++)
          {
          i1=i+cases[thecase][0];
          j1=j+cases[thecase][1];
          k1=k+cases[thecase][2];
          if (i1>=0 && i1<=2 && j1>=0 && j1<=2 && k1>=0 && k1<=2)
            { /* if in 3x3x3 cube */
            vv=i1+j1*3+k1*9;
            if (vv!=13) elist[h][++count]=vv;
            }
          }
        elist[h][0]=count;                                                              
        h++;
        }                       

  h=0;
  for (k=0;k<3;k++)  /* for each spot in 3x3x3 cube */
    for(j=0;j<3;j++)
      for(i=0;i<3;i++)
        if (!((i==0 || i==2)&&(j==0 || j==2)&&(k==0 || k==2)))
        if ((i!=1)||(j!=1)||(k!=1))
          {
          cases1826[h][0]=i-1;
          cases1826[h][1]=j-1;
          cases1826[h][2]=k-1;
          h++;
          }
  h=0;
  for (k=0;k<3;k++)  /* for each spot in 3x3x3 cube */
    for(j=0;j<3;j++)
      for(i=0;i<3;i++)
        {
        count=0; /* count of nbrs */
        for (thecase=0;thecase<18;thecase++)
          {
          i1=i+cases1826[thecase][0];
          j1=j+cases1826[thecase][1];
          k1=k+cases1826[thecase][2];
          if (i1>=0 && i1<=2 && j1>=0 && j1<=2 && k1>=0 && k1<=2)
            { /* if in 3x3x3 cube */
            vv=i1+j1*3+k1*9;
            if (vv!=13) elist18[h][++count]=vv;
            }
          }
        elist18[h][0]=count;                                                              
        h++;
        }                       

  h=0;
  for (k=0;k<3;k++)  /* for each spot in 3x3x3 cube */
    for(j=0;j<3;j++)
      for(i=0;i<3;i++)
        if ((i!=1)||(j!=1)||(k!=1))
          {
          cases1826[h][0]=i-1;
          cases1826[h][1]=j-1;
          cases1826[h][2]=k-1;
          h++;
          }
  h=0;
  for (k=0;k<3;k++)  /* for each spot in 3x3x3 cube */
    for(j=0;j<3;j++)
      for(i=0;i<3;i++)
        {
        count=0; /* count of nbrs */
        for (thecase=0;thecase<26;thecase++)
          {
          i1=i+cases1826[thecase][0];
          j1=j+cases1826[thecase][1];
          k1=k+cases1826[thecase][2];
          if (i1>=0 && i1<=2 && j1>=0 && j1<=2 && k1>=0 && k1<=2)
            { /* if in 3x3x3 cube */
            vv=i1+j1*3+k1*9;
            if (vv!=13) elist26[h][++count]=vv;
            }
          }
        elist26[h][0]=count;                                                              
        h++;
        } 

  }

 
 
static void process_row(int axis,float *start)
  {
  register mwSize len;
  register float *p, *p2, *p3, *p_end, pv, p2v, **jp, **j2p, **j_end, *x2p;
  register float x, x0, x2, *recip, *square, dx;

  dx = g_deltax[axis];
  if (dx < 0.f) dx = -dx;
  j_end = g_j;
  p = start;
  p_end = start + g_axis_len[axis];
  do {
    if (*p != FLT_MAX) *j_end++ = p;
  } while (++p != p_end);
  if (j_end == g_j) return;
  jp = j2p = g_j;
  x2p = g_x;
  *g_x = -FLT_MAX;
  if (++jp != j_end) {
    pv = *(p = *jp);
    p2v = *(p2 = *j2p);
    x2 = -FLT_MAX;
    x0 = dx*(p - start);
    square = g_square[axis];
    recip = g_recip[axis];
    while (1) {
      len = p - p2;
      x = x0 + (pv - p2v - square[len])*recip[len];
      if (x > x2) {
        *++j2p = p;
        *++x2p = x;
        if (++jp == j_end) break;
        p2 = p;
        p2v = pv;
        x2 = x;
        pv = *(p = *jp);
        x0 = dx*(p - start);
      } else {
        p2v = *(p2 = *--j2p);
        x2 = *--x2p;
      }
    }
  }

  len = g_axis_len[axis];
  p = p_end = g_tmp_row + len;
  p3 = start + len;
  while (len--) {
    x = dx*len;
    while (*x2p > x) {
      j2p--;
      x2p--;
    }
    p2 = *j2p;
    x = dx*(--p3 - p2);
    *--p = *p2 + x*x;
  }
  p2 = start;
  while (p != p_end) *p2++ = *p++;
  return;
  }

static void recursive_add_dist_squared(int axis,float *start)
  {
  mwSize len;
  float *p, *p_end, *p2, *p2_end, *p3;

  if (axis == 0) {
    process_row(0, start);
    return;
  }

  len = g_stride[axis];
  p = start;
  p_end = p + len;
  p2_end = g_tmp + g_axis_len[axis];
  while (p != p_end) {
    p2 = g_tmp;
    p3 = p;
    while (p2 != p2_end) {
      *p2++ = *p3;
      p3 += len;
    }
    process_row(axis, g_tmp);
    p2 = g_tmp;
    p3 = p++;
    while (p2 != p2_end) {
      *p3 = *p2++;
      p3 += len;
    }
  }

  p = start;
  len = g_axis_len[axis];
  while (len--) {
    recursive_add_dist_squared(axis-1, p);
    p += g_stride[axis];
  }
  return;
  }



static int dist_squared(
     int rank,
     mwSize *axis_len,
     float *deltax,
     register char *inimage,
     register char inobject,
     float *outdist_squared)
  {
  mwSize max_axis_len, data_len;
  register mwSize len;
  register float *p, *p2, ftmp, ftmp2;
  int i;

  if(!(rank >= 0)) 
    {error_msg("Error during distance transform.\n",__LINE__);return(1);}
  if (rank == 0) {
    *outdist_squared = (*inimage == inobject)? 0.f : FLT_MAX;
    return 0;
  }

  g_stride = (mwSize *)Gcalloc(rank,sizeof(mwSize),0);
  if (g_stride==NULL) {error_msg("Memory error.\n",__LINE__);return(1);}

  max_axis_len = 2;
  data_len = 1;
  for (i=0; i<rank; i++) {
    if (!(axis_len[i] > 1)) return(1);
    if (!(deltax[i] != 0.f)) return(1);
    if (max_axis_len < axis_len[i]) max_axis_len = axis_len[i];
    g_stride[i] = data_len;
    data_len *= axis_len[i];
  }

  g_tmp = (float *)Gcalloc(max_axis_len,sizeof(float),0);
  if (g_tmp==NULL) {error_msg("Memory error.\n",__LINE__);return(1);}

  g_tmp_row = (float *)Gcalloc(max_axis_len,sizeof(float),0);
  if (g_tmp_row==NULL) {error_msg("Memory error.\n",__LINE__);return(1);}

  g_j = (float **)Gcalloc(max_axis_len,sizeof(float *),0);
  if (g_j==NULL) {error_msg("Memory error.\n",__LINE__);return(1);}

  g_x = (float *)Gcalloc(max_axis_len,sizeof(float),0);
  if (g_x==NULL) {error_msg("Memory error.\n",__LINE__);return(1);}

  p = outdist_squared;
  len = g_stride[rank - 1];
  while (len--) {
    *p++ = (*inimage++ == inobject)? 0.f : FLT_MAX;
  }
  p2 = outdist_squared;
  ftmp = deltax[rank - 1];
  if (ftmp < 0.f) ftmp = -ftmp;
  len = data_len - g_stride[rank - 1];
  while (len--) {
    if (*inimage++ == inobject) *p = 0.f;
    else if ((*p = *p2) != FLT_MAX) *p += ftmp;
    p++;
    p2++;
  }
  if (rank > 1) {
    g_axis_len = axis_len;
    g_deltax = deltax;

    g_recip = (float **)Gcalloc(rank - 1,sizeof(float *),0);
    if (g_recip==NULL) {error_msg("Memory error.\n",__LINE__);return(1);}
    g_square = (float **)Gcalloc(rank - 1,sizeof(float *),0);
    if (g_square==NULL) {error_msg("Memory error.\n",__LINE__);return(1);}

    for (i=0; i<rank-1; i++) {g_recip[i]=g_square[i]=NULL;}

    for (i=0; i<rank-1; i++) {

      g_recip[i] = (float *)Gcalloc(axis_len[i] + 1,sizeof(float),0);
      if (g_recip[i]==NULL) {error_msg("Memory error.\n",__LINE__);return(1);}
      g_square[i] = (float *)Gcalloc(axis_len[i] + 1,sizeof(float),0);
      if (g_square[i]==NULL) {error_msg("Memory error.\n",__LINE__);return(1);}

      len = axis_len[i];
      do {
        ftmp2 = len*deltax[i];
        if (ftmp2 < 0.f) ftmp2 = -ftmp2;
        g_recip[i][len] = 0.5f/ftmp2;
        g_square[i][len] = ftmp2*ftmp2;
      } while (--len);
    }
    while (p2 != outdist_squared) {
      len = g_stride[rank - 1];
      while (len--) {
        --p2;
        if (*--p != FLT_MAX) {
          ftmp2 = *p + ftmp;
          *p *= *p;
          if (*p2 > ftmp2) *p2 = ftmp2;
        }
      }
      recursive_add_dist_squared(rank-2, p);
    }
    len = g_stride[rank - 1];
    while (len--) {
      if (*--p != FLT_MAX) *p *= *p;
    }
    recursive_add_dist_squared(rank-2, p);
    for (i=0; i<rank-1; i++) {
      Gfree(g_recip[i]);
      Gfree(g_square[i]);
    }
    Gfree(g_recip);
    Gfree(g_square);
  } else {
    len = data_len - 1;
    while (len--) {
      --p2;
      if (*--p != FLT_MAX) {
        ftmp2 = *p + ftmp;
        *p *= *p;
        if (*p2 > ftmp2) *p2 = ftmp2;
      }
    }
    if (*--p != FLT_MAX) *p *= *p;
  }
  Gfree(g_stride);
  Gfree(g_tmp);
  Gfree(g_tmp_row);
  Gfree(g_j);
  Gfree(g_x);
  return 0;
  }
 
 
static int sub2ind(int horiz, int vert, int depth, int img_horiz, int img_vert)
  {  /* convert subscripts to linear index */
  return(horiz+(vert+depth*img_vert)*img_horiz); /* factoring to save a mult... why not? */
  }


static int get_cc(unsigned char * zpic, int *que, int * status, int *dims, int *ac, int connectivity)
  {
  int i,*offs,nbrs6[6],nbrs18[18],*g_counts;
  int j,k,h,vox,c,w,groups,que_pos,que_len;  
  int autocrop[3][2],woffsh,w0;

  vox=dims[0]*dims[1]*dims[2];

  nbrs6[0]=-1;             nbrs6[1]=1;
  nbrs6[2]=-dims[0];       nbrs6[3]=dims[0];
  nbrs6[4]=-dims[0]*dims[1]; nbrs6[5]=-nbrs6[4];

  for (k=h=0;k<3;k++)
    for(j=0;j<3;j++)
      for(i=0;i<3;i++)
        if (!((i==0 || i==2)&&(j==0 || j==2)&&(k==0 || k==2)))
          if ((i!=1)||(j!=1)||(k!=1))  /* calculate offsets for all 18 neighboring voxels */
            nbrs18[h++]=sub2ind(i,j,k,dims[0],dims[1])-sub2ind(1,1,1,dims[0],dims[1]);

  if (connectivity==18) {c=18; offs=nbrs18;}
  else {c=6; offs=nbrs6;}

  for(i=0;i<vox;i++) {status[i]=0; zpic[i]=(zpic[i]>0);}

  w=groups=0;
  for(k=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++)
      for(i=0;i<dims[0];i++)
      {
      if (zpic[w]==1) /* not processed yet */
        {
        groups++;  /* new group */
        status[w]=groups;
        zpic[w]=2;  /* 0=background, 1=not on que, 2=on que or done with */
        que_len=0;
        for(h=0;h<c;h++)  /* add nbrs to que */
          if (zpic[woffsh=w+offs[h]]==1) zpic[que[que_len++]=woffsh]=2;

        que_pos=0;
        while (que_pos<que_len)
          {
          status[w0=que[que_pos++]]=groups;  /* get a voxel from que */
          for(h=0;h<c;h++)  /* add nbrs to que */
            if (zpic[woffsh=w0+offs[h]]==1) zpic[que[que_len++]=woffsh]=2;
          }
        }
      w++;
      } /* for loop through voxels */

  g_counts=(int *)Gcalloc(groups+1,sizeof(int),0);
  if (g_counts==NULL) {error_msg("Memory error.\n",__LINE__);return(1);}

  for(i=0;i<vox;i++) g_counts[status[i]]++;
  w=h=0;
  for (i=1;i<=groups;i++) if (g_counts[i]>w) {w=g_counts[i]; h=i;}
  for(i=0;i<vox;i++) {zpic[i]=(status[i]==h); status[i]=0;}

  for(k=0;k<3;k++) {autocrop[k][0]=dims[k]; autocrop[k][1]=0;}
  for(k=h=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++)
      for(i=0;i<dims[0];i++)
        {
        if (zpic[h])
          {
          if (autocrop[0][0]>i) autocrop[0][0]=i; if (autocrop[0][1]<i) autocrop[0][1]=i;
          if (autocrop[1][0]>j) autocrop[1][0]=j; if (autocrop[1][1]<j) autocrop[1][1]=j;
          if (autocrop[2][0]>k) autocrop[2][0]=k; if (autocrop[2][1]<k) autocrop[2][1]=k;
          }
        h++;
        }
  h=0;
  for(i=0;i<3;i++)
    for(j=0;j<2;j++)
      ac[h++]=autocrop[i][j];

  Gfree(g_counts);
  return(0);
  }


static int set_up(genus0parameters * g0)
  {
  int i,j,h,k,totlen,hz,*pad,ac[6];
  unsigned short value, * input;
  float *m,minvoxelsize;

  /* set some global variables */
  verbose=(g0->verbose!=0);
  if (verbose) print_msg("Setting up...\n");

  if((input=g0->input)==NULL)
    {error_msg("No input volume.\n",__LINE__);return(1);}    

  if ((g0->dims[0]<=0)||(g0->dims[1]<=0)||(g0->dims[2]<=0))
    {error_msg("Bad input volume dimensions.\n",__LINE__);return(1);}    

  calloc_list=(void**)basic_calloc(max_calloced,sizeof(void*));
  if (calloc_list==NULL)
    {error_msg("Memory error.\n",__LINE__);return(1);}

  persist=(int*)basic_calloc(max_calloced,sizeof(int));
  if (persist==NULL)
    {error_msg("Memory error.\n",__LINE__);return(1);}

  if ((g0->connectivity!=6)&&(g0->connectivity!=18)) g0->connectivity=6;
  connectivity=g0->connectivity;
  invconnectivity=(g0->connectivity==6)?18:6;

  pad=g0->pad;
  for(i=0;i<3;i++) if (pad[i]<2) pad[i]=2;

  g0->biggest_component=(g0->biggest_component!=0);
  cut_loops=g0->cut_loops=(g0->cut_loops!=0);
  g0->return_surface=(g0->return_surface!=0);
  g0->return_adjusted_label_map=(g0->return_adjusted_label_map!=0);

  value=g0->value; 

  /* get cropping limits */
  for(k=0;k<3;k++) {autocrop[k][0]=(g0->dims)[k]; autocrop[k][1]=0;}
  for(k=h=0;k<(g0->dims[2]);k++)
    for(j=0;j<(g0->dims[1]);j++)
      for(i=0;i<(g0->dims[0]);i++)
        {
        if (input[h]>=value)
          {
          if (autocrop[0][0]>i) autocrop[0][0]=i; if (autocrop[0][1]<i) autocrop[0][1]=i;
          if (autocrop[1][0]>j) autocrop[1][0]=j; if (autocrop[1][1]<j) autocrop[1][1]=j;
          if (autocrop[2][0]>k) autocrop[2][0]=k; if (autocrop[2][1]<k) autocrop[2][1]=k;
          }
        h++;
        }

  if ((autocrop[0][0]>autocrop[0][1]) || (autocrop[1][0]>autocrop[1][1]) ||
      (autocrop[2][0]>autocrop[2][1])) 
    {error_msg("No data in volume matches specified value.\n",__LINE__);return(1);}

  /* calculate cropped dimensions, and total length */
  img_horiz= autocrop[0][1]-autocrop[0][0]+1+pad[0]*2;
  img_vert = autocrop[1][1]-autocrop[1][0]+1+pad[1]*2;
  img_depth= autocrop[2][1]-autocrop[2][0]+1+pad[2]*2;
  totlen=img_horiz*img_vert*img_depth;
  paddeddims[0]=img_horiz; paddeddims[1]=img_vert; paddeddims[2]=img_depth;

  if ((zpic=(unsigned char *)Gcalloc(totlen,sizeof(unsigned char),0))==NULL)
    {error_msg("Memory error.\n",__LINE__);return(1);}

  for(k=0;k<(img_depth-pad[2]*2);k++)
    for(j=0;j<(img_vert-pad[1]*2);j++)
      for(i=0;i<(img_horiz-pad[0]*2);i++)
        {
        h=sub2ind(pad[0]+i,pad[1]+j,pad[2]+k,img_horiz,img_vert);
        hz=sub2ind(autocrop[0][0]+i,autocrop[1][0]+j,autocrop[2][0]+k,g0->dims[0],g0->dims[1]);
        zpic[h]=(int)(input[hz]>=value);
        }

  que=NULL;
  status=NULL;

  if (g0->connected_component)
    {
    que_size=totlen;
    if ((que=(int *)Gcalloc(que_size,sizeof(int),0))==NULL) 
      {error_msg("Memory error.\n",__LINE__);return(1);}

    /* initialize voxel status */
    if ((status=(int *)Gcalloc(totlen,sizeof(int),0))==NULL) 
      {error_msg("Memory error.\n",__LINE__);return(1);}

    /* zpic is binary. return connected component */
    if(get_cc(zpic,que,status,paddeddims,ac,g0->connectivity))
      {error_msg("Connected component error.\n",__LINE__);return(1);}

    /* now we can crop zpic even more */
    autocrop[0][0]+=ac[0]-pad[0]; 
    autocrop[0][1]=autocrop[0][0]+(ac[1]-ac[0]); 

    autocrop[1][0]+=ac[2]-pad[1]; 
    autocrop[1][1]=autocrop[1][0]+(ac[3]-ac[2]);

    autocrop[2][0]+=ac[4]-pad[2]; 
    autocrop[2][1]=autocrop[2][0]+(ac[5]-ac[4]); 

    img_horiz= autocrop[0][1]-autocrop[0][0]+1+pad[0]*2;
    img_vert = autocrop[1][1]-autocrop[1][0]+1+pad[1]*2;
    img_depth= autocrop[2][1]-autocrop[2][0]+1+pad[2]*2;
    totlen=img_horiz*img_vert*img_depth;

    h=0;
    for(k=0;k<img_depth;k++)
      for(j=0;j<img_vert;j++)
        for(i=0;i<img_horiz;i++)
          {
          hz=sub2ind(ac[0]-pad[0]+i,ac[2]-pad[1]+j,ac[4]-pad[2]+k,paddeddims[0],paddeddims[1]);
          zpic[h++]=zpic[hz];
          }
    paddeddims[0]=img_horiz;  paddeddims[1]=img_vert; paddeddims[2]=img_depth;
    }

  if ((m=g0->ijk2ras)==NULL) 
    voxelsize[0]=voxelsize[1]=voxelsize[2]=1.0;
  else
    for(i=0;i<3;i++) 
      voxelsize[i]=sqrt(m[0+i]*m[0+i]+m[4+i]*m[4+i]+m[8+i]*m[8+i]);
  for(i=0;i<3;i++) if ((g0->extraijkscale)[i]<=0.0) (g0->extraijkscale)[i]=1.0;
  for(i=0;i<3;i++) voxelsize[i]*=(g0->extraijkscale)[i];
  minvoxelsize=voxelsize[0];
  if (minvoxelsize>voxelsize[1]) minvoxelsize=voxelsize[1];
  if (minvoxelsize>voxelsize[2]) minvoxelsize=voxelsize[2];
  if (minvoxelsize<=0.0) minvoxelsize=1.0;

  /* calculate a bunch of offsets to face-sharing voxel neighbors */
  nbrs[0]=-1;                    nbrs[1]=-nbrs[0];
  nbrs[2]=-img_horiz;            nbrs[3]=-nbrs[2];
  nbrs[4]=-(img_horiz*img_vert); nbrs[5]=-nbrs[4];

  for (k=h=0;k<3;k++)
    for(j=0;j<3;j++)
      for(i=0;i<3;i++)   /* calculate offsets for all 27 voxels in neighborhood */
        offs[h++]=sub2ind(i,j,k,img_horiz,img_vert)-sub2ind(1,1,1,img_horiz,img_vert);

  for (k=h=0;k<3;k++)
    for(j=0;j<3;j++)
      for(i=0;i<3;i++)
        if (!((i==0 || i==2)&&(j==0 || j==2)&&(k==0 || k==2)))
          if ((i!=1)||(j!=1)||(k!=1))  /* calculate offsets for all 18 neighboring voxels */
            nbrs18[h++]=sub2ind(i,j,k,img_horiz,img_vert)-sub2ind(1,1,1,img_horiz,img_vert);

  for (k=h=0;k<3;k++)
    for(j=0;j<3;j++)
      for(i=0;i<3;i++)
         if ((i!=1)||(j!=1)||(k!=1)) /* calculate offsets for 26 neighboring voxels */
          nbrs26[h++]=sub2ind(i,j,k,img_horiz,img_vert)-sub2ind(1,1,1,img_horiz,img_vert);


  /* even more offsets */
  calc_elist();
  for(i=0;i<27;i++) pass[i]=1;

  if (!(g0->any_genus))
    {
    que_size=totlen;
    if (que==NULL) /* might have already calloced it above ... */
      if ((que=(int *)Gcalloc(que_size,sizeof(int),0))==NULL) 
        {error_msg("Memory error.\n",__LINE__);return(1);}

    if ((fzpic=(float *)Gcalloc(totlen,sizeof(float),0))==NULL)
      {error_msg("Memory error.\n",__LINE__);return(1);}

    paddeddims[0]=img_horiz; paddeddims[1]=img_vert; paddeddims[2]=img_depth;

    paddeddims0[0]=img_horiz; paddeddims0[1]=img_vert; paddeddims0[2]=img_depth;
    /* calculate distance transform */
    /* sets fzpic to dist from {zpic=(1-cut_loops)} */
    if (dist_squared(3,paddeddims0,voxelsize,(char*)zpic,(char)(1-g0->cut_loops),fzpic))
      {error_msg("Memory error.\n",__LINE__);return(1);}
      
    fzpicmax=0.0;
    for(i=0;i<totlen;i++)
      {
      fzpic[i]=sqrt(fzpic[i]);
      if (fzpic[i]>fzpicmax) fzpicmax=fzpic[i];
      }

    maxlevels=(int)(fzpicmax/minvoxelsize+2.5);

    /* Formula for cm_size.  certain to be <= totlen + 12 */
    cm_size=totlen+12;
    if ((cm=(int *)Gcalloc(cm_size,sizeof(int),0))==NULL) 
      {error_msg("Memory error.\n",__LINE__);return(1);}

    /* initialize voxel status */
    if (status==NULL) /* might have already calloced it above ... */
      if ((status=(int *)Gcalloc(totlen,sizeof(int),0))==NULL) 
        {error_msg("Memory error.\n",__LINE__);return(1);}

    for (i=0;i<totlen;i++) 
      status[i]=(zpic[i]==g0->cut_loops)?0:3; /* 3 if not wanted component */

    /* set status of double boundary to 4 */
    for(i=0;i<img_horiz;i++)
      for(j=0;j<img_vert;j++)
        {
        status[sub2ind(i,j,0,img_horiz,img_vert)]=4;
        status[sub2ind(i,j,img_depth-1,img_horiz,img_vert)]=4;
        }
    for(i=0;i<img_horiz;i++)
      for(j=0;j<img_depth;j++)
        {
        status[sub2ind(i,0,j,img_horiz,img_vert)]=4;
        status[sub2ind(i,img_vert-1,j,img_horiz,img_vert)]=4;
        }
    for(i=0;i<img_vert;i++)
      for(j=0;j<img_depth;j++)
        {
        status[sub2ind(0,i,j,img_horiz,img_vert)]=4;
        status[sub2ind(img_horiz-1,i,j,img_horiz,img_vert)]=4;
        }

    /* set fzpic of single boundary to fzpicmax */
    for(i=1;i<img_horiz-1;i++)
      for(j=1;j<img_vert-1;j++)
        {
        fzpic[sub2ind(i,j,1,img_horiz,img_vert)]=fzpicmax;
        fzpic[sub2ind(i,j,img_depth-2,img_horiz,img_vert)]=fzpicmax;
        }
    for(i=1;i<img_horiz-1;i++)
      for(j=1;j<img_depth-1;j++)
        {
        fzpic[sub2ind(i,1,j,img_horiz,img_vert)]=fzpicmax;
        fzpic[sub2ind(i,img_vert-2,j,img_horiz,img_vert)]=fzpicmax;
        }
    for(i=1;i<img_vert-1;i++)
      for(j=1;j<img_depth-1;j++)
        {
        fzpic[sub2ind(1,i,j,img_horiz,img_vert)]=fzpicmax;
        fzpic[sub2ind(img_horiz-2,i,j,img_horiz,img_vert)]=fzpicmax;
        }

    } /* if (!(g0->any_genus)) */
  else
    {
    if (que!=NULL) Gfree(que); /* don't need the que if not doing topology correction */
    }
  return(0); /* no error */
  }
 
 
static int truecm(int st) /* return true component, set cm[st] to true comp */
  {
  int s0,s1;
  if (cm[st]!=st)
    {
    s0=st;
    while (cm[st]!=st) st=cm[st];
    while (cm[s0]!=st) {s1=cm[s0]; cm[s0]=st; s0=s1;}
    }
  return(st);
  }


static int truecmvx(int vx) /* return true component of voxel */
  {
  return(status[vx]=truecm(status[vx]));
  }


static int test18(int qqp, int *nc)
  {
  int elQqpj,stqqpn,i,j,ec_count=0,ec[27],found_another=1;
  int st[19],Que_len,Que_pos,Que[27],Qqp;
  static int idx[18]={1,3,4,5,7,9,10,11,12,14,15,16,17,19,21,22,23,25};

  for(i=0;i<27;i++) ec[i]=0;

  for(i=0;i<18;i++) /* for each 18 neighbor */
    if (status[qqp+nbrs18[i]]>10) /* if nbr is in an Mcubes component */
      if (!ec[idx[i]]) /* if not assigned an edge component */
        {
        /* create new edge component */
        stqqpn=truecmvx(qqp+nbrs18[i]);
        st[++ec_count]=stqqpn; /* remember the status associated with the new component */

        /* find all in the edge component */

        Que_pos=0; Que_len=1; Que[0]=idx[i]; /* add it to Que */
        ec[idx[i]]=ec_count;
        while (Que_pos<Que_len)
          {
          Qqp=Que[Que_pos];
          for (j=1;j<=elist18[Qqp][0];j++) /* for each nbr of Qqp */
           {
           elQqpj=elist18[Qqp][j];
           if ((status[qqp+offs[elQqpj]]>10) && (!ec[elQqpj]))
              ec[Que[Que_len++]=elQqpj]=ec_count; /* add nbr to que */
           }
          Que_pos++;
          }

        /* compare with previous components */
        for(j=1;j<ec_count;j++)
          {
          if (st[j]==stqqpn) /* bad news if components are the same */
            {
            found_another=0;
            i=18;
            break;
            }
          }
        }

  *nc=(ec_count>0)?(st[1]):0;
  return(found_another);
  }


static int test6(int qqp, int *nc)
  {
  int elQqpj,stqqpn,i,j,ec_count=0,ec[27],found_another=1;
  int st[7],Que_len,Que_pos,Que[27],Qqp;
  static int idx[6]={12,14,10,16,4,22};

  for(i=0;i<27;i++) ec[i]=0;

  pass[1]=pass[3]=pass[5]=pass[7]=pass[9]=pass[11]=pass[15]=pass[17]=
  pass[19]=pass[21]=pass[23]=pass[25]=0;
  if(status[qqp+offs[ 4]]>10) pass[ 1]=pass[ 3]=pass[ 5]=pass[ 7]=1;
  if(status[qqp+offs[10]]>10) pass[ 1]=pass[ 9]=pass[11]=pass[19]=1;
  if(status[qqp+offs[12]]>10) pass[ 3]=pass[ 9]=pass[15]=pass[21]=1;
  if(status[qqp+offs[14]]>10) pass[ 5]=pass[11]=pass[17]=pass[23]=1;
  if(status[qqp+offs[16]]>10) pass[ 7]=pass[15]=pass[17]=pass[25]=1;
  if(status[qqp+offs[22]]>10) pass[19]=pass[21]=pass[23]=pass[25]=1;

  for(i=0;i<6;i++) /* for each neighbor */
    if (status[qqp+nbrs[i]]>10) /* if nbr is in an Mcubes component */
      if (!ec[idx[i]]) /* if not assigned an edge component */
        {
        /* create new edge component */

        stqqpn=truecmvx(qqp+nbrs[i]);
        st[++ec_count]=stqqpn; /* remember the status associated with the new component */

        /* find all in the edge component */

        Que_pos=0; Que_len=1; Que[0]=idx[i]; /* add it to Que */
        ec[idx[i]]=ec_count;
        while (Que_pos<Que_len)
          {
          Qqp=Que[Que_pos];
          for (j=1;j<=elist[Qqp][0];j++) /* for each nbr of Qqp */
            {
            elQqpj=elist[Qqp][j];
            if ((status[qqp+offs[elQqpj]]>10) && (!ec[elQqpj]) && pass[elQqpj])
              ec[Que[Que_len++]=elQqpj]=ec_count; /* add nbr to que */
            }
          Que_pos++;
          }

        /* compare with previous components */
        for(j=1;j<ec_count;j++)
          {
          if (st[j]==stqqpn) /* bad news if components are the same */
            {
            found_another=0;
            i=6;
            break;
            }
          }
        }

  *nc=(ec_count>0)?(st[1]):0;
  return(found_another);
  }
 
 
static int cmtostat(void)
  {
  int j,i,*cmremap=NULL,*ccount=NULL,totlen;
  char msg[200];

  if ((cmremap=(int *)Gcalloc(comp_count+1,sizeof(int),0))==NULL)
    {error_msg("Memory error.\n",__LINE__);return(1);}
  if ((ccount=(int *)Gcalloc(comp_count+1,sizeof(int),0))==NULL)
    {error_msg("Memory error.\n",__LINE__);return(1);}

  for(i=11;i<=comp_count;i++) ccount[truecm(i)]++;

  j=10;
  for(i=11;i<=comp_count;i++)
    if (ccount[i])
      {
      j++;
      cmremap[i]=j;
      }
  totlen=img_horiz*img_vert*img_depth;
  for (i=0;i<totlen;i++) if (status[i]>10) status[i]=cmremap[cm[status[i]]];

  if (verbose&&0)
    {
    sprintf(msg,"Components reduced from %u to %u.\n",comp_count-10,j-10);
    print_msg(msg);
    }
  comp_count=j;
  for(i=11;i<=comp_count;i++) cm[i]=i;
  Gfree(cmremap);
  Gfree(ccount);
  return(0);
  }
 
 
static void find_component(int level)
  {
  int i,qqp,qqpni,vox,nc;
  int found_another,*nbrs0,totlen;
  int (*test)(int, int *);
  float flevel;
  int theconnectivity;  

  if (cut_loops) theconnectivity=connectivity;
  else theconnectivity=invconnectivity;

  totlen=img_horiz*img_vert*img_depth;
  nbrs0=nbrs; test=test6;
  if (theconnectivity==18) {nbrs0=nbrs18; test=test18;}

  flevel=(level-1.0)/(maxlevels-1.0)*fzpicmax;

  for(vox=0;vox<totlen;vox++)
    {

    /* if level good and voxel unprocessed */
    if ((fzpic[vox]>=flevel) && (status[vox]==0)) 
      { 
      que_len=que_pos=0;
      que[que_len]=vox; /* add it to que */
      que_len++; if (que_len==que_size) que_len=0;
      status[vox]=2; /* mark as on que */
      while (que_pos!=que_len)
        {
        qqp=que[que_pos];
        /* check if can add */
        nc=0;
        found_another=test(qqp,&nc);
        /* if you can, add it, and combine components if needed */
        if (found_another)
          {
          if(nc==0) /* if no neighboring component */
            {
            comp_count++;
            cm[comp_count]=comp_count;
            status[qqp]=comp_count;
            }
          else
            {
            status[qqp]=nc; /* there was this true component nc */
            for(i=0;i<theconnectivity;i++)
              {
              qqpni=qqp+nbrs0[i]; /* for each nbr */
              if (status[qqpni]>10) /* if part of a component */
                if (truecmvx(qqpni)!=nc) /* if not the same component as qqp */
                  cm[status[qqpni]]=status[qqpni]=nc;
              }
            }
          /* add nbrs to que, if level ok, etc. */
          for(i=0;i<theconnectivity;i++)
            {
            qqpni=qqp+nbrs0[i];
            if ((fzpic[qqpni]>=flevel) && (status[qqpni]==0)) 
              {  
              que[que_len]=qqpni; /* add it to que */
              que_len++; if (que_len==que_size) que_len=0;
              status[qqpni]=2; /* mark as on que */
              }
            } /* for i */
          } /* if found another */
        else
          {
          status[qqp]=0;
          }
        /* move to next point in que */
        que_pos++; if (que_pos==que_size) que_pos=0;
        } /* while que not empty */
      } /* end if >=level and status==0 */
    } /* end for vox */
  } /* end find_component */
 

static int GetSurf(unsigned char * J, unsigned char val, int * dims, int connectivity,
  int ** Tris, float ** Verts, int * Tri_count, int * Vert_count, genus0parameters *g0)
  {
  unsigned char *status,*cidx,*pidx;
  int img_vert,iv1,img_horiz,img_depth,cellplane,cellplanem[3];
  int i,j,k,u,v,w,*v_idx,tri_count,vert_count,ic,cells,tc[256],vc[256];
  int ih1,id1,stat,pc,vert_count1,tri_count1,tcm[3],tsvw;
  int offs[8],coffs[12],ctype[12],*ts,hep[3],planeidx[12],planeoff[2];
  int hep2[3]={5,6,11},*tris,vert_count_times2;
  unsigned char pows[8]={64, 128, 16, 32, 4, 8, 1, 2}; /* 2^(6 7 4 5 2 3 0 1) */ 
  float half[3][3]={{0.5,1.0,1.0},{1.0,0.5,1.0},{1.0,1.0,0.5}},*verts;

  #include "tricases.h"

  *Tris=NULL;
  *Verts=NULL;
  *Tri_count=*Vert_count=0;

  for(i=0;i<256;i++)
    for(j=15;j<19;j++)
      tricases[i][j]=-1;

  if (connectivity>6) /* flip table up and down */
    for(i=0;i<128;i++)
      for(j=0;j<19;j++)
        {
        k=tricases[i][j]; tricases[i][j]=tricases[255-i][j]; tricases[255-i][j]=k;
        }
  if (connectivity==26) /* modify table for 26 nbr connectivity */
    for(i=0;i<19;i++)
      {
      tricases[65  /*190*/][i]=altcases[0][i];
      tricases[130 /*125*/][i]=altcases[1][i];
      tricases[40  /*215*/][i]=altcases[2][i];
      tricases[20  /*235*/][i]=altcases[3][i];
      }
  if (connectivity>6) /* switch triangle orientation if we've flipped */
    for(i=0;i<256;i++)
      for(j=0;j<18;j+=3)
        {
        k=tricases[i][j]; tricases[i][j]=tricases[i][j+1]; tricases[i][j+1]=k;
        }

  img_horiz=dims[0]; /* cols, fastest in memory */   
  img_vert=dims[1]; /* rows */
  img_depth=dims[2]; /* planes, slowest in mem*/

  ih1=img_horiz-1;
  iv1=img_vert-1;
  id1=img_depth-1;

  cellplane=ih1*iv1;
  cells=id1*cellplane;
  cellplanem[0]=0; cellplanem[1]=cellplane; cellplanem[2]=(cellplane<<1);

  v_idx=(int*)Gcalloc(cellplane*3*2,sizeof(int),0);
  if (v_idx==NULL) {error_msg("Memory error.\n",__LINE__);return(1);}

  status=(unsigned char*)Gcalloc(cells,sizeof(unsigned char),0);
  if (status==NULL) {error_msg("Memory error.\n",__LINE__);return(1);}

  for(i=0;i<256;i++) /* calculate number of tris, new verts for each case */
    {
    hep[0]=hep[1]=hep[2]=0;
    k=0;
    for(j=0;j<19;j++)
      {
      if (tricases[i][j]>=0) k++;
      for(w=0;w<3;w++) if (tricases[i][j]==hep2[w]) hep[w]=1;
      }
    tc[i]=k/3;
    vc[i]=hep[0]+hep[1]+hep[2];
    }

  /* calculate offsets to adjacent cells */
  offs[0]=0; offs[1]=1; offs[2]=img_horiz; offs[3]=ih1;
  for(i=4;i<8;i++) offs[i]=offs[i-4]+cellplane;

  /* calculate offsets to point into v_idx */
  /* v_idx holds vertex number data for two cell planes */
  /* cell edge midpoints 4..11 are in the current plane */
  /* 0..3 are in the previous plane */
  for(i=0;i<4;i++) planeidx[i]=0;
  for(i=4;i<12;i++) planeidx[i]=1;
  coffs[0]=-ih1; ctype[0]=1; /* cell edge midpoint type '6' */
  coffs[1]=0; ctype[1]=0;    /* cell edge midpoint type '5' */
  coffs[2]=0; ctype[2]=1;
  coffs[3]=-1; ctype[3]=0;
  for(i=4;i<8;i++) {coffs[i]=coffs[i-4]; ctype[i]=ctype[i-4];}
  coffs[8]=-ih1-1; ctype[8]=2; /* cell edge midpoint type '11' */
  coffs[9]=-ih1; ctype[9]=2;
  coffs[10]=-1; ctype[10]=2;
  coffs[11]=0; ctype[11]=2;
  /* end of setup */

  /* calculate status of each cell */
  pidx=J+1+img_horiz+(img_horiz*img_vert); /* begin at point (1,1,1) */
  cidx=status; /* begin at cell (0,0,0) */
  /* for each interior point */
  for (k=img_depth-2;k>=1;k--,pidx+=(img_horiz<<1),cidx+=ih1)
    for (j=img_vert-2;j>=1;j--,pidx+=2,cidx++)
      for (i=img_horiz-2;i>=1;i--,pidx++,cidx++)
        if (*pidx==val)
          for (w=0;w<8;w++) /* update cells which have pidx as a corner */
            *(cidx+offs[w])+=pows[w];

  /* quick count of the number of tris, number of verts needed */
  for (tri_count=vert_count=i=0;i<cells;i++)
    {
    tri_count+=tc[stat=status[i]];
    vert_count+=vc[stat];
    }
  tcm[0]=0; tcm[1]=tri_count; tcm[2]=(tri_count<<1);

  vert_count_times2=vert_count*2; /* just twice vert_count */

  /* allocate memory for tris and verts */
  tris=(int*)Gcalloc(tri_count*3,sizeof(int),g0->return_surface);
  if (tris==NULL) {error_msg("Memory error.\n",__LINE__);return(1);}

  verts=(float*)Gcalloc(vert_count*3,sizeof(float),g0->return_surface);
  if (verts==NULL) {error_msg("Memory error.\n",__LINE__);return(1);}

  *Tris=tris;
  *Verts=verts;
  *Vert_count=vert_count;
  *Tri_count=tri_count;

  /* extract the surface... */
  ic=0; /* number of cell within cell volume */
  vert_count1=tri_count1=0;
  for (k=0;k<id1;k++) /* for each cell plane */
    {
    planeoff[1]=(k%2)*3*cellplane; /* offset to current cell plane */
    planeoff[0]=3*cellplane-planeoff[1];  /* offset to previous plane */
    pc=0; /* number of cell within cell plane */
    for (j=0;j<iv1;j++)    /* for each cell */
      for (i=0;i<ih1;i++)
        {
        if (tc[stat=status[ic]]) /* no tris implies no verts */
          {
          ts=tricases[stat]; /* index into case table */

          /* create vertex numbers for this cell, if needed */
          hep[0]=hep[1]=hep[2]=0; /* looking for edge midpoints 5,6,11 */
          for (w=0;w<15;w++)
            {
            if(ts[w]<0) /* negative indicates end of vert list */
              {
              w=15; /* crude escape of w loop */
              }
            else
              {
              for(u=0;u<3;u++) /* look for 5,6,11 */
                if ((ts[w]==hep2[u]) && (hep[u]==0))
                  {
                  /* create a vertex number for the cell edge midpoint */
                  v_idx[planeoff[1]+pc+cellplanem[u]]=vert_count1;
                  hep[u]=1; /* be sure not to count it twice */
                  /* calc xyz coords of vert */
                  verts[vert_count1]= ((float)i) + half[u][1] ; 
                  verts[vert_count1+vert_count]= ((float)j) + half[u][0] ;
                  verts[vert_count1+vert_count_times2]= ((float)k) + half[u][2] ;
                  vert_count1++;
                  }
              }
            } /* end for create vertices */

          /* create triangles */
          u=3*tc[stat];
          for (w=0;w<u;w+=3) /* for each triangle */
            {
            for(v=0;v<3;v++) /* for each vert of triangle */
              {
              tsvw=ts[v+w];
              tris[tri_count1+tcm[v]]=v_idx[planeoff[planeidx[tsvw]]+pc
                +coffs[tsvw]+cellplanem[ctype[tsvw]]] ; /* whew! */
              }
            tri_count1++;
            } /* end for create triangles */

          } /* end if tc not zero */
        ic++;
        pc++;
        } /* end for i */
    } /* end for k */

  if (status!=NULL) Gfree(status);
  if (v_idx!=NULL) Gfree(v_idx);

  return(0); /* return with success status */

  } /*end!*/
 
 
static int big_component(int * Tris, float * Verts, int * Vert_count, int * Tri_count)
  {
  int ov,ot,tri_count,vert_count,*v_count,*tp,*v_ran,*v_idx,u0,v0,i,j,*tris,w;
  int *que,que_pos,que_len,*component,comp_count,*c_count,tc[3],k,qqp,vt,*tr;
  float *verts;

  ot=tri_count=*Tri_count;
  ov=vert_count=*Vert_count;

  tris=Tris;
  verts=Verts;

  if (verbose) print_msg("Extracting component with largest number of vertices...\n");

  if ((v_count=(int *)Gcalloc(vert_count,sizeof(int),0))==NULL) 
    {error_msg("Memory error.\n",__LINE__);return(1);}  

  tp=tris;
  for(i=0;i<tri_count;i++)
    {
    v_count[*(tp++)]++;
    v_count[*(tp++)]++;
    v_count[*(tp++)]++;
    }

  if ((v_ran=(int *)Gcalloc(vert_count+1,sizeof(int),0))==NULL)
    {error_msg("Memory error.\n",__LINE__);return(1);}  

  v_ran[0]=0;
  for(i=0;i<vert_count;i++)   
    {
    v_ran[i+1]=v_ran[i]+v_count[i];
    v_count[i]=0;
    }

  if ((v_idx=(int *)Gcalloc(v_ran[vert_count],sizeof(int),0))==NULL)
    {error_msg("Memory error.\n",__LINE__);return(1);}  

  tp=tris;
  for (j=0;j<3;j++)
    for(i=0;i<tri_count;i++)
      {
      v0=*(tp++); /* a vertex of triangle i */
      u0=v_ran[v0]+v_count[v0]; /* offset into v_idx */
      v_idx[u0]=i; /* triangle number */
      v_count[v0]++;
      }

  component=v_count;
  for(i=0;i<vert_count;i++) component[i]=0;

  if ((que=(int *)Gcalloc(vert_count,sizeof(int),0))==NULL)
    {error_msg("Memory error.\n",__LINE__);return(1);}  

  tc[0]=0; tc[1]=tri_count; tc[2]=tri_count*2;
  comp_count=0;
  for (vt=0;vt<vert_count;vt++)
    {
    if (component[vt]==0) /* if not processed */
      {
      comp_count++; /* increase component count */
      que_len=que_pos=0; /* reset que */
      component[que[que_len++]=vt]=comp_count; /* add vt to start of que */    
      while(que_pos<que_len) /* while que not empty */
        {
        qqp=que[que_pos];
        /* add nbrs to que */
        for(i=v_ran[qqp];i<v_ran[qqp+1];i++)
          {
          tr=tris+v_idx[i]; /* triangle containing qqp */
          for(j=0;j<3;j++)
            {
            k=tr[tc[j]]; /* a vertex nbr of qqp */
            if (component[k]==0) component[que[que_len++]=k]=comp_count;
            }
          }
        que_pos++;
        }    
      }
    }

  if (comp_count>1)
    {
    if ((c_count=(int *)Gcalloc(1+comp_count,sizeof(int),0))==NULL)
      {error_msg("Memory error.\n",__LINE__);return(1);}  

    for(i=0;i<vert_count;i++) c_count[component[i]]++;

    j=1; k=c_count[j];
    for(i=1;i<=comp_count;i++) if (k<c_count[i]) {j=i; k=c_count[i];}
    Gfree(c_count);

    k=j; /* k is number of largest component */
    j=0;
    for (vt=0;vt<vert_count;vt++)
      {
      if (component[vt]==k)
        {
        que[vt]=j; /* create remapping index */
        verts[j]=verts[vt]; /* move vertex data */
        verts[j+vert_count]=verts[vt+vert_count];
        verts[j+vert_count*2]=verts[vt+vert_count*2];
        j++;
        }
      else que[vt]=0;
      }

    w=0;
    for(vt=0;vt<tri_count;vt++)
      {
      tr=tris+vt;
      if (component[tr[0]]==k)
        {
        tris[w]=      que[tr[    0]];
        tris[w+tc[1]]=que[tr[tc[1]]];
        tris[w+tc[2]]=que[tr[tc[2]]];
        w++;
        }
      }

    tri_count=w;
    vert_count=j;    
    }

  if (tri_count!=ot)
    {
    for(i=0;i<tri_count;i++) tris[tri_count+i]=tris[ot+i];
    for(i=0;i<tri_count;i++) tris[(tri_count<<1)+i]=tris[(ot<<1)+i];   
    }

  if (vert_count!=ov)
    {
    for(i=0;i<vert_count;i++) verts[i+vert_count]=verts[i+ov]; 
    for(i=0;i<vert_count;i++) verts[i+(vert_count<<1)]=verts[i+(ov<<1)];
    }

  *Vert_count=vert_count;
  *Tri_count=tri_count;

  Gfree(que);
  Gfree(v_count);
  Gfree(v_ran);
  Gfree(v_idx);

  return(0);

  }
 
 
static int save_image(genus0parameters * g0)
  {
  int i,j,k,totlen,h,h1,sti;
  int *pad,dims[3],origlen,vo[3],v2[3],pos[3];
  char msg[200];
  unsigned short *output,zp;
  float *m,*verts,hv[3];

  pad=g0->pad;
  totlen=img_horiz*img_vert*img_depth;

  if (!(g0->any_genus)) /* if we made topological corrections */
    {
    j=0;
    for(i=0;i<totlen;i++)
      {
      zp=zpic[i];  /* original zpic value */
      sti=(status[i]<=10); /* not part of a component */
      zpic[i]=(g0->cut_loops)?(1-sti):sti; /* switch to complement */
      if ((zp==g0->cut_loops) && (status[i]<=10) &&(status[i]!=4))
        {
        status[i]=1; /* we wanted this one */
        j++;
        }
      else status[i]=0;
      }
    if (verbose) {sprintf(msg,"Made %d adjustments.\n",j);print_msg(msg);}
    Gfree(que);
    Gfree(cm);
    Gfree(fzpic);
    }

  /* Get the surface ! */
  dims[0]=img_horiz; dims[1]=img_vert; dims[2]=img_depth;
  if(GetSurf(zpic,1,dims,g0->connectivity,
    &(g0->triangles), &(g0->vertices),
    &(g0->tri_count),&(g0->vert_count),g0)) return(1);
  if (g0->biggest_component)
    if ((g0->tri_count>0)&&(g0->vert_count>0))
      if(big_component(g0->triangles,g0->vertices,&(g0->vert_count),&(g0->tri_count)))
        {error_msg("Error getting surface components.\n",__LINE__);return(1);}

  output=g0->output;
  if (g0->return_adjusted_label_map) /* they want a new label map back */
    {
    origlen=(g0->dims[0])*(g0->dims[1])*(g0->dims[2]);
    if (output==NULL) /* they didn't allocate.  So we need to */
      {
      if ((output=(unsigned short *)Gcalloc(origlen,sizeof(unsigned short),1))==NULL)
        {error_msg("Memory error.\n",__LINE__);return(1);}
      g0->calloced_output=1;
      }
    g0->output=output;
    for(i=0;i<origlen;i++) output[i]=(g0->input)[i]; /* just copy original input */
    if (!(g0->any_genus))
      {  
      for(k=0;k<(img_depth-pad[2]*2);k++)
        for(j=0;j<(img_vert-pad[1]*2);j++)
          for(i=0;i<(img_horiz-pad[0]*2);i++)
            {
            h1=sub2ind(pad[0]+i,pad[1]+j,pad[2]+k,img_horiz,img_vert);
            h=sub2ind(autocrop[0][0]+i,autocrop[1][0]+j,autocrop[2][0]+k,g0->dims[0],g0->dims[1]);
            if (status[h1]) output[h]=g0->alt_value;
            }
      }
    }
  vo[0]=0; vo[1]=g0->vert_count; vo[2]=(g0->vert_count)*2;
  v2[0]=1; v2[1]=dims[0]; v2[2]=dims[0]*dims[1];
  verts=g0->vertices;

  if (g0->return_adjusted_label_map) /* they want a new label map back */
    {
    if (g0->contour_value==65535) g0->contour_value=g0->alt_value;
    if (g0->alt_contour_value==65535) g0->alt_contour_value=g0->alt_value;
    for(h=0;h<g0->vert_count;h++)
      {
      j=0; /* j will be direction in which 0.5 position occurs */
      for (i=0;i<3;i++) 
        {
        pos[i]=(int)verts[h+vo[i]];
        if ((verts[h+vo[i]]-pos[i])==0.5) j=i; /* must happen exactly once */
        }
      i=sub2ind(pos[0],pos[1],pos[2],dims[0],dims[1]); /* index into zpic */
      if (zpic[i]==0) pos[j]++; /* subscript into boundary voxel of zpic */
      for (i=0;i<3;i++) pos[i]+=autocrop[i][0]-pad[i]; /* subscript into boundary voxel of input */
      i=sub2ind(pos[0],pos[1],pos[2],g0->dims[0],g0->dims[1]); /* index into input */
      output[i]=((g0->input)[i]==g0->value)?(g0->contour_value):(g0->alt_contour_value);
      }
    }

  /* need to translate verts back into input space.  zpic has been cropped */
  for(h=0;h<g0->vert_count;h++) for (i=0;i<3;i++) verts[h+vo[i]]+=autocrop[i][0]-pad[i];

  /* multiply verts by ijk2ras matrix */
  if (g0->ijk2ras!=NULL)
    {
    m=g0->ijk2ras;
    for(h=0;h<g0->vert_count;h++)
      {
      hv[0]=(g0->vertices)[h]; 
      hv[1]=(g0->vertices)[h+(g0->vert_count)]; 
      hv[2]=(g0->vertices)[h+((g0->vert_count)<<1)];
      (g0->vertices)[h]                      =m[0]*hv[0]+m[1]*hv[1]+m[ 2]*hv[2]+m[3];
      (g0->vertices)[h+g0->vert_count]       =m[4]*hv[0]+m[5]*hv[1]+m[ 6]*hv[2]+m[7];
      (g0->vertices)[h+((g0->vert_count)<<1)]=m[8]*hv[0]+m[9]*hv[1]+m[10]*hv[2]+m[11];
      }
    }

  Gfree(zpic);
  if (!(g0->any_genus)) Gfree(status);
  if (!(g0->return_surface)) /* don't return surface if they didn't want it */
    {
    Gfree(g0->vertices); g0->vertices=NULL;
    Gfree(g0->triangles); g0->triangles=NULL;
    }

  if (verbose) 
    {
    sprintf(msg,"Vertices: %d  Triangles: %d\n",g0->vert_count,g0->tri_count);
    print_msg(msg);
    }
  return(0);
  }

extern void genus0init(genus0parameters * g0)
  {
  g0->input=NULL;
  g0->return_adjusted_label_map=1;
  g0->return_surface=1;
  (g0->dims)[0]=(g0->dims)[1]=(g0->dims)[2]=0;
  g0->alt_value=g0->value=0;
  g0->contour_value=65535;
  g0->alt_contour_value=65535;
  g0->cut_loops=0;
  g0->pad[0]=g0->pad[1]=g0->pad[2]=2;
  g0->connectivity=6;
  g0->any_genus=0;
  g0->biggest_component=1;
  g0->connected_component=1;
  g0->ijk2ras=NULL;
  (g0->extraijkscale)[0]=(g0->extraijkscale)[1]=(g0->extraijkscale)[2]=1.0;
  g0->verbose=0;

  g0->output=NULL;
  g0->vert_count=0;
  g0->vertices=NULL;
  g0->tri_count=0;
  g0->triangles=NULL;
  g0->calloced_output=0; /* private */
  }
 
 

extern int genus0(genus0parameters * g0)
  {
  int thistenth,lasttenth,level;

  if (set_up(g0)) {Gfree_all(0);return(1);}
  if (!(g0->any_genus))
    {
    if (verbose) printf("Starting main process...\n");
    lasttenth=0;
    for(level=maxlevels;level>=1;level--) /* [maxlevels...1] */
      {
      find_component(level);
      thistenth=(int)(((float)(maxlevels-level))/(maxlevels-1.0)*10.0+0.5);
      if (thistenth!=lasttenth)
        {
        lasttenth=thistenth;
        if (verbose) printf("Done with %d percent.\n",thistenth*10);
        if (thistenth<10)
          if (cmtostat()) {Gfree_all(0);return(1);}
        }
      }
    if(cmtostat()) {Gfree_all(0);return(1);}
    }
  if (save_image(g0)) {Gfree_all(0);return(1);}
  Gfree_all(1); /* isn't really be needed, if we've freed our calloc's */
  return(0); /* normal, error free return */
  }

extern void genus0destruct(genus0parameters * g0)
  {
  if (g0->vertices!=NULL) 
    {
    basic_free(g0->vertices);
    g0->vertices=NULL;
    }
  if (g0->triangles!=NULL) 
    {
    basic_free(g0->triangles);
    g0->triangles=NULL;
    }
  if (g0->calloced_output)
    if (g0->output!=NULL)
      {
      basic_free(g0->output);
      g0->output=NULL;
      g0->calloced_output=0;
      }
  }
