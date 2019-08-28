
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <ctype.h>

#define     MAX_LINE_LENGTH 256
#define     MAX_INTEGER 2147483647
#define     MAX_FLOAT   1.e32

typedef struct {
  char  model[MAX_LINE_LENGTH]; /* model-name header*/
  char  **uheader; /* user header */
  char  **pheader; /* project header (remark: leading dataset-project-headers are stored in the related dataset!) */
  int   u;         /* number of user headers */
  int   p;         /* number of project headers */
  int   n;         /* number of nodes */
  int   e;         /* number of elements  */
  int   f;         /* number of faces */
  int   g;         /* number of edges */
  int   t;         /* number of texts */
  int   sets;      /* sets (groups) of entities */
  int   mats;      /* materials   */
  int   amps;      /* amplitudes  */
  int   l;         /* number of loadcases (Datasets) */
  int   b;         /* number of nodeBlocks */
  int   c;         /* number of 'cuts' over all nodeBlocks (block-to-block interfaces for isaac) */
  int   etype[100];/* number of elements of a certain type */
  int   nmax;      /* maximum nodenumber */
  int   nmin;      /* minimum nodenumber */
  int   emax;      /* maximum elemnumber */
  int   emin;      /* minimum elemnumber */
  int   orignmax;  /* max-node-nr-of-original-nodes (w/o nodes for drawing purposes) */
  int   orign;     /* nr-of-original-nodes (w/o nodes for drawing purposes) */
  int   olc;       /* nr-of-original-loadcases (w/o cgx generated datasets (lc)) */
  int   noffs;     /* node-nr offset, defined with asgn */
  int   eoffs;     /* elem-nr offset, defined with asgn */
} Summen;


typedef struct {
  int   nr;              /*   external node-nr (node[node-indx].nr) */
  int   indx;            /*   node-index (node[ext-node-nr].indx)   */
  char  pflag;           /*   1 if used for display purposes    */
                         /*  -1 if the node is deleted          */
                         /*   0 default                         */
  double nx;             /*   coordinates  node[ext-node-nr].nx */
  double ny;
  double nz;
} Nodes;

typedef struct {
  int nr;                /* external element-nr */
  //  int indx;              /* -index (elem[external elem-nr].indx)   */
  int type;              /* element type (1:Hexa8)  */
  int group;
  int mat;
  int attr;              /* -1: unstructured mesh tr3u (-2 for mesh with libGLu tr3g ) */
                         /*  0: default           */
                         /*  1: reduced integration he8r he20r */
                         /*  2: incompatible modes he8i */
                         /*  3: DASHPOTA be2d */
                         /*  4: plane strain (CPE) tr3e tr6e qu4e qu8e */
                         /*  5: plane stress (CPS) tr3s */
                         /*  6: axisymmetric  (CAX) tr3c */
                         /*  7: fluid he8f */
                         /*  8: tet10m */
  int nod[27];
  double **side;         /* side[Nr.][x|y|z]== normal vector */
} Elements;

typedef struct {
  char  **pheader;    /* project header */
  int   npheader;              /* number of headers */
  char  **compName;
  char  **icname;
  char  name[MAX_LINE_LENGTH];
  char  dataset_name[MAX_LINE_LENGTH];
  char  dataset_text[MAX_LINE_LENGTH];
  char  analysis_name[MAX_LINE_LENGTH];
  float value;
  char  filename[MAX_LINE_LENGTH];
  FILE *handle;
  fpos_t *fileptr;
  int   loaded;       /* if data are stored:1 else: 0 */
  int format_flag;
  int analysis_type;
  int step_number;
  int ncomps;
  int irtype;
  int *menu;
  int *ictype;
  int *icind1;
  int *icind2;
  int *iexist;
  float **dat;        /* node related data */
  float ***edat;      /* element related data */
  float *max;         /* maximum datum */
  float *min;         /* minimum datum */
  int *nmax;          /* node with maximum datum */
  int *nmin;          /* node with minimum datum */
} Datasets;

void freeDatasets(Datasets *lcase, int nr);
int readfrd(char *datin, Summen *anz, Nodes **nptr, Elements **eptr, Datasets **lptr, int read_mode );
int readfrdblock( int lc, Summen *anz,   Nodes     *node, Datasets *lcase );
int stoi(char *string, int a, int b);
double stof(char *string, int a, int b);
void stos(char *string, int a, int b, char *puffer);
int compare (char *str1, char *str2, int length);
int frecord( FILE *handle1,  char *string);
int elemChecker(int sum_e, Nodes *node, Elements *elem);
void v_result( double *A, double *B, double *C );
void v_prod( double *A, double *B, double *C );
double v_betrag(double *a);
int strsplt( char *rec_str, char breakchar, char ***ptr);
