/*  testSerial.c  */

#include "../Bridge.h"

void Factor ( ) ;
void MatMul ( ) ;
void Solve ( ) ;

/*--------------------------------------------------------------------*/

void main ( int argc, char *argv[] )
/*
   ----------------------------------------------------------
   read in Harwell-Boeing matrices, use serial factor, solve,
   and multiply routines based on spooles, invoke eigensolver

   created  -- 98mar31 jcp
   modified -- 98dec18, cca
   ----------------------------------------------------------
*/
{
Bridge    bridge ;
char      *inFileName_A, *inFileName_B, *outFileName, 
          *parmFileName, *type ;
char      buffer[20], pbtype[4], which[4] ;
double    lftend, rhtend, center, shfscl, t1, t2 ;
double    c__1 = 1.0, c__4 = 4.0, tolact = 2.309970868130169e-11 ;
double    eigval[1000], sigma[2];
double    *evec;
int       error, fstevl, lfinit, lstevl, mxbksz, msglvl, ncol, ndiscd,
          neig, neigvl, nfound, nnonzeros, nrhs, nrow, prbtyp, rc, 
          retc, rfinit, seed, warnng ;
int       c__5 = 5, output = 6 ;
int       *lanczos_wksp;
InpMtx    *inpmtxA, *inpmtxB ;
FILE      *msgFile, *parmFile;

/*--------------------------------------------------------------------*/

if ( argc != 7 ) {
   fprintf(stdout, 
  "\n\n usage : %s msglvl msgFile parmFile seed inFileA inFileB"
  "\n    msglvl   -- message level"
  "\n    msgFile  -- message file"
  "\n    parmFile -- input parameters file"
  "\n    seed     -- random number seed, used for ordering"
  "\n    inFileA -- stiffness matrix in Harwell-Boeing format"
  "\n    inFileB -- mass matrix in Harwell-Boeing format"
  "\n               used for prbtyp = 1 or 2"
  "\n", argv[0]) ;
   return ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(argv[2], "a")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n",
           argv[0], argv[2]) ;
   exit(-1) ;
}
parmFileName = argv[3] ;
seed         = atoi(argv[4]) ;
inFileName_A = argv[5] ;
inFileName_B = argv[6] ;
fprintf(msgFile, 
        "\n %s "
        "\n msglvl         -- %d" 
        "\n msgFile        -- %s" 
        "\n parmFile       -- %s" 
        "\n seed           -- %d" 
        "\n stiffness file -- %s" 
        "\n mass file      -- %s" 
        "\n",
        argv[0], msglvl, argv[2], parmFileName, seed, 
        inFileName_A, inFileName_B) ;
fflush(msgFile) ;
/*
   ---------------------------------------------
   read in the Harwell-Boeing matrix information
   ---------------------------------------------
*/
if ( strcmp(inFileName_A, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
MARKTIME(t1) ;
readHB_info (inFileName_A, &nrow, &ncol, &nnonzeros, &type, &nrhs) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : read in header information for A",
        t2 - t1) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   read in eigenvalue problem data
   neigvl -- # of desired eigenvalues
   which  -- which eigenvalues to compute
     'l' or 'L' lowest (smallest magnitude)
     'h' or 'H' highest (largest magnitude)
     'n' or 'N' nearest to central value
     'c' or 'C' nearest to central value
     'a' or 'A' all eigenvalues in interval
   pbtype -- type of problem
     'v' or 'V' generalized symmetric problem (K,M)
                with M positive semidefinite (vibration problem)
     'b' or 'B' generalized symmetric problem (K,K_s)
                with K positive semidefinite
                with K_s posibly indefinite (buckling problem)
     'o' or 'O' ordinary symmetric eigenproblem
   lfinit -- if true, lftend is restriction on lower bound of 
             eigenvalues. if false, no restriction on lower bound
   lftend -- left endpoint of interval
   rfinit -- if true, rhtend is restriction on upper bound of
             eigenvalues.  if false, no restriction on upper bound
   rhtend -- right endpoint of interval
   center -- center of interval
   mxbksz -- upper bound on block size for Lanczos recurrence
   shfscl -- shift scaling parameter, an estimate on the magnitude
             of the smallest nonzero eigenvalues
   ---------------------------------------------------------------
*/
MARKTIME(t1) ;
parmFile = fopen(parmFileName, "r");
fscanf(parmFile, "%d %s %s %d %le %d %le %le %d %le", 
       &neigvl, which, pbtype, &lfinit, &lftend, 
       &rfinit, &rhtend, &center, &mxbksz, &shfscl) ;
fclose(parmFile);
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : read in eigenvalue problem data",
        t2 - t1) ;
/*
   ----------------------------------------
   check and set the problem type parameter
   ----------------------------------------
*/
switch ( pbtype[1] ) {
case 'v' : case 'V' : prbtyp = 1 ; break ;
case 'b' : case 'B' : prbtyp = 2 ; break ;
case 'o' : case 'O' : prbtyp = 3 ; break ;
default :
   fprintf(stderr, "\n invalid problem type %s", pbtype) ;
   exit(-1) ;
}
/*
   ----------------------------
   Initialize Lanczos workspace
   ----------------------------
*/
MARKTIME(t1) ;
lanczos_init_ ( &lanczos_wksp ) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : initialize lanczos workspace", 
        t2 - t1) ;
/*
   ----------------------------------
   initialize communication structure
   ----------------------------------
*/
MARKTIME(t1) ;
lanczos_set_parm( &lanczos_wksp, "order-of-problem",   &nrow,   &retc );
lanczos_set_parm( &lanczos_wksp, "accuracy-tolerance", &tolact, &retc );
lanczos_set_parm( &lanczos_wksp, "max-block-size",     &mxbksz, &retc );
lanczos_set_parm( &lanczos_wksp, "shift-scale",        &shfscl, &retc );
lanczos_set_parm( &lanczos_wksp, "message_level",      &msglvl, &retc );
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : init lanczos communication structure", 
        t2 - t1) ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------
   create the InpMtx objects for matrix A and B
   ---------------------------------------------
*/
if ( strcmp(inFileName_A, "none") == 0 ) {
   fprintf(msgFile, "\n no file to read from") ;
   exit(0) ;
}
MARKTIME(t1) ;
inpmtxA = InpMtx_new() ;
InpMtx_readFromHBfile ( inpmtxA, inFileName_A ) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : read in A", t2 - t1) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n InpMtx A object after loading") ;
   InpMtx_writeForHumanEye(inpmtxA, msgFile) ;
   fflush(msgFile) ;
}
MARKTIME(t1) ;
lanczos_set_parm( &lanczos_wksp, "matrix-type", &c__1, &retc );
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : set A's parameters", t2 - t1) ;
if ( prbtyp != 3 ) {
   if ( strcmp(inFileName_B, "none") == 0 ) {
      fprintf(msgFile, "\n no file to read from") ;
      exit(0) ;
   }
   MARKTIME(t1) ;
   inpmtxB = InpMtx_new() ;
   InpMtx_readFromHBfile ( inpmtxB, inFileName_B ) ;
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %8.3f : read in B", t2 - t1) ;
} else {
   MARKTIME(t1) ;
   inpmtxB = NULL ;
   lanczos_set_parm( &lanczos_wksp, "matrix-type", &c__4, &retc );
   MARKTIME(t2) ;
   fprintf(msgFile, "\n CPU %8.3f : set B's parameters", t2 - t1) ;
}
if ( msglvl > 2  && prbtyp != 3 ) {
   fprintf(msgFile, "\n\n InpMtx B object after loading") ;
   InpMtx_writeForHumanEye(inpmtxB, msgFile) ;
   fflush(msgFile) ;
 }
/*
   -----------------------------
   set up the solver environment
   -----------------------------
*/
MARKTIME(t1) ;
rc = Setup((void *) &bridge, &prbtyp, &nrow, &mxbksz, inpmtxA, inpmtxB,
           &seed, &msglvl, msgFile) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : set up solver environment", t2 - t1) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n fatal error %d from Setup()", rc) ;
   exit(-1) ;
}
/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------
   invoke eigensolver
   nfound -- # of eigenvalues found and kept
   ndisc  -- # of additional eigenvalues discarded
   -----------------------------------------------
*/
MARKTIME(t1) ;
lanczos_run(&neigvl, &which[1] , &pbtype[1], &lfinit, &lftend, 
	    &rfinit, &rhtend, &center, &lanczos_wksp, &bridge, &nfound,
	    &ndiscd, &warnng, &error, Factor, MatMul, Solve ) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : time for lanczos run", t2 - t1) ;
/*
   -------------------------
   get eigenvalues and print
   -------------------------
*/
MARKTIME(t1) ;
neig   = nfound + ndiscd ;
lstevl = nfound ;
lanczos_eigenvalues (&lanczos_wksp, eigval, &neig, &retc);
fstevl = 1 ;
if ( nfound == 0 ) fstevl = -1 ;
if ( ndiscd > 0 ) lstevl = -ndiscd ;
hdslp5_ ("computed eigenvalues returned by hdserl",
         &neig, eigval, &output, 39L ) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : get and print eigenvalues ", t2 - t1) ;
/*
   -------------------------
   get eigenvectors and print
   -------------------------
*/
/*
MARKTIME(t1) ;
neig = min ( 50, nrow );
Lncz_ALLOCATE(evec, double, nrow, retc);

for ( i = 1 ; i <= nfound ; i++ ) {
   lanczos_eigenvector ( &lanczos_wksp, &i, &i, newToOld,
                        evec, &nrow, &retc ) ;
   hdslp5_ ( "computed eigenvector returned by hdserc",
             &neig, evec, &output, 39L ) ;
}
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : get and print eigenvectors ", t2 - t1) ;
*/
/*
   ------------------------
   free the working storage
   ------------------------
*/
MARKTIME(t1) ;
lanczos_free( &lanczos_wksp ) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : free lanczos workspace ", t2 - t1) ;
MARKTIME(t1) ;
rc = Cleanup(&bridge) ;
MARKTIME(t2) ;
fprintf(msgFile, "\n CPU %8.3f : free solver workspace ", t2 - t1) ;
if ( rc != 1 ) {
   fprintf(stderr, "\n error return %d from Cleanup()", rc) ;
   exit(-1) ;
}
fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return ; }

