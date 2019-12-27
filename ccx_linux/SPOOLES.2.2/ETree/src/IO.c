/*  IO.c  */

#include "../ETree.h"

static const char *suffixb = ".etreeb" ;
static const char *suffixf = ".etreef" ;

/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   purpose -- to read in an ETree object from a file

   input --

      fn -- filename, must be *.etreeb or *.etreef

   return value -- 1 if success, 0 if failure

   created -- 95nov15, cca
   -------------------------------------------------
*/
int
ETree_readFromFile ( 
   ETree   *etree, 
   char    *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || fn == NULL ) {
   fprintf(stderr, "\n error in ETree_readFromFile(%p,%s)"
           "\n bad input\n", etree, fn) ;
   return(0) ;
}
/*
   -------------
   read the file
   -------------
*/
fnlength = strlen(fn) ;
sulength = strlen(suffixb) ;
if ( fnlength > sulength ) {
   if ( strcmp(&fn[fnlength-sulength], suffixb) == 0 ) {
      if ( (fp = fopen(fn, "rb")) == NULL ) {
         fprintf(stderr, "\n error in ETree_readFromFile(%p,%s)"
                 "\n unable to open file %s", etree, fn, fn) ;
         rc = 0 ;
      } else {
         rc = ETree_readFromBinaryFile(etree, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "r")) == NULL ) {
         fprintf(stderr, "\n error in ETree_readFromFile(%p,%s)"
                 "\n unable to open file %s", etree, fn, fn) ;
         rc = 0 ;
      } else {
         rc = ETree_readFromFormattedFile(etree, fp) ;
         fclose(fp) ;
      }
   } else {
      fprintf(stderr, "\n error in ETree_readFromFile(%p,%s)"
              "\n bad ETree file name %s,"
              "\n must end in %s (binary) or %s (formatted)\n",
              etree, fn, fn, suffixb, suffixf) ;
      rc = 0 ;
   }
} else {
   fprintf(stderr, "\n error in ETree_readFromFile(%p,%s)"
       "\n bad ETree file name %s,"
       "\n must end in %s (binary) or %s (formatted)\n",
       etree, fn, fn, suffixb, suffixf) ;
   rc = 0 ;
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   purpose -- to read an ETree object from a formatted file

   return value -- 1 if success, 0 if failure

   created -- 95nov15, cca
   --------------------------------------------------------
*/
int
ETree_readFromFormattedFile ( 
   ETree   *etree, 
   FILE    *fp 
) {
int    rc ;
int    itemp[2] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in ETree_readFromFormattedFile(%p,%p)"
           "\n bad input\n", etree, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
ETree_clearData(etree) ;
/*
   ---------------------------
   initialize the ETree object
   ---------------------------
*/
ETree_init1(etree, 0, 0) ;
/*
   -----------------------------
   read in the two scalar fields
   -----------------------------
*/
if ( (rc = IVfscanf(fp, 2, itemp)) != 2 ) {
   fprintf(stderr, "\n error in ETree_readFromFormattedFile(%p,%p)"
           "\n %d items of %d read\n", etree, fp, rc, 2) ;
   return(0) ;
}
etree->nfront = itemp[0] ;
etree->nvtx   = itemp[1] ;
/*
   -----------------------
   read in the Tree object
   -----------------------
*/
Tree_readFromFormattedFile(etree->tree, fp) ;
/*
   ------------------------------
   read in the nodwghts IV object
   ------------------------------
*/
IV_readFromFormattedFile(etree->nodwghtsIV, fp) ;
/*
   ------------------------------
   read in the bndwghts IV object
   ------------------------------
*/
IV_readFromFormattedFile(etree->bndwghtsIV, fp) ;
/*
   --------------------------------
   read in the vtxToFront IV object
   --------------------------------
*/
IV_readFromFormattedFile(etree->vtxToFrontIV, fp) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   purpose -- to read an ETree object from a binary file

   return value -- 1 if success, 0  if failure

   created -- 95nov15, cca
   ----------------------------------------------------
*/
int
ETree_readFromBinaryFile ( 
   ETree    *etree, 
   FILE   *fp 
) {
int    rc ;
int    itemp[2] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || fp == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_readFromBinaryFile(%p,%p)"
           "\n bad input\n", etree, fp) ;
   return(0) ;
}
/*
   ---------------------
   clear the data fields
   ---------------------
*/
ETree_clearData(etree) ;
/*
   ---------------------------
   initialize the ETree object
   ---------------------------
*/
ETree_init1(etree, 0, 0) ;
/*
   -----------------------------
   read in the two scalar fields
   -----------------------------
*/
if ( (rc = fread((void *) itemp, sizeof(int), 2, fp)) != 2 ) {
   fprintf(stderr, "\n error in ETree_readFromBinaryFile(%p,%p)"
           "\n itemp(2) : %d items of %d read\n", etree, fp, rc, 2) ;
  return(0) ;
}
etree->nfront = itemp[0] ;
etree->nvtx   = itemp[1] ;
/*
   -----------------------
   read in the Tree object
   -----------------------
*/
Tree_readFromBinaryFile(etree->tree, fp) ;
/*
   ------------------------------
   read in the nodwghts IV object
   ------------------------------
*/
IV_readFromBinaryFile(etree->nodwghtsIV, fp) ;
/*
   ------------------------------
   read in the bndwghts IV object
   ------------------------------
*/
IV_readFromBinaryFile(etree->bndwghtsIV, fp) ;
/*
   --------------------------------
   read in the vtxToFront IV object
   --------------------------------
*/
IV_readFromBinaryFile(etree->vtxToFrontIV, fp) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   purpose -- to write an ETree object to a file

   input --

      fn -- filename
        *.etreeb -- binary
        *.etreef -- formatted
        anything else -- for human eye

   return value -- 1 if success, 0 otherwise

   created -- 95nov15, cca
   --------------------------------------------
*/
int
ETree_writeToFile ( 
   ETree   *etree, 
   char   *fn 
) {
FILE   *fp ;
int    fnlength, rc, sulength ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || fn == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_writeToFile(%p,%s)"
    "\n bad input\n", etree, fn) ; 
}
/*
   ------------------
   write out the file
   ------------------
*/
fnlength = strlen(fn) ;
sulength = strlen(suffixb) ;
if ( fnlength > sulength ) {
   if ( strcmp(&fn[fnlength-sulength], suffixb) == 0 ) {
      if ( (fp = fopen(fn, "wb")) == NULL ) {
         fprintf(stderr, "\n error in ETree_writeToFile(%p,%s)"
                 "\n unable to open file %s", etree, fn, fn) ;
         rc = 0 ;
      } else {
         rc = ETree_writeToBinaryFile(etree, fp) ;
         fclose(fp) ;
      }
   } else if ( strcmp(&fn[fnlength-sulength], suffixf) == 0 ) {
      if ( (fp = fopen(fn, "w")) == NULL ) {
         fprintf(stderr, "\n error in ETree_writeToFile(%p,%s)"
                 "\n unable to open file %s", etree, fn, fn) ;
         rc = 0 ;
      } else {
         rc = ETree_writeToFormattedFile(etree, fp) ;
         fclose(fp) ;
      }
   } else {
      if ( (fp = fopen(fn, "a")) == NULL ) {
         fprintf(stderr, "\n error in ETree_writeToFile(%p,%s)"
                 "\n unable to open file %s", etree, fn, fn) ;
         rc = 0 ;
      } else {
         rc = ETree_writeForHumanEye(etree, fp) ;
         fclose(fp) ;
      }
   }
} else {
   if ( (fp = fopen(fn, "a")) == NULL ) {
      fprintf(stderr, "\n error in ETree_writeToFile(%p,%s)"
              "\n unable to open file %s", etree, fn, fn) ;
      rc = 0 ;
   } else {
      rc = ETree_writeForHumanEye(etree, fp) ;
      fclose(fp) ;
   }
}
return(rc) ; }

/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------------
   purpose -- to write an ETree object to a formatted file

   return value -- 1 if success, 0 otherwise

   created -- 95nov15, cca
   ------------------------------------------------------
*/
int
ETree_writeToFormattedFile ( 
   ETree   *etree, 
   FILE    *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || fp == NULL || etree->tree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_writeToFormattedFile(%p,%p)"
           "\n bad input\n", etree, fp) ;
   exit(-1) ;
}
/*
   ---------------------------
   write the two scalar fields
   ---------------------------
*/
rc = fprintf(fp, "\n %d %d", etree->nfront, etree->nvtx) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in ETree_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from first fprintf\n", etree, fp, rc) ;
   return(0) ;
}
/*
   ---------------------------------
   write the Tree object to the file 
   ---------------------------------
*/
rc = Tree_writeToFormattedFile(etree->tree, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in ETree_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from writing Tree to file\n", 
           etree, fp, rc) ;
   return(0) ;
}
/*
   ----------------------------------------
   write the nodwghts IV object to the file 
   ----------------------------------------
*/
rc = IV_writeToFormattedFile(etree->nodwghtsIV, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in ETree_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from writing nodwghtsIV to file\n", 
           etree, fp, rc) ;
   return(0) ;
}
/*
   ----------------------------------------
   write the bndwghts IV object to the file 
   ----------------------------------------
*/
rc = IV_writeToFormattedFile(etree->bndwghtsIV, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in ETree_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from writing bndwghtsIV to file\n", 
           etree, fp, rc) ;
   return(0) ;
}
/*
   ------------------------------------------
   write the vtxToFront IV object to the file 
   ------------------------------------------
*/
rc = IV_writeToFormattedFile(etree->vtxToFrontIV, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in ETree_writeToFormattedFile(%p,%p)"
           "\n rc = %d, return from writing vtxToFrontIV to file\n", 
           etree, fp, rc) ;
   return(0) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to write an ETree object to a binary file

   return value -- 1 if success, 0 otherwise

   created -- 95nov15, cca
   ---------------------------------------------------
*/
int
ETree_writeToBinaryFile ( 
   ETree    *etree, 
   FILE   *fp 
) {
int   rc ;
int   itemp[2] ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || fp == NULL || etree->tree == NULL ) {
   fprintf(stderr, "\n fatal error in ETree_writeToBinaryFile(%p,%p)"
           "\n bad input\n", etree, fp) ;
   exit(-1) ;
}
/*
   ---------------------------
   write the two scalar fields
   ---------------------------
*/
itemp[0] = etree->nfront ;
itemp[1] = etree->nvtx   ;
rc = fwrite((void *) itemp, sizeof(int), 2, fp) ;
if ( rc != 2 ) {
   fprintf(stderr, "\n error in ETree_writeToBinaryFile(%p,%p)"
           "\n %d of %d scalar items written\n", etree, fp, rc, 2) ;
   return(0) ;
}
/*
   ---------------------------------
   write the Tree object to the file 
   ---------------------------------
*/
rc = Tree_writeToBinaryFile(etree->tree, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in ETree_writeToBinaryFile(%p,%p)"
           "\n rc = %d, return from writing Tree to file\n", 
           etree, fp, rc) ;
   return(0) ;
}
/*
   ----------------------------------------
   write the nodwghts IV object to the file 
   ----------------------------------------
*/
rc = IV_writeToBinaryFile(etree->nodwghtsIV, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in ETree_writeToBinaryFile(%p,%p)"
           "\n rc = %d, return from writing nodwghtsIV to file\n", 
           etree, fp, rc) ;
   return(0) ;
}
/*
   ----------------------------------------
   write the bndwghts IV object to the file 
   ----------------------------------------
*/
rc = IV_writeToBinaryFile(etree->bndwghtsIV, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in ETree_writeToBinaryFile(%p,%p)"
           "\n rc = %d, return from writing bndwghtsIV to file\n", 
           etree, fp, rc) ;
   return(0) ;
}
/*
   ------------------------------------------
   write the vtxToFront IV object to the file 
   ------------------------------------------
*/
rc = IV_writeToBinaryFile(etree->vtxToFrontIV, fp) ;
if ( rc < 0 ) {
   fprintf(stderr, "\n fatal error in ETree_writeToBinaryFile(%p,%p)"
           "\n rc = %d, return from writing vtxToFrontIV to file\n", 
           etree, fp, rc) ;
   return(0) ;
}

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------
   purpose -- to write an ETree object for a human eye

   return value -- 1 if success, 0 otherwise

   created -- 95nov15, cca
   ---------------------------------------------------
*/
int
ETree_writeForHumanEye ( 
   ETree    *etree, 
   FILE   *fp 
) {
int   nfront, rc, v ;
int   *bndwghts, *fch, *nodwghts, *par, *sib ;

if ( etree == NULL || fp == NULL || (nfront = etree->nfront) <= 0 ) {
   fprintf(stderr, "\n fatal error in ETree_writeForHumanEye(%p,%p)"
           "\n bad input\n", etree, fp) ;
   exit(-1) ;
}
if ( (rc = ETree_writeStats(etree, fp)) == 0 ) {
   fprintf(stderr, "\n fatal error in ETree_writeForHumanEye(%p,%p)"
           "\n rc = %d, return from ETree_writeStats(%p,%p)\n",
           etree, fp, rc, etree, fp) ;
   return(0) ;
}
par = etree->tree->par ;
fch = etree->tree->fch ;
sib = etree->tree->sib ;
nodwghts = IV_entries(etree->nodwghtsIV) ;
bndwghts = IV_entries(etree->bndwghtsIV) ;
fprintf(fp, 
        "\n front    parent   fchild   sibling   nodwght   bndwght") ;
for ( v = 0 ; v < nfront ; v++ ) {
   fprintf(fp, "\n %5d %9d %9d %9d %9d %9d ", 
          v, par[v], fch[v], sib[v], nodwghts[v], bndwghts[v]) ;
}
fflush(fp) ;
fprintf(fp, "\n\n vtxToFront IV object") ;
IV_writeForHumanEye(etree->vtxToFrontIV, fp) ;
fflush(fp) ;

return(1) ; }

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------------------------
   purpose -- to write out the statistics for the ETree object

   return value -- 1 if success, 0 otherwise

   created -- 95nov15, cca
   -----------------------------------------------------------
*/
int
ETree_writeStats ( 
   ETree    *etree, 
   FILE   *fp 
) {
int   rc ;
/*
   ---------------
   check the input
   ---------------
*/
if ( etree == NULL || fp == NULL ) {
   fprintf(stderr, "\n error in ETree_writeStats(%p,%p)"
           "\n bad input\n", etree, fp) ;
   exit(-1) ;
}
rc = fprintf(fp, 
      "\n ETree : etree object, %d fronts, %d vertices, takes %d bytes",
      etree->nfront, etree->nvtx, ETree_sizeOf(etree)) ;
if ( rc < 0 ) { 
   fprintf(stderr, "\n fatal error in ETree_writeStats(%p,%p)"
           "\n rc = %d, return from fprintf\n", etree, fp, rc) ;
   return(0) ;
}
return(1) ; }

/*--------------------------------------------------------------------*/
