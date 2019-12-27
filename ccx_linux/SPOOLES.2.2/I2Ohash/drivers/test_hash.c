/*  test_hash.c  */

#include "../I2Ohash.h"
#include "../../Drand.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] )
/*
   -------------------------------------------
   generate random (key1, key2, value) entries
   and place them into a hash table.

   created -- 98jan28, cca
   -------------------------------------------
*/
{
Drand     *drand ;
FILE      *msgFile ;
I2Ohash   *hashtbl ;
int       grow, ii, key1, key2, maxkey, msglvl, nent, seed, size ;

if ( argc != 8 ) {
   fprintf(stdout, 
      "\n\n %% usage : %s msglvl msgFile size grow maxkey nent seed"
      "\n %%    msglvl  -- message level"
      "\n %%    msgFile -- message file"
      "\n %%    size    -- number of lists in the hash table"
      "\n %%    grow    -- growth for list items"
      "\n %%    maxkey  -- maximum key value"
      "\n %%    nent    -- number of entries to insert"
      "\n %%    seed    -- random number seed"
      "\n", argv[0]) ;
   return(0) ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else if ( (msgFile = fopen(argv[2], "a")) == NULL ) {
   fprintf(stderr, "\n fatal error in %s"
           "\n unable to open file %s\n",
           argv[0], argv[2]) ;
   return(-1) ;
}
size   = atoi(argv[3]) ;
grow   = atoi(argv[4]) ;
maxkey = atoi(argv[5]) ;
nent   = atoi(argv[6]) ;
seed   = atoi(argv[7]) ;
/*
fprintf(msgFile, 
        "\n %% %s "
        "\n %% msglvl  -- %d" 
        "\n %% msgFile -- %s" 
        "\n %% size    -- %d" 
        "\n %% grow    -- %d" 
        "\n %% maxkey  -- %d" 
        "\n %% nent    -- %d" 
        "\n",
        argv[0], msglvl, argv[2], size, grow, maxkey, nent) ;
fflush(msgFile) ;
*/
/*
   --------------------------------
   initialize the hash table object
   --------------------------------
*/
hashtbl = I2Ohash_new() ;
I2Ohash_init(hashtbl, size, nent/2, grow) ;
drand = Drand_new() ;
Drand_setSeed(drand, seed) ;
Drand_setUniform(drand, 0, maxkey) ;
/*
   ------------------------------------------------------------
   generate (key1, key2, value) pairs and insert into the table
   ------------------------------------------------------------
*/
for ( ii = 0 ; ii < nent ; ii++ ) {
   key1 = Drand_value(drand) ;
   key2 = Drand_value(drand) ;
   if ( msglvl > 2 ) {
      fprintf(msgFile, "\n inserting <%d,%d,NULL> into hash table",
              key1, key2) ;
   }
   I2Ohash_insert(hashtbl, key1, key2, NULL) ;
   if ( msglvl > 1 ) {
      fprintf(msgFile, "\n %d entries, measure = %.3f", 
              ii+1, I2Ohash_measure(hashtbl)) ;
   }
   if ( msglvl > 3 ) {
      I2Ohash_writeForHumanEye(hashtbl, msgFile) ;
   }
}
if ( msglvl > 0 ) {
   fprintf(msgFile, "\n  %12d %12d %12.3f",
           size, nent, I2Ohash_measure(hashtbl)) ;
}

fprintf(msgFile, "\n") ;
fclose(msgFile) ;

return(1) ; }

/*--------------------------------------------------------------------*/
