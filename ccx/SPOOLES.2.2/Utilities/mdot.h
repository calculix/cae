/*  mdot.h  */

/*--------------------------------------------------------------------*/
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
      sums[2] = row0[*] * col2[*]
      sums[3] = row1[*] * col0[*]
      sums[4] = row1[*] * col1[*]
      sums[5] = row1[*] * col2[*]
      sums[6] = row2[*] * col0[*]
      sums[7] = row2[*] * col1[*]
      sums[8] = row2[*] * col2[*]

   created -- 96oct20, cca
   -----------------------------------------
*/
void
mdot3x3 (
   double    sums[],
   int       n,
   double    row0[],
   double    row1[],
   double    row2[],
   double    col0[],
   double    col1[],
   double    col2[]
) ;
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
      sums[2] = row0[*] * col2[*]
      sums[3] = row1[*] * col0[*]
      sums[4] = row1[*] * col1[*]
      sums[5] = row1[*] * col2[*]

   created -- 96oct20, cca
   -----------------------------------------
*/
void
mdot2x3 (
   double    sums[],
   int       n,
   double    row0[],
   double    row1[],
   double    col0[],
   double    col1[],
   double    col2[]
) ;
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
      sums[2] = row0[*] * col2[*]

   created -- 96oct20, cca
   -----------------------------------------
*/
void
mdot1x3 (
   double    sums[],
   int       n,
   double    row0[],
   double    col0[],
   double    col1[],
   double    col2[]
) ;
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
      sums[2] = row1[*] * col0[*]
      sums[3] = row1[*] * col1[*]
      sums[4] = row2[*] * col0[*]
      sums[5] = row2[*] * col1[*]

   created -- 96oct20, cca
   -----------------------------------------
*/
void
mdot3x2 (
   double    sums[],
   int       n,
   double    row0[],
   double    row1[],
   double    row2[],
   double    col0[],
   double    col1[]
) ;
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]
      sums[3] = row1[*] * col0[*]
      sums[4] = row1[*] * col1[*]

   created -- 96oct20, cca
   -----------------------------------------
*/
void
mdot2x2 (
   double    sums[],
   int       n,
   double    row0[],
   double    row1[],
   double    col0[],
   double    col1[]
) ;
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[1] = row0[*] * col1[*]

   created -- 96oct20, cca
   -----------------------------------------
*/
void
mdot1x2 (
   double    sums[],
   int       n,
   double    row0[],
   double    col0[],
   double    col1[]
) ;
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[1] = row1[*] * col0[*]
      sums[2] = row2[*] * col0[*]

   created -- 96oct20, cca
   -----------------------------------------
*/
void
mdot3x1 (
   double    sums[],
   int       n,
   double    row0[],
   double    row1[],
   double    row2[],
   double    col0[]
) ;
/*
   -----------------------------------------
   purpose -- compute a multiple dot product

      sums[0] = row0[*] * col0[*]
      sums[1] = row1[*] * col0[*]

   created -- 96oct20, cca
   -----------------------------------------
*/
void
mdot2x1 (
   double    sums[],
   int       n,
   double    row0[],
   double    row1[],
   double    col0[]
) ;
/*
   -----------------------------------------
   purpose -- compute a dot product

      sums[0] = row0[*] * col0[*]

   created -- 96oct20, cca
   -----------------------------------------
*/
void
mdot1x1 (
   double    sums[],
   int       n,
   double    row0[],
   double    col0[]
) ;
/*--------------------------------------------------------------------*/
