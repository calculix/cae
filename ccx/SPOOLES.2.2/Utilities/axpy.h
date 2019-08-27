/*  axpy.h  */

/*
   -----------------------------------------------------------------
   y0[] = y0[] + alpha[0] * x0[] + alpha[1] * x1[] + alpha[2] * x2[]
   y1[] = y1[] + alpha[3] * x0[] + alpha[4] * x1[] + alpha[5] * x2[]
   y2[] = y2[] + alpha[6] * x0[] + alpha[7] * x1[] + alpha[8] * x2[]

   created -- 98dec10, cca
   -----------------------------------------------------------------
*/
void
DVaxpy33 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   alpha[],
   double   x0[],
   double   x1[],
   double   x2[]
) ;
/*
   -----------------------------------------------
   y0[] = y0[] + alpha[0] * x0[] + alpha[1] * x1[] 
   y1[] = y1[] + alpha[2] * x0[] + alpha[3] * x1[] 
   y2[] = y2[] + alpha[4] * x0[] + alpha[5] * x1[] 

   created -- 98dec10, cca
   -----------------------------------------------
*/
void
DVaxpy32 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   alpha[],
   double   x0[],
   double   x1[]
) ;
/*
   -----------------------------
   y0[] = y0[] + alpha[0] * x0[]
   y1[] = y1[] + alpha[1] * x0[]
   y2[] = y2[] + alpha[2] * x0[]

   created -- 98dec10, cca
   -----------------------------
*/
void
DVaxpy31 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   alpha[],
   double   x0[]
) ;
/*
   -----------------------------------------------------------------
   y0[] = y0[] + alpha[0] * x0[] + alpha[1] * x1[] + alpha[2] * x2[]
   y1[] = y1[] + alpha[3] * x0[] + alpha[4] * x1[] + alpha[5] * x2[]

   created -- 98dec10, cca
   -----------------------------------------------------------------
*/
void
DVaxpy23 (
   int      n,
   double   y0[],
   double   y1[],
   double   alpha[],
   double   x0[],
   double   x1[],
   double   x2[]
) ;
/*
   -----------------------------------------------
   y0[] = y0[] + alpha[0] * x0[] + alpha[1] * x1[]
   y1[] = y1[] + alpha[2] * x0[] + alpha[3] * x1[]

   created -- 98dec10, cca
   -----------------------------------------------
*/
void
DVaxpy22 (
   int      n,
   double   y0[],
   double   y1[],
   double   alpha[],
   double   x0[],
   double   x1[]
) ;
/*
   -----------------------------
   y0[] = y0[] + alpha[0] * x0[]
   y1[] = y1[] + alpha[1] * x0[]

   created -- 98dec10, cca
   -----------------------------
*/
void
DVaxpy21 (
   int      n,
   double   y0[],
   double   y1[],
   double   alpha[],
   double   x0[]
) ;
/*
   -----------------------------------------------------------------
   y0[] = y0[] + alpha[0] * x0[] + alpha[1] * x1[] + alpha[2] * x2[]

   created -- 98dec10, cca
   -----------------------------------------------------------------
*/
void
DVaxpy13 (
   int      n,
   double   y0[],
   double   alpha[],
   double   x0[],
   double   x1[],
   double   x2[]
) ;
/*
   -----------------------------------------------
   y0[] = y0[] + alpha[0] * x0[] + alpha[1] * x1[]

   created -- 98dec10, cca
   -----------------------------------------------
*/
void
DVaxpy12 (
   int      n,
   double   y0[],
   double   alpha[],
   double   x0[],
   double   x1[]
) ;
/*
   -----------------------------
   y0[] = y0[] + alpha[0] * x0[]

   created -- 98dec10, cca
   -----------------------------
*/
void
DVaxpy11 (
   int      n,
   double   y0[],
   double   alpha[],
   double   x0[]
) ;
/*
   -----------------------------------------------------------------
   y0[] = y0[] + alpha[0:1] * x0[]
               + alpha[2:3] * x1[] + alpha[4:5] * x2[]
   y1[] = y1[] + alpha[6:7] * x0[]
               + alpha[8:9] * x1[] + alpha[10:11] * x2[]
   y2[] = y2[] + alpha[12:13] * x0[]
               + alpha[14:15] * x1[] + alpha[16:17] * x2[]

   created -- 98dec10, cca
   -----------------------------------------------------------------
*/
void
ZVaxpy33 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   alpha[],
   double   x0[],
   double   x1[],
   double   x2[]
) ;
/*
   -----------------------------------------------------
   y0[] = y0[] + alpha[0:1] * x0[] + alpha[2:3] * x1[]
   y1[] = y1[] + alpha[4:5] * x0[] + alpha[6:7] * x1[]
   y2[] = y2[] + alpha[8:9] * x0[] + alpha[10:11] * x1[]

   created -- 98dec10, cca
   -----------------------------------------------------
*/
void
ZVaxpy32 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   alpha[],
   double   x0[],
   double   x1[]
) ;
/*
   -------------------------------
   y0[] = y0[] + alpha[0:1] * x0[]
   y1[] = y1[] + alpha[2:3] * x0[]
   y2[] = y2[] + alpha[4:5] * x0[]

   created -- 98dec10, cca
   -------------------------------
*/
void
ZVaxpy31 (
   int      n,
   double   y0[],
   double   y1[],
   double   y2[],
   double   alpha[],
   double   x0[]
) ;
/*
   -----------------------------------------------------------------
   y0[] = y0[] + alpha[0:1] * x0[]
               + alpha[2:3] * x1[] + alpha[4:5] * x2[]
   y1[] = y1[] + alpha[6:7] * x0[]
               + alpha[8:9] * x1[] + alpha[10:11] * x2[]

   created -- 98dec10, cca
   -----------------------------------------------------------------
*/
void
ZVaxpy23 (
   int      n,
   double   y0[],
   double   y1[],
   double   alpha[],
   double   x0[],
   double   x1[],
   double   x2[]
) ;
/*
   -----------------------------------------------------
   y0[] = y0[] + alpha[0:1] * x0[] + alpha[2:3] * x1[]
   y1[] = y1[] + alpha[4:5] * x0[] + alpha[6:7] * x1[]

   created -- 98dec10, cca
   -----------------------------------------------------
*/
void
ZVaxpy22 (
   int      n,
   double   y0[],
   double   y1[],
   double   alpha[],
   double   x0[],
   double   x1[]
) ;
/*
   -------------------------------
   y0[] = y0[] + alpha[0:1] * x0[]
   y1[] = y1[] + alpha[2:3] * x0[]

   created -- 98dec10, cca
   -------------------------------
*/
void
ZVaxpy21 (
   int      n,
   double   y0[],
   double   y1[],
   double   alpha[],
   double   x0[]
) ;
/*
   ---------------------------------------------------
   y0[] = y0[] + alpha[0:1] * x0[]
               + alpha[2:3] * x1[] + alpha[4:5] * x2[]

   created -- 98dec10, cca
   ---------------------------------------------------
*/
void
ZVaxpy13 (
   int      n,
   double   y0[],
   double   alpha[],
   double   x0[],
   double   x1[],
   double   x2[]
) ;
/*
   -----------------------------------------------------
   y0[] = y0[] + alpha[0:1] * x0[] + alpha[2:3] * x1[]

   created -- 98dec10, cca
   -----------------------------------------------------
*/
void
ZVaxpy12 (
   int      n,
   double   y0[],
   double   alpha[],
   double   x0[],
   double   x1[]
) ;
/*
   -------------------------------
   y0[] = y0[] + alpha[0:1] * x0[]

   created -- 98dec10, cca
   -------------------------------
*/
void
ZVaxpy11 (
   int      n,
   double   y0[],
   double   alpha[],
   double   x0[]
) ;
