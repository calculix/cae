/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

#include "readfrd.h"


/* liefert double aus string von position a bis b */
double stof(char *string, int a, int b)
{
  register int    n, i;
  static char  puffer[MAX_LINE_LENGTH];

  n=-1;
  for (i=a-1; i<b; i++)
    {
    n++;
    if ((i>=MAX_LINE_LENGTH)||(n>=MAX_LINE_LENGTH)) break;
    puffer[n] = string[i];
  }
  puffer[n+1] = '\0';
  return( atof( puffer ) );
}

