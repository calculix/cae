/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2019 Guido Dhondt                     */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                                                                       */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"
#include "mortar.h"
/** \brief inserts a new nonzero matrix position into the sorted sparse matrix data structure
 *
 * @param [in,out] irowp		field saving the row numbers
 * @param [in] i1   		row number (FORTRAN convention)
 * @param [in] i2   		column number (FORTRAN convention)
 * @param [in,out] ifree	position of next free space in row,imast1,bd
 * @param [in,out] nzs_		size of row,imast1,bd
 * @param [in] contribution	value to be saved in bd
 * @param [in,out] bdp		values of sparse matrix  
**/
void insertas_ws(ITG **irowp, ITG *i1,
		 ITG *i2, ITG *ifree, ITG *nzs_, double *contribution, double **bdp){
  
  /*   inserts a new nonzero matrix position into the data structure
       the structure is not assumed to be symmetric 
       i1: row number (FORTRAN convention) 
       i2: column number (FORTRAN convention) */
  
  ITG  idof1,idof2,*irow=NULL,*mast1=NULL,nzs_old,i;
  double *bd=NULL;
  
  irow=*irowp;   
  bd=*bdp;
  
  idof1 = *i1;
  idof2 = *i2;
  
  if(*ifree>*nzs_){
    //      printf("Insertas RENEW ifree = %d,nzs = %d\n",*ifree,*nzs_);
    //      *nzs_=(ITG )(1.1**nzs_);
    nzs_old=*nzs_;
    *nzs_=(ITG)(1.5**nzs_+10);
    RENEW(irow,ITG ,*nzs_);
    for(i=nzs_old;i<*nzs_;i++){irow[i]=0;}
    if (irow==NULL) printf("WARNING !!!!\n");
    RENEW(bd,double,*nzs_);
    for(i=nzs_old;i<*nzs_;i++){bd[i]=0.;}
  }
  irow[*ifree-1]=idof1;
  bd[*ifree-1]=*contribution;
  ++*ifree;
  //    printf("ifree %d\n",*ifree);
  
  *irowp=irow;
  *bdp=bd;
  
  return;
  
}
