/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2019 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "CalculiX.h"
#include "mortar.h"
/** 
* \brief function preparing for mortar output in frd-file
* Author: Saskia Sitzmann
 * @param [in,out] ne		number of elements
 * @param [in] nslavs		number of slave nodes
 * @param [in] mi		(1) max # of integration points per element (2) max degree of freedom per element
 * @param [in] nk		number of nodes 
 * @param [in,out] nkon		size of kon
 * @param [in,out] stxp		usually: (1:6,k,i) buckling stress of element i in integration point k here: used to transmit contact variables
  * @param [in] cdisp		mortar contact output CDISP, CSTRESS
 * @param [in,out] fn		internal forces
 * @param [in] cfs		contat forces
 * @param [in] cfm		not used any more
**/
void mortar_prefrd(ITG *ne,ITG *nslavs, ITG *mi, ITG *nk, ITG *nkon, 
		   double **stxp, double *cdisp,
		   double *fn, double *cfs, double *cfm){

  ITG i,k,l,mt=mi[1]+1;
  
  double *stx=NULL;
  
  /* save CDISP,CPRESS in stx */
  stx=*stxp;	      
  RENEW(stx,double,6*mi[0]*(*ne+*nslavs));      
  for(k=0;k<*nslavs;k++){
    for(l=0;l<6;l++){
      stx[6*mi[0]*(*ne+k)+l]=cdisp[6*k+l];
      //printf("node %d deg %d v %e",k,l,cdisp[6*k+l]);
    }     
  }
  /* add slave nodes to ne and nkon */
  *ne+=*nslavs;
  *nkon+=*nslavs;    
  /* add contact forces to reaction force */
  for(i=0;i<mt**nk;i++){
    fn[i]=fn[i]+cfs[i];      
  }
  *stxp=stx;
  
  return;
}
