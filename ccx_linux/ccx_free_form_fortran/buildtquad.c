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
#include <time.h>
#include "CalculiX.h"
#include "mortar.h"
/**
 * Calculates the transformation matrices \f$ T \f$ and \f$ T^{-1} \f$ for quad-quad or quad-lin method
 * see phd-thesis Sitzmann Chapter 4.1
 * Author: Saskia Sitzmann
 * 
 * @param [in] ntie		number of contraints
 * @param [in] ipkon		pointer into field kon...
 * @param [in] kon 		.. for element i storing the connectivity list of elem. in succ. order 
 * @param [in] nk 		number of nodes
 * @param [in] lakon		(i) label for element i 
 * @param [in] nslavnode	(i)pointer into field isalvnode for contact tie i 
 * @param [in] itiefac		pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i 
 * @param [in] tieset           (i) name of tie i 
 * @param [in] islavnode	field storing the nodes of the slave surface
 * @param [in] islavsurf	islavsurf(1,i) slaveface i islavsurf(2,i) # integration points generated before looking at face i  
 * @param [out] irowtlocp		field containing row numbers of autloc
 * @param [out] jqtloc	        pointer into field irowtloc
 * @param [out] autlocp		transformation matrix \f$ T[p,q]\f$ for slave nodes \f$ p,q \f$ 
 * @param [out] irowtlocinvp	field containing row numbers of autlocinv
 * @param [out] jqtlocinv	pointer into field irowtlocinv
 * @param [out] autlocinvp	transformation matrix \f$ T^{-1}[p,q]\f$ for slave nodes \f$ p,q \f$ 
 * @param [in]  iflagdualquad   flag indicating what mortar contact is used (=1 quad-lin, =2 quad-quad, =3 PG quad-lin, =4 PG quad-quad)
*/
void buildtquad(ITG *ntie,ITG *ipkon,ITG *kon,ITG *nk,char *lakon,
		ITG *nslavnode,ITG *itiefac,char *tieset,
		ITG *islavnode,ITG *islavsurf,
		ITG **irowtlocp,ITG *jqtloc,double **autlocp,
		ITG **irowtlocinvp,ITG *jqtlocinv,double **autlocinvp,
		ITG *iflagdualquad){  
  
  ITG i,j,l,nodesf,nodem,istart,icounter,ndim,ifree,ifree2,
    nzstloc,nzstlocinv,*idcontr1=NULL,*idcontr2=NULL,debug,
    *mast1=NULL,*mast2=NULL,*irowtloc=NULL,*irowtlocinv=NULL;
  
  double contribution,*dcontr=NULL,*autloc=NULL,*autlocinv=NULL;
  
  irowtloc=*irowtlocp; autloc=*autlocp;
  irowtlocinv=*irowtlocinvp; autlocinv=*autlocinvp;
  
  debug=0;
  
  /** built T and T^-1 **/
  
  //debug=1;
  nzstloc=3*nslavnode[*ntie];
  nzstlocinv=3*nslavnode[*ntie];
  NNEW(mast1,ITG,nzstloc);
  RENEW(autloc,double,nzstloc);
  RENEW(irowtloc,ITG,nzstloc);
  NNEW(mast2,ITG,nzstlocinv);
  RENEW(autlocinv,double,nzstloc);
  RENEW(irowtlocinv,ITG,nzstloc);	
  ifree=1;ifree2=1;
  NNEW(dcontr,double,16);
  NNEW(idcontr1,ITG,16);
  NNEW(idcontr2,ITG,16);
  for(i=0;i<*ntie;i++){	    
    if(tieset[i*(81*3)+80]=='C'){				
      for(l=itiefac[2*i];l<=itiefac[2*i+1];l++){
	if(debug==1)printf("bdfill face %" ITGFORMAT "  \n",l);
	if(*iflagdualquad==2 || *iflagdualquad==4){
	  FORTRAN(createtele,(ipkon,kon,lakon,islavsurf,
			      dcontr,idcontr1,idcontr2,&icounter,&l));
	}else{
	  FORTRAN(createtele_lin,(ipkon,kon,lakon,islavsurf,
				  dcontr,idcontr1,idcontr2,&icounter,&l));	     
	}
	
	for(j=0;j<icounter;j++){
	  contribution=dcontr[j];
	  nodesf=idcontr1[j];				
	  nodem=idcontr2[j];				
	  if(debug==1)printf("\tbdfill: T nodesf %" ITGFORMAT " nodem %"
			     ITGFORMAT " c %e \n",nodesf,nodem,contribution);                                			        
	  /// T
	  insertas(&irowtloc,&mast1,&nodesf,&nodem,&ifree,&nzstloc,
		   &contribution,&autloc);				
	}
	if(*iflagdualquad==2 || *iflagdualquad==4){
	  FORTRAN(createteleinv,(ipkon,kon,lakon,islavsurf,
				 dcontr,idcontr1,idcontr2,&icounter,&l));
	}else{
	  FORTRAN(createteleinv_lin,(ipkon,kon,lakon,islavsurf,
				     dcontr,idcontr1,idcontr2,&icounter,&l));
	}
	for(j=0;j<icounter;j++){
	  contribution=dcontr[j];
	  nodesf=idcontr1[j];				
	  nodem=idcontr2[j];				
	  if(debug==1)printf("\tbdfill: Tinv nodesf %" ITGFORMAT " nodem %"
			     ITGFORMAT " c %e \n",nodesf,nodem,contribution);  
	  /// T^-1
	  insertas(&irowtlocinv,&mast2,&nodesf,&nodem,&ifree2,&nzstlocinv,
		   &contribution,&autlocinv);
	}	      
      }
    }
  }
  SFREE(dcontr);SFREE(idcontr1);SFREE(idcontr2);
    
  nzstloc=ifree-1;
  printf(" buildquad: \tsize autloc %" ITGFORMAT " \n",ifree);
  ndim=*nk;
  matrixsort(autloc,mast1,irowtloc,jqtloc,&nzstloc,&ndim);
  
  /*Calculation of the real contribution*/
  
  icounter=0;
  for (i=0;i<*nk;i++){
    if(jqtloc[i]!=jqtloc[i+1]){
      //printf("\t bdfill: column %" ITGFORMAT " row %" ITGFORMAT " v %e\n",i,irowtloc[jqtloc[i]-1],autloc[jqtloc[i]-1]);
      irowtloc[icounter]=irowtloc[jqtloc[i]-1];
      autloc[icounter]=autloc[jqtloc[i]-1];
      icounter++;
      istart=icounter;
      for (j=jqtloc[i];j<jqtloc[i+1]-1;j++){
	if (irowtloc[j]==irowtloc[icounter-1]){
	  //printf("\t bdfill:same column %" ITGFORMAT " row %" ITGFORMAT " v %e\n",i,irowtloc[j],autloc[j]);
	  autloc[icounter-1]=autloc[j];   
	}else{
	  //printf("\t bdfill: column %" ITGFORMAT " row %" ITGFORMAT " v %e\n",i,irowtloc[j],autloc[j]);
	  irowtloc[icounter]=irowtloc[j];
	  autloc[icounter]=autloc[j];
	  icounter++;
	  }
      }
    }else{ istart=icounter+1;}   
    jqtloc[i]=istart;
  }
  jqtloc[*nk]=icounter+1; 
  printf(" buildquad: \tsize autloc %" ITGFORMAT " \n",icounter);
  RENEW(irowtloc,ITG,icounter+1);
  RENEW(autloc,double,icounter+1); 
  SFREE(mast1);	
  //FORTRAN(stop,());
  
  
  nzstlocinv=ifree2-1; 
  ndim=*nk;
  matrixsort(autlocinv,mast2,irowtlocinv,jqtlocinv,&nzstlocinv,&ndim);  
 
  /*Calculation of the real contribution*/
  
  icounter=0;
  for (i=0;i<*nk;i++){
    if(jqtlocinv[i]!=jqtlocinv[i+1]){
      //printf("\t bdfill: column %" ITGFORMAT " row %" ITGFORMAT " v %e\n",i,irowtlocinv[jqtlocinv[i]-1],autlocinv[jqtlocinv[i]-1]);
      irowtlocinv[icounter]=irowtlocinv[jqtlocinv[i]-1];
      autlocinv[icounter]=autlocinv[jqtlocinv[i]-1];
      icounter++;
      istart=icounter;
      for (j=jqtlocinv[i];j<jqtlocinv[i+1]-1;j++){
	if (irowtlocinv[j]==irowtlocinv[icounter-1]){
	  //printf("\t bdfill:same column %" ITGFORMAT " row %" ITGFORMAT " v %e\n",i,irowtlocinv[j],autlocinv[j]);
	  autlocinv[icounter-1]=autlocinv[j];   
	}else{
	  //printf("\t bdfill: column %" ITGFORMAT " row %" ITGFORMAT " v %e\n",i,irowtlocinv[j],autlocinv[j]);
	  irowtlocinv[icounter]=irowtlocinv[j];
	  autlocinv[icounter]=autlocinv[j];
	  icounter++;
	}
      }
    }else{ istart=icounter+1;}   
    jqtlocinv[i]=istart;
  }
  jqtlocinv[*nk]=icounter+1;
  
  printf(" buildquad: \tsize autlocinv %" ITGFORMAT " \n\n",icounter);
  RENEW(irowtlocinv,ITG,icounter+1);
  RENEW(autlocinv,double,icounter+1); 
  SFREE(mast2);
  
  *irowtlocp=irowtloc;*autlocp=autloc;
  *irowtlocinvp=irowtlocinv;*autlocinvp=autlocinv;
  //FORTRAN(stop,());
  return;
}	    

  
