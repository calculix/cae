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
 * \brief transform SPCs/MPCs for quadratic slave elements needed for dual mortar contact
 * Author: Saskia Sitzmann
 * @todo: debug error with mid edge nodes
 *
 * @param [in] nboun            number of SPCs
 * @param [in] ndirboun		(i) direction of SPC i 
 * @param [in] nodeboun         (i) node of SPC i
 * @param [in] xboun            (i) value of SPC i
 * @param [in] nmpc		number of mpcs
 * @param [in] ipompc           (i) pointer to nodempc and coeffmpc for MPC i
 * @param [in] nodempc          nodes and directions of MPCs
 * @param [in] coefmpc          coefficients of MPCs
 * @param [in] labmpc		MPCs labels 
 * @param [in] ikboun           sorted dofs idof=8*(node-1)+dir for SPCs
 * @param [in] ilboun           SPC numbers for sorted dofs
 * @param [in] ikmpc 		sorted dofs idof=8*(node-1)+dir for MPCs
 * @param [in] ilmpc		SPC numbers for sorted dofs
 * @param [out] nboun2           number of transformed SPCs
 * @param [out] ndirboun2p	(i) direction of transformed SPC i 
 * @param [out] nodeboun2p      (i) node of transformed SPC i
 * @param [out] xboun2p         (i) value of transformed SPC i
 * @param [out] nmpc2		number of transformed mpcs
 * @param [out] ipompc2p        (i) pointer to nodempc and coeffmpc for transformed MPC i
 * @param [out] nodempc2p       nodes and directions of transformed MPCs
 * @param [out] coefmpc2p       coefficients of transformed MPCs
 * @param [out] labmpc2p	transformed MPCs labels 
 * @param [out] ikboun2p        sorted dofs idof=8*(node-1)+dir for transformed SPCs
 * @param [out] ilboun2p        transformed SPC numbers for sorted dofs
 * @param [out] ikmpc2p 	sorted dofs idof=8*(node-1)+dir for transformed MPCs
 * @param [out] ilmpc2p		transformed SPC numbers for sorted dofs 
 * @param [in] irowtlocinv	field containing row numbers of autlocinv
 * @param [in] jqtlocinv	pointer into field irowtlocinv
 * @param [in] autlocinv	transformation matrix \f$ T^{-1}[p,q]\f$ for slave nodes \f$ p,q \f$  
 * @param [in] nk		number of nodes 
 * @param [out] nk2		number or generated points needed for transformed SPCs 
 * @param [in]  iflagdualquad   flag indicating what mortar contact is used (=1 quad-lin, =2 quad-quad, =3 PG quad-lin, =4 PG quad-quad) 
 * @param [in] ntie		number of contraints 
 * @param [in] tieset           (1,i) name of tie constraint (2,i) dependent surface (3,i) independent surface 
 * @param [in] itiefac 		pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i
 * @param [in] islavsurf	islavsurf(1,i) slaveface i islavsurf(2,i) # integration points generated before looking at face i 
 * @param [in] lakon		(i) label for element i
 * @param [in] ipkon		pointer into field kon...
 * @param [in] kon 		.. for element i storing the connectivity list of elem. in succ. order 
 * @param [in] mt		maximum degrees of freedoms per node  
 * @param [in] memmpc_	size of nodempc/coefmpc 
 * @param [in] nodeforc		point force, node
 * @param [in] ndirforc		point force, dir
 * @param [in] xforc		point force, value
 * @param [in] nforc		number of point forces 
 * @param [out] nodeforc2p	transformed point force, node
 * @param [out] ndirforc2p	transformed point force, dir
 * @param [out] xforc2p		transformed point force, value
 * @param [out] nforc2		number of transformed point forces  
**/
void transformspcsmpcs_quad(ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
			    ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
			    ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
			    ITG *nboun2,ITG **ndirboun2p,ITG **nodeboun2p,double **xboun2p,
			    ITG *nmpc2,ITG **ipompc2p,ITG **nodempc2p,double **coefmpc2p,char **labmpc2p,
			    ITG **ikboun2p,ITG **ilboun2p,ITG **ikmpc2p,ITG **ilmpc2p,
			    ITG *irowtlocinv, ITG *jqtlocinv,double *autlocinv, 
			    ITG *nk,ITG *nk2,ITG *iflagdualquad,
			    ITG *ntie, char *tieset, ITG *itiefac,ITG *islavsurf,
			    char *lakon,ITG *ipkon,ITG *kon,ITG *mt,
			    ITG *memmpc_,
			    ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
			    ITG **nodeforc2p,ITG **ndirforc2p,double **xforc2p,ITG *nforc2){

  ITG i,ii,j,jj,jj2,k,l,ndimboun2,ndimmpc2a,ndimmpc2b,node,nodedep,nodeind,ndir,ndirdep,ndirind,
    ist,index,ifree,idof,id,nhelp,nelems,ifaces,jfaces,idummy,nope,dimnk2,mpcfree2,memmpc_2,
    nopes,nodevdep,nodewdep,nodevind,nodewind,nodes[8],konl[20],ifac,newnode,
    *ndirboun2=NULL, *nodeboun2=NULL, *ipompc2=NULL, *nodempc2=NULL,ndimforc2,
    *ikboun2=NULL, *ilboun2=NULL, *ikmpc2=NULL, *ilmpc2=NULL,
    *nodeforc2=NULL,*ndirforc2=NULL;
  
  
  double alpha,a1,b1,fixed_disp,coeffdep,coeffind,fixed_forc,
    *xboun2=NULL, *coefmpc2=NULL, *xforc2=NULL; 
  
  char *labmpc2=NULL;
  
  ndirboun2=*ndirboun2p; nodeboun2=*nodeboun2p; xboun2=*xboun2p;
  ipompc2=*ipompc2p; nodempc2=*nodempc2p; coefmpc2=*coefmpc2p;
  ikboun2=*ikboun2p; ilboun2=*ilboun2p; ikmpc2=*ikmpc2p; ilmpc2=*ilmpc2p;  
  labmpc2=*labmpc2p;nodeforc2=*nodeforc2p; ndirforc2=*ndirforc2p;
  xforc2=*xforc2p;
  // initialize
  
  dimnk2=10;
  *nk2=0;
  *nboun2=0;
  ndimboun2=*nboun;
  RENEW(ndirboun2,ITG, *nboun);
  RENEW(nodeboun2,ITG, *nboun);
  RENEW(xboun2,double,*nboun);
  *nmpc2=0;
  ndimmpc2a=*nmpc;
  ndimmpc2b=3**memmpc_;
  RENEW(ipompc2,ITG, ndimmpc2a);
  RENEW(nodempc2,ITG, 3*ndimmpc2b);
  RENEW(coefmpc2,double,ndimmpc2b);
  RENEW(ikboun2,ITG, *nboun);
  RENEW(ilboun2,ITG, *nboun);
  RENEW(ikmpc2,ITG, *nmpc);
  RENEW(ilmpc2,ITG, *nmpc);
  RENEW(labmpc2,char,20**nmpc+1);
  *nforc2=0; 
  ndimforc2=*nforc;
  RENEW(nodeforc2,ITG,2**nforc);
  RENEW(ndirforc2,ITG,*nforc);
  RENEW(xforc2,double,*nforc);  
  
  ifree=1;
  a1=0.0;
  b1=0.0;
  if(*iflagdualquad==1 || *iflagdualquad==3){
    alpha=1.0/2.0;
    b1=1.0;
    a1=alpha;
  }else{
    alpha=1.0/5.0;
    b1=(1-2*alpha);
    a1=alpha; 
  }
  /*
    for(i=0;i<*nboun;i++){
    xboun2[i]=xboun[i];
    ndirboun2[i]=ndirboun[i];
    nodeboun2[i]=nodeboun[i];
    ikboun2[i]=ikboun[i];
    ilboun2[i]=ilboun[i];  
    }
    *nboun2=*nboun;
    for(i=0;i<*nmpc;i++){ 
    ist=ipompc[i];
    ipompc2[i]=ipompc[i];
    ilmpc2[i]=ilmpc[i];
    ikmpc2[i]=ikmpc[i];	
    nodempc2[0+(ist-1)*3]=nodempc[0+(ist-1)*3];
    nodempc2[1+(ist-1)*3]=nodempc[1+(ist-1)*3];
    nodempc2[2+(ist-1)*3]=nodempc[2+(ist-1)*3];
    coefmpc2[ist-1]=coefmpc[ist-1];
    index=nodempc[2+(ist-1)*3];
    if(index!=0){
    do{
    nodempc2[0+(index-1)*3]=nodempc[0+(index-1)*3];
    nodempc2[1+(index-1)*3]=nodempc[1+(index-1)*3];
    nodempc2[2+(index-1)*3]=nodempc[2+(index-1)*3];
    coefmpc2[index-1]=coefmpc[index-1];
    index=nodempc[2+(index-1)*3];
    if(index==0){
    break;
    }
    }while(1);
    }
    }
    
    *ndirboun2p=ndirboun2; *nodeboun2p=nodeboun2; *xboun2p=xboun2;
    *ipompc2p=ipompc2; *nodempc2p=nodempc2; *coefmpc2p=coefmpc2;
    *ikboun2p=ikboun2; *ilboun2p=ilboun2; *ikmpc2p=ikmpc2; *ilmpc2p=ilmpc2;
   
    *labmpc2p=labmpc2;
    return;
  */
  
  printf(" transformspcsmpcs_quad: nboun %" ITGFORMAT "\t",*nboun);
  // loop over all SPCs
  for(i=0;i<*nboun;i++){
    fixed_disp=xboun[i];
    ndir=ndirboun[i];
    node=nodeboun[i];
    idof=8*(node-1)+ndir;
    //printf("SPC %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " v %e ifree %" ITGFORMAT " size %" ITGFORMAT "\n",i+1,node,ndir,fixed_disp,*nboun2+1,ndimboun2);
    nhelp=0;
    nodevdep=0;
    nodewdep=0;
    // @todo: set this part active again!
    /*
    if(jqtlocinv[node]-jqtlocinv[node-1]>0){
      //check for mid node
      do{
	for(j=0;j<*ntie;j++){	    	    
	  if(tieset[j*(81*3)+80]=='C'){  
	    for(l=itiefac[2*j];l<=itiefac[2*j+1];l++){
	      ifaces = islavsurf[2*(l-1)+0];
	      nelems = (ITG)(ifaces/10);
	      jfaces = ifaces - nelems*10;
	      FORTRAN(getnumberofnodes,(&nelems,&jfaces,lakon,&nope,&nopes,&idummy)); 
	      for(jj=0;jj<nope;jj++){
		konl[jj]=kon[ipkon[nelems-1]+jj];
	      }
	      for(jj=0;jj<nopes;jj++){
		jj2=jj+1;
		ifac=FORTRAN(getiface,(&jj2,&jfaces,&nope));
		nodes[jj]=konl[ifac-1]; 
	      }
	      ii=-1;
	      for(jj=0;jj<nopes;jj++){
		if(nodes[jj]==node){ii=jj;}
	      }
	      if(ii>-1){
		break;
	      }
	    }
	  }
	}
	break;
      }while(1);
      if((ii>2 && nopes==6)||(ii>3 && nopes==8)){
	nhelp=3;
	if(nopes==6){
	  if(ii==5){
	    nodevdep=nodes[2];
	    nodewdep=nodes[0];
	  }else{
	    nodevdep=nodes[ii-3];
	    nodewdep=nodes[ii-2];	
	  }
	}else{
	  if(ii==7){
	    nodevdep=nodes[3];
	    nodewdep=nodes[0];	    
	  }else{
	    nodevdep=nodes[ii-4];
	    nodewdep=nodes[ii-3];	    
	  }
	}	
      }
    }
    */
    if(nhelp==3 && ndirboun[i]<4){// slave node & mid node for quadratic elements
      //  spc is transformed to mpc
      // create new node
      *nk2=*nk2+1;
      newnode=*nk+*nk2;
      // add spc's
      for(jj=0;jj<4;jj++){
	if(jj==ndir){
	  xboun2[*nboun2]=fixed_disp;
	}else{
	  xboun2[*nboun2]=0.0;   
	}
	ndirboun2[*nboun2]=jj;
	nodeboun2[*nboun2]=newnode;
	idof=8*(newnode-1)+jj;
	FORTRAN(nident,(ikboun2,&idof,nboun2,&id));
	//printf("\t id %" ITGFORMAT " \n",id);
	for( j=*nboun2;j>id;j--){
	  ikboun2[j]=ikboun2[j-1];
	  ilboun2[j]=ilboun2[j-1];
	  //printf("\t\t %" ITGFORMAT " to %" ITGFORMAT "\n",j-1,j);
	}
	//printf("\t place at %" ITGFORMAT "\n",id);
	ikboun2[id]=idof;
	ilboun2[id]=*nboun2+1;    
	*nboun2=*nboun2+1;
	if(*nboun2+1>ndimboun2){
	  ndimboun2=ndimboun2*1.5+1;
	  RENEW(xboun2,double,ndimboun2);
	  RENEW(ndirboun2,ITG, ndimboun2);
	  RENEW(nodeboun2,ITG, ndimboun2);
	  RENEW(ikboun2,ITG, ndimboun2);
	  RENEW(ilboun2,ITG, ndimboun2);
	}
	//printf("SPC2 %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " v %e ifree %" ITGFORMAT " size %" ITGFORMAT "\n",i+1,nodeboun2[*nboun2-1],ndirboun2[*nboun2-1],xboun2[*nboun2-1],*nboun2,ndimboun2);
      }
      // add mpc
      if(*nmpc2+1>ndimmpc2a){
	ndimmpc2a=ndimmpc2a*1.5+1;
	RENEW(ipompc2,ITG, ndimmpc2a);
	RENEW(ikmpc2,ITG, ndimmpc2a);
	RENEW(ilmpc2,ITG, ndimmpc2a);
	RENEW(labmpc2,char,20*ndimmpc2a+1);
      } 
      if(ifree+4>ndimmpc2b){
	ndimmpc2b=ndimmpc2b*1.5+4;
	RENEW(coefmpc2,double,ndimmpc2b);
	RENEW(nodempc2,ITG, 3*ndimmpc2b);
      }	
      for(jj=0;jj<20;jj++){labmpc2[20*(*nmpc2)+jj]=' ';}
      ipompc2[*nmpc2]=ifree;
      nodempc2[0+(ifree-1)*3]=node;
      nodempc2[1+(ifree-1)*3]=ndir;
      nodempc2[2+(ifree-1)*3]=ifree+1;
      coefmpc2[ifree-1]=(b1);
      //printf("\t ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " \n",ipompc2[*nmpc2],nodempc2[0+(ifree-1)*3],nodempc2[1+(ifree-1)*3],coefmpc2[ifree-1],nodempc2[2+(ifree-1)*3]);
      ifree++;
      nodempc2[0+(ifree-1)*3]=nodevdep;
      nodempc2[1+(ifree-1)*3]=ndir;
      nodempc2[2+(ifree-1)*3]=ifree+1;
      coefmpc2[ifree-1]=(a1);
      //printf("\t ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " \n",ipompc2[*nmpc2],nodempc2[0+(ifree-1)*3],nodempc2[1+(ifree-1)*3],coefmpc2[ifree-1],nodempc2[2+(ifree-1)*3]);
      ifree++;
      nodempc2[0+(ifree-1)*3]=nodewdep;
      nodempc2[1+(ifree-1)*3]=ndir;
      nodempc2[2+(ifree-1)*3]=ifree+1;
      coefmpc2[ifree-1]=(a1);
      //printf("\t ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " \n",ipompc2[*nmpc2],nodempc2[0+(ifree-1)*3],nodempc2[1+(ifree-1)*3],coefmpc2[ifree-1],nodempc2[2+(ifree-1)*3]);
      ifree++;       
      nodempc2[0+(ifree-1)*3]=newnode;
      nodempc2[1+(ifree-1)*3]=ndir;
      nodempc2[2+(ifree-1)*3]=0;
      coefmpc2[ifree-1]=-1.0;
      //printf("\t ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " \n",ipompc2[*nmpc2],nodempc2[0+(ifree-1)*3],nodempc2[1+(ifree-1)*3],coefmpc2[ifree-1],nodempc2[2+(ifree-1)*3]);
      ifree++;      
      idof=8*(node-1)+ndir;
      FORTRAN(nident,(ikmpc2,&idof,nmpc2,&id));
      //printf("\t id %" ITGFORMAT " \n",id);
      for( j=*nmpc2;j>id;j--){
	ikmpc2[j]=ikmpc2[j-1];
	ilmpc2[j]=ilmpc2[j-1];
	//printf("\t\t %" ITGFORMAT " to %" ITGFORMAT "\n",j-1,j);
      }
      //printf("\t place at %" ITGFORMAT "\n",id);
      ikmpc2[id]=idof;
      ilmpc2[id]=*nmpc2+1;    
      *nmpc2=*nmpc2+1; 
    }else{ // no slave node or corner slave node
      //nothing to do 
      xboun2[*nboun2]=fixed_disp;
      ndirboun2[*nboun2]=ndir;
      nodeboun2[*nboun2]=node;
      FORTRAN(nident,(ikboun2,&idof,nboun2,&id));
      //printf("\t id %" ITGFORMAT " \n",id);
      for( j=*nboun2;j>id;j--){
	ikboun2[j]=ikboun2[j-1];
	ilboun2[j]=ilboun2[j-1];
	//printf("\t\t %" ITGFORMAT " to %" ITGFORMAT "\n",j-1,j);
      }
      //printf("\t place at %" ITGFORMAT "\n",id);
      ikboun2[id]=idof;
      ilboun2[id]=*nboun2+1;    
      *nboun2=*nboun2+1;
      if(*nboun2+1>ndimboun2){
	ndimboun2=ndimboun2*1.5+1;
	RENEW(xboun2,double,ndimboun2);
	RENEW(ndirboun2,ITG, ndimboun2);
	RENEW(nodeboun2,ITG, ndimboun2);
	RENEW(ikboun2,ITG, ndimboun2);
	RENEW(ilboun2,ITG, ndimboun2);
      }
      //printf("SPC2 %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " v %e ifree %" ITGFORMAT " size %" ITGFORMAT "\n",i+1,nodeboun2[*nboun2-1],ndirboun2[*nboun2-1],xboun2[*nboun2-1],*nboun2,ndimboun2);
    }  
  }
  printf("nboun2 %" ITGFORMAT "\n",*nboun2);
  RENEW(xboun2,double,*nboun2);
  RENEW(ndirboun2,ITG, *nboun2);
  RENEW(nodeboun2,ITG, *nboun2);
  
  // loop over all MPCs
  printf(" transformspcsmpcs_quad: nmpc %" ITGFORMAT "\t",*nmpc);
  for(i=0;i<*nmpc;i++){ 
    ist=ipompc[i];
    nodedep=nodempc[0+(ist-1)*3];
    ndirdep=nodempc[1+(ist-1)*3];
    coeffdep=coefmpc[ist-1];
    nhelp=0;
    /*
    if(jqtlocinv[nodedep]-jqtlocinv[nodedep-1]>0){
      //check for mid node
      ii=-1;
      do{
	for(j=0;j<*ntie;j++){	    	    
	  if(tieset[j*(81*3)+80]=='C'){  
	    for(l=itiefac[2*j];l<=itiefac[2*j+1];l++){
	      ifaces = islavsurf[2*(l-1)+0];
	      nelems = (ITG) (ifaces/10);
	      jfaces = ifaces - nelems*10;
	      FORTRAN(getnumberofnodes,(&nelems,&jfaces,lakon,&nope,&nopes,&idummy)); 
	      for(jj=0;jj<nope;jj++){
		konl[jj]=kon[ipkon[nelems-1]+jj];
	      }
	      //if(nodedep==24707||nodedep==24713||nodedep==24710||nodedep==24715)printf("face %" ITGFORMAT " nodes ",l);
	      for(jj=0;jj<nopes;jj++){
		jj2=jj+1;
		ifac=FORTRAN(getiface,(&jj2,&jfaces,&nope));
		nodes[jj]=konl[ifac-1]; 
	        //if(nodedep==24707||nodedep==24713||nodedep==24710||nodedep==24715)printf(" %" ITGFORMAT " ",nodes[jj]);
	      }
	      //if(nodedep==24707||nodedep==24713||nodedep==24710||nodedep==24715)printf("\n");
	      ii=-1;
	      for(jj=0;jj<nopes;jj++){
		if(nodes[jj]==nodedep){ii=jj;}
	      }
	      //if(nodedep==24707||nodedep==24713||nodedep==24710||nodedep==24715)printf("\t ii %" ITGFORMAT "\n",ii);
	      if(ii>-1){
		break;
	      }
	    }
	    
	  }
	  if(ii>-1){
	    break;
	  }
	}
	break;
      }while(1);
      //if(nodedep==24716||nodedep==24713||nodedep==24710||nodedep==24715)printf("\t ii %" ITGFORMAT "\n",ii);
      if((ii>2 && nopes==6)||(ii>3 && nopes==8)){
	nhelp=3;
	if(nopes==6){
	  if(ii==5){
	    nodevdep=nodes[2];
	    nodewdep=nodes[0];
	  }else{
	    nodevdep=nodes[ii-3];
	    nodewdep=nodes[ii-2];	
	  }
	}else{
	  if(ii==7){
	    nodevdep=nodes[3];
	    nodewdep=nodes[0];	    
	  }else{
	    nodevdep=nodes[ii-4];
	    nodewdep=nodes[ii-3];	    
	  }
	}	
      }
    }
    */
    //if(nodedep==24716||nodedep==24713||nodedep==24710||nodedep==24715){printf("node %" ITGFORMAT " midnode %" ITGFORMAT " \n",nodedep,nhelp);}
    //printf("MPC %" ITGFORMAT " ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " nmpc2 %" ITGFORMAT " size %" ITGFORMAT "\n",i+1,ist,nodedep,ndirdep,coeffdep,nodempc[2+(ist-1)*3],*nmpc2+1,ndimmpc2a);
    if(*nmpc2+1>ndimmpc2a){
      ndimmpc2a=ndimmpc2a*1.5+1;
      RENEW(ipompc2,ITG, ndimmpc2a);
      RENEW(ikmpc2,ITG, ndimmpc2a);
      RENEW(ilmpc2,ITG, ndimmpc2a);
      RENEW(labmpc2,char,20*ndimmpc2a+1);
    }    
    for(jj=0;jj<20;jj++){labmpc2[20*(*nmpc2)+jj]=labmpc[20*i+jj];}
    if(nhelp==3 && ndirdep<4){// slave node & mid node for quadratic elements
      // find adjacent nodes
      if(ifree+3>ndimmpc2b){
	ndimmpc2b=ndimmpc2b*1.5+3;
	RENEW(coefmpc2,double,ndimmpc2b);
	RENEW(nodempc2,ITG, 3*ndimmpc2b);
      }	
      ipompc2[*nmpc2]=ifree;
      nodempc2[0+(ifree-1)*3]=nodedep;
      nodempc2[1+(ifree-1)*3]=ndirdep;
      nodempc2[2+(ifree-1)*3]=ifree+1;
      coefmpc2[ifree-1]=coeffdep*(b1);
      //printf("\t ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " \n",ipompc2[*nmpc2],nodempc2[0+(ifree-1)*3],nodempc2[1+(ifree-1)*3],coefmpc2[ifree-1],nodempc2[2+(ifree-1)*3]);
      ifree++;
      nodempc2[0+(ifree-1)*3]=nodevdep;
      nodempc2[1+(ifree-1)*3]=ndirdep;
      nodempc2[2+(ifree-1)*3]=ifree+1;
      coefmpc2[ifree-1]=coeffdep*(a1);
      //printf("\t ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " \n",ipompc2[*nmpc2],nodempc2[0+(ifree-1)*3],nodempc2[1+(ifree-1)*3],coefmpc2[ifree-1],nodempc2[2+(ifree-1)*3]);
      ifree++;
      nodempc2[0+(ifree-1)*3]=nodewdep;
      nodempc2[1+(ifree-1)*3]=ndirdep;
      nodempc2[2+(ifree-1)*3]=ifree+1;
      coefmpc2[ifree-1]=coeffdep*(a1);
      //printf("\t ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " \n",ipompc2[*nmpc2],nodempc2[0+(ifree-1)*3],nodempc2[1+(ifree-1)*3],coefmpc2[ifree-1],nodempc2[2+(ifree-1)*3]);
      ifree++;
    }else{ // no slave node or corner slave node or mpc not for displacements 
      ipompc2[*nmpc2]=ifree;
      nodempc2[0+(ifree-1)*3]=nodedep;
      nodempc2[1+(ifree-1)*3]=ndirdep;
      nodempc2[2+(ifree-1)*3]=ifree+1;
      coefmpc2[ifree-1]=coeffdep;
      //printf("\t ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " \n",ipompc2[*nmpc2],nodempc2[0+(ifree-1)*3],nodempc2[1+(ifree-1)*3],coefmpc2[ifree-1],nodempc2[2+(ifree-1)*3]);
      ifree++;
      if(ifree>ndimmpc2b){
	ndimmpc2b=ndimmpc2b*1.5+1;
	RENEW(coefmpc2,double,ndimmpc2b);
	RENEW(nodempc2,ITG, 3*ndimmpc2b);
      }
    }
    index=nodempc[2+(ist-1)*3];
    if(index!=0){
      do{
	nodeind=nodempc[0+(index-1)*3];
	ndirind=nodempc[1+(index-1)*3];
	coeffind=coefmpc[index-1];
	nhelp=0;
	/*
	if(jqtlocinv[nodeind]-jqtlocinv[nodeind-1]>0){
	  //check for mid node
	  ii=-1;
	  do{
	    for(j=0;j<*ntie;j++){	    	    
	      if(tieset[j*(81*3)+80]=='C'){  
		for(l=itiefac[2*j];l<=itiefac[2*j+1];l++){
		  ifaces = islavsurf[2*(l-1)+0];
		  nelems = (ITG) (ifaces/10);
		  jfaces = ifaces - nelems*10;
		  //printf("\telem %" ITGFORMAT " iface %" ITGFORMAT " jface %" ITGFORMAT "\n",nelems,ifaces,jfaces);
		  FORTRAN(getnumberofnodes,(&nelems,&jfaces,lakon,&nope,&nopes,&idummy)); 
		  for(jj=0;jj<nope;jj++){
		    konl[jj]=kon[ipkon[nelems-1]+jj];
		  }
		  for(jj=0;jj<nopes;jj++){
		    jj2=jj+1;
		    ifac=FORTRAN(getiface,(&jj2,&jfaces,&nope));
		    nodes[jj]=konl[ifac-1]; 
		    //printf(" nodes(%" ITGFORMAT ",%" ITGFORMAT ")=%" ITGFORMAT " ",ifac,jj,nodes[jj]);
		  }
		  ii=-1;
		  for(jj=0;jj<nopes;jj++){
		    if(nodes[jj]==nodeind){ii=jj;}
		  }
		  //printf("\n\t ii %" ITGFORMAT "\n",ii);
		  if(ii>-1){
		    break;
		  }
		}
	      }
	      if(ii>-1){
		break;
	      } 
	    }
	    break;
	  }while(1);
	  if((ii>2 && nopes==6)||(ii>3 && nopes==8)){
	    nhelp=3;
	    if(nopes==6){
	      if(ii==5){
		nodevind=nodes[2];
		nodewind=nodes[0];
	      }else{
		nodevind=nodes[ii-3];
		nodewind=nodes[ii-2];	
	      }
	    }else{
	      if(ii==7){
		nodevind=nodes[3];
		nodewind=nodes[0];	    
	      }else{
		nodevind=nodes[ii-4];
		nodewind=nodes[ii-3];	    
	      }
	    }	
	  }
	}
	*/
	//printf("\t index %" ITGFORMAT " nodeind %" ITGFORMAT " dirind %" ITGFORMAT " cind %e index %" ITGFORMAT " ifree %" ITGFORMAT " size %" ITGFORMAT "\n",index,nodeind,ndirind,coeffind,nodempc[2+(index-1)*3],ifree,ndimmpc2b);
	if(nhelp==3 && ndirdep<4){// slave node & mid node for quadratic elements
	  // find adjacent nodes	
	  if(ifree+3>ndimmpc2b){
	    ndimmpc2b=ndimmpc2b*1.5+3;
	    RENEW(coefmpc2,double,ndimmpc2b);
	    RENEW(nodempc2,ITG, 3*ndimmpc2b);
	  }
	  nodempc2[0+(ifree-1)*3]=nodeind;
	  nodempc2[1+(ifree-1)*3]=ndirind;
	  nodempc2[2+(ifree-1)*3]=ifree+1;   
	  coefmpc2[ifree-1]=coeffind*(b1);
	  ifree++;
	  //printf("\t ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " \n",ipompc2[*nmpc2],nodempc2[0+(ifree-2)*3],nodempc2[1+(ifree-2)*3],coefmpc2[ifree-2],nodempc2[2+(ifree-2)*3]);
	  nodempc2[0+(ifree-1)*3]=nodevind;
	  nodempc2[1+(ifree-1)*3]=ndirind;
	  nodempc2[2+(ifree-1)*3]=ifree+1;
	  coefmpc2[ifree-1]=coeffind*(a1);
	  ifree++;
	  //printf("\t ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " \n",ipompc2[*nmpc2],nodempc2[0+(ifree-2)*3],nodempc2[1+(ifree-2)*3],coefmpc2[ifree-2],nodempc2[2+(ifree-2)*3]);
	  nodempc2[0+(ifree-1)*3]=nodewind;
	  nodempc2[1+(ifree-1)*3]=ndirind;
	  coefmpc2[ifree-1]=coeffind*(a1);
	  ifree++;	    	    
	}else{// no slave node or corner slave node or mpc not for displacements   
	  nodempc2[0+(ifree-1)*3]=nodeind;
	  nodempc2[1+(ifree-1)*3]=ndirind;
	  coefmpc2[ifree-1]=coeffind;
	  ifree++;
	  if(ifree>ndimmpc2b){
	    ndimmpc2b=ndimmpc2b*1.5+1;
	    RENEW(coefmpc2,double,ndimmpc2b);
	    RENEW(nodempc2,ITG, 3*ndimmpc2b);
	  }	    
	  
	}	
	index=nodempc[2+(index-1)*3];
	if(index==0){
	  nodempc2[2+(ifree-2)*3]=0;
	  //printf("\t ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " \n",ipompc2[*nmpc2],nodempc2[0+(ifree-2)*3],nodempc2[1+(ifree-2)*3],coefmpc2[ifree-2],nodempc2[2+(ifree-2)*3]);
	  break;
	}else{
	  nodempc2[2+(ifree-2)*3]=ifree;
	  //printf("\t ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " \n",ipompc2[*nmpc2],nodempc2[0+(ifree-2)*3],nodempc2[1+(ifree-2)*3],coefmpc2[ifree-2],nodempc2[2+(ifree-2)*3]);
	}
      }while(1);	
    }
    idof=8*(nodedep-1)+ndirdep;
    FORTRAN(nident,(ikmpc2,&idof,nmpc2,&id));
    //printf("\t id %" ITGFORMAT " \n",id);
    for( j=*nmpc2;j>id;j--){
      ikmpc2[j]=ikmpc2[j-1];
      ilmpc2[j]=ilmpc2[j-1];
      //printf("\t\t %" ITGFORMAT " to %" ITGFORMAT "\n",j-1,j);
    }
    //printf("\t place at %" ITGFORMAT "\n",id);
    ikmpc2[id]=idof;
    ilmpc2[id]=*nmpc2+1;    
    *nmpc2=*nmpc2+1;
    if(*nmpc2+1>ndimmpc2a){
      ndimmpc2a=ndimmpc2a*1.5+1;
      RENEW(ipompc2,ITG, ndimmpc2a);
      RENEW(ikmpc2,ITG, ndimmpc2a);
      RENEW(ilmpc2,ITG, ndimmpc2a);
      RENEW(labmpc2,char,20*ndimmpc2a+1);
    } 
  }
  RENEW(ipompc2,ITG, *nmpc2);
  RENEW(ikmpc2,ITG, *nmpc2);
  RENEW(ilmpc2,ITG, *nmpc2);
  RENEW(labmpc2,char,20**nmpc2+1);
  RENEW(coefmpc2,double,(ifree-1));
  RENEW(nodempc2,ITG, 3*(ifree-1));
  //nodempc2[3*(ifree-1)+2]=0;
  printf("nmpc2 %" ITGFORMAT "\n",*nmpc2);
  
  printf(" transformspcsmpcs_quad: nforc %" ITGFORMAT "\t",*nforc);
  // loop over all point forces
  for(i=0;i<*nforc;i++){
    fixed_forc=xforc[i];
    ndir=ndirforc[i];
    node=nodeforc[2*i];
    idof=8*(node-1)+ndir;
    //printf("Node Forc %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " v %e ifree %" ITGFORMAT " size %" ITGFORMAT "\n",i+1,node,ndir,fixed_forc,*nforc2+1,ndimforc2);
    nhelp=0;
    nodevdep=0;
    nodewdep=0;
    // @todo: set this part active again!
    /*
    if(jqtlocinv[node]-jqtlocinv[node-1]>0){
      //check for mid node
      do{
	for(j=0;j<*ntie;j++){	    	    
	  if(tieset[j*(81*3)+80]=='C'){  
	    for(l=itiefac[2*j];l<=itiefac[2*j+1];l++){
	      ifaces = islavsurf[2*(l-1)+0];
	      nelems = (ITG)(ifaces/10);
	      jfaces = ifaces - nelems*10;
	      FORTRAN(getnumberofnodes,(&nelems,&jfaces,lakon,&nope,&nopes,&idummy)); 
	      for(jj=0;jj<nope;jj++){
		konl[jj]=kon[ipkon[nelems-1]+jj];
	      }
	      for(jj=0;jj<nopes;jj++){
		jj2=jj+1;
		ifac=FORTRAN(getiface,(&jj2,&jfaces,&nope));
		nodes[jj]=konl[ifac-1]; 
	      }
	      ii=-1;
	      for(jj=0;jj<nopes;jj++){
		if(nodes[jj]==node){ii=jj;}
	      }
	      if(ii>-1){
		break;
	      }
	    }
	  }
	}
	break;
      }while(1);
      if((ii>2 && nopes==6)||(ii>3 && nopes==8)){
	nhelp=3;
	if(nopes==6){
	  if(ii==5){
	    nodevdep=nodes[2];
	    nodewdep=nodes[0];
	  }else{
	    nodevdep=nodes[ii-3];
	    nodewdep=nodes[ii-2];	
	  }
	}else{
	  if(ii==7){
	    nodevdep=nodes[3];
	    nodewdep=nodes[0];	    
	  }else{
	    nodevdep=nodes[ii-4];
	    nodewdep=nodes[ii-3];	    
	  }
	}	
      }
    }
    */
    if(nhelp==3){
      if(*nforc2+4>ndimforc2){
	ndimforc2=ndimforc2*1.5+3;
	RENEW(xforc2,double,ndimforc2);
	RENEW(ndirforc2,ITG, ndimforc2);
	RENEW(nodeforc2,ITG, 2*ndimforc2);
      }
      //mid node
      xforc2[*nforc2]=fixed_forc*b1;
      ndirforc2[*nforc2]=ndir;
      nodeforc2[2**nforc2]=node;
      //printf("Node Forc %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " v %e ifree %" ITGFORMAT " size %" ITGFORMAT "\n",i+1,nodeforc2[2**nforc2],ndirforc2[*nforc2],xforc2[*nforc2],*nforc2+1,ndimforc2);
      *nforc2=*nforc2+1;     
      //corner node 1
      xforc2[*nforc2]=fixed_forc*a1;
      ndirforc2[*nforc2]=ndir;
      nodeforc2[2**nforc2]=nodevdep;
      //printf("Node Forc %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " v %e ifree %" ITGFORMAT " size %" ITGFORMAT "\n",i+1,nodeforc2[2**nforc2],ndirforc2[*nforc2],xforc2[*nforc2],*nforc2+1,ndimforc2);
      *nforc2=*nforc2+1;
      //corner node 2
      xforc2[*nforc2]=fixed_forc*a1;
      ndirforc2[*nforc2]=ndir;
      nodeforc2[2**nforc2]=nodewdep;
      //printf("Node Forc %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " v %e ifree %" ITGFORMAT " size %" ITGFORMAT "\n",i+1,nodeforc2[2**nforc2],ndirforc2[*nforc2],xforc2[*nforc2],*nforc2+1,ndimforc2);
      *nforc2=*nforc2+1;      
      
    }else{
      //nothing to do 
      xforc2[*nforc2]=fixed_forc;
      ndirforc2[*nforc2]=ndir;
      nodeforc2[2**nforc2]=node;
      *nforc2=*nforc2+1;
      if(*nforc2+1>ndimforc2){
	ndimforc2=ndimforc2*1.5+1;
	RENEW(xforc2,double,ndimforc2);
	RENEW(ndirforc2,ITG, ndimforc2);
	RENEW(nodeforc2,ITG, 2*ndimforc2);
      }      
    }
  }
   printf("nforc2 %" ITGFORMAT "\n\n",*nforc2); 
  
  
  /*for(i=0;i<*nboun2;i++){
    fixed_disp=xboun2[i];
    ndir=ndirboun2[i];
    node=nodeboun2[i];
    idof=8*(node-1)+ndir;
    printf("SPC2 %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " v %e ifree %" ITGFORMAT " size %" ITGFORMAT "\n",i+1,node,ndir,fixed_disp,*nboun2+1,ndimboun2);
    }*/
  /*
    for(i=0;i<*nmpc2;i++){ 
    ist=ipompc2[i];
    nodedep=nodempc2[0+(ist-1)*3];
    ndirdep=nodempc2[1+(ist-1)*3];
    coeffdep=coefmpc2[ist-1];
    printf("MPC2 %" ITGFORMAT " ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " nmpc2 %" ITGFORMAT " size %" ITGFORMAT "\n",i+1,ist,nodedep,ndirdep,coeffdep,nodempc2[2+(ist-1)*3],*nmpc2+1,ndimmpc2a);
    index=nodempc2[2+(ist-1)*3];
    if(index!=0){
    do{
    nodeind=nodempc2[0+(index-1)*3];
    ndirind=nodempc2[1+(index-1)*3];
    coeffind=coefmpc2[index-1];
    printf("\t nodeind %" ITGFORMAT " dirind %" ITGFORMAT " cind %e index %" ITGFORMAT " ifree %" ITGFORMAT " size %" ITGFORMAT "\n",nodeind,ndirind,coeffind,nodempc2[2+(index-1)*3],ifree,ndimmpc2b);
    index=nodempc2[2+(index-1)*3];
    if(index==0){
    break;
    }
    }while(1);
    }
    }
  */
  
  /*for(i=0;i<*nmpc2;i++){ 
    ist=ipompc2[i];
    nodedep=nodempc2[0+(ist-1)*3];
    ndirdep=nodempc2[1+(ist-1)*3];
    coeffdep=coefmpc2[ist-1];
    if(nodedep==580||nodedep==543||nodedep==577||nodedep==565){
      printf("MPC2 %" ITGFORMAT " ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " ilmpc2 %" ITGFORMAT " ikmpc2 %" ITGFORMAT " ilmpc %" ITGFORMAT " ikmpc %" ITGFORMAT "\n",i+1,ist,nodedep,ndirdep,coeffdep,nodempc2[2+(ist-1)*3],ilmpc2[i],ikmpc2[i],ilmpc[i],ikmpc[i]);
    }
    index=nodempc2[2+(ist-1)*3];
    if(index!=0){
      do{
	nodeind=nodempc2[0+(index-1)*3];
	ndirind=nodempc2[1+(index-1)*3];
	coeffind=coefmpc2[index-1];
	if(nodedep==580||nodedep==543||nodedep==577||nodedep==565){
	  printf("\t nodeind %" ITGFORMAT " dirind %" ITGFORMAT " cind %e index %" ITGFORMAT " ifree %" ITGFORMAT " size %" ITGFORMAT "\n",nodeind,ndirind,coeffind,nodempc2[2+(index-1)*3],ifree,ndimmpc2b);
	}
	index=nodempc2[2+(index-1)*3];
	if(index==0){
	  break;
	}
      }while(1);
    }
  }*/
  //FORTRAN(stop,());
  
  //decascade mpcs
  mpcfree2=0;
  memmpc_2=ifree-1;
  decascade_mortar(nmpc2,ipompc2,&nodempc2,&coefmpc2,
	         ikmpc2,ilmpc2,&memmpc_2,&mpcfree2);
 //printf("\n\n");
 //i=2;
 
  /*for(i=0;i<*nmpc2;i++){ 
    ist=ipompc2[i];
    nodedep=nodempc2[0+(ist-1)*3];
    ndirdep=nodempc2[1+(ist-1)*3];
    coeffdep=coefmpc2[ist-1];
    if(nodedep==356||nodedep==381||nodedep==380||nodedep==565){
      printf("MPC3 %" ITGFORMAT " ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " nmpc2 %" ITGFORMAT " size %" ITGFORMAT "\n",i+1,ist,nodedep,ndirdep,coeffdep,nodempc2[2+(ist-1)*3],*nmpc2+1,ndimmpc2a);
    
    index=nodempc2[2+(ist-1)*3];
    if(index!=0){
      do{
	nodeind=nodempc2[0+(index-1)*3];
	ndirind=nodempc2[1+(index-1)*3];
	coeffind=coefmpc2[index-1];
	//if(nodedep==356||nodedep==381||nodedep==380||nodedep==565){
	  printf("\t nodeind %" ITGFORMAT " dirind %" ITGFORMAT " cind %e index %" ITGFORMAT " ifree %" ITGFORMAT " size %" ITGFORMAT "\n",nodeind,ndirind,coeffind,nodempc2[2+(index-1)*3],ifree,ndimmpc2b);
	//}
	index=nodempc2[2+(index-1)*3];
	if(index==0){
	  break;
	}
      }while(1);
    }
    }
  }*/
  //printf("\n\n");
  /*for(i=0;i<*nmpc;i++){ 
    //i=2;
    ist=ipompc[i];
    nodedep=nodempc[0+(ist-1)*3];
    ndirdep=nodempc[1+(ist-1)*3];
    coeffdep=coefmpc[ist-1];
    if(nodedep==356||nodedep==381||nodedep==380||nodedep==565){
      printf("MPC %" ITGFORMAT " ist %" ITGFORMAT " node %" ITGFORMAT " dir %" ITGFORMAT " c %e index %" ITGFORMAT " labmpc %c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c\n",i+1,ist,nodedep,ndirdep,coeffdep,nodempc[2+(ist-1)*3],labmpc[20*i+0],labmpc[20*i+1],labmpc[20*i+2],labmpc[20*i+3],labmpc[20*i+4],labmpc[20*i+5],labmpc[20*i+6],labmpc[20*i+7],labmpc[20*i+8],labmpc[20*i+9],labmpc[20*i+10],labmpc[20*i+11],labmpc[20*i+12],labmpc[20*i+13],labmpc[20*i+14],labmpc[20*i+15],labmpc[20*i+16],labmpc[20*i+17],labmpc[20*i+18],labmpc[20*i+19]);
    }
    index=nodempc[2+(ist-1)*3];
    if(index!=0){
      do{
	if(nodedep==356||nodedep==381||nodedep==380||nodedep==565){
	  printf("\t index %" ITGFORMAT " nodeind %" ITGFORMAT " dirind %" ITGFORMAT " cind %e index %" ITGFORMAT " \n",index,nodempc[0+(index-1)*3],nodempc[1+(index-1)*3],coefmpc[index-1],nodempc[2+(index-1)*3]);
	}
	index=nodempc[2+(index-1)*3];
	if(index==0){
	  break;
	}
      }while(1);
    }
  }*/
    /*for(i=0;i<*nboun2;i++){ 
     printf("nodeboun %" ITGFORMAT " dirbounn %" ITGFORMAT " xboun %e nodeboun2 %" ITGFORMAT " dirbounn2 %" ITGFORMAT " xboun2 %e \n",nodeboun[i],ndirboun[i],xboun[i],nodeboun2[i],ndirboun2[i],xboun2[i]); 
    }*/
    //printf("nk2 %" ITGFORMAT " \n",*nk2);
  *ndirboun2p=ndirboun2; *nodeboun2p=nodeboun2; *xboun2p=xboun2;
  *ipompc2p=ipompc2; *nodempc2p=nodempc2; *coefmpc2p=coefmpc2;
  *ikboun2p=ikboun2; *ilboun2p=ilboun2; *ikmpc2p=ikmpc2; *ilmpc2p=ilmpc2;
  *labmpc2p=labmpc2;
  *nodeforc2p=nodeforc2;*ndirforc2p=ndirforc2;
  *xforc2p=xforc2;  
  //FORTRAN(stop,());
  return;
}
