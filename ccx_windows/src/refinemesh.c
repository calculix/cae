/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2020 Guido Dhondt                          */

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

void refinemesh(ITG *nk,ITG *ne,double *co,ITG *ipkon,ITG *kon,
		double *v,double *veold,double *stn,double *een,
		double *emn,double *epn,double *enern,double *qfn,
		double *errn,char *filab,ITG *mi,char *lakon,
                char *jobnamec,ITG *istartset,ITG *iendset,
		ITG *ialset,char *set,ITG *nset,char *matname,
		ITG *ithermal,char *output,ITG *nmat){

  /* refinement of tetrahedral meshes as a function of a field
     variable */

  ITG nktet_,netet_,*kontet=NULL,*ifatet=NULL,*ifac=NULL,*itetfa=NULL,
    *ipofa=NULL,*ipoeln=NULL,*ieln=NULL,ifreetet=1,ifreefa=1,ielem,
    ifreeln=1,*iedg=NULL,*ipoed=NULL,*ipoeled=NULL,*ieled=NULL,
    *iedtet=NULL,ifreeed=1,ifreele=1,*integerglob=NULL,*iexternfa=NULL,
    *iexternedg=NULL,*iexternnode=NULL,*n=NULL,nnewnodes,newsize,
    *ibasenewnodes=NULL,i,k,*irand1=NULL,*irand2=NULL,iext,iedge,
    kflag,node,ibase,nktet,*nh=NULL,ij,nloop,*iedgnewnodes=NULL,
    *kontetor=NULL,iquad=0,*iedgext=NULL,nexternedg,*iedgmid=NULL,
    *ifacext=NULL,nexternfa,*iexternfaor=NULL,*itreated=NULL,
    *ialsete=NULL,nexternel,*iedgextfa=NULL,*ifacexted=NULL,
    *ilist=NULL,*isharp=NULL,*idimsh=NULL,ifreenn=1,*iponn=NULL,
    *inn=NULL,*n1newnodes=NULL,*n2newnodes=NULL,*jfix=NULL,*number=NULL,
    *iparentel=NULL;

  //  unsigned ITG unsiint,seed=184389;
  unsigned long unsiint,seed=184389;

  double *cotet=NULL,*planfa=NULL,*bc=NULL,*cg=NULL, *h=NULL,dcg,
    *doubleglob=NULL,*r=NULL,*d=NULL,*conewnodes=NULL,dmin,*quality=NULL,
    randomnumber,*conewnodesorig=NULL,cotetorig[3],height,
    *hnewnodes=NULL;

  nktet=*nk;
  nktet_=*nk+1000;
  netet_=2**ne+1000;
    
  NNEW(cotet,double,3*nktet_);
  memcpy(&cotet[0],&co[0],sizeof(double)*3**nk);
  NNEW(h,double,nktet_);
  NNEW(jfix,ITG,nktet_);
  NNEW(kontet,ITG,4*netet_);
  NNEW(kontetor,ITG,6*netet_);

  /* catalogueing the tetrahedral elements */

  NNEW(ifatet,ITG,4*netet_);
  NNEW(ifac,ITG,16*netet_);
  NNEW(itetfa,ITG,8*netet_);
  NNEW(planfa,double,16*netet_);
  NNEW(ipofa,ITG,nktet_);
  NNEW(bc,double,4*netet_);
  NNEW(cg,double,3*netet_);
  NNEW(ipoeln,ITG,nktet_);
  NNEW(ieln,ITG,8*netet_);
  NNEW(iparentel,ITG,netet_);

  FORTRAN(cattet,(kontet,&netet_,ifac,ne,ipkon,kon,ifatet,&ifreetet,
		  bc,itetfa,&ifreefa,planfa,ipofa,cotet,cg,ipoeln,
		  ieln,&ifreeln,lakon,kontetor,&iquad,istartset,iendset,
		  ialset,set,nset,filab,jfix,iparentel,jobnamec));
  
  /* catalogueing the edges of the tetrahedral elements */

  NNEW(iedg,ITG,18*netet_);
  NNEW(ipoed,ITG,nktet_);
  NNEW(ipoeled,ITG,6*netet_);
  NNEW(ieled,ITG,12*netet_);
  NNEW(iedtet,ITG,6*netet_);

  FORTRAN(catedges_refine,(&netet_,iedg,kontet,ipoed,&ifreeed,iedtet,
			   ipoeled,ieled,&ifreele));

  /* determining the external faces, edges and nodes */

  NNEW(iexternedg,ITG,6*netet_);
  NNEW(iexternfa,ITG,4*netet_);
  NNEW(iexternnode,ITG,nktet_);

  FORTRAN(determineextern,(ifac,itetfa,iedg,ipoed,iexternedg,
			   iexternfa,iexternnode,&nktet_,ipofa));

  /* store the nodes of the external edges in the unrefined mesh
     in iedgext; for edge i is iexternedg(i) a pointer to iedgext */

  NNEW(iedgext,ITG,3*(ifreeed-1));

  FORTRAN(midexternaledges,(iexternedg,&nexternedg,iedgext,
			    &ifreeed,ieled,ipoeled,iedg,iedtet,
			    kontetor));

  RENEW(iedgext,ITG,3*nexternedg);

  /* store the nodes of the external faces in the unrefined mesh
     in ifacext; for face i is iexternfac(i) a pointer to ifacext */

  NNEW(iedgextfa,ITG,2*nexternedg);
  NNEW(ifacext,ITG,6*(ifreefa-1));
  NNEW(ifacexted,ITG,3*(ifreefa-1));
  NNEW(iexternfaor,ITG,ifreefa-1);
  memcpy(iexternfaor,iexternfa,sizeof(ITG)*(ifreefa-1));
  NNEW(ialsete,ITG,ifreefa-1);

  FORTRAN(midexternalfaces,(iexternfa,&nexternfa,ifacext,
			    &ifreefa,itetfa,ifac,kontet,
			    kontetor,ialsete,&nexternel,
			    iedgextfa,ifacexted,ipoed,iedg,
			    iexternedg));

  RENEW(ifacext,ITG,6*nexternfa);
  RENEW(ifacexted,ITG,3*nexternfa);
  RENEW(ialsete,ITG,nexternel);
  SFREE(kontetor);

  NNEW(isharp,ITG,nexternedg);
  FORTRAN(checksharp,(&nexternedg,iedgextfa,cotet,ifacext,isharp));

  //  ij=9;
  NNEW(quality,double,nktet_);
  /*  FORTRAN(writemeshinp,(kontet,&netet_,cotet,&nktet,&ij,ipoed,iedg,
      iexternedg,quality));*/

  nloop=5;
  for(ij=0;ij<nloop;ij++){

    /* calculating the target element size in the nodes of the
       original mesh (master mesh) */

    NNEW(d,double,6*netet_);

    if(ij==0){
      NNEW(nh,ITG,nktet_);
      FORTRAN(calculateh,(nk,v,veold,stn,een,emn,epn,enern,qfn,errn,h,
			  filab,mi,d,nh,&dmin,ipoed,iedg,cotet,jfix));
      SFREE(nh);

      /* storing the master mesh and the target size field for 
	 subsequent interpolation purposes */

      getlocalresults(&integerglob,&doubleglob,&nktet_,cotet,h,&netet_,
		      kontet,ifatet,planfa,kontetor);

    }else{
      FORTRAN(updategeodata,(&nktet,&netet_,h,d,&dmin,ipoed,iedg,cotet,
			     planfa,bc,cg,kontet,ifac,ipofa,doubleglob,
			     integerglob,ipoeln));
    }

    /* treating the surface edges */

    iext=1;

    /* deciding which edges have to be split and in how many
       parts */

    NNEW(n,ITG,6*netet_);
    NNEW(r,double,6*netet_);

    FORTRAN(edgedivide,(&nnewnodes,&nktet_,ipoed,iexternedg,iedg,
			d,h,n,r,&iext,jfix));

    /* determining new nodes in the edges to be divided */

    NNEW(conewnodes,double,3*nnewnodes);
    NNEW(ibasenewnodes,ITG,nnewnodes);
    NNEW(iedgnewnodes,ITG,nnewnodes);
    NNEW(n1newnodes,ITG,nnewnodes);
    NNEW(n2newnodes,ITG,nnewnodes);
    NNEW(hnewnodes,double,nnewnodes);
	
    FORTRAN(newnodes,(&nktet_,ipoed,n,iedg,h,d,r,
		      conewnodes,cotet,ibasenewnodes,ipoeled,
		      ieled,doubleglob,integerglob,&nnewnodes,
		      iedgnewnodes,hnewnodes,n1newnodes,n2newnodes));

    NNEW(conewnodesorig,double,3*nnewnodes);
    memcpy(conewnodesorig,conewnodes,sizeof(double)*(3*nnewnodes));
	
    SFREE(n);SFREE(r);SFREE(d);

    /* adding random contributions to the new nodes towords
       the center of gravity of the base element */

    if(1.e-4*dmin<1.e-5){
      dmin=1.e-4*dmin;
    }else{
      dmin=1.e-5;
    }
    sgenrand(seed);
	
    for(i=0;i<nnewnodes;i++){

      /* generating a random number between 0.5 and 1 */

      unsiint=genrand();
      randomnumber=(double)unsiint/(unsigned long)0xffffffff;

      //	    printf("randomnumber=%f\n",randomnumber);

      randomnumber=0.5*(1.+randomnumber);

      /* moving the node to be inserted towards the center of
	 gravity of its basis element */

      dcg=sqrt((conewnodes[3*i]-cg[3*ibasenewnodes[i]-3])*
	       (conewnodes[3*i]-cg[3*ibasenewnodes[i]-3])+
	       (conewnodes[3*i+1]-cg[3*ibasenewnodes[i]-2])*
	       (conewnodes[3*i+1]-cg[3*ibasenewnodes[i]-2])+
	       (conewnodes[3*i+2]-cg[3*ibasenewnodes[i]-1])*
	       (conewnodes[3*i+2]-cg[3*ibasenewnodes[i]-1]));
      conewnodes[3*i]=conewnodes[3*i]+dmin*randomnumber*
	(cg[3*ibasenewnodes[i]-3]-conewnodes[3*i])/dcg;
      conewnodes[3*i+1]=conewnodes[3*i+1]+dmin*randomnumber*
	(cg[3*ibasenewnodes[i]-2]-conewnodes[3*i+1])/dcg;
      conewnodes[3*i+2]=conewnodes[3*i+2]+dmin*randomnumber*
	(cg[3*ibasenewnodes[i]-1]-conewnodes[3*i+2])/dcg;
    }

    /* rearranging the node order in an aleatoric way */

    sgenrand(seed);
    NNEW(irand1,ITG,nnewnodes);
    NNEW(irand2,ITG,nnewnodes);
    for(i=0;i<nnewnodes;i++){
      unsiint=genrand();
      irand1[i]=unsiint-2147483000-648;
      irand2[i]=i;
    }
    kflag=2;
    FORTRAN(isortii,(irand1,irand2,&nnewnodes,&kflag));
    SFREE(irand1);

    /* inserting the nodes one by one */

    for(i=0;i<nnewnodes;i++){
      node=irand2[i]+1;

      /* check whether the edge consisting of nodes n1 and n2 (n1<n2)
         still exists */

      FORTRAN(checkexiedge,(&n1newnodes[node-1],&n2newnodes[node-1],
			    ipoed,iedg,&node));
      if(node==0) continue;
      if(ibasenewnodes[node-1]==0) continue;

      /* finding the element to which the node belongs */

      ibase=ibasenewnodes[node-1];
      iedge=iedgnewnodes[node-1];
      height=hnewnodes[node-1];
      h[nktet]=height;
      ibasenewnodes[node-1]=0;

      /* adding the node to cotet */

      nktet++;
      cotet[3*nktet-3]=conewnodes[3*node-3];
      cotet[3*nktet-2]=conewnodes[3*node-2];
      cotet[3*nktet-1]=conewnodes[3*node-1];

      cotetorig[0]=conewnodesorig[3*node-3];
      cotetorig[1]=conewnodesorig[3*node-2];
      cotetorig[2]=conewnodesorig[3*node-1];

      node=nktet;

      /* finding the cavity and remeshing the cavity */
  
      FORTRAN(cavityext_refine,(kontet,ifatet,&ifreetet,bc,ifac,itetfa,
				&ifreefa,planfa,ipofa,cotet,&ibase,&node,
				iexternfa,ipoed,iedg,&ifreeed,ipoeled,ieled,
				&ifreele,&nktet,&netet_,ibasenewnodes,
				conewnodes,ipoeln,ieln,&ifreeln,&nnewnodes,
				iexternedg,iedtet,cotetorig,&iedge,iexternnode,
				cg,&height,iparentel));
      //	    printf("after cavityext in refinemesh.c\n");
	      
      /* check whether the fields for the tetrahedra information
	 have to be enlarged */
	      
      if(ifreetet>netet_-1000){
	newsize=(ITG)(1.1*netet_);
	kontet[4*netet_-1]=netet_+1;
		
	RENEW(kontet,ITG,4*newsize);
	RENEW(ifatet,ITG,4*newsize);
	RENEW(bc,double,4*newsize);
	RENEW(cg,double,3*newsize);
	RENEW(ifac,ITG,16*newsize);
	RENEW(iparentel,ITG,newsize);
		
	RENEW(iexternfa,ITG,4*newsize);
	for(k=4*netet_;k<4*newsize;k++)iexternfa[k]=0;
		
	RENEW(itetfa,ITG,8*newsize);
	RENEW(planfa,double,16*newsize);
	RENEW(ieln,ITG,8*newsize);
	RENEW(iedg,ITG,18*newsize);
		
	RENEW(iexternedg,ITG,6*newsize);
	for(k=6*netet_;k<6*newsize;k++)iexternedg[k]=0;
		
	RENEW(ipoeled,ITG,6*newsize);
	for(k=6*netet_;k<6*newsize;k++)ipoeled[k]=0;
		
	RENEW(ieled,ITG,12*newsize);
	RENEW(iedtet,ITG,6*newsize);
		
	FORTRAN(reinit_refine,(kontet,ifac,ieln,&netet_,&newsize,
			       ifatet,itetfa,iedg,ieled));
	netet_=newsize;
      }
	
      /* check whether the fields for the node information
	 have to be enlarged */
	
      if(nktet>nktet_-100){
	newsize=(ITG)(1.1*nktet_);
	RENEW(cotet,double,3*newsize);
	for(i=3*nktet_;i<3*newsize;i++)cotet[i]=0.;
	RENEW(h,double,newsize);
	RENEW(jfix,ITG,newsize);
	for(i=nktet_;i<newsize;i++)jfix[i]=0;
		
	RENEW(ipoeln,ITG,newsize);
	for(i=nktet_;i<newsize;i++)ipoeln[i]=0;
		
	RENEW(ipofa,ITG,newsize);
	for(i=nktet_;i<newsize;i++)ipofa[i]=0;
		
	RENEW(ipoed,ITG,newsize);
	for(i=nktet_;i<newsize;i++)ipoed[i]=0;
		
	RENEW(iexternnode,ITG,newsize);
	for(i=nktet_;i<newsize;i++)iexternnode[i]=0;
		
	nktet_=newsize;
      }
    }

    SFREE(conewnodes);SFREE(ibasenewnodes);SFREE(iedgnewnodes);
    SFREE(irand2);SFREE(conewnodesorig);SFREE(n1newnodes);
    SFREE(hnewnodes);SFREE(n2newnodes);

    /* treating the subsurface edges */

    iext=0;

    NNEW(d,double,6*netet_);
    FORTRAN(calculated,(&nktet,d,&dmin,ipoed,iedg,cotet));

    /* deciding which edges have to be split and in how many
       parts */

    NNEW(n,ITG,6*netet_);
    NNEW(r,double,6*netet_);

    FORTRAN(edgedivide,(&nnewnodes,&nktet_,ipoed,iexternedg,iedg,
			d,h,n,r,&iext,jfix));

    /* determining new nodes in the edges to be divided */

    NNEW(conewnodes,double,3*nnewnodes);
    NNEW(ibasenewnodes,ITG,nnewnodes);
    NNEW(iedgnewnodes,ITG,nnewnodes);
    NNEW(n1newnodes,ITG,nnewnodes);
    NNEW(n2newnodes,ITG,nnewnodes);
    NNEW(hnewnodes,double,nnewnodes);
	
    FORTRAN(newnodes,(&nktet_,ipoed,n,iedg,h,d,r,
		      conewnodes,cotet,ibasenewnodes,ipoeled,
		      ieled,doubleglob,integerglob,&nnewnodes,
		      iedgnewnodes,hnewnodes,n1newnodes,n2newnodes));
	
    SFREE(n);SFREE(r);SFREE(d);

    /* adding random contributions to the new nodes towords
       the center of gravity of the base element */

    if(1.e-4*dmin<1.e-5){
      dmin=1.e-4*dmin;
    }else{
      dmin=1.e-5;
    }
    //	mtseed(&iseed);
    sgenrand(seed);
	
    for(i=0;i<nnewnodes;i++){

      /* generating a random number between 0.5 and 1 */

      unsiint=genrand();
      randomnumber=(double)unsiint/(unsigned long)0xffffffff;

      //	    printf("randomnumber=%f\n",randomnumber);

      randomnumber=0.5*(1.+randomnumber);

      /* moving the node to be inserted towards the center of
	 gravity of its basis element */

      dcg=sqrt((conewnodes[3*i]-cg[3*ibasenewnodes[i]-3])*
	       (conewnodes[3*i]-cg[3*ibasenewnodes[i]-3])+
	       (conewnodes[3*i+1]-cg[3*ibasenewnodes[i]-2])*
	       (conewnodes[3*i+1]-cg[3*ibasenewnodes[i]-2])+
	       (conewnodes[3*i+2]-cg[3*ibasenewnodes[i]-1])*
	       (conewnodes[3*i+2]-cg[3*ibasenewnodes[i]-1]));
      conewnodes[3*i]=conewnodes[3*i]+dmin*randomnumber*
	(cg[3*ibasenewnodes[i]-3]-conewnodes[3*i])/dcg;
      conewnodes[3*i+1]=conewnodes[3*i+1]+dmin*randomnumber*
	(cg[3*ibasenewnodes[i]-2]-conewnodes[3*i+1])/dcg;
      conewnodes[3*i+2]=conewnodes[3*i+2]+dmin*randomnumber*
	(cg[3*ibasenewnodes[i]-1]-conewnodes[3*i+2])/dcg;
    }

    /* rearranging the node order in an aleatoric way */

    sgenrand(seed);
    NNEW(irand1,ITG,nnewnodes);
    NNEW(irand2,ITG,nnewnodes);
    for(i=0;i<nnewnodes;i++){
      unsiint=genrand();
      irand1[i]=unsiint-2147483000-648;
      irand2[i]=i;
    }
    kflag=2;
    FORTRAN(isortii,(irand1,irand2,&nnewnodes,&kflag));
    SFREE(irand1);

    /* inserting the nodes one by one */

    for(i=0;i<nnewnodes;i++){
      node=irand2[i]+1;

      /* check whether the edge consisting of nodes n1 and n2 (n1<n2)
         still exists */

      FORTRAN(checkexiedge,(&n1newnodes[node-1],&n2newnodes[node-1],
			    ipoed,iedg,&node));
      if(node==0) continue;
      if(ibasenewnodes[node-1]==0) continue;

      /* finding the element to which the node belongs */

      ibase=ibasenewnodes[node-1];
      height=hnewnodes[node-1];
      h[nktet]=height;
      ibasenewnodes[node-1]=0;

      /* adding the node to cotet */

      nktet++;
      cotet[3*nktet-3]=conewnodes[3*node-3];
      cotet[3*nktet-2]=conewnodes[3*node-2];
      cotet[3*nktet-1]=conewnodes[3*node-1];

      node=nktet;

      /* finding the cavity and remeshing the cavity */
  
      FORTRAN(cavity_refine,(kontet,ifatet,&ifreetet,bc,ifac,itetfa,
			     &ifreefa,planfa,ipofa,cotet,&ibase,&node,iexternfa,
			     ipoed,iedg,&ifreeed,ipoeled,ieled,&ifreele,&nktet,
			     &netet_,ibasenewnodes,conewnodes,ipoeln,ieln,
			     &ifreeln,&nnewnodes,iexternedg,iedtet,cg,
			     &height,iparentel));
	      
      /* check whether the fields for the tetrahedra information
	 have to be enlarged */
	      
      if(ifreetet>netet_-1000){
	newsize=(ITG)(1.1*netet_);
	kontet[4*netet_-1]=netet_+1;
		
	RENEW(kontet,ITG,4*newsize);
	RENEW(ifatet,ITG,4*newsize);
	RENEW(bc,double,4*newsize);
	RENEW(cg,double,3*newsize);
	RENEW(ifac,ITG,16*newsize);
	RENEW(iparentel,ITG,newsize);
		
	RENEW(iexternfa,ITG,4*newsize);
	for(k=4*netet_;k<4*newsize;k++)iexternfa[k]=0;
		
	RENEW(itetfa,ITG,8*newsize);
	RENEW(planfa,double,16*newsize);
	RENEW(ieln,ITG,8*newsize);
	RENEW(iedg,ITG,18*newsize);
		
	RENEW(iexternedg,ITG,6*newsize);
	for(k=6*netet_;k<6*newsize;k++)iexternedg[k]=0;
		
	RENEW(ipoeled,ITG,6*newsize);
	for(k=6*netet_;k<6*newsize;k++)ipoeled[k]=0;
		
	RENEW(ieled,ITG,12*newsize);
	RENEW(iedtet,ITG,6*newsize);
		
	FORTRAN(reinit_refine,(kontet,ifac,ieln,&netet_,&newsize,
			       ifatet,itetfa,iedg,ieled));
	netet_=newsize;
      }
	
      /* check whether the fields for the node information
	 have to be enlarged */
	
      if(nktet>nktet_-100){
	newsize=(ITG)(1.1*nktet_);
	RENEW(cotet,double,3*newsize);
	for(i=3*nktet_;i<3*newsize;i++)cotet[i]=0.;
	RENEW(h,double,newsize);
	RENEW(jfix,ITG,newsize);
	for(i=nktet_;i<newsize;i++)jfix[i]=0;
		
	RENEW(ipoeln,ITG,newsize);
	for(i=nktet_;i<newsize;i++)ipoeln[i]=0;
		
	RENEW(ipofa,ITG,newsize);
	for(i=nktet_;i<newsize;i++)ipofa[i]=0;
		
	RENEW(ipoed,ITG,newsize);
	for(i=nktet_;i<newsize;i++)ipoed[i]=0;
		
	RENEW(iexternnode,ITG,newsize);
	for(i=nktet_;i<newsize;i++)iexternnode[i]=0;
		
	nktet_=newsize;
      }
    }

    SFREE(conewnodes);SFREE(ibasenewnodes);SFREE(iedgnewnodes);
    SFREE(irand2);SFREE(hnewnodes);SFREE(n1newnodes);SFREE(n2newnodes);

    /* catalogueing the nodes according to the sharpness of the
       attached edges; determining the neighboring nodes to be 
       used in the smoothing */

    NNEW(iponn,ITG,nktet_);
    NNEW(inn,ITG,12*netet_);
    NNEW(idimsh,ITG,nktet_);
    NNEW(quality,double,netet_);

    FORTRAN(catnodes,(&ifreenn,inn,iponn,iedg,ipoed,&nktet_,iexternnode,
		      idimsh,isharp,iexternedg));

    RENEW(inn,ITG,2*ifreenn);

    /* smoothing the mesh */
    
    FORTRAN(smoothing,(inn,iponn,&nktet,iexternnode,&netet_,kontet,
		       cotet,ipoeln,ieln,h,quality,jfix));

    SFREE(iponn);SFREE(inn);SFREE(quality);SFREE(idimsh);
    ifreenn=1;
	
  }  /* enddo ij loop */
  
  /* catalogueing the nodes according to the sharpness of the
     attached edges; determine the neighboring nodes to be 
     used in the smoothing */

  NNEW(quality,double,netet_);
  NNEW(iponn,ITG,nktet_);
  NNEW(inn,ITG,12*netet_);
  NNEW(idimsh,ITG,nktet_);

  FORTRAN(catnodes,(&ifreenn,inn,iponn,iedg,ipoed,&nktet_,iexternnode,
		    idimsh,isharp,iexternedg));

  RENEW(inn,ITG,2*ifreenn);

  /* projecting the vertex nodes */
  
  NNEW(itreated,ITG,nktet);
  NNEW(ilist,ITG,nexternfa);

  FORTRAN(projectvertexnodes,(ipoed,iexternedg,iedgext,cotet,&nktet,iedg,
			      iexternfa,ifacext,itreated,ilist,isharp,ipofa,
			      ifac,iedgextfa,ifacexted,co,idimsh));

  /* determining the quality of the elements */

  ielem=0;
  FORTRAN(meshquality,(&netet_,kontet,cotet,quality,&ielem));

  /* remove the slivers */

  FORTRAN(removesliver,(&netet_,kontet,iexternnode,iedtet,iexternedg,quality));

  /* smoothing the mesh */

  FORTRAN(smoothing,(inn,iponn,&nktet,iexternnode,&netet_,kontet,
		     cotet,ipoeln,ieln,h,quality,jfix));

  SFREE(iponn);SFREE(inn);

  /* determining the quality of the elements */


  // NEXT LINES ARE NOT NEEDED? QUALITY IS NOT USED?
  //  ielem=0;
  //  FORTRAN(meshquality,(&netet_,kontet,cotet,quality,&ielem));

  /* projecting the vertex nodes */

  FORTRAN(projectvertexnodes,(ipoed,iexternedg,iedgext,cotet,&nktet,iedg,
			      iexternfa,ifacext,itreated,ilist,isharp,ipofa,
			      ifac,iedgextfa,ifacexted,co,idimsh));

  /* determining the quality of the elements */

  ielem=0;
  FORTRAN(meshquality,(&netet_,kontet,cotet,quality,&ielem));

  /* remove the slivers */

  FORTRAN(removesliver,(&netet_,kontet,iexternnode,iedtet,iexternedg,quality));
  
  /* generating middle nodes if needed (maximum of (6 * number of
     elements) = (maximum number of edges)) 
     projecting the middle nodes onto the surface */

  if(iquad==1){
    RENEW(cotet,double,3*nktet_+18*netet_);
    for(i=3*nktet_;i<3*nktet+18*netet_;i++)cotet[i]=0.;
  }

  NNEW(iedgmid,ITG,6*netet_);

  /* next line was set to force linear elements
     should be removed in the final version */

  //   iquad=0;
  FORTRAN(projectnodes,(&nktet_,ipoed,iedgmid,iexternedg,iedgext,cotet,&nktet,
			iedg,&iquad,iexternfa,ifacext,itreated,ilist,isharp,
			ipofa,ifac,iedgextfa,ifacexted,jfix,co,idimsh));

  SFREE(ialsete);SFREE(ilist);SFREE(isharp);SFREE(idimsh);

  /* store the refined mesh in input format */
  
  //  nktet=nktet+(*nk);
  NNEW(number,ITG,nktet);
  
  FORTRAN(writerefinemesh,(kontet,&netet_,cotet,&nktet,jobnamec,
			   ipkon,kon,lakon,&iquad,iedtet,iedgmid,
			   ne,number,jfix,iparentel,nk));

  SFREE(number);

  /* store the refined mesh of the part of the mesh which was refined 
     to frd-file */
  
  writenewmesh(&nktet,&netet_,cotet,&iquad,kontet,iedgmid,iedtet,mi,
	       matname,ithermal,jobnamec,output,nmat);
  
  SFREE(itreated);
   
  SFREE(iedgmid);

  SFREE(iexternedg);SFREE(iexternfa);SFREE(iexternnode);SFREE(iedgext);
  SFREE(ifacext);SFREE(ifacexted);SFREE(iexternfaor);SFREE(iedgextfa);

  SFREE(integerglob);SFREE(doubleglob);

  SFREE(iedg);SFREE(ipoed);SFREE(ipoeled);SFREE(ieled);SFREE(iedtet);

  SFREE(cotet);SFREE(kontet);SFREE(ifatet);SFREE(ifac);SFREE(itetfa);
  SFREE(planfa);SFREE(ipofa);SFREE(bc);SFREE(cg);SFREE(ipoeln);
  SFREE(ieln);SFREE(quality);SFREE(jfix);SFREE(h);SFREE(iparentel);

  return;
  
}
