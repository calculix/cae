/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2018 Guido Dhondt                          */

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
                char *jobnamec){

     /* refinement of tetrahedral meshes as a function of a field
        variable */

    ITG nktet_,netet_,*kontet=NULL,*ifatet=NULL,*ifac=NULL,*itetfa=NULL,
	*ipofa=NULL,*ipoeln=NULL,*ieln=NULL,ifreetet=1,ifreefa=1,
	ifreeln=1,*iedg=NULL,*ipoed=NULL,*ipoeled=NULL,*ieled=NULL,
	*iedtet=NULL,ifreeed=1,ifreele=1,*integerglob=NULL,*iexternfa=NULL,
	*iexternedg=NULL,*iexternnode=NULL,*n=NULL,nnewnodes,newsize,
        *ibasenewnodes=NULL,i,j,k,*irand1=NULL,*irand2=NULL,iext,iedge,
	kflag,node,ibase,nktet,*nh=NULL,ij,nloop,*iedgnewnodes=NULL,
        *kontetor=NULL,iquad=0,*iedgext=NULL,nexternedg,*iedgmid=NULL,
        *ifacext=NULL,nexternfa,*iexternfaor=NULL,*itreated=NULL,
        *ialsete=NULL,nexternel,*iedgextfa=NULL,*ifacexted=NULL,
        *ilist=NULL,*isharp=NULL;

//    ITG iseed=184389;
//    unsigned ITG unsiint;

    unsigned ITG unsiint,seed=184389;

    double *cotet=NULL,*planfa=NULL,*bc=NULL,*cg=NULL, *h=NULL,dcg,
	*doubleglob=NULL,*r=NULL,*d=NULL,*conewnodes=NULL,dmin,
	randomnumber,*conewnodesorig=NULL,cotetorig[3],height,
        *hnewnodes=NULL;

    nktet=*nk;
    nktet_=*nk+1000;
    netet_=*ne+1000;
    
    NNEW(cotet,double,3*nktet_);
    memcpy(&cotet[0],&co[0],sizeof(double)*3**nk);
    NNEW(h,double,nktet_);
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

    FORTRAN(cattet,(kontet,&netet_,ifac,ne,ipkon,kon,ifatet,&ifreetet,
		    bc,itetfa,&ifreefa,planfa,ipofa,cotet,cg,ipoeln,
		    ieln,&ifreeln,lakon,kontetor,&iquad));

    /* catalogueing the edges of the tetrahedral elements */

    NNEW(iedg,ITG,18*netet_);
    NNEW(ipoed,ITG,nktet_);
    NNEW(ipoeled,ITG,6*netet_);
    NNEW(ieled,ITG,12*netet_);
    NNEW(iedtet,ITG,6*netet_);

    FORTRAN(catedges,(&netet_,iedg,kontet,ipoed,&ifreeed,iedtet,
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

//    SFREE(kontetor);

    nloop=3;
    for(ij=0;ij<nloop;ij++){

        /* calculating the target element size in the nodes of the
           original mesh (master mesh) */

	NNEW(d,double,6*netet_);

	if(ij==0){
	    NNEW(nh,ITG,nktet_);
	    FORTRAN(calculateh,(nk,v,veold,stn,een,emn,epn,enern,qfn,errn,h,
				filab,mi,d,nh,&dmin,ipoed,iedg,cotet));
	    SFREE(nh);

            /* storing the master mesh and the target size field for 
               subsequent interpolation purposes */

	    getlocalresults(&integerglob,&doubleglob,&nktet_,cotet,h,&netet_,
			    kontet,ifatet,planfa,kontetor);

	    SFREE(kontetor);
	}else{
	    FORTRAN(calculated,(&nktet,d,&dmin,ipoed,iedg,cotet));
	}

        /* treating the surface edges */

	iext=1;

        /* deciding which edges have to be split and in how many
           parts */

	NNEW(n,ITG,6*netet_);
	NNEW(r,double,6*netet_);

	FORTRAN(edgedivide,(&nnewnodes,&nktet_,ipoed,iexternedg,iedg,
			    d,h,n,r,&iext));

        /* determining new nodes in the edges to be divided */

	NNEW(conewnodes,double,3*nnewnodes);
	NNEW(ibasenewnodes,ITG,nnewnodes);
	NNEW(iedgnewnodes,ITG,nnewnodes);
	NNEW(hnewnodes,double,nnewnodes);
	
	FORTRAN(newnodes,(&nktet_,ipoed,n,iedg,h,d,r,
			  conewnodes,cotet,ibasenewnodes,ipoeled,
			  ieled,doubleglob,integerglob,&nnewnodes,
			  iedgnewnodes,hnewnodes,iexternedg));

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
//	mtseed(&iseed);
	sgenrand(seed);
	
	for(i=0;i<nnewnodes;i++){

	    /* generateing a random number between 0.5 and 1 */

//	    unsiint=randmt();
//	    randomnumber=((unsiint-2147483000-648)/2147483648./2.);

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

//	mtseed(&iseed);
	sgenrand(seed);
	NNEW(irand1,ITG,nnewnodes);
	NNEW(irand2,ITG,nnewnodes);
	for(i=0;i<nnewnodes;i++){
//	    unsiint=randmt();
	    unsiint=genrand();
	    irand1[i]=unsiint-2147483000-648;
	    irand2[i]=i;
	}
	kflag=2;
	FORTRAN(isortii,(irand1,irand2,&nnewnodes,&kflag));
	SFREE(irand1);

        /* inserting the nodes one by one */

//	nnewnodes=4;
	for(i=0;i<nnewnodes;i++){
	    node=irand2[i]+1;

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
  
	    FORTRAN(cavityext,(kontet,ifatet,&ifreetet,bc,ifac,itetfa,
			    &ifreefa,planfa,ipofa,cotet,&ibase,&node,iexternfa,
			    ipoed,iedg,&ifreeed,ipoeled,ieled,&ifreele,&nktet,
			    &netet_,ibasenewnodes,conewnodes,ipoeln,ieln,
			    &ifreeln,&nnewnodes,iexternedg,iedtet,
			    cotetorig,&iedge,iexternnode,cg,&height));
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
		
		FORTRAN(reinit2,(kontet,ifac,ieln,&netet_,&newsize,
				 ifatet,itetfa,iedg,ieled));
		netet_=newsize;
	    }
	
	    /* check whether the fields for the node information
	       have to be enlarged */
	
	    if(nktet>nktet_-100){
		newsize=(ITG)(1.1*nktet_);
		RENEW(cotet,double,3*newsize);
		RENEW(h,double,newsize);
		
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
	SFREE(irand2);SFREE(conewnodesorig);
	SFREE(hnewnodes);

        /* treating the subsurface edges */

	iext=0;

	NNEW(d,double,6*netet_);
	FORTRAN(calculated,(&nktet,d,&dmin,ipoed,iedg,cotet));

        /* deciding which edges have to be split and in how many
           parts */

	NNEW(n,ITG,6*netet_);
	NNEW(r,double,6*netet_);

	FORTRAN(edgedivide,(&nnewnodes,&nktet_,ipoed,iexternedg,iedg,
			    d,h,n,r,&iext));

        /* determining new nodes in the edges to be divided */

	NNEW(conewnodes,double,3*nnewnodes);
	NNEW(ibasenewnodes,ITG,nnewnodes);
	NNEW(iedgnewnodes,ITG,nnewnodes);
	NNEW(hnewnodes,double,nnewnodes);
	
	FORTRAN(newnodes,(&nktet_,ipoed,n,iedg,h,d,r,
			  conewnodes,cotet,ibasenewnodes,ipoeled,
			  ieled,doubleglob,integerglob,&nnewnodes,
			  iedgnewnodes,hnewnodes,iexternedg));
	
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

	    /* generateing a random number between 0.5 and 1 */

//	    unsiint=randmt();
//	    randomnumber=((unsiint-2147483000-648)/2147483648./2.);
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

//	mtseed(&iseed);
	sgenrand(seed);
	NNEW(irand1,ITG,nnewnodes);
	NNEW(irand2,ITG,nnewnodes);
	for(i=0;i<nnewnodes;i++){
//	    unsiint=randmt();
	    unsiint=genrand();
	    irand1[i]=unsiint-2147483000-648;
	    irand2[i]=i;
	}
	kflag=2;
	FORTRAN(isortii,(irand1,irand2,&nnewnodes,&kflag));
	SFREE(irand1);

        /* inserting the nodes one by one */

/*	if(ij==0){
	    nnewnodes=1;
	}else{
	    nnewnodes=0;
	    }*/
//		nnewnodes=0;
	for(i=0;i<nnewnodes;i++){
	    node=irand2[i]+1;

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
  
	    FORTRAN(cavity,(kontet,ifatet,&ifreetet,bc,ifac,itetfa,
			    &ifreefa,planfa,ipofa,cotet,&ibase,&node,iexternfa,
			    ipoed,iedg,&ifreeed,ipoeled,ieled,&ifreele,&nktet,
			    &netet_,ibasenewnodes,conewnodes,ipoeln,ieln,
			    &ifreeln,&nnewnodes,iexternedg,iedtet,cg,
                            &height));
//	    printf("after cavity in refinemesh.c\n");
	      
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
		
		FORTRAN(reinit2,(kontet,ifac,ieln,&netet_,&newsize,
				 ifatet,itetfa,iedg,ieled));
		netet_=newsize;
	    }
	
	    /* check whether the fields for the node information
	       have to be enlarged */
	
	    if(nktet>nktet_-100){
		newsize=(ITG)(1.1*nktet_);
		RENEW(cotet,double,3*newsize);
		RENEW(h,double,newsize);
		
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
	SFREE(irand2);SFREE(hnewnodes);

	FORTRAN(writemeshinp,(kontet,&netet_,cotet,&nktet,&ij));
	
    }  /* enddo ij loop */

    /* projecting the newly generated nodes onto the external surface;
       generating middle nodes if needed (maximum of (6 * number of
       elements) = (maximum number of edges)) */

    RENEW(cotet,double,3*nktet_+18*netet_);
    SFREE(h);
    NNEW(iedgmid,ITG,6*netet_);
    NNEW(itreated,ITG,nktet);
    NNEW(ilist,ITG,nexternfa);
    NNEW(isharp,ITG,nexternedg);

    /* next line was set to force linear elements
       should be removed in the final version */

//    iquad=0;
    FORTRAN(projectnodes,(&nktet_,ipoed,iedgmid,iexternedg,
			  iedgext,cotet,&nktet,iedg,&iquad,
			  integerglob,doubleglob,iexternfa,
			  ifacext,itreated,ialsete,&nexternel,
                          ilist,isharp,ipofa,ifac,iedgextfa,
                          ifacexted,&nexternedg));

    SFREE(ialsete);SFREE(ilist);SFREE(isharp);

    FORTRAN(writerefinemesh,(kontet,&netet_,cotet,&nktet,jobnamec,
			     ipkon,kon,lakon,&iquad,iedtet,iedgmid,
		             itreated,ne));

    SFREE(itreated);
   
    SFREE(iedgmid);

    SFREE(iexternedg);SFREE(iexternfa);SFREE(iexternnode);SFREE(iedgext);
    SFREE(ifacext);SFREE(ifacexted);SFREE(iexternfaor);

    SFREE(integerglob);SFREE(doubleglob);

    SFREE(iedg);SFREE(ipoed);SFREE(ipoeled);SFREE(ieled);SFREE(iedtet);

    SFREE(cotet);SFREE(kontet);SFREE(ifatet);SFREE(ifac);SFREE(itetfa);
    SFREE(planfa);SFREE(ipofa);SFREE(bc);SFREE(cg);SFREE(ipoeln);
    SFREE(ieln);

  return;
  
}
