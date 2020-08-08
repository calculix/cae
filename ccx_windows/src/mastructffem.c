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

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

void mastructffem(ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
		  ITG *nodeboun,ITG *ndirboun,ITG *nboun,ITG *ipompc,
		  ITG *nodempc,ITG *nmpc,ITG *nactdoh,ITG *icolt,
		  ITG *icolv,ITG *icolp,ITG *icolk,ITG *jqt,ITG *jqv,ITG *jqp,
		  ITG *jqk,ITG **mast1p,ITG **irowtp,ITG **irowvp,ITG **irowpp,
		  ITG **irowkp,ITG *isolver,ITG *neqt,ITG *neqv,ITG *neqp,
		  ITG *neqk,ITG *ikmpc,ITG *ilmpc,ITG *ipointer,
		  ITG *nzst,ITG *nzsv,ITG *nzsp,ITG *nzsk,
		  ITG *ithermal,ITG *ikboun,ITG *ilboun,ITG *turbulent,
		  ITG *nactdok,ITG *ifreestream,ITG *nfreestream,
		  ITG *isolidface,ITG *nsolidface,ITG *nzs,ITG *iexplicit,
		  ITG *ielmat,ITG *inomat,char *labmpc){

  ITG i,j,k,l,jj,ll,id,index,jdof1,jdof2,idof1,idof2,mpc1,mpc2,id1,id2,
    ist1,ist2,node1,node2,isubtract,nmast,ifree,istart,istartold,idir,
    index1,index2,m,node,nzs_,ist,kflag,indexe,nope,isize,*mast1=NULL,
    *irowt=NULL,*irowv=NULL,*irowp=NULL,*irowk=NULL,fluid,imaterial;

  /* the indices in the comments follow FORTRAN convention, i.e. the
     fields start with 1 */

  mast1=*mast1p;
  irowt=*irowtp;irowv=*irowvp;irowp=*irowpp;irowk=*irowkp;

  kflag=2;
  
  /* initialisation of nactdoh */
  
  for(i=0;i<5**nk;++i){nactdoh[i]=0;}
  
  /* determining the active degrees of freedom due to elements */
  
  for(i=0;i<*ne;++i){
      
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[8*i],"F")!=0) continue;
    indexe=ipkon[i];
    if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
    else if (strcmp1(&lakon[8*i+3],"4")==0)nope=4;
    else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
    else continue;
      
    for(j=0;j<nope;++j){
      node=kon[indexe+j]-1;
      for(k=0;k<5;++k){
	nactdoh[5*node+k]=1;
	//	      inomat[node]=ielmat[i];
      }
    }
  }
  
  /* determining the active degrees of freedom due to mpc's */
  
  for(i=0;i<*nmpc;++i){

    /* check whether fluid MPC: a MPC is a fluid MPC if it
       contains at least one node belonging to a fluid
       element */

    index=ipompc[i]-1;
    fluid=0;
    do{
      node=nodempc[3*index]-1;

      /* check whether first velocity DOF is active */

      if(nactdoh[5*node+1]!=0){
	fluid=1;
	imaterial=inomat[node];
	break;
      }
      index=nodempc[3*index+2];
      if(index==0) break;
      index--;
    }while(1);

    if(fluid){
      index=ipompc[i]-1;
      do{
	node=nodempc[3*index]-1;
	idir=nodempc[3*index+1];
	inomat[node]=imaterial;
	/* next line may be necessary for pressure and semi-implicit scheme */
	/*	      nactdoh[5*node+idir]=1;*/
	index=nodempc[3*index+2];
	if(index==0) break;
	index--;
      }while(1);
    }
  }
  
  /* subtracting pressure SPCs and MPCs (only for the semi-implicit scheme) */

  if(*iexplicit==0){
    for(i=0;i<*nboun;++i){
      if(ndirboun[i]!=4){continue;}
      nactdoh[5*(nodeboun[i]-1)+ndirboun[i]]=0;
    }

    for(i=0;i<*nmpc;++i){
      index=ipompc[i]-1;
      if(nodempc[3*index+1]>4){continue;}
      if(nodempc[3*index+1]==4){
	nactdoh[5*(nodempc[3*index]-1)+4]=0;
      }
    }
  }      

  /* subtracting velocity=0 SPCs: setting the dof to a negative
     number, the absolute value of which is twice the SPC number */

  /*  for(i=0;i<*nboun;++i){
    if((ndirboun[i]<1)||(ndirboun[i]>3)) continue;
    nactdoh[5*(nodeboun[i]-1)+ndirboun[i]]=-2*(i+1);
    }*/
  
  /* numbering the active degrees of freedom */
  
  *neqt=0;*neqv=0;*neqp=0;
  for(i=0;i<*nk;++i){
    if(*ithermal>1){
      if(nactdoh[5*i]>0){
	++(*neqt);
	nactdoh[5*i]=*neqt;
      }
    }
    for(j=1;j<4;++j){
      if(nactdoh[5*i+j]>0){
	++(*neqv);
	nactdoh[5*i+j]=*neqv;
      }
    }
    if(nactdoh[5*i+4]>0){
      ++(*neqp);
      nactdoh[5*i+4]=*neqp;
    }
  }
  if(*ithermal>1) printf("neqttt=%d\n",*neqt);
  printf("neqvvv=%d\n",*neqv);
  printf("neqppp=%d\n",*neqp);
  
  /* determining the turbulence degrees of freedom */

  if(*turbulent!=0){

    /* initialisation of nactdok */
      
    for(i=0;i<*nk;++i){nactdok[i]=0;}
      
    /* determining the turbulence degrees of freedom due to elements */
      
    for(i=0;i<*ne;++i){
	  
      if(ipkon[i]<0) continue;
      if(strcmp1(&lakon[8*i],"F")!=0) continue;
      indexe=ipkon[i];
      if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
      else if (strcmp1(&lakon[8*i+3],"4")==0)nope=4;
      else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
      else continue;
	  
      for(j=0;j<nope;++j){
	node=kon[indexe+j]-1;
	nactdok[node]=1;
      }
    }
  
    /* determining the active degrees of freedom due to mpc's */

    /* if the MPC is a cyclic pressure MPC it is
       also applied to the turbulence parameters */
  
    for(i=0;i<*nmpc;++i){
      if(strcmp1(&labmpc[20*i],"CYCLIC")!=0) continue;
      index=ipompc[i]-1;
      do{
	node=nodempc[3*index]-1;
	if(inomat[node]==0) break;
	idir=nodempc[3*index+1];
	if(idir!=4) break;
	/*	      nactdok[node]=1;*/
	index=nodempc[3*index+2];
	if(index==0) break;
	index--;
      }while(1);
    }
      
    /* numbering the active degrees of freedom */
      
    *neqk=0;
    for(i=0;i<*nk;++i){
      if(nactdok[i]!=0){
	++(*neqk);
	nactdok[i]=*neqk;
      }
    }
  }    
  
  /* temperature entries*/
  
  if(*ithermal>1){

    ifree=0;
    nzs_=*nzs;
      
    /* determining the position of each nonzero matrix element
	 
       mast1(ipointer(i)) = first nonzero row in column i
       irow(ipointer(i))  points to further nonzero elements in 
       column i */
      
    for(i=0;i<3**nk;++i){ipointer[i]=0;}
      
    for(i=0;i<*ne;++i){
	  
      if(ipkon[i]<0) continue;
      if(strcmp1(&lakon[8*i],"F")!=0) continue;
      indexe=ipkon[i];
      if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
      else if (strcmp1(&lakon[8*i+3],"4")==0)nope=4;
      else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
      else continue;
	  
      for(jj=0;jj<nope;++jj){
	      
	j=jj;
	      
	node1=kon[indexe+j];
	jdof1=nactdoh[5*(node1-1)];
	      
	for(ll=jj;ll<nope;++ll){
		  
	  l=ll;
		  
	  node2=kon[indexe+l];
	  jdof2=nactdoh[5*(node2-1)];
		  
	  /* check whether one of the DOF belongs to a SPC or MPC */
		  
	  if((jdof1>0)&&(jdof2>0)){
	    insertcbs(ipointer,&mast1,&irowt,&jdof1,&jdof2,&ifree,&nzs_);
	  }
	  else if((jdof1>0)||(jdof2>0)){
		      
	    /* idof1: genuine DOF
	       idof2: nominal DOF of the SPC/MPC */
		      
	    if(jdof1<=0){
	      idof1=jdof2;
	      idof2=jdof1;}
	    else{
	      idof1=jdof1;
	      idof2=jdof2;}
		      
	    if(*nmpc>0){

	      if(idof2!=2*(idof2/2)){
			      
		/* regular DOF / MPC */
			      
		id=(-idof2+1)/2;
		ist=ipompc[id-1];
		index=nodempc[3*ist-1];
		if(index==0) continue;
		while(1){
		  idof2=nactdoh[5*(nodempc[3*index-3]-1)+nodempc[3*index-2]];
		  if(idof2>0){
		    insertcbs(ipointer,&mast1,&irowt,&idof1,&idof2,&ifree,&nzs_);
		  }
		  index=nodempc[3*index-1];
		  if(index==0) break;
		}
		continue;
	      }
	    }
		      
	  }
		  
	  else{
	    idof1=jdof1;
	    idof2=jdof2;
	    mpc1=0;
	    mpc2=0;
	    if(*nmpc>0){
	      if(idof1!=2*(idof1/2)) mpc1=1;
	      if(idof2!=2*(idof2/2)) mpc2=1;
	    }
	    if((mpc1==1)&&(mpc2==1)){
	      id1=(-idof1+1)/2;
	      id2=(-idof2+1)/2;
	      if(id1==id2){
			      
		/* MPC id1 / MPC id1 */
			      
		ist=ipompc[id1-1];
		index1=nodempc[3*ist-1];
		if(index1==0) continue;
		while(1){
		  idof1=nactdoh[5*(nodempc[3*index1-3]-1)+nodempc[3*index1-2]];
		  index2=index1;
		  while(1){
		    idof2=nactdoh[5*(nodempc[3*index2-3]-1)+nodempc[3*index2-2]];
		    if((idof1>0)&&(idof2>0)){
		      insertcbs(ipointer,&mast1,&irowt,&idof1,&idof2,&ifree,&nzs_);}
		    index2=nodempc[3*index2-1];
		    if(index2==0) break;
		  }
		  index1=nodempc[3*index1-1];
		  if(index1==0) break;
		}
	      }
			  
	      else{
			      
		/* MPC id1 /MPC id2 */
			      
		ist1=ipompc[id1-1];
		index1=nodempc[3*ist1-1];
		if(index1==0) continue;
		while(1){
		  idof1=nactdoh[5*(nodempc[3*index1-3]-1)+nodempc[3*index1-2]];
		  ist2=ipompc[id2-1];
		  index2=nodempc[3*ist2-1];
		  if(index2==0){
		    index1=nodempc[3*index1-1];
		    if(index1==0){break;}
		    else{continue;}
		  }
		  while(1){
		    idof2=nactdoh[5*(nodempc[3*index2-3]-1)+nodempc[3*index2-2]];
		    if((idof1>0)&&(idof2>0)){
		      insertcbs(ipointer,&mast1,&irowt,&idof1,&idof2,&ifree,&nzs_);}
		    index2=nodempc[3*index2-1];
		    if(index2==0) break;
		  }
		  index1=nodempc[3*index1-1];
		  if(index1==0) break;
		}
	      }
	    }
	  }
	}
      }
    }
      
    for(i=0;i<*neqt;++i){
      if(ipointer[i]==0){
	if(i>=*neqt) continue;
	printf("*ERROR in mastructf: zero column\n");
	FORTRAN(stop,());
      }
      istart=ipointer[i];
      while(1){
	istartold=istart;
	istart=irowt[istart-1];
	irowt[istartold-1]=i+1;
	if(istart==0) break;
      }
    }
      
    /* defining icolt and jqt */
      
    if(*neqt==0){
      printf("\n*WARNING in mastructf: no degrees of freedom in the temperature matrix\n\n");
    }
      
    nmast=ifree;
      
    /* summary */
      
    printf(" number of temperature equations\n");
    printf(" %d\n",*neqt);
    printf(" number of nonzero temperature matrix elements\n");
    printf(" %d\n",nmast);
    printf("\n");
      
    /* changing the meaning of icolt,jqt,mast1,irowt:
	 
       - irowt is going to contain the row numbers of the SUBdiagonal
       nonzero's, column per column
       - mast1 contains the column numbers
       - icolt(i)=# SUBdiagonal nonzero's in column i
       - jqt(i)= location in field irow of the first SUBdiagonal
       nonzero in column i
	 
    */
      
    /* switching from a SUPERdiagonal inventory to a SUBdiagonal one */
      
    FORTRAN(isortii,(mast1,irowt,&nmast,&kflag));
      
    /* filtering out the diagonal elements and generating icolt and jqt */
      
    isubtract=0;
    for(i=0;i<*neqt;++i){icolt[i]=0;}
    k=0;
    for(i=0;i<nmast;++i){
      if(mast1[i]==irowt[i]){++isubtract;}
      else{
	mast1[i-isubtract]=mast1[i];
	irowt[i-isubtract]=irowt[i];
	if(k!=mast1[i]){
	  for(l=k;l<mast1[i];++l){jqt[l]=i+1-isubtract;}
	  k=mast1[i];
	}
	++icolt[k-1];
      }
    }
    nmast=nmast-isubtract;
    for(l=k;l<*neqt+1;++l){jqt[l]=nmast+1;}
      
    for(i=0;i<*neqt;++i){
      if(jqt[i+1]-jqt[i]>0){
	isize=jqt[i+1]-jqt[i];
	FORTRAN(isortii,(&irowt[jqt[i]-1],&mast1[jqt[i]-1],&isize,&kflag));
      }
    }
      
    if(*neqt==0){*nzst=0;}
    else{*nzst=jqt[*neqt]-1;}
      
  }
  
  /* velocity entries */
  
  ifree=0;
  nzs_=*nzs;
  RENEW(mast1,ITG,nzs_);
  for(i=0;i<nzs_;i++){mast1[i]=0;}
  
  /* determining the position of each nonzero matrix element
     
     mast1(ipointer(i)) = first nonzero row in column i
     irowv(ipointer(i))  points to further nonzero elements in 
     column i */
  
  for(i=0;i<3**nk;++i){ipointer[i]=0;}
  
  for(i=0;i<*ne;++i){
      
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[8*i],"F")!=0) continue;
    indexe=ipkon[i];
    if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
    else if (strcmp1(&lakon[8*i+3],"4")==0)nope=4;
    else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
    else continue;
      
    for(jj=0;jj<3*nope;++jj){
	  
      j=jj/3;
      k=jj-3*j;
	  
      node1=kon[indexe+j];
      jdof1=nactdoh[5*(node1-1)+k+1];
	  
      for(ll=jj;ll<3*nope;++ll){
	      
	l=ll/3;
	m=ll-3*l;
	      
	node2=kon[indexe+l];
	jdof2=nactdoh[5*(node2-1)+m+1];
	      
	/* check whether one of the DOF belongs to a SPC or MPC */
	      
	if((jdof1>0)&&(jdof2>0)){
	  insertcbs(ipointer,&mast1,&irowv,&jdof1,&jdof2,&ifree,&nzs_);
	}
	else if((jdof1>0)||(jdof2>0)){
		  
	  /* idof1: genuine DOF
	     idof2: nominal DOF of the SPC/MPC */
		  
	  if(jdof1<=0){
	    idof1=jdof2;
	    idof2=jdof1;}
	  else{
	    idof1=jdof1;
	    idof2=jdof2;}
		  
	  if(*nmpc>0){

	    if(idof2!=2*(idof2/2)){
			  
	      /* regular DOF / MPC */
			  
	      id=(-idof2+1)/2;
	      ist=ipompc[id-1];
	      index=nodempc[3*ist-1];
	      if(index==0) continue;
	      while(1){
		idof2=nactdoh[5*(nodempc[3*index-3]-1)+nodempc[3*index-2]];
		if(idof2>0){
		  insertcbs(ipointer,&mast1,&irowv,&idof1,&idof2,&ifree,&nzs_);
		}
		index=nodempc[3*index-1];
		if(index==0) break;
	      }
	      continue;
	    }
	  }
	}
	      
	else{
	  idof1=jdof1;
	  idof2=jdof2;
	  mpc1=0;
	  mpc2=0;
	  if(*nmpc>0){
	    if(idof1!=2*(idof1/2)) mpc1=1;
	    if(idof2!=2*(idof2/2)) mpc2=1;
	  }
	  if((mpc1==1)&&(mpc2==1)){
	    id1=(-idof1+1)/2;
	    id2=(-idof2+1)/2;
	    if(id1==id2){
			  
	      /* MPC id1 / MPC id1 */
			  
	      ist=ipompc[id1-1];
	      index1=nodempc[3*ist-1];
	      if(index1==0) continue;
	      while(1){
		idof1=nactdoh[5*(nodempc[3*index1-3]-1)+nodempc[3*index1-2]];
		index2=index1;
		while(1){
		  idof2=nactdoh[5*(nodempc[3*index2-3]-1)+nodempc[3*index2-2]];
		  if((idof1>0)&&(idof2>0)){
		    insertcbs(ipointer,&mast1,&irowv,&idof1,&idof2,&ifree,&nzs_);}
		  index2=nodempc[3*index2-1];
		  if(index2==0) break;
		}
		index1=nodempc[3*index1-1];
		if(index1==0) break;
	      }
	    }
		      
	    else{
			  
	      /* MPC id1 /MPC id2 */
			  
	      ist1=ipompc[id1-1];
	      index1=nodempc[3*ist1-1];
	      if(index1==0) continue;
	      while(1){
		idof1=nactdoh[5*(nodempc[3*index1-3]-1)+nodempc[3*index1-2]];
		ist2=ipompc[id2-1];
		index2=nodempc[3*ist2-1];
		if(index2==0){
		  index1=nodempc[3*index1-1];
		  if(index1==0){break;}
		  else{continue;}
		}
		while(1){
		  idof2=nactdoh[5*(nodempc[3*index2-3]-1)+nodempc[3*index2-2]];
		  if((idof1>0)&&(idof2>0)){
		    insertcbs(ipointer,&mast1,&irowv,&idof1,&idof2,&ifree,&nzs_);}
		  index2=nodempc[3*index2-1];
		  if(index2==0) break;
		}
		index1=nodempc[3*index1-1];
		if(index1==0) break;
	      }
	    }
	  }
	}
      }
    }
  }
  
  for(i=0;i<*neqv;++i){
    if(ipointer[i]==0){
      if(i>=*neqv) continue;
      printf("*ERROR in mastructf: zero column in the velocity matrix\n");
      printf("       DOF %d\n",i);
      FORTRAN(stop,());
    }
    istart=ipointer[i];
    while(1){
      istartold=istart;
      istart=irowv[istart-1];
      irowv[istartold-1]=i+1;
      if(istart==0) break;
    }
  }
  
  /* defining icolv and jqv */
  
  if(*neqv==0){
    printf("\n*WARNING in mastructf: no degrees of freedom in the velocity matrix\n\n");
  }
  
  nmast=ifree;
  
  /* summary */
  
  printf(" number of velocity equations\n");
  printf(" %d\n",*neqv);
  printf(" number of nonzero velocity matrix elements\n");
  printf(" %d\n",nmast);
  printf("\n");
  
  /* changing the meaning of icolv,jqv,mast1,irowv:
     
     - irowv is going to contain the row numbers of the SUBdiagonal
     nonzero's, column per column
     - mast1 contains the column numbers
     - icolv(i)=# SUBdiagonal nonzero's in column i
     - jqv(i)= location in field irow of the first SUBdiagonal
     nonzero in column i
  
  */
  
  /* switching from a SUPERdiagonal inventory to a SUBdiagonal one */
  
  FORTRAN(isortii,(mast1,irowv,&nmast,&kflag));
  
  /* filtering out the diagonal elements and generating icolv and jqv */
  
  isubtract=0;
  for(i=0;i<*neqv;++i){icolv[i]=0;}
  k=0;
  for(i=0;i<nmast;++i){
    if(mast1[i]==irowv[i]){++isubtract;}
    else{
      mast1[i-isubtract]=mast1[i];
      irowv[i-isubtract]=irowv[i];
      if(k!=mast1[i]){
	for(l=k;l<mast1[i];++l){jqv[l]=i+1-isubtract;}
	k=mast1[i];
      }
      ++icolv[k-1];
    }
  }
  nmast=nmast-isubtract;
  for(l=k;l<*neqv+1;++l){jqv[l]=nmast+1;}
  
  for(i=0;i<*neqv;++i){
    if(jqv[i+1]-jqv[i]>0){
      isize=jqv[i+1]-jqv[i];
      FORTRAN(isortii,(&irowv[jqv[i]-1],&mast1[jqv[i]-1],&isize,&kflag));
    }
  }
  
  if(*neqv==0){*nzsv=0;}
  else{*nzsv=jqv[*neqv]-1;}
  
  /* pressure entries */
  
  ifree=0;
  nzs_=*nzs;
  RENEW(mast1,ITG,nzs_);
  for(i=0;i<nzs_;i++){mast1[i]=0;}
  
  /* determining the position of each nonzero matrix element
     
     mast1(ipointer(i)) = first nonzero row in column i
     irowp(ipointer(i))  points to further nonzero elements in 
     column i */
  
  for(i=0;i<3**nk;++i){ipointer[i]=0;}
  
  for(i=0;i<*ne;++i){
      
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[8*i],"F")!=0) continue;
    indexe=ipkon[i];
    if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
    else if (strcmp1(&lakon[8*i+3],"4")==0)nope=4;
    else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
    else continue;
      
    for(jj=0;jj<nope;++jj){
	  
      j=jj;
	  
      node1=kon[indexe+j];
      jdof1=nactdoh[5*(node1-1)+4];
	  
      for(ll=jj;ll<nope;++ll){
	      
	l=ll;
	      
	node2=kon[indexe+l];
	jdof2=nactdoh[5*(node2-1)+4];
	      
	/* check whether one of the DOF belongs to a SPC or MPC */
	      
	if((jdof1>0)&&(jdof2>0)){
	  insertcbs(ipointer,&mast1,&irowp,&jdof1,&jdof2,&ifree,&nzs_);
	}
	else if((jdof1>0)||(jdof2>0)){
		  
	  /* idof1: genuine DOF
	     idof2: nominal DOF of the SPC/MPC */
		  
	  if(jdof1==0){
	    idof1=jdof2;
	    idof2=jdof1;}
	  else{
	    idof1=jdof1;
	    idof2=jdof2;}
		  
	  if(*nmpc>0){


	    if(idof2!=2*(idof2/2)){
			  
	      /* regular DOF / MPC */
			  
	      id=(-idof2+1)/2;
	      ist=ipompc[id-1];
	      index=nodempc[3*ist-1];
	      if(index==0) continue;
	      while(1){
		idof2=nactdoh[5*(nodempc[3*index-3]-1)+4];
		if(idof2>0){
		  insertcbs(ipointer,&mast1,&irowp,&idof1,&idof2,&ifree,&nzs_);
		}
		index=nodempc[3*index-1];
		if(index==0) break;
	      }
	      continue;
	    }
	  }
	}
	      
	else{
	  idof1=jdof1;
	  idof2=jdof2;
	  mpc1=0;
	  mpc2=0;
	  if(*nmpc>0){
	    if(idof1!=2*(idof1/2)) mpc1=1;
	    if(idof2!=2*(idof2/2)) mpc2=1;
	  }
	  if((mpc1==1)&&(mpc2==1)){
	    id1=(-idof1+1)/2;
	    id2=(-idof2+1)/2;
	    if(id1==id2){
			  
	      /* MPC id1 / MPC id1 */
			  
	      ist=ipompc[id1-1];
	      index1=nodempc[3*ist-1];
	      if(index1==0) continue;
	      while(1){
		idof1=nactdoh[5*(nodempc[3*index1-3]-1)+4];
		index2=index1;
		while(1){
		  idof2=nactdoh[5*(nodempc[3*index2-3]-1)+4];
		  if((idof1>0)&&(idof2>0)){
		    insertcbs(ipointer,&mast1,&irowp,&idof1,&idof2,&ifree,&nzs_);}
		  index2=nodempc[3*index2-1];
		  if(index2==0) break;
		}
		index1=nodempc[3*index1-1];
		if(index1==0) break;
	      }
	    }
		      
	    else{
			  
	      /* MPC id1 /MPC id2 */
			  
	      ist1=ipompc[id1-1];
	      index1=nodempc[3*ist1-1];
	      if(index1==0) continue;
	      while(1){
		idof1=nactdoh[5*(nodempc[3*index1-3]-1)+4];
		ist2=ipompc[id2-1];
		index2=nodempc[3*ist2-1];
		if(index2==0){
		  index1=nodempc[3*index1-1];
		  if(index1==0){break;}
		  else{continue;}
		}
		while(1){
		  idof2=nactdoh[5*(nodempc[3*index2-3]-1)+4];
		  if((idof1>0)&&(idof2>0)){
		    insertcbs(ipointer,&mast1,&irowp,&idof1,&idof2,&ifree,&nzs_);}
		  index2=nodempc[3*index2-1];
		  if(index2==0) break;
		}
		index1=nodempc[3*index1-1];
		if(index1==0) break;
	      }
	    }
	  }
	}
      }
    }
  }
    
  for(i=0;i<*neqp;++i){
    if(ipointer[i]==0){
      if(i>=*neqp) continue;
      printf("*ERROR in mastructf: zero column\n");
      FORTRAN(stop,());
    }
    istart=ipointer[i];
    while(1){
      istartold=istart;
      istart=irowp[istart-1];
      irowp[istartold-1]=i+1;
      if(istart==0) break;
    }
  }
  
  /* defining icolp and jqp */
  
  if(*neqp==0){
    printf("\n*WARNING in matructf: no degrees of freedom in the pressure matrix\n\n");
  }
  
  nmast=ifree;
  
  /* summary */
  
  printf(" number of pressure equations\n");
  printf(" %d\n",*neqp);
  printf(" number of nonzero pressure matrix elements\n");
  printf(" %d\n",nmast);
  printf("\n");
  
  /* changing the meaning of icolp,jqp,mast1,irowp:
     
     - irowp is going to contain the row numbers of the SUBdiagonal
     nonzero's, column per column
     - mast1 contains the column numbers
     - icolp(i)=# SUBdiagonal nonzero's in column i
     - jqp(i)= location in field irow of the first SUBdiagonal
     nonzero in column i
  
  */
  
  /* switching from a SUPERdiagonal inventory to a SUBdiagonal one */
  
  FORTRAN(isortii,(mast1,irowp,&nmast,&kflag));
  
  /* filtering out the diagonal elements and generating icolp and jqp */
  
  isubtract=0;
  for(i=0;i<*neqp;++i){icolp[i]=0;}
  k=0;
  for(i=0;i<nmast;++i){
    if(mast1[i]==irowp[i]){++isubtract;}
    else{
      mast1[i-isubtract]=mast1[i];
      irowp[i-isubtract]=irowp[i];
      if(k!=mast1[i]){
	for(l=k;l<mast1[i];++l){jqp[l]=i+1-isubtract;}
	k=mast1[i];
      }
      ++icolp[k-1];
    }
  }
  nmast=nmast-isubtract;
  for(l=k;l<*neqp+1;++l){jqp[l]=nmast+1;}
  
  for(i=0;i<*neqp;++i){
    if(jqp[i+1]-jqp[i]>0){
      isize=jqp[i+1]-jqp[i];
      FORTRAN(isortii,(&irowp[jqp[i]-1],&mast1[jqp[i]-1],&isize,&kflag));
    }
  }
  
  if(*neqp==0){*nzsp=0;}
  else{*nzsp=jqp[*neqp]-1;}
  
  /* turbulence entries */
  
  if(*turbulent!=0){
      
    ifree=0;
    nzs_=*nzs;
    RENEW(mast1,ITG,nzs_);
    for(i=0;i<nzs_;i++){mast1[i]=0;}
      
    /* determining the position of each nonzero matrix element
	 
       mast1(ipointer(i)) = first nonzero row in column i
       irowk(ipointer(i))  points to further nonzero elements in 
       column i */
      
    for(i=0;i<3**nk;++i){ipointer[i]=0;}
      
    for(i=0;i<*ne;++i){
	  
      if(ipkon[i]<0) continue;
      if(strcmp1(&lakon[8*i],"F")!=0) continue;
      indexe=ipkon[i];
      if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
      else if (strcmp1(&lakon[8*i+3],"4")==0)nope=4;
      else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
      else continue;
	  
      for(jj=0;jj<nope;++jj){
	      
	j=jj;
	      
	node1=kon[indexe+j];
	jdof1=nactdok[node1-1];
	      
	for(ll=jj;ll<nope;++ll){
		  
	  l=ll;
		  
	  node2=kon[indexe+l];
	  jdof2=nactdok[node2-1];
		  
	  /* check whether one of the DOF belongs to a SPC or MPC */
		  
	  if((jdof1>0)&&(jdof2>0)){
	    insertcbs(ipointer,&mast1,&irowk,&jdof1,&jdof2,&ifree,&nzs_);
	  }
	  else if((jdof1>0)||(jdof2>0)){
		      
	    /* idof1: genuine DOF
	       idof2: nominal DOF of the SPC/MPC */
		      
	    if(jdof1<=0){
	      idof1=jdof2;
	      idof2=jdof1;}
	    else{
	      idof1=jdof1;
	      idof2=jdof2;}
		      
	    if(*nmpc>0){

	      if(idof2!=2*(idof2/2)){
			      
		/* regular DOF / MPC */
			      
		id=(-idof2+1)/2;
		ist=ipompc[id-1];
		index=nodempc[3*ist-1];
		if(index==0) continue;
		while(1){
		  idof2=nactdok[nodempc[3*index-3]-1];
		  if(idof2>0){
		    insertcbs(ipointer,&mast1,&irowk,&idof1,&idof2,&ifree,&nzs_);
		  }
		  index=nodempc[3*index-1];
		  if(index==0) break;
		}
		continue;
	      }
	    }
	  }
		  
	  else{
	    idof1=jdof1;
	    idof2=jdof2;
	    mpc1=0;
	    mpc2=0;
	    if(*nmpc>0){
	      if(idof1!=2*(idof1/2)) mpc1=1;
	      if(idof2!=2*(idof2/2)) mpc2=1;
	    }
	    if((mpc1==1)&&(mpc2==1)){
	      id1=(-idof1+1)/2;
	      id2=(-idof2+1)/2;
	      if(id1==id2){
			      
		/* MPC id1 / MPC id1 */
			      
		ist=ipompc[id1-1];
		index1=nodempc[3*ist-1];
		if(index1==0) continue;
		while(1){
		  idof1=nactdok[nodempc[3*index1-3]-1];
		  index2=index1;
		  while(1){
		    idof2=nactdok[nodempc[3*index2-3]-1];
		    if((idof1>0)&&(idof2>0)){
		      insertcbs(ipointer,&mast1,&irowk,&idof1,&idof2,&ifree,&nzs_);}
		    index2=nodempc[3*index2-1];
		    if(index2==0) break;
		  }
		  index1=nodempc[3*index1-1];
		  if(index1==0) break;
		}
	      }
			  
	      else{
			      
		/* MPC id1 /MPC id2 */
			      
		ist1=ipompc[id1-1];
		index1=nodempc[3*ist1-1];
		if(index1==0) continue;
		while(1){
		  idof1=nactdok[nodempc[3*index1-3]-1];
		  ist2=ipompc[id2-1];
		  index2=nodempc[3*ist2-1];
		  if(index2==0){
		    index1=nodempc[3*index1-1];
		    if(index1==0){break;}
		    else{continue;}
		  }
		  while(1){
		    idof2=nactdok[nodempc[3*index2-3]-1];
		    if((idof1>0)&&(idof2>0)){
		      insertcbs(ipointer,&mast1,&irowk,&idof1,&idof2,&ifree,&nzs_);}
		    index2=nodempc[3*index2-1];
		    if(index2==0) break;
		  }
		  index1=nodempc[3*index1-1];
		  if(index1==0) break;
		}
	      }
	    }
	  }
	}
      }
    }
      
    for(i=0;i<*neqk;++i){
      if(ipointer[i]==0){
	if(i>=*neqk) continue;
	printf("*ERROR in mastructf: zero column\n");
	FORTRAN(stop,());
      }
      istart=ipointer[i];
      while(1){
	istartold=istart;
	istart=irowk[istart-1];
	irowk[istartold-1]=i+1;
	if(istart==0) break;
      }
    }
      
    /* defining icolk and jqk */
      
    if(*neqk==0){
      printf("\n*WARNING in matructf: no degrees of freedom in the turbulence matrix\n\n");
    }
      
    nmast=ifree;
      
    /* summary */
      
    printf(" number of turbulence equations\n");
    printf(" %d\n",*neqk);
    printf(" number of nonzero turbulence matrix elements\n");
    printf(" %d\n",nmast);
    printf("\n");
      
    /* changing the meaning of icolk,jqk,mast1,irowk:
	 
       - irowk is going to contain the row numbers of the SUBdiagonal
       nonzero's, column per column
       - mast1 contains the column numbers
       - icolk(i)=# SUBdiagonal nonzero's in column i
       - jqk(i)= location in field irow of the first SUBdiagonal
       nonzero in column i
      
    */
      
    /* switching from a SUPERdiagonal inventory to a SUBdiagonal one */
      
    FORTRAN(isortii,(mast1,irowk,&nmast,&kflag));
      
    /* filtering out the diagonal elements and generating icolk and jqk */
      
    isubtract=0;
    for(i=0;i<*neqk;++i){icolk[i]=0;}
    k=0;
    for(i=0;i<nmast;++i){
      if(mast1[i]==irowk[i]){++isubtract;}
      else{
	mast1[i-isubtract]=mast1[i];
	irowk[i-isubtract]=irowk[i];
	if(k!=mast1[i]){
	  for(l=k;l<mast1[i];++l){jqk[l]=i+1-isubtract;}
	  k=mast1[i];
	}
	++icolk[k-1];
      }
    }
    nmast=nmast-isubtract;
    for(l=k;l<*neqk+1;++l){jqk[l]=nmast+1;}
      
    for(i=0;i<*neqk;++i){
      if(jqk[i+1]-jqk[i]>0){
	isize=jqk[i+1]-jqk[i];
	FORTRAN(isortii,(&irowk[jqk[i]-1],&mast1[jqk[i]-1],&isize,&kflag));
      }
    }
      
    if(*neqk==0){*nzsk=0;}
    else{*nzsk=jqk[*neqk]-1;}

  }

  *mast1p=mast1;
  *irowtp=irowt;*irowvp=irowv;*irowpp=irowp;*irowkp=irowk;

  return;

}
