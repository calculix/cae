/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2019 Guido Dhondt                     */

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

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define abs(a) (((a) < (0)) ? (-a) : (a))
/**  transforming the system back to the standard basis functions \f$ \tilde{u}\rightarrow u \f$ and updating the active
 *       and inactive sets and the Langrange Multipliers (LM) 
 *       see Sitzmann Algorithm 2, p.71
 * 
 * Author: Saskia Sitzmann
 *
 * @param [in] 	bhat		intermediate right hand side 
 * @param [in] 	adc		intermediate system matrix, diagonal terms
 * @param [in] 	auc		intermediate system matrix, nondiagonal terms
 * @param [in] 	jqc             pointer to irowc
 * @param [in] 	irowc		row numbers of auc
 * @param [in] 	neq             number of active degrees of freedom
 * @param [in] 	gap		(i) gap for slave node i
 * @param [in,out] b		in: differenzial displacement out:real displacement
 * @param [in,out] islavact	(i) indicates, if slave node i is active (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, =0 inactive node, =1 sticky node, =2 slipping/active node) 
 * @param [out] irowddinv	field containing row numbers of auddinv
 * @param [out] jqddinv		pointer into field irowddinv
 * @param [out] auddinv		coupling matrix \f$ \tilde{D}^{-1}_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [in] irowtloc		field containing row numbers of autloc
 * @param [in] jqtloc	        pointer into field irowtloc
 * @param [in] autloc		transformation matrix \f$ T[p,q]\f$ for slave nodes \f$ p,q \f$ 
 * @param [in] irowtlocinv	field containing row numbers of autlocinv
 * @param [in] jqtlocinv	pointer into field irowtlocinv
 * @param [in] autlocinv	transformation matrix \f$ T^{-1}[p,q]\f$ for slave nodes \f$ p,q \f$  
 * @param [in] ntie		number of contraints 
 * @param [in] nslavnode	(i)pointer into field isalvnode for contact tie i 
 * @param [in] islavnode	field storing the nodes of the slave surface 
 * @param [in] nmastnode	(i)pointer into field imastnode for contact tie i
 * @param [in] imastnode	field storing the nodes of the master surfaces 
 * @param [in] slavnor		slave normals
 * @param [in] slavtan		slave tangents 
 * @param [in] nactdof 		(i,j) actual degree of freedom for direction i of node j
 * @param [out]	iflagact	flag indicating if semi-smooth Newton has converged
 * @param [in,out] cstress	current Lagrange multiplier 
 * @param [in] cstressini	Lagrange multiplier at start of the increment
 * @param [in] mi		(1) max # of integration points per element (2) max degree of freedom per element 
 * @param [out] cdisp		vector saving contact variables for frd-output
 * @param [out] f_cs            contact forces for active degrees of freedom
 * @param [out] f_cm            not used any more
 * @param [in] iit		iteration number of Newton-Raphson iteration 
 * @param [in] iinc		increment number
 * @param [in] vold 		displacements
 * @param [in] vini 		displacements at the start of the increment
 * @param [out] bp		current friction bound
 * @param [in] nk		number of nodes
 * @param [in] nboun2            number of transformed SPCs
 * @param [in] ndirboun2	(i) direction of transformed SPC i 
 * @param [in] nodeboun2        (i) node of transformed SPC i
 * @param [in] xboun2           (i) value of transformed SPC i
 * @param [in] nmpc2		number of transformed mpcs
 * @param [in] ipompc2          (i) pointer to nodempc and coeffmpc for transformed MPC i
 * @param [in] nodempc2         nodes and directions of transformed MPCs
 * @param [in] coefmpc2         coefficients of transformed MPCs
 * @param [in] ikboun2          sorted dofs idof=8*(node-1)+dir for transformed SPCs
 * @param [in] ilboun2          transformed SPC numbers for sorted dofs
 * @param [in] ikmpc2 		sorted dofs idof=8*(node-1)+dir for transformed MPCs
 * @param [in] ilmpc2		transformed SPC numbers for sorted dofs 
 * @param [in] nmpc		number of mpcs
 * @param [in] ipompc           (i) pointer to nodempc and coeffmpc for MPC i
 * @param [in] nodempc          nodes and directions of MPCs
 * @param [in] coefmpc          coefficients of MPCs
 * @param [in] ikboun           sorted dofs idof=8*(node-1)+dir for SPCs
 * @param [in] ilboun           SPC numbers for sorted dofs
 * @param [in] ikmpc 		sorted dofs idof=8*(node-1)+dir for MPCs
 * @param [in] ilmpc		SPC numbers for sorted dofs
 * @param [in] nslavspc2	(2*i) pointer to islavspc2...
 * @param [in] islavspc2         ... which stores transformed SPCs for slave node i
 * @param [in] nsspc2            number of transformed SPC for slave nodes
 * @param [in] nslavmpc2	(2*i) pointer to islavmpc2...
 * @param [in] islavmpc2	... which stores transformed MPCs for slave node i
 * @param [in] nsmpc2		number of transformed MPC for slave nodes 
 * @param [in] nmastspc		(2*i) pointer to imastspc...
 * @param [in] imastspc         ... which stores SPCs for master node i
 * @param [in] nmspc            number of SPC for master nodes
 * @param [in] nmastmpc		(2*i) pointer to imastmpc...
 * @param [in] imastmpc		... which stores MPCs for master node i
 * @param [in] nmmpc		number of MPC for master nodes 
 * @param [in] tieset           (1,i) name of tie constraint (2,i) dependent surface (3,i) independent surface 
 * @param [in] elcon		material parameters
 * @param [in] tietol		(1,i) tie tolerance (2,i) contant interaction material definition
 * @param [in] ncmat_		maximum number of elastic material constants 
 * @param [in] ntmat_           maximum number of temperature data points for any material
 * @param [in] plicon		isotropic hardening curve or points for pressure-overclosure=tabular
 * @param [in] nplicon          isotropic hardening curve. 
 * @param [in] npmat_		maximum number of data points for plicon
 * @param [in] nelcon           (1,i) number of elastic constants for material i (2,i) number of temperature points 
 * @param [in] dtime		delta time
 * @param [out] cfs 		contact force 
 * @param [out] cfm 		not used any more
 * @param [in] islavnodeinv     (i) slave node index for node i
 * @param [out] Bd		coupling matrix \f$ B_d[p,q]=\int \psi_p \phi_q dS \f$, \f$ p \in S, q \in M \f$ 
 * @param [out] irowb		field containing row numbers of Bd
 * @param [out] jqb		pointer into field irowb
 * @param [out] Dd		coupling matrix \f$ D_d[p,q]=\int \psi_p \phi_q dS \f$, \f$ p,q \in S \f$ 
 * @param [out] irowd		field containing row numbers of Dd
 * @param [out] jqd		pointer into field irowd
 * @param [out] Ddtil		coupling matrix \f$ \tilde{D}_d[p,q]=\int \psi_p \tilde{\phi}_q dS \f$, \f$ p,q \in S \f$ 
 * @param [out] irowdtil	field containing row numbers of Ddtil
 * @param [out] jqdtil		pointer into field irowdtil 
 * @param [out] Bpgd		Petrov-Galerkin coupling matrix \f$ B_d^{PG}[p,q]=\int \tilde{\phi}_p \phi_q dS \f$, \f$ p \in S, q \in M \f$ 
 * @param [out] irowbpg		field containing row numbers of Bpgd
 * @param [out] jqbpg		pointer into field irowbpg
 * @param [out] Dpgd		Petrov-Galerkin coupling matrix \f$ D_d[p,q]=\int \tilde{\phi}_p \phi_q dS \f$, \f$ p,q \in S \f$ 
 * @param [out] irowdpg		field containing row numbers of Dpgd
 * @param [out] jqdpg		pointer into field irowdpg 
 * @param [in] lambdaiwan       Lagrange multiplier splitted to Iwan elements
 * @param [in] lambdaiwanini    Lagrange multiplier splitted to Iwan elements at start of increment
 * @param [in] nmethod		analysis method 
 * @param [in] bet		parameter used in alpha-method
 * @param [in]  iflagdualquad   flag indicating what mortar contact is used (=1 quad-lin, =2 quad-quad, =3 PG quad-lin, =4 PG quad-quad)
 * @param [in] ithermal         thermal method
 * @param [in] iperturb		geometrical method
 * @param [in] labmpc           labels of MPCs
 * @param [in] labmpc2		labels of transformed MPCs
 * @param [in] cam              dummy variable here
 * @param [in] veold		velocity
 * @param [in] accold		acceleration
 * @param [in] gam		parameter for alpha method
 * @param [in] nk2		number or generated points needed for transformed SPCs
*/

void contactstress2(double *bhat, double *adc, double *auc,ITG *jqc, 
		    ITG *irowc, ITG *neq, double *gap,
		    double *b, ITG *islavact,
		    ITG *irowddinv, ITG *jqddinv,double *auddinv,
		    ITG *irowtloc, ITG *jqtloc,double *autloc,  
		    ITG *irowtlocinv, ITG *jqtlocinv,double *autlocinv, 
		    ITG *ntie, ITG *nslavnode,
		    ITG *islavnode, ITG *nmastnode, ITG *imastnode, 
		    double *slavnor,double *slavtan,
		    ITG *nactdof,ITG *iflagact,double *cstress,
		    double *cstressini, ITG *mi,
		    double *cdisp, double *f_cs, double *f_cm, ITG *iit, 
		    ITG *iinc,
		    double *vold,double *vini,double* bp,ITG *nk,
		    ITG *nboun2,ITG *ndirboun2,ITG *nodeboun2,double *xboun2,
		    ITG *nmpc2,ITG *ipompc2,ITG *nodempc2,double *coefmpc2,
		    ITG *ikboun2,ITG *ilboun2,ITG *ikmpc2,ITG *ilmpc2,
		    ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,
		    ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
		    ITG *nslavspc2,ITG *islavspc2,ITG *nsspc2,ITG *nslavmpc2,
		    ITG *islavmpc2,ITG *nsmpc2,
		    ITG *nmastspc,ITG *imastspc,ITG *nmspc,ITG *nmastmpc,
		    ITG *imastmpc,ITG *nmmpc,
		    char *tieset,
		    double  *elcon, double *tietol,ITG *ncmat_,ITG *ntmat_,
		    double *plicon,ITG *nplicon, ITG *npmat_,ITG *nelcon,
		    double *dtime,double *cfs, double *cfm,ITG *islavnodeinv,
		    double *Bd, ITG *irowb,ITG *jqb,
		    double *Dd, ITG *irowd,ITG *jqd,
		    double *Ddtil, ITG *irowdtil,ITG *jqdtil,
		    double *Bdtil, ITG *irowbtil,ITG *jqbtil,
		    double *Bpgd, ITG *irowbpg,ITG *jqbpg,
		    double *Dpgd, ITG *irowdpg,ITG *jqdpg,
		    double *lambdaiwan,double *lambdaiwanini, ITG *nmethod, 
		    double *bet,
		    ITG *iflagdualquad,ITG *ithermal,ITG *iperturb,
		    char *labmpc,char *labmpc2,double *cam, double *veold, 
		    double *accold, double *gam,ITG *nk2,
		    double *cfsini,double *cfstil,double *plkcon, 
		    ITG *nplkcon,char *filab,double *f,double *fn,
                    double *qa,ITG *nprint,char *prlab,double *xforc,
                    ITG *nforc){
  
  ITG i,ii,j,l,jj,k,kk,idof1,idof2,idof3,nodes,mt=mi[1]+1,nstick=0,nslip=0,
      ninacti=0,nnogap=0,nolm=0,ndiverg,nhelp,
      debug,idof,node2,dirind,dirdep,index,ist,idofs,islavk2,id,dim,node1,
      imastk2,yielded,number,iout,num_cpus=1,
      keepset,derivmode,regmode,node_max_ncf_n,node_max_ncf_t[2],ndof,
      regmodet,niwan,calcul_fn,calcul_f,intpointvarm,calcul_qa,calcul_cauchy,
      iener,ikin,intpointvart;
  
  double aux,stressnormal,ddispnormal,*unitmatrix=NULL,atau2,stressnormalo,
      contribution,*veold2=NULL,*accold2=NULL,
      constant=1.E10,constantn=1.E10,constantt=1.E10,stresst[2],stressinit[2],
      stressinit2[2],stresstildet[2],stresstildet2[2],ddispt[2],
      disp_totalnormal,disp_tildet[2],
      uholdt1,uholdt2,bpold,nw_t=0.0,*du=NULL,f_csn, f_cmn,dd[3], ch,dhs,
      coefdep,disp_t[2],disp_tildeto[2],scal,
      w_t[3],*bhat2=NULL,epsilon=1.e-5*dtime[0],nut,nstresst,n[3],t[6],
      *fmpc=NULL,lm_t1_av,lm_t2_av,*rc=NULL,heatflux,betac,rheat[3],hheat,
      f_cs_tot[4], f_cm_tot[4], *vectornull=NULL,*u_oldt=NULL,resreg[2],
      rslip[6],ltslip[6],ltu[2],resreg2[2],*tmean=NULL,conductance,
      *vtilt=NULL,*dtemp=NULL,
      *cstress2=NULL,*cstresstil=NULL,*cstressini2=NULL,*u_old=NULL,aninvloc,
      atauinvloc,gnc,gtc[2],*cold=NULL,*cold2=NULL,uturn_t[2],disp_t2[2],
      ln_old,ncf_n, max_ncf_n,ncf_t[2], max_ncf_t[2],mu,mumax,p0,beta,sb,rb,
      reversed,nstresstildet,nstresstildet2,nstressinit,reversed2,*b2=NULL,
      alpha, *rs=NULL,*rm=NULL;
  
  clock_t debut;
  clock_t fin;
  debug=0;
  alpha = 1-2*sqrt(*bet);
  keepset=0;
  mumax=-1.0;
  max_ncf_n=-1.0; max_ncf_t[0]=-1.0; max_ncf_t[1]=-1.0;
  node_max_ncf_n=0;node_max_ncf_t[0]=0;node_max_ncf_t[1]=0;
  lm_t1_av=0;lm_t2_av=0;
  *iflagact=0;
  NNEW(vectornull,double,neq[1]);
  NNEW(cstress2,double,neq[1]);  
  NNEW(bhat2,double,mt*(*nk+*nk2));
  NNEW(b2,double,mt*(*nk+*nk2));
  NNEW(cold,double,3*nslavnode[*ntie]);
  NNEW(cold2,double,3*nslavnode[*ntie]);
  NNEW(cstresstil,double,mt*nslavnode[*ntie]);
  NNEW(cstressini2,double,mt*nslavnode[*ntie]);
  NNEW(du,double,3*nslavnode[*ntie]);
  NNEW(u_old,double,3*nslavnode[*ntie]);
  NNEW(u_oldt,double,3*nslavnode[*ntie]);
  if(ithermal[0]==3){NNEW(tmean,double,nslavnode[*ntie]);NNEW(dtemp,double,nslavnode[*ntie]);}
  ndiverg=14;
  debut = clock();
  
  /* get cstress2 (Sitzmann, Equation (4.14)) */
  FORTRAN(opnonsym,(&neq[1],&aux,b,vectornull,adc,auc,jqc,irowc));
  for(i=0;i<neq[1];i++){f_cs[i]=bhat[i]-vectornull[i];}


  for(i=0;i<neq[1];i++){vectornull[i]=0.0;}	
  FORTRAN(opnonsymt,(&neq[1],&aux,f_cs,cstress2,vectornull,auddinv,jqddinv,irowddinv));
  for(i=0;i<neq[1];i++){f_cs[i]=0.0;}
  
  NNEW(unitmatrix,double,neq[1]);
  // mechanical part
  for(i=0;i<neq[0];i++){
    unitmatrix[i]=1.;vectornull[i]=0.0;
    if(*nmethod == 4){
      cstress2[i]*=1.0/(1.0+alpha);
    }
  }
  // thermal part
  for(i=neq[0];i<neq[1];i++){
    unitmatrix[i]=1.;vectornull[i]=0.0;
    //if(*nmethod == 4){
    //  cstress2[i]*=1.0/(1.0+alpha);
    //}
  }  
  for(i=0;i<mt**nk;i++){bhat2[i]=0.0;}
  iout=1;
  /* fill in missing results from SPCs/MPCs */
  /** @TODO: changes for nmethod==4 , voldtil for incremental mpcs and spc on slave side**/
  
  FORTRAN(resultsini_mortar,(nk,b2,ithermal,iperturb,
			     nactdof,&iout,vold,b,nodeboun2,ndirboun2,
			     xboun2,nboun2,ipompc2,nodempc2,coefmpc2,labmpc2,nmpc2,nmethod,cam,neq,
			     veold,accold,bet,gam,dtime,mi,vini));
/** @TODO: Fix for du in xboun2: dxboun2()=xboun2(i)-vold(ndirboun2(i),nodeboun2(i))**/
    /*NNEW(dxboun2,double,nboun2);
    for(i=0; i<nboun2;i++){
        if(ndirboun[i]<=mi[1]){
	 dxboun2[i]=xboun2[i]-vold[ndirboun2[i],nodeboun2[i]]
	}
    }
    FORTRAN(resultsini,(nk,b2,ithermal,filab,iperturb,f,fn,
       nactdof,&iout,qa,vold,b,nodeboun2,ndirboun2,
       dxboun2,nboun2,ipompc2,nodempc2,coefmpc2,labmpc2,nmpc2,nmethod,cam,neq,
       veold,accold,bet,gam,dtime,mi,vini,nprint,prlab,
       &intpointvarm,&calcul_fn,&calcul_f,&calcul_qa,&calcul_cauchy,&iener,
       &ikin,&intpointvart,xforc,nforc));
      subroutine resultsini(nk,v,ithermal,filab,iperturb,f,fn,
     &  nactdof,iout,qa,vold,b,nodeboun,ndirboun,
     &  xboun,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,
     &  veold,accold,bet,gam,dtime,mi,vini,nprint,prlab,
     &  intpointvarm,calcul_fn,calcul_f,calcul_qa,calcul_cauchy,
     &  ikin,intpointvart,typeboun)       
    FREE(dxboun2);
    */
  /* transformation delta tildeu^q->delta u:
     u=T tildeu^q*/
  for(i=0;i<*nk;i++){
    nodes=i+1;
    for(jj=jqtloc[nodes-1]-1;jj<jqtloc[nodes-1+1]-1;jj++){
      for(l=0;l<3;l++){
	idof1=mt*nodes-3+l;
	idof2=mt*irowtloc[jj]-3+l;
	bhat2[idof2]+=autloc[jj]*b2[idof1];
      }
      if(ithermal[0]==3){
	idof1=mt*nodes-4;
	idof2=mt*irowtloc[jj]-4;
	bhat2[idof2]+=autloc[jj]*b2[idof1];	
      }
      if(debug==1){
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*(i+1)-4],bhat2[mt*(i+1)-4],vold[mt*(i+1)-4]);
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*(i+1)-3],bhat2[mt*(i+1)-3],vold[mt*(i+1)-3]); 
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*(i+1)-2],bhat2[mt*(i+1)-2],vold[mt*(i+1)-2]);
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*(i+1)-1],bhat2[mt*(i+1)-1],vold[mt*(i+1)-1]);
      } 
    }	
  }
 /*if(*nk>55068){
	 i=55068;
	for(k=0;k<3;k++){
	 if(nactdof[mt*i-3+k]>0){printf("node %d dof %d b %e \n",i,k+1,b[nactdof[mt*i-3+k]-1]);}
	} 
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*i-3],bhat2[mt*i-3],vold[mt*i-3]); 
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*i-2],bhat2[mt*i-2],vold[mt*i-2]);
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*i-1],bhat2[mt*i-1],vold[mt*i-1]);  
 }*/
 	/* if(*nk>580){
	   //mid
	 i=381;
	for(k=0;k<3;k++){
	 if(nactdof[mt*i-3+k]>0){printf("node %d dof %d b %e \n",i,k+1,b[nactdof[mt*i-3+k]-1]);}
	} 
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*i-3],bhat2[mt*i-3],vold[mt*i-3]); 
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*i-2],bhat2[mt*i-2],vold[mt*i-2]);
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*i-1],bhat2[mt*i-1],vold[mt*i-1]); 
      	 i=1553;
	for(k=0;k<3;k++){
	 if(nactdof[mt*i-3+k]>0){printf("node %d dof %d b %e \n",i,k+1,b[nactdof[mt*i-3+k]-1]);}
	} 
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*i-3],bhat2[mt*i-3],vold[mt*i-3]); 
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*i-2],bhat2[mt*i-2],vold[mt*i-2]);
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*i-1],bhat2[mt*i-1],vold[mt*i-1]);
      printf("\n\n");
      //corner
	 i=356;
	for(k=0;k<3;k++){
	 if(nactdof[mt*i-3+k]>0){printf("node %d dof %d b %e \n",i,k+1,b[nactdof[mt*i-3+k]-1]);}
	} 
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*i-3],bhat2[mt*i-3],vold[mt*i-3]); 
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*i-2],bhat2[mt*i-2],vold[mt*i-2]);
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*i-1],bhat2[mt*i-1],vold[mt*i-1]); 
      	 i=1550;
	for(k=0;k<3;k++){
	 if(nactdof[mt*i-3+k]>0){printf("node %d dof %d b %e \n",i,k+1,b[nactdof[mt*i-3+k]-1]);}
	} 
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*i-3],bhat2[mt*i-3],vold[mt*i-3]); 
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*i-2],bhat2[mt*i-2],vold[mt*i-2]);
      printf("node %d b2 %e du %e u %e \n",i,b2[mt*i-1],bhat2[mt*i-1],vold[mt*i-1]);
	 } */
  /* overwrite b with untransformed delta u*/
  for( i=0; i<*ntie; i++){      	
    if(tieset[i*(81*3)+80]=='C'){
      nhelp=0;
      for(j=nslavnode[i]; j<nslavnode[i+1]; j++){
	nodes=islavnode[j];
	for(l=0;l<3;l++){
	  b2[mt*nodes-3+l]=bhat2[mt*nodes-3+l];
	  idof1=nactdof[mt*nodes-3+l]-1;
	  if(idof1>-1){	      
	    if(*nmethod==4){
	      b[idof1]=(b2[mt*nodes-3+l])/(*bet**dtime**dtime);
	    }else{
	      b[idof1]=b2[mt*nodes-3+l];
	    }
	  }
	}
	if(ithermal[0]==3){
          b2[mt*nodes-4]=bhat2[mt*nodes-4];
	  idof1=nactdof[mt*nodes-4]-1;
	  if(idof1>-1){	      
	      b[idof1]=b2[mt*nodes-4];
	  }	  
	}
	if(islavact[j]>-1){nhelp++;}
      }
      ndiverg=max(ndiverg,(nhelp/100)+*ntie);
    }
  }
  
  /* calculate hatu=D u^S+ B u^M for update in semi-smooth Newton, see Sitzmann Chapter 3.4. */
  if(*iflagdualquad>2){
    /* Petrov-Galerkin formulation*/  
    /* get du^hat */ 
    /* get uhat_k-1 */
    for(i=0;i<*nk;i++){
      nodes=i+1;
      for(jj=jqdpg[nodes-1]-1;jj<jqdpg[nodes-1+1]-1;jj++){
	for(k=0;k<3;k++){
	  du[(islavnodeinv[irowdpg[jj]-1]-1)*3+k]+=Dpgd[jj]*b2[mt*nodes-3+k];
	  u_oldt[(islavnodeinv[irowdpg[jj]-1]-1)*3+k]+=Dpgd[jj]*(vold[mt*(nodes)-3+k]-vini[mt*(nodes)-3+k]);
	  u_old[(islavnodeinv[irowdpg[jj]-1]-1)*3+k]+=Dpgd[jj]*(vold[mt*(nodes)-3+k]);
	}
      }	    
    }
    for(i=0;i<*nk;i++){
      nodes=i+1;
      for(jj=jqbpg[nodes-1]-1;jj<jqbpg[nodes-1+1]-1;jj++){
	for(k=0;k<3;k++){
	  du[(islavnodeinv[irowbpg[jj]-1]-1)*3+k]+=Bpgd[jj]*b2[mt*nodes-3+k];
	  u_oldt[(islavnodeinv[irowbpg[jj]-1]-1)*3+k]+=Bpgd[jj]*(vold[mt*(nodes)-3+k]-vini[mt*(nodes)-3+k]);
	  u_old[(islavnodeinv[irowbpg[jj]-1]-1)*3+k]+=Bpgd[jj]*(vold[mt*(nodes)-3]+k);
	}
      }	  
    } 
  }else{
    /* normal formulation */  
    /* get du^hat */ 
    /* get uhat_k-1 */
    for(i=0;i<*nk;i++){
      nodes=i+1;
      for(jj=jqd[nodes-1]-1;jj<jqd[nodes-1+1]-1;jj++){
	for(k=0;k<3;k++){
	  du[(islavnodeinv[irowd[jj]-1]-1)*3+k]+=Dd[jj]*b2[mt*nodes-3+k];
	  u_oldt[(islavnodeinv[irowd[jj]-1]-1)*3+k]+=Dd[jj]*(vold[mt*(nodes)-3+k]-vini[mt*(nodes)-3+k]);
	  u_old[(islavnodeinv[irowd[jj]-1]-1)*3+k]+=Dd[jj]*(vold[mt*(nodes)-3+k]);
	}
	if(ithermal[0]==3){
	  tmean[islavnodeinv[irowd[jj]-1]-1]+=Dd[jj]*(b2[mt*nodes-4]+vold[mt*(nodes)-4])*0.5;
	  dtemp[islavnodeinv[irowd[jj]-1]-1]+=Dd[jj]*(b2[mt*nodes-4]+vold[mt*(nodes)-4]);
	}
      }	    
    }
    for(i=0;i<*nk;i++){
      nodes=i+1;
      for(jj=jqb[nodes-1]-1;jj<jqb[nodes-1+1]-1;jj++){
	for(k=0;k<3;k++){
	  du[(islavnodeinv[irowb[jj]-1]-1)*3+k]+=Bd[jj]*b2[mt*nodes-3+k];	
	  u_oldt[(islavnodeinv[irowb[jj]-1]-1)*3+k]+=Bd[jj]*(vold[mt*(nodes)-3+k]-vini[mt*(nodes)-3+k]);
	  u_old[(islavnodeinv[irowb[jj]-1]-1)*3+k]+=Bd[jj]*(vold[mt*(nodes)-3+k]);
	}
	if(ithermal[0]==3){
	  tmean[islavnodeinv[irowb[jj]-1]-1]-=Bd[jj]*(b2[mt*nodes-4]+vold[mt*(nodes)-4])*0.5;
	  dtemp[islavnodeinv[irowb[jj]-1]-1]+=Bd[jj]*(b2[mt*nodes-4]+vold[mt*(nodes)-4]);
	}	
      }	  
    } 
  }
  /* get lambda_scaled */
  
  for (i=0;i<*ntie;i++){
    if(tieset[i*(81*3)+80]=='C'){      
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	stressnormalo=cstress[mt*j]*slavnor[3*j]+cstress[mt*j+1]*slavnor[3*j+1]+cstress[mt*j+2]*slavnor[3*j+2];
	cold[3*j+0]=cstress[mt*j];
	cold[3*j+1]=cstress[mt*j+1];
	cold[3*j+2]=cstress[mt*j+2];       
	nodes=islavnode[j];	    	
	idof1=nactdof[mt*nodes-3]-1;	
	idof2=nactdof[mt*nodes-2]-1;	
	idof3=nactdof[mt*nodes-1]-1; 
	ndof=0;
	if(idof1>-1){ndof++;}if(idof2>-1){ndof++;}if(idof3>-1){ndof++;}
	for(k=0;k<3;k++){	  
	  idof=nactdof[mt*nodes-3+k]-1;	  
	  if(idof>-1) {	  	  	    	  	    
	    cstress[mt*j+k]=cstress2[idof];	  	       	 	  
	  }else{	    	  	    
	    cstress[mt*j+k]=0.0;	     	  	    
	    for(jj=nslavmpc2[2*(j)];jj<nslavmpc2[2*(j)+1];jj++){             	    	      
	      ist=islavmpc2[2*jj];            	    	      
	      dirdep=nodempc2[3*(ist-1)+1];            	    	      
	      coefdep=coefmpc2[ist-1];            	    	      
	      index=nodempc2[3*(ist-1)+2];	     	    	      
	      node1=nodempc2[3*(ist-1)];	     	    	      
	      if(dirdep==k+1){               	      		
		while(index!=0){               				  
		  dirind=nodempc2[3*(index-1)+1];	       				  
		  node2=nodempc2[3*(index-1)];	       				  
		  ch=0.0;
		  if(node2>*nk){
		    idof=-1;
		  }else{
		    idof=nactdof[mt*node2-3+(dirind-1)]-1;
		  }
		  if(idof>-1){	       		 		    
		    ch=cstress2[idof];
		  } 	       
		  cstress[mt*j+dirdep-1]=cstress[mt*j+dirdep-1]-coefmpc2[index-1]*ch/coefdep;	       		             		  
		  index=nodempc2[3*(index-1)+2];	       	      		
		}	      	    	      
	      }	     	  	    
	    }	     	     	  	   		  
	  }		
	}
	if(ithermal[0]==3){
	  idof=nactdof[mt*nodes-4]-1;	  
	  if(idof>-1) {	  	  	    	  	    
	    cstress[mt*j+3]=cstress2[idof];	  	       	 	  
	  }	
	}
      }
    }
  }
  
  /** generate cstressini2, cstresstil **/
  for (i=0;i<*ntie;i++){      	
    if(tieset[i*(81*3)+80]=='C'){      	
      for(j=nslavnode[i];j<nslavnode[i+1];j++){	    
	nodes=islavnode[j];
	for(jj=jqdtil[nodes-1]-1;jj<jqdtil[nodes-1+1]-1;jj++){
	  for(l=0;l<3;l++){
	    cstresstil[(islavnodeinv[irowdtil[jj]-1]-1)*mt+l]+=Ddtil[jj]*cstress[(islavnodeinv[nodes-1]-1)*mt+l];
	    cstressini2[(islavnodeinv[irowdtil[jj]-1]-1)*mt+l]+=Ddtil[jj]*cstressini[(islavnodeinv[nodes-1]-1)*mt+l];
	    cold2[(islavnodeinv[irowdtil[jj]-1]-1)*3+l]+=Ddtil[jj]*cold[(islavnodeinv[nodes-1]-1)*3+l];
	  }
	  if(ithermal[0]==3){
	    cstresstil[(islavnodeinv[irowdtil[jj]-1]-1)*mt+3]+=Ddtil[jj]*cstress[(islavnodeinv[nodes-1]-1)*mt+3];
	    cstressini2[(islavnodeinv[irowdtil[jj]-1]-1)*mt+3]+=Ddtil[jj]*cstressini[(islavnodeinv[nodes-1]-1)*mt+3];	    
	  }
	}	  
      }	
    }    
  } 
  
  /* update slave nodes according to semi-smooth Newton */
  for (i=0;i<*ntie;i++){
    if(tieset[i*(81*3)+80]=='C'){      
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	nodes=islavnode[j];
	idof1=nactdof[mt*nodes-3]-1;	
	idof2=nactdof[mt*nodes-2]-1;	
	idof3=nactdof[mt*nodes-1]-1; 
	ndof=0;
	if(idof1>-1){ndof++;}if(idof2>-1){ndof++;}if(idof3>-1){ndof++;}
	FORTRAN(getcontactparams,(&mu,&regmode,&regmodet,&aninvloc,&atauinvloc,&p0,&beta,tietol,elcon,&i,ncmat_,ntmat_,&niwan));
	mumax=max(mu,mumax);
        stressnormalo=cold2[3*j+0]*slavnor[3*j]+cold2[3*j+1]*slavnor[3*j+1]+cold2[3*j+2]*slavnor[3*j+2];
	stressnormal=cstresstil[mt*j+0]*slavnor[3*j]+cstresstil[mt*j+1]*slavnor[3*j+1]+cstresstil[mt*j+2]*slavnor[3*j+2];	
	stresst[0]=cstresstil[mt*j+0]*slavtan[6*j]+cstresstil[mt*j+1]*slavtan[6*j+1]+cstresstil[mt*j+2]*slavtan[6*j+2];	
	stresst[1]=cstresstil[mt*j+0]*slavtan[6*j+3]+cstresstil[mt*j+1]*slavtan[6*j+4]+cstresstil[mt*j+2]*slavtan[6*j+5];
	stressinit[0]=cstressini2[mt*j+0]*slavtan[6*j]+cstressini2[mt*j+1]*slavtan[6*j+1]+cstressini2[mt*j+2]*slavtan[6*j+2];	
	stressinit[1]=cstressini2[mt*j+0]*slavtan[6*j+3]+cstressini2[mt*j+1]*slavtan[6*j+4]+cstressini2[mt*j+2]*slavtan[6*j+5];  
	stresstildet[0]=(stresst[0]-stressinit[0]);
	stresstildet[1]=(stresst[1]-stressinit[1]);
	nstresstildet=sqrt(stresstildet[0]*stresstildet[0]+stresstildet[1]*stresstildet[1]);
	nstressinit=sqrt(stressinit[0]*stressinit[0]+stressinit[1]*stressinit[1]);
        if(ithermal[0]==3){heatflux=cstresstil[mt*j+3];}
	ddispnormal=du[3*j+0]*slavnor[3*j]+du[3*j+1]*slavnor[3*j+1]+du[3*j+2]*slavnor[3*j+2];	
	ddispt[0]=du[3*j+0]*slavtan[6*j]+du[3*j+1]*slavtan[6*j+1]+du[3*j+2]*slavtan[6*j+2];	
	ddispt[1]=du[3*j+0]*slavtan[6*j+3]+du[3*j+1]*slavtan[6*j+4]+du[3*j+2]*slavtan[6*j+5];			
	//disp_totalnormal=(du[3*j+0]+u_old[j*3])*slavnor[3*j]
	//  +(du[3*j+1]+u_old[j*3+1])*slavnor[3*j+1]
	//  +(du[3*j+2]+u_old[j*3+2])*slavnor[3*j+2];
	disp_totalnormal=(du[3*j+0]+u_oldt[j*3])*slavnor[3*j]
	  +(du[3*j+1]+u_oldt[j*3+1])*slavnor[3*j+1]
	  +(du[3*j+2]+u_oldt[j*3+2])*slavnor[3*j+2];	  
	disp_tildet[0]=(du[3*j+0]+u_oldt[j*3])*slavtan[6*j]
	  +(du[3*j+1]+u_oldt[j*3+1])*slavtan[6*j+1]
	  +(du[3*j+2]+u_oldt[j*3+2])*slavtan[6*j+2];
	disp_tildet[1]=(du[3*j+0]+u_oldt[j*3])*slavtan[6*j+3]
	  +(du[3*j+1]+u_oldt[j*3+1])*slavtan[6*j+4]
	  +(du[3*j+2]+u_oldt[j*3+2])*slavtan[6*j+5];
	disp_tildeto[0]=(u_oldt[j*3])*slavtan[6*j]
	  +(u_oldt[j*3+1])*slavtan[6*j+1]
	  +(u_oldt[j*3+2])*slavtan[6*j+2];
	disp_tildeto[1]=(u_oldt[j*3])*slavtan[6*j+3]
	  +(u_oldt[j*3+1])*slavtan[6*j+4]
	  +(u_oldt[j*3+2])*slavtan[6*j+5];
	disp_t[0]=(du[3*j+0]+u_old[j*3])*slavtan[6*j]
	  +(du[3*j+1]+u_old[j*3+1])*slavtan[6*j+1]
	  +(du[3*j+2]+u_old[j*3+2])*slavtan[6*j+2];
	disp_t[1]=(du[3*j+0]+u_old[j*3])*slavtan[6*j+3]
	  +(du[3*j+1]+u_old[j*3+1])*slavtan[6*j+4]
	  +(du[3*j+2]+u_old[j*3+2])*slavtan[6*j+5];		
	uholdt1=u_old[j*3]*slavtan[6*j]
	  +u_old[j*3+1]*slavtan[6*j+1]
	  +u_old[j*3+2]*slavtan[6*j+2];
	uholdt2=u_old[j*3]*slavtan[6*j+3]
	  +u_old[j*3+1]*slavtan[6*j+4]
	  +u_old[j*3+2]*slavtan[6*j+5];
	ln_old=cold2[3*j+0]*slavnor[3*j]+cold2[3*j+1]*slavnor[3*j+1]+cold2[3*j+2]*slavnor[3*j+2];
	nut=sqrt(disp_tildet[0]*disp_tildet[0]+disp_tildet[1]*disp_tildet[1]);
	nstresst=sqrt(stresst[0]*stresst[0]+stresst[1]*stresst[1]);
	bpold=bp[j];
	
	double gnc_old,dgnc_old;
	derivmode=0;
	if(islavact[j]>-1){scal=Ddtil[jqdtil[nodes-1]-1];}else{scal=0.0;}
	FORTRAN(regularization_gn_c,(&ln_old,&derivmode,&regmode,&gnc_old,&aninvloc,&p0,&beta,elcon,nelcon,&i,ntmat_,
				     plicon,nplicon,npmat_,ncmat_,tietol,&scal));
	derivmode=1;
	FORTRAN(regularization_gn_c,(&ln_old,&derivmode,&regmode,&dgnc_old,&aninvloc,&p0,&beta,elcon,nelcon,&i,ntmat_,
				     plicon,nplicon,npmat_,ncmat_,tietol,&scal));
	derivmode=0;
	FORTRAN(regularization_gn_c,(&stressnormal,&derivmode,&regmode,&gnc,&aninvloc,&p0,&beta,elcon,nelcon,&i,ntmat_,
				     plicon,nplicon,npmat_,ncmat_,tietol,&scal));
	
	derivmode=0;
	FORTRAN(regularization_gt_c,(stresstildet,&derivmode,&regmodet,gtc,&atauinvloc,&p0,&beta,elcon,nelcon,&i,ntmat_,
				     plicon,nplicon,npmat_,ncmat_,tietol));
	constantt=min(constant,1.0/atauinvloc);	
	constantn=constant;
	//if(*nmethod==4 && aninvloc>1.e-10){constantn=1.E8;}
	
/*	if(ithermal[0]==3){
         FORTRAN(heat_conduction_contact,(&stressnormal,vtilt,&mu,&betac,resreg,&derivmode,&regmode,
          &debug,&tmean[j],&scal,elcon,nelcon,&i,ntmat_,plkcon,nplkcon,npmat_,ncmat_,tietol,
          n,n,t,t,&conductance,rheat,&hheat));	  
	  }*/
	
	if(mu>1.E-10){ 	  
	  bp[j]=mu*(stressnormal+constantn*(ddispnormal-gap[j]-gnc));	
	}else{	  
	  bp[j]=(stressnormal+constantn*(ddispnormal-gap[j]-gnc));	
	}
	yielded=0;
	resreg[0]=0.0;resreg[1]=0.0;
	if(regmodet==2 && islavact[j]>0){	  
	  jj=j+1;
	  atau2=1.0/atauinvloc;
	  n[0]=slavnor[3*j];n[1]=slavnor[3*j+1];n[2]=slavnor[3*j+2];
	  for(k=0;k<6;k++){t[k]=slavtan[6*j+k];}
	  derivmode=0;
	  // update lambdaiwan
	  regmode=2;
	  FORTRAN(regularization_slip_iwan,(&stressnormal,
					    disp_tildeto,&bpold,&atau2,resreg2,&derivmode,&regmode,lambdaiwan,
					    lambdaiwanini,&jj,n,t,&mu,rslip,ltslip,ltu,&yielded,&nodes,&debug,&niwan,ddispt));
	  bpold=max(0.0,bp[j]);	
	  // calc resreg
	  regmode=1;
	  FORTRAN(regularization_slip_iwan,(&stressnormal,
					    disp_tildet,&bpold,&atau2,resreg,&derivmode,&regmode,lambdaiwan,
					    lambdaiwanini,&jj,n,t,&mu,rslip,ltslip,ltu,&yielded,&nodes,&debug,&niwan,ddispt));
	}
	disp_t2[0]=stressinit[0]+constantt*disp_tildet[0];
	disp_t2[1]=stressinit[1]+constantt*disp_tildet[1];
	
	w_t[0]=stresst[0]+constantt*(disp_tildet[0]-gtc[0]);	
	w_t[1]=stresst[1]+constantt*(disp_tildet[1]-gtc[1]);		
	nw_t=sqrt(w_t[0]*w_t[0]+w_t[1]*w_t[1]);  
	ncf_t[0]=0.0;ncf_t[1]=0.0;ncf_n=0.0;
	/* evaluate non-linear complementary functions (Sitzmann, Equation (3.54),(3.65)) */ 
	if(mu>1.E-10){
	  ncf_n=mu*stressnormal-max(0.0,bp[j]);
	  if(regmodet==1){
	    if(islavact[j]==0){
	      ncf_t[0]=stresst[0];
	      ncf_t[1]=stresst[1];
	    }else if(islavact[j]==1){
	      ncf_t[0]=disp_tildet[0]-gtc[0];
	      ncf_t[1]=disp_tildet[1]-gtc[1];
	    }else if(islavact[j]==2){
	      ncf_t[0]=stresst[0]-bp[j]*w_t[0]/nw_t;
	      ncf_t[1]=stresst[1]-bp[j]*w_t[1]/nw_t;
	    }else{
	      ncf_t[0]=0.0;ncf_t[1]=0.0; 
	    }
	  }else{
	    if(islavact[j]==2){ 
	      ncf_t[0]=stresst[0]-resreg[0];
	      ncf_t[1]=stresst[1]-resreg[1];
	    }else{
	      ncf_t[0]=stresst[0];
	      ncf_t[1]=stresst[1];     
	    }
	  }
	}else{
	  ncf_n=stressnormal-max(0.0,bp[j]);
	  ncf_t[0]=0.0;ncf_t[1]=0.0;
	}  
	if( ncf_n<0.0)ncf_n=-ncf_n;
	if( ncf_t[0]<0.0)ncf_t[0]=-ncf_t[0];
	if( ncf_t[1]<0.0)ncf_t[1]=-ncf_t[1];
	if(ncf_n>max_ncf_n && ndof>0 && islavact[j]>-1){max_ncf_n=ncf_n;node_max_ncf_n=nodes;}
	if(ncf_t[0]>max_ncf_t[0] && ndof>0 && islavact[j]>-1){max_ncf_t[0]=ncf_t[0];node_max_ncf_t[0]=nodes;}
	if(ncf_t[1]>max_ncf_t[1] && ndof>0 && islavact[j]>-1){max_ncf_t[1]=ncf_t[1];node_max_ncf_t[1]=nodes;}	
	//debug=1;  
	if(debug==1  ){  
	  printf("u(%" ITGFORMAT "): %" ITGFORMAT " %" ITGFORMAT " %" ITGFORMAT " act %" ITGFORMAT " %" ITGFORMAT "\n", nodes, idof1+1,idof2+1, idof3+1,islavact[j],yielded);
	  if(regmodet==2){
	    printf("\t resreg2-lmt %e %e  \n",resreg2[0]-stresst[0],resreg2[1]-stresst[1]);
	  }
	  printf("\t resreg %e %e bp %e \n",resreg[0],resreg[1],bp[j]);
	  printf("\t lmtilde %e %e utilde %e %e gtc %e %e \n",stresstildet[0],stresstildet[1],disp_tildet[0],disp_tildet[1],gtc[0],gtc[1]);
	  printf("\t w_t2 %e %e \n",disp_t2[0],disp_t2[1]);
	  printf("\t uhat(%" ITGFORMAT ") : %e %e %e \n", nodes, du[3*j+0],du[3*j+1], du[3*j+2]);
	  printf("\t u2(%" ITGFORMAT "): %e  \n", nodes,ddispnormal);
	  printf("\t u_tot2(%" ITGFORMAT "): %e %e %e \n", nodes, disp_totalnormal,disp_t[0],disp_t[1]);
	  printf("\t cstress(%" ITGFORMAT "): %e %e %e, actif: %" ITGFORMAT "  \n",nodes,cstress[mt*j+0],cstress[mt*j+1],cstress[mt*j+2],islavact[j]);
	  printf("\t cstresstil(%" ITGFORMAT "): %e %e %e\n",nodes,cstresstil[mt*j+0],cstresstil[mt*j+1],cstresstil[mt*j+2]);
	  printf("\t lm(%" ITGFORMAT ")     : %e %e %e \n",nodes, stressnormal, stresst[0],stresst[1]);		
	  printf("\t newton eq. n: %e = %e diff %e \n",ddispnormal-dgnc_old*stressnormal, gap[j]+gnc_old-dgnc_old*ln_old, ddispnormal-dgnc_old*stressnormal-(gap[j]+gnc_old-dgnc_old*ln_old));
	  printf("\t u_n %e -dgnc_o*ln %e =gap %e +gnc_o %e -dgnc_o*lo %e   \n",ddispnormal,dgnc_old*stressnormal,gap[j],gnc_old, dgnc_old*ln_old);
	  printf("\t ln %e bp/c %e disp-gap %e gnc(ln) %e \n",stressnormal,ddispnormal-gap[j]-gnc,ddispnormal-gap[j],gnc);
	  if(mu>1.E-10){
	    printf("\t uh_t1= %e uh_t2= %e nuh_t= %e \n",disp_tildet[0],disp_tildet[1],sqrt(disp_tildet[0]*disp_tildet[0]+disp_tildet[1]*disp_tildet[1])) ;		
	    printf("\t w_t1= %e w_t2= %e nw_t= %e \n",w_t[0],w_t[1],nw_t);
	    printf("\t ncf_n %e  ncf_t %e %e\n",ncf_n,ncf_t[0],ncf_t[1]);		
	    if(islavact[j]==1){
	      printf("stick: nlmt %e < bp %e and utildetau-lmtildetau %e %e =0? \n",sqrt((stresst[0])*(stresst[0])+(stresst[1])*(stresst[1])),bp[j],disp_tildet[0]-gtc[0],disp_tildet[1]-gtc[1]);  
	    }else if(islavact[j]==2){
	      printf("slip: nlmt %e = bp %e and utildetau-lmtildetau %e %e >0? \n",sqrt((stresst[0])*(stresst[0])+(stresst[1])*(stresst[1])),bp[j],disp_tildet[0]-gtc[0],disp_tildet[1]-gtc[1]);    
	    }
	  }
	  if(ithermal[0]==3){
	   printf("\t heat flux %e tmean %e dtemp %e conductance %e newton %e = %e diff %e\n",heatflux, tmean[j],dtemp[j],conductance, heatflux, dtemp[j]*conductance,heatflux-dtemp[j]*conductance); 
	  }
	}
	debug=0;
	
	if(keepset==0){
	  if(mu>1.E-10){ //contact tie with friction	
	    if (islavact[j]>-1){	
	      if(regmodet==1){
		/*** Active/Inactive set condition cf: Sitzmann Equation (3.58),(3.70) ***/		    	    
		if(bp[j]>-1E-10 && nw_t<(bp[j])){				      
		  nstick++;			
		  if (islavact[j]!=1) {*iflagact = 1;}				
		  islavact[j]=1;
		  cdisp[6*j]=gap[j]-ddispnormal;		  
		  cdisp[6*j+1]=disp_t[0];
		  cdisp[6*j+2]=disp_t[1];							      
		  cdisp[6*j+3]=stressnormal;				      
		  cdisp[6*j+4]=stresst[0];				      
		  cdisp[6*j+5]=stresst[1];
		  lm_t1_av=lm_t1_av+abs(stresst[0]);
		  lm_t2_av=lm_t2_av+abs(stresst[1]);
		}else if(bp[j]>-1E-10 && nw_t>=(bp[j])){				      
		  nslip++;						      
		  if (islavact[j]!=2) {*iflagact = 1;}	
		  islavact[j]=2;
		  cdisp[6*j]=gap[j]-ddispnormal;		  
		  cdisp[6*j+1]=disp_t[0];
		  cdisp[6*j+2]=disp_t[1];						
		  cdisp[6*j+3]=stressnormal;				      
		  cdisp[6*j+4]=stresst[0];				      
		  cdisp[6*j+5]=stresst[1];
		  lm_t1_av=lm_t1_av+abs(stresst[0]);
		  lm_t2_av=lm_t2_av+abs(stresst[1]);		  
		}else{				      
		  if (islavact[j]>0){ *iflagact = 1;}	
		  ninacti++;				      		  			
		  islavact[j]=0;				      
		  cstress[mt*j+0]=0.;				      
		  cstress[mt*j+1]=0.;				      
		  cstress[mt*j+2]=0.;
		  if(ithermal[0]==3){cstress[mt*j+3]=0.;}
		  cdisp[6*j]=gap[j]-ddispnormal;		        	      
		  cdisp[6*j+1]=0.;		        	      
		  cdisp[6*j+2]=0.;		        	      
		  cdisp[6*j+3]=0.;		        	      
		  cdisp[6*j+4]=0.;		        	      
		  cdisp[6*j+5]=0.;
		  if(idof1>-1)cstress2[idof1]=0.0;
		  if(idof2>-1)cstress2[idof2]=0.0;
		  if(idof3>-1)cstress2[idof3]=0.0;
		}
	      }else{
		if(bp[j]>-1E-10){
		  nslip++;						      
		  if (islavact[j]!=2) {*iflagact = 1;}	
		  islavact[j]=2;
		  cdisp[6*j]=gap[j]-ddispnormal;		  
		  cdisp[6*j+1]=disp_t[0];
		  cdisp[6*j+2]=disp_t[1];						
		  cdisp[6*j+3]=stressnormal;				      
		  cdisp[6*j+4]=stresst[0];				      
		  cdisp[6*j+5]=stresst[1];
		  
		}else{
		  if (islavact[j]>0){ *iflagact = 1;}	
		  ninacti++;				      		  			
		  islavact[j]=0;				      
		  cstress[mt*j+0]=0.;				      
		  cstress[mt*j+1]=0.;				      
		  cstress[mt*j+2]=0.;
		  if(ithermal[0]==3){cstress[mt*j+3]=0.;}
		  cdisp[6*j]=gap[j]-ddispnormal;		        	      
		  cdisp[6*j+1]=0.;		        	      
		  cdisp[6*j+2]=0.;		        	      
		  cdisp[6*j+3]=0.;		        	      
		  cdisp[6*j+4]=0.;		        	      
		  cdisp[6*j+5]=0.;
		  if(idof1>-1)cstress2[idof1]=0.0;
		  if(idof2>-1)cstress2[idof2]=0.0;
		  if(idof3>-1)cstress2[idof3]=0.0;
		  if(regmodet==2){
		    for(k=0;k<3*niwan;k++){  
		      lambdaiwan[3*niwan*j+k]=0.0;	      
		    }
		  }
		}
	      }
	    }else{
	      /* nodes without LM contribution */
	      if(islavact[j]==-1){	      
		nnogap++;		    
	      }else if(islavact[j]==-2){	      
		nolm++;	    
	      }	            
	      cstress[mt*j+0]=0.;		    
	      cstress[mt*j+1]=0.;		    
	      cstress[mt*j+2]=0.;
	      if(ithermal[0]==3){cstress[mt*j+3]=0.;}
	      cdisp[6*j]=0.;		    	    
	      cdisp[6*j+1]=0.;		    	    
	      cdisp[6*j+2]=0.;		    	    
	      cdisp[6*j+3]=0.;		    	    
	      cdisp[6*j+4]=0.;		    	    
	      cdisp[6*j+5]=0.;
	      if(idof1>-1)cstress2[idof1]=0.0;
	      if(idof2>-1)cstress2[idof2]=0.0;
	      if(idof3>-1)cstress2[idof3]=0.0;
	      if(regmodet==2){
		for(k=0;k<3*niwan;k++){	
		  lambdaiwan[3*niwan*j+k]=0.0;	      
		}
	      }
	    }   	
	  }else{ //no friction            	  
	    if (islavact[j]>-1){  		    	     
	      /*** Active/Inactive set condition cf: Sitzmann Equation (3.58) ***/		    	    
	      if(bp[j]>-1E-10){	      
		nslip++;				      
		if (islavact[j]!=2) {*iflagact = 1;}
		islavact[j]=2;
		cdisp[6*j]=gap[j]-ddispnormal;
		cdisp[6*j+1]=disp_t[0];
		cdisp[6*j+2]=disp_t[1];						
		cdisp[6*j+3]=stressnormal;			
		cdisp[6*j+4]=stresst[0];			
		cdisp[6*j+5]=stresst[1];			    	    
	      }else{				          
		if (islavact[j]!=0){ *iflagact = 1;}
		ninacti++;				      
		islavact[j]=0;
		cstress[mt*j+0]=0.;		  				      
		cstress[mt*j+1]=0.;				      
		cstress[mt*j+2]=0.;
		if(ithermal[0]==3){cstress[mt*j+3]=0.;}
		cdisp[6*j]=ddispnormal-gap[j];		        	      
		cdisp[6*j+1]=0.;	      
		cdisp[6*j+2]=0;						        	      
		cdisp[6*j+3]=0.;		        	      
		cdisp[6*j+4]=0.;		        	      
		cdisp[6*j+5]=0.;
		if(idof1>-1)cstress2[idof1]=0.0;
		if(idof2>-1)cstress2[idof2]=0.0;
		if(idof3>-1)cstress2[idof3]=0.0;
	      }			  
	    }else{
	      /* nodes without LM contribution */
	      if(islavact[j]==-1){	      
		nnogap++;		    
	      }else{	      
		nolm++;	    
	      }	            	    
	      cstress[mt*j+0]=0.;		    	    
	      cstress[mt*j+1]=0.;		    	    
	      cstress[mt*j+2]=0.;
	      if(ithermal[0]==3){cstress[mt*j+3]=0.;}
	      cdisp[6*j]=0.;			    	    
	      cdisp[6*j+1]=0.;            
	      cdisp[6*j+2]=0.0;			    		    	    
	      cdisp[6*j+3]=0.;		    	    
	      cdisp[6*j+4]=0.;		    	    
	      cdisp[6*j+5]=0.;
	      if(idof1>-1)cstress2[idof1]=0.0;
	      if(idof2>-1)cstress2[idof2]=0.0;
	      if(idof3>-1)cstress2[idof3]=0.0;
	    }			    	
	  }	
	}else{
	  /* in case active set is fixed */
	  if(islavact[j]>0){				     
	    cdisp[6*j]=gap[j]-ddispnormal;
	    cdisp[6*j+1]=disp_t[0];
	    cdisp[6*j+2]=disp_t[1];	      			
	    cdisp[6*j+3]=stressnormal;			
	    cdisp[6*j+4]=stresst[0];			
	    cdisp[6*j+5]=stresst[1];		   
	    
	  }else{
	    cstress[mt*j+0]=0.;			
	    cstress[mt*j+1]=0.;			
	    cstress[mt*j+2]=0.;
	    if(ithermal[0]==3){cstress[mt*j+3]=0.;}
	    cdisp[6*j]=ddispnormal-gap[j];		        
	    cdisp[6*j+1]=0.;
	    cdisp[6*j+2]=0;						        
	    cdisp[6*j+3]=0.;		        
	    cdisp[6*j+4]=0.;		        
	    cdisp[6*j+5]=0.;
	    if(idof1>-1)cstress2[idof1]=0.0;
	    if(idof2>-1)cstress2[idof2]=0.0;
	    if(idof3>-1)cstress2[idof3]=0.0;
	  }	  	
	} 
	/* update weighted dual gap */
	gap[j]=gap[j]-ddispnormal;      
      }       
    }   
  }    
  if(keepset==1)printf("contactstress_fric2: keep active set!!!\n");     
  
  printf("\n max_ncf_n %e node %" ITGFORMAT "\n",max_ncf_n,node_max_ncf_n);
  printf(" max_ncf_t %e %e node %" ITGFORMAT " av %e\n",max_ncf_t[0],max_ncf_t[0]/((lm_t1_av+0.001)/(nstick+nslip+0.001)),node_max_ncf_t[0],((lm_t1_av+0.001)/(nstick+nslip+0.001)));
  printf(" max_ncf_t %e %e node %" ITGFORMAT " av %e\n",max_ncf_t[1],max_ncf_t[1]/((lm_t2_av+0.001)/(nstick+nslip+0.001)),node_max_ncf_t[1],((lm_t2_av+0.001)/(nstick+nslip+0.001))); 
  if(keepset==0){
    /* relative convergence critera for semi-smooth Newton */
    if(max_ncf_n>1.e-3  ){*iflagact=1;}
    if((max_ncf_t[0]/((lm_t1_av+0.001)/(nstick+nslip+0.001))>9.e-4 ||max_ncf_t[1]/((lm_t2_av+0.001)/(nstick+nslip+0.001))>9.e-4) && mumax>1.E-10 ){*iflagact=1;} 
  }
  
  SFREE(unitmatrix);
  /* calculating the  contact forces 
     Ph.D. Thesis Stefan Hartmann eqn. (6.26) */
  // fill cstress for nogap nodes
  for(i=0;i<mt**nk;i++){cfs[i]=0.0;}
  if(*nmethod==4){for(i=0;i<mt**nk;i++){cfstil[i]=0.0;}}
  
  // get contact forces
  for(i=0;i<*nk;i++){
    for(j=jqd[i]-1;j<jqd[i+1]-1;j++){
      for(l=0;l<3;l++){
	cfs[mt*(i+1)-3+l]+=Dd[j]*cstress[mt*(islavnodeinv[irowd[j]-1]-1)+l];
      }
      if(ithermal[0]==3){
	cfs[mt*(i+1)-4]+=Dd[j]*cstress[mt*(islavnodeinv[irowd[j]-1]-1)+3];
      }
    }
  }
  for(i=0;i<*nk;i++){
    for(j=jqb[i]-1;j<jqb[i+1]-1;j++){
      for(l=0;l<3;l++){
	cfs[mt*(i+1)-3+l]+=Bd[j]*cstress[mt*(islavnodeinv[irowb[j]-1]-1)+l];
      }
      if(ithermal[0]==3){
	cfs[mt*(i+1)-4]+=Bd[j]*cstress[mt*(islavnodeinv[irowb[j]-1]-1)+3];
      }      
    }
  }
  if(*nmethod==4){
    for(i=0;i<*nk;i++){
      for(j=jqdtil[i]-1;j<jqdtil[i+1]-1;j++){
	for(l=0;l<3;l++){
	  cfstil[mt*(i+1)-3+l]+=Ddtil[j]*cstress[mt*(islavnodeinv[irowdtil[j]-1]-1)+l];
	}
      }
    }
    for(i=0;i<*nk;i++){
      for(j=jqbtil[i]-1;j<jqbtil[i+1]-1;j++){
	for(l=0;l<3;l++){
	  cfstil[mt*(i+1)-3+l]+=Bdtil[j]*cstress[mt*(islavnodeinv[irowbtil[j]-1]-1)+l];
	}
	
      }
    }
    
    for(i=0;i<mt**nk;i++){cfsini[i]=0.0;}
      for(i=0;i<*nk;i++){
      for(j=jqd[i]-1;j<jqd[i+1]-1;j++){
      for(l=0;l<3;l++){
      //cfs[3*(islavnodeinv[i]-1)+l]+=Dd[j]*cstress[3*(islavnodeinv[irowd[j]-1]-1)+l];
      cfsini[mt*(i+1)-3+l]+=Dd[j]*cstressini[mt*(islavnodeinv[irowd[j]-1]-1)+l];
      }
      }
      }
      for(i=0;i<*nk;i++){
      for(j=jqb[i]-1;j<jqb[i+1]-1;j++){
      //prinft("\t r %" ITGFORMAT " c %" ITGFORMAT " %" ITGFORMAT " %" ITGFORMAT "  \n",i+1,irowb[j],-islavnodeinv[i],islavnodeinv[irowd[j]-1])
      for(l=0;l<3;l++){
      //cfm[3*((-islavnodeinv[i])-1)+l]+=Bd[j]*cstress[3*(islavnodeinv[irowb[j]-1]-1)+l];
      cfsini[mt*(i+1)-3+l]+=Bd[j]*cstressini[mt*(islavnodeinv[irowb[j]-1]-1)+l];
	}
	}
	} 
    
  }
  
  NNEW(fmpc,double,*nmpc);
  NNEW(rc,double,*nk*mt);
  calcul_fn=1;
  calcul_f=1;
  for(i=0;i<*nk;i++){
    for(l=0;l<3;l++){      
    if(*nmethod==4){
      /* modification for dynamic calculations */
      /* Sitzmann Equation (5.12) */
      rc[mt*(i+1)-3+l]=cfs[mt*(i+1)-3+l]*(1+alpha)-cfsini[mt*(i+1)-3+l]*(alpha);
    }else{
      rc[mt*(i+1)-3+l]=cfs[mt*(i+1)-3+l];
    }
    }
    if(ithermal[0]==3){
      rc[mt*(i+1)-4]=cfs[mt*(i+1)-4];
      //printf("rc %e\n",rc[mt*(i+1)-4]);
    }
  }
  resultsforc(nk,f_cs,rc,nactdof,ipompc,nodempc,
	      coefmpc,labmpc,nmpc,mi,fmpc,&calcul_fn,
	      &calcul_f,&num_cpus);
  SFREE(fmpc);
  
  /* print total contact force per contact tie for debugging purpose */
  for( i=0; i<*ntie; i++){
    if(tieset[i*(81*3)+80]=='C'){
      for(jj=0;jj<3;jj++){	
	f_cs_tot[jj]=0.0;	
	f_cm_tot[jj]=0.0;     
      }
      if(ithermal[0]==3){
      	f_cs_tot[3]=0.0;	
	f_cm_tot[3]=0.0; 
      }
      for(j=nslavnode[i]; j<nslavnode[i+1]; j++){	  
	nodes = islavnode[j];	  
	for(l=0;l<3;l++){	    
	  idof1=nactdof[mt*nodes-3+l]-1;
	  //cfs[3*j+l]=0.0;
	  if(idof1>-1){
	    f_cs_tot[l]+=f_cs[idof1];
	    f_cs_tot[l]+=f_cm[idof1];
	  }
	}
	if(ithermal[0]==3){
	  idof1=nactdof[mt*nodes-4]-1;
	  //cfs[3*j+l]=0.0;
	  if(idof1>-1){
	    f_cs_tot[3]+=f_cs[idof1];
	    f_cs_tot[3]+=f_cm[idof1];
	  }	  
	}
      }			
      for(j=nmastnode[i]; j<nmastnode[i+1]; j++){	  
	nodes = imastnode[j];	    
	for(l=0;l<3;l++){		
	  idof1=nactdof[mt*nodes-3+l]-1;
	  if(idof1>-1){		 
	    f_cm_tot[l]+=f_cs[idof1];
	    
	  }	    
	}
	if(ithermal[0]==3){
	  idof1=nactdof[mt*nodes-4]-1;
	  //cfs[3*j+l]=0.0;
	  if(idof1>-1){
	    f_cm_tot[3]+=f_cs[idof1];
	  }	  
	}	
      }
      f_csn=sqrt(f_cs_tot[0]*f_cs_tot[0]+f_cs_tot[1]*f_cs_tot[1]+f_cs_tot[2]*f_cs_tot[2]);      
      f_cmn=sqrt(f_cm_tot[0]*f_cm_tot[0]+f_cm_tot[1]*f_cm_tot[1]+f_cm_tot[2]*f_cm_tot[2]); 
      if(ithermal[0]==3){
      printf("\n tie %" ITGFORMAT " slave contact force : %e %e %e and norm: %e heatflux %e \n",i+1,f_cs_tot[0],f_cs_tot[1], f_cs_tot[2], f_csn,f_cs_tot[3]);      
      printf(" tie %" ITGFORMAT " master contact force: %e %e %e and norm: %e heatflux %e \n",i+1,f_cm_tot[0],f_cm_tot[1], f_cm_tot[2], f_cmn,f_cm_tot[3]);	
      }else{
      printf("\n tie %" ITGFORMAT " slave contact force : %e %e %e and norm: %e \n",i+1,f_cs_tot[0],f_cs_tot[1], f_cs_tot[2], f_csn);      
      printf(" tie %" ITGFORMAT " master contact force: %e %e %e and norm: %e \n",i+1,f_cm_tot[0],f_cm_tot[1], f_cm_tot[2], f_cmn );
      }
    } 	      
  }
  
  /* transform contact traction for frd-output I(LM)=T*LM */
  for( i=0; i<*ntie; i++){
    if(tieset[i*(81*3)+80]=='C'){
      for(j=nslavnode[i]; j<nslavnode[i+1]; j++){
	for(jj=0;jj<mt;jj++){	
	  cstresstil[mt*j+jj]=0.0;	    
	}
      }
    }
  }
  for(i=0;i<*nk;i++){
    nodes=i+1;
    for(jj=jqtloc[nodes-1]-1;jj<jqtloc[nodes-1+1]-1;jj++){
      for(l=0;l<3;l++){
	idof1=mt*(islavnodeinv[nodes-1]-1)+l;
	idof2=mt*(islavnodeinv[irowtloc[jj]-1]-1)+l;
	cstresstil[idof2]+=autloc[jj]*cstress[idof1];
      }
    }	
  }
  for( i=0; i<*ntie; i++){      	
    if(tieset[i*(81*3)+80]=='C'){	        
      for(j=nslavnode[i]; j<nslavnode[i+1]; j++){
	stressnormal=cstresstil[mt*j+0]*slavnor[3*j]+cstresstil[mt*j+1]*slavnor[3*j+1]+cstresstil[mt*j+2]*slavnor[3*j+2];	
	stresst[0]=cstresstil[mt*j+0]*slavtan[6*j]+cstresstil[mt*j+1]*slavtan[6*j+1]+cstresstil[mt*j+2]*slavtan[6*j+2];	
	stresst[1]=cstresstil[mt*j+0]*slavtan[6*j+3]+cstresstil[mt*j+1]*slavtan[6*j+4]+cstresstil[mt*j+2]*slavtan[6*j+5];
	cdisp[6*j+3]=stressnormal;			
	cdisp[6*j+4]=stresst[0];			
	cdisp[6*j+5]=stresst[1];
      }
    }
  }
  
  /* needed for adaptive time stepping */
  if(*iit>ndiverg){*iflagact=0;} 
  printf("\n contacstress : N_stick : %" ITGFORMAT "\t N_slip : %" ITGFORMAT "\tN_Inactiv : %" ITGFORMAT "\t N_nogap : %" ITGFORMAT "  N_nolm : %" ITGFORMAT "\n Flag = %" ITGFORMAT "\n",nstick,nslip,ninacti,nnogap,nolm,*iflagact);
  fin= clock();
  printf(" contactstress : ndiverg %" ITGFORMAT " \n",ndiverg);
  printf(" contactstress : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);
  
  number=1;
  dim=neq[1];   
  
  
  SFREE(bhat2);
  SFREE(b2);
  SFREE(cstress2);
  SFREE(u_old);
  SFREE(u_oldt);
  SFREE(du);
  SFREE(cstresstil);
  SFREE(cstressini2);
  SFREE(cold);
  SFREE(cold2);
  SFREE(vectornull);
  SFREE(rc);
  if(ithermal[0]==3){SFREE(tmean);SFREE(dtemp);}
  return;
}
