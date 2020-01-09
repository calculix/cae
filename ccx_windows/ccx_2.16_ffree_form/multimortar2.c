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
#include <time.h>
#include <string.h> 
#include "CalculiX.h"
#include "mortar.h"

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
/**
 * \brief embedding the contact conditions into the matrix system
 *        and condening the Lagrange multiplier
 *        see phd-thesis Sitzmann equation (4.15) for quad-quad/quad-lin method  or (4.24) for PG quad-lin method
 * Author: Saskia Sitzmann
 *
 * @param [in,out] aup		nondiagonal entries of matrix \f$ K \f$
 * @param [in,out] ad		diagonal entries of matrix \f$ K \f$
 * @param [in,out] irowp		field containing row number for entries in au
 * @param [in,out] jq		(i) first element in au belonging to column i	
 * @param [in,out] nzs		size of au 
 * @param [out] aucp		nondiagonal entries of intermediate matrix \f$ K_c \f$
 * @param [out] adc		diagonal entries of intermediate matrix \f$ K_c \f$
 * @param [out] irowcp		field containing row number for entries in auc
 * @param [out] jqc		(i) first element in auc belonging to column i	
 * @param [out] nzsc		size of auc 
 * @param [in] aubd		coupling matrix \f$ B_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [in] irowbd		field containing row numbers of aubd
 * @param [in] jqbd		pointer into field irowbd
 * @param [in] aubdtil		matrix \f$ \tilde{D}^{-1}\tilde{B}_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [in] irowbdtil	field containing row numbers of aubd
 * @param [in] jqbdtil		pointer into field irowbdtil
 * @param [in] aubdtil2		coupling matrix \f$ \tilde{D}$ and $\tilde{B}^2_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [in] irowbdtil2	field containing row numbers of aubdtil2
 * @param [in] jqbdtil2		pointer into field irowbdtil2
 * @param [in] irowdd		field containing row numbers of audd
 * @param [in] jqdd		pointer into field irowdd
 * @param [in] audd		coupling matrix \f$ D_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [in] irowddtil2	field containing row numbers of audd
 * @param [in] jqddtil2		pointer into field irowdd
 * @param [in] auddtil2		matrix \f$ Id_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [in] irowddinv	field containing row numbers of auddinv
 * @param [in] jqddinv		pointer into field irowddinv
 * @param [in] auddinv		coupling matrix \f$ \tilde{D}^{-1}_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [in] Bd		coupling matrix \f$ B_d[p,q]=\int \psi_p \phi_q dS \f$, \f$ p \in S, q \in M \f$ 
 * @param [in] irowb		field containing row numbers of Bd
 * @param [in] jqb		pointer into field irowb
 * @param [in] Dd		coupling matrix \f$ D_d[p,q]=\int \psi_p \phi_q dS \f$, \f$ p,q \in S \f$ 
 * @param [in] irowd		field containing row numbers of Dd
 * @param [in] jqd		pointer into field irowd
 * @param [in] Ddtil		coupling matrix \f$ \tilde{D}_d[p,q]=\int \psi_p \tilde{\phi}_q dS \f$, \f$ p,q \in S \f$ 
 * @param [in] irowdtil	field containing row numbers of Ddtil
 * @param [in] jqdtil		pointer into field irowdtil 
 * @param [in] neq		(0) # of mechanical equations (1) sum of mechanical and thermal equations (2) neq(1+ # of single point contraints)
 * @param [in,out] b		right hand side
 * @param [out] bhat		intermediate right hand side
 * @param [in] islavnode	field storing the nodes of the slave surface
 * @param [in] imastnode	field storing the nodes of the master surfaces
 * @param [in] nslavnode	(i)pointer into field isalvnode for contact tie i 
 * @param [in] nmastnode	(i)pointer into field imastnode for contact tie i
 * @param [in] islavact		(i) indicates, if slave node i is active (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, =0 inactive node, =1 sticky node, =2 slipping/active node) 
 * @param [in] islavactdof      (i)=10*slavenodenumber+direction for active dof i
 * @param [in] gap		(i) \f$ g_i= <g, \Psi_i> \f$ for node i on slave surface
 * @param [in] slavnor		slave normal
 * @param [in] slavtan		slave tangent  
 * @param [in] vold		displacement of nodes
 * @param [in] vini		displacement at start of the increment
 * @param [in] cstress		current Lagrange multiplier 
 * @param [in] cstressini	Lagrange multiplier at start of the increment
 * @param [in] bp_old		old friction bounds 
 * @param [in] nactdof 		(i,j) actual degree of freedom for direction i of node j
 * @param [in] ntie		number of ties
 * @param [in] mi		(1) max # of integration points per element (2) max degree of freedom per element
 * @param [in] nk		number of nodes
 * @param [in] nboun            number of SPCs
 * @param [in] ndirboun		(i) direction of SPC i 
 * @param [in] nodeboun         (i) node of SPC i
 * @param [in] xboun            (i) value of SPC i
 * @param [in] nmpc		number of mpcs
 * @param [in] ipompc           (i) pointer to nodempc and coeffmpc for MPC i
 * @param [in] nodempc          nodes and directions of MPCs
 * @param [in] coefmpc          coefficients of MPCs
 * @param [in] ikboun           sorted dofs idof=8*(node-1)+dir for SPCs
 * @param [in] ilboun           SPC numbers for sorted dofs
 * @param [in] ikmpc 		sorted dofs idof=8*(node-1)+dir for MPCs
 * @param [in] ilmpc		SPC numbers for sorted dofs
 * @param [in] nslavspc		(2*i) pointer to islavspc...
 * @param [in] islavspc         ... which stores SPCs for slave node i
 * @param [in] nsspc            number of SPC for slave nodes
 * @param [in] nslavmpc		(2*i) pointer to islavmpc...
 * @param [in] islavmpc		... which stores MPCs for slave node i
 * @param [in] nsmpc		number of MPC for slave nodes
 * @param [in] nmastspc		(2*i) pointer to imastspc...
 * @param [in] imastspc         ... which stores SPCs for master node i
 * @param [in] nmspc            number of SPC for master nodes
 * @param [in] nmastmpc		(2*i) pointer to imastmpc...
 * @param [in] imastmpc		... which stores MPCs for master node i
 * @param [in] nmmpc		number of MPC for master nodes
 * @param [in] tieset           (1,i) name of tie constraint (2,i) dependent surface (3,i) independent surface
 * @param [in] islavactdoftie   (i)=tie number for active dof i
 * @param [in] nelcon           (1,i) number of elastic constants for material i (2,i) number of temperature points  
 * @param [in] elcon		material parameters 
 * @param [in] tietol		(1,i) tie tolerance (2,i) contant interaction material definition
 * @param [in] ncmat_		maximum number of elastic material constants 
 * @param [in] ntmat_           maximum number of temperature data points for any material
 * @param [in] plicon		isotropic hardening curve or points for pressure-overclosure=tabular
 * @param [in] nplicon          isotropic hardening curve. 
 * @param [in] npmat_		maximum number of data points for plicon
 * @param [in] dtime		delta time 
 * @param [in] irowtloc		field containing row numbers of autloc
 * @param [in] jqtloc	        pointer into field irowtloc
 * @param [in] autloc		transformation matrix \f$ T[p,q]\f$ for slave nodes \f$ p,q \f$ 
 * @param [in] irowtlocinv	field containing row numbers of autlocinv
 * @param [in] jqtlocinv	pointer into field irowtlocinv
 * @param [in] autlocinv	transformation matrix \f$ T^{-1}[p,q]\f$ for slave nodes \f$ p,q \f$  
 * @param [in] islavnodeinv     (i) slave node index for node i 
 * @param [in] lambdaiwan       Lagrange multiplier splitted to Iwan elements
 * @param [in] lambdaiwanini    Lagrange multiplier splitted to Iwan elements at start of increment
 * @param [in] iit		iteration number of Newton-Raphson iteration
 * @param [in] nmethod		analysis method
 * @param [in] bet		parameter used in alpha-method 
*/

void multimortar2(double **aup, double *ad, ITG **irowp, ITG *jq, ITG *nzs,
		  double **aucp, double *adc, ITG **irowcp, ITG *jqc, ITG *nzsc,
		  double *aubd,ITG *irowbd,ITG *jqbd,
		  double *aubdtil,ITG *irowbdtil,ITG *jqbdtil,
		  double *aubdtil2,ITG *irowbdtil2,ITG *jqbdtil2,
		  ITG *irowdd, ITG *jqdd, double *audd,
		  ITG *irowddtil2, ITG *jqddtil2, double *auddtil2,
		  ITG *irowddinv, ITG *jqddinv, double *auddinv,
		  double *Bd, ITG *irowb,ITG *jqb,
		  double *Dd, ITG *irowd,ITG *jqd,
		  double *Ddtil, ITG *irowdtil,ITG *jqdtil,
		  ITG *neq,double *b, double *bhat,ITG *islavnode,ITG *imastnode,
		  ITG *nslavnode,ITG *nmastnode,
		  ITG *islavact,ITG *islavactdof,
		  double *gap,
		  double *slavnor, double *slavtan,
		  double *vold,double *vini, double *cstress, double *cstressini,
		  double *bp_old,ITG *nactdof, ITG *ntie, ITG *mi,ITG *nk,
		  ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
		  ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,
		  ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
		  ITG *nslavspc,ITG *islavspc,ITG *nsspc,ITG *nslavmpc,ITG *islavmpc,ITG *nsmpc,
		  ITG *nmastspc,ITG *imastspc,ITG *nmspc,ITG *nmastmpc,ITG *imastmpc,ITG *nmmpc,
		  char *tieset,
		  ITG *islavactdoftie,ITG *nelcon, double  *elcon, double *tietol,ITG *ncmat_,ITG *ntmat_,
		  double *plicon,ITG *nplicon, ITG *npmat_,double *dtime,
		  ITG *irowtloc, ITG *jqtloc,double *autloc,  
		  ITG *irowtlocinv, ITG *jqtlocinv,double *autlocinv,
		  ITG *islavnodeinv,double *lambdaiwan,double *lambdaiwanini,ITG *iit, ITG *nmethod, double *bet, ITG *ithermal,
		  double *plkcon, ITG *nplkcon
		  ){ 
  
  
  ITG i,j,k,l,ll,m,icol,mt=mi[1]+1,nodesf,nodem,row_ln,row_lm,row_li,row_la,kflag=2,
    nmtrue,nstrue,nntrue,debug,dofm,dofs,jrow,islavnodeentry,impclack,
    idof1,idof2,idof3,node,jrow2,nslavs,nmasts,dim,
    nzs_nn,nzs_nm,nzs_ni,nzs_na,
    nzs_mn,nzs_mm,nzs_mi,nzs_ma,
    nzs_in,nzs_im,nzs_ii,nzs_ia,
    nzs_an,nzs_am,nzs_ai,nzs_aa,
    nzs_bd1,nzs_bdtil2,nzs_ddtil2i,nzs_ddtil2a,nzsbdtil,nzsbdtil2,
    nzs_mmf,nzs_iif,nzs_aaf,nzs_intmn,nzs_intmm,nzs_intmi,nzs_intma,
    nzs_mntil,nzs_mmtil,nzs_mitil,nzs_matil,nzsdd,nzs_intan1,nzs_intan2,
    nzs_intam1,nzs_intam2,nzs_intai1,nzs_intai2,nzs_intaa1,nzs_intaa2,
    nzs_intin1,nzs_intin2,nzsddinv,
    nzs_intim1,nzs_intim2,nzs_intii1,nzs_intii2,nzs_intia1,nzs_intia2,
    nzs_ba,nzsddtil2,nzs_diiinv,nzs_diainv,nzs_daiinv,nzs_daainv,nzs_da,nzs_di,
    nzs_dan,nzs_dai,nzs_dam,nzs_daa,
    nzs_din,nzs_dii,nzs_dim,nzs_dia,
    ifree_nn,ifree_nm,ifree_ni,ifree_na,
    ifree_mn,ifree_mm,ifree_mi,ifree_ma,
    ifree_in,ifree_im,ifree_ii,ifree_ia,
    ifree_an,ifree_am,ifree_ai,ifree_aa,
    ifree, numb,	
    *irow=NULL,*irow_iid=NULL,
    *irow_nn=NULL,*irow_nm=NULL,*irow_ni=NULL,*irow_na=NULL,
    *irow_mn=NULL,*irow_mm=NULL,*irow_mi=NULL,*irow_ma=NULL,
    *irow_in=NULL,*irow_im=NULL,*irow_ii=NULL,*irow_ia=NULL,
    *irow_an=NULL,*irow_am=NULL,*irow_ai=NULL,*irow_aa=NULL,
    *jq_iid=NULL,
    *jq_nn=NULL,*jq_nm=NULL,*jq_ni=NULL,*jq_na=NULL,
    *jq_mn=NULL,*jq_mm=NULL,*jq_mi=NULL,*jq_ma=NULL,
    *jq_in=NULL,*jq_im=NULL,*jq_ii=NULL,*jq_ia=NULL,
    *jq_an=NULL,*jq_am=NULL,*jq_ai=NULL,*jq_aa=NULL,
    *irowc=NULL,*irow_bdtil2=NULL,*jq_bdtil2=NULL,
    *irow_ddtil2i=NULL,*jq_ddtil2i=NULL,*irow_ddtil2a=NULL,*jq_ddtil2a=NULL,
    *irow_bd1=NULL,*jq_bd1=NULL,
    *irow_aad=NULL,*irow_mmd=NULL,*jq_aad=NULL,*jq_mmd=NULL,
    *irow_mmf=NULL,*jq_mmf=NULL,*irow_mmtil=NULL,*jq_mmtil=NULL,
    *irow_mntil=NULL,*jq_mntil=NULL,*irow_mitil=NULL,*jq_mitil=NULL,
    *irow_matil=NULL,*jq_matil=NULL,
    *irow_intmn=NULL,*jq_intmn=NULL,*irow_intmm=NULL,*jq_intmm=NULL,
    *irow_intmi=NULL,*jq_intmi=NULL,*irow_intma=NULL,*jq_intma=NULL,
    *irow_iif=NULL,*jq_iif=NULL,*irow_aaf=NULL,*jq_aaf=NULL,
    *irow_ba=NULL, *jq_ba=NULL,*irow_diainv=NULL,*jq_diainv=NULL,*irow_daainv=NULL,*jq_daainv=NULL,
    *irow_daiinv=NULL,*jq_daiinv=NULL,*irow_diiinv=NULL,*jq_diiinv=NULL,*irow_di=NULL,*jq_di=NULL,
    *irow_da=NULL,*jq_da=NULL,
    *irow_intan1=NULL,*jq_intan1=NULL,*irow_intan2=NULL,*jq_intan2=NULL,
    *irow_intam1=NULL,*jq_intam1=NULL,*irow_intam2=NULL,*jq_intam2=NULL,
    *irow_intai1=NULL,*jq_intai1=NULL,*irow_intai2=NULL,*jq_intai2=NULL,
    *irow_intaa1=NULL,*jq_intaa1=NULL,*irow_intaa2=NULL,*jq_intaa2=NULL,
    *irow_intin1=NULL,*jq_intin1=NULL,*irow_intin2=NULL,*jq_intin2=NULL,
    *irow_intim1=NULL,*jq_intim1=NULL,*irow_intim2=NULL,*jq_intim2=NULL,
    *irow_intii1=NULL,*jq_intii1=NULL,*irow_intii2=NULL,*jq_intii2=NULL,
    *irow_intia1=NULL,*jq_intia1=NULL,*irow_intia2=NULL,*jq_intia2=NULL,
    *irow_dan=NULL,*jq_dan=NULL,*irow_dam=NULL,*jq_dam=NULL,
    *irow_dai=NULL,*jq_dai=NULL,*irow_daa=NULL,*jq_daa=NULL,
    *irow_din=NULL,*jq_din=NULL,*irow_dim=NULL,*jq_dim=NULL,
    *irow_dii=NULL,*jq_dii=NULL,*irow_dia=NULL,*jq_dia=NULL,
    *irow_antil=NULL,*jq_antil=NULL,*irow_amtil=NULL,*jq_amtil=NULL,
    *irow_aitil=NULL,*jq_aitil=NULL,*irow_aatil=NULL,*jq_aatil=NULL,
    *irow_t=NULL,*jq_t=NULL,nzs_t=*nzs,*mast1=NULL,
    *l_flag=NULL,*n_flag=NULL,*m_flag=NULL,*a_flag=NULL,*i_flag=NULL,number=1,
    iact,*n_flagr=NULL,*m_flagr=NULL,*i_flagr=NULL,*a_flagr=NULL,
    nzsddtil2a,nzsddtil2i,nzs_sym;
  
  double help,e1,e2,e3,
    *auc=NULL,*au=NULL,*ad_nn=NULL,*ad_mm=NULL,
    *ad_ii=NULL,*ad_aa,*au_bd1=NULL,*au_bdtil2=NULL,*au_ddtil2i=NULL,*au_ddtil2a=NULL,*au_ssd=NULL,
    *au_nn=NULL,*au_nm=NULL,*au_ni=NULL,*au_na=NULL,
    *au_mn=NULL,*au_mm=NULL,*au_mi=NULL,*au_ma=NULL,
    *au_in=NULL,*au_im=NULL,*au_ii=NULL,*au_ia=NULL,
    *au_an=NULL,*au_am=NULL,*au_ai=NULL,*au_aa=NULL,
    *au_mmf=NULL,*au_mmtil=NULL,*au_mntil=NULL,
    *au_mitil=NULL,*au_matil=NULL,*au_aaf=NULL,*au_iif=NULL,
    *au_intmn=NULL,*au_intmm=NULL,*au_intmi=NULL,*au_intma=NULL,
    *au_ba=NULL, *au_daainv=NULL,*au_diiinv=NULL,*au_diainv=NULL,*au_daiinv=NULL,*au_di=NULL,*au_da=NULL,
    *au_intan1=NULL,*au_intan2=NULL,*au_intam1=NULL,*au_intam2=NULL,
    *au_intai1=NULL,*au_intai2=NULL,*au_intaa1=NULL,*au_intaa2=NULL,
    *au_intin1=NULL,*au_intin2=NULL,*au_intim1=NULL,*au_intim2=NULL,
    *au_intii1=NULL,*au_intii2=NULL,*au_intia1=NULL,*au_intia2=NULL,
    *au_dan=NULL,*au_dam=NULL,*au_dai=NULL,*au_daa=NULL,
    *au_din=NULL,*au_dim=NULL,*au_dii=NULL,*au_dia=NULL,
    *au_antil=NULL,*au_amtil=NULL,*au_aitil=NULL,*au_aatil=NULL,
    *au_t=NULL,	
    *f_a=NULL,*f_m=NULL,*f_i=NULL,*f_inta1=NULL,*f_inta2=NULL,*v_r=NULL,*f_da=NULL,*f_di=NULL,*f_atil=NULL;
  
  
  clock_t debut;
  clock_t fin; 
  irowc = *irowcp; auc=*aucp;
  irow = *irowp; au=*aup;
  nslavs=nslavnode[*ntie];
  nmasts=nmastnode[*ntie];                
  debug=0;
  /* save full K matrix without contact conditions
     needed for calculation of Lagrange multiplier in contactstress2 */
  /* Au is symmetric compute the non_symmeric whole au_matrix auc=au+au^T **/		
  NNEW(au_t,double,nzs_t);	
  NNEW(irow_t,ITG, nzs_t);	
  NNEW(jq_t,ITG, neq[1]+1);	
  dim=neq[1];
  transpose(au,jq,irow,&dim,
	    au_t,jq_t,irow_t);
  
  
  add_rect(au,irow,jq,neq[1],neq[1],
	   au_t,irow_t,jq_t,neq[1],neq[1],
	   &auc,&irowc,jqc,nzsc);
  for(i=0;i<neq[1];i++){adc[i]=ad[i];}
  *nzsc=jqc[neq[1]]-1;
  RENEW(auc,double,*nzsc);
  RENEW(irowc,ITG, *nzsc);
  
  
  /* Flag to produce the bijection between local and global dof **/  
  /*
   *  l_flag[i]=1 for M Master dof
   *  l_flag[i]=2 for I Slave dof
   *  l_flag[i]=3 for A Slave dof
   *  l_flag[i]=0 for N rest of the dof
   *  
   *  n_flag contains local N_row number
   *  m_flag contains local M_row number
   *  i_flag contains local I_row number
   *  a_flag contains local A_row number
  **/ 	
  NNEW(l_flag,ITG, neq[1]);	
  NNEW(n_flag,ITG, neq[1]);	
  NNEW(m_flag,ITG, neq[1]);	
  NNEW(i_flag,ITG, neq[1]);
  NNEW(a_flag,ITG, neq[1]);
  
  /* Fill l_flag*/       
  //Slave
  for (i=0;i<*ntie;i++){
    if(tieset[i*(81*3)+80]=='C'){      
      for(j=nslavnode[i];j<nslavnode[i+1];j++){         
	nodesf=islavnode[j];
	// right now temperature degrees of freedom are set to N, since l_flag[k-1]==0
	for (l=0;l<3;l++){ //test of dof	    
	  k = nactdof[mt*(nodesf)-3+l];	                              
	  if(k>0&&islavact[j]==0) l_flag[k-1]=2;
	  if(k>0&&islavact[j]>0) {l_flag[k-1]=3;iact++;}
	  // no-gap and no-LM nodes must be treated as master nodes
	  if(k>0&&(islavact[j]==-1 ||islavact[j]==-2)) l_flag[k-1]=1;
	  // nodes with no integration points can be treated N-nodes
	  if(k>0&&islavact[j]==-3)l_flag[k-1]=0;
	  //printf("slavedof %" ITGFORMAT " node %" ITGFORMAT " act %" ITGFORMAT " flag %" ITGFORMAT "\n",k,nodesf,islavact[j],l_flag[k-1]);
	}
	if (ithermal[0]==3){ //coupled thermo-mechanical calculation    
	  k = nactdof[mt*(nodesf)-4];	                              
	  if(k>0&&islavact[j]==0) l_flag[k-1]=2;
	  if(k>0&&islavact[j]>0) {l_flag[k-1]=3;iact++;}
	  // no-gap and no-LM nodes must be treated as master nodes
	  if(k>0&&(islavact[j]==-1 ||islavact[j]==-2)) l_flag[k-1]=1;
	  // nodes with no integration points can be treated N-nodes
	  if(k>0&&islavact[j]==-3)l_flag[k-1]=0;
	  //printf("slavedof %" ITGFORMAT " node %" ITGFORMAT " act %" ITGFORMAT " flag %" ITGFORMAT "\n",k,nodesf,islavact[j],l_flag[k-1]);
	}	
      }
    }
  }
  //Master
  for (i=0;i<*ntie;i++){
    if(tieset[i*(81*3)+80]=='C'){      
      for(j=nmastnode[i];j<nmastnode[i+1];j++){        
	nodem=imastnode[j];	 
	for (l=0;l<3;l++){ //test of dof	    
	  k = nactdof[mt*(nodem)-3+l];            
	  if(k>0){if(jqbdtil[k]-jqbdtil[k-1]>0){ l_flag[k-1]=1;}else{ l_flag[k-1]=0;}}       
	}
	if (ithermal[0]==3){ //coupled thermo-mechanical calculation 
	  k = nactdof[mt*(nodem)-4];            
	  if(k>0){if(jqbdtil[k]-jqbdtil[k-1]>0){ l_flag[k-1]=1;}else{ l_flag[k-1]=0;}}	
	}
      }
    }
  }
  /*** Fill the local row ***/        
  row_ln=0;        
  row_lm=0;        
  row_la=0;
  row_li=0;
  /* Stock of the diagonale */		
  NNEW(ad_nn,double,neq[1]);	
  NNEW(ad_mm,double,neq[1]);
  NNEW(irow_mmd,ITG, neq[1]);	
  NNEW(jq_mmd,ITG, neq[1]);
  NNEW(ad_aa,double,neq[1]);	
  NNEW(irow_aad,ITG, neq[1]);	
  NNEW(jq_aad,ITG, neq[1]);
  NNEW(ad_ii,double,neq[1]);	
  NNEW(irow_iid,ITG, neq[1]);	
  NNEW(jq_iid,ITG, neq[1]);
  
  /**** For the construction of Bhat, the new rhs***/
  
  NNEW(f_a,double,neq[1]);
  NNEW(f_i,double,neq[1]);
  NNEW(f_m,double,neq[1]);	        
  for (i=0;i<neq[1];i++){ 
    bhat[i]=b[i];
    switch(l_flag[i]){           
    case 0 : 	{//N			                       
      ad_nn[row_ln]=ad[i];                        			
      n_flag[i]=(++row_ln);                        
      break;
    }           
    case 1 : 	{//M                         
      f_m[row_lm]=b[i];			
      m_flag[i]=(++row_lm);				
      ad_mm[row_lm-1]=ad[i];			
      jq_mmd[row_lm-1]=row_lm;			
      irow_mmd[row_lm-1]=row_lm;                        
      break;
    }    
    case 2 : 	{//I			                       
      ad_ii[row_li]=ad[i];                        
      f_i[row_li]=b[i];
      i_flag[i]=(++row_li); 
      jq_iid[row_li-1]=row_li;			
      irow_iid[row_li-1]=row_li;
      break;
    }     
    case 3 : 	{//A	
      f_a[row_la]=b[i];			
      a_flag[i]=(++row_la);                        
      ad_aa[row_la-1]=ad[i];			
      jq_aad[row_la-1]=row_la;			
      irow_aad[row_la-1]=row_la;			                       
      break;
    }          
    default :    {
      printf("error l_flag\n");                        
      break;
    } 	         
    }		
  }
  jq_iid[row_li]=row_li+1;
  jq_aad[row_la]=row_la+1;	
  jq_mmd[row_lm]=row_lm+1;		
  RENEW(f_a,double,row_la);
  RENEW(f_i,double,row_li);
  RENEW(f_m,double,row_lm);	
  RENEW(ad_nn,double,row_ln);	
  RENEW(ad_mm,double,row_lm);	
  RENEW(irow_mmd,ITG, row_lm);	
  RENEW(jq_mmd,ITG, row_lm+1);	
  RENEW(ad_aa,double,row_la);	
  RENEW(irow_aad,ITG, row_la);	
  RENEW(jq_aad,ITG, row_la+1);
  RENEW(ad_ii,double,row_li);	
  RENEW(irow_iid,ITG, row_li);	
  RENEW(jq_iid,ITG, row_li+1);  
  
  /* extract submatrices **/
  
  NNEW(au_nn,double,*nzs/8);        
  NNEW(au_nm,double,*nzs/8);
  NNEW(au_ni,double,*nzs/8);
  NNEW(au_na,double,*nzs/8);
  NNEW(au_mn,double,*nzs/8);                
  NNEW(au_mm,double,*nzs/8);  
  NNEW(au_mi,double,*nzs/8);
  NNEW(au_ma,double,*nzs/8);  
  NNEW(au_in,double,*nzs/8);
  NNEW(au_im,double,*nzs/8);
  NNEW(au_ii,double,*nzs/8);
  NNEW(au_ia,double,*nzs/8);
  NNEW(au_an,double,*nzs/8);
  NNEW(au_am,double,*nzs/8);  
  NNEW(au_ai,double,*nzs/8);
  NNEW(au_aa,double,*nzs/8);
  
  NNEW(irow_nm,ITG, *nzs/8); 
  NNEW(irow_nn,ITG, *nzs/8);
  NNEW(irow_ni,ITG, *nzs/8);
  NNEW(irow_na,ITG, *nzs/8);  
  NNEW(irow_mm,ITG, *nzs/8); 
  NNEW(irow_mn,ITG, *nzs/8);
  NNEW(irow_mi,ITG, *nzs/8);
  NNEW(irow_ma,ITG, *nzs/8); 
  NNEW(irow_im,ITG, *nzs/8); 
  NNEW(irow_in,ITG, *nzs/8);
  NNEW(irow_ii,ITG, *nzs/8);
  NNEW(irow_ia,ITG, *nzs/8);  
  NNEW(irow_an,ITG, *nzs/8);               
  NNEW(irow_am,ITG, *nzs/8);        
  NNEW(irow_ai,ITG, *nzs/8); 
  NNEW(irow_aa,ITG, *nzs/8);
  
  NNEW(jq_nn,ITG, row_ln+1);      
  NNEW(jq_nm,ITG, row_lm+1); 
  NNEW(jq_ni,ITG, row_li+1); 
  NNEW(jq_na,ITG, row_la+1);   
  NNEW(jq_mn,ITG, row_ln+1);      
  NNEW(jq_mm,ITG, row_lm+1); 
  NNEW(jq_mi,ITG, row_li+1); 
  NNEW(jq_ma,ITG, row_la+1);   
  NNEW(jq_in,ITG, row_ln+1);      
  NNEW(jq_im,ITG, row_lm+1); 
  NNEW(jq_ii,ITG, row_li+1); 
  NNEW(jq_ia,ITG, row_la+1);    
  NNEW(jq_an,ITG, row_ln+1);               
  NNEW(jq_am,ITG, row_lm+1); 
  NNEW(jq_ai,ITG, row_li+1);
  NNEW(jq_aa,ITG, row_la+1); 
  
  jq_nn[0]=1;  
  jq_nm[0]=1; 
  jq_ni[0]=1; 
  jq_na[0]=1;
  jq_mn[0]=1;      
  jq_mm[0]=1;   
  jq_mi[0]=1;  
  jq_ma[0]=1;  
  jq_in[0]=1;  
  jq_im[0]=1;
  jq_ii[0]=1;
  jq_ia[0]=1;
  jq_an[0]=1;               
  jq_am[0]=1;   
  jq_ai[0]=1;  
  jq_aa[0]=1;   
  
  nzs_nn=*nzs/8.0;
  nzs_nm=*nzs/8.0;
  nzs_ni=*nzs/8.0;
  nzs_na=*nzs/8.0;
  nzs_mn=*nzs/8.0;
  nzs_mm=*nzs/8.0;
  nzs_mi=*nzs/8.0;
  nzs_ma=*nzs/8.0;
  nzs_in=*nzs/8.0;
  nzs_im=*nzs/8.0;
  nzs_ii=*nzs/8.0;
  nzs_ia=*nzs/8.0;
  nzs_an=*nzs/8.0;
  nzs_am=*nzs/8.0;
  nzs_ai=*nzs/8.0;
  nzs_aa=*nzs/8.0;
  
  ifree_nn=1;
  ifree_nm=1;
  ifree_ni=1;
  ifree_na=1;
  ifree_mn=1;
  ifree_mm=1;
  ifree_mi=1;
  ifree_ma=1;
  ifree_in=1;
  ifree_im=1;
  ifree_ii=1;
  ifree_ia=1;
  ifree_an=1;
  ifree_am=1;
  ifree_ai=1;
  ifree_aa=1;  
  
  for(j=0;j<neq[1];j++){	   
    m=j+1;
    for(i=jqc[j]-1;i<jqc[j+1]-1;i++){             
      switch(l_flag[j]){               	
      case 0 : {// matrices XN 
	m=n_flag[j]+1;
	switch(l_flag[irowc[i]-1]){                          
	case 0 : {//NN                                   	
	//printf("multimortar2 101 %e\n",auc[i]);
	  insertas_ws(&irow_nn,&(n_flag[irowc[i]-1]),&m,&ifree_nn,
		      &nzs_nn,&auc[i],&au_nn);	  
	  break;
	}                          
	case 1 : //MN                                    	 
	//printf("multimortar2 102 %e\n",auc[i]);
	  insertas_ws(&irow_mn,&(m_flag[irowc[i]-1]),&m,&ifree_mn,
		      &nzs_mn,&auc[i],&au_mn);
	  break;                           
	case 2 : //IN                                    	  
	//printf("multimortar2 103 %e\n",auc[i]);
	  insertas_ws(&irow_in,&(i_flag[irowc[i]-1]),&m,&ifree_in,
		      &nzs_in,&auc[i],&au_in);
	  break;
	  
	case 3: //AN
	//printf("multimortar2 104 %e\n",auc[i]);
	  insertas_ws(&irow_an,&(a_flag[irowc[i]-1]),&m,&ifree_an,
		      &nzs_an,&auc[i],&au_an);
	  break;
	default : 	
	  break;                        
	}                       
	break;
      }	    
      case 1 : {// matrices XM  
	m=m_flag[j]+1;
	switch(l_flag[irowc[i]-1]){ 
        case 0 : {//NM                                   	
	//printf("multimortar2 105 %e\n",auc[i]);
	  insertas_ws(&irow_nm,&(n_flag[irowc[i]-1]),&m,&ifree_nm,
		      &nzs_nm,&auc[i],&au_nm);
	  break;
	}   
	case 1 : {//MM                                   	 
	//printf("multimortar2 106 %e\n",auc[i]);
	  insertas_ws(&irow_mm,&(m_flag[irowc[i]-1]),&m,&ifree_mm,
		      &nzs_mm,&auc[i],&au_mm);
	  break;
	}                           
	case 2 : //IM                                    	 
	//printf("multimortar2 107 %e\n",auc[i]);
	  insertas_ws(&irow_im,&(i_flag[irowc[i]-1]),&m,&ifree_im,
		      &nzs_im,&auc[i],&au_im);
	  break;
	case 3 : {//AM                                   	
	//printf("multimortar2 108 %e\n",auc[i]);
	  insertas_ws(&irow_am,&(a_flag[irowc[i]-1]),&m,&ifree_am,
		      &nzs_am,&auc[i],&au_am);
	  break;
	} 	  
	default : 	
	  break;	  
	}
	break;
      }
      case 2 : {// matrices XI 
	m=i_flag[j]+1;
	switch(l_flag[irowc[i]-1]){                           
        case 0 : {//NI                                   	 
	//printf("multimortar2 109 %e\n",auc[i]);
	  insertas_ws(&irow_ni,&(n_flag[irowc[i]-1]),&m,&ifree_ni,
		      &nzs_ni,&auc[i],&au_ni);
	  break;
	}   
	case 1 : {//MI                                   	
	//printf("multimortar2 110 %e\n",auc[i]);
	  insertas_ws(&irow_mi,&(m_flag[irowc[i]-1]),&m,&ifree_mi,
		      &nzs_mi,&auc[i],&au_mi);
	  break;
	}                           
	case 2 : //II                                   	
	//printf("multimortar2 111 %e\n",auc[i]);
	  insertas_ws(&irow_ii,&(i_flag[irowc[i]-1]),&m,&ifree_ii,
		      &nzs_ii,&auc[i],&au_ii);
	  break;
	case 3 : {//AI                                   	
	//printf("multimortar2 112 %e\n",auc[i]);
	  insertas_ws(&irow_ai,&(a_flag[irowc[i]-1]),&m,&ifree_ai,
		      &nzs_ai,&auc[i],&au_ai);
	  break;
	} 	  
	default : 	
	  break;                     
	}
	break;
      }
      case 3 : {// matrices XA 
	m=a_flag[j]+1;
	switch(l_flag[irowc[i]-1]){                           
        case 0 : {//NA                                   	
	//printf("multimortar2 113 %e\n",auc[i]);
	  insertas_ws(&irow_na,&(n_flag[irowc[i]-1]),&m,&ifree_na,
		      &nzs_na,&auc[i],&au_na);
	  break;
	}   
	case 1 : {//MA                                   	
	//printf("multimortar2 114 %e\n",auc[i]);
	  insertas_ws(&irow_ma,&(m_flag[irowc[i]-1]),&m,&ifree_ma,
		      &nzs_ma,&auc[i],&au_ma);
	  break;
	}                           
	case 2 : //IA                                  	
	//printf("multimortar2 115 %e\n",auc[i]);
	  insertas_ws(&irow_ia,&(i_flag[irowc[i]-1]),&m,&ifree_ia,
		      &nzs_ia,&auc[i],&au_ia);
	  break;
	case 3 : {//AA                                   	
	//printf("multimortar2 116 %e\n",auc[i]);
	  insertas_ws(&irow_aa,&(a_flag[irowc[i]-1]),&m,&ifree_aa,
		      &nzs_aa,&auc[i],&au_aa);
	  break;
	} 	  
	default : 	
	  break;                     
	}
	break;
      }      
      default : 
	break;             
      } // end switch column lflag           
    } // end loop nzs column
    
    /* actualisation of the jq **/
    
    if (n_flag[j]!=0) jq_nn[n_flag[j]]=ifree_nn;		
    if (n_flag[j]!=0)jq_mn[n_flag[j]]=ifree_mn;
    if (n_flag[j]!=0)jq_in[n_flag[j]]=ifree_in;
    if (n_flag[j]!=0)jq_an[n_flag[j]]=ifree_an;    
    if (m_flag[j]!=0) jq_nm[m_flag[j]]=ifree_nm;		
    if (m_flag[j]!=0)jq_mm[m_flag[j]]=ifree_mm;
    if (m_flag[j]!=0)jq_im[m_flag[j]]=ifree_im;
    if (m_flag[j]!=0)jq_am[m_flag[j]]=ifree_am;   
    if (i_flag[j]!=0) jq_ni[i_flag[j]]=ifree_ni;		
    if (i_flag[j]!=0)jq_mi[i_flag[j]]=ifree_mi;
    if (i_flag[j]!=0)jq_ii[i_flag[j]]=ifree_ii;
    if (i_flag[j]!=0)jq_ai[i_flag[j]]=ifree_ai;    
    if (a_flag[j]!=0) jq_na[a_flag[j]]=ifree_na;		
    if (a_flag[j]!=0)jq_ma[a_flag[j]]=ifree_ma;
    if (a_flag[j]!=0)jq_ia[a_flag[j]]=ifree_ia;
    if (a_flag[j]!=0)jq_aa[a_flag[j]]=ifree_aa;   
  } // end loop over the global column 
  /*****************************/        
  
  
  /****************************************************************/   
  printf("\tmm2: N %" ITGFORMAT " M %" ITGFORMAT " I %" ITGFORMAT " A %" ITGFORMAT "  nzs %" ITGFORMAT " \n",row_ln,row_lm,row_li,row_la,*nzs);  
  nstrue=3*nslavs;   
  nmtrue=3*nmasts;
  nntrue=neq[1]-nstrue-nmtrue; 	   
  
  if ((ifree_nn-1)!=0){ 
    nzs_nn=ifree_nn-1;
    RENEW(au_nn,double,nzs_nn);           
    RENEW(irow_nn,ITG, nzs_nn);
  }else{
    nzs_nn=ifree_nn-1;
    //printf("A_NN null\n");           
    RENEW(au_nn,double,1);           
    RENEW(irow_nn,ITG, 1);		
  }
  if ((ifree_nm-1)!=0){ 
    nzs_nm=ifree_nm-1;
    RENEW(au_nm,double,nzs_nm);           
    RENEW(irow_nm,ITG, nzs_nm);    
  }else{
    nzs_nm=ifree_nm-1;
    //printf("A_NM null\n");           
    RENEW(au_nm,double,1);           
    RENEW(irow_nm,ITG, 1);		
  }
  if ((ifree_ni-1)!=0){ 
    nzs_ni=ifree_ni-1;
    RENEW(au_ni,double,nzs_ni);           
    RENEW(irow_ni,ITG, nzs_ni);   
  }else{
    nzs_ni=ifree_ni-1;
    //printf("A_NI null\n");           
    RENEW(au_ni,double,1);           
    RENEW(irow_ni,ITG, 1);		
  }
  if ((ifree_na-1)!=0){
    nzs_na=ifree_na-1;
    RENEW(au_na,double,nzs_na);           
    RENEW(irow_na,ITG, nzs_na);       	   
  }else{
    nzs_na=ifree_na-1;
    //printf("A_NA null\n");           
    RENEW(au_na,double,1);           
    RENEW(irow_na,ITG, 1);		
  }  
  if ((ifree_mn-1)!=0){
    nzs_mn=ifree_mn-1;
    RENEW(au_mn,double,nzs_mn);           
    RENEW(irow_mn,ITG, nzs_mn);       	   
  }else{
    nzs_mn=ifree_mn-1;
    //printf("A_MN null\n");           
    RENEW(au_mn,double,1);           
    RENEW(irow_mn,ITG, 1);		
  }
  if ((ifree_mm-1)!=0){
    nzs_mm=ifree_mm-1;
    RENEW(au_mm,double,nzs_mm);           
    RENEW(irow_mm,ITG, nzs_mm);       	   
  }else{
    nzs_mm=ifree_mm-1;
    //printf("A_mM null\n");           
    RENEW(au_mm,double,1);           
    RENEW(irow_mm,ITG, 1);		
  }
  if ((ifree_mi-1)!=0){
    nzs_mi=ifree_mi-1;
    RENEW(au_mi,double,nzs_mi);           
    RENEW(irow_mi,ITG, nzs_mi);       	   
  }else{
    nzs_mi=ifree_mi-1;
    //printf("A_MI null\n");           
    RENEW(au_mi,double,1);           
    RENEW(irow_mi,ITG, 1);		
  }
  if ((ifree_ma-1)!=0){
    nzs_ma=ifree_ma-1;
    RENEW(au_ma,double,nzs_ma);           
    RENEW(irow_ma,ITG, nzs_ma);       	   
  }else{
    nzs_ma=ifree_ma-1;
    //printf("A_MA null\n");           
    RENEW(au_ma,double,1);           
    RENEW(irow_ma,ITG, 1);		
  } 
  if ((ifree_in-1)!=0){ 
    nzs_in=ifree_in-1;
    RENEW(au_in,double,nzs_in);           
    RENEW(irow_in,ITG, nzs_in);       	   
  }else{	
    nzs_in=ifree_in-1;
    //printf("A_IN null\n");           
    RENEW(au_in,double,1);           
    RENEW(irow_in,ITG, 1);		
  }
  if ((ifree_im-1)!=0){ 
    nzs_im=ifree_im-1;
    RENEW(au_im,double,nzs_im);           
    RENEW(irow_im,ITG, nzs_im);       	   
  }else{
    nzs_im=ifree_im-1;
    //printf("A_IM null\n");           
    RENEW(au_im,double,1);           
    RENEW(irow_im,ITG, 1);		
  }
  if ((ifree_ii-1)!=0){
    nzs_ii=ifree_ii-1;
    RENEW(au_ii,double,nzs_ii);           
    RENEW(irow_ii,ITG, nzs_ii);       	   
  }else{	
    nzs_ii=ifree_ii-1;
    //printf("A_II null\n");           
    RENEW(au_ii,double,1);           
    RENEW(irow_ii,ITG, 1);		
  }
  if ((ifree_ia-1)!=0){
    nzs_ia=ifree_ia-1;
    RENEW(au_ia,double,nzs_ia);           
    RENEW(irow_ia,ITG, nzs_ia);       	   
  }else{
    nzs_ia=ifree_ia-1;
    //printf("A_IA null\n");           
    RENEW(au_ia,double,1);           
    RENEW(irow_ia,ITG, 1);		
  }  
  if ((ifree_an-1)!=0){
    nzs_an=ifree_an-1;
    RENEW(au_an,double,nzs_an);           
    RENEW(irow_an,ITG, nzs_an);       	   
  }else{
    nzs_an=ifree_an-1;
    //printf("A_AN null\n");           
    RENEW(au_an,double,1);           
    RENEW(irow_an,ITG, 1);		
  }
  if ((ifree_am-1)!=0){
    nzs_am=ifree_am-1;
    RENEW(au_am,double,nzs_am);           
    RENEW(irow_am,ITG, nzs_am);       	   
  }else{
    nzs_am=ifree_am-1;
    //printf("A_AM null\n");           
    RENEW(au_am,double,1);           
    RENEW(irow_am,ITG, 1);		
  }
  if ((ifree_ai-1)!=0){
    nzs_ai=ifree_ai-1;
    RENEW(au_ai,double,nzs_ai);           
    RENEW(irow_ai,ITG, nzs_ai);       	   
  }else{
    nzs_ai=ifree_ai-1;
    //printf("A_AI null\n");           
    RENEW(au_ai,double,1);           
    RENEW(irow_ai,ITG, 1);		
  }
  if ((ifree_aa-1)!=0){
    nzs_aa=ifree_aa-1;
    RENEW(au_aa,double,nzs_aa);           
    RENEW(irow_aa,ITG, nzs_aa);       	   
  }else{
    nzs_aa=ifree_aa-1;
    //printf("A_AA null\n");           
    RENEW(au_aa,double,1);           
    RENEW(irow_aa,ITG, 1);		
  }  
  
  printf("\tmm2: MN %" ITGFORMAT " %" ITGFORMAT " IN %" ITGFORMAT " %" ITGFORMAT " AN %" ITGFORMAT " %" ITGFORMAT " \n",nzs_mn,nzs_nm,nzs_in,nzs_ni,nzs_an,nzs_na);
  printf("\tmm2: IM %" ITGFORMAT " %" ITGFORMAT " AM %" ITGFORMAT " %" ITGFORMAT " AI %" ITGFORMAT " %" ITGFORMAT " \n",nzs_im,nzs_mi,nzs_am,nzs_ma,nzs_ai,nzs_ia);
  
  /* generate D_d, D_d^-1,B_a ,B_d^til in local dofs **/
  
  nzsbdtil2=jqbdtil[neq[1]]-1;
  NNEW(au_bd1,double,nzsbdtil2);	  
  NNEW(irow_bd1,ITG, nzsbdtil2);	  
  NNEW(jq_bd1,ITG, row_lm+1);	
  nzs_bd1=0;	  
  jq_bd1[0]=1;
  nzsbdtil2=jqbdtil2[neq[1]]-1;
  NNEW(au_bdtil2,double,nzsbdtil2);	  
  NNEW(irow_bdtil2,ITG, nzsbdtil2);	  
  NNEW(jq_bdtil2,ITG, row_lm+1);	
  nzs_bdtil2=0;	  
  jq_bdtil2[0]=1;
  impclack=0;
  nzsddtil2a=jqbdtil2[neq[1]]-1;
  NNEW(au_ddtil2a,double,nzsddtil2a);	  
  NNEW(irow_ddtil2a,ITG, nzsddtil2a);	  
  NNEW(jq_ddtil2a,ITG, row_la+1);	
  nzs_ddtil2a=0;	  
  jq_ddtil2a[0]=1;
  nzsddtil2i=jqbdtil2[neq[1]]-1;
  NNEW(au_ddtil2i,double,nzsddtil2i);	  
  NNEW(irow_ddtil2i,ITG, nzsddtil2i);	  
  NNEW(jq_ddtil2i,ITG, row_li+1);	
  nzs_ddtil2i=0;	  
  jq_ddtil2i[0]=1;
  for(j=0;j<neq[1];j++){	        	    
    for (i=jqbdtil2[j]-1;i<jqbdtil2[j+1]-1;i++){	      
      switch(l_flag[j]){		
      case 1 : //Matrix XM			
	switch(l_flag[irowbdtil2[i]-1]){			  
	case 3 : //AM here Bdtild				    
	  irow_bdtil2[nzs_bdtil2]=a_flag[irowbdtil2[i]-1];
	  au_bdtil2[nzs_bdtil2]=aubdtil2[i]; // Needed for slave part
	  nzs_bdtil2++;
	  break;
	default : 
	  break;				    			
	}
	break;
      case 0 :
	impclack++;
	break;	
      case 3 :
	//impclack++;
	break;
      default :
	break;	      
      }	    
    }	       
    if (m_flag[j]!=0) jq_bdtil2[m_flag[j]]=nzs_bdtil2+1;
  }
  impclack=0;
  for(j=0;j<neq[1];j++){	        	    
    for (i=jqbdtil2[j]-1;i<jqbdtil2[j+1]-1;i++){	      
      switch(l_flag[j]){		
      case 3 : //Matrix XA			
	switch(l_flag[irowbdtil2[i]-1]){			  
	case 3 : //AA here Ddtild				    
	  irow_ddtil2a[nzs_ddtil2a]=a_flag[irowbdtil2[i]-1];
	  au_ddtil2a[nzs_ddtil2a]=aubdtil2[i]; // Needed for slave part
	  nzs_ddtil2a++;
	  break;
	default : 
	  break;				    			
	}
	break;
      case 0 :
	impclack++;
	break;	
      default :
	break;	      
      }	    
    }	       
    if (a_flag[j]!=0) jq_ddtil2a[a_flag[j]]=nzs_ddtil2a+1;
  }
  impclack=0;
  for(j=0;j<neq[1];j++){	        	    
    for (i=jqbdtil2[j]-1;i<jqbdtil2[j+1]-1;i++){	      
      switch(l_flag[j]){		
      case 2 : //Matrix XI			
	switch(l_flag[irowbdtil2[i]-1]){			  
	case 3 : //AI here Ddtild				    
	  irow_ddtil2i[nzs_ddtil2i]=a_flag[irowbdtil2[i]-1];
	  au_ddtil2i[nzs_ddtil2i]=aubdtil2[i]; // Needed for slave part
	  nzs_ddtil2i++;
	  break;
	default : 
	  break;				    			
	}
	break;
      case 0 :
	impclack++;
	break;	
      default :
	break;	      
      }	    
    }	       
    if (i_flag[j]!=0) jq_ddtil2i[i_flag[j]]=nzs_ddtil2i+1;
  }
  for(j=0;j<neq[1];j++){	        	    
    for (i=jqbdtil[j]-1;i<jqbdtil[j+1]-1;i++){	      
      switch(l_flag[j]){		
      case 1 : //Matrix XM			
	switch(l_flag[irowbdtil[i]-1]){			  
	case 3 : //AM here Bdtild				    
	  irow_bd1[nzs_bd1]=a_flag[irowbdtil[i]-1];
	  au_bd1[nzs_bd1]=-aubdtil[i]; // Needed for master part
	  nzs_bd1++;
	  break;			  
	default : 
	  break;				    			
	}
	break;
      default : 	       
	break;	      
      }	    
    }	    
    if (m_flag[j]!=0) jq_bd1[m_flag[j]]=nzs_bd1+1;
  }
  RENEW(irow_bd1,ITG, nzs_bd1);	 
  RENEW(au_bd1,double,nzs_bd1);
  RENEW(irow_bdtil2,ITG, nzs_bdtil2);	 
  RENEW(au_bdtil2,double,nzs_bdtil2);
  RENEW(irow_ddtil2i,ITG, nzs_ddtil2i);	 
  RENEW(au_ddtil2i,double,nzs_ddtil2i);
  RENEW(irow_ddtil2a,ITG, nzs_ddtil2a);	 
  RENEW(au_ddtil2a,double,nzs_ddtil2a);
  printf("\tmm2: nzs_ddtil2a %" ITGFORMAT " mpclack %" ITGFORMAT " \n",nzs_ddtil2a,impclack);
  printf("\tmm2: nzs_ddtil2i %" ITGFORMAT " mpclack %" ITGFORMAT " \n",nzs_ddtil2i,impclack);
  printf("\tmm2: nzs_bdtil %" ITGFORMAT " mpclack %" ITGFORMAT " \n",nzs_bdtil2,impclack);
  printf("\tmm2: nzs_bd1 %" ITGFORMAT "\n",nzs_bd1);
  /*for (i=0;i<row_lm;i++){    
    if(jq_bdtil2[i+1]-jq_bdtil2[i]>0){   
    numb=jq_bdtil2[i+1]-jq_bdtil2[i];    
    FORTRAN(isortid,(&irow_bdtil2[jq_bdtil2[i]-1],&au_bdtil2[jq_bdtil2[i]-1],&numb,&kflag));    
    }  
    }*/



    /* use Ddtil2=Id instead of Ddtil^-1 **/
    /* this is a adaption of an older version, were the matrix and 
      the right hand side was multiplied by auddinv before embeding
      the contact conditions. don't get confused by names 
      au_diiinv,au_diainv,au_daiinv,au_daainv */
  /*nzsddtil2=jqddtil2[neq[1]]-1;
  NNEW(au_diainv,double,nzsddtil2);	  
  NNEW(irow_diainv,ITG, nzsddtil2);	  
  NNEW(jq_diainv,ITG, row_la+1);
  NNEW(au_daainv,double,nzsddtil2);	  
  NNEW(irow_daainv,ITG, nzsddtil2);	  
  NNEW(jq_daainv,ITG, row_la+1);
  nzs_diainv=0;	  
  jq_diainv[0]=1;
  nzs_daainv=0;	  
  jq_daainv[0]=1;
  NNEW(au_daiinv,double,nzsddtil2);	  
  NNEW(irow_daiinv,ITG, nzsddtil2);	  
  NNEW(jq_daiinv,ITG, row_li+1);
  NNEW(au_diiinv,double,nzsddtil2);	  
  NNEW(irow_diiinv,ITG, nzsddtil2);	  
  NNEW(jq_diiinv,ITG, row_li+1);
  nzs_daiinv=0;	  
  jq_daiinv[0]=1;
  nzs_diiinv=0;	  
  jq_diiinv[0]=1;
  for(j=0;j<neq[1];j++){	        	    
    for (i=jqddtil2[j]-1;i<jqddtil2[j+1]-1;i++){	      
      switch(l_flag[j]){
      case 2 : //Matrix XI			
	switch(l_flag[irowddtil2[i]-1]){
	case 2 ://II here D_iiinv
	  irow_diiinv[nzs_diiinv]=i_flag[irowddtil2[i]-1];				    
	  au_diiinv[nzs_diiinv]=auddtil2[i];				    
	  nzs_diiinv++;	
	  break;
	case 3 : //AI here D_aiinv
	  // NOW EMPTY
	  irow_daiinv[nzs_daiinv]=a_flag[irowddtil2[i]-1];				    
	  au_daiinv[nzs_daiinv]=auddtil2[i];				    
	  nzs_daiinv++;				    
	  break;			  
	default : 
	  break;				    			
	}
      break;
      
      case 3 : //Matrix XA			
	switch(l_flag[irowddtil2[i]-1]){
	case 2 ://IA here D_iainv
	  // NOW EMPTY
	  irow_diainv[nzs_diainv]=i_flag[irowddtil2[i]-1];				    
	  au_diainv[nzs_diainv]=auddtil2[i];				    
	  nzs_diainv++;	
	  break;
	case 3 : //AA here D_aainv				    
	  irow_daainv[nzs_daainv]=a_flag[irowddtil2[i]-1];				    
	  au_daainv[nzs_daainv]=auddtil2[i];				    
	  nzs_daainv++;				    
	  break;			  
	default : 
	  break;				    			
	}			  
      default : 	       
	break;	      
      }	    
    }	        
    if (i_flag[j]!=0) jq_daiinv[i_flag[j]]=nzs_daiinv+1;
    if (i_flag[j]!=0) jq_diiinv[i_flag[j]]=nzs_diiinv+1;
    if (a_flag[j]!=0) jq_daainv[a_flag[j]]=nzs_daainv+1;
    if (a_flag[j]!=0) jq_diainv[a_flag[j]]=nzs_diainv+1;
  }
  
  RENEW(irow_diainv,ITG, nzs_diainv);	 
  RENEW(au_diainv,double,nzs_diainv);
  RENEW(irow_diiinv,ITG, nzs_diiinv);	 
  RENEW(au_diiinv,double,nzs_diiinv);
  RENEW(irow_daainv,ITG, nzs_daainv);	 
  RENEW(au_daainv,double,nzs_daainv);
  RENEW(irow_daiinv,ITG, nzs_daiinv);	 
  RENEW(au_daiinv,double,nzs_daiinv); 
  printf("\tmm2: dia %" ITGFORMAT " daa %" ITGFORMAT " dai %" ITGFORMAT " dii %" ITGFORMAT " \n",nzs_diainv,nzs_daainv,nzs_daiinv,nzs_diiinv);
  */
  /* add diagonals in K **/ 
  printf("\tmm2: add diagonals\n");
  nzs_mmf=nzs_mm+row_lm;
  NNEW(jq_mmf,ITG, row_lm+1);	
  NNEW(au_mmf,double,nzs_mmf);	
  NNEW(irow_mmf,ITG, nzs_mmf);
  add_rect(ad_mm,irow_mmd,jq_mmd,row_lm,row_lm,
	   au_mm,irow_mm,jq_mm,row_lm,row_lm,
	   &au_mmf,&irow_mmf,jq_mmf,&nzs_mmf);
  
  nzs_iif=nzs_ii+row_li;
  NNEW(jq_iif,ITG, row_li+1);	
  NNEW(au_iif,double,nzs_iif);	
  NNEW(irow_iif,ITG, nzs_iif);	   
  add_rect(ad_ii,irow_iid,jq_iid,row_li,row_li,
	   au_ii,irow_ii,jq_ii,row_li,row_li,
	   &au_iif,&irow_iif,jq_iif,&nzs_iif);
  
  nzs_aaf=nzs_aa+row_la;
  NNEW(jq_aaf,ITG, row_la+1);	
  NNEW(au_aaf,double,nzs_aaf);	
  NNEW(irow_aaf,ITG, nzs_aaf);	   
  add_rect(ad_aa,irow_aad,jq_aad,row_la,row_la,
	   au_aa,irow_aa,jq_aa,row_la,row_la,
	   &au_aaf,&irow_aaf,jq_aaf,&nzs_aaf);
  printf("\tmm2: NN %" ITGFORMAT " MM %" ITGFORMAT " II %" ITGFORMAT " AA %" ITGFORMAT "\n",nzs_nn,nzs_mmf,nzs_iif,nzs_aaf);
  SFREE(au_ii);SFREE(irow_ii);SFREE(jq_ii);
  SFREE(au_aa);SFREE(irow_aa);SFREE(jq_aa);
  SFREE(au_mm);SFREE(irow_mm);SFREE(jq_mm);
  
  
  
  /* Local => Global topology **/	  
  NNEW(n_flagr,ITG, row_ln);	
  NNEW(m_flagr,ITG, row_lm);	
  NNEW(i_flagr,ITG, row_li);
  NNEW(a_flagr,ITG, row_la);  
  for(j=0;j<neq[1];j++){		
    switch(l_flag[j]){			
    case 0 : {n_flagr[n_flag[j]-1]=j+1;break;} //N			
    case 1 : {m_flagr[m_flag[j]-1]=j+1;break;} //M			
    case 2 : {i_flagr[i_flag[j]-1]=j+1;break;} //I
    case 3 : {a_flagr[a_flag[j]-1]=j+1;break;} //A
    default : break;		
    }	
  }    
  
  /* K_MX -> K_MX^til: K_MX^til=K_MX-Btil^T*K_AX **/
  printf("\tmm2: alter K_MX\n");
  /* K_MN^til **/
  nzs_intmn=nzs_mn+1;
  NNEW(jq_intmn,ITG, row_ln+1);        
  NNEW(au_intmn,double,nzs_intmn);	
  NNEW(irow_intmn,ITG, nzs_intmn);        	
  multi_rect(au_bd1,irow_bd1,jq_bd1,row_la,row_lm,
	     au_an,irow_an,jq_an,row_la,row_ln,
	     &au_intmn,&irow_intmn,jq_intmn, &nzs_intmn);
  
  nzs_mntil=nzs_mn+1;	
  NNEW(jq_mntil,ITG, row_ln+1);	
  NNEW(au_mntil,double,nzs_mntil);	
  NNEW(irow_mntil,ITG, nzs_mntil);		
  add_rect(au_intmn,irow_intmn,jq_intmn,row_lm,row_ln,
	   au_mn,irow_mn,jq_mn,row_lm,row_ln,
	   &au_mntil,&irow_mntil,jq_mntil,&nzs_mntil);
  
  
  /* K_MM^til **/
  nzs_intmm=nzs_mmf;
  NNEW(jq_intmm,ITG, row_lm+1);        
  NNEW(au_intmm,double,nzs_intmm);	
  NNEW(irow_intmm,ITG, nzs_intmm);        	
  multi_rect(au_bd1,irow_bd1,jq_bd1,row_la,row_lm,
	     au_am,irow_am,jq_am,row_la,row_lm,
	     &au_intmm,&irow_intmm,jq_intmm, &nzs_intmm);
  
  nzs_mmtil=nzs_mmf;
  NNEW(jq_mmtil,ITG, row_lm+1);	
  NNEW(au_mmtil,double,nzs_mmtil);	
  NNEW(irow_mmtil,ITG, nzs_mmtil);		
  add_rect(au_intmm,irow_intmm,jq_intmm,row_lm,row_lm,
	   au_mmf,irow_mmf,jq_mmf,row_lm,row_lm,
	   &au_mmtil,&irow_mmtil,jq_mmtil,&nzs_mmtil);
  
  /* K_MI^til **/
  nzs_intmi=max(nzs_mi,nzs_bdtil2);
  NNEW(jq_intmi,ITG, row_li+1);        
  NNEW(au_intmi,double,nzs_intmi);	
  NNEW(irow_intmi,ITG, nzs_intmi);        	
  multi_rect(au_bd1,irow_bd1,jq_bd1,row_la,row_lm,
	     au_ai,irow_ai,jq_ai,row_la,row_li,
	     &au_intmi,&irow_intmi,jq_intmi, &nzs_intmi);
  
  nzs_mitil=max(nzs_mi,nzs_intmi);	     
  NNEW(jq_mitil,ITG, row_li+1);	
  NNEW(au_mitil,double,nzs_mitil);	
  NNEW(irow_mitil,ITG, nzs_mitil);		
  add_rect(au_intmi,irow_intmi,jq_intmi,row_lm,row_li,
	   au_mi,irow_mi,jq_mi,row_lm,row_li,
	   &au_mitil,&irow_mitil,jq_mitil,&nzs_mitil);
  
  /* K_MA^til **/
  
  nzs_intma=max(nzs_ma,nzs_bdtil2);
  NNEW(jq_intma,ITG, row_la+1);        
  NNEW(au_intma,double,nzs_intma);	
  NNEW(irow_intma,ITG, nzs_intma);        	
  multi_rect(au_bd1,irow_bd1,jq_bd1,row_la,row_lm,
	     au_aaf,irow_aaf,jq_aaf,row_la,row_la,
	     &au_intma,&irow_intma,jq_intma, &nzs_intma);
  
  nzs_matil=max(nzs_ma,nzs_intma);
  NNEW(jq_matil,ITG, row_la+1);	
  NNEW(au_matil,double,nzs_matil);	
  NNEW(irow_matil,ITG, nzs_matil);		
  add_rect(au_intma,irow_intma,jq_intma,row_lm,row_la,
	   au_ma,irow_ma,jq_ma,row_lm,row_la,
	   &au_matil,&irow_matil,jq_matil,&nzs_matil); 
  printf("\tmm2: MNtil %" ITGFORMAT " %" ITGFORMAT " MMtil %" ITGFORMAT " %" ITGFORMAT " MItil %" ITGFORMAT " %" ITGFORMAT " MAtil %" ITGFORMAT " %" ITGFORMAT " \n",nzs_intmn,nzs_mntil,nzs_intmm,nzs_mmtil,nzs_intmi,nzs_mitil,nzs_intma,nzs_matil);	   
  SFREE(au_mmf);SFREE(irow_mmf);SFREE(jq_mmf);
  SFREE(au_mn);SFREE(irow_mn);SFREE(jq_mn);
  SFREE(au_mi);SFREE(irow_mi);SFREE(jq_mi);
  SFREE(au_ma);SFREE(irow_ma);SFREE(jq_ma);
  
  
  /* call trafoNT **/
  NNEW(jq_antil,ITG, row_ln+1);
  NNEW(au_antil,double,jq_an[row_ln]-1);
  NNEW(irow_antil,ITG, jq_an[row_ln]-1);
  NNEW(jq_amtil,ITG, row_lm+1);
  NNEW(au_amtil,double,jq_am[row_lm]-1);
  NNEW(irow_amtil,ITG, jq_am[row_lm]-1);
  NNEW(jq_aitil,ITG, row_li+1);
  NNEW(au_aitil,double,jq_ai[row_li]-1);
  NNEW(irow_aitil,ITG, jq_ai[row_li]-1);
  NNEW(jq_aatil,ITG, row_la+1);
  NNEW(au_aatil,double,jq_aaf[row_la]-1);
  NNEW(irow_aatil,ITG, jq_aaf[row_la]-1);
  NNEW(f_atil,double,row_la);
  trafontmortar2(neq,nzs,islavactdof,islavact,nslavnode,nmastnode,
		 f_a,f_atil,
		 au_an, irow_an, jq_an,
		 au_am, irow_am, jq_am,
		 au_ai,irow_ai, jq_ai,
		 au_aaf, irow_aaf, jq_aaf,
		 &au_antil, &irow_antil, jq_antil,
		 &au_amtil, &irow_amtil, jq_amtil,
		 &au_aitil, &irow_aitil, jq_aitil,
		 &au_aatil, &irow_aatil, jq_aatil,
		 gap,
		 Bd,irowb,jqb,
		 Dd,irowd,jqd,
		 Ddtil,irowdtil,jqdtil,
		 au_bdtil2,irow_bdtil2,jq_bdtil2,
		 au_ddtil2i,irow_ddtil2i,jq_ddtil2i,
		 au_ddtil2a,irow_ddtil2a,jq_ddtil2a,
		 m_flagr,i_flagr,a_flagr,a_flag,i_flag,m_flag,
		 &row_ln,&row_lm,&row_li,&row_la,
		 slavnor,slavtan,
		 vold,vini,cstress,cstressini,
		 bp_old,nactdof,islavnode,imastnode,ntie,mi,nk,
		 nboun,ndirboun,nodeboun,xboun,
		 nmpc,ipompc,nodempc,coefmpc,
		 ikboun,ilboun,ikmpc,ilmpc,
		 nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
		 nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
		 tieset,
		 islavactdoftie,nelcon,elcon,tietol,ncmat_,ntmat_,
		 plicon,nplicon,npmat_,dtime,
		 irowtloc, jqtloc, autloc,
	         irowtlocinv, jqtlocinv, autlocinv,
		 islavnodeinv,lambdaiwan,lambdaiwanini,iit,nmethod,bet,
		 ithermal,plkcon,nplkcon
		 );
		 
  SFREE(au_an);SFREE(irow_an);SFREE(jq_an);
  SFREE(au_am);SFREE(irow_am);SFREE(jq_am);
  SFREE(au_ai);SFREE(irow_ai);SFREE(jq_ai);
  SFREE(au_aaf);SFREE(irow_aaf);SFREE(jq_aaf);	 
  
  
  /* Calculation of f_mhat = f_m + Bd_T*f_a **/	
  
  multi_rectv(au_bd1,irow_bd1,jq_bd1,row_la,row_lm,f_a,&v_r); //Bd_T*f_a        
  //printf("row_lm %" ITGFORMAT " \n",row_lm);	
  for (i=0;i<row_lm;i++) {b[m_flagr[i]-1]=f_m[i]+v_r[i];} // Just Master part is concerned	
  SFREE(v_r);  
  for (i=0;i<row_li;i++) {b[i_flagr[i]-1]=f_i[i];}
  for (i=0;i<row_la;i++) {b[a_flagr[i]-1]=f_atil[i];} 
  //SFREE(f_da);

  /* up to this point nzs is still the number of nonzero off-diagonals in the
     lower triangle of the symmetric matrix (= upper triangle transpose) */

  nzs_sym=*nzs;
  
  /* ASSEMBLE **/
  *nzs=*nzsc;	
  ifree=1;	
  
  RENEW(au,double,*nzs);

  /* CHANGE on 20190329: */
  
  for(j=nzs_sym;j<*nzs;j++){au[j]=0.;}

  /* END CHANGE */
  
  RENEW(irow,ITG,*nzs+neq[1]);

  /* CHANGE on 20190329: */
  
  for(j=nzs_sym;j<*nzs;j++){irow[j]=0;}

  /* END CHANGE */
  
  NNEW(mast1,ITG,*nzs);
    ITG loc_n,loc_i,loc_a,loc_m;	
  for(j=0;j<neq[1];j++){ad[j]=0.0;}
  for (j=0;j<neq[1];j++){		
    m=j+1;
    //printf("\tassemble dof %" ITGFORMAT " %" ITGFORMAT " ifree %" ITGFORMAT " \n",m,l_flag[j],ifree);
    switch(l_flag[j]){			
    case 0 : // insert Matrices NN, MN, IN, AN				    	
      loc_n=n_flag[j];				    	
      for(l=jq_nn[loc_n-1]-1;l<jq_nn[loc_n]-1;l++){ //NN
	//printf("multimortar2 1 %e\n",au_nn[l]);
	insertas(&irow,&mast1,&n_flagr[irow_nn[l]-1],&m,&ifree,nzs,&au_nn[l],&au);					 			  		
      }
      for(l=jq_mntil[loc_n-1]-1;l<jq_mntil[loc_n]-1;l++){ //MN					 
	//printf("multimortar2 2 %e\n",au_mntil[l]);
	insertas(&irow,&mast1,&m_flagr[irow_mntil[l]-1],&m,&ifree,nzs,&au_mntil[l],&au);					 				 			  		
      }		
      for(l=jq_in[loc_n-1]-1;l<jq_in[loc_n]-1;l++){ //IN					 
	//printf("multimortar2 3 %e\n",au_in[l]);
	insertas(&irow,&mast1,&i_flagr[irow_in[l]-1],&m,&ifree,nzs,&au_in[l],&au);													  		
      }
      for(l=jq_antil[loc_n-1]-1;l<jq_antil[loc_n]-1;l++){ //AN					 
	//printf("multimortar2 4 %e\n",au_antil[l]);
	insertas(&irow,&mast1,&a_flagr[irow_antil[l]-1],&m,&ifree,nzs,&au_antil[l],&au);													  		
      }
      ad[j]=ad_nn[loc_n-1];				    
      break;			
    case 1 : // insert Matrices NM, MM, IM, AM				    	
      loc_m=m_flag[j];	
      for(l=jq_nm[loc_m-1]-1;l<jq_nm[loc_m]-1;l++){ //NM					 
	//printf("multimortar2 5 %e\n",au_nm[l]);
	insertas(&irow,&mast1,&n_flagr[irow_nm[l]-1],&m,&ifree,nzs,&au_nm[l],&au);					 			  		
      }
      for(l=jq_mmtil[loc_m-1]-1;l<jq_mmtil[loc_m]-1;l++){ //MM	
	//printf("\tcdof %" ITGFORMAT " lcdof %" ITGFORMAT " rdof %" ITGFORMAT " lrdof %" ITGFORMAT " v %e slav %" ITGFORMAT " %" ITGFORMAT " \n",m,loc_m,m_flagr[irow_mmtil[l]-1],irow_mmtil[l],au_mmtil[l],islavactdof[j],islavactdof[m_flagr[irow_mmtil[l]-1]-1]);
        if ((m==m_flagr[irow_mmtil[l]-1])){ //diagonal of mm 						
	  ad[j]=au_mmtil[l];					
	}else{					 
	//printf("multimortar2 6 %e\n",au_mmtil[l]);
	  insertas(&irow,&mast1,&m_flagr[irow_mmtil[l]-1],&m,&ifree,nzs,&au_mmtil[l],&au);					  
	}			  		
      }		
      for(l=jq_im[loc_m-1]-1;l<jq_im[loc_m]-1;l++){ //IM					 
	//printf("multimortar2 7 %e\n",au_im[l]);
	insertas(&irow,&mast1,&i_flagr[irow_im[l]-1],&m,&ifree,nzs,&au_im[l],&au);													  		
      }
      for(l=jq_amtil[loc_m-1]-1;l<jq_amtil[loc_m]-1;l++){ //AM					 
	//printf("multimortar2 8 %e\n",au_amtil[l]);
	insertas(&irow,&mast1,&a_flagr[irow_amtil[l]-1],&m,&ifree,nzs,&au_amtil[l],&au);													  		
      }    
      break;				
    case 2 : // insert Matrix NI, MI, II, AI				    	
      loc_i=i_flag[j];
      for(l=jq_ni[loc_i-1]-1;l<jq_ni[loc_i]-1;l++){ //NI					 
	//printf("multimortar2 9 %e\n",au_ni[l]);
	insertas(&irow,&mast1,&n_flagr[irow_ni[l]-1],&m,&ifree,nzs,&au_ni[l],&au);					 			  		
      }
      for(l=jq_mitil[loc_i-1]-1;l<jq_mitil[loc_i]-1;l++){ //MI					 				 
	//printf("multimortar2 10 %e\n",au_mitil[l]);
	insertas(&irow,&mast1,&m_flagr[irow_mitil[l]-1],&m,&ifree,nzs,&au_mitil[l],&au);					  
      }		
      for(l=jq_iif[loc_i-1]-1;l<jq_iif[loc_i]-1;l++){ //II					 
	if (m==i_flagr[irow_iif[l]-1]){							
	  ad[j]=au_iif[l];					
	}else{					 
	//printf("multimortar2 11 %e\n",au_iif[l]);
	  insertas(&irow,&mast1,&i_flagr[irow_iif[l]-1],&m,&ifree,nzs,&au_iif[l],&au);					 
	}													  		
      }
      for(l=jq_aitil[loc_i-1]-1;l<jq_aitil[loc_i]-1;l++){ //AI					 
	//printf("multimortar2 12 %e\n",au_aitil[l]);
	insertas(&irow,&mast1,&a_flagr[irow_aitil[l]-1],&m,&ifree,nzs,&au_aitil[l],&au);													  		
      }  
      break;
    case 3 : // insert Matrix NA, MA, IA, AA				    	
      loc_a=a_flag[j];
      for(l=jq_na[loc_a-1]-1;l<jq_na[loc_a]-1;l++){ //NA					 
	//printf("multimortar2 13 %e\n",au_na[l]);
	insertas(&irow,&mast1,&n_flagr[irow_na[l]-1],&m,&ifree,nzs,&au_na[l],&au);					 			  		
      }
      for(l=jq_matil[loc_a-1]-1;l<jq_matil[loc_a]-1;l++){ //MA					 				 
	//printf("multimortar2 14 %e\n",au_matil[l]);
	insertas(&irow,&mast1,&m_flagr[irow_matil[l]-1],&m,&ifree,nzs,&au_matil[l],&au);					  
      }		
      for(l=jq_ia[loc_a-1]-1;l<jq_ia[loc_a]-1;l++){ //IA					 				 
	//printf("multimortar2 15 %e\n",au_ia[l]);
	insertas(&irow,&mast1,&i_flagr[irow_ia[l]-1],&m,&ifree,nzs,&au_ia[l],&au);					 													  		
      }
      for(l=jq_aatil[loc_a-1]-1;l<jq_aatil[loc_a]-1;l++){ //AA	
	//printf("\tcdof %" ITGFORMAT " lcdof %" ITGFORMAT " rdof %" ITGFORMAT " lrdof %" ITGFORMAT " v %e slav %" ITGFORMAT " %" ITGFORMAT " \n",m,loc_a,a_flagr[irow_aatil[l]-1],irow_aatil[l],au_aatil[l],islavactdof[j],islavactdof[a_flagr[irow_aatil[l]-1]-1]);
	if (m==a_flagr[irow_aatil[l]-1]){
	  //printf("\t diag! \n");
	  ad[j]=au_aatil[l];					
	}else{					 
	//printf("multimortar2 16 %e\n",au_aatil[l]);
	  insertas(&irow,&mast1,&a_flagr[irow_aatil[l]-1],&m,&ifree,nzs,&au_aatil[l],&au);					 
	}	
      }  
      break;      
    default :   
      break;		
    }		
  }
  // sort pro column 
  //   for(j=0;j<*nzs;j++){printf("multimortar2.c %d au=%f\n",j,au[j]);}
  *nzs=ifree-1;
  //   for(j=0;j<452;j++){printf("multimortar2.c %d au=%f\n",j,au[j]);}
  RENEW(au,double,*nzs);	
  RENEW(irow,ITG, *nzs);	
  //printf("trafoNT2: nzs %" ITGFORMAT "\n",*nzs);
  /*****************************************************/
  debut=clock();  
  dim=neq[1];
  matrixsort(au,mast1,irow,jq,nzs,&dim);
  SFREE(mast1);
  fin=clock();
  printf("\tmm2 tri fin : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);
  /*for (i=0;i<row_lm;i++){    
    for(j=jq_bdtil2[i]-1;j<jq_bdtil2[i+1]-1;j++){
      printf(" (%" ITGFORMAT ",%" ITGFORMAT ") %e ",a_flagr[irow_bdtil2[j]-1],m_flagr[i],au_bdtil2[j]);
    }
    printf("\n");
  } */
  
  /*********** Free the intermediate matrices ********/	
  SFREE(au_nn);SFREE(irow_nn);SFREE(jq_nn);		
  //SFREE(au_din);SFREE(irow_din);SFREE(jq_din);
  SFREE(au_antil);SFREE(irow_antil);SFREE(jq_antil);
  SFREE(au_nm);SFREE(irow_nm);SFREE(jq_nm);		
  //SFREE(au_dim);SFREE(irow_dim);SFREE(jq_dim);
  SFREE(au_amtil);SFREE(irow_amtil);SFREE(jq_amtil);
  SFREE(au_ni);SFREE(irow_ni);SFREE(jq_ni);	  
  SFREE(au_aitil);SFREE(irow_aitil);SFREE(jq_aitil);
  SFREE(au_na);SFREE(irow_na);SFREE(jq_na);	 	
  SFREE(au_dia);SFREE(irow_dia);SFREE(jq_dia);
  SFREE(au_aatil);SFREE(irow_aatil);SFREE(jq_aatil);
  SFREE(au_bd1);SFREE(irow_bd1);SFREE(jq_bd1);
  SFREE(au_bdtil2);SFREE(irow_bdtil2);SFREE(jq_bdtil2);	
  SFREE(au_ddtil2i);SFREE(irow_ddtil2i);SFREE(jq_ddtil2i);	
  SFREE(au_ddtil2a);SFREE(irow_ddtil2a);SFREE(jq_ddtil2a);	
  SFREE(ad_ii);SFREE(irow_iid);SFREE(jq_iid);
  SFREE(ad_aa);SFREE(irow_aad);SFREE(jq_aad);
  SFREE(ad_mm);SFREE(irow_mmd);SFREE(jq_mmd);	
  //SFREE(au_dii);SFREE(irow_dii);SFREE(jq_dii);
  SFREE(au_intmn);SFREE(irow_intmn);SFREE(jq_intmn);
  SFREE(au_intmm);SFREE(irow_intmm);SFREE(jq_intmm);
  SFREE(au_intmi);SFREE(irow_intmi);SFREE(jq_intmi);
  SFREE(au_intma);SFREE(irow_intma);SFREE(jq_intma);
  SFREE(au_mntil);SFREE(irow_mntil);SFREE(jq_mntil);
  SFREE(au_mmtil);SFREE(irow_mmtil);SFREE(jq_mmtil);
  SFREE(au_mitil);SFREE(irow_mitil);SFREE(jq_mitil);
  SFREE(au_matil);SFREE(irow_matil);SFREE(jq_matil);
  //SFREE(au_da);SFREE(irow_da);SFREE(jq_da);
  //SFREE(au_di);SFREE(irow_di);SFREE(jq_di);
  SFREE(au_ba);SFREE(irow_ba);SFREE(jq_ba);
  SFREE(f_a); SFREE(f_i);
  SFREE(f_m); SFREE(f_atil);
  SFREE(ad_nn);
  SFREE(au_in);SFREE(irow_in);SFREE(jq_in);
  SFREE(au_im);SFREE(irow_im);SFREE(jq_im);
  SFREE(au_ia);SFREE(irow_ia);SFREE(jq_ia);
  SFREE(au_iif);SFREE(irow_iif);SFREE(jq_iif);
  /*SFREE(au_diiinv);SFREE(jq_diiinv);SFREE(irow_diiinv);
  SFREE(au_diainv);SFREE(jq_diainv);SFREE(irow_diainv);
  SFREE(au_daiinv);SFREE(jq_daiinv);SFREE(irow_daiinv);
  SFREE(au_daainv);SFREE(jq_daainv);SFREE(irow_daainv);*/
  SFREE(au_t);SFREE(jq_t);SFREE(irow_t);
  SFREE(l_flag);SFREE(n_flag);SFREE(i_flag);SFREE(a_flag);SFREE(m_flag);	
  SFREE(n_flagr);SFREE(i_flagr);SFREE(a_flagr);SFREE(m_flagr);
  SFREE(f_di);
  /*************************/	
  /*END transmit the new stiffness matrix*/
  //*nzs=jq[neq[1]]-1;
  //  RENEW(au,double,*nzs);	
  //  RENEW(irow,ITG, *nzs);	
  printf("\tnzsc %" ITGFORMAT " nzs %" ITGFORMAT " \n",*nzsc,*nzs);
   debug=0;
   //for(j=0;j<452;j++){printf("multimortar2.c %d, irow=%d\n",j,irow[j]);}
   //   for(j=0;j<452;j++){printf("multimortar2.c %d au=%f\n",j,au[j]);}
    
  if(debug==1){	
    number=3;		
    FORTRAN(writematrix,(auc,adc,irowc,jqc,&neq[1],&number));	    	
    
    number=4;		
    FORTRAN(writematrix,(au,ad,irow,jq,&neq[1],&number));		
    printf("\n");	
    number=5;		
    FORTRAN(writevector,(b,&neq[1],&number));		
    number=6;		
    FORTRAN(writevector,(bhat,&neq[1],&number));
  }

  *irowcp = irowc; *aucp=auc;
  *irowp = irow; *aup=au;
  //FORTRAN(stop,());
  return;
}
