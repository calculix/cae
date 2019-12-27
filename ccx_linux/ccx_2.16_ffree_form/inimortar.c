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

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
/** \brief  function initializing mortar contact at the start of nonlingeo.c
 *   getting contact parameters, allocating needed fields, get results from last step,
 *   determine used mortar method, (transform and cataloque SPCs/MPCs)
 *   and generate islavelinv and  islavnodeinv
 *
 *   Author: Saskia Sitzmann
 *
 * @param [in,out] enerp		?
 * @param [in] mi		(1) max # of integration points per element (2) max degree of freedom per element 
 * @param [in] ne 		number of elements
 * @param [in] nslavs		number of slave nodes
 * @param [in] nk 		number of nodes
 * @param [in] nener 		flag indicating wheter ener is allocated
 * @param [in] ipkonp		pointer into field kon...
 * @param [in] lakonp		(i) label for element i
 * @param [in] konp 		.. for element i storing the connectivity list of elem. in succ. order
 * @param [in] nkon 		size kon
 * @param [in] maxprevcontel	maximum number of previous contact elements
 * @param [in,out] xstatep	state variables
 * @param [in] nstate_		maximum number of state variables 
 * @param [out] islavactdoftiep   (i)=tie number for active dof i
 * @param [out] bpp		friction bounds 
 * @param [out] islavactp	(i) indicates, if slave node i is active (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, =0 inactive node, =1 sticky node, =2 slipping/active node)  
 * @param [out] gapp		(i) gap for node i on slave surface
 * @param [out] slavnorp		slave normal
 * @param [out] slavtanp		slave tangent 
 * @param [out] cdispp		vector saving contact variables for frd-output 
 * @param [out] cstressp	current Lagrange multiplier 
 * @param [out] cfsp 		contact force 
 * @param [out] cfmp 		not used any more 
 * @param [out] bpinip		friction bounds at start of the increment
 * @param [out] islavactinip	islavact at the start of the increment
 * @param [out] cstressinip	Lagrange multiplier at start of the increment
 * @param [out] niwan		number of iwan elements
 * @param [in]  ntie		number of ties
 * @param [in] tieset           (1,i) name of tie constraint (2,i) dependent surface (3,i) independent surface  
 * @param [in] nslavnode	(i)pointer into field isalvnode for contact tie i 
 * @param [in] islavnode	field storing the nodes of the slave surface 
 * @param [out] islavnodeinvp     (i) slave node index for node i
 * @param [out] islavelinvp       (i)==0 if there is no slave node in the element, >0 otherwise
 * @param [out] pslavdualp	(:,i)coefficients \f$ \alpha_{ij}\f$, \f$ 1,j=1,..8\f$ for dual shape functions for face i
 * @param [out] pslavdualpgp	(:,i)coefficients \f$ \alpha_{ij}\f$, \f$ 1,j=1,..8\f$ for Petrov-Galerkin shape functions for face i 
 * @param [out] autlocp		transformation matrix \f$ T[p,q]\f$ for slave nodes \f$ p,q \f$
 * @param [out] irowtlocp		field containing row numbers of autloc
 * @param [out] jqtlocp	        pointer into field irowtloc
 * @param [out] autlocinvp	transformation matrix \f$ T^{-1}[p,q]\f$ for slave nodes \f$ p,q \f$ 
 * @param [out] irowtlocinvp	field containing row numbers of autlocinv
 * @param [out] jqtlocinvp	pointer into field irowtlocinv
 * @param [out] Bdp		coupling matrix \f$ B_d[p,q]=\int \psi_p \phi_q dS \f$, \f$ p \in S, q \in M \f$ 
 * @param [out] irowbp		field containing row numbers of Bd
 * @param [out] jqbp		pointer into field irowb
 * @param [out] Bdhelpp		coupling matrix \f$ Bhelp_d[p,q]=\tilde{D}^{-1}\tilde{B}\f$, \f$ p \in S, q \in M \f$ 
 * @param [out] irowbhelpp	field containing row numbers of Bdhelp
 * @param [out] jqbhelpp		pointer into field irowbhelp
 * @param [out] Ddp		coupling matrix \f$ D_d[p,q]=\int \psi_p \phi_q dS \f$, \f$ p,q \in S \f$ 
 * @param [out] irowdp		field containing row numbers of Dd
 * @param [out] jqdp		pointer into field irowd
 * @param [out] Ddtilp		coupling matrix \f$ \tilde{D}_d[p,q]=\int \psi_p \tilde{\phi}_q dS \f$, \f$ p,q \in S \f$ 
 * @param [out] irowdtilp	field containing row numbers of Ddtil
 * @param [out] jqdtilp		pointer into field irowdtil 
 * @param [out] Bdtilp		coupling matrix \f$ \tilde{B}_d[p,q]=\int \psi_p \tilde{\phi}_q dS \f$, \f$ p \in S, q \in M \f$ 
 * @param [out] irowbtilp	field containing row numbers of Bdtil
 * @param [out] jqbtilp		pointer into field irowbtil
 * @param [out] Bpgdp		Petrov-Galerkin coupling matrix \f$ B_d^{PG}[p,q]=\int \tilde{\phi}_p \phi_q dS \f$, \f$ p \in S, q \in M \f$ 
 * @param [out] irowbpgp		field containing row numbers of Bpgd
 * @param [out] jqbpgp		pointer into field irowbpg
 * @param [out] Dpgdp		Petrov-Galerkin coupling matrix \f$ D_d[p,q]=\int \tilde{\phi}_p \phi_q dS \f$, \f$ p,q \in S \f$ 
 * @param [out] irowdpgp		field containing row numbers of Dpgd
 * @param [out] jqdpgp		pointer into field irowdpg
 * @param [out] Dpgdtilp		transformed Petrov-Galerkin coupling matrix \f$ D_d[p,q]=\int \tilde{\phi}_p \tilde{\phi}_q dS \f$, \f$ p,q \in S \f$ 
 * @param [out] irowdpgtilp	field containing row numbers of Dpgdtil
 * @param [out] jqdpgtilp	pointer into field irowdpgtil
 * @param [out] Bpgdtilp		transformed Petrov-Galerkin coupling matrix \f$ B_d^{PG}[p,q]=\int \tilde{\phi}_p \tilde{\phi}_q dS \f$, \f$ p \in S, q \in M \f$ 
 * @param [out] irowbpgtilp	field containing row numbers of Bpgdtil
 * @param [out] jqbpgtilp	pointer into field irowbpgtil
 * @param [in]  iflagdualquad   flag indicating what mortar contact is used (=1 quad-lin, =2 quad-quad, =3 PG quad-lin, =4 PG quad-quad)
 * @param [in] itiefac 		pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i
 * @param [in] islavsurf	islavsurf(1,i) slaveface i islavsurf(2,i) # integration points generated before looking at face i 
 * @param [in] nboun            number of SPCs
 * @param [in] ndirboun		(i) direction of SPC i 
 * @param [in] nodeboun         (i) node of SPC i
 * @param [in] xboun            (i) value of SPC i
 * @param [in] nmpc		number of mpcs
 * @param [in] ipompc           (i) pointer to nodempc and coeffmpc for MPC i
 * @param [in] nodempc          nodes and directions of MPCs
 * @param [in] coefmpc          coefficients of MPCs
 * @param [in] labmpc		MPC labels 
 * @param [in] ikboun           sorted dofs idof=8*(node-1)+dir for SPCs
 * @param [in] ilboun           SPC numbers for sorted dofs
 * @param [in] ikmpc 		sorted dofs idof=8*(node-1)+dir for MPCs
 * @param [in] ilmpc		SPC numbers for sorted dofs
 * @param [out] nboun2            number of transformed SPCs
 * @param [out] ndirboun2p	(i) direction of transformed SPC i 
 * @param [out] nodeboun2p       (i) node of transformed SPC i
 * @param [out] xboun2p          (i) value of transformed SPC i
 * @param [out] nmpc2		number of transformed mpcs
 * @param [out] ipompc2p         (i) pointer to nodempc and coeffmpc for transformed MPC i
 * @param [out] nodempc2p        nodes and directions of transformed MPCs
 * @param [out] coefmpc2p        coefficients of transformed MPCs
 * @param [out] labmpc2p	transformed MPC labels 
 * @param [out] ikboun2p         sorted dofs idof=8*(node-1)+dir for transformed SPCs
 * @param [out] ilboun2p         transformed SPC numbers for sorted dofs
 * @param [out] ikmpc2p 		sorted dofs idof=8*(node-1)+dir for transformed MPCs
 * @param [out] ilmpc2p		transformed SPC numbers for sorted dofs
 * @param [out] nslavspcp	(2*i) pointer to islavspc...
 * @param [out] islavspcp         ... which stores SPCs for slave node i
 * @param [out] nslavmpcp	(2*i) pointer to islavmpc...
 * @param [out] islavmpcp	... which stores MPCs for slave node i
 * @param [out] nslavspc2p	(2*i) pointer to islavspc2...
 * @param [out] islavspc2p       ... which stores transformed SPCs for slave node i
 * @param [out] nslavmpc2p	(2*i) pointer to islavmpc2...
 * @param [out] islavmpc2p	... which stores transformed MPCs for slave node i
 * @param [out] nmastspcp	(2*i) pointer to imastspc...
 * @param [out] imastspcp         ... which stores SPCs for master node i
 * @param [out] nmastmpcp	(2*i) pointer to imastmpc...
 * @param [out] imastmpcp	... which stores MPCs for master node i 
 * @param [out] nsspc            number of SPC for slave nodes 
 * @param [out] nsspc2          number of transformed SPC for slave nodes
 * @param [out] nsmpc		number of MPC for slave nodes
 * @param [out] nsmpc2		number of transformed MPC for slave nodes 
 * @param [in] imastnode	field storing the nodes of the master surfaces
 * @param [in] nmastnode	(i)pointer into field imastnode for contact tie i 
 * @param [out] nmspc            number of SPC for master nodes
 * @param [out] nmmpc		number of MPC for master nodes 
 * @param [in] iponoels         (i) pointer to inoels
 * @param [in] inoels           (3,i) element number, local node number and pointer to another entry
 * @param [in] tietol		(1,i) tie tolerance (2,i) contant interaction material definition
 * @param [in] elcon		material parameters
 * @param [in] ncmat_		maximum number of elastic material constants 
 * @param [in] ntmat_           maximum number of temperature data points for any material 
 * @param [in] nasym		flag indicating whether matrix is symmetric or not
 * @param [in] iflag_fric	flag indicating whether iwan model is used
 * @param [out] lambdaiwanp       Lagrange multiplier splitted to Iwan elements
 * @param [out] lambdaiwaninip    Lagrange multiplier splitted to Iwan elements at start of increment 
 * @param [out] nk2		number or generated points needed for transformed SPCs
 * @param [in] vold 		displacements
 * @param [in] nset		number of sets 
 * @param [in] set		(i) name of set i
 * @param [in] mortar		param indicating what contact formulation is used (=0 NTS penalty, =1 GPTS penalty , >1 STS mortar)
 * @param [in] memmpc_		size of nodempc/coefmpc
 * @param [in] ielmatp           (j,i) material number of layer j
 * @param [in] ielorienp		(j,i) orientation number of layer j 
 * @param [in] norien		number of orientations
**/
void inimortar(double **enerp, ITG *mi, ITG *ne ,ITG *nslavs,ITG *nk,ITG *nener,
	       ITG **ipkonp, char **lakonp, ITG **konp, ITG *nkon,
	       ITG *maxprevcontel, double **xstatep, ITG *nstate_,
	       ITG **islavactdoftiep, double **bpp, ITG **islavactp,
	       double **gapp,double **slavnorp, double **slavtanp,
	       double **cdispp,
	       double **cstressp, double **cfsp, double **cfmp,
	       double **cfsinip,double **cfsinitilp,double **cfstilp,
	       double **bpinip, ITG **islavactinip, double **cstressinip,
	       ITG *niwan,ITG *ntie, char *tieset,
	       ITG *nslavnode, ITG *islavnode,
	       ITG **islavnodeinvp, ITG **islavelinvp,double **pslavdualp,
	       double **pslavdualpgp,
	       double **autlocp,ITG **irowtlocp,ITG **jqtlocp,	
	       double **autlocinvp,ITG **irowtlocinvp,ITG **jqtlocinvp,
	       double **Bdp,ITG **irowbp,ITG **jqbp,
	       double **Bdhelpp,ITG **irowbhelpp,ITG **jqbhelpp,
	       double **Ddp,ITG **irowdp,ITG **jqdp, 
	       double **Ddtilp,ITG **irowdtilp,ITG **jqdtilp,
	       double **Bdtilp,ITG **irowbtilp,ITG **jqbtilp,
	       double **Bpgdp,ITG **irowbpgp,ITG **jqbpgp,
	       double **Dpgdp,ITG **irowdpgp,ITG **jqdpgp, 
	       double **Dpgdtilp,ITG **irowdpgtilp,ITG **jqdpgtilp,
	       double **Bpgdtilp,ITG **irowbpgtilp,ITG **jqbpgtilp,
	       ITG *iflagdualquad,ITG *itiefac, ITG *islavsurf,
	       ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
	       ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
	       ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
	       ITG *nboun2,ITG **ndirboun2p,ITG **nodeboun2p,double **xboun2p,
	       ITG *nmpc2,ITG **ipompc2p, ITG **nodempc2p,double **coefmpc2p,
	       char **labmpc2p,
	       ITG **ikboun2p,ITG **ilboun2p,ITG **ikmpc2p,ITG **ilmpc2p,
	       ITG **nslavspcp,ITG **islavspcp,ITG **nslavmpcp,ITG **islavmpcp,
	       ITG **nslavspc2p,ITG **islavspc2p, ITG **nslavmpc2p,
	       ITG **islavmpc2p,
	       ITG **nmastspcp,ITG **imastspcp,ITG **nmastmpcp,ITG **imastmpcp,
	       ITG **nmastmpc2p,ITG **imastmpc2p,ITG *nmmpc2,
	       ITG *nsspc, ITG *nsspc2, ITG *nsmpc, ITG *nsmpc2,
	       ITG *imastnode,ITG *nmastnode,ITG *nmspc, ITG *nmmpc,
	       ITG *iponoels,ITG *inoels,
	       double *tietol,double *elcon,ITG *ncmat_,ITG *ntmat_,ITG *nasym,
	       ITG *iflag_fric,double **lambdaiwanp,double **lambdaiwaninip,
	       ITG *nk2, double *vold,ITG *nset,char *set,ITG *mortar,
	       ITG *memmpc_,
	       ITG **ielmatp,ITG **ielorienp,ITG *norien,ITG *nmethod,
	       ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
	       ITG **nodeforc2p,ITG **ndirforc2p,double **xforc2p,ITG *nforc2){

    char *lakon=NULL,*labmpc2=NULL; 
    
  ITG k,i,j,node,mt=mi[1]+1,regmode,regmodet,*ipkon=NULL,*kon=NULL,
      *islavactdoftie=NULL,*islavact=NULL,
      *islavactini=NULL,*islavnodeinv=NULL,*islavelinv=NULL,*irowtloc=NULL,
      *jqtloc=NULL,
      *irowtlocinv=NULL,*jqtlocinv=NULL,*irowbhelp=NULL,*jqbhelp=NULL,
      *irowb=NULL,*jqb=NULL,
      *irowd=NULL,*jqd=NULL,*irowdtil=NULL,*jqdtil=NULL,
      *irowbtil=NULL,*jqbtil=NULL,*irowbpg=NULL,*jqbpg=NULL,
      *irowdpg=NULL,*jqdpg=NULL,*irowdpgtil=NULL,*jqdpgtil=NULL,
      *irowbpgtil=NULL,*jqbpgtil=NULL,*ndirboun2=NULL,*nodeboun2=NULL,
      *ipompc2=NULL,*nodempc2=NULL,
      *ikboun2=NULL,*ilboun2=NULL,*ikmpc2=NULL,*ilmpc2=NULL,
      *nslavspc=NULL,*islavspc=NULL,*nslavmpc=NULL,*islavmpc=NULL,
      *nslavspc2=NULL,*islavspc2=NULL,*nslavmpc2=NULL,*islavmpc2=NULL,
      *nmastspc=NULL,*imastspc=NULL,*nmastmpc=NULL,*imastmpc=NULL,
      *nmastmpc2=NULL,*imastmpc2=NULL,*ielmat=NULL,*ielorien=NULL,
      *nodeforc2=NULL,*ndirforc2=NULL;
  
  double mu,fkninv,fktauinv,p0,beta,*ener=NULL,*xstate=NULL,*bp=NULL,
      *gap=NULL,*slavnor=NULL,*slavtan=NULL,
      *cdisp=NULL,*cstress=NULL,*cfs=NULL,*cfsini=NULL,*cfsinitil=NULL,
      *cfstil=NULL,*cfm=NULL,*bpini=NULL,*cstressini=NULL,*pslavdual=NULL,
      *pslavdualpg=NULL,*autloc=NULL,
      *autlocinv=NULL,*Bd=NULL,*Bdhelp=NULL,*Dd=NULL,*Ddtil=NULL,*Bdtil=NULL,
      *xboun2=NULL,*coefmpc2=NULL,*lambdaiwan=NULL,
      *lambdaiwanini=NULL,*Bpgd=NULL,*Dpgd=NULL,*Dpgdtil=NULL,*Bpgdtil=NULL,
      *xforc2=NULL; 
  
  ener=*enerp;ipkon=*ipkonp;lakon=*lakonp;kon=*konp;xstate=*xstatep;
  islavactdoftie=*islavactdoftiep;bp=*bpp;islavact=*islavactp;gap=*gapp;
  slavnor=*slavnorp;slavtan=*slavtanp;cdisp=*cdispp;cstress=*cstressp;
  cfs=*cfsp;cfsini=*cfsinip;cfsinitil=*cfsinitilp;cfstil=*cfstilp;
  cfm=*cfmp;bpini=*bpinip;islavactini=*islavactinip;cstressini=*cstressinip;
  islavnodeinv=*islavnodeinvp;islavelinv=*islavelinvp;pslavdual=*pslavdualp;
  pslavdualpg=*pslavdualpgp;
  autloc=*autlocp;irowtloc=*irowtlocp;jqtloc=*jqtlocp;	
  autlocinv=*autlocinvp;irowtlocinv=*irowtlocinvp;jqtlocinv=*jqtlocinvp;
  Bd=*Bdp;irowb=*irowbp;jqb=*jqbp;
  Bdhelp=*Bdhelpp;irowbhelp=*irowbhelpp;jqbhelp=*jqbhelpp;
  Dd=*Ddp;irowd=*irowdp;jqd=*jqdp;
  Ddtil=*Ddtilp;irowdtil=*irowdtilp;jqdtil=*jqdtilp;
  Bdtil=*Bdtilp;irowbtil=*irowbtilp;jqbtil=*jqbtilp;
  Bpgd=*Bpgdp;irowbpg=*irowbpgp;jqbpg=*jqbpgp;
  Dpgd=*Dpgdp;irowdpg=*irowdpgp;jqdpg=*jqdpgp;
  Dpgdtil=*Dpgdtilp;irowdpgtil=*irowdpgtilp;jqdpgtil=*jqdpgtilp;
  Bpgdtil=*Bpgdtilp;irowbpgtil=*irowbpgtilp;jqbpgtil=*jqbpgtilp;  
  ndirboun2=*ndirboun2p;nodeboun2=*nodeboun2p;xboun2=*xboun2p;
  ipompc2=*ipompc2p;nodempc2=*nodempc2p;coefmpc2=*coefmpc2p;
  labmpc2=*labmpc2p;ikboun2=*ikboun2p;ilboun2=*ilboun2p;ikmpc2=*ikmpc2p;
  ilmpc2=*ilmpc2p;
  nslavspc=*nslavspcp;islavspc=*islavspcp;nslavmpc=*nslavmpcp;
  islavmpc=*islavmpcp;
  nslavspc2=*nslavspc2p;islavspc2=*islavspc2p;nslavmpc2=*nslavmpc2p;
  islavmpc2=*islavmpc2p;
  nmastspc=*nmastspcp;imastspc=*imastspcp;nmastmpc=*nmastmpcp;
  imastmpc=*imastmpcp;nmastmpc2=*nmastmpc2p;imastmpc2=*imastmpc2p;
  lambdaiwan=*lambdaiwanp;lambdaiwanini=*lambdaiwaninip;
  ielmat=*ielmatp;ielorien=*ielorienp;
  nodeforc2=*nodeforc2p; ndirforc2=*ndirforc2p;
  xforc2=*xforc2p;

  
  fflush(stdout); 
  RENEW(ielmat,ITG,mi[2]*(*ne+*nslavs));
  for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++) ielmat[k]=1;
  if(*norien>0){      
    RENEW(ielorien,ITG,mi[2]*(*ne+*nslavs));
    for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++) ielorien[k]=0;
  }  
  
  if(*nener==1){
    RENEW(ener,double,mi[0]*(*ne+*nslavs)*2);
  }
  RENEW(ipkon,ITG,*ne+*nslavs);
  RENEW(lakon,char,8*(*ne+*nslavs));
  
  /* adding one element per slave node, similar to 
     spring elements;
     needed for output in frd-format of CDISP and CSTRES */
  
  RENEW(kon,ITG,*nkon+*nslavs);
  if((*maxprevcontel==0)&&(*nslavs!=0)){
    RENEW(xstate,double,*nstate_*mi[0]*(*ne+*nslavs));
    for(k=*nstate_*mi[0]**ne;k<*nstate_*mi[0]*(*ne+*nslavs);k++){
      xstate[k]=0.;
    }
  }      
  
  for(k=0;k<*nslavs;k++){
    ipkon[*ne+k]=*nkon+k;
    kon[*nkon+k]=islavnode[k];
    strcpy1(&lakon[8*(*ne+k)]," S    C0",8);
  }
  NNEW(islavactdoftie,ITG,nslavnode[*ntie]);
  NNEW(bp,double,nslavnode[*ntie]);
  NNEW(islavact,ITG,nslavnode[*ntie]);
  NNEW(gap,double,nslavnode[*ntie]);
  NNEW(slavnor,double,3*nslavnode[*ntie]);
  NNEW(slavtan,double,6*nslavnode[*ntie]);
  NNEW(cdisp,double,6*nslavnode[*ntie]);
  
  /* allocation of temperary fields: stores the structure
     of the stiffness matrix without mortar contact */
  
  NNEW(cstress,double,mt*nslavnode[*ntie]);
  NNEW(cfs,double,mt**nk);
  if(*nmethod==4){
   NNEW(cfsini,double,mt**nk);
   NNEW(cfstil,double,mt**nk);
   NNEW(cfsinitil,double,mt**nk);
  }
  /* fields for cutback  */
  
  NNEW(bpini,double,nslavnode[*ntie]);
  NNEW(islavactini,ITG,nslavnode[*ntie]);  
  NNEW(cstressini,double,mt*nslavnode[*ntie]);
  
  /* check for friction */
  
  *iflag_fric=0;
  printf("\t***mortar***\n");
  *niwan=1;
  for (i=0;i<*ntie;i++){
    if(tieset[i*(81*3)+80]=='C'){
      FORTRAN(getcontactparams,(&mu,&regmode,&regmodet,&fkninv,&fktauinv,
				&p0,&beta,tietol,elcon,&i,ncmat_,ntmat_,&k));
      printf("\ttie %" ITGFORMAT " mu %e stsl %e  niwan %" ITGFORMAT
	     " regmode %" ITGFORMAT " %" ITGFORMAT " kinv %e p0 %e beta %e \n",
	     i+1,mu,1.0/fktauinv,k,regmode,regmodet,fkninv, p0, beta );
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	islavactdoftie[j]=i;
      }
      if(mu>1.e-10 && regmodet==2){*iflag_fric=1;*niwan=max(*niwan,k);}
    }
  }	
  printf("\tiflag_firc %" ITGFORMAT " niwan %" ITGFORMAT
	 " \n \t************\n",*iflag_fric,*niwan);
  if(*iflag_fric==1){
    NNEW(lambdaiwan,double,3**niwan*nslavnode[*ntie]);
    NNEW(lambdaiwanini,double,3**niwan*nslavnode[*ntie]);
  }
  
  /* get results from last step */
  
  if(*maxprevcontel!=0){
    for (i=0;i<*ntie;i++){
      if(tieset[i*(81*3)+80]=='C'){
	if(*nstate_*mi[0]>0){
	  for(j=nslavnode[i];j<nslavnode[i+1];j++){
	    for(k=0;k<3;k++){
	      cstress[mt*j+k]=xstate[*nstate_*mi[0]*(*ne+j)+k];
	      xstate[*nstate_*mi[0]*(*ne+j)+k]=0.;		
	    }
	    if(xstate[*nstate_*mi[0]*(*ne+j)+3]>1.0 && xstate[*nstate_*mi[0]*(*ne+j)+3]<2.0){
	      islavact[j]=1;
	    }else if(xstate[*nstate_*mi[0]*(*ne+j)+3]>2.0 && xstate[*nstate_*mi[0]*(*ne+j)+3]<3.0){
	      islavact[j]=2;
	    }
	    xstate[*nstate_*mi[0]*(*ne+j)+3]=0.;
	    if(*iflag_fric==1){
	      for(k=0;k<3**niwan;k++){
		lambdaiwan[3**niwan*j+k]=xstate[*nstate_*mi[0]*(*ne+j)+4+k];
		lambdaiwanini[3**niwan*j+k]=xstate[*nstate_*mi[0]*(*ne+j)+4+k];
		xstate[*nstate_*mi[0]*(*ne+j)+4+k]=0.0;
	      }
	    }  
	  }
	}
      }
    }	
  }
  NNEW(islavnodeinv,ITG,*nk);
  NNEW(islavelinv,ITG,*ne);
  for( i=0; i<*ntie; i++){    
    if(tieset[i*(81*3)+80]=='C'){ 
      for(j=nmastnode[i]; j<nmastnode[i+1]; j++){
	node=imastnode[j];
	islavnodeinv[node-1]=-(j+1);
      }				
    }
  }    
  for( i=0; i<*ntie; i++){    
    if(tieset[i*(81*3)+80]=='C'){ 
      for(j=nslavnode[i]; j<nslavnode[i+1]; j++){
	node=islavnode[j];
	islavnodeinv[node-1]=j+1;
      }			
    }
  }
  /* coeffs for dual basis functions */
  NNEW(pslavdual,double,64*itiefac[2**ntie-1]);
  
  /*  T and T^-1 and coupling matrices in nodes  */
  NNEW(autloc,double, 3*nslavnode[*ntie]);
  NNEW(irowtloc,ITG,3*nslavnode[*ntie]);
  NNEW(jqtloc,ITG, *nk+1);	
  NNEW(autlocinv,double, 3*nslavnode[*ntie]);
  NNEW(irowtlocinv,ITG, 3*nslavnode[*ntie]);
  NNEW(jqtlocinv,ITG, *nk+1);
  NNEW(Bd,double,1);
  NNEW(irowb,ITG,1);
  NNEW(jqb,ITG,*nk+1);
  NNEW(Bdhelp,double,1);
  NNEW(irowbhelp,ITG,1);
  NNEW(jqbhelp,ITG,*nk+1);
  NNEW(Dd,double,1);
  NNEW(irowd,ITG,1);
  NNEW(jqd,ITG,*nk+1);
  NNEW(Ddtil,double,1);
  NNEW(irowdtil,ITG,1);
  NNEW(jqdtil,ITG,*nk+1);
  NNEW(Bdtil,double,1);
  NNEW(irowbtil,ITG,1);
  NNEW(jqbtil,ITG,*nk+1);
  NNEW(Bpgd,double,1);
  NNEW(irowbpg,ITG,1);
  NNEW(jqbpg,ITG,*nk+1);
  NNEW(Dpgd,double,1);
  NNEW(irowdpg,ITG,1);
  NNEW(jqdpg,ITG,*nk+1);
  NNEW(Dpgdtil,double,1);
  NNEW(irowdpgtil,ITG,1);
  NNEW(jqdpgtil,ITG,*nk+1);
  NNEW(Bpgdtil,double,1);
  NNEW(irowbpgtil,ITG,1);
  NNEW(jqbpgtil,ITG,*nk+1);
  
  /* iflagdualquad==1 : linear-quadratic (mortar==3: linmortar)
   * iflagdualquad==2 : quadratric-quadratic (mortar==2: mortar)
   * iflagdualquad==3 : petrov-galerkin lin-quad (mortar==4: pglinmortar)
   * iflagdualquad==4 : petrov-galerkin quad-quad (mortar==5: pgmortar) */
  
  if(*mortar==2){
    *iflagdualquad=2;
  }else if(*mortar==3){
    *iflagdualquad=1;
  }else if(*mortar==4){
    *iflagdualquad=3;
  }else if(*mortar==5){
    *iflagdualquad=4;
  }else{
    *iflagdualquad=2;
  }
  
  if(*iflagdualquad==3 ||*iflagdualquad==4){
  NNEW(pslavdualpg,double,64*itiefac[2**ntie-1]);
  }else{
  NNEW(pslavdualpg,double,64); 
  }
  
  buildtquad(ntie,ipkon,kon,nk,lakon,nslavnode,itiefac,tieset,
	     islavnode,islavsurf,&irowtloc,jqtloc,&autloc,
	     &irowtlocinv,jqtlocinv,&autlocinv,iflagdualquad);
  
  /* checking for SPC's and MPC's on slave and master surface */
  
  NNEW(ndirboun2,ITG, 1);
  NNEW(nodeboun2,ITG, 1);
  NNEW(xboun2,double,1);
  NNEW(ipompc2,ITG, 1);
  NNEW(nodempc2,ITG, 1);
  NNEW(coefmpc2,double,1);
  NNEW(ikboun2,ITG, 1);
  NNEW(ilboun2,ITG, 1);
  NNEW(ikmpc2,ITG, 1);
  NNEW(ilmpc2,ITG, 1);
  NNEW(labmpc2,char,20*1);
  NNEW(nodeforc2,ITG, 2);
  NNEW(ndirforc2,ITG, 1);
  NNEW(xforc2,double, 1);

  NNEW(nslavspc,ITG,2*nslavnode[*ntie]);
  NNEW(islavspc,ITG,2**nboun);
  NNEW(nslavmpc,ITG,2*nslavnode[*ntie]);
  NNEW(islavmpc,ITG,2**nmpc);
  NNEW(nslavspc2,ITG,2*nslavnode[*ntie]);
  NNEW(islavspc2,ITG,2**nboun);
  NNEW(nslavmpc2,ITG,2*nslavnode[*ntie]);
  NNEW(islavmpc2,ITG,2**nmpc);
  NNEW(nmastspc,ITG,2*nmastnode[*ntie]);
  NNEW(imastspc,ITG,2**nboun);
  NNEW(nmastmpc,ITG,2*nmastnode[*ntie]);
  NNEW(imastmpc,ITG,2**nmpc); 
  NNEW(nmastmpc2,ITG,2*nmastnode[*ntie]);
  NNEW(imastmpc2,ITG,2**nmpc);  
  
  FORTRAN(genislavelinv,(islavelinv,jqtloc,lakon,ipkon,kon,ne,nasym));
  
  *enerp=ener;*ipkonp=ipkon;*lakonp=lakon;*konp=kon;*xstatep=xstate;
  *islavactdoftiep=islavactdoftie;*bpp=bp;*islavactp=islavact;*gapp=gap;
  *slavnorp=slavnor;*slavtanp=slavtan;*cdispp=cdisp;*cstressp=cstress;
  *cfsp=cfs;*cfsinip=cfsini;*cfsinitilp=cfsinitil;*cfstilp=cfstil;
  *cfmp=cfm;*bpinip=bpini;*islavactinip=islavactini;*cstressinip=cstressini;
  *islavnodeinvp=islavnodeinv;*islavelinvp=islavelinv;*pslavdualp=pslavdual;
  *pslavdualpgp=pslavdualpg;
  *autlocp=autloc;*irowtlocp=irowtloc;*jqtlocp=jqtloc;	
  *autlocinvp=autlocinv;*irowtlocinvp=irowtlocinv;*jqtlocinvp=jqtlocinv;
  *Bdp=Bd;*irowbp=irowb;*jqbp=jqb;
  *Bdhelpp=Bdhelp;*irowbhelpp=irowbhelp;*jqbhelpp=jqbhelp;
  *Ddp=Dd;*irowdp=irowd;*jqdp=jqd;
  *Ddtilp=Ddtil;*irowdtilp=irowdtil;*jqdtilp=jqdtil;
  *Bdtilp=Bdtil;*irowbtilp=irowbtil;*jqbtilp=jqbtil;
  *Bpgdp=Bpgd;*irowbpgp=irowbpg;*jqbpgp=jqbpg;
  *Dpgdp=Dpgd;*irowdpgp=irowdpg;*jqdpgp=jqdpg;
  *Dpgdtilp=Dpgdtil;*irowdpgtilp=irowdpgtil;*jqdpgtilp=jqdpgtil;
  *Bpgdtilp=Bpgdtil;*irowbpgtilp=irowbpgtil;*jqbpgtilp=jqbpgtil;
  *ndirboun2p=ndirboun2;*nodeboun2p=nodeboun2;*xboun2p=xboun2;
  *ipompc2p=ipompc2;*nodempc2p=nodempc2;*coefmpc2p=coefmpc2;
  *labmpc2p=labmpc2;*ikboun2p=ikboun2;*ilboun2p=ilboun2;*ikmpc2p=ikmpc2;
  *ilmpc2p=ilmpc2;
  *nslavspcp=nslavspc;*islavspcp=islavspc;*nslavmpcp=nslavmpc;
  *islavmpcp=islavmpc;
  *nslavspc2p=nslavspc2;*islavspc2p=islavspc2;*nslavmpc2p=nslavmpc2;
  *islavmpc2p=islavmpc2;
  *nmastspcp=nmastspc;*imastspcp=imastspc;*nmastmpcp=nmastmpc;
  *imastmpcp=imastmpc;*nmastmpc2p=nmastmpc2;*imastmpc2p=imastmpc2;
  *lambdaiwanp=lambdaiwan;*lambdaiwaninip=lambdaiwanini;
  *ielmatp=ielmat;*ielorienp=ielorien;
  *nodeforc2p=nodeforc2;*ndirforc2p=ndirforc2;
  *xforc2p=xforc2;
  
  return;
}
