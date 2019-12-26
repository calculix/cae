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

/**
 * \brief function to include contact conditions with the dual mortar method in the transformed system
 *        see phd-thesis Sitzmann, Algorithm 2, p.71
 *
 * Author: Saskia Sitzmann
 * 
 * @param [in] ncont		number of triangles in triagulation of master side
 * @param [in] ntie		number of ties
 * @param [in] tieset           (1,i) name of tie constraint (2,i) dependent surface (3,i) independent surface
 * @param [in] nset		number of sets
 * @param [in] set		(i) name of set i 
 * @param [in] istartset        (i) pointer to ialset containing the first set member
 * @param [in] iendset          (i) pointer to ialset containing the last set member
 * @param [in] ialset           set members
 * @param [in] itietri          (1,i) first triangle in field koncont for contact contraint i (2,i) last one
 * @param [in] lakon		(i) label for element i
 * @param [in] ipkon		pointer into field kon...
 * @param [in] kon 		.. for element i storing the connectivity list of elem. in succ. order 
 * @param [in] koncont		connectivity for master surface triangulation
 * @param [in] ne 		number of elements
 * @param [in] cg		center of gravity for slave faces
 * @param [in] straight		plane equations for master surface triangulation
 * @param [in] co		coordinates of nodes
 * @param [in] vold		displacement of nodes 
 * @param [in] ielmat           (j,i) material number of layer j
 * @param [in] elcon		material parameters
 * @param [in] istep		step number
 * @param [in] iinc		increment number
 * @param [in] iit		iteration number of Newton-Raphson iteration
 * @param [in] ncmat_		maximum number of elastic material constants 
 * @param [in] ntmat_           maximum number of temperature data points for any material
 * @param [in] ne0              number of elements without contact elements
 * @param [in] vini 		displacements at the start of the increment
 * @param [in] nmethod		analysis method
 * @param [in] neq		number of active degrees of freedom
 * @param [in] nzs 		number of nonzero,nondiagonal entries of matrix K 
 * @param [in] nactdof 		(i,j) actual degree of freedom for direction i of node j
 * @param [in] itiefac 		pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i
 * @param [in] islavsurf	islavsurf(1,i) slaveface i islavsurf(2,i) # integration points generated before looking at face i 
 * @param [in] islavnode	field storing the nodes of the slave surface
 * @param [in] imastnode	field storing the nodes of the master surfaces
 * @param [in] nslavnode	(i)pointer into field isalvnode for contact tie i 
 * @param [in] nmastnode	(i)pointer into field imastnode for contact tie i 
 * @param [in,out] ad		diagonal terms of system matrix K
 * @param [in,out] aup           system matrix K
 * @param [in,out] b 		right hand side of linear system
 * @param [in,out] irowp   	rows for system matrix K
 * @param [in,out] icol 	colums for systme matrix K
 * @param [in,out] jq 		pointer to irow
 * @param [in] imastop		connectivty of master triangulation jumping from triangle to triange via common edge
 * @param [in] iponoels         (i) pointer to inoels
 * @param [in] inoels           (3,i) element number, local node number and pointer to another entry
 * @param [out] nzsc		number of nonzero,nondiagonal entries of intermediate system matrix
 * @param [out] aucp		intermediate system matrix
 * @param [out] adc             intermediate system matrix, diagonal terms
 * @param [out] irowcp           rows for intermediate system matrix
 * @param [out] jqc		pointer to irowc
 * @param [in,out] islavact	(i) indicates, if slave node i is active (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, =0 inactive node, =1 sticky node, =2 slipping/active node) 
 * @param [in,out] gap		(i) gap for node i on slave surface
 * @param [in,out] slavnor		slave normal
 * @param [in,out] slavtan		slave tangent 
 * @param [out] bhat            intermediate right hand side 
 * @param [out] irowbdp		field containing row numbers of aubd
 * @param [out] jqbd		pointer into field irowbd
 * @param [out] aubdp		coupling matrix \f$ B_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [out] irowbdtilp	field containing row numbers of aubd
 * @param [out] jqbdtil		pointer into field irowbdtil
 * @param [out] aubdtilp		matrix \f$ \tilde{D}^{-1}\tilde{B}_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [out] irowbdtil2p	field containing row numbers of aubdtil2
 * @param [out] jqbdtil2	pointer into field irowbdtil2
 * @param [out] aubdtil2p	coupling matrix \f$ \tilde{D}$ and $\tilde{B}^2_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [out] irowddp		field containing row numbers of audd
 * @param [out] jqdd		pointer into field irowdd
 * @param [out] auddp		coupling matrix \f$ D_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [out] irowddtilp	field containing row numbers of audd
 * @param [out] jqddtil		pointer into field irowdd
 * @param [out] auddtilp		coupling matrix \f$ \tilde{D}_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [out] irowddtil2p	field containing row numbers of audd
 * @param [out] jqddtil2	pointer into field irowdd
 * @param [out] auddtil2p	matrix \f$ Id_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [out] irowddinvp	field containing row numbers of auddinv
 * @param [out] jqddinv		pointer into field irowddinv
 * @param [out] auddinvp		coupling matrix \f$ \tilde{D}^{-1}_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [in] irowtloc		field containing row numbers of autloc
 * @param [in] jqtloc	        pointer into field irowtloc
 * @param [in] autloc		transformation matrix \f$ T[p,q]\f$ for slave nodes \f$ p,q \f$ 
 * @param [in] irowtlocinv	field containing row numbers of autlocinv
 * @param [in] jqtlocinv	pointer into field irowtlocinv
 * @param [in] autlocinv	transformation matrix \f$ T^{-1}[p,q]\f$ for slave nodes \f$ p,q \f$  
 * @param [in] mi		(1) max # of integration points per element (2) max degree of freedom per element
 * @param [in] ipe		(i) pointer to ime for node i 
 * @param [in] ime              ... cataloging the edges with node i
 * @param [in] tietol		(1,i) tie tolerance (2,i) contant interaction material definition
 * @param [in] iflagact         here: flag indicating if coupling matrices should be updated every iteration or only once per increment (==0)
 * @param [in] cstress		current Lagrange multiplier 
 * @param [in] cstressini	Lagrange multiplier at start of the increment
 * @param [in] bp_old		old friction bounds
 * @param [in] iflag_fric	flag indicating if iwan friction model is used
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
 * @param [in] nboun2            number of transformed SPCs
 * @param [in] ndirboun2		(i) direction of transformed SPC i 
 * @param [in] nodeboun2         (i) node of transformed SPC i
 * @param [in] xboun2            (i) value of transformed SPC i
 * @param [in] nmpc2		number of transformed mpcs
 * @param [in] ipompc2           (i) pointer to nodempc and coeffmpc for transformed MPC i
 * @param [in] nodempc2          nodes and directions of transformed MPCs
 * @param [in] coefmpc2          coefficients of transformed MPCs
 * @param [in] ikboun2           sorted dofs idof=8*(node-1)+dir for transformed SPCs
 * @param [in] ilboun2           transformed SPC numbers for sorted dofs
 * @param [in] ikmpc2 		sorted dofs idof=8*(node-1)+dir for transformed MPCs
 * @param [in] ilmpc2		transformed SPC numbers for sorted dofs
 * @param [in] nslavspc		(2*i) pointer to islavspc...
 * @param [in] islavspc         ... which stores SPCs for slave node i
 * @param [in] nsspc            number of SPC for slave nodes
 * @param [in] nslavmpc		(2*i) pointer to islavmpc...
 * @param [in] islavmpc		... which stores MPCs for slave node i
 * @param [in] nsmpc		number of MPC for slave nodes
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
 * @param [in] pslavdual	(:,i)coefficients \f$ \alpha_{ij}\f$, \f$ 1,j=1,..8\f$ for dual shape functions for face i
 * @param [in] pslavdualpg	(:,i)coefficients \f$ \alpha_{ij}\f$, \f$ 1,j=1,..8\f$ for Petrov-Galerkin shape functions for face i 
 * @param [in] islavactdof      (i)=10*slavenodenumber+direction for active dof i
 * @param [in] islavactdoftie   (i)=tie number for active dof i
 * @param [in] plicon		isotropic hardening curve or points for pressure-overclosure=tabular
 * @param [in] nplicon          isotropic hardening curve. 
 * @param [in] npmat_		maximum number of data points for plicon
 * @param [in] nelcon           (1,i) number of elastic constants for material i (2,i) number of temperature points 
 * @param [in] dtime		delta time
 * @param [in] islavnodeinv     (i) slave node index for node i
 * @param [out] Bdp		coupling matrix \f$ B_d[p,q]=\int \psi_p \phi_q dS \f$, \f$ p \in S, q \in M \f$ 
 * @param [out] irowbp		field containing row numbers of Bd
 * @param [out] jqb		pointer into field irowb
 * @param [out] Bdhelpp		coupling matrix \f$ Bhelp_d[p,q]=\tilde{D}^{-1}\tilde{B}\f$, \f$ p \in S, q \in M \f$ 
 * @param [out] irowbhelpp	field containing row numbers of Bdhelp
 * @param [out] jqbhelp		pointer into field irowbhelp
 * @param [out] Ddp		coupling matrix \f$ D_d[p,q]=\int \psi_p \phi_q dS \f$, \f$ p,q \in S \f$ 
 * @param [out] irowdp		field containing row numbers of Dd
 * @param [out] jqd		pointer into field irowd
 * @param [out] Ddtilp		coupling matrix \f$ \tilde{D}_d[p,q]=\int \psi_p \tilde{\phi}_q dS \f$, \f$ p,q \in S \f$ 
 * @param [out] irowdtilp	field containing row numbers of Ddtil
 * @param [out] jqdtil		pointer into field irowdtil 
 * @param [out] Bdtilp		coupling matrix \f$ \tilde{B}_d[p,q]=\int \psi_p \tilde{\phi}_q dS \f$, \f$ p \in S, q \in M \f$ 
 * @param [out] irowbtilp	field containing row numbers of Bdtil
 * @param [out] jqbtil		pointer into field irowbtil
 * @param [out] Bpgdp		Petrov-Galerkin coupling matrix \f$ B_d^{PG}[p,q]=\int \tilde{\phi}_p \phi_q dS \f$, \f$ p \in S, q \in M \f$ 
 * @param [out] irowbpgp		field containing row numbers of Bpgd
 * @param [out] jqbpg		pointer into field irowbpg
 * @param [out] Dpgdp		Petrov-Galerkin coupling matrix \f$ D_d[p,q]=\int \tilde{\phi}_p \phi_q dS \f$, \f$ p,q \in S \f$ 
 * @param [out] irowdpgp		field containing row numbers of Dpgd
 * @param [out] jqdpg		pointer into field irowdpg
 * @param [out] Dpgdtilp		transformed Petrov-Galerkin coupling matrix \f$ D_d[p,q]=\int \tilde{\phi}_p \tilde{\phi}_q dS \f$, \f$ p,q \in S \f$ 
 * @param [out] irowdpgtilp	field containing row numbers of Dpgdtil
 * @param [out] jqdpgtil	pointer into field irowdpgtil
 * @param [out] Bpgdtilp		transformed Petrov-Galerkin coupling matrix \f$ B_d^{PG}[p,q]=\int \tilde{\phi}_p \tilde{\phi}_q dS \f$, \f$ p \in S, q \in M \f$ 
 * @param [out] irowbpgtilp	field containing row numbers of Bpgdtil
 * @param [out] jqbpgtil	pointer into field irowbpgtil
 * @param [in] lambdaiwan       Lagrange multiplier splitted to Iwan elements
 * @param [in] lambdaiwanini    Lagrange multiplier splitted to Iwan elements at start of increment
 * @param [in] bet		parameter used in alpha-method
 * @param [in]  iflagdualquad   flag indicating what mortar contact is used (=1 quad-lin, =2 quad-quad, =3 PG quad-lin, =4 PG quad-quad)
 * @param [in]  labmpc2
 * @param [in,out]  cfsinitil \f$ \tilde{\Phi}_{c,j}\f$ contact forces from last increment, needed for dynamic calculations  
 * @param [in]  reltime       relative step time, needed for shrink
 */
void contactmortar(ITG *ncont, ITG *ntie, char *tieset, ITG *nset, char *set,
		   ITG *istartset, ITG *iendset, ITG *ialset, ITG *itietri,
		   char *lakon, ITG *ipkon, ITG *kon, ITG *koncont, ITG *ne,
		   double *cg, double *straight, double *co,
		   double *vold, ITG *ielmat,double *elcon,
		   ITG *istep,ITG *iinc,ITG *iit,ITG *ncmat_,ITG *ntmat_,
		   ITG *ne0, double *vini,
		   ITG *nmethod,ITG *neq, ITG *nzs, ITG *nactdof, ITG *itiefac,
		   ITG *islavsurf, ITG *islavnode, ITG *imastnode,
		   ITG *nslavnode, ITG *nmastnode, double *ad,
		   double **aup, double *b, ITG **irowp, ITG *icol, ITG *jq, ITG *imastop,
		   ITG *iponoels, ITG *inoels, ITG *nzsc, double **aucp,
		   double *adc, ITG **irowcp, ITG *jqc, ITG *islavact,
		   double *gap,
		   double *slavnor,double *slavtan, 
		   double *bhat,
		   ITG **irowbdp, ITG *jqbd, double **aubdp,
		   ITG **irowbdtilp, ITG *jqbdtil,double **aubdtilp,
		   ITG **irowbdtil2p, ITG *jqbdtil2,double **aubdtil2p,
		   ITG **irowddp, ITG *jqdd,double **auddp, 
		   ITG **irowddtilp, ITG *jqddtil,double **auddtilp, 
		   ITG **irowddtil2p, ITG *jqddtil2,double **auddtil2p,
		   ITG **irowddinvp, ITG *jqddinv,double **auddinvp,
		   ITG *irowtloc, ITG *jqtloc,double *autloc,   
		   ITG *irowtlocinv, ITG *jqtlocinv,double *autlocinv,    
		   ITG *mi,ITG *ipe, ITG *ime,double *tietol,ITG *iflagact,double *cstress,
		   double *cstressini,double *bp_old,ITG *iflag_fric, ITG *nk, 
		   ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
		   ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,
		   ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
		   ITG *nboun2,ITG *ndirboun2,ITG *nodeboun2,double *xboun2,
		   ITG *nmpc2,ITG *ipompc2,ITG *nodempc2,double *coefmpc2,
		   ITG *ikboun2,ITG *ilboun2,ITG *ikmpc2,ITG *ilmpc2,
		   ITG *nslavspc,ITG *islavspc,ITG *nsspc,ITG *nslavmpc,ITG *islavmpc,ITG *nsmpc,
		   ITG *nslavspc2,ITG *islavspc2,ITG *nsspc2,ITG *nslavmpc2,ITG *islavmpc2,ITG *nsmpc2,
		   ITG *nmastspc,ITG *imastspc,ITG *nmspc,ITG *nmastmpc,ITG *imastmpc,ITG *nmmpc,
		   ITG *nmastmpc2,ITG *imastmpc2,ITG *nmmpc2,
		   double *pslavdual,double *pslavdualpg,
		   ITG *islavactdof,
		   ITG *islavactdoftie,
		   double *plicon,ITG *nplicon, ITG *npmat_,ITG *nelcon, double *dtime,
		   ITG *islavnodeinv,
		   double **Bdp, ITG **irowbp,ITG *jqb,
		   double **Bdhelpp, ITG **irowbhelpp,ITG *jqbhelp,
		   double **Ddp, ITG **irowdp,ITG *jqd,
		   double **Ddtilp, ITG **irowdtilp,ITG *jqdtil,
		   double **Bdtilp, ITG **irowbtilp,ITG *jqbtil,
		   double **Bpgdp, ITG **irowbpgp,ITG *jqbpg,
		   double **Dpgdp, ITG **irowdpgp,ITG *jqdpg,
		   double **Dpgdtilp, ITG **irowdpgtilp,ITG *jqdpgtil,
		   double **Bpgdtilp, ITG **irowbpgtilp,ITG *jqbpgtil,
		   double *lambdaiwan,double *lambdaiwanini, double *bet, ITG *iflagdualquad,
		   char *labmpc2, double *cfsinitil,double *reltime, ITG *ithermal,
		   double *plkcon, ITG *nplkcon){
  
  ITG i,j,k,numb,ntrimax,*nx=NULL,*ny=NULL,*nz=NULL,nintpoint=0,calcul_fn,calcul_f,
    nzsbd,*irowbd=NULL,*irowdd=NULL,*irowddinv=NULL,*irowddtil=NULL,*irowbdtil=NULL,
    nzs2,niwan,*irowddtil2=NULL,*irowbdtil2=NULL,
    l,nstart,kflag,ntri,ii,number,regmode,derivmode,regmodet=1,
    *irowc=NULL,*imastsurf=NULL,imastk2,num_cpus=1,
    *irow=NULL,*irowqdt=NULL, *irowb=NULL,*irowbhelp=NULL, *irowd=NULL,*irowdtil=NULL, *irowbtil=NULL,
    *irowbpg=NULL, *irowdpg=NULL,*irowdpgtil=NULL, *irowbpgtil=NULL,
    debug,nacti, ninacti, nnogap,nstick,nnolm,nnoslav,nzsbdtil,nzsbdtil2;
  
  double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL,*aubd=NULL,*cstresstil=NULL,scal,
    *audd=NULL, *auddtil=NULL, *auddtil2=NULL, *auddinv=NULL,
    *auc=NULL, *pmastsurf=NULL,*auqdt=NULL,*gapmints=NULL,
    *au=NULL,*pslavsurf=NULL,*aubdtil=NULL,*aubdtil2=NULL,
    *Bd=NULL,*Bdhelp=NULL, *Dd=NULL,*Ddtil=NULL, *Bdtil=NULL,
    *Bpgd=NULL, *Dpgd=NULL,*Dpgdtil=NULL, *Bpgdtil=NULL,
    *areaslav=NULL,*anull=NULL,mu,fkninv,fktauinv,p0,beta,*rs=NULL,*rsb=NULL, alpha,*fmpc=NULL;
  
  double aninvloc,gnc,u_nold,u_told[2],nu_told,lnold,lt[2],nltold,constant=1.e10,atau2,sb,rb,resreg[2];
  
  double *u_old=NULL,*u_oldt=NULL;
  ITG mt=mi[1]+1,nodes,id,islavk2,dim,idof1,idof2,jj,kk;
  
  
  clock_t debut;
  clock_t fin;
  irow = *irowp; au=*aup; auc=*aucp;irowc=*irowcp;
  aubd=*aubdp; irowbd=*irowbdp;
  aubdtil=*aubdtilp; irowbdtil=*irowbdtilp;
  aubdtil2=*aubdtil2p; irowbdtil2=*irowbdtil2p;
  irowdd = *irowddp; audd=*auddp; 
  irowddtil = *irowddtilp; auddtil=*auddtilp; 
  irowddtil2 = *irowddtil2p; auddtil2=*auddtil2p;
  irowddinv = *irowddinvp; auddinv=*auddinvp; 
  irowb = *irowbp; Bd =*Bdp;irowbhelp = *irowbhelpp; Bdhelp =*Bdhelpp;
  irowd = *irowdp; Dd =*Ddp;
  irowdtil = *irowdtilp; Ddtil =*Ddtilp;
  irowbtil = *irowbtilp; Bdtil =*Bdtilp;
  irowbpg = *irowbpgp; Bpgd =*Bpgdp;
  irowdpg = *irowdpgp; Dpgd =*Dpgdp;
  irowdpgtil = *irowdpgtilp; Dpgdtil =*Dpgdtilp;
  irowbpgtil = *irowbpgtilp; Bpgdtil =*Bpgdtilp;
  
  
  debug=0;
  NNEW(rs,double,(mt)**nk);
  NNEW(rsb,double,neq[1]);
  printf(" contactmortar: start\n");
  /* coupling the slave degrees of freedom with the corresponding
     slave nodes and doing the same for the master nodes */
  FORTRAN(genislavactdof,(ntie,tieset,neq,nactdof,nslavnode,islavact,
			  nmastnode,imastnode,islavactdof,
			  islavnode,mi,ithermal));
  
  /* right now iflagact is 1 in the first iteration of every increment and 0 for all subsequent iterations.
     Thus the update of the normal and tangentials as well as the segmentation of the contact surface nedded
     for the calculation of the coupling matrices is done only once per increment, since a combined fix-point
     Newton approach in implemented, see phd-thesis Saskia Sitzmann, Chapter 3 introduction  */			      
  if(*iflagact==0){	
    /* update the location of the center of gravity of 
       the master triangles and the coefficients of their
       bounding planes needed for the search algorithm in slavintmortar->neartriangle */
    FORTRAN(updatecont,(koncont,ncont,co,vold,
			cg,straight,mi));		
    /* determining the size of the auxiliary fields 
       (needed for the master triangle search for any
       given location on the slave faces) */	
    ntrimax=0;	
    for(i=0;i<*ntie;i++){	    
      if(itietri[2*i+1]-itietri[2*i]+1>ntrimax)		
	ntrimax=itietri[2*i+1]-itietri[2*i]+1;  	
    }
    //if ((*iinc==1)&&(*iit==1)){ 
    /* For the first step, first increment, first iteration an initial guess for 
       the active set is generated analogous to node-to-surface penalty */
    if ((*iinc==1)&&(*iit==1)&&(*istep==1)){	    
      debut=clock();	    
      NNEW(xo,double,ntrimax);	    
      NNEW(yo,double,ntrimax);	    
      NNEW(zo,double,ntrimax);	    
      NNEW(x,double,ntrimax);	    
      NNEW(y,double,ntrimax);	    
      NNEW(z,double,ntrimax);	   
      NNEW(nx,ITG,ntrimax);	   
      NNEW(ny,ITG,ntrimax);	    
      NNEW(nz,ITG,ntrimax);	    
      NNEW(areaslav,double,itiefac[2*(*ntie-1)+1]);	    
      ITG ifree=0;	    
      FORTRAN(genfirstactif,(tieset,ntie,itietri,ne,ipkon,kon,lakon,
			     cg,straight,koncont,
			     co,vold,xo,yo,zo,x,y,z,nx,ny,nz,ielmat,elcon,istep,
			     iinc,iit,ncmat_,ntmat_,ne0,vini,nmethod,mi,
			     imastop,nslavnode,islavnode,islavsurf,itiefac,areaslav,iponoels,
			     inoels,set,nset,istartset,iendset,ialset,islavact,&ifree,
			     tietol));
      printf("\tFrist Active Set : %" ITGFORMAT " nodes\n",ifree);	    
      fin= clock();	    
      printf("\tgenfirstactiv : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);	    
      SFREE(xo);SFREE(yo);SFREE(zo);SFREE(x);SFREE(y);SFREE(z);SFREE(nx);	    
      SFREE(ny);SFREE(nz);	    
      SFREE(areaslav);	
    }
    fflush(stdout);	
    nacti=0; ninacti=0; nnogap=0;nstick=0;nnolm=0;nnoslav=0;
    for (i=0;i<*ntie;i++){	    	    
      if(tieset[i*(81*3)+80]=='C'){	   	      
	for(j=nslavnode[i];j<nslavnode[i+1];j++){	
	  if(islavact[j]<0){islavact[j]=-3;}
	  if(islavact[j]==2){nacti++;}
	  if(islavact[j]==0){ninacti++;}
	  if(islavact[j]==-3){nnoslav++;}	    
	}	    
      }	
    }
    printf("\tcm: N_Activ: %" ITGFORMAT "\t N_stick: %" ITGFORMAT "\tN_Inactiv: %" ITGFORMAT "\t N_nogap: %" ITGFORMAT "\t N_nolm %" ITGFORMAT "\n",nacti,nstick,ninacti,nnogap,nnolm);		
    NNEW(xo,double,ntrimax);	
    NNEW(yo,double,ntrimax);	
    NNEW(zo,double,ntrimax);	
    NNEW(x,double,ntrimax);	
    NNEW(y,double,ntrimax);	
    NNEW(z,double,ntrimax);	
    NNEW(nx,ITG,ntrimax);	
    NNEW(ny,ITG,ntrimax);	
    NNEW(nz,ITG,ntrimax);	
    /* calculating the normals,tangents in the nodes of the slave
       surface */	
    debut=clock();	
    FORTRAN(gencontrel,(tieset,ntie,itietri,ipkon,kon,
			lakon,set,cg,straight,
			koncont,co,vold,nset,
			islavsurf,itiefac,
			islavnode,nslavnode,slavnor,slavtan,mi));
    fin= clock();
    printf("\tgencontrel : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);  
    
    /* Calculating the location of the matched slave/master
       integration points, see phd-thesis Saskia Sitzmann, Appendix A */
    debut=clock();    
    NNEW(imastsurf,ITG,66);	
    NNEW(gapmints,double,66);	
    NNEW(pmastsurf,double,132);	
    NNEW(pslavsurf,double,198);	
    islavsurf[1]=0;	
    //#pragma omp for
    for(i=0;i<*ntie;i++){	    
      ii=i+1;	    
      if(tieset[i*(81*3)+80]=='C'){		
	nstart=itietri[2*i]-1;		
	ntri=(itietri[2*i+1]-nstart);		
	for(j=0;j<ntri;j++){		    
	  xo[j]=cg[(nstart+j)*3];		    
	  x[j]=xo[j];		   
	  nx[j]=j+1;		    
	  yo[j]=cg[(nstart+j)*3+1];		    
	  y[j]=yo[j];		    
	  ny[j]=j+1;		    
	  zo[j]=cg[(nstart+j)*3+2];		    
	  z[j]=zo[j];		    
	  nz[j]=j+1;		
	}
	kflag=2;		
	FORTRAN(dsort,(x,nx,&ntri,&kflag));		
	FORTRAN(dsort,(y,ny,&ntri,&kflag));		
	FORTRAN(dsort,(z,nz,&ntri,&kflag));		
	for(l=itiefac[2*i];l<=itiefac[2*i+1];l++){		    
	  RENEW(imastsurf,ITG,nintpoint+ntri*8*7);		    
	  RENEW(gapmints,double,nintpoint+ntri*8*7);		    
	  RENEW(pmastsurf,double,2*(nintpoint+ntri*8*7));		    
	  RENEW(pslavsurf,double,3*(nintpoint+ntri*8*7));		    
	  FORTRAN(slavintmortar,(tieset,ntie,itietri,ipkon,kon,
				 lakon,set,cg,straight,&nintpoint,
				 koncont,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,nset,
				 iinc,iit,
				 islavsurf,imastsurf,pmastsurf,itiefac,
				 islavnode,nslavnode,slavnor,slavtan,imastop,gapmints,
				 islavact,mi,ncont,ipe,ime,pslavsurf,&ii,&l,&ntri,tietol,
				 reltime,nmethod));
	}	    
      }	
    }
    fin= clock();	
    printf("\tslavintmortar : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);	
    printf("\tnumber of slave integration points = %" ITGFORMAT "\n",nintpoint);	
    if (nintpoint!=0){	    
      RENEW(imastsurf,ITG,nintpoint);	
    }else{	    
      RENEW(imastsurf,ITG,1);	
    }	
    if (nintpoint!=0){	    
      RENEW(gapmints,double,nintpoint);	
    }else{	    
      RENEW(gapmints,double,1);	
    }
    if (nintpoint!=0){    
      RENEW(pmastsurf,double,2*nintpoint);	
    }else{	    
      RENEW(pmastsurf,double,2);	
    }		
    if (nintpoint!=0){	    
      RENEW(pslavsurf,double,3*nintpoint);	
    }else{	    
      RENEW(pslavsurf,double,3);
    }
    SFREE(xo);SFREE(yo);SFREE(zo);SFREE(x);SFREE(y);SFREE(z);SFREE(nx);     
    SFREE(ny);SFREE(nz);
    
    /* check SPC's and MPC's on slave nodes for compability and set all slave nodes involed in SPCs/MPCs to no-LM nodes */
    debut=clock();
    FORTRAN(checkspcmpc,(lakon,ipkon,kon,ntie,tieset,
			 islavnode,
			 imastnode,nslavnode,nmastnode,
			 slavnor,islavact,
			 nboun,ndirboun,nodeboun,xboun,
			 nmpc,ipompc,nodempc,coefmpc,
			 ikboun,ilboun,ikmpc,ilmpc,
			 nboun2,ndirboun2,nodeboun2,xboun2,
			 nmpc2,ipompc2,nodempc2,coefmpc2,			 
			 ikboun2,ilboun2,ikmpc2,ilmpc2,
			 nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
			 nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc));
    fin= clock();
    printf("\tcheckspcmpc : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);		
    nacti=0; ninacti=0; nnogap=0;nstick=0;nnolm=0;nnoslav=0;	
    for (i=0;i<*ntie;i++){	
      if(tieset[i*(81*3)+80]=='C'){	    
	for(j=nslavnode[i];j<nslavnode[i+1];j++){		
	  if(islavact[j]==2){nacti++;}
	    if(islavact[j]==1){nstick++;}		
	    if(islavact[j]==0){ninacti++;}
	    if(islavact[j]==-1){nnogap++;}
	    if(islavact[j]==-2){nnolm++;}
	    if(islavact[j]==-3){nnoslav++;}
	}	 
      }	 	
    }
    printf("\tcm: N_Activ: %" ITGFORMAT "\t N_stick: %" ITGFORMAT "\tN_Inactiv: %" ITGFORMAT "\t N_nogap: %" ITGFORMAT "\t N_nolm: %" ITGFORMAT "\n",nacti,nstick,ninacti,nnogap,nnolm);	      
    
    /* calculating the coeffs of dual basis functions (Sitzmann, Chapter 3.3.) and redistribute contributions of nogap
       and noLM nodes other slave nodes (Sitzmann, Chapter 4.3.)*/		
    debut=clock();	
    FORTRAN(gendualcoeffs,(tieset,ntie,ipkon,kon,
			   lakon,
			   co,vold,
			   iinc,iit,islavact,
			   islavsurf,itiefac,
			   islavnode,nslavnode,
			   mi,pslavsurf,pslavdual,pslavdualpg,iflagdualquad));
    fin= clock();	
    printf("\tgendualcoeffs : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);
    
    /* calculate all mortar coupling matrices as well as the dual gap via the segmentation of the slave surface */
    nzsbd = 6*nslavnode[*ntie];
    nzsbdtil=6*nslavnode[*ntie];
    debut=clock();
    bdfill(&irowbd, jqbd, &aubd, &nzsbd,
	   &irowbdtil, jqbdtil, &aubdtil,&nzsbdtil,
	   &irowbdtil2, jqbdtil2, &aubdtil2,&nzsbdtil2,
	   &irowdd, jqdd, &audd,
	   &irowddtil, jqddtil, &auddtil,
	   &irowddtil2, jqddtil2, &auddtil2,
	   &irowddinv, jqddinv, &auddinv, 	     
	   irowtloc, jqtloc, autloc, 
	   irowtlocinv, jqtlocinv, autlocinv, 
	   ntie,
	   ipkon, kon, lakon, nslavnode, nmastnode, imastnode, islavnode, 
	   islavsurf, imastsurf, pmastsurf, itiefac,tieset, neq, nactdof,co,vold,
	   iponoels, inoels,mi,gapmints,gap,pslavsurf,pslavdual,pslavdualpg,&nintpoint,slavnor,nk,
	   nboun,ndirboun,nodeboun,xboun,
	   nmpc,ipompc,nodempc,coefmpc,
	   ikboun,ilboun,ikmpc,ilmpc,
	   nboun2,ndirboun2,nodeboun2,xboun2,
	   nmpc2,ipompc2,nodempc2,coefmpc2,
	   ikboun2,ilboun2,ikmpc2,ilmpc2,
	   nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
	   nslavspc2,islavspc2,nsspc2,nslavmpc2,islavmpc2,nsmpc2,
	   nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
	   nmastmpc2,imastmpc2,nmmpc2,
	   iit, iinc,islavactdof,islavact,islavnodeinv,
	   &Bd,&irowb,jqb,
	   &Bdhelp,&irowbhelp,jqbhelp,
	   &Dd,&irowd,jqd,
	   &Ddtil,&irowdtil,jqdtil,
	   &Bdtil,&irowbtil,jqbtil,
	   &Bpgd,&irowbpg,jqbpg,
	   &Dpgd,&irowdpg,jqdpg,
	   &Dpgdtil,&irowdpgtil,jqdpgtil,
	   &Bpgdtil,&irowbpgtil,jqbpgtil,
	   iflagdualquad,ithermal); 
    fin= clock();
    printf("\tbdfill : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);	
    SFREE(imastsurf);SFREE(pmastsurf);SFREE(gapmints);SFREE(pslavsurf);       
    fflush(stdout);   
  }
  
  
  /* get uhat_k-1 for first increment and first iteration**/
  if(*iit==1){     
    NNEW(u_old,double,3*nslavnode[*ntie]); 
    NNEW(u_oldt,double,3*nslavnode[*ntie]);
    NNEW(cstresstil,double,mt*nslavnode[*ntie]);
    for (i=0;i<*ntie;i++){      	
      if(tieset[i*(81*3)+80]=='C'){      	
	for(j=nslavnode[i];j<nslavnode[i+1];j++){	    
	  nodes=islavnode[j];
	  for(jj=jqd[nodes-1]-1;jj<jqd[nodes-1+1]-1;jj++){      
	    u_oldt[(islavnodeinv[irowd[jj]-1]-1)*3]+=Dd[jj]*(vold[mt*(nodes)-3]-vini[mt*(nodes)-3]);
	    u_oldt[(islavnodeinv[irowd[jj]-1]-1)*3+1]+=Dd[jj]*(vold[mt*(nodes)-2]-vini[mt*(nodes)-2]);
	    u_oldt[(islavnodeinv[irowd[jj]-1]-1)*3+2]+=Dd[jj]*(vold[mt*(nodes)-1]-vini[mt*(nodes)-1]);
	    u_old[(islavnodeinv[irowd[jj]-1]-1)*3]+=Dd[jj]*(vold[mt*(nodes)-3]);
	    u_old[(islavnodeinv[irowd[jj]-1]-1)*3+1]+=Dd[jj]*(vold[mt*(nodes)-2]);
	    u_old[(islavnodeinv[irowd[jj]-1]-1)*3+2]+=Dd[jj]*(vold[mt*(nodes)-1]);
	  }	  
	}	
      }    
    }
    for (i=0;i<*ntie;i++){      	
      if(tieset[i*(81*3)+80]=='C'){      	
	for(j=nmastnode[i];j<nmastnode[i+1];j++){	    
	  nodes=imastnode[j];
	  for(jj=jqb[nodes-1]-1;jj<jqb[nodes-1+1]-1;jj++){      
	    u_oldt[(islavnodeinv[irowb[jj]-1]-1)*3]+=Bd[jj]*(vold[mt*(nodes)-3]-vini[mt*(nodes)-3]);
	    u_oldt[(islavnodeinv[irowb[jj]-1]-1)*3+1]+=Bd[jj]*(vold[mt*(nodes)-2]-vini[mt*(nodes)-2]);
	    u_oldt[(islavnodeinv[irowb[jj]-1]-1)*3+2]+=Bd[jj]*(vold[mt*(nodes)-1]-vini[mt*(nodes)-1]);
	    u_old[(islavnodeinv[irowb[jj]-1]-1)*3]+=Bd[jj]*(vold[mt*(nodes)-3]);
	    u_old[(islavnodeinv[irowb[jj]-1]-1)*3+1]+=Bd[jj]*(vold[mt*(nodes)-2]);
	    u_old[(islavnodeinv[irowb[jj]-1]-1)*3+2]+=Bd[jj]*(vold[mt*(nodes)-1]);
	  }	  
	}	
      }    
    }
    for (i=0;i<*ntie;i++){      	
      if(tieset[i*(81*3)+80]=='C'){      	
	for(j=nslavnode[i];j<nslavnode[i+1];j++){	    
	  nodes=islavnode[j];
	  for(jj=jqdtil[nodes-1]-1;jj<jqdtil[nodes-1+1]-1;jj++){     
	    cstresstil[(islavnodeinv[irowdtil[jj]-1]-1)*(mt)]+=Ddtil[jj]*cstress[(islavnodeinv[nodes-1]-1)*(mt)+0];
	    cstresstil[(islavnodeinv[irowdtil[jj]-1]-1)*(mt)+1]+=Ddtil[jj]*cstress[(islavnodeinv[nodes-1]-1)*(mt)+1];
	    cstresstil[(islavnodeinv[irowdtil[jj]-1]-1)*(mt)+2]+=Ddtil[jj]*cstress[(islavnodeinv[nodes-1]-1)*(mt)+2];
	  }	  
	}	
      }    
    }     
  }
  
  nacti=0; ninacti=0; nnogap=0;nstick=0;nnolm=0;nnoslav=0;		
  for (i=0;i<*ntie;i++){	  
    if(tieset[i*(81*3)+80]=='C'){	    
      for(j=nslavnode[i];j<nslavnode[i+1];j++){	      
	/* adjust active set for first iteration of first increment **/	        
	if((*iit==1)){		 
	  u_told[0]=u_old[(j)*3+0]*slavtan[(j*6)+0]+u_old[(j)*3+1]*slavtan[(j*6)+1]+u_old[(j)*3+2]*slavtan[(j*6)+2]; 	
	  u_told[1]=u_old[(j)*3+0]*slavtan[(j*6)+3]+u_old[(j)*3+1]*slavtan[(j*6)+4]+u_old[(j)*3+2]*slavtan[(j*6)+5];		 
	  nu_told=sqrt(u_told[0]*u_told[0]+u_told[1]*u_told[1]);
	  lnold=cstresstil[(j)*(mt)+0]*slavnor[(j*3)+0]+cstresstil[(j)*(mt)+1]*slavnor[(j*3)+1]+cstresstil[(j)*(mt)+2]*slavnor[(j*3)+2];
	  lt[0]=cstresstil[(j)*(mt)+0]*slavtan[(j*6)+0]+cstresstil[(j)*(mt)+1]*slavtan[(j*6)+1]+cstresstil[(j)*(mt)+2]*slavtan[(j*6)+2]; 	
	  lt[1]=cstresstil[(j)*(mt)+0]*slavtan[(j*6)+3]+cstresstil[(j)*(mt)+1]*slavtan[(j*6)+4]+cstresstil[(j)*(mt)+2]*slavtan[(j*6)+5];		 
	  nltold=sqrt(lt[0]*lt[0]+lt[1]*lt[1]);
	  lnold=cstresstil[(j)*(mt)+0]*slavnor[(j*3)+0]+cstresstil[(j)*(mt)+1]*slavnor[(j*3)+1]+cstresstil[(j)*(mt)+2]*slavnor[(j*3)+2];
	  FORTRAN(getcontactparams,(&mu,&regmode,&regmodet,&fkninv,&fktauinv,&p0,&beta,tietol,elcon,&i,ncmat_,ntmat_,&niwan));
	  derivmode=0;
	  if(islavact[j]>-1){scal=Ddtil[jqdtil[islavnode[j]-1]-1];}else{scal=0.0;}
	  FORTRAN(regularization_gn_c,(&lnold,&derivmode,&regmode,&gnc,&aninvloc,&p0,&beta,elcon,nelcon,&i,ntmat_,
				       plicon,nplicon,npmat_,ncmat_,tietol,&scal));
	  if(mu>1.E-10){     
	    bp_old[j]=mu*(lnold);
	    
	  }else{
	    bp_old[j]=(lnold);
	  }
	  jj=j+1;
	  nltold=sqrt((lt[0])*(lt[0])+(lt[1])*(lt[1]));
	  // in case of NO friction node must set "slip"
	  if(mu>1.E-10 ){
	    if(*iinc==1){
	      if(islavact[j]==0 && gap[j]<1.e-9 ) {islavact[j]=1;}	    
	      if(islavact[j]>0 && bp_old[j]<1.e-14){bp_old[j]=1;}
	      if(regmodet==1){
		if(islavact[j]>0 && nltold <1.e-5){islavact[j]=1;}// first step
	      }else{
		if(islavact[j]==1){islavact[j]=2;}
	      }
	    }
	  }else{
	    if(*iinc==1){
	      if(gap[j]>1E-10 && islavact[j]>0 ){islavact[j]=0; bp_old[j]=0.0;}
	      if(gap[j]<1E-10 && islavact[j]==0 ){islavact[j]=2;}
	      if(islavact[j]==1){islavact[j]=2;}
	    }
	  }	 
	}	
	if(islavact[j]==2){nacti++;}
	if(islavact[j]==1){nstick++;}		
	if(islavact[j]==0){ninacti++;}
	if(islavact[j]==-1){nnogap++;}
	if(islavact[j]==-2){nnolm++;}	
	if(islavact[j]==-3){nnoslav++;}
      }	  
    }	
  }
  
  if(*iinc==1 && *iit==1 && *istep==1 ){
    // set initial value for bp_old in first iteration of first increment of first step 
    for (i=0;i<*ntie;i++){	 
      if(tieset[i*(81*3)+80]=='C'){	 
	for(j=nslavnode[i];j<nslavnode[i+1];j++){	  
	  bp_old[j]=1.0;  	 
	}	 
      }       
    }   
  }
  if(*iit==1){SFREE(u_old);SFREE(u_oldt);SFREE(cstresstil);}
  printf("\tcm: N_Activ: %" ITGFORMAT "\t N_stick: %" ITGFORMAT "\tN_Inactiv: %" ITGFORMAT "\t N_nogap: %" ITGFORMAT "\t N_nolm: %" ITGFORMAT " N_noslav: %" ITGFORMAT "\n",nacti,nstick,ninacti,nnogap,nnolm,nnoslav);       
  
  /* modifying the stiffnes matrix K with the coupling matrices; the
     expanded (symmetric) matrix is described in asymmetric form by
     the fields auc, adc, irowc, jqc and nzsc, bhat */       
  nzsbd=jqbd[neq[1]]-1;		
  *nzsc = nzs[1];	
  debut=clock();		
  nzs2=nzs[1];
  
  
  alpha = 1-2*sqrt(*bet);
  /* alter mechanical part of residuum for dynamic calculations, see Sitzmann Chapter 5.1., Equation (5.18) */
  if(*nmethod==4){
    //if(*iinc==1 && *iit==1){
      for(i=0;i<mt**nk;i++){cfsinitil[i]=0.0;}
      for(i=0;i<*nk;i++){
	for(j=jqdtil[i]-1;j<jqdtil[i+1]-1;j++){
	  for(l=0;l<3;l++){
	    cfsinitil[mt*(i+1)-3+l]+=Ddtil[j]*cstressini[mt*(islavnodeinv[irowdtil[j]-1]-1)+l];
	  }
	}
      }
      for(i=0;i<*nk;i++){
	for(j=jqbtil[i]-1;j<jqbtil[i+1]-1;j++){
	  for(l=0;l<3;l++){
	    cfsinitil[mt*(i+1)-3+l]+=Bdtil[j]*cstressini[mt*(islavnodeinv[irowbtil[j]-1]-1)+l];
	  }
	}    
      }
    //}
    for(i=0;i<*nk*mt;i++){
      rs[i]=cfsinitil[i]*(alpha);
    }
    NNEW(fmpc,double,*nmpc);
    calcul_fn=1;
    calcul_f=1;
    resultsforc(nk,rsb,rs,nactdof,ipompc2,nodempc2,
		coefmpc2,labmpc2,nmpc2,mi,fmpc,&calcul_fn,
		&calcul_f,&num_cpus);
    SFREE(fmpc);	    
    for(i=0;i<neq[1];i++){b[i]+=rsb[i];}    
  }	    

  /* modifying the stiffnes matrix K with the coupling matrices; embedding of the contact conditions
     the expanded (symmetric) matrix is described in asymmetric form by
     the fields auc, adc, irowc, jqc and nzsc, bhat */
  /* k needed in semi-smooth Newton in tangential direction:
     k=1 stick is assumed in all nodes, k=2 stick or slip is assummed 
     according to active set entry*/
  if(*iit==1 && *iinc==1){k=1;}else{k=2;}
  if(*iflagdualquad>2){
    /* Petrov-Galerkin formulation (Sitzmann, Chapter 4.2.) */
    multimortar2(&au, ad, &irow, jq, &nzs2,
		 &auc, adc, &irowc, jqc, nzsc,
		 aubd,irowbd,jqbd,
		 aubdtil,irowbdtil,jqbdtil,
		 aubdtil2,irowbdtil2,jqbdtil2,
		 irowdd, jqdd, audd,
		 irowddtil2, jqddtil2, auddtil2,
		 irowddinv, jqddinv, auddinv,
		 Bpgd,irowbpg,jqbpg,
		 Dpgd,irowdpg,jqdpg,
		 Ddtil,irowdtil,jqdtil,
		 neq,b,bhat,islavnode,imastnode,nslavnode,nmastnode,
		 islavact,islavactdof,
		 gap,
		 slavnor,slavtan,
		 vold,vini,cstress,cstressini,
		 bp_old,nactdof,ntie,mi,nk,
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
		 islavnodeinv,lambdaiwan,lambdaiwanini,&k,nmethod,bet,
		 ithermal,plkcon,nplkcon);   
    
  }else{
    /* normal formulation (Sitzmann, Chapter 4.1.)*/
    multimortar2(&au, ad, &irow, jq, &nzs2,
		 &auc, adc, &irowc, jqc, nzsc,
		 aubd,irowbd,jqbd,
		 aubdtil,irowbdtil,jqbdtil,
		 aubdtil2,irowbdtil2,jqbdtil2,
		 irowdd, jqdd, audd,
		 irowddtil2, jqddtil2, auddtil2,
		 irowddinv, jqddinv, auddinv,
		 Bd,irowb,jqb,
		 Dd,irowd,jqd,
		 Ddtil,irowdtil,jqdtil,
		 neq,b,bhat,islavnode,imastnode,nslavnode,nmastnode,
		 islavact,islavactdof,
		 gap,
		 slavnor,slavtan,
		 vold,vini,cstress,cstressini,
		 bp_old,nactdof,ntie,mi,nk,
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
		 islavnodeinv,lambdaiwan,lambdaiwanini,&k,nmethod,bet,
		 ithermal,plkcon,nplkcon);
  }
 
  fin= clock();	
  printf("\tmultimortar : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);
  number=10;	
  nzs[0]=jq[neq[1]]-1;  
  nzs[1]=jq[neq[1]]-1;
  
  fflush(stdout); 
  debug=1;
  
  /* calculating icol and icolc (needed for SPOOLES) */   
  for(i=0; i<neq[1]; i++){	
    icol[i] = jq[i+1]-jq[i];
  } 
  
  /* nzlc is the number of the rightmost column with 
     nonzero off-diagonal terms */    
  number=10;    
  
  *irowp = irow; *aup=au; 
  *aucp=auc; *irowcp=irowc;
  *aubdp=aubd; *irowbdp=irowbd;
  *aubdtilp=aubdtil; *irowbdtilp=irowbdtil;
  *aubdtil2p=aubdtil2; *irowbdtil2p=irowbdtil2;
  *auddp=audd; *irowddp=irowdd;
  *auddinvp=auddinv; *irowddinvp=irowddinv;    
  *auddtilp=auddtil; *irowddtilp=irowddtil;
  *auddtil2p=auddtil2; *irowddtil2p=irowddtil2;
  *Bdp = Bd; *irowbp= irowb;
  *Bdhelpp = Bdhelp; *irowbhelpp= irowbhelp;
  *Ddp = Dd; *irowdp= irowd;
  *Ddtilp = Ddtil; *irowdtilp= irowdtil;
  *Bdtilp = Bdtil; *irowbtilp= irowbtil;
  *Bpgdp = Bpgd; *irowbpgp= irowbpg;
  *Dpgdp = Dpgd; *irowdpgp= irowdpg;
  *Dpgdtilp = Dpgdtil; *irowdpgtilp= irowdpgtil;
  *Bpgdtilp = Bpgdtil; *irowbpgtilp= irowbpgtil;
  
  SFREE(rs);SFREE(rsb);
  fflush(stdout); 
  printf(" contactmortar: end\n\n"); 
  //FORTRAN(stop,());
  return;
}
