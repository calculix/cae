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
#ifdef SPOOLES
   #include "spooles.h"
#endif
#ifdef SGI
   #include "sgi.h"
#endif
#ifdef TAUCS
   #include "tau.h"
#endif
#ifdef PARDISO
   #include "pardiso.h"
#endif

/**
 * \brief function called before solver transforming the SPCs/MPCs, the matrix and the right hand side for quadratic elements (see Phd-thesis Sitzmann,Chapter 4)
 * Author: Saskia Sitzmann
 *
 * @param [out] iflagact	flag indicating, whether the coupling matrices have to be recalculated
 * @param [in]	ismallsliding	flag indicating, whether coupling matrices are callculated in every iteration (==0) or only in the first iteration (==1)
 * @param [in]  nzs		size of stiffness matrix K
 * @param [out] nzsc2		number of nonzero,nondiagonal entries of intermediate system matrix
 * @param [out] auc2p		intermediate system matrix
 * @param [out] adc2p           intermediate system matrix, diagonal terms
 * @param [out] irowc2p         rows for intermediate system matrix
 * @param [out] icolc2p		columns for intermediate system matrix
 * @param [out] jqc2p		pointer to irowc 
 * @param [out] aubdp		coupling matrix \f$ B_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms 
 * @param [out] irowbdp		field containing row numbers of aubd
 * @param [out] jqbdp		pointer into field irowbd
 * @param [out] aubdtilp	matrix \f$ \tilde{D}^{-1}\tilde{B}_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [out] irowbdtilp	field containing row numbers of aubd
 * @param [out] jqbdtilp	pointer into field irowbdtil 
 * @param [out] aubdtil2p	coupling matrix \f$ \tilde{D}$ and $\tilde{B}^2_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms 
 * @param [out] irowbdtil2p	field containing row numbers of aubdtil2
 * @param [out] jqbdtil2p	pointer into field irowbdtil2
 * @param [out] auddp		coupling matrix \f$ D_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [out] irowddp		field containing row numbers of audd
 * @param [out] jqddp		pointer into field irowdd
 * @param [out] auddtilp	coupling matrix \f$ \tilde{D}_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [out] irowddtilp	field containing row numbers of audd
 * @param [out] jqddtilp	pointer into field irowdd
 * @param [out] auddtil2p	matrix \f$ Id_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [out] irowddtil2p	field containing row numbers of audd
 * @param [out] jqddtil2p	pointer into field irowdd
 * @param [out] auddinvp	coupling matrix \f$ \tilde{D}^{-1}_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
 * @param [out] irowddinvp	field containing row numbers of auddinv
 * @param [out] jqddinvp	pointer into field irowddinv
 * @param [out] jqtempp		field storing the untransformed stiffness matrix representation
 * @param [out] irowtempp	field storing the untransformed stiffness matrix representation
 * @param [out] icoltempp	field storing the untransformed stiffness matrix representation
 * @param [out] nzstemp		field storing the untransformed stiffness matrix size
 * @param [in] iit		iteration number of Newton-Raphson iteration 
 * @param [out] slavnor		slave normal
 * @param [out] slavtan		slave tangent 
 * @param [in,out] icol		stiffness matrix representation, columns
 * @param [in,out] irow		stiffness matrix representation, rows
 * @param [in,out] jq	        stiffness matrix representation, column ointer to irow
 * @param [in] ikboun           sorted dofs idof=8*(node-1)+dir for SPCs
 * @param [in] ilboun           SPC numbers for sorted dofs
 * @param [in] ikmpc 		sorted dofs idof=8*(node-1)+dir for MPCs
 * @param [in] ilmpc		SPC numbers for sorted dofs 
 * @param [out] nboun2          number of transformed SPCs
 * @param [out] ndirboun2p	(i) direction of transformed SPC i 
 * @param [out] nodeboun2p      (i) node of transformed SPC i
 * @param [out] xboun2p         (i) value of transformed SPC i
 * @param [out] nmpc2		number of transformed mpcs
 * @param [out] ipompc2p        (i) pointer to nodempc and coeffmpc for transformed MPC i
 * @param [out] nodempc2p       i and directions of transformed MPCs
 * @param [out] coefmpc2p       coefficients of transformed MPCs
 * @param [out] labmpc2p 	transformed mpc labels
 * @param [out] ikboun2p        sorted dofs idof=8*(node-1)+dir for transformed SPCs
 * @param [out] ilboun2p        transformed SPC numbers for sorted dofs
 * @param [out] ikmpc2p 	sorted dofs idof=8*(node-1)+dir for transformed MPCs
 * @param [out] ilmpc2p		transformed SPC numbers for sorted dofs
 * @param [out] nslavspcp	(2*i) pointer to islavspc...
 * @param [out] islavspcp       ... which stores SPCs for slave node i
 * @param [out] nsspc            number of SPC for slave i
 * @param [out] nslavmpcp	(2*i) pointer to islavmpc...
 * @param [out] islavmpcp	... which stores MPCs for slave node i
 * @param [out] nsmpc		number of MPC for slave i
 * @param [out] nslavspc2p	(2*i) pointer to islavspc2...
 * @param [out] islavspc2p      ... which stores transformed SPCs for slave node i
 * @param [out] nsspc2          number of transformed SPC for slave i
 * @param [out] nslavmpc2p	(2*i) pointer to islavmpc2...
 * @param [out] islavmpc2p	... which stores transformed MPCs for slave node i
 * @param [out] nsmpc2		number of transformed MPC for slave i 
 * @param [in] imastnode	field storing the i of the master surfaces
 * @param [in] nmastnode	(i)pointer into field imastnode for contact tie i 
 * @param [out] nmastspcp	(2*i) pointer to imastspc...
 * @param [out] imastspcp        ... which stores SPCs for master node i
 * @param [out] nmspc            number of SPC for master i
 * @param [out] nmastmpcp	(2*i) pointer to imastmpc...
 * @param [out] imastmpcp	... which stores MPCs for master node i
 * @param [out] nmmpc		number of MPC for master i 
 * @param [in] co		coordinates of i 
 * @param [in] nk		number of i
 * @param [in] kon 		.. for element i storing the connectivity list of elem. in succ. order 
 * @param [in] ipkon		pointer into field kon... 
 * @param [in] lakon		(i) label for element i
 * @param [in] ne		number of elements
 * @param [in] stn		?
 * @param [in] elcon		material parameters
 * @param [in] nelcon           (1,i) number of elastic constants for material i (2,i) number of temperature points 
 * @param [in] rhcon		(0,j,i) temperature (1,j,i) density at the density temperature point j of material i
 * @param [in] nrhcon		(i) number of temperature data points for the density material i
 * @param [in] alcon		(0,j,i) temperature, (k,j,i) expansion coefficient k at expansion temperature point j of material i
 * @param [in] nalcon		(1,i) number of expansion constants  (2,i) number of temperature data points for expansion coefficients for material i
 * @param [in] alzero		used in material data
 * @param [in] ielmat           (j,i) material number of layer j
 * @param [in] ielorien		(j,i) orientation number of layer j
 * @param [in] norien		number of orientations
 * @param [in] orab		(7,*) description of local coordinate system
 * @param [in] ntmat_		maximum number of temperature data points for any material
 * @param [in] t0		needed for spring elements
 * @param [in] t1		needed for spring elements
 * @param [in] ithermal		>0 thermal effects are taken into account 
 * @param [in] prestr		prestress values /plastic inital strain ???
 * @param [in] iprestr		flag indicating whether prestress is defined 
 * @param [in] filab		?
 * @param [in] eme		mechanical strain for element quadrature point-wise
 * @param [in] emn		?
 * @param [in] een		?
 * @param [in] iperturb		geometrical method
 * @param [in] f		internal forces
 * @param [in] nactdof 		(i,j) actual degree of freedom for direction i of node j
 * @param [in] iout		flag indicating what to calculate in results
 * @param [in] qa
 * @param [in] vold		displacements in i
 * @param [in] b		right hand side
 * @param [in] nboun            number of SPCs
 * @param [in] ndirboun		(i) direction of SPC i 
 * @param [in] nodeboun         (i) node of SPC i
 * @param [in] xbounact         (i) current value of SPC i 
 * @param [in] xboun            (i) value of SPC i
 * @param [in] nmpc		number of mpcs
 * @param [in] ipompc           (i) pointer to nodempc and coeffmpc for MPC i
 * @param [in] nodempc          i and directions of MPCs
 * @param [in] labmpc		MPC labels
 * @param [in] coefmpc          coefficients of MPCs
 * @param [in] nmethod		analysis method 
 * @param [in] 	neq             number of active degrees of freedom
 * @param [in] veold		velocities in i
 * @param [in] accold 		accelerations in i
 * @param [in] dtime		real time increment size
 * @param [in] time		real time size of all previous increments including the present increment 
 * @param [in] ttime		real time size of all previous steps
 * @param [in] plicon		isotropic hardening curve 
 * @param [in] nplicon          pointer to isotropic hardening curve 
 * @param [in] plkcon		kinematic hardening curve
 * @param [in] nplkcon		pointer to kinematik hardening curve
 * @param [in] xstateini	state variables at start if the increment
 * @param [in] xstiff		(1:27,k,i) stiffness matrix of element i in integration point k for last increment
 * @param [in] xstate		current state variables
 * @param [in] npmat_		maximum number of data points for plicon
 * @param [in] matname		(i) name of material i
 * @param [in] mi		(1) max # of integration points per element (2) max degree of freedom per element
 * @param [in] ielas		used in mechmodel   
 * @param [in] icmd 		used in mechmodel  
 * @param [in] ncmat_		maximum number of elastic material constants 
 * @param [in] nstate_		maximum number of state variables
 * @param [in] stiini		(1:6,k,i) 2nd order stress of element i in integration point k at start of the increment
 * @param [in] vini		displacements at the start of the increment
 * @param [in] ener		internal energy and kinetic energy
 * @param [in] enern		?
 * @param [in] emeini		?
 * @param [in] xstaten		?
 * @param [in] eei		mechanical strain ?
 * @param [in] enerini		internal engery at the start od the increment
 * @param [in] cocon		used for thermal calculation, @see materialdata_th
 * @param [in] ncocon		used for thermal calculation, @see materialdata_th 
 * @param [in] set		(i) name of set i 
 * @param [in] nset		number of sets 
 * @param [in] istartset        (i) pointer to ialset containing the first set member
 * @param [in] iendset          (i) pointer to ialset containing the last set member
 * @param [in] ialset           set members
 * @param [in] nprint		?
 * @param [in] prlab		?
 * @param [in] prset		?
 * @param [in] qfx		qflux in integration points
 * @param [in] qfn		?
 * @param [in] trab		?
 * @param [in] inotr		?
 * @param [in] ntrans		? 
 * @param [in] nelemload	element faces of which are loaded
 * @param [in] nload		number of facial distributed loads 
 * @param [in] istep		step number
 * @param [in] iinc		increment number
 * @param [in] springarea	used in spring elements
 * @param [in] reltime		theta+dtheta
 * @param [in] ne0		?
 * @param [in] xforc		pointforce, value
 * @param [in] nforc		number of point forces
 * @param [in] thicke		thickness matrix
 * @param [in] shcon		used for thermal calculation, @see materialdata_th
 * @param [in] nshcon		used for thermal calculation, @see materialdata_th
 * @param [in] sideload		load label
 * @param [in] xload		concentrated load 
 * @param [in] xloadold		concentrated load old 
 * @param [in] icfd		?
 * @param [in] inomat		?
 * @param [in] islavelinv       (i)==0 if there is no slave node in the element, >0 otherwise
 * @param [in] islavsurf	islavsurf(1,i) slaveface i islavsurf(2,i) # integration points generated before looking at face i 
 * @param [in] iponoels		?
 * @param [in] inoels		?
 * @param [out] fexttil		transformed external forces
 * @param [out] ftil		transformed internal forces
 * @param [in]  fextinitil	transformed external forces at start of the increment
 * @param [in]  finitil		transformed internal forces at start of the increment
 * @param [in] mortar		param indicating what contact formulation is used (=0 NTS penalty, =1 GPTS penalty , >1 STS mortar)
 * @param [in] nslavnode	(i)pointer into field isalvnode for contact tie i 
 * @param [in] islavnode	field storing the i of the slave surface 
 * @param [in] nslavs		number of slave i
 * @param [in] ntie		number of ties
 * @param [in] autloc		transformation matrix \f$ T[p,q]\f$ for slave i \f$ p,q \f$  
 * @param [in] irowtloc		field containing row numbers of autloc
 * @param [in] jqtloc	        pointer into field irowtloc
 * @param [in] autlocinv	transformation matrix \f$ T^{-1}[p,q]\f$ for slave i \f$ p,q \f$ 
 * @param [in] irowtlocinv	field containing row numbers of autlocinv
 * @param [in] jqtlocinv	pointer into field irowtlocinv
 * @param [in] nk2		number or generated points needed for transformed SPCs 
 * @param [in]  iflagdualquad   flag indicating what mortar contact is used (=1 quad-lin, =2 quad-quad, =3 PG quad-lin, =4 PG quad-quad) 
 * @param [in] tieset           (1,i) name of tie constraint (2,i) dependent surface (3,i) independent surface 
 * @param [in] itiefac 		pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i
 * @param [in] rhsi		flag indicating whether to calculate the right hand side 
 * @param [in,out] au		stiffness matrix K, non-diagonal terms
 * @param [in,out] ad		stiffness matrix K, dialgonal terms 
 * @param [out] f_cmp		not used any more
 * @param [out] f_csp		contact force for active degrees of freedom
 * @param [in] t1act
 * @param [in] cam
 * @param [in] bet		parameter for alpha method
 * @param [in] gam		parameter for alpha method
 * @param [in] epn
 * @param [in] xloadact		magnitude of load
 * @param [in] nodeforc		point force, node
 * @param [in] ndirforc		point force, dir
 * @param [in] xforcact		pointforce, value
 * @param [in] xbodyact		body forces
 * @param [in] ipobody		pointer to xbody...	
 * @param [in] nbody		number of mechanical body loads
 * @param [in] cgr		gravity forces
 * @param [out] nzl		highes column numbeer with non-zero entry
 * @param [in] sti		(1:6,k,i) 2nd order stress of element i in integration point k 
 * @param [in] iexpl		flag indicating whether explicit (=1) or implicit (=0) method is used 
 * @param [in] mass		flag indicating whether to calculate the mass matrix
 * @param [in] buckling		flag indicating whether to calculate the buckling stiffness
 * @param [in] stiffness	flag indicating whether to calculate the stiffness matrix
 * @param [in] intscheme	flag indicating what integration scheme to use
 * @param [in] physcon		used for thermal calculation
 * @param [in] coriolis		flag indicating whether to calculate the coriolis matrix
 * @param [in] ibody 		used for assinging body forces
 * @param [in] integerglob	used for submodel
 * @param [in] doubleglob	used for submodel
 * @param [in] nasym		flag indicating whether matrix is symmetric or not
 * @param [in] alpham		parameter for damping
 * @param [in]	betam		parameter for damping
 * @param [in]	auxtil2		auxilary field
 * @param [in] pslavsurf	field storing  position xil, etal and weight for integration point on slave side
 * @param [in] pmastsurf 	field storing position and etal for integration points on master side
 * @param [in] clearini		used for spring elements
 * @param [in] ielprop		used for beam elements
 * @param [in] prop		used for beam elements
 * @param [in] islavact		(i) indicates, if slave node i is active (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, =0 inactive node, =1 sticky node, =2 slipping/active node) 
 * @param [in] cdn		?
 * @param [in] memmpc_		size of nodempc, coeffmpc
 * @param [in] cvtilini		C*v at start of the increment
 * @param [in] cvtil		C*v
 * @param [in] idamping		flag indicating whether damping is used
 * @param [in] ilin		flag indicating wheter fist iteration is calculated linear geometrically 
 * @param [in] iperturb_sav	saved iperturb values 
 * @param [out] nodeforc2p	transformed point force, node
 * @param [out] ndirforc2p	transformed point force, dir
 * @param [out] xforc2p		transformed point force, value
 * @param [out] nforc2		number of transformed point forces   
**/
void precontact_mortar(ITG *iflagact,ITG *ismallsliding,ITG *nzs,ITG *nzsc2,
                       double **auc2p,double **adc2p,ITG **irowc2p, ITG **icolc2p,ITG **jqc2p,
                       double **aubdp,ITG **irowbdp,ITG **jqbdp,
                       double **aubdtilp,ITG **irowbdtilp,ITG **jqbdtilp,
                       double **aubdtil2p,ITG **irowbdtil2p,ITG **jqbdtil2p,
                       double **auddp,ITG **irowddp,ITG **jqddp,
                       double **auddtilp,ITG **irowddtilp,ITG **jqddtilp,
                       double **auddtil2p,ITG **irowddtil2p,ITG **jqddtil2p,
                       double **auddinvp,ITG **irowddinvp,ITG **jqddinvp,
                       ITG **jqtempp,ITG **irowtempp, ITG **icoltempp,ITG *nzstemp,
                       ITG *iit,double *slavnor, double *slavtan,
                       ITG *icol, ITG *irow, ITG *jq,
                       ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
                       ITG *nboun2,ITG **ndirboun2p,ITG **nodeboun2p,double **xboun2p,
                       ITG *nmpc2,ITG **ipompc2p, ITG **nodempc2p,double **coefmpc2p,char **labmpc2p,
                       ITG **ikboun2p,ITG **ilboun2p,ITG **ikmpc2p,ITG **ilmpc2p,
                       ITG **nslavspcp,ITG **islavspcp,ITG **nslavmpcp,ITG **islavmpcp,
                       ITG **nslavspc2p,ITG **islavspc2p, ITG **nslavmpc2p,ITG **islavmpc2p,
                       ITG **nmastspcp,ITG **imastspcp,ITG **nmastmpcp,ITG **imastmpcp,
                       ITG **nmastmpc2p,ITG **imastmpc2p,ITG *nmmpc2,
                       ITG *nsspc, ITG *nsspc2, ITG *nsmpc, ITG *nsmpc2,
                       ITG *imastnode,ITG *nmastnode,ITG *nmspc, ITG *nmmpc,		       
                       double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
                       ITG *ne,double *stn, 
                       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
                       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
                       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
                       double *t0,double *t1,ITG *ithermal,double *prestr, 
                       ITG *iprestr,char *filab,double *eme,double *emn,
                       double *een,ITG *iperturb,double *f,ITG *nactdof,
                       ITG *iout,double *qa,
                       double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
                       double *xbounact,double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,
                       double *coefmpc,char *labmpc,ITG *nmpc,ITG *nmethod, 
                       ITG *neq,double *veold,double *accold,
                       double *dtime,double *time,
                       double *ttime,double *plicon,
                       ITG *nplicon,double *plkcon,ITG *nplkcon,
                       double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
                       char *matname,ITG *mi,ITG *ielas,
                       ITG *icmd,ITG *ncmat_,ITG *nstate_,double *stiini,
                       double *vini,double *ener,
                       double *enern,double *emeini,double *xstaten,double *eei,
                       double *enerini,double *cocon,ITG *ncocon,char *set, 
                       ITG *nset,ITG *istartset,
                       ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
                       char *prset,double *qfx,double *qfn,double *trab,
                       ITG *inotr,ITG *ntrans,ITG *nelemload,
                       ITG *nload,ITG *istep,ITG *iinc,
                       double *springarea,double *reltime,ITG *ne0,double *xforc,
                       ITG *nforc,double *thicke,
                       double *shcon,ITG *nshcon,char *sideload,double *xload,
                       double *xloadold,ITG *icfd,ITG *inomat,
                       ITG *islavelinv,ITG *islavsurf,
                       ITG *iponoels,ITG *inoels,                       
                       double *fexttil, double *ftil,double *fextinitil, double *finitil,
                       ITG *mortar, ITG *nslavnode,ITG *islavnode,ITG *nslavs, ITG *ntie,
                       double *autloc,ITG *irowtloc,ITG *jqtloc,
                       double *autlocinv,ITG *irowtlocinv,ITG *jqtlocinv,
                       ITG *nk2, ITG *iflagdualquad,
                       char *tieset,ITG *itiefac  ,ITG *rhsi,
                       double *au,double *ad, double **f_cmp, double **f_csp,
                       double *t1act,double *cam,double *bet, double *gam,double *epn,
                       double *xloadact,ITG *nodeforc,ITG *ndirforc,double *xforcact,
                       double *xbodyact,ITG *ipobody, ITG *nbody,double *cgr,
                       ITG *nzl,double *sti,ITG *iexpl,ITG *mass,ITG *buckling,ITG *stiffness,
                       ITG *intscheme,double *physcon,ITG *coriolis,ITG *ibody,
                       ITG *integerglob,double *doubleglob,ITG *nasym,
                       double *alpham,double *betam,double *auxtil2,
                       double *pslavsurf,double *pmastsurf,
                       double *clearini,ITG *ielprop,double *prop,
                       ITG *islavact,double *cdn,ITG *memmpc_,
                       double *cvinitil,double *cvtil,ITG *idamping,
                       ITG *ilin, ITG *iperturb_sav,
                       double *adbtil, double *aubtil,double *adb, double *aub,
                       ITG **nodeforc2p,ITG **ndirforc2p,double **xforc2p,ITG *nforc2,
                       ITG *itietri,double *cg, double *straight,ITG *koncont,
		       double *energyini,
       			double *energy,ITG *kscale,ITG *iponoel,ITG *inoel,ITG *nener,
       			char *orname,ITG *network,
		       char *typeboun,ITG *num_cpus
                       ){
  
  ITG im,i,ii,j,jj,k,l,mt=mi[1]+1,node,idof1,idof2,jfaces,nelems,ifaces,nope,nopes,idummy,
    *ndirboun2=NULL,*nodeboun2=NULL,*ipompc2=NULL,*nodempc2=NULL,nodes[8],konl[20],jj2,ifac,
    *ikboun2=NULL,*ilboun2=NULL,*ikmpc2=NULL,*ilmpc2=NULL,
    *nslavspc=NULL,*islavspc=NULL,*nslavmpc=NULL,*islavmpc=NULL,
    *nslavspc2=NULL,*islavspc2=NULL,*nslavmpc2=NULL,*islavmpc2=NULL,
    *nmastspc=NULL,*imastspc=NULL,*nmastmpc=NULL,*imastmpc2=NULL,*nmastmpc2=NULL,*imastmpc=NULL,
    *irowc2=NULL,*icolc2=NULL,*jqc2=NULL,*irowbd=NULL,*jqbd=NULL,*irowbdtil=NULL,*jqbdtil=NULL,
    *irowbdtil2=NULL,*jqbdtil2=NULL,*irowdd=NULL,*jqdd=NULL,*irowddtil=NULL,*jqddtil=NULL,
    *irowddinv=NULL,*jqddinv=NULL,*irowddtil2=NULL,*jqddtil2=NULL,
    *irowtemp=NULL, *icoltemp=NULL,*jqtemp=NULL,*inum=NULL,*icoltil=NULL,
    *irowtil=NULL,*jqtil=NULL,*irowbtil=NULL,*jqbtil=NULL,
    *nodeforc2=NULL,*ndirforc2=NULL;
  
  double scal1,alpha,*xboun2=NULL,*coefmpc2=NULL,*auc2=NULL,*adc2=NULL,*aubd=NULL,*aubdtil=NULL,*aubdtil2=NULL,
    *audd=NULL,*auddtil=NULL,*auddinv=NULL, *auddtil2=NULL,*v=NULL,*stx=NULL,*fn=NULL,*autil=NULL,
    *adtil=NULL,*fmpc2=NULL,*voldtil=NULL,*vinitil=NULL,*accoldtil=NULL,*btil=NULL,
    *adctil=NULL,*auctil=NULL,*veoldtil=NULL,*volddummy=NULL,*vectornull=NULL,*f_cs=NULL,*f_cm=NULL,
    *xforc2=NULL, *fnext=NULL,*veolddummy=NULL,*accolddummy=NULL;
  
  char *labmpc2=NULL;
  
  clock_t debut;
  clock_t fin; 
    
  alpha=1-2*sqrt(*bet);
  
  ndirboun2=*ndirboun2p;nodeboun2=*nodeboun2p;xboun2=*xboun2p;
  ipompc2=*ipompc2p;nodempc2=*nodempc2p;coefmpc2=*coefmpc2p;
  labmpc2=*labmpc2p;ikboun2=*ikboun2p;ilboun2=*ilboun2p;ikmpc2=*ikmpc2p;ilmpc2=*ilmpc2p;
  nslavspc=*nslavspcp;islavspc=*islavspcp;nslavmpc=*nslavmpcp;islavmpc=*islavmpcp;
  nslavspc2=*nslavspc2p;islavspc2=*islavspc2p;nslavmpc2=*nslavmpc2p;islavmpc2=*islavmpc2p;
  nmastspc=*nmastspcp;imastspc=*imastspcp;nmastmpc=*nmastmpcp;imastmpc=*imastmpcp;nmastmpc2=*nmastmpc2p;imastmpc2=*imastmpc2p;
  auc2=*auc2p;adc2=*adc2p;irowc2=*irowc2p;icolc2=*icolc2p;jqc2=*jqc2p;
  aubd=*aubdp;irowbd=*irowbdp;jqbd=*jqbdp;
  aubdtil=*aubdtilp;irowbdtil=*irowbdtilp;jqbdtil=*jqbdtilp;
  aubdtil2=*aubdtil2p;irowbdtil2=*irowbdtil2p;jqbdtil2=*jqbdtil2p;
  audd=*auddp;irowdd=*irowddp;jqdd=*jqddp;
  auddtil=*auddtilp;irowddtil=*irowddtilp;jqddtil=*jqddtilp;
  auddtil2=*auddtil2p;irowddtil2=*irowddtil2p;jqddtil2=*jqddtil2p;
  auddinv=*auddinvp;irowddinv=*irowddinvp;jqddinv=*jqddinvp;
  irowtemp=*irowtempp;icoltemp=*icoltempp;jqtemp=*jqtempp;
  f_cs=*f_csp;f_cm=*f_cmp;
  ndirforc2=*ndirforc2p;nodeforc2=*nodeforc2p;xforc2=*xforc2p;
  
  fflush(stdout);
  
  NNEW(f_cs,double,neq[1]);
  NNEW(f_cm,double,neq[1]);
  
  // check for coupled thermo-mechanical calculation
  if(*ithermal>1){
    printf("\tprecontactmortar: coupled thermo-mechanical calculations NOT supported yet!\n \tPlease use surface-to-surface penalty contact instead.\n\n STOP!\n");
    fflush(stdout);
    //FORTRAN(stop,());
  }
  
  // fix for linear calculation in first iteration of first increment
  if((*nmethod!=4)&&(*ilin==1)&&(*iit==1)&&(*iinc==1)){  
    *ielas=1;  
    iperturb[0]=-1;  
    iperturb[1]=0;	  
  }  
  
  /* small sliding is autometically set active due to combined fix-point Newton approach 
     do NOT change this unless the additional derivates neglected here have been implemented*/
  *iflagact=0;
  *ismallsliding=1;
  *nzsc2=nzs[1];
  NNEW(auc2,double,*nzsc2);
  NNEW(adc2,double,neq[1]);
  NNEW(irowc2,ITG,*nzsc2);
  NNEW(icolc2,ITG,neq[1]);
  NNEW(jqc2,ITG,neq[1]+1); 
  if(*iit==1 || *ismallsliding==0){
    NNEW(aubd,double, 6*nslavnode[*ntie]);
    NNEW(irowbd,ITG, 6*nslavnode[*ntie]);
    NNEW(jqbd,ITG, neq[1]+1);
    NNEW(aubdtil,double, 6*nslavnode[*ntie]);
    NNEW(irowbdtil,ITG, 6*nslavnode[*ntie]);
    NNEW(jqbdtil,ITG, neq[1]+1);
    NNEW(aubdtil2,double, 6*nslavnode[*ntie]);
    NNEW(irowbdtil2,ITG, 6*nslavnode[*ntie]);
    NNEW(jqbdtil2,ITG, neq[1]+1);
    NNEW(audd,double, 3*nslavnode[*ntie]);
    NNEW(irowdd,ITG, 3*nslavnode[*ntie]);
    NNEW(jqdd,ITG, neq[1]+1);
    NNEW(auddtil,double, 3*nslavnode[*ntie]);
    NNEW(irowddtil,ITG, 3*nslavnode[*ntie]);
    NNEW(jqddtil,ITG, neq[1]+1);
    NNEW(auddtil2,double, 3*nslavnode[*ntie]);
    NNEW(irowddtil2,ITG, 3*nslavnode[*ntie]);
    NNEW(jqddtil2,ITG, neq[1]+1);
    NNEW(auddinv,double,3*nslavnode[*ntie]);
    NNEW(irowddinv,ITG, 3*nslavnode[*ntie]);
    NNEW(jqddinv,ITG, neq[1]+1);	    
  }
  
  //printf("vold(%" ITGFORMAT ")= %f %f %f \n",103298,vold[mt*(103298)-3+0],vold[mt*(103298)-3+1],vold[mt*(103298)-3+2]);
  
  /* storing the original stiffness matrix */
  NNEW(jqtemp,ITG,neq[1]+1);
  NNEW(irowtemp,ITG,nzs[1]);
  NNEW(icoltemp,ITG,neq[1]);
  for(i=0;i<3;i++){nzstemp[i]=nzs[i];}
  for (i=0;i<neq[1];i++){jqtemp[i]=jq[i];icoltemp[i]=icol[i];}
  jqtemp[neq[1]]=jq[neq[1]];
  for (i=0;i<nzs[1];i++){irowtemp[i]=irow[i];}
  
  if(*iit==1 || *ismallsliding==0){
    DMEMSET(slavnor,0,3*nslavnode[*ntie],0.);
    DMEMSET(slavtan,0,6*nslavnode[*ntie],0.);
  }
  
  // setting iflagact=1 before calling contactmortar invokes combined fix-point Newton approach
  if(*iit>1 && *ismallsliding==1){*iflagact=1;}
  
  /* transform SPCs/MPCs in case of quadratic finite elements
     Caution: there is still a problem with MPCs on slave mid nodes, avoid this!!! */
  //wird immer aufgerufen, da sich xbounact geaendert haben kann
  transformspcsmpcs_quad(nboun,ndirboun,nodeboun,xbounact,
			 nmpc,ipompc,nodempc,coefmpc,labmpc,
			 ikboun,ilboun,ikmpc,ilmpc,
			 nboun2,&ndirboun2,&nodeboun2,&xboun2,
			 nmpc2,&ipompc2,&nodempc2,&coefmpc2,&labmpc2,
			 &ikboun2,&ilboun2,&ikmpc2,&ilmpc2,
			 irowtlocinv,jqtlocinv,autlocinv, 
			 nk,nk2,iflagdualquad,
			 ntie,tieset,itiefac,islavsurf,
			 lakon,ipkon,kon,&mt,memmpc_,
			 nodeforc,ndirforc,xforcact,nforc,
			 &nodeforc2,&ndirforc2,&xforc2,nforc2);
  //SFREE(islavspc);SFREE(islavmpc);SFREE(islavspc2);SFREE(islavmpc2);SFREE(imastspc);SFREE(imastmpc);
  RENEW(islavspc,ITG,2**nboun);
  RENEW(islavmpc,ITG,2**nmpc);
  RENEW(islavspc2,ITG,2**nboun2);
  RENEW(islavmpc2,ITG,2**nmpc2);
  RENEW(imastspc,ITG,2**nboun);
  RENEW(imastmpc,ITG,2**nmpc);
  RENEW(imastmpc2,ITG,2**nmpc);
  /* cataloque SPCs/MPCs */
  FORTRAN(conttiemortar,(lakon,ipkon,kon,ntie,tieset,nset,set,
			 itiefac,islavsurf,islavnode,imastnode,nslavnode,nmastnode,
			 nslavs,iponoels,inoels,nk,
			 nboun,ndirboun,nodeboun,xbounact,
			 nmpc,ipompc,nodempc,coefmpc,
			 ikboun,ilboun,ikmpc,ilmpc,
			 nboun2,ndirboun2,nodeboun2,xboun2,
			 nmpc2,ipompc2,nodempc2,coefmpc2,
			 ikboun2,ilboun2,ikmpc2,ilmpc2,
			 nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
			 nslavspc2,islavspc2,nsspc2,nslavmpc2,islavmpc2,nsmpc2,
			 nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
			 nmastmpc2,imastmpc2,nmmpc2));       
  RENEW(islavspc,ITG,2**nsspc+1);
  RENEW(islavmpc,ITG,2**nsmpc+1);
  RENEW(islavspc2,ITG,2**nsspc2+1);
  RENEW(islavmpc2,ITG,2**nsmpc2+1);
  RENEW(imastspc,ITG,2**nmspc+1);
  RENEW(imastmpc,ITG,2**nmmpc+1);
  RENEW(imastmpc2,ITG,2**nmmpc2+1);
  
  /* Check for additional MPCs on slave mid nodes */
    if(*iit==1 && *iinc==1){
      
    FORTRAN(gencontrel,(tieset,ntie,itietri,ipkon,kon,
			lakon,set,cg,straight,
			koncont,co,vold,nset,
			islavsurf,itiefac,
			islavnode,nslavnode,slavnor,slavtan,mi));      
      // call checkspsmpc
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
			 
        
      for (i=0;i<*ntie;i++){  
	if(tieset[i*(81*3)+80]=='C'){      
	  for(j=nslavnode[i];j<nslavnode[i+1];j++){     
	    node=islavnode[j];
	    int checkformidnode=0;
	    for(jj=nslavmpc[2*(j)];jj<nslavmpc[2*(j)+1];jj++){
//	      printf("%"ITGFORMAT" %"ITGFORMAT" \n",islavmpc[2*jj],islavmpc[2*jj+1]);
	      if(islavmpc[2*jj+1]==-2){
		// MPC cannot be identified
		checkformidnode=1;
		printf(" check for mid node %"ITGFORMAT"\n",node);
	      }	      
	    }
     // @todo: set this part active again!
    if(checkformidnode==1){
      //check for mid node
      do{
	    	      
	    for(l=itiefac[2*i];l<=itiefac[2*i+1];l++){
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

	break;
      }while(1);
      if((ii>2 && nopes==6)||(ii>3 && nopes==8)){
	// mid node found with extra not supported mpc ->error
	printf(" precontactmortar: Problem with slave mid node  \n\n");
	printf(" *ERROR: Slave mid node %"ITGFORMAT" has additional MPC which is not a directional blocking MPC or a 1-to-1 cyclic symmetry MPC.  \n",node);
	printf("\t\t This is not supported yet!!!!!!!!!!!!!\n");
	fflush(stdout);
	FORTRAN(stop,());
      }
    }        
	  }
	}
      }
  if(*iit==1 || *ismallsliding==0){
    DMEMSET(slavnor,0,3*nslavnode[*ntie],0.);
    DMEMSET(slavtan,0,6*nslavnode[*ntie],0.);
  } 
    printf(" precontactmortar: MPC-check OK\n\n");
    }
  /* fix for quadratic FE */
  debut=clock();
  NNEW(v,double,mt**nk);
  NNEW(stx,double,6*mi[0]**ne);
  NNEW(fn,double,mt*(*nk+*nk2));
  NNEW(fmpc2,double,*nmpc2);
  //lsprintf("nmpc2 %"ITGFORMAT" \n",*nmpc2);
  //printf("vold(%" ITGFORMAT ")= %f %f %f \n",103298,vold[mt*(103298)-3+0],vold[mt*(103298)-3+1],vold[mt*(103298)-3+2]);
  memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
  *iout=-1;
  for(i=0;i<neq[1];i++){ftil[i]=0.0;fexttil[i]=0.0;}
  NNEW(inum,ITG,*nk);
  
  //printf("vold(%" ITGFORMAT ")= %f %f %f \n",103298,vold[mt*(103298)-3+0],vold[mt*(103298)-3+1],vold[mt*(103298)-3+2]);
  /// @TODO: for coupled calculations calling results_dstil causes strange behaviour, debug this for quadratic elements
  if(ithermal[0]==3){

  }else{
  printf(" precontactmortar: call results_dstil\n");
  results_dstil(co,nk,kon,ipkon,lakon,ne,
  		v,stn,inum,stx,elcon,nelcon,
		rhcon,nrhcon,alcon,nalcon,alzero,
		ielmat,ielorien,norien,orab,ntmat_,
		t0,
		t1act,ithermal,prestr,iprestr,filab,
		eme,emn,
		een,iperturb,ftil,fn,nactdof,iout,
		qa,vold,b,nodeboun,ndirboun,
		xbounact,nboun,ipompc,nodempc,coefmpc,
		labmpc,nmpc,nmethod,cam,&neq[1],veold,
		accold,bet,gam,dtime,time,ttime,plicon,nplicon,plkcon,
		nplkcon,xstateini,xstiff,xstate,npmat_,
		epn,matname,mi,ielas,icmd,ncmat_,
		nstate_,
		stiini,vini,ikboun,ilboun,ener,
		enern,emeini,xstaten,eei,enerini,
		cocon,ncocon,set,nset,istartset,
		iendset,
	    	ialset,nprint,prlab,prset,qfx,qfn,
		trab,
		inotr,ntrans,fmpc2,nelemload,nload,
		ikmpc,ilmpc,
		istep,iinc,springarea,reltime,ne0,
		thicke,
		shcon,nshcon,sideload,xloadact,
		xloadold,icfd,inomat,pslavsurf,
		pmastsurf,mortar,islavact,cdn,
		islavnode,nslavnode,ntie,clearini,
        	islavsurf,ielprop,prop,energyini,
       		energy,kscale,iponoel,inoel,nener,
       		orname,network,ipobody,xbodyact,ibody,
       		typeboun,	
		islavelinv,autloc,irowtloc,jqtloc,
		nboun2,ndirboun2,nodeboun2,xboun2,
		nmpc2,ipompc2,nodempc2,coefmpc2,labmpc2,
		ikboun2,ilboun2,ikmpc2,ilmpc2);
		
    printf(" precontactmortar: results_dstil finished\n");
    fflush(stdout);
  }
  
  SFREE(v);SFREE(stx);SFREE(fn);SFREE(inum);SFREE(fmpc2);
  *iout=0;	    
  
  *rhsi=1;
  NNEW(adtil,double,neq[1]);
  NNEW(autil,double,nzs[1]);
  //  for(i=0;i<nzs[1];i++) printf("autil=%d %d %e\n",i,nzs[1],autil[i]);
  NNEW(irowtil,ITG,nzs[1]);
  NNEW(jqtil,ITG,neq[1]+1);
  NNEW(icoltil,ITG,neq[1]);
  //if(*nmethod==4){
  //  NNEW(adbtil,double,neq[1]);
  //  NNEW(aubtil,double,nzs[1]);
    //mass[0] = 1;
  //}
  for(i=0;i<neq[1]+1;i++){jqtil[i]=jq[i];}
  for(i=0;i<neq[1];i++){icoltil[i]=icol[i];}
  for(i=0;i<nzs[1];i++){irowtil[i]=irow[i];}
  printf(" precontactmortar: call mafillsmmain_dstil\n");
  fflush(stdout);
  mafillsmmain_dstil(co,nk,kon,ipkon,lakon,ne,nodeboun2,ndirboun2,xboun2,nboun2,
          ipompc2,nodempc2,coefmpc2,nmpc2,nodeforc2,ndirforc2,xforc2,
          nforc2,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
          nbody,cgr,adtil,autil,fexttil,nactdof,icol,jqtil,irowtil,neq,nzl,
          nmethod,ikmpc2,ilmpc2,ikboun2,ilboun2,
          elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
          ielmat,ielorien,norien,orab,ntmat_,
          t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
          nzs,stx,adbtil,aubtil,iexpl,plicon,nplicon,plkcon,nplkcon,
          xstiff,npmat_,dtime,matname,mi,
          ncmat_,mass,stiffness,buckling,rhsi,intscheme,
          physcon,shcon,nshcon,cocon,ncocon,ttime,time,istep,iinc,
          coriolis,ibody,xloadold,reltime,veold,springarea,nstate_,
          xstateini,xstate,thicke,integerglob,doubleglob,
          tieset,istartset,iendset,ialset,ntie,nasym,pslavsurf,
          pmastsurf,mortar,clearini,ielprop,prop,ne0,fnext,kscale,
	  iponoel,inoel,network,ntrans,inotr,trab,
          nslavnode,islavnode,islavsurf,islavelinv,
	  autloc,irowtloc,jqtloc);  
    printf(" precontactmortar: mafillsmmain_dstil finished\n");
    fflush(stdout);

		
  /* overwrite K matrix with Ktilde */
    /* 
    int number=3;		
    FORTRAN(writematrix,(au,ad,irow,jq,&neq[1],&number));
    
        number=4;		
    FORTRAN(writematrix,(autil,adtil,irow,jq,&neq[1],&number));  
    */
  
  
  for(i=0;i<neq[1];i++){
    ad[i]=adtil[i];
    if(*nmethod==4){
      //adb[i]=adbtil[i];
    }
  }
  for(i=0;i<nzs[1];i++){
    au[i]=autil[i];
    //   printf("autil=%d %d %e\n",i,nzs[1],au[i]);
    if(*nmethod==4){
      //aub[i]=aubtil[i];
    }
  }
  
  
  /* transform vold,vini,accold for dynamic calculations */ 
  NNEW(voldtil,double,mt**nk);
  NNEW(veoldtil,double,mt**nk);
  NNEW(vinitil,double,mt**nk);
  NNEW(accoldtil,double,mt**nk);
  if(*nmethod==4){
    for(i=0;i<*nk;i++){
      node=i+1;
      if(jqtlocinv[node-1+1]-jqtlocinv[node-1]>0){
	for(jj=jqtlocinv[node-1]-1;jj<jqtlocinv[node-1+1]-1;jj++){
	  for(l=0;l<3;l++){
	    idof1=mt*node-3+l;
	    idof2=mt*irowtlocinv[jj]-3+l;
	    accoldtil[idof2]+=autlocinv[jj]*accold[idof1];
	    voldtil[idof2]+=autlocinv[jj]*vold[idof1];
	    veoldtil[idof2]+=autlocinv[jj]*veold[idof1];
	    vinitil[idof2]+=autlocinv[jj]*vini[idof1];
	  }
	  if(ithermal[0]>1){
	    idof1=mt*(node-1);
	    idof2=mt*(irowtlocinv[jj]-1);
	    accoldtil[idof2]+=autlocinv[jj]*accold[idof1];
	    voldtil[idof2]+=autlocinv[jj]*vold[idof1];
	    veoldtil[idof2]+=autlocinv[jj]*veold[idof1];
	    vinitil[idof2]+=autlocinv[jj]*vini[idof1];	    
	  }
	}	
      }else{
	for(l=0;l<3;l++){
	  idof1=mt*node-3+l;
	  accoldtil[idof1]=accold[idof1];
	  voldtil[idof1]=vold[idof1];
	  veoldtil[idof1]=veold[idof1];
	  vinitil[idof1]=vini[idof1];
	}
	if(ithermal[0]>1){
	  idof1=mt*(node-1);
	  accoldtil[idof1]=accold[idof1];
	  voldtil[idof1]=vold[idof1];
	  veoldtil[idof1]=veold[idof1];
	  vinitil[idof1]=vini[idof1];    
	}	
      }
    }	   
  }
  if(((*alpham>0.)||(*betam>0.))&&(*iexpl<=1)){
    NNEW(adctil,double,neq[1]);
    NNEW(auctil,double,nzs[1]);
    for(k=0;k<neq[0];k++){adctil[k]=*alpham*adbtil[k]+*betam*adtil[k];}
    for(k=0;k<nzs[0];k++){auctil[k]=*alpham*aubtil[k]+*betam*autil[k];}
  }
  NNEW(btil,double,neq[1]);
  //if(*nmethod!=4){
  /* overwrite r with rtilde */
  printf(" precontactmortar: call calcresidual\n");
  calcresidual(nmethod,neq,btil,fexttil,ftil,iexpl,nactdof,auxtil2,voldtil,
	       vinitil,dtime,accoldtil,nk,adbtil,aubtil,jqtil,irowtil,nzl,
	       &alpha,fextinitil,finitil,islavnode,nslavnode,mortar,ntie,f_cm,
	       f_cs,mi,nzs,nasym,idamping,veoldtil,adctil,auctil,cvinitil,cvtil,
	       alpham,num_cpus);
  printf(" precontactmortar: calcresidual finished\n");
  fflush(stdout);  
    /*
 	number=3;		
    FORTRAN(writevector,(b, &neq[1],&number));
 	number=4;		
    FORTRAN(writevector,(btil, &neq[1],&number)); 
    */   
      	 /*i=381;
	for(k=0;k<3;k++){
	 if(nactdof[mt*i-3+k]>0){printf("node %d dof %d f %e fext %e  b %e \n",i,k+1,ftil[nactdof[mt*i-3+k]-1],fexttil[nactdof[mt*i-3+k]-1],btil[nactdof[mt*i-3+k]-1]);}
	} 
	 i=1553;
	for(k=0;k<3;k++){
	 if(nactdof[mt*i-3+k]>0){printf("node %d dof %d f %e fext %e  b %e \n",i,k+1,ftil[nactdof[mt*i-3+k]-1],fexttil[nactdof[mt*i-3+k]-1],btil[nactdof[mt*i-3+k]-1]);}
	} */
	
  ///@TODO: for coupled thermo-mechanical calculations b!=btil!!! debug this for quadratic elements
  if(ithermal[0]==3){
    
  }else{  
    for(k=0;k<neq[1];k++){b[k]=btil[k];}
  }
  //}
  SFREE(btil);
  
  fin= clock();

  SFREE(voldtil);SFREE(vinitil);SFREE(accoldtil);SFREE(veoldtil); 
  printf(" mafillsm_dstil : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);
  /* calculation of M* for dynamic calculations */ 
  
  if(*iexpl<=1){
    if(*nmethod==4){		    
      /* mechanical part */		    
      if(*ithermal!=2){
	scal1=*bet**dtime**dtime*(1.+alpha);
	for(k=0;k<neq[0];++k){
	  ad[k]=adbtil[k]+scal1*ad[k];
	}
	for(k=0;k<nzs[0];++k){
	  au[k]=aubtil[k]+scal1*au[k];
	}
	/* upper triangle of asymmetric matrix */
	if(*nasym>0){
	  for(k=nzs[2];k<nzs[2]+nzs[0];++k){
	    au[k]=aubtil[k]+scal1*au[k];
	  }
	}
	/* damping */
	
	if((*alpham>0.)||(*betam>0.)){
	  scal1=*gam**dtime*(1.+alpha);
	  for(k=0;k<neq[0];++k){
	    ad[k]+=scal1*adctil[k];
	  }
	  for(k=0;k<nzs[0];++k){
	    au[k]+=scal1*auctil[k];
	  }
	  
	  /* upper triangle of asymmetric matrix */
	  
	  if(*nasym>0){
	    for(k=nzs[2];k<nzs[2]+nzs[0];++k){
	      au[k]+=scal1*auctil[k];
	    }
	  }
	}
	
      }
      /* thermal part */		    
      if(*ithermal>1){
	for(k=neq[0];k<neq[1];++k){
	  ad[k]=adbtil[k]/(*dtime)+ad[k];
			}
	for(k=nzs[0];k<nzs[1];++k){
	  au[k]=aubtil[k]/(*dtime)+au[k];
			}
	/* upper triangle of asymmetric matrix */
	if(*nasym>0){
	  for(k=nzs[2]+nzs[0];k<nzs[2]+nzs[1];++k){
	    au[k]=aubtil[k]/(*dtime)+au[k];
	  }
	}
      }
    }	
  }
  /*	  int number=3;		
    FORTRAN(writematrix,(au,ad,irow,jq,&neq[1],&number));*/
    /*
        number=4;		
    FORTRAN(writematrix,(aub,adb,irow,jq,&neq[1],&number));
    
            number=6;		
    FORTRAN(writematrix,(aubtil,adbtil,irow,jq,&neq[1],&number));
     number=5;		
    FORTRAN(writevector,(b,&neq[1],&number));   */

 
       /*if(*nmethod==4 && *nk>668){
	i=668;
	for(k=0;k<3;k++){
	 if(nactdof[mt*i-3+k]>0){printf("node %d dof %d ad %e \n",i,k+1,ad[nactdof[mt*i-3+k]-1]);}
	}
      }*/ 
  //if(*nmethod==4){SFREE(aubtil);SFREE(adbtil);mass[0] = 0;}
  SFREE(adtil);SFREE(autil);SFREE(irowtil);SFREE(jqtil);SFREE(icoltil);
  
  /* update vold due to spcs to get gap right for rigid body movements */         
  if(*iinc==1 && *iit==1 &&*nmethod!=4){       
    NNEW(v,double,mt**nk);	       
    NNEW(volddummy,double,mt**nk);
    NNEW(veolddummy,double,mt**nk);
    NNEW(accolddummy,double,mt**nk);	       
    for(k=0;k<mt**nk;k++){volddummy[k]=0.0;veolddummy[k]=0.0;accolddummy[k]=0.0;}           
    memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);	       
    NNEW(vectornull,double,neq[1]);	       
    *iout=-1;         
    FORTRAN(resultsini_mortar,(nk,v,ithermal,iperturb,
			     nactdof,iout,volddummy,vectornull,nodeboun,ndirboun,
			     xbounact,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,
			     veold,accold,bet,gam,dtime,mi,vini));
   /* @TODO: Check if this is working right with resultsini instead of resultsini_mortar !!! */
   /* FORTRAN(resultsini,(nk,v,ithermal,filab,iperturb,f,fn,
       nactdof,iout,qa,volddummy,vectornull,nodeboun2,ndirboun2,
       xboun2,nboun2,ipompc2,nodempc2,coefmpc2,labmpc2,nmpc2,nmethod,cam,neq,
       veolddummy,accolddummy,bet,gam,dtime,mi,vini,nprint,prlab,
       &intpointvarm,&calcul_fn,&calcul_f,&calcul_qa,&calcul_cauchy,
       &ikin,&intpointvart,typeboun));
     */  
    memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);	     	
    SFREE(v);SFREE(vectornull);SFREE(volddummy);SFREE(veolddummy);SFREE(accolddummy);	    
  }
  if(*nmethod==4){
    	  /*for(i=0;i<neq[0];i++){
	   printf("i %d finini %e fextini %e  \n",i,finitil[i],fextinitil[i]); 
	  }
  	  for(i=0;i<neq[0];i++){
	   printf("i %d fint %e fext %e b %e \n",i,ftil[i],fexttil[i],b[i]); 
	  }*/
  }
  
  *ielas=0;
  *iout=0;
  iperturb[0]=iperturb_sav[0];
  iperturb[1]=iperturb_sav[1];
  
  *ndirboun2p=ndirboun2;*nodeboun2p=nodeboun2;*xboun2p=xboun2;
  *ipompc2p=ipompc2;*nodempc2p=nodempc2;*coefmpc2p=coefmpc2;
  *labmpc2p=labmpc2;*ikboun2p=ikboun2;*ilboun2p=ilboun2;*ikmpc2p=ikmpc2;*ilmpc2p=ilmpc2;
  *nslavspcp=nslavspc;*islavspcp=islavspc;*nslavmpcp=nslavmpc;*islavmpcp=islavmpc;
  *nslavspc2p=nslavspc2;*islavspc2p=islavspc2;*nslavmpc2p=nslavmpc2;*islavmpc2p=islavmpc2;
  *nmastspcp=nmastspc;*imastspcp=imastspc;*nmastmpcp=nmastmpc;*nmastmpc2p=nmastmpc2;*imastmpcp=imastmpc;*imastmpc2p=imastmpc2;   
  *auc2p=auc2;*adc2p=adc2;*irowc2p=irowc2;*icolc2p=icolc2;*jqc2p=jqc2;
  *aubdp=aubd;*irowbdp=irowbd;*jqbdp=jqbd;
  *aubdtilp=aubdtil;*irowbdtilp=irowbdtil;*jqbdtilp=jqbdtil;
  *aubdtil2p=aubdtil2;*irowbdtil2p=irowbdtil2;*jqbdtil2p=jqbdtil2;
  *auddp=audd;*irowddp=irowdd;*jqddp=jqdd;
  *auddtilp=auddtil;*irowddtilp=irowddtil;*jqddtilp=jqddtil;
  *auddtil2p=auddtil2;*irowddtil2p=irowddtil2;*jqddtil2p=jqddtil2;
  *auddinvp=auddinv;*irowddinvp=irowddinv;*jqddinvp=jqddinv;
  *irowtempp=irowtemp;*icoltempp=icoltemp;*jqtempp=jqtemp;
  *f_csp=f_cs;*f_cmp=f_cm;
  *nodeforc2p=nodeforc2;*ndirforc2p=ndirforc2;
  *xforc2p=xforc2;
  
  return;
}
