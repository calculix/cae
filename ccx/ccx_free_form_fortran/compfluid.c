/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2018 Guido Dhondt                     */

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

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"
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

static ITG num_cpus;

void compfluid(double **cop, ITG *nk, ITG **ipkonfp, ITG *konf, char **lakonfp,
    char **sidefacep, ITG *ifreestream, 
    ITG *nfreestream, ITG *isolidsurf, ITG *neighsolidsurf,
    ITG *nsolidsurf, ITG *nshcon, double *shcon,
    ITG *nrhcon, double *rhcon, double **voldp, ITG *ntmat_,ITG *nodeboun, 
    ITG *ndirboun, ITG *nboun, ITG *ipompc,ITG *nodempc, ITG *nmpc,
    ITG *ikmpc, ITG *ilmpc, ITG *ithermal, ITG *ikboun, ITG *ilboun,
    ITG *iturbulent, ITG *isolver, ITG *iexpl, double *ttime,
    double *time, double *dtime, ITG *nodeforc,ITG *ndirforc,double *xforc,
    ITG *nforc, ITG *nelemload, char *sideload, double *xload,ITG *nload,
    double *xbody,ITG *ipobody,ITG *nbody, ITG *ielmatf, char *matname,
    ITG *mi, ITG *ncmat_, double *physcon, ITG *istep, ITG *iinc,
    ITG *ibody, double *xloadold, double *xboun,
    double *coefmpc, ITG *nmethod, double *xforcold, double *xforcact,
    ITG *iamforc,ITG *iamload, double *xbodyold, double *xbodyact,
    double *t1old, double *t1, double *t1act, ITG *iamt1, double *amta,
    ITG *namta, ITG *nam, double *ampli, double *xbounold, double *xbounact,
    ITG *iamboun, ITG *itg, ITG *ntg, char *amname, double *t0, 
    ITG **nelemfacep,
    ITG *nface, double *cocon, ITG *ncocon, double *xloadact, double *tper,
    ITG *jmax, ITG *jout, char *set, ITG *nset, ITG *istartset,
    ITG *iendset, ITG *ialset, char *prset, char *prlab, ITG *nprint,
    double *trab, ITG *inotr, ITG *ntrans, char *filab, char *labmpc, 
    double *sti, ITG *norien, double *orab, char *jobnamef,char *tieset,
    ITG *ntie, ITG *mcs, ITG *ics, double *cs, ITG *nkon, ITG *mpcfree,
    ITG *memmpc_,double *fmpc,ITG *nef,ITG **inomatp,double *qfx,
    ITG *neifa,ITG *neiel,ITG *ielfa,ITG *ifaext,double *vfa,double *vel,
    ITG *ipnei,ITG *nflnei,ITG *nfaext,char *typeboun,ITG *neij,
    double *tincf,ITG *nactdoh,ITG *nactdohinv,ITG *ielorienf,char*jobnamec,
    ITG *ifatie,ITG *nstate_,double *xstate,char *orname,ITG *kon,
    double *ctrl,ITG *kode,double *velo,double *veloo,ITG *initial){

    /* main computational fluid dynamics routine */
  
    char cflag[1],*lakonf=NULL,*sideface=NULL,fncvg[132]="",*lakon=NULL;

  ITG *ipointer=NULL,*mast1=NULL,*irow=NULL,*icol=NULL,*jq=NULL,
      nzs=20000000,compressible,*ifabou=NULL,*ja=NULL,
      nfabou,im,iflag,*ipkon=NULL,*ielprop=NULL,*ielmat=NULL,
      *ipkonf=NULL,*nelemface=NULL,last=0,icyclic,*iau6=NULL,
      *inomat=NULL,ithermalref,*integerglob=NULL,iincf,
      iconvergence=0,i,*inum=NULL,iitf,ifreefa,*neielcp=NULL,
      *iponofa=NULL,*inofa=NULL,*ia=NULL,*ielpropf=NULL,
      icent=0,isti=0,iqfx=0,nfield,ndim,iorienglob,force=0,icfd=1,
      imach=0,ikappa=0,iit,iatleastonepressurebc,iturb=0,
      *inoel=NULL,*iponoel=NULL,icounter,ischeme=1,isimplec=0,
      iitt,iitg,iitp,iexceed,*iam=NULL,*jam=NULL,*iamorig=NULL,
      nz_num,*nestart=NULL,*ineighblock=NULL,*neighblock=NULL,
      nneighblock,iittf,iittfo,ip,ierrmax;

  ITG nelt,isym,itol,itmax,iunit,lrgw,*igwk=NULL,ligw,ierr,*iwork=NULL,iter,
      nsave,lenw,leniw;

  double *umfa=NULL,reltime,*doubleglob=NULL,relnormmin,
      *co=NULL,*vold=NULL,*coel=NULL,*cosa=NULL,*gradvel=NULL,*gradvfa=NULL,
      *xxn=NULL,*xxi=NULL,*xle=NULL,*xlen=NULL,*xlet=NULL,timef,dtimef,
      *cofa=NULL,*area=NULL,*xrlfa=NULL,reltimef,ttimef,*hcfa=NULL,*cvel=NULL,
      *au=NULL,*ad=NULL,*b=NULL,*volume=NULL,*body=NULL,*dy=NULL,
      *advfa=NULL,*ap=NULL,*bp=NULL,*xxj=NULL,*gradkel=NULL,*gradoel=NULL,
      *v=NULL,*cosb=NULL,dmin,tincfguess,
      *hel=NULL,*hfa=NULL,*auv=NULL,*adv=NULL,*bv=NULL,*sel=NULL,*gamma=NULL,
      *gradtfa=NULL,*gradtel=NULL,*umel=NULL,*cvfa=NULL,*gradpel=NULL,
      *eei=NULL,*ener=NULL,*thicke=NULL,*eme=NULL,c[9],*gradkfa=NULL,
      ptimef,*stn=NULL,*qfn=NULL,*hcel=NULL,*aua=NULL,a1,a2,a3,beta,
      *prop=NULL,*dp=NULL,*xxni=NULL,*xxnj=NULL,*xxicn=NULL,*xturb=NULL,
      *xmach=NULL,*xkappa=NULL,urelax,*flux=NULL,velnormo[5],velnorm[5],
      relnormt,relnormv,relnormp=0,relnormmax=1.e30,*temp=NULL,*auv6=NULL,
      *adv6=NULL,*auv3=NULL,*bv3=NULL,*vela=NULL,*velaa=NULL,*yy=NULL,
      *gradofa=NULL,betam=0.1,*gradpfa=NULL,*gradpcel=NULL,*gradpcfa=NULL,
      *am=NULL,*f1=NULL,*of2=NULL;

  double tol,*rgwk=NULL,err,*sb=NULL,*sx=NULL,*rwork=NULL,*rf=NULL;

  FILE *f3;

  co=*cop;
  ipkonf=*ipkonfp;lakonf=*lakonfp;
  nelemface=*nelemfacep;sideface=*sidefacep;
  vold=*voldp;inomat=*inomatp;

#ifdef SGI
  ITG token;
#endif
	  
  strcpy(fncvg,jobnamec);
  strcat(fncvg,".fcv");

  if((f3=fopen(fncvg,"w"))==NULL){
//  if((f3=fopen("fluidconvergence","w"))==NULL){
      printf("*ERROR in compfluid: cannot open cvg file for writing...");
      exit(0);
  }
  fprintf(f3,"temperature    velocity    pressure\n\n");

  urelax=.2;

  /* relative time at the end of the mechanical increment */

  reltime=(*time)/(*tper);

  /* open frd-file for fluids */

  FORTRAN(openfilefluid,(jobnamef));

  /* variables for multithreading procedure */

  ITG sys_cpus;
  char *env,*envloc,*envsys;
      
  num_cpus=0;
  sys_cpus=0;
  
  /* explicit user declaration prevails */
  
  envsys=getenv("NUMBER_OF_CPUS");
  if(envsys){
      sys_cpus=atoi(envsys);
      if(sys_cpus<0) sys_cpus=0;
  }
  
  /* automatic detection of available number of processors */
  
  if(sys_cpus==0){
      sys_cpus = getSystemCPUs();
      if(sys_cpus<1) sys_cpus=1;
  }
  
  /* local declaration prevails, if strictly positive */
  
  envloc = getenv("CCX_NPROC_CFD");
  if(envloc){
      num_cpus=atoi(envloc);
      if(num_cpus<0){
	  num_cpus=0;
      }else if(num_cpus>sys_cpus){
	  num_cpus=sys_cpus;
      }
  }
  
  /* else global declaration, if any, applies */
  
  env = getenv("OMP_NUM_THREADS");
  if(num_cpus==0){
      if (env)
	  num_cpus = atoi(env);
      if (num_cpus < 1) {
	  num_cpus=1;
      }else if(num_cpus>sys_cpus){
	  num_cpus=sys_cpus;
      }
  }
  
// next line is to be inserted in a similar way for all other paralell parts
  
  if(*nef<num_cpus) num_cpus=*nef;
  
  printf(" Using up to %" ITGFORMAT " cpu(s) for CFD.\n", num_cpus);
  
  pthread_t tid[num_cpus];
  
  /*  *iexpl==0:  structure:implicit, fluid:incompressible
      *iexpl==1:  structure:implicit, fluid:compressible
      *iexpl==2:  structure:explicit, fluid:incompressible
      *iexpl==3:  structure:explicit, fluid:compressible */

  if((*iexpl==1)||(*iexpl==3)){
      compressible=1;
  }else{
      compressible=0;
  }

  /* ischeme=1: ud
     ischeme=2: modified smart */

  ischeme=(ITG)floor(ctrl[47]);
  if(compressible==1){
      if(ischeme==1){
	  printf(" CFD scheme: upwind difference\n");
      }else{
	  printf(" CFD scheme: modified smart\n");
      }
  }else{
      printf(" CFD scheme: Gamma\n");
  }

  /* isimplec=0: simple
     isimplec=1: simplec */

  isimplec=(ITG)floor(ctrl[48]);
  if(compressible==1){
      if(isimplec==0){
	  printf(" algorithm: simple\n\n");
      }else{
	  printf(" algorithm: simplec\n\n");
      }
  }

  /* iitt: max number of transient iterations
     iitg: max number of geometric iterations (extrapol_*.f)
     iitp: max number of pressure iterations */

  iitt=(ITG)floor(ctrl[49]);
  iitg=(ITG)floor(ctrl[50]);
  iitp=(ITG)floor(ctrl[51]);

  /* if initial conditions are specified for the temperature, 
     it is assumed that the temperature is an unknown */

  ithermalref=*ithermal;
  if(*ithermal==1){
    *ithermal=2;
  }

  /* determining the matrix structure */
  
  NNEW(ipointer,ITG,*nef);
  NNEW(mast1,ITG,nzs);
  NNEW(irow,ITG,nzs);
  NNEW(icol,ITG,*nef);
  NNEW(jq,ITG,*nef+1);

  mastructf(nk,konf,ipkonf,lakonf,nef,icol,jq,&mast1,&irow,
	    isolver,ipointer,&nzs,ipnei,neiel,mi);

  SFREE(ipointer);SFREE(mast1);
 
  NNEW(iau6,ITG,6**nef);
  FORTRAN(create_iau6,(nef,ipnei,neiel,jq,irow,&nzs,iau6,lakonf));
		  
  if(compressible==0){
      NNEW(ia,ITG,nzs+*nef);
      NNEW(ja,ITG,*nef+1);
      NNEW(aua,double,nzs+*nef);
      FORTRAN(preconvert2slapcol,(irow,ia,jq,ja,&nzs,nef));
  }

  /* calculation geometric data */

  NNEW(coel,double,3**nef);
  NNEW(volume,double,*nef);
  NNEW(cosa,double,*nflnei);
  NNEW(cosb,double,*nflnei);
  NNEW(xxn,double,3**nflnei);
  NNEW(xxi,double,3**nflnei);
  NNEW(xxj,double,3**nflnei);
  NNEW(xxni,double,3**nflnei);
  NNEW(xxicn,double,3**nflnei);
  NNEW(xxnj,double,3**nflnei);
  NNEW(xle,double,*nflnei);
  NNEW(xlen,double,*nflnei);
  NNEW(xlet,double,*nflnei);
  NNEW(cofa,double,3**nface);
  NNEW(area,double,*nface);
  NNEW(xrlfa,double,3**nface);
  NNEW(rf,double,3**nface);

  /* closest distance to a node from a solid surface (dy) 
     distance from any element center to a solid surface (yy) */

  if(*iturbulent>0){
      NNEW(dy,double,*nsolidsurf);
      if(*iturbulent>2){
	  NNEW(yy,double,*nef);
      }
  }

  FORTRAN(initialcfd,(nef,ipkonf,konf,lakonf,co,coel,cofa,nface,
	  ielfa,area,ipnei,neiel,xxn,xxi,xle,xlen,xlet,xrlfa,cosa,
	  volume,neifa,xxj,cosb,&dmin,ifatie,cs,tieset,&icyclic,c,
	  neij,physcon,isolidsurf,nsolidsurf,dy,xxni,xxnj,xxicn,
	  nflnei,iturbulent,rf,yy,vel,velo,veloo));

  /* creating iam from neiel:
     1) get rid of zero neighbors
     2) get rid of cyclic symmetry neighbors
     3) for parallel execution: get rid of neighbors belonging
     to a different block */

  NNEW(jam,ITG,*nef+1);
  NNEW(iam,ITG,*nflnei+*nef);
  NNEW(iamorig,ITG,*nflnei+*nef);
  NNEW(nestart,ITG,num_cpus+1);
  NNEW(ineighblock,ITG,num_cpus+1);
  NNEW(neighblock,ITG,3**nflnei);
  
  FORTRAN(createblock,(nef,ipnei,neiel,iam,jam,iamorig,nflnei,
		     &nz_num,&num_cpus,nestart,ineighblock,
                     neighblock,&icyclic,neifa,ifatie,&nneighblock));
  RENEW(iam,ITG,nz_num);
  RENEW(iamorig,ITG,nz_num);
  RENEW(neighblock,ITG,3*nneighblock);
  NNEW(am,double,nz_num);

  /* storing pointers to the boundary conditions in ielfa */

  NNEW(ifabou,ITG,7**nfaext);
  FORTRAN(applyboun,(ifaext,nfaext,ielfa,ikboun,ilboun,
       nboun,typeboun,nelemload,nload,sideload,isolidsurf,nsolidsurf,
       ifabou,&nfabou,nface,nodeboun,ndirboun,ikmpc,ilmpc,labmpc,nmpc,
       nactdohinv,&compressible,&iatleastonepressurebc,ipkonf,kon,konf));
  RENEW(ifabou,ITG,nfabou);

  /* catalogueing the nodes for output purposes (interpolation at
     the nodes */
  
  NNEW(iponofa,ITG,*nk);
  NNEW(inofa,ITG,2**nface*4);

  FORTRAN(cataloguenodes,(iponofa,inofa,&ifreefa,ielfa,ifabou,ipkonf,
			  konf,lakonf,nface,nk));

  RENEW(inofa,ITG,2*ifreefa);

  /* material properties for athermal calculations 
     = calculation for which no initial thermal conditions
     were defined */

  NNEW(umfa,double,*nface);
  NNEW(umel,double,*nef);
  
  if(*ithermal==0){
      
      /* athermal incompressible calculations */

      /* calculating the dynamic viscosity at the element centers */
      
      FORTRAN(calcumel,(nef,vel,shcon,nshcon,ielmatf,ntmat_,
			    ithermal,mi,umel));

  }


  if(*ithermal!=0){
      NNEW(hcfa,double,*nface);
      NNEW(cvel,double,*nef);
      NNEW(cvfa,double,*nface);
  }

  if(*nbody>0) NNEW(body,double,4**nef);

  /* v is a auxiliary field: set to zero for the calls to
     tempload */

  NNEW(v,double,5**nk);

  /* next section is for stationary calculations */
  
  if(*nmethod==1){
      
      /* boundary conditions at the end of the mechanical
	 increment */
      
      FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
	     xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
             xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
             namta,nam,ampli,time,&reltime,ttime,dtime,ithermal,nmethod,
             xbounold,xboun,xbounact,iamboun,nboun,
             nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
	     co,v,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
	     ntrans,trab,inotr,vold,integerglob,doubleglob,tieset,istartset,
             iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,
             ipobody,iponoel,inoel,ipkon,kon,ielprop,prop,ielmat,
             shcon,nshcon,rhcon,nrhcon,cocon,ncocon,ntmat_,lakon));

      /* body forces (gravity, centrifugal and Coriolis forces */

      if(*nbody>0){
	  FORTRAN(inicalcbody,(nef,body,ipobody,ibody,xbody,coel,vel,lakonf,
                            nactdohinv,&icent));
      }
  }

  /* extrapolating the velocity from the elements centers to the face
     centers, thereby taking the boundary conditions into account */

  NNEW(gradvel,double,9**nef);
  NNEW(gradvfa,double,9**nface);

  extrapol_velmain(nface,ielfa,xrlfa,vel,vfa,
       ifabou,xbounact,ipnei,nef,&icyclic,c,ifatie,xxn,gradvel,
       gradvfa,neifa,rf,area,volume,xle,xxi,xxj,xlet,
       coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
       &iitg,&num_cpus);
/*      printf("after extrapol_velmain\n");
      for(i=0;i<*nface;i++){
	  printf("i=%d,cofa[0]=%f,vfa[1]=%f,%f\n",i,cofa[3*i+1],vfa[8*i+1],vfa[8*i+1]/cofa[3*i+1]);
	  }*/

  /* extrapolation of the pressure at the element centers
     to the face centers */

  NNEW(gradpel,double,3**nef);
  NNEW(gradpfa,double,3**nface);

  extrapol_pelmain(nface,ielfa,xrlfa,vel,vfa,
       ifabou,xbounact,nef,gradpel,gradpfa,neifa,rf,area,volume,
       xle,xxi,&icyclic,xxn,ipnei,ifatie,
       coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
       &iitg,c,&num_cpus);

  /* generate fields for the pressure correction gradients */

  NNEW(gradpcel,double,3**nef);
  NNEW(gradpcfa,double,3**nface);

  /* extrapolation of the temperature at the element centers
     to the face centers */

  if(*ithermal>0){

      NNEW(gradtel,double,3**nef);
      NNEW(gradtfa,double,3**nface);

/*	  printf("before extrapol_telmain: t\n");
	  for(i=0;i<*nef;i++){
	      printf("i=%d,co[0]=%f,vel[0]=%f %f\n",i,coel[3*i+1],vel[i],vel[i]/(1.+(1.-coel[3*i+1])*coel[3*i+1]/2.));
	      }*/

      extrapol_telmain(nface,ielfa,xrlfa,vel,vfa,
       ifabou,xbounact,nef,gradtel,gradtfa,neifa,rf,area,volume,
       xle,xxi,&icyclic,xxn,ipnei,ifatie,xload,xlet,xxj,
       coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
       &iitg,c,&num_cpus);

/*	  printf("after extrapol_telmain: t\n");
	      for(i=0;i<*nface;i++){
		  printf("i=%d,cofa[0]=%f,vfa[0]=%f %f\n",i,cofa[3*i+1],vfa[8*i],vfa[8*i]/(1.+(1.-cofa[3*i+1])*cofa[3*i+1]/2.));
		  }*/
	  
      /* calculating the heat conduction at the face centers */
	  
      FORTRAN(calchcfa,(nface,vfa,cocon,ncocon,ielmatf,ntmat_,
			mi,ielfa,hcfa));

      if(compressible==0){
	  
	  /* calculating the specific heat at constant volume at the 
	     face centers (secant value) */
	  
	  FORTRAN(calccvfa,(nface,vfa,shcon,nshcon,ielmatf,ntmat_,
			    mi,ielfa,cvfa,physcon));
      }else{

	  /* calculating the specific heat at constant volume at the 
	     face centers (secant value) */
	  
	  FORTRAN(calccvfacomp,(nface,vfa,shcon,nshcon,ielmatf,ntmat_,
				mi,ielfa,cvfa,physcon));
      }
  }

//  NNEW(flux,double,6**nef);
  NNEW(flux,double,*nflnei);
      
  if(compressible==0){

      /* calculating the density at the element centers */
      
      FORTRAN(calcrhoel,(nef,vel,rhcon,nrhcon,ielmatf,ntmat_,
			 ithermal,mi));
      
      /* calculating the density at the face centers */
      
      FORTRAN(calcrhofa,(nface,vfa,rhcon,nrhcon,ielmatf,ntmat_,
			 ithermal,mi,ielfa));

  }else{

      /* calculating the density at the element centers */
      
      calcrhoelcompmain(nef,vel,shcon,ielmatf,ntmat_,mi,
			&num_cpus);

/*	  printf("after calcrhoelcompmain: rho\n");
	  for(i=0;i<*nef;i++){
	      printf("i=%d,co[0]=%f,vel[5]=%f,%f\n",i,coel[3*i+1],vel[5**nef+i],vel[5**nef+i]*(1.+coel[3*i+1]*(1.-coel[3*i+1])/2.)*0.285714286);
	      }*/
      
      /* calculating the density at the face centers */
      
      if(ischeme==1){
	  calcrhofacomp_udmain(nface,vfa,shcon,ielmatf,ntmat_,
			       mi,ielfa,ipnei,vel,nef,flux,
			       &num_cpus,xxi,xle,gradpel,gradtel,
                               neij);
/*	      printf("after calcrhofacomp_udmain: rho\n");
	      for(i=0;i<*nface;i++){
		  printf("i=%d,cofa[0]=%f,vfa[5]=%f,%f\n",i,cofa[3*i+1],vfa[8*i+5],vfa[8*i+5]*(1.+cofa[3*i+1]*(1.-cofa[3*i+1])/2.)*0.285714286);
		  }*/
      }else{
	  calcrhofacomp_mod_smartmain(nface,vfa,shcon,ielmatf,ntmat_,
				      mi,ielfa,ipnei,vel,nef,flux,
				      gradpel,gradtel,xxj,xlet,
				      &num_cpus);
      }

  }
  
  /* calculating the initial mass flux */

  FORTRAN(calcinitialflux,(area,vfa,xxn,ipnei,nef,neifa,lakonf,flux));

  /* calculating the dynamic viscosity at the face centers */
  
  FORTRAN(calcumfa,(nface,vfa,shcon,nshcon,ielmatf,ntmat_,
		    ithermal,mi,ielfa,umfa));

  /* extrapolation of the turbulence variables at the element centers
     to the face centers */

  if(*iturbulent>0){

      NNEW(gradkel,double,3**nef);
      NNEW(gradkfa,double,3**nface);
      NNEW(gradoel,double,3**nef);
      NNEW(gradofa,double,3**nface);

      DMEMSET(vel,7**nef,8**nef,1.);

      extrapol_kelmain(nface,ielfa,xrlfa,vel,vfa,
	  ifabou,xbounact,nef,gradkel,gradkfa,neifa,rf,area,volume,
	  xle,xxi,&icyclic,xxn,ipnei,ifatie,xlet,xxj,
	  coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
	  umfa,physcon,&iitg,c,&num_cpus);

      extrapol_oelmain(nface,ielfa,xrlfa,vel,vfa,
	  ifabou,xbounact,nef,gradoel,gradofa,neifa,rf,area,volume,
	  xle,xxi,&icyclic,xxn,ipnei,ifatie,xlet,xxj,
	  coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
	  umfa,physcon,&iitg,c,dy,&num_cpus);

      if(*iturbulent>2){
	  NNEW(f1,double,*nef);
	  if(*iturbulent>3) NNEW(of2,double,*nef);
      }

  }

  /* calculating the time increment */

  FORTRAN(calcguesstincf,(nface,&dmin,vfa,umfa,cvfa,hcfa,ithermal,&tincfguess,
                          &compressible));

  /* start of the major loop */

  NNEW(advfa,double,*nface);
  NNEW(hfa,double,3**nface);

  NNEW(ap,double,*nface);
  NNEW(bp,double,*nface);

  NNEW(au,double,*nflnei+*nef);
  NNEW(ad,double,*nef);
  NNEW(b,double,*nef);

  NNEW(auv,double,*nflnei+*nef);

  NNEW(bv,double,3**nef);
  NNEW(hel,double,3**nef);
  NNEW(sel,double,3**nef);

  NNEW(rwork,double,*nef);

  NNEW(inum,ITG,*nk);

//  NNEW(velo,double,8**nef);
//  NNEW(veloo,double,8**nef);

  /* initializing velo and veloo */

  if(*initial==1){
      memcpy(&veloo[0],&vel[0],sizeof(double)*8**nef);
      memcpy(&velo[0],&vel[0],sizeof(double)*8**nef);
  }

  /* check output requests */

  if((strcmp1(&filab[1914],"MACH")==0)||
     (strcmp1(&filab[3132],"PTF")==0)||
     (strcmp1(&filab[3219],"TTF")==0)){
      imach=1;
  }
  
  if((strcmp1(&filab[3132],"PTF")==0)||
     (strcmp1(&filab[3219],"TTF")==0)){
      ikappa=1;
  }
  
  if(strcmp1(&filab[2088],"TURB")==0){
      iturb=1;
  }
  
  for(i=0;i<*nprint;i++){
      if(imach==0){
	  if((strcmp1(&prlab[6*i],"MACH")==0)||
	     (strcmp1(&prlab[6*i],"PTF")==0)||
	     (strcmp1(&prlab[6*i],"TTF")==0)){
	      imach=1;
	  }
      }
      if(ikappa==0){
	  if((strcmp1(&prlab[6*i],"PTF")==0)||
	     (strcmp1(&prlab[6*i],"TTF")==0)){
	      ikappa=1;
	  }
      }
      if(iturb==0){
	  if(strcmp1(&prlab[6*i],"TURB")==0){
	      iturb=1;
	  }
      }
  }

  iincf=0;

  if(*tincf<=0.) *tincf=tincfguess;
  printf("time increment for the CFD-calculations = %e\n\n",*tincf);

  ttimef=*ttime;
  timef=*time-*dtime;
  dtimef=*tincf;

  if(compressible==0){
      a1=1.5/dtimef;
      a2=-2./dtimef;
      a3=0.5/dtimef;
  }else{
      a1=1./dtimef;
      a2=-a1;
      a3=0.;
  }

  iittfo=iitt;

  NNEW(temp,double,8**nef);
  NNEW(gamma,double,*nface);

  icounter=0;

  do{

      iincf++;

      printf("fluid increment = %d\n",iincf);

//      if(iincf>5000) ischeme=2;

//      relnormmin=1.e30;
      ierrmax=0;

      timef+=dtimef;
      if((*time<timef)&&(*nmethod==4)){
	  dtimef-=(timef-*time);
	  timef=*time;
	  last=1;
	  beta=dtimef/(*tincf);
	  a1=(2.+beta)/(1.+beta)/dtimef;
	  a2=-(1.+beta)/beta/dtimef;
	  a3=1./(beta*(1.+beta))/dtimef;
      }

      /* starting iterations till convergence of the fluid increment */

      iittf=0;
      
      do{
	  iittf++;

	  printf("      iteration = %d\n",iittf);

      /* conditions for transient calculations */
     
      if(*nmethod==4){

          /* boundary conditions at end of fluid increment */

	  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
             xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
             xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
             namta,nam,ampli,&timef,&reltimef,&ttimef,&dtimef,ithermal,nmethod,
             xbounold,xboun,xbounact,iamboun,nboun,
             nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
	     co,v,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
             ntrans,trab,inotr,vold,integerglob,doubleglob,tieset,istartset,
             iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,
             ipobody,iponoel,inoel,ipkon,kon,ielprop,prop,ielmat,
             shcon,nshcon,rhcon,nrhcon,cocon,ncocon,ntmat_,lakon));

	  /* body forces (gravity, centrifugal and Coriolis forces) */
      
	  if(*nbody>0){
	      FORTRAN(calcbody,(nef,body,ipobody,ibody,xbody,coel,vel,lakonf,
				nactdohinv));
	  }

      }else if(icent==1){
	  
	  /* body forces (gravity, centrifugal and Coriolis forces;
             only if centrifugal forces are active => the ensuing
             Coriolis forces depend on the actual velocity) */
	  
	  FORTRAN(calcbody,(nef,body,ipobody,ibody,xbody,coel,vel,lakonf,
				nactdohinv));
      }

      /* updating of the material properties */

      if(*ithermal>0){
	  

	  if(compressible==0){

	      /* calculating material data 
		 density (elements+faces)
		 heat capacity at constant volume (elements+faces)
		 dynamic viscosity (elements+faces)
		 heat conduction (faces) */

	      FORTRAN(materialdata_cfd,(nef,vel,shcon,nshcon,ielmatf,
		      ntmat_,mi,cvel,vfa,cocon,ncocon,physcon,cvfa,
		      ithermal,nface,umel,umfa,ielfa,hcfa,rhcon,nrhcon));

	  }else{

	      /* calculating material data 
		 heat capacity at constant volume (elements+faces)
		 dynamic viscosity (elements+faces)
		 heat conduction (faces) */

	      materialdata_cfd_compmain(nef,vel,shcon,nshcon,ielmatf,
					ntmat_,mi,cvel,vfa,cocon,ncocon,
					physcon,cvfa,ithermal,nface,umel,
					umfa,ielfa,hcfa,&num_cpus);
	  }

      }

      /* filling the lhs and rhs's for the balance of momentum
	 equations */
      
      if(icounter==0){
	  DMEMSET(auv,0,*nflnei+*nef,0.);
	  DMEMSET(bv,0,3**nef,0.);
      }
	  
      if(compressible==0){
      
              /* calculate gamma (Ph.D. Thesis Jasak) */

	  calcgammavmain(nface,ielfa,vel,gradvel,gamma,xlet,xxj,ipnei,
			 &betam,nef,flux,&num_cpus);
	      
	  mafillvmain(nef,ipnei,neifa,neiel,vfa,xxn,area,
		       auv,&auv[*nflnei],jq,irow,&nzs,bv,vel,cosa,umfa,
		       xlet,xle,gradvfa,xxi,
		       body,volume,ielfa,lakonf,ifabou,nbody,
		       &dtimef,velo,veloo,sel,xrlfa,gamma,xxj,nactdohinv,&a1,
		       &a2,&a3,flux,&icyclic,c,ifatie,iau6,xxni,xxnj,
		       iturbulent,gradvel,of2,yy,umel);
	  
      }else{
		  
	  /* calculating the density at the element centers */
	  
	  calcrhoelcompmain(nef,vel,shcon,ielmatf,ntmat_,mi,
			    &num_cpus);

/*	  printf("after calcrhoelcompmain: rho\n");
	  for(i=0;i<*nef;i++){
	      printf("i=%d,co[0]=%f,vel[5]=%f,%f\n",i,coel[3*i+1],vel[5**nef+i],vel[5**nef+i]*(1.+coel[3*i+1]*(1.-coel[3*i+1])/2.)*0.285714286);
	      }*/

	  /* calculating the density at the face centers
		 (gamma method) */
	
	  if(ischeme==1){
	      calcrhofacomp_udmain(nface,vfa,shcon,ielmatf,ntmat_,
				   mi,ielfa,ipnei,vel,nef,flux,
				   &num_cpus,xxi,xle,gradpel,
                                   gradtel,neij);
/*	      printf("after calcrhofacomp_udmain: rho\n");
	      for(i=0;i<*nface;i++){
		  printf("i=%d,cofa[0]=%f,vfa[5]=%f,%f\n",i,cofa[3*i+1],vfa[8*i+5],vfa[8*i+5]*(1.+cofa[3*i+1]*(1.-cofa[3*i+1])/2.)*0.285714286);
		  }*/
	  }else{
	      calcrhofacomp_mod_smartmain(nface,vfa,shcon,ielmatf,ntmat_,
					  mi,ielfa,ipnei,vel,nef,flux,
					  gradpel,gradtel,xxj,xlet,
					  &num_cpus);
	  }
/*	      printf("before hrv_udmain: flux\n");
	      for(i=0;i<*nflnei;i++){
		  printf("i=%d,flux[i]=%f\n",i,flux[i]);
		  }*/
      
	  /* convection scheme */

	  if(ischeme==1){
	      hrv_udmain(nface,ielfa,vel,ipnei,nef,flux,vfa,&num_cpus,
			 xxi,xle,gradvel,neij);
	  }else{
	      hrv_mod_smartmain(nface,ielfa,vel,gradvel,xlet,xxj,ipnei,
				nef,flux,vfa,&num_cpus,gamma);
	  }
/*	  printf("before mafillvcompmain: vx\n");
	  for(i=0;i<*nface;i++){
	      printf("i=%d,cofa[0]=%f,vfa[1]=%f,%f\n",i,cofa[3*i+1],vfa[8*i+1],vfa[8*i+1]/cofa[3*i+1]);
	      }*/
/*	  printf("before mafillvcompmain: gradvel\n");
	  for(i=0;i<1080;i++){
	      printf("i=%d,gradvel=%f\n",i,gradvel[i]);
	      }
	  printf("before mafillvcompmain: gradvfa\n");
	  for(i=0;i<4518;i++){
	      printf("i=%d,gradvfa=%f\n",i,gradvfa[i]);
	      }*/
	      
	  mafillvcompmain(nef,ipnei,neifa,neiel,vfa,xxn,area,
		       auv,&auv[*nflnei],jq,irow,&nzs,bv,vel,cosa,umfa,
                       xlet,xle,gradvfa,xxi,
		       body,volume,ielfa,lakonf,ifabou,nbody,
		       &dtimef,velo,veloo,sel,xrlfa,gamma,xxj,nactdohinv,&a1,
		       &a2,&a3,flux,&icyclic,c,ifatie,iau6,xxni,xxnj,
		       iturbulent,gradvel,of2,yy,umel);
      }

      for(i=0;i<*nef;i++){rwork[i]=1./auv[*nflnei+i];}

      /* underrelaxation of the velocity only for compressible
         simple scheme */

      if((compressible==1)&&(isimplec==0)){
	  memcpy(&temp[*nef],&vel[*nef],sizeof(double)*3**nef);
      }

      /* reordering the lhs (getting rid of zeros)  */

      reorderlhsmain(auv,am,iamorig,&nz_num,&num_cpus);

      /* modifying the rhs (taking common contributions
         between the cpu-blocks into account) */

      FORTRAN(reorderrhs,(auv,&bv[0],&vel[*nef],neighblock,&nneighblock));
      FORTRAN(reorderrhs,(auv,&bv[*nef],&vel[2**nef],neighblock,&nneighblock));
      FORTRAN(reorderrhs,(auv,&bv[2**nef],&vel[3**nef],neighblock,&nneighblock));

      /* calculation of the momentum residual */

      calcresvfluidmain(nef,am,&bv[0],&auv[*nflnei],iam,jam,
			&vel[*nef],&relnormv,nestart,&num_cpus);

      dgmresmain(nef,&bv[0],&vel[*nef],&nz_num,iam,jam,am,
		      &isym,&itol,&tol,&itmax,&iter,
                      &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
		      &ligw,rwork,iwork,nestart,&num_cpus);
      if(ierr>1){
	  printf("*WARNING in compfluid: error message from predgmres (v_x)=%d\n",ierr);
	  if(ierr>ierrmax)ierrmax=ierr;
      }

      dgmresmain(nef,&bv[*nef],&vel[2**nef],&nz_num,iam,jam,am,
		      &isym,&itol,&tol,&itmax,&iter,
                      &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
                      &ligw,rwork,iwork,nestart,&num_cpus);
      if(ierr>1){
	  printf("*WARNING in compfluid: error message from predgmres (v_y)=%d\n",ierr);
	  if(ierr>ierrmax)ierrmax=ierr;
      }

      dgmresmain(nef,&bv[2**nef],&vel[3**nef],&nz_num,iam,jam,am,
		      &isym,&itol,&tol,&itmax,&iter,
                      &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
                      &ligw,rwork,iwork,nestart,&num_cpus);
      if(ierr>1){
	  printf("*WARNING in compfluid: error message from predgmres (v_z)=%d\n",ierr);
	  if(ierr>ierrmax)ierrmax=ierr;
      }

      /* underrelaxation of the velocity only for compressible
         simple scheme */

/*      printf("mafillvcompmain: vx original\n");
      for(i=0;i<*nef;i++){
	  printf("i=%d,co[0]=%f,vel[1]=%f,%f\n",i,coel[3*i+1],temp[*nef+i],temp[*nef+i]/coel[3*i+1]);
	  }

      printf("mafillvcompmain: vx new\n");
      for(i=0;i<*nef;i++){
	  printf("i=%d,co[0]=%f,vel[1]=%f,%f\n",i,coel[3*i+1],vel[*nef+i],vel[*nef+i]/coel[3*i+1]);
	  }*/

      if((compressible==1)&&(isimplec==0)){
	  for(i=*nef;i<4**nef;i++){vel[i]=0.8*vel[i]+0.2*temp[i];}
      }

/*      printf("after mafillvcompmain: vx mix\n");
      for(i=0;i<*nef;i++){
	  printf("i=%d,co[0]=%f,vel[1]=%f,%f\n",i,coel[3*i+1],vel[*nef+i],vel[*nef+i]/coel[3*i+1]);
	  }*/

      for(iitf=0;iitf<1;iitf++){
	  
	  /* generating vol/ad and v* at the face centers (advfa and hfa) */

	  if((compressible==0)||(isimplec==0)){
	      extrapolate_d_v_simplemain(nface,ielfa,xrlfa,&auv[*nflnei],advfa,hfa,
					 &icyclic,c,ifatie,vel,nef,volume,&num_cpus);
	  }else{
	      extrapolate_d_v_simplecmain(nface,ielfa,xrlfa,&auv[*nflnei],advfa,hfa,
					  &icyclic,c,ifatie,vel,nef,volume,auv,
					  ipnei,&num_cpus);
	  }
	  
          /* calculate the mass flow based on the newly calculated
             velocity */
/*	      printf("before calcfluxmain: flux\n");
	      for(i=0;i<*nflnei;i++){
		  printf("i=%d,flux[i]=%f\n",i,flux[i]);
		  }*/

	  calcfluxmain(area,vfa,xxn,ipnei,nef,neifa,flux,xxj,gradpfa,xlet,xle,vel,
		       advfa,ielfa,neiel,ifabou,hfa,&num_cpus);
/*	      printf("after calcfluxmain: flux\n");
	      for(i=0;i<*nflnei;i++){
		  printf("i=%d,flux[i]=%f\n",i,flux[i]);
		  }*/

	  /* calculating the lhs and rhs of the equation system to determine
	     p' (balance of mass) */

	  if(compressible==0){

	      /* incompressible media */

              /* temporarily store the pressure in temp */

	      memcpy(&temp[4**nef],&vel[4**nef],sizeof(double)**nef);
		  
	      /* first iteration: calculating both lhs and rhs */
		  
	      DMEMSET(ad,0,*nef,0.);
	      DMEMSET(au,0,nzs,0.);
	      DMEMSET(b,0,*nef,0.);
	      
	      mafillpmain(nef,lakonf,ipnei,neifa,neiel,vfa,area,
				   advfa,xlet,cosa,volume,au,ad,jq,irow,ap,
				   ielfa,ifabou,xle,b,xxn,nef,
				   &nzs,hfa,gradpel,bp,xxi,neij,xlen,cosb,
			           &iatleastonepressurebc,iau6,xxicn,flux);

	      FORTRAN(convert2slapcol,(au,ad,jq,&nzs,nef,aua));
	      
	      nelt=nzs+*nef;
	      isym=1;

	      /* next line was changed from 10 to 3 on 22.12.2016 */

	      nsave=3;
	      itol=0;
	      tol=1.e-6;

	      /* next line was changed from 110 to 10 on 22.12.2016 */

	      itmax=10;
	      iunit=0;
	      lenw=131+17**nef+2*nelt;
	      NNEW(rgwk,double,lenw);
	      leniw=32+4**nef+2*nelt;
	      NNEW(igwk,ITG,leniw);

              /* initial guess: 0 */

	      DMEMSET(vel,4**nef,5**nef,0.);
      
	      FORTRAN(dslugm,(nef,&b[0],&vel[4**nef],&nelt,ia,ja,aua,
			      &isym,&nsave,&itol,&tol,&itmax,&iter,
			      &err,&ierr,&iunit,rgwk,&lenw,igwk,&leniw));
	      SFREE(rgwk);SFREE(igwk);

              /* non-orthogonal pressure correction */

	      /* calculate the p' gradient at the 
		 face centers */

	      for(ip=0;ip<iitp;ip++){

	      iflag=1;
	      extrapol_dpelmain(nface,ielfa,xrlfa,vel,vfa,
                     ifabou,xbounact,nef,gradpcel,gradpcfa,neifa,rf,area,volume,
                     xle,xxi,&icyclic,xxn,ipnei,ifatie,
                     coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,
		     nactdoh,&iflag,xxj,xlet,c,&num_cpus);
	      
	      /* update the right hand side (taking skewness of
		 elements into account) */
	      
	      DMEMSET(b,0,*nef,0.);
	      rhspmain(nef,lakonf,ipnei,neifa,neiel,vfa,area,
		       advfa,xlet,cosa,volume,au,ad,jq,irow,ap,ielfa,ifabou,xle,
		       b,xxn,nef,&nzs,hfa,gradpel,bp,xxi,neij,xlen,
		       &iatleastonepressurebc,xxicn,flux,xxnj,gradpcfa,cosb);
	      
	      /* calculate a better pressure correction p' */
	      
	      nelt=nzs+*nef;
	      isym=1;
	      nsave=3;
	      itol=0;
	      tol=1.e-6;
	      itmax=10;
	      iunit=0;
	      lenw=131+17**nef+2*nelt;
	      NNEW(rgwk,double,lenw);
	      leniw=32+4**nef+2*nelt;
	      NNEW(igwk,ITG,leniw);

              /* initial guess: 0 */
	      
	      DMEMSET(vel,4**nef,5**nef,0.);
	      
	      FORTRAN(dslugm,(nef,&b[0],&vel[4**nef],&nelt,ia,ja,aua,
			      &isym,&nsave,&itol,&tol,&itmax,&iter,
			      &err,&ierr,&iunit,rgwk,&lenw,igwk,&leniw));
			      SFREE(rgwk);SFREE(igwk);

	      }// end loop iitp
	      
	  }else{
	      
	      /* compressible media */

              /* temporarily store the pressure in temp */

	      memcpy(&temp[4**nef],&vel[4**nef],sizeof(double)**nef);
		  
	      DMEMSET(au,0,*nflnei+*nef,0.);
	      DMEMSET(b,0,*nef,0.);

/*	  printf("b before (0): p\n");
	  for(i=0;i<*nef;i++){
	      printf("i=%d,co[0]=%f,b=%f\n",i,coel[3*i+1],b[i]);
	      }*/
	      
	      mafillpcompmain(nef,lakonf,ipnei,neifa,neiel,vfa,area,
			  advfa,xlet,cosa,volume,au,&au[*nflnei],jq,irow,ap,
			  ielfa,ifabou,xle,b,xxn,nef,
			  &nzs,hfa,gradpel,bp,xxi,neij,xlen,cosb,
			  ielmatf,mi,&a1,&a2,&a3,velo,veloo,&dtimef,shcon,
			  ntmat_,vel,nactdohinv,xrlfa,flux,iau6,xxicn,
                          gamma);

/*	  printf("b after (0): p\n");
	  for(i=0;i<*nef;i++){
	      printf("i=%d,co[0]=%f,b=%f\n",i,coel[3*i+1],b[i]);
	      }*/

	      for(i=0;i<*nef;i++){rwork[i]=1./au[*nflnei+i];}

              /* initial guess: 0 */

	      DMEMSET(vel,4**nef,5**nef,0.);

	      /* reordering the lhs (getting rid of zeros) */
	      
	      reorderlhsmain(au,am,iamorig,&nz_num,&num_cpus);

              /* no change of rhs (reorderrhs) since initial guess is zero */

/*	  printf("b (0): p\n");
	  for(i=0;i<*nef;i++){
	      printf("i=%d,co[0]=%f,b=%f\n",i,coel[3*i+1],b[i]);
	      }*/

	      dgmresmain(nef,&b[0],&vel[4**nef],&nz_num,iam,jam,am,
				 &isym,&itol,&tol,&itmax,&iter,
				 &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
				 &ligw,rwork,iwork,nestart,&num_cpus);

	      if(ierr>1){
		  printf("*WARNING in compfluid: error message from predgmres (p)=%d\n",ierr);
		  if(ierr>ierrmax)ierrmax=ierr;
	      }

/*	  printf("pressure correction (0): p\n");
	  for(i=0;i<*nef;i++){
	      printf("i=%d,co[0]=%f,vel[5]=%f\n",i,coel[3*i+1],velo[4**nef+i]);
	      }

	  printf("pressure correction (1): p\n");
	  for(i=0;i<*nef;i++){
	      printf("i=%d,co[0]=%f,vel[5]=%f\n",i,coel[3*i+1],vel[4**nef+i]);
	      }*/

              /* non-orthogonal pressure correction */

	      for(ip=0;ip<iitp;ip++){

	      /* calculate the p' gradient at the 
		 face centers */

	      iflag=1;
	      extrapol_dpelmain(nface,ielfa,xrlfa,vel,vfa,
                     ifabou,xbounact,nef,gradpcel,gradpcfa,neifa,rf,area,volume,
                     xle,xxi,&icyclic,xxn,ipnei,ifatie,
                     coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,
		     nactdoh,&iflag,xxj,xlet,c,&num_cpus);
	      
	      /* update the right hand side (taking skewness of
		 elements into account) */

	      DMEMSET(b,0,*nef,0.);
	      
	      rhspcompmain(nef,lakonf,ipnei,neifa,neiel,vfa,area,
			  advfa,xlet,cosa,volume,au,&au[*nflnei],jq,irow,ap,
			  ielfa,ifabou,xle,b,xxn,nef,
			  &nzs,hfa,gradpel,bp,xxi,neij,xlen,cosb,
			  ielmatf,mi,&a1,&a2,&a3,velo,veloo,&dtimef,shcon,
			  ntmat_,vel,nactdohinv,xrlfa,flux,iau6,xxicn,
			   gamma,xxnj,gradpcfa);

	      for(i=0;i<*nef;i++){rwork[i]=1./au[*nflnei+i];}

	      /* reordering the lhs (getting rid of zeros) */
	      
	      reorderlhsmain(au,am,iamorig,&nz_num,&num_cpus);

	      /* modifying the rhs (taking common contributions
		 between the cpu-blocks into account) */
	      
	      FORTRAN(reorderrhs,(au,&b[0],&vel[4**nef],neighblock,&nneighblock));

              /* initial guess: 0 */

	      DMEMSET(vel,4**nef,5**nef,0.);

              /* calculation of the mass conservation residual */

	      FORTRAN(calcrespfluid,(nef,&b[0],&auv[*nflnei],
                            &temp[4**nef],&relnormp));

	      dgmresmain(nef,&b[0],&vel[4**nef],&nz_num,iam,jam,am,
				 &isym,&itol,&tol,&itmax,&iter,
				 &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
				 &ligw,rwork,iwork,nestart,&num_cpus);

	      if(ierr>1){
		  printf("*WARNING in compfluid: error message from predgmres (p)=%d\n",ierr);
		  if(ierr>ierrmax)ierrmax=ierr;
	      }

/*	  printf("pressure correction (2): p\n");
	  for(i=0;i<*nef;i++){
	      printf("i=%d,co[0]=%f,vel[5]=%f\n",i,coel[3*i+1],vel[4**nef+i]);
	      }*/

	      }   // end loop iitp
	      
	  }

          /* calculate the p' gradient at the 
             element centers */

	  iflag=0;
          extrapol_dpelmain(nface,ielfa,xrlfa,vel,vfa,
                  ifabou,xbounact,nef,gradpcel,gradpcfa,neifa,rf,area,volume,
                  xle,xxi,&icyclic,xxn,ipnei,ifatie,
                  coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,
		  nactdoh,&iflag,xxj,xlet,c,&num_cpus);
	  
	  /* correction of the velocity v* at the element centers due
             to the pressure change in order to get v** */

//	  if((compressible==0)||(isimplec==0)){
	      correctvelmain(&auv[*nflnei],nef,volume,gradpcel,vel,&num_cpus);
/*	  }else{
	      FORTRAN(correctvel_simplec,(&auv[*nflnei],nef,volume,gradpcel,vel,
                                      ipnei,auv));
				      }*/

/*      printf("after correctvelmain: vx \n");
      for(i=0;i<*nef;i++){
	  printf("i=%d,co[0]=%f,vel[1]=%f,%f\n",i,coel[3*i+1],temp[*nef+i],temp[*nef+i]/coel[3*i+1]);
	  }*/
	  
      }

      /* correcting the flux to get mf** */
      
      if(compressible==0){
	  correctfluxmain(nef,ipnei,neifa,neiel,flux,vfa,advfa,area,
			  vel,xlet,ielfa,xle,ifabou,&num_cpus);
      }else{

	  /* correcting the density at the element centers */
      
/*	  FORTRAN(correctrhoelcomp,(nef,vel,shcon,ielmatf,ntmat_,
	  mi));*/
      
	  /* correcting the density at the face centers */
      
/*	  FORTRAN(correctrhofacomp,(nface,vfa,shcon,ielmatf,ntmat_,
	      mi,ielfa,ipnei,vel,nef,flux,gradpel,gradtel,xxj,
	      &betam,xlet));*/
		  
/*	  FORTRAN(calcrhofacomp,(nface,vfa,shcon,ielmatf,ntmat_,
		      mi,ielfa,ipnei,vel,nef,flux,gradpel,gradtel,xxj,
		      &betam,xlet));*/

          /* correcting the flux at the faces */

	  correctfluxcompmain(nef,ipnei,neifa,neiel,flux,vfa,advfa,area,
			      vel,xlet,ielfa,xle,ifabou,ielmatf,mi,shcon,
			      ntmat_,&num_cpus);
/*	      printf("after correctfluxcompmain: flux\n");
	      for(i=0;i<*nflnei;i++){
		  printf("i=%d,flux[i]=%f\n",i,flux[i]);
		  }*/
      }

      /* correcting the pressure to get p* */

      /* underrelaxation always except for the compressible simplec
         scheme */

      if((compressible==1)&&(isimplec==1)){
	  for(i=4**nef;i<5**nef;i++){vel[i]=vel[i]+temp[i];}
      }else{
	  for(i=4**nef;i<5**nef;i++){vel[i]=0.2*vel[i]+temp[i];}
      }

      extrapol_pelmain(nface,ielfa,xrlfa,vel,vfa,
	      ifabou,xbounact,nef,gradpel,gradpfa,neifa,rf,area,volume,
	      xle,xxi,&icyclic,xxn,ipnei,ifatie,
	      coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
	      &iitg,c,&num_cpus);
/*	      printf("after extrapol_pelmain: p\n");
	      for(i=0;i<*nface;i++){
		  printf("i=%d,cofa[0]=%f,vfa[4]=%f\n",i,cofa[3*i+1],vfa[8*i+4]);
		  }*/

      if(*ithermal>0){

          /* calculating the lhs and rhs of the energy equation */

	  DMEMSET(ad,0,*nef,0.);
	  DMEMSET(au,0,*nflnei+*nef,0.);
	  DMEMSET(b,0,*nef,0.);

	  if(compressible==0){
      
             /* calculate gamma (Ph.D. Thesis Jasak) */

	      calcgammatmain(nface,ielfa,vel,gradtel,gamma,xlet,xxj,ipnei,
			     &betam,nef,flux,&num_cpus);

	      mafilltmain(nef,ipnei,neifa,neiel,vfa,xxn,area,
	       au,&au[*nflnei],jq,irow,&nzs,b,vel,umel,xlet,xle,gradtfa,xxi,
	       body,volume,ielfa,lakonf,ifabou,nbody,nef,
	       &dtimef,velo,veloo,cvfa,hcfa,cvel,gradvel,xload,gamma,
	       xrlfa,xxj,nactdohinv,&a1,&a2,&a3,flux,iau6,xxni,xxnj,
	       iturbulent,of2);

	  }else{

              /* underrelaxation of the temperature only for the
                 compressible simple scheme */

	      if(isimplec==0){
		  memcpy(&temp[0],&vel[0],sizeof(double)**nef);
	      }
 
              /* convective scheme */


/*	  printf("before hrt_udmain: t\n");
	      for(i=0;i<*nface;i++){
		  printf("i=%d,cofa[0]=%f,vfa[0]=%f %f\n",i,cofa[3*i+1],vfa[8*i],vfa[8*i]/(1.+(1.-cofa[3*i+1])*cofa[3*i+1]/2.));
		  }*/
//	      if(ischeme==1){
	      hrt_udmain(nface,ielfa,vel,ipnei,nef,flux,vfa,&num_cpus,
			 xxi,xle,gradtel,neij);
/*	      }else{
		  FORTRAN(hrt_mod_smart,(nface,ielfa,vel,gradtel,gamma,
                                    xlet,xxn,xxj,ipnei,&betam,nef,flux,vfa));
				    }*/

/*	  printf("after hrt_udmain: t\n");
	      for(i=0;i<*nface;i++){
		  printf("i=%d,cofa[0]=%f,vfa[0]=%f %f\n",i,cofa[3*i+1],vfa[8*i],vfa[8*i]/(1.+(1.-cofa[3*i+1])*cofa[3*i+1]/2.));
		  }*/

	      mafilltcompmain(nef,ipnei,neifa,neiel,vfa,xxn,area,
	       au,&au[*nflnei],jq,irow,&nzs,b,vel,umel,xlet,xle,gradtfa,xxi,
	       body,volume,ielfa,lakonf,ifabou,nbody,nef,
	       &dtimef,velo,veloo,cvfa,hcfa,cvel,gradvel,xload,gamma,
	       xrlfa,xxj,nactdohinv,&a1,&a2,&a3,flux,iau6,xxni,xxnj,
	       iturbulent,of2);
	  }

	  for(i=0;i<*nef;i++){rwork[i]=1./au[*nflnei+i];}

	  /* reordering the lhs (getting rid of zeros) */
	  
	  reorderlhsmain(au,am,iamorig,&nz_num,&num_cpus);

	  /* modifying the rhs (taking common contributions
	     between the cpu-blocks into account) */
	  
	  FORTRAN(reorderrhs,(au,&b[0],&vel[0],neighblock,&nneighblock));

          /* calculation of the energy residual */

	  calcrestfluidmain(nef,am,&b[0],&auv[*nflnei],iam,jam,
                            &vel[0],&relnormt,nestart,&num_cpus);

	  dgmresmain(nef,&b[0],&vel[0],&nz_num,iam,jam,am,
			     &isym,&itol,&tol,&itmax,&iter,
			     &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
			     &ligw,rwork,iwork,nestart,&num_cpus);


	  if(ierr>1){
	      printf("*WARNING in compfluid: error message from predgmres (T)=%d\n",ierr);
	      if(ierr>ierrmax)ierrmax=ierr;
	  }

	  /* underrelaxation of the temperature only for compressible
	     simple scheme */

/*	  printf("temperature new: t\n");
	  for(i=0;i<*nef;i++){
	      printf("i=%d,co[0]=%f,vel[0]=%f %f\n",i,coel[3*i+1],vel[i],vel[i]/(1.+(1.-coel[3*i+1])*coel[3*i+1]/2.));
	      }*/

	  if((compressible==1)&&(isimplec==0)){
	      for(i=0;i<*nef;i++){vel[i]=0.8*vel[i]+0.2*temp[i];}
	  }

/*	  printf("temperature mix: t\n");
	  for(i=0;i<*nef;i++){
	      printf("i=%d,co[0]=%f,vel[0]=%f %f\n",i,coel[3*i+1],vel[i],vel[i]/(1.+(1.-coel[3*i+1])*coel[3*i+1]/2.));
	      }*/

	  /* extrapolation of the temperature at the element centers
	     to the face centers */

	  extrapol_telmain(nface,ielfa,xrlfa,vel,vfa,
		  ifabou,xbounact,nef,gradtel,gradtfa,neifa,rf,area,volume,
		  xle,xxi,&icyclic,xxn,ipnei,ifatie,xload,xlet,xxj,
		  coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
		  &iitg,c,&num_cpus);

/*	  printf("after extrapol_telmain: t\n");
	      for(i=0;i<*nface;i++){
		  printf("i=%d,cofa[0]=%f,vfa[0]=%f %f\n",i,cofa[3*i+1],vfa[8*i],vfa[8*i]/(1.+(1.-cofa[3*i+1])*cofa[3*i+1]/2.));
		  }*/

          /* recalculating the density  */
      
	  if(compressible==0){
	      
	      /* calculating the density at the element centers */
	      
	      FORTRAN(calcrhoel,(nef,vel,rhcon,nrhcon,ielmatf,ntmat_,
				 ithermal,mi));
	      
	      /* calculating the density at the face centers */
	      
	      FORTRAN(calcrhofa,(nface,vfa,rhcon,nrhcon,ielmatf,ntmat_,
				 ithermal,mi,ielfa));

	  }

      }
      
      if(*iturbulent>0){

          /* calculating the lhs and rhs of the k-equation */

	  DMEMSET(au,0,*nflnei+*nef,0.);
	  DMEMSET(b,0,*nef,0.);

	  if(compressible==0){
      
              /* calculate gamma (Ph.D. Thesis Jasak) */

	      FORTRAN(calcgammak,(nface,ielfa,vel,gradkel,gamma,xlet,xxn,xxj,
				  ipnei,&betam,nef,flux));

	      mafillkmain(nef,ipnei,neifa,neiel,vfa,xxn,area,
	       au,&au[*nflnei],jq,irow,&nzs,b,vel,umfa,xlet,xle,gradkfa,xxi,
	       body,volume,ielfa,lakonf,ifabou,nbody,nef,
	       &dtimef,velo,veloo,cvfa,hcfa,cvel,gradvel,xload,gamma,
	       xrlfa,xxj,nactdohinv,&a1,&a2,&a3,flux,iau6,xxni,xxnj,
	       iturbulent,f1,of2,yy,umel,gradkel,gradoel);
	  }else{

              /* underrelaxation of the temperature only for the
                 compressible simple scheme */

	      if(isimplec==0){
		  memcpy(&temp[5**nef],&vel[5**nef],sizeof(double)**nef);
	      }
 
              /* convective scheme */

//	      if(ischeme==1){
	          hrk_udmain(nface,ielfa,vel,ipnei,nef,flux,vfa,&num_cpus);
/*	      }else{
		  FORTRAN(hrk_mod_smart,(nface,ielfa,vel,gradtel,gamma,
                                    xlet,xxn,xxj,ipnei,&betam,nef,flux,vfa));
				    }*/

	      mafillkcompmain(nef,ipnei,neifa,neiel,vfa,xxn,area,
	       au,&au[*nflnei],jq,irow,&nzs,b,vel,umfa,xlet,xle,gradkfa,xxi,
	       body,volume,ielfa,lakonf,ifabou,nbody,nef,
	       &dtimef,velo,veloo,cvfa,hcfa,cvel,gradvel,xload,
	       xrlfa,xxj,nactdohinv,&a1,&a2,&a3,flux,iau6,xxni,xxnj,
	       iturbulent,f1,of2,yy,umel,gradkel,gradoel);
	  }

	  for(i=0;i<*nef;i++){rwork[i]=1./au[*nflnei+i];}

	  /* reordering the lhs (getting rid of zeros) */
	  
	  reorderlhsmain(au,am,iamorig,&nz_num,&num_cpus);

	  /* modifying the rhs (taking common contributions
	     between the cpu-blocks into account) */
	  
	  FORTRAN(reorderrhs,(au,&b[0],&vel[6**nef],neighblock,&nneighblock));

	  dgmresmain(nef,&b[0],&vel[5**nef],&nz_num,iam,jam,am,
			     &isym,&itol,&tol,&itmax,&iter,
			     &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
			     &ligw,rwork,iwork,nestart,&num_cpus);
	  if(ierr>1){
	      printf("*WARNING in compfluid: error message from predgmres (k)=%d\n",ierr);
	      if(ierr>ierrmax)ierrmax=ierr;
	  }

	  /* underrelaxation of the temperature only for compressible
	     simple scheme */

	  if((compressible==1)&&(isimplec==0)){
	      for(i=5**nef;i<6**nef;i++){vel[i]=0.8*vel[i]+0.2*temp[i];}
	  }

          /* calculating the lhs and rhs of the omega-equation */

	  DMEMSET(au,0,*nflnei+*nef,0.);
	  DMEMSET(b,0,*nef,0.);

	  if(compressible==0){
      
	      /* calculate gamma (Ph.D. Thesis Jasak) */

	      FORTRAN(calcgammao,(nface,ielfa,vel,gradoel,gamma,xlet,xxn,xxj,
				  ipnei,&betam,nef,flux));

	      mafillomain(nef,ipnei,neifa,neiel,vfa,xxn,area,
	       au,&au[*nflnei],jq,irow,&nzs,b,vel,umfa,xlet,xle,gradofa,xxi,
	       body,volume,ielfa,lakonf,ifabou,nbody,nef,
	       &dtimef,velo,veloo,cvfa,hcfa,cvel,gradvel,xload,gamma,
	       xrlfa,xxj,nactdohinv,&a1,&a2,&a3,flux,iau6,xxni,xxnj,
	       iturbulent,f1,of2,gradkel,gradoel);
	  }else{

              /* underrelaxation of the temperature only for the
                 compressible simple scheme */

	      if(isimplec==0){
		  memcpy(&temp[6**nef],&vel[6**nef],sizeof(double)**nef);
	      }
 
              /* convective scheme */

//	      if(ischeme==1){
	          hro_udmain(nface,ielfa,vel,ipnei,nef,flux,vfa,&num_cpus);
/*	      }else{
		  FORTRAN(hro_mod_smart,(nface,ielfa,vel,gradtel,gamma,
                                    xlet,xxn,xxj,ipnei,&betam,nef,flux,vfa));
				    }*/

	      mafillocompmain(nef,ipnei,neifa,neiel,vfa,xxn,area,
	       au,&au[*nflnei],jq,irow,&nzs,b,vel,umfa,xlet,xle,gradofa,xxi,
	       body,volume,ielfa,lakonf,ifabou,nbody,nef,
	       &dtimef,velo,veloo,cvfa,hcfa,cvel,gradvel,xload,
	       xrlfa,xxj,nactdohinv,&a1,&a2,&a3,flux,iau6,xxni,xxnj,
	       iturbulent,f1,of2,gradkel,gradoel);
	  }

	  for(i=0;i<*nef;i++){rwork[i]=1./au[*nflnei+i];}

	  /* reordering the lhs (getting rid of zeros) */
	  
	  reorderlhsmain(au,am,iamorig,&nz_num,&num_cpus);

	  /* modifying the rhs (taking common contributions
	     between the cpu-blocks into account) */
	  
	  FORTRAN(reorderrhs,(au,&b[0],&vel[7**nef],neighblock,&nneighblock));

	  dgmresmain(nef,&b[0],&vel[7**nef],&nz_num,iam,jam,am,
			     &isym,&itol,&tol,&itmax,&iter,
			     &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
			     &ligw,rwork,iwork,nestart,&num_cpus);
	  if(ierr>1){
	      printf("*WARNING in compfluid: error message from predgmres (om)=%d\n",ierr);
	      if(ierr>ierrmax)ierrmax=ierr;
	  }

	  /* underrelaxation of the temperature only for compressible
	     simple scheme */

	  if((compressible==1)&&(isimplec==0)){
	      for(i=6**nef;i<7**nef;i++){vel[i]=0.8*vel[i]+0.2*temp[i];}
	  }

	  /* extrapolation of the turbulence variables at the element centers
	     to the face centers */

	  extrapol_kelmain(nface,ielfa,xrlfa,vel,vfa,
		  ifabou,xbounact,nef,gradkel,gradkfa,neifa,rf,area,volume,
                  xle,xxi,&icyclic,xxn,ipnei,ifatie,xlet,xxj,
                  coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
		  umfa,physcon,&iitg,c,&num_cpus);

	  extrapol_oelmain(nface,ielfa,xrlfa,vel,vfa,
		  ifabou,xbounact,nef,gradoel,gradofa,neifa,rf,area,volume,
                  xle,xxi,&icyclic,xxn,ipnei,ifatie,xlet,xxj,
                  coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
		  umfa,physcon,&iitg,c,dy,&num_cpus);

      }

      /* extrapolating the velocity from the elements centers to the face
	 centers, thereby taking the boundary conditions into account */


/*      printf("before extrapol_velmain2\n");
      for(i=0;i<*nef;i++){
	  printf("i=%d,co[0]=%f,vel[1]=%f,%f\n",i,coel[3*i+1],vel[*nef+i],vel[*nef+i]/coel[3*i+1]);
	  }*/
      extrapol_velmain(nface,ielfa,xrlfa,vel,vfa,
              ifabou,xbounact,ipnei,nef,&icyclic,c,ifatie,xxn,gradvel,
              gradvfa,neifa,rf,area,volume,xle,xxi,xxj,xlet,
	      coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
	      &iitg,&num_cpus);
/*      printf("after extrapol_velmain2\n");
      for(i=0;i<*nface;i++){
	  printf("i=%d,cofa[0]=%f,vfa[1]=%f,%f\n",i,cofa[3*i+1],vfa[8*i+1],vfa[8*i+1]/cofa[3*i+1]);
	  }*/

//      FORTRAN(writevfa,(vfa,nface,nactdohinv,ielfa));

      /* end subiterations */
      
/*      for(i=0;i<5;i++){velnorm[i]=0;}
      FORTRAN(norm,(vel,velnorm,nef));

      relnormt=0.;
      relnormv=0.;
      relnormp=0.;
      relnormmax=0.;

      if(*ithermal!=0){
	  if(velnorm[0]/(*nef)>1.e-10){
	      relnormt=fabs(velnorm[0]-velnormo[0])/(velnorm[0]);
	      if(relnormt>relnormmax) relnormmax=relnormt;
	  }
      }
      if((velnorm[1]+velnorm[2]+velnorm[3])/(*nef)>1.e-10){
	  relnormv=fabs(velnorm[1]+velnorm[2]+velnorm[3]-velnormo[1]-velnormo[2]-velnormo[3])/(velnorm[1]+velnorm[2]+velnorm[3]);
	  if(relnormv>relnormmax) relnormmax=relnormv;
      }
      if(velnorm[4]/(*nef)>1.e-10){
	  relnormp=fabs(velnorm[4]-velnormo[4])/(velnorm[4]);
	  if(relnormp>relnormmax) relnormmax=relnormp;
	  }*/
      relnormmax=relnormt;
      if(relnormv>relnormmax){relnormmax=relnormv;}
      if(relnormp>relnormmax){relnormmax=relnormp;}
      
      fprintf(f3,"%d %d %11.4e %11.4e %11.4e %11.4e\n",iincf,iittf,dtimef,relnormt,relnormv,relnormp);

//      memcpy(velnormo,velnorm,sizeof(double)*5);

      if(*nmethod==1){

          /* steady state flow:
	     calculate the velocity only once in each increment */

	  if(relnormmax<1.e-10) iconvergence=1;
//	  if(relnormmax<1.e-5) iconvergence=1;
	  if(compressible==0){
	      break;
	  }else if((relnormmax<1.e-3)||(iittf==iitt)){

              /* compressible steady state flow */

/*	      if(((iittf>2*iittfo)&&(iittf>10))||
	      (ierrmax>0)){*/
		  
		  /* divergence */
		  
/*		  timef-=dtimef;
		  memcpy(&vel[0],&velo[0],sizeof(double)*8**nef);
		  memcpy(&velo[0],&veloo[0],sizeof(double)*8**nef);
		  dtimef*=0.25;
		  printf("divergence: recalculated increment with reduced time increment  %e\n",dtimef);
		  a1=1./dtimef;
		  a2=-a1;
		  iittf=iittfo;
		  iincf--;
		  }else if((iittfo<2)&&(iittf<2)){*/
		  
		  /* fast convergence */
		  
/*		  dtimef*=1.01;
		  printf("increased time increment to %e\n",dtimef);
		  a1=1./dtimef;
		  a2=-a1;
	      }
	      iittfo=iittf;*/
	      break;
	  }
      }else{

          /* transient flow:
             calculate the velocity repeatedly in each increment */

	  if((relnormmax<1.e-5)||(iittf==iitt)){

              /* dynamic change of time increment for transient
                 compressible flow */

	      if(compressible!=0){

		  /* compressible flow */

		  if((iittfo<3)&&(iittf<3)){

		      /* fast convergence */

		      dtimef*=1.05;
		      printf("increased time increment to %e\n",dtimef);
		      a1=1./dtimef;
		      a2=-a1;
		  }else if(iittf>iitt-1){

		      /* divergence */

		      timef-=dtimef;
		      memcpy(&vel[0],&velo[0],sizeof(double)*8**nef);
		      memcpy(&velo[0],&veloo[0],sizeof(double)*8**nef);
		      dtimef*=0.25;
		      printf("divergence: recalculated increment with reduced time increment  %e\n",dtimef);
		      a1=1./dtimef;
		      a2=-a1;
//                  iittf=iittfo????
		  }else if((iittfo>10)&&(iittf>10)){

		      /* slow convergence */

		      dtimef*=0.95;
		      printf("decreased time increment to %e\n",dtimef);
		      a1=1./dtimef;
		      a2=-a1;
		  }
		  iittfo=iittf;
	      }
	      break;
	  }
      }

      }while(1);
      
      if(((iincf/jout[1])*jout[1]==iincf)||(iconvergence==1)||
	 (iincf==jmax[1])){

	  /* calculating the stress and the heat flow at the
             integration points, if requested */

	  if((strcmp1(&filab[3306],"SF  ")==0)||
             (strcmp1(&filab[3480],"SVF ")==0))isti=1;
          if(strcmp1(&filab[3393],"HFLF")==0)iqfx=1;
	  for(i=0;i<*nprint;i++){
	      if(strcmp1(&prlab[6*i],"SVF")==0) isti=1;
	      if(strcmp1(&prlab[6*i],"HFLF")==0)iqfx=1;
	  }

          /* calculating the heat conduction at the element centers */

	  if(iqfx==1){
	      NNEW(hcel,double,*nef);
	      FORTRAN(calchcel,(vel,cocon,ncocon,ielmatf,ntmat_,mi,
			       hcel,nef));
	  }

	  /* calculating the stress and/or the heat flux at the
             element centers */

	  if((isti==1)||(iqfx==1)){
	      FORTRAN(calcstressheatflux,(sti,umel,gradvel,qfx,hcel,
					  gradtel,nef,&isti,&iqfx,mi));
	      if(iqfx==1)SFREE(hcel);
	  }
 
          /* extrapolating the stresses */

	  if((strcmp1(&filab[3306],"SF  ")==0)||
             (strcmp1(&filab[3480],"SVF ")==0)){
	      nfield=6;
	      ndim=6;
	      if((*norien>0)&&
                 ((strcmp1(&filab[3311],"L")==0)||(strcmp1(&filab[3485],"L")==0))){
		  iorienglob=1;
	      }else{
		  iorienglob=0;
	      }
	      strcpy1(&cflag[0],&filab[2962],1);
	      NNEW(stn,double,6**nk);
	      FORTRAN(extrapolate,(sti,stn,ipkonf,inum,konf,lakonf,
		      &nfield,nk,nef,mi,&ndim,orab,ielorienf,co,&iorienglob,
		      cflag,vold,&force,ielmatf,thicke,ielpropf,prop));
	  }

	  /* extrapolating the heat flow */

	  
	  if(strcmp1(&filab[3393],"HFLF")==0){
	      nfield=3;
	      ndim=3;
	      if((*norien>0)&&(strcmp1(&filab[3398],"L")==0)){
		  iorienglob=1;
	      }else{
		  iorienglob=0;
	      }
	      strcpy1(&cflag[0],&filab[3049],1);
	      NNEW(qfn,double,3**nk);
	      FORTRAN(extrapolate,(qfx,qfn,ipkonf,inum,konf,lakonf,
		      &nfield,nk,nef,mi,&ndim,orab,ielorienf,co,&iorienglob,
		      cflag,vold,&force,ielmatf,thicke,ielpropf,prop));
	  }
	  
	  /* extrapolating the facial values of the static temperature 
             and/or the velocity and/or the static pressure to the nodes */

	  if(imach){NNEW(xmach,double,*nk);}
	  if(ikappa){NNEW(xkappa,double,*nk);}
	  if(iturb){NNEW(xturb,double,2**nk);}

/*	  printf("before extrapolatefluid: t\n");
	      for(i=0;i<*nface;i++){
		  printf("i=%d,cofa[0]=%f,vfa[0]=%f %f\n",i,cofa[3*i+1],vfa[8*i],vfa[8*i]/(1.+(1.-cofa[3*i+1])*cofa[3*i+1]/2.));
		  }*/

	  FORTRAN(extrapolatefluid,(nk,iponofa,inofa,inum,vfa,vold,ielfa,
                  ithermal,&imach,&ikappa,xmach,xkappa,shcon,nshcon,ntmat_,
		  ielmatf,physcon,mi,&iturb,xturb,gradtfa,gradvfa,gradpfa,
		  gradkfa,gradofa,co,cofa,ifabou));

/*	  printf("after extrapolatefluid: t\n");
	      for(i=0;i<*nk;i++){
		  printf("i=%d,co[0]=%f,vold[0]=%f %f\n",i,co[3*i+1],vold[5*i],vold[5*i]/(1.+(1.-co[3*i+1])*co[3*i+1]/2.));
		  }*/

          /* storing the results in dat-format */

	  ptimef=ttimef+timef;
	  FORTRAN(printoutfluid,(set,nset,istartset,iendset,ialset,nprint,
				 prlab,prset,ipkonf,lakonf,sti,eei,
                                 xstate,ener,mi,nstate_,co,konf,qfx,
                                 &ptimef,trab,inotr,ntrans,orab,ielorienf,
                                 norien,vold,ielmatf,
                                 thicke,eme,xturb,physcon,nactdoh,
                                 ielpropf,prop,xkappa,xmach,ithermal,
                                 orname));

          /* thermal flux and drag: storage in dat-format */

	  FORTRAN(printoutface,(co,rhcon,nrhcon,ntmat_,vold,shcon,nshcon,
		  cocon,ncocon,&compressible,istartset,iendset,ipkonf,
		  lakonf,konf,
		  ialset,prset,&ptimef,nset,set,nprint,prlab,ielmatf,mi,
		  ithermal,nactdoh,&icfd,time,stn));
	  
	  /* storing the results in frd-format */
	  

/*      printf("before frdfluid: vx mix\n");
      for(i=0;i<*nef;i++){
	  printf("i=%d,co[0]=%f,vel[1]=%f,%f\n",i,coel[3*i+1],vel[*nef+i],vel[*nef+i]/coel[3*i+1]);
	  }*/
	  FORTRAN(frdfluid,(co,nk,konf,ipkonf,lakonf,nef,vold,kode,&timef,ielmatf,
			    matname,filab,inum,ntrans,inotr,trab,mi,istep,
                            stn,qfn,nactdohinv,xmach,xkappa,physcon,xturb,
                            coel,vel,cofa,vfa,nface));

//	  FORTRAN(writevfa,(vfa,nface,nactdohinv,ielfa));

	  if((strcmp1(&filab[3306],"SF  ")==0)||
             (strcmp1(&filab[3480],"SVF ")==0)){SFREE(stn);}
	  if(strcmp1(&filab[3393],"HFLF")==0){SFREE(qfn);}

	  if(imach){SFREE(xmach);}
	  if(ikappa){SFREE(xkappa);}
	  if(iturb){SFREE(xturb);}

      }
      
      if(iincf==jmax[1]){
	  printf("*INFO: maximum number of fluid increments reached\n\n");
	  fclose(f3);
	  break;
//	  FORTRAN(stop,());
      }
      if(last==1){
	  printf("*INFO: mechanical time increment reached: time=%e\n\n",*dtime);
	  fclose(f3);
//	  FORTRAN(stop,());
	  break;
      }
      if(iconvergence==1){
	  printf("*INFO: steady state reached\n\n");
	  fclose(f3);
//	  FORTRAN(stop,());
	  break;
      }
      
      
      memcpy(&veloo[0],&velo[0],sizeof(double)*8**nef);
      memcpy(&velo[0],&vel[0],sizeof(double)*8**nef);
      
  }while(1);
  
  FORTRAN(closefilefluid,());

  SFREE(flux);

  if(compressible==0){SFREE(ia);SFREE(ja);SFREE(aua);}

  SFREE(irow);SFREE(icol);SFREE(jq);SFREE(iau6);
  
  SFREE(coel);SFREE(cosa);SFREE(xxn);SFREE(xxi);SFREE(xle);SFREE(xlen);
  SFREE(xlet);SFREE(cofa);SFREE(area);SFREE(xrlfa);SFREE(volume);
  SFREE(cosb);SFREE(xxni);SFREE(xxnj);SFREE(xxicn);SFREE(xxj);
  SFREE(rf);
  if(*iturbulent>0){
      SFREE(dy);
      if(*iturbulent>2) SFREE(yy);
  }

  SFREE(ifabou);SFREE(umfa);SFREE(umel);

  SFREE(gradvel);SFREE(gradvfa);SFREE(au);SFREE(ad);SFREE(b);SFREE(advfa);
  SFREE(ap);SFREE(bp);SFREE(gradpel);SFREE(rwork);SFREE(gradpfa);
  SFREE(gradpcel);SFREE(gradpcfa);
  SFREE(hfa);SFREE(hel);SFREE(adv);SFREE(bv);SFREE(sel);
  SFREE(auv);

  if(*ithermal>0){
      SFREE(gradtel);SFREE(gradtfa);SFREE(hcfa);SFREE(cvel);SFREE(cvfa);
  }

  if(*iturbulent>0){
      SFREE(gradkel);SFREE(gradkfa);SFREE(gradoel);SFREE(gradofa);
      if(*iturbulent>2){
	  SFREE(f1);
	  if(*iturbulent>3) SFREE(of2);
      }
  }

  SFREE(inum);SFREE(v);

  SFREE(iponofa);SFREE(inofa);

  if(*nbody>0) SFREE(body);

  *ithermal=ithermalref;

  SFREE(temp);SFREE(gamma);

  SFREE(iam);SFREE(jam);SFREE(iamorig);SFREE(am);
  SFREE(nestart);SFREE(ineighblock);SFREE(neighblock);

  return;
  
} 
