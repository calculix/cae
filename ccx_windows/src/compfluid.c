/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2020 Guido Dhondt                     */

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

void compfluid(double **cop,ITG *nk,ITG **ipkonfp,ITG *konf,char **lakonfp,
	       char **sidefacep,ITG *ifreestream,ITG *nfreestream,
	       ITG *isolidsurf,ITG *neighsolidsurf,ITG *nsolidsurf,
	       ITG *nshcon,double *shcon,ITG *nrhcon,double *rhcon,
	       double **voldp,ITG *ntmat_,ITG *nodeboun,ITG *ndirboun,
	       ITG *nboun,ITG *ipompc,ITG *nodempc,ITG *nmpc,ITG *ikmpc,
	       ITG *ilmpc,ITG *ithermal,ITG *ikboun,ITG *ilboun,
	       ITG *iturbulent,ITG *isolver,ITG *iexpl,double *ttime,
	       double *time,double *dtime,ITG *nodeforc,ITG *ndirforc,
	       double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,double *xbody,ITG *ipobody,ITG *nbody,
	       ITG *ielmatf,char *matname,ITG *mi,ITG *ncmat_,double *physcon,
	       ITG *istep,ITG *iinc,ITG *ibody,double *xloadold,double *xboun,
	       double *coefmpc,ITG *nmethod,double *xforcold,double *xforcact,
	       ITG *iamforc,ITG *iamload,double *xbodyold,double *xbodyact,
	       double *t1old,double *t1,double *t1act,ITG *iamt1,double *amta,
	       ITG *namta,ITG *nam,double *ampli,double *xbounold,
	       double *xbounact,ITG *iamboun,ITG *itg,ITG *ntg,char *amname,
	       double *t0,ITG **nelemfacep,ITG *nface,double *cocon,
	       ITG *ncocon,double *xloadact,double *tper,ITG *jmax,ITG *jout,
	       char *set,ITG *nset,ITG *istartset,ITG *iendset,ITG *ialset,
	       char *prset,char *prlab,ITG *nprint,double *trab,ITG *inotr,
	       ITG *ntrans,char *filab,char *labmpc,double *sti,ITG *norien,
	       double *orab,char *jobnamef,char *tieset,ITG *ntie,ITG *mcs,
	       ITG *ics,double *cs,ITG *nkon,ITG *mpcfree,ITG *memmpc_,
	       double *fmpc,ITG *nef,ITG **inomatp,double *qfx,ITG *neifa,
	       ITG *neiel,ITG *ielfa,ITG *ifaext,double *vfa,double *vel,
	       ITG *ipnei,ITG *nflnei,ITG *nfaext,char *typeboun,ITG *neij,
	       double *tincf,ITG *nactdoh,ITG *nactdohinv,ITG *ielorienf,
	       char*jobnamec,ITG *ifatie,ITG *nstate_,double *xstate,
	       char *orname,ITG *kon,double *ctrl,ITG *kode,double *velo,
	       double *veloo,ITG *initial){

  /* main computational fluid dynamics routine */
  
  char cflag[1],*lakonf=NULL,fncvg[132]="",*lakon=NULL;

  ITG *ipointer=NULL,*mast1=NULL,*irow=NULL,*icol=NULL,*jq=NULL,nzs=20000000,
    compressible,*ifabou=NULL,*ja=NULL,*ikf=NULL,nfabou,im,iflag,*ipkon=NULL,
    *ielprop=NULL,*ielmat=NULL,*ipkonf=NULL,last=0,icyclic,*iau6=NULL,
    ifixdtimef=0,ithermalref,*integerglob=NULL,iincf,ipower=64,*konl=NULL,
    iconvergence=0,i,*inum=NULL,j,k,ifreefa,isiz,*ipofano=NULL,*ifano=NULL,
    *ia=NULL,*ielpropf=NULL,icent=0,isti=0,iqfx=0,nfield,ndim,iorienglob,
    force=0,icfd=1,imach=0,ikappa=0,iatleastonepressurebc,iturb=0,*inoel=NULL,
    *iponoel=NULL,icounter,ischeme=1,isimplec=0,iitf,iitg,iitp,*iam=NULL,
    *jam=NULL,*iamorig=NULL,nz_num,*nestart=NULL,*ineighblock=NULL,
    *neighblock=NULL,nneighblock,iit,iito,ip,ierrmax,ncfd,iitpt,*inlet=NULL,
    *ipgradfa=NULL,*ipbount=NULL,*ipbounv1=NULL,*ipbounv2=NULL,
    *ipbounv3=NULL,*ipbounp=NULL,*ibount=NULL,*ibounv1=NULL,*ibounv2=NULL,
    *ibounv3=NULL,*ibounp=NULL,nbt,nbv1,nbv2,nbv3,nbp,nbpt,nbpv1,nbpv2,nbpv3,
    nbpp,nkf,*iponofa=NULL,*inofa=NULL,symmetryflag=2,inputformat=4;

  ITG nelt,isym,itol,itmax,iunit,lrgw,*igwk=NULL,ligw,ierr,*iwork=NULL,iter,
    nsave,lenw,leniw;

  double *umfa=NULL,reltime,*doubleglob=NULL,dtimefold,
    *co=NULL,*vold=NULL,*coel=NULL,*cosa=NULL,*gradvel=NULL,*gradvfa=NULL,
    *xxn=NULL,*xxi=NULL,*xle=NULL,*xlen=NULL,*xlet=NULL,timef,dtimef,
    *cofa=NULL,*area=NULL,*xrlfa=NULL,reltimef,ttimef,*hcfa=NULL,*cvel=NULL,
    *au=NULL,*ad=NULL,*b=NULL,*volume=NULL,*body=NULL,*dy=NULL,
    *advfa=NULL,*ap=NULL,*bp=NULL,*xxj=NULL,*gradkel=NULL,*gradoel=NULL,
    *cosb=NULL,hmin,tincfguess,*h=NULL,*gradelsh=NULL,
    *hel=NULL,*hfa=NULL,*auv=NULL,*adv=NULL,*bv=NULL,*sel=NULL,*gamma=NULL,
    *gradtfa=NULL,*gradtel=NULL,*umel=NULL,*cvfa=NULL,*gradpel=NULL,
    *eei=NULL,*ener=NULL,*thicke=NULL,*eme=NULL,c[9],*gradkfa=NULL,
    ptimef,*stn=NULL,*qfn=NULL,*hcel=NULL,*aua=NULL,a1,a2,a3,beta=0.,
    *prop=NULL,*xxni=NULL,*xxnj=NULL,*xxicn=NULL,*xturb=NULL,
    *xmach=NULL,*xkappa=NULL,*flux=NULL,*sc=NULL,*gradfash=NULL,
    relnormt,relnormv,relnormp=0,relnormmax=1.e30,*temp=NULL,*yy=NULL,
    *gradofa=NULL,betam=0.1,*gradpfa=NULL,*gradpcel=NULL,*gradpcfa=NULL,
    *am=NULL,*f1=NULL,*of2=NULL,*xxna=NULL,*ale=NULL,*alet=NULL,
    *ratio=NULL,dgmrestol=1.e-12,*bvcp=NULL,*bcp=NULL,*xbount=NULL,
    *xbounv1=NULL,*xbounv2=NULL,*xbounv3=NULL,*xbounp=NULL,*adb=NULL,
    *aub=NULL,sigma=0;

  double tol,*rgwk=NULL,err,*sb=NULL,*sx=NULL,*rwork=NULL,*rf=NULL;

  FILE *f3;

  co=*cop;ipkonf=*ipkonfp;lakonf=*lakonfp;vold=*voldp;

#ifdef SGI
  ITG token;
#endif
  
  ncfd=(ITG)physcon[9];
  //  printf("ncfd=%d\n",ncfd);
	  
  strcpy(fncvg,jobnamec);
  strcat(fncvg,".fcv");

  if((f3=fopen(fncvg,"w"))==NULL){
    //  if((f3=fopen("fluidconvergence","w"))==NULL){
    printf("*ERROR in compfluid: cannot open cvg file for writing...");
    exit(0);
  }
  fprintf(f3,"temperature    velocity    pressure\n\n");

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
    
  /* local declaration prevails,if strictly positive */
    
  envloc = getenv("CCX_NPROC_CFD");
  if(envloc){
    num_cpus=atoi(envloc);
    if(num_cpus<0){
      num_cpus=0;
    }else if(num_cpus>sys_cpus){
      num_cpus=sys_cpus;
    }
  }
    
  /* else global declaration,if any,applies */
    
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
    
  printf(" Using up to %" ITGFORMAT " cpu(s) for CFD.\n",num_cpus);
    
  //  pthread_t tid[num_cpus];
    
  /*  *iexpl==0:  structure:implicit,fluid:incompressible
   *iexpl==1:  structure:implicit,fluid:compressible
   *iexpl==2:  structure:explicit,fluid:incompressible
   *iexpl==3:  structure:explicit,fluid:compressible */
    
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
    
  /* iitf: max number of transient iterations
     iitg: max number of geometric iterations (extrapol_*.f)
     iitp: max number of pressure iterations */
    
  iitf=(ITG)floor(ctrl[49]);
  iitg=(ITG)floor(ctrl[50]);
  iitp=(ITG)floor(ctrl[51]);
  iitpt=(ITG)floor(ctrl[52]);
    
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
  NNEW(h,double,*nflnei);
  NNEW(sc,double,*nef);
  DMEMSET(sc,0,*nef,1.);
  NNEW(cosa,double,*nflnei);
  NNEW(cosb,double,*nflnei);
  NNEW(xxn,double,3**nflnei);
  NNEW(xxna,double,3**nflnei);
  NNEW(xxi,double,3**nflnei);
  NNEW(xxj,double,3**nflnei);
  NNEW(xxni,double,3**nflnei);
  NNEW(xxicn,double,3**nflnei);
  NNEW(xxnj,double,3**nflnei);
  NNEW(xle,double,*nflnei);
  NNEW(ale,double,*nflnei);
  NNEW(xlen,double,*nflnei);
  NNEW(xlet,double,*nflnei);
  NNEW(alet,double,*nflnei);
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
    
  FORTRAN(geocfd,(nef,ipkonf,konf,lakonf,co,coel,cofa,nface,
		  ielfa,area,ipnei,neiel,xxn,xxi,xle,xlen,xlet,xrlfa,cosa,
		  volume,neifa,xxj,cosb,&hmin,ifatie,cs,tieset,&icyclic,c,
		  neij,physcon,isolidsurf,nsolidsurf,dy,xxni,xxnj,xxicn,
		  nflnei,iturbulent,rf,yy,vel,velo,veloo,xxna,ale,alet,h));

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
    
  NNEW(inlet,ITG,*nface);
  NNEW(ifabou,ITG,7**nfaext);
  FORTRAN(applyboun,(ifaext,nfaext,ielfa,ikboun,ilboun,
		     nboun,typeboun,nelemload,nload,sideload,isolidsurf,
		     nsolidsurf,ifabou,&nfabou,nface,nodeboun,ndirboun,ikmpc,
		     ilmpc,labmpc,nmpc,nactdohinv,&compressible,
		     &iatleastonepressurebc,ipkonf,kon,konf,inlet));
  RENEW(ifabou,ITG,nfabou);
    
  /* catalogueing the nodes for output purposes (interpolation at
     the nodes */
    
  NNEW(ipofano,ITG,*nk);
  NNEW(ifano,ITG,2**nface*4);
    
  FORTRAN(cataloguenodes,(ipofano,ifano,&ifreefa,ielfa,ifabou,ipkonf,
			  konf,lakonf,nface,nk));
    
  RENEW(ifano,ITG,2*ifreefa);

  /* determining the coefficients to get from the element center values:
     - the nodal values (by extra/interpolation of the element center values)
     - the facial values (by interpolation from the nodal values)
     - the center gradient values (by use of the shape functions using the
     nodal values)
     - the facial gradient values (by use of the shape functions using the
     nodal values) */
  
  /* elem2node(coel,nef,&ncfd,&num_cpus,ipkonf,konf,lakonf,co,&hmin,ipnei,neiel,
	    xxn,nkon,nk,neifa,area,&ikf,&konl,&ratio,&gradelsh,&gradfash,nflnei,
	    &ipgradfa,&ipbount,&ipbounv1,&ipbounv2,&ipbounv3,&ipbounp,
	    &ibount,&ibounv1,&ibounv2,&ibounv3,&ibounp,&xbount,&xbounv1,
	    &xbounv2,&xbounv3,&xbounp,&nbt,&nbv1,&nbv2,&nbv3,&nbp,ithermal,
	    ipofano,ifano,&nbpt,&nbpv1,&nbpv2,&nbpv3,&nbpp,ielfa,ifabou,
	    &ifreefa,&nkf,&iponofa,&inofa,nface);*/
  
  //  SFREE(ikf);SFREE(konl);SFREE(ratio);
  //  SFREE(gradfash);SFREE(gradelsh);SFREE(ipgradfa);
    
  /* material properties for athermal calculations 
     = calculation for which no initial thermal conditions
     were defined */
    
  NNEW(umfa,double,*nface);
  NNEW(umel,double,*nef);
    
  //  if((*ithermal==0)||(*iturbulent>0)){
  if(*ithermal==0){
	
    /* athermal incompressible calculations */
	
    /* calculating the dynamic viscosity at the element centers */
	
    FORTRAN(calcumel,(nef,vel,shcon,nshcon,ielmatf,ntmat_,
		      ithermal,mi,umel));
	
  }
    
    
  if(*ithermal>0){
    NNEW(hcfa,double,*nface);
    NNEW(cvel,double,*nef);
    NNEW(cvfa,double,*nface);
  }
    
  if(*nbody>0) NNEW(body,double,4**nef);
    
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
	co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
	ntrans,trab,inotr,vold,integerglob,doubleglob,tieset,istartset,
	iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,
	ipobody,iponoel,inoel,ipkon,kon,ielprop,prop,ielmat,
	shcon,nshcon,rhcon,nrhcon,cocon,ncocon,ntmat_,lakon));
	
    /* body forces (gravity,centrifugal and Coriolis forces */
	
    if(*nbody>0){
      FORTRAN(inicalcbody,(nef,body,ipobody,ibody,xbody,coel,vel,lakonf,
			   nactdohinv,&icent));
    }
  }
    
  /* extrapolating the velocity from the elements centers to the face
     centers,thereby taking the boundary conditions into account */
    
  NNEW(gradvel,double,9**nef);
  NNEW(gradvfa,double,9**nface);
    
  extrapol_velmain(nface,ielfa,xrlfa,vel,vfa,
		   ifabou,xbounact,ipnei,nef,&icyclic,c,ifatie,xxn,gradvel,
		   gradvfa,neifa,rf,area,volume,xle,xxi,xxj,xlet,
		   coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
		   &iitg,&num_cpus,&compressible,xxna,&ncfd,cofa);
    
  /* extrapolation of the pressure at the element centers
     to the face centers */
    
  NNEW(gradpel,double,3**nef);
  NNEW(gradpfa,double,3**nface);
    
  extrapol_pelmain(nface,ielfa,xrlfa,vel,vfa,
		   ifabou,xbounact,nef,gradpel,gradpfa,neifa,rf,area,volume,
		   xle,xxi,&icyclic,xxn,ipnei,ifatie,
		   coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
		   &iitg,c,&num_cpus,&compressible,xxna,&ncfd);
    
  /* generate fields for the pressure correction gradients */
    
  NNEW(gradpcel,double,3**nef);
  NNEW(gradpcfa,double,3**nface);
    
  /* extrapolation of the temperature at the element centers
     to the face centers */
    
  if(*ithermal>0){
	
    NNEW(gradtel,double,3**nef);
    NNEW(gradtfa,double,3**nface);
	
    extrapol_telmain(nface,ielfa,xrlfa,vel,vfa,
		     ifabou,xbounact,nef,gradtel,gradtfa,neifa,rf,area,volume,
		     xle,xxi,&icyclic,xxn,ipnei,ifatie,xload,xlet,xxj,
		     coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
		     &iitg,c,&num_cpus,&compressible,xxna,&ncfd);
	
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
	
    /* calculating the density at the face centers */
	
    if(ischeme==1){
      hrr_udmain(nface,vfa,shcon,ielmatf,ntmat_,
		 mi,ielfa,ipnei,vel,nef,flux,
		 &num_cpus,xxi,xle,gradpel,gradtel,
		 neij);
    }else{
      hrr_mod_smartmain(nface,vfa,shcon,ielmatf,ntmat_,
			mi,ielfa,ipnei,vel,nef,flux,
			gradpel,gradtel,xxj,xlet,
			&num_cpus);
    }
	
  }
    
  /* calculating the initial mass flux */
    
  FORTRAN(calcinitialflux,(area,vfa,xxna,ipnei,nef,neifa,lakonf,flux));

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
		     umfa,physcon,&iitg,c,&num_cpus,&compressible,xxna,&ncfd,
		     inlet);

    extrapol_oelmain(nface,ielfa,xrlfa,vel,vfa,
		     ifabou,xbounact,nef,gradoel,gradofa,neifa,rf,area,volume,
		     xle,xxi,&icyclic,xxn,ipnei,ifatie,xlet,xxj,
		     coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
		     umfa,physcon,&iitg,c,dy,&num_cpus,&compressible,xxna,&ncfd,
		     inlet);

    if(*iturbulent>2){
      NNEW(f1,double,*nef);
      if(*iturbulent>3) NNEW(of2,double,*nef);
    }
	
  }

  /* calculating the time increment */
    
  FORTRAN(initincf,(nface,&hmin,vfa,umfa,cvfa,hcfa,ithermal,&dtimef,
		    &compressible));
    
  /* start of the major loop */
    
  NNEW(advfa,double,*nface);
  NNEW(hfa,double,3**nface);
    
  NNEW(ap,double,*nface);
  NNEW(bp,double,*nface);
    
  NNEW(au,double,*nflnei+*nef);
  NNEW(ad,double,*nef);
  NNEW(b,double,*nef);
  NNEW(bcp,double,*nef);
    
  NNEW(auv,double,*nflnei+*nef);
    
  NNEW(bv,double,3**nef);
  NNEW(bvcp,double,3**nef);
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
    
  if((strcmp1(&filab[2088],"TURB")==0)&&(*iturbulent>0)){
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

  /* if the user has specified a fixed fluid time increment,use it */
    
  if(*tincf>0.){
    dtimef=*tincf;
    ifixdtimef=1;
  }
    
  printf("time increment for the CFD-calculations = %e\n\n",dtimef);
    
  ttimef=*ttime;
  timef=*time-*dtime;
    
  if(compressible==0){
    a1=1.5/dtimef;
    a2=-2./dtimef;
    a3=0.5/dtimef;
  }else{
    a1=1./dtimef;
    a2=-a1;
    a3=0.;
  }
    
  iito=iitf;
    
  NNEW(temp,double,8**nef);
  NNEW(gamma,double,*nface);
    
  icounter=0;
    
  do{
	
    iincf++;
	
    //        printf("fluid increment = %d\n",iincf);
	
    if((iincf/ipower)*ipower==iincf){
      printf("fluid increment = %d\n",iincf);
	    
      /* reevaluating the time increment size 
	 only for steady state compressible calculations*/
	    
      if((*nmethod==1)&&(ifixdtimef==0)){
		
	if(*ithermal>0){
	  NNEW(hcel,double,*nef);
	  FORTRAN(calchcel,(vel,cocon,ncocon,ielmatf,ntmat_,mi,
			    hcel,nef));
	}
		
	dtimefold=dtimef;
	FORTRAN(newtincf,(ithermal,&dtimef,&compressible,vel,
			  hcel,umel,cvel,h,sc,iturbulent,ipkonf,
			  nmethod,nef,lakonf,xxn,ipnei));
		
	if(compressible==1){
	  a1=1./dtimef;
	  a2=-a1;
	}else{
	  beta=dtimef/dtimefold;
	  a1=(2.+beta)/(1.+beta)/dtimef;
	  a2=-(1.+beta)/beta/dtimef;
	  a3=1./(beta*(1.+beta))/dtimef;
	}
	if(*ithermal>0) SFREE(hcel);
      }
      ipower*=2;
    }
	
    ierrmax=0;
	
    timef+=dtimef;
    if((*time<timef)&&(*nmethod==4)){
      dtimefold=dtimef;
      dtimef-=(timef-*time);
      timef=*time;
      last=1;
      beta=dtimef/dtimefold;
      a1=(2.+beta)/(1.+beta)/dtimef;
      a2=-(1.+beta)/beta/dtimef;
      a3=1./(beta*(1.+beta))/dtimef;
    }
	
    /* starting iterations till convergence of the fluid increment */
	
    iit=0;
	
    do{
      iit++;
	    
      //	  printf("      iteration = %d\n",iit);
	    
      /* conditions for transient calculations */
	    
      if(*nmethod==4){
		
	/* boundary conditions at end of fluid increment */
		
	FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
			  xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
			  xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
			  namta,nam,ampli,&timef,&reltimef,&ttimef,&dtimef,ithermal,nmethod,
			  xbounold,xboun,xbounact,iamboun,nboun,
			  nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
			  co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
			  ntrans,trab,inotr,vold,integerglob,doubleglob,tieset,istartset,
			  iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,
			  ipobody,iponoel,inoel,ipkon,kon,ielprop,prop,ielmat,
			  shcon,nshcon,rhcon,nrhcon,cocon,ncocon,ntmat_,lakon));
		
	/* body forces (gravity,centrifugal and Coriolis forces) */
		
	if(*nbody>0){
	  FORTRAN(calcbody,(nef,body,ipobody,ibody,xbody,coel,vel,lakonf,
			    nactdohinv));
	}

      }else if(icent==1){
		
	/* body forces (gravity,centrifugal and Coriolis forces;
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
		    alet,ale,gradvfa,xxi,
		    body,volume,ielfa,lakonf,ifabou,nbody,
		    &dtimef,velo,veloo,sel,xrlfa,gamma,xxj,nactdohinv,&a1,
		    &a2,&a3,flux,&icyclic,c,ifatie,iau6,xxna,xxnj,
		    iturbulent,gradvel,of2,yy,umel,&ncfd,inlet,sc);
		
      }else{
		
	/* convection scheme */
		
	if(ischeme==1){
	  hrv_udmain(nface,ielfa,vel,ipnei,nef,flux,vfa,&num_cpus,
		     xxi,xle,gradvel,neij);
	}else{
	  hrv_mod_smartmain(nface,ielfa,vel,gradvel,xlet,xxj,ipnei,
			    nef,flux,vfa,&num_cpus);
	}

	//	printf("mafillvcompmain\n");
	mafillvcompmain(nef,ipnei,neifa,neiel,vfa,xxn,area,
			auv,&auv[*nflnei],jq,irow,&nzs,bv,vel,cosa,umfa,
			alet,ale,gradvfa,xxi,
			body,volume,ielfa,lakonf,ifabou,nbody,
			&dtimef,velo,veloo,sel,xrlfa,gamma,xxj,nactdohinv,&a1,
			&a2,&a3,flux,&icyclic,c,ifatie,iau6,xxna,xxnj,
			iturbulent,gradvel,of2,yy,umel,&ncfd,inlet,sc);
      }
	    
      for(i=0;i<*nef;i++){rwork[i]=1./auv[*nflnei+i];}
	    
      /* underrelaxation of the velocity only for compressible
	 simple scheme */
	    
      if((compressible==1)&&(isimplec==0)){
	memcpy(&temp[*nef],&vel[*nef],sizeof(double)*ncfd**nef);
      }
	    
      /* reordering the lhs (getting rid of zeros)  */
	    
      reorderlhsmain(auv,am,iamorig,&nz_num,&num_cpus);
	    
      /* modifying the rhs (taking common contributions
	 between the cpu-blocks into account) */
	    
      memcpy(&bvcp[0],&bv[0],sizeof(double)*3**nef);
      for(i=0;i<ncfd;i++){
	FORTRAN(reorderrhs,(auv,&bv[i**nef],&vel[(i+1)**nef],neighblock,&nneighblock));
      }
	    
      /* calculation of the momentum residual */
	    
      calcresvfluidmain(nef,am,&bv[0],&auv[*nflnei],iam,jam,
			&vel[*nef],&relnormv,nestart,&num_cpus,
                        &ncfd);
	    
      for(k=1;k<num_cpus+1;k++){
	if(k>1){
	  memcpy(&bv[0],&bvcp[0],sizeof(double)*3**nef);
	}
	for(i=0;i<ncfd;i++){
	  if(k>1){
	    FORTRAN(reorderrhs,(auv,&bv[i**nef],&vel[(i+1)**nef],neighblock,
				&nneighblock));
	  }
	  dgmresmain(nef,&bv[i**nef],&vel[(i+1)**nef],&nz_num,iam,jam,am,
		     &isym,&itol,&tol,&itmax,&iter,
		     &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
		     &ligw,rwork,iwork,nestart,&num_cpus,&dgmrestol);
	  
	  if(ierr>1){
	    printf("*WARNING in compfluid: error message from predgmres (v_%d)=%d\n",i+1,ierr);
	    if(ierr>ierrmax)ierrmax=ierr;
	  }
	}
      }
      
      /*      for(i=0;i<ncfd;i++){
	spooles(&auv[*nflnei],am,adb,aub,&sigma,&bv[i**nef],iam,jam,
		nef,&nz_num,&symmetryflag,&inputformat,&nz_num);
	memcpy(&vel[(i+1)**nef],&bv[i**nef],sizeof(double)**nef);
	}*/
	    
      /* underrelaxation of the velocity only for compressible
	 simple scheme */
	    
      if((compressible==1)&&(isimplec==0)){
	for(i=*nef;i<(ncfd+1)**nef;i++){vel[i]=0.8*vel[i]+0.2*temp[i];}
      }

      //      printf("iitpt=%d,iitp=%d\n",iitpt,iitp);
      
      for(j=0;j<iitpt;j++){
	//	printf("j=%d,iitpt=%d\n",j,iitpt);
		
	/* generating vol/ad and v* at the face centers (advfa and hfa) */
		
	if((compressible==0)||(isimplec==0)){
	  extrapolate_d_v_simplemain(nface,ielfa,xrlfa,&auv[*nflnei],advfa,hfa,
				     &icyclic,c,ifatie,vel,nef,volume,&num_cpus,
				     &ncfd);
	}else{
	  extrapolate_d_v_simplecmain(nface,ielfa,xrlfa,&auv[*nflnei],advfa,hfa,
				      &icyclic,c,ifatie,vel,nef,volume,auv,
				      ipnei,&num_cpus,&ncfd);
	}
		
	/* calculate the mass flow based on the newly calculated
	   velocity */
		
	calcfluxmain(area,vfa,xxna,ipnei,nef,neifa,flux,xxj,gradpfa,xlet,xle,vel,
		     advfa,ielfa,neiel,ifabou,hfa,&num_cpus);
		
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
			      ifabou,xbounact,nef,gradpcel,gradpcfa,neifa,rf,
			      area,volume,
			      xle,xxi,&icyclic,xxn,ipnei,ifatie,
			      coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,
			      nactdoh,&iflag,xxj,xlet,c,&num_cpus,xxna,&ncfd,
			      &ip);
			
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
		    
	  /* calculate the p' gradient at the 
	     element centers */

	  /* 16.04.2020: should be iflag=1????? */

	  
	  iflag=0;
	  extrapol_dpelmain(nface,ielfa,xrlfa,vel,vfa,
			    ifabou,xbounact,nef,gradpcel,gradpcfa,neifa,rf,area,volume,
			    xle,xxi,&icyclic,xxn,ipnei,ifatie,
			    coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,
			    nactdoh,&iflag,xxj,xlet,c,&num_cpus,xxna,&ncfd,&ip);
		    
	  /* correction of the velocity v* at the element centers due
	     to the pressure change in order to get v** */
		    
	  correctvelmain(&auv[*nflnei],nef,volume,gradpcel,vel,&num_cpus);
		    
	  /* extrapolating the velocity from the elements centers to the face
	     centers,thereby taking the boundary conditions into account */
		    
	  extrapol_velmain(nface,ielfa,xrlfa,vel,vfa,
			   ifabou,xbounact,ipnei,nef,&icyclic,c,ifatie,xxn,gradvel,
			   gradvfa,neifa,rf,area,volume,xle,xxi,xxj,xlet,
			   coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
			   &iitg,&num_cpus,&compressible,xxna,&ncfd,cofa);
		    
	  /* correcting the flux to get mf** */
		    
	  correctfluxmain(nef,ipnei,neifa,neiel,flux,vfa,advfa,area,
			  vel,alet,ielfa,ale,ifabou,&num_cpus,xxnj,
			  gradpcfa);
		    
	  /* correcting the pressure to get p* */
		    
	  /* underrelaxation always except for the compressible simplec
	     scheme */
		    
	  for(i=4**nef;i<5**nef;i++){vel[i]=0.2*vel[i]+temp[i];}
		    
	  extrapol_pelmain(nface,ielfa,xrlfa,vel,vfa,
			   ifabou,xbounact,nef,gradpel,gradpfa,neifa,rf,area,volume,
			   xle,xxi,&icyclic,xxn,ipnei,ifatie,
			   coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
			   &iitg,c,&num_cpus,&compressible,xxna,&ncfd);
	      
	}else{
		    
	  /* compressible media */
		    
	  /* temporarily store the pressure in temp */
		    
	  memcpy(&temp[4**nef],&vel[4**nef],sizeof(double)**nef);
		    
	  DMEMSET(au,0,*nflnei+*nef,0.);
	  DMEMSET(b,0,*nef,0.);
		    
	  //	printf("mafillpcompmain\n");
	  mafillpcompmain(nef,lakonf,ipnei,neifa,neiel,vfa,area,
			  advfa,xlet,cosa,volume,au,&au[*nflnei],jq,irow,ap,
			  ielfa,ifabou,xle,b,xxn,nef,
			  &nzs,hfa,gradpel,bp,xxi,neij,xlen,cosb,
			  ielmatf,mi,&a1,&a2,&a3,velo,veloo,&dtimef,shcon,
			  ntmat_,vel,nactdohinv,xrlfa,flux,iau6,xxicn,
			  gamma,inlet);
		    
	  for(i=0;i<*nef;i++){rwork[i]=1./au[*nflnei+i];}
		    
	  /* initial guess: 0 */
		    
	  DMEMSET(vel,4**nef,5**nef,0.);
		    
	  /* reordering the lhs (getting rid of zeros) */
		    
	  reorderlhsmain(au,am,iamorig,&nz_num,&num_cpus);
		    
	  /* calculation of the mass conservation residual */
		    
	  FORTRAN(calcrespfluid,(nef,&b[0],&au[*nflnei],
				 &temp[4**nef],&relnormp));
		    
	  /* no change of rhs (reorderrhs) since initial guess is zero */
		    
	  for(k=1;k<num_cpus+1;k++){
	    if(k==1){
	      memcpy(&bcp[0],&b[0],sizeof(double)**nef);
	    }else{
	      memcpy(&b[0],&bcp[0],sizeof(double)**nef);
	      FORTRAN(reorderrhs,(au,&b[0],&vel[4**nef],neighblock,
				  &nneighblock));
	    }
	    dgmresmain(nef,&b[0],&vel[4**nef],&nz_num,iam,jam,am,
		       &isym,&itol,&tol,&itmax,&iter,
		       &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
		       &ligw,rwork,iwork,nestart,&num_cpus,&dgmrestol);
		    
	    if(ierr>1){
	      printf("*WARNING in compfluid: error message from predgmres (p)=%d\n",ierr);
	      if(ierr>ierrmax)ierrmax=ierr;
	    }
	  }
      
	  /*	  spooles(&au[*nflnei],am,adb,aub,&sigma,&b[0],iam,jam,
		  nef,&nz_num,&symmetryflag,&inputformat,&nz_num);
		  memcpy(&vel[4**nef],&b[0],sizeof(double)**nef);*/
	  
	  /* non-orthogonal pressure correction */
		    
	  for(ip=0;ip<iitp;ip++){
			
	    /* calculate the p' gradient at the 
	       face centers */
			
	    iflag=1;
	    extrapol_dpelmain(nface,ielfa,xrlfa,vel,vfa,
			      ifabou,xbounact,nef,gradpcel,gradpcfa,neifa,rf,area,volume,
			      xle,xxi,&icyclic,xxn,ipnei,ifatie,
			      coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,
			      nactdoh,&iflag,xxj,xlet,c,&num_cpus,xxna,&ncfd,&ip);
	      
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
			
	    memcpy(&bcp[0],&b[0],sizeof(double)**nef);
	    FORTRAN(reorderrhs,(au,&b[0],&vel[4**nef],neighblock,&nneighblock));
			
	    /* initial guess: 0 */
			
	    DMEMSET(vel,4**nef,5**nef,0.);
			
	    /* calculation of the mass conservation residual */
			
	    FORTRAN(calcrespfluid,(nef,&b[0],&au[*nflnei],
				   &temp[4**nef],&relnormp));

	    for(k=1;k<num_cpus+1;k++){
	      if(k>1){
		memcpy(&b[0],&bcp[0],sizeof(double)**nef);
		FORTRAN(reorderrhs,(au,&b[0],&vel[4**nef],neighblock,
				    &nneighblock));
	      }
	      dgmresmain(nef,&b[0],&vel[4**nef],&nz_num,iam,jam,am,
			 &isym,&itol,&tol,&itmax,&iter,
			 &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
			 &ligw,rwork,iwork,nestart,&num_cpus,&dgmrestol);
			
	      if(ierr>1){
		printf("*WARNING in compfluid: error message from predgmres (p)=%d\n",ierr);
		if(ierr>ierrmax)ierrmax=ierr;
	      }
	    }
      
	    /*	    spooles(&au[*nflnei],am,adb,aub,&sigma,&b[0],iam,jam,
		    nef,&nz_num,&symmetryflag,&inputformat,&nz_num);
		    memcpy(&vel[4**nef],&b[0],sizeof(double)**nef);*/
			
	  }   // end loop iitp
		    
	  /* calculate the p' gradient at the 
	     element centers */

	  /* 16.04.2020: should be iflag=1????? */
		    
	  iflag=0;
	  extrapol_dpelmain(nface,ielfa,xrlfa,vel,vfa,
			    ifabou,xbounact,nef,gradpcel,gradpcfa,neifa,rf,area,volume,
			    xle,xxi,&icyclic,xxn,ipnei,ifatie,
			    coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,
			    nactdoh,&iflag,xxj,xlet,c,&num_cpus,xxna,&ncfd,&ip);
		    
	  /* correction of the velocity v* at the element centers due
	     to the pressure change in order to get v** */
		    
	  //	  if((compressible==0)||(isimplec==0)){
	  correctvelmain(&auv[*nflnei],nef,volume,gradpcel,vel,&num_cpus);
	  /*	  }else{
		  FORTRAN(correctvel_simplec,(&auv[*nflnei],nef,volume,gradpcel,vel,
		  ipnei,auv));
		  }*/
		    
	  /* extrapolating the velocity from the elements centers to the face
	     centers,thereby taking the boundary conditions into account */
		    
	  extrapol_velmain(nface,ielfa,xrlfa,vel,vfa,
			   ifabou,xbounact,ipnei,nef,&icyclic,c,ifatie,xxn,gradvel,
			   gradvfa,neifa,rf,area,volume,xle,xxi,xxj,xlet,
			   coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
			   &iitg,&num_cpus,&compressible,xxna,&ncfd,cofa);
		    
	  /* correcting the flux to get mf** */
		    
	  correctfluxcompmain(nef,ipnei,neifa,neiel,flux,vfa,advfa,area,
			      vel,alet,ielfa,ale,ifabou,ielmatf,mi,shcon,
			      ntmat_,&num_cpus,xxnj,gradpcfa,inlet);
		    
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
			   &iitg,c,&num_cpus,&compressible,xxna,&ncfd);
		    
	  /* calculating the lhs and rhs of the energy equation */
		    
	  DMEMSET(ad,0,*nef,0.);
	  DMEMSET(au,0,*nflnei+*nef,0.);
	  DMEMSET(b,0,*nef,0.);
		    
	  /* underrelaxation of the temperature only for the
	     compressible simple scheme */
		    
	  if(isimplec==0){
	    memcpy(&temp[0],&vel[0],sizeof(double)**nef);
	  }
		    
	  /* convective scheme */
		    
	  //	      if(ischeme==1){
	  hrt_udmain(nface,ielfa,vel,ipnei,nef,flux,vfa,&num_cpus,
		     xxi,xle,gradtel,neij);
	  /*	      }else{
		      FORTRAN(hrt_mod_smart,(nface,ielfa,vel,gradtel,gamma,
		      xlet,xxn,xxj,ipnei,&betam,nef,flux,vfa));
		      }*/
		    
	  //	printf("mafilltcompmain\n");
	  mafilltcompmain(nef,ipnei,neifa,neiel,vfa,xxn,area,
			  au,&au[*nflnei],jq,irow,&nzs,b,vel,umel,alet,ale,gradtfa,xxi,
			  body,volume,ielfa,lakonf,ifabou,nbody,nef,
			  &dtimef,velo,veloo,cvfa,hcfa,cvel,gradvel,xload,gamma,
			  xrlfa,xxj,nactdohinv,&a1,&a2,&a3,flux,iau6,xxni,xxnj,
			  iturbulent,of2,sc);
		    
	  for(i=0;i<*nef;i++){rwork[i]=1./au[*nflnei+i];}
	      
	  /* reordering the lhs (getting rid of zeros) */
		    
	  reorderlhsmain(au,am,iamorig,&nz_num,&num_cpus);
		    
	  /* modifying the rhs (taking common contributions
	     between the cpu-blocks into account) */
		    
	  memcpy(&bcp[0],&b[0],sizeof(double)**nef);
	  FORTRAN(reorderrhs,(au,&b[0],&vel[0],neighblock,&nneighblock));
		    
	  /* calculation of the energy residual */
		    
	  calcrestfluidmain(nef,am,&b[0],&au[*nflnei],iam,jam,
                            &vel[0],&relnormt,nestart,&num_cpus);

	  for(k=1;k<num_cpus+1;k++){
	    if(k>1){
	      memcpy(&b[0],&bcp[0],sizeof(double)**nef);
	      FORTRAN(reorderrhs,(au,&b[0],&vel[0],neighblock,&nneighblock));
	    }
	    dgmresmain(nef,&b[0],&vel[0],&nz_num,iam,jam,am,
		       &isym,&itol,&tol,&itmax,&iter,
		       &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
		       &ligw,rwork,iwork,nestart,&num_cpus,&dgmrestol);
		    
	    if(ierr>1){
	      printf("*WARNING in compfluid: error message from predgmres (T)=%d\n",ierr);
	      if(ierr>ierrmax)ierrmax=ierr;
	    }
	  }
      
	  /*	  spooles(&au[*nflnei],am,adb,aub,&sigma,&b[0],iam,jam,
		  nef,&nz_num,&symmetryflag,&inputformat,&nz_num);
		  memcpy(&vel[0],&b[0],sizeof(double)**nef);*/

	  /* underrelaxation of the temperature only for compressible
	     simple scheme */
		    
	  if(isimplec==0){
	    for(i=0;i<*nef;i++){vel[i]=0.8*vel[i]+0.2*temp[i];}
	  }
		    
	  /* extrapolation of the temperature at the element centers
	     to the face centers */
		    
	  extrapol_telmain(nface,ielfa,xrlfa,vel,vfa,
			   ifabou,xbounact,nef,gradtel,gradtfa,neifa,rf,area,volume,
			   xle,xxi,&icyclic,xxn,ipnei,ifatie,xload,xlet,xxj,
			   coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
			   &iitg,c,&num_cpus,&compressible,xxna,&ncfd);
		    
	  /* recalculating the density  */
		    
	  /* calculating the density at the element centers */
		    
	  calcrhoelcompmain(nef,vel,shcon,ielmatf,ntmat_,mi,
			    &num_cpus);
		    
	  /* calculating the density at the face centers */
		    
	  if(ischeme==1){
	    hrr_udmain(nface,vfa,shcon,ielmatf,ntmat_,
		       mi,ielfa,ipnei,vel,nef,flux,
		       &num_cpus,xxi,xle,gradpel,
		       gradtel,neij);
	  }else{
	    hrr_mod_smartmain(nface,vfa,shcon,ielmatf,ntmat_,
			      mi,ielfa,ipnei,vel,nef,flux,
			      gradpel,gradtel,xxj,xlet,
			      &num_cpus);
	  }
		    
	}
		
      }

      if((*ithermal>0)&&(compressible==0)){

	/* calculating the lhs and rhs of the energy equation */
		
	DMEMSET(ad,0,*nef,0.);
	DMEMSET(au,0,*nflnei+*nef,0.);
	DMEMSET(b,0,*nef,0.);
		
	/* calculate gamma (Ph.D. Thesis Jasak) */

	calcgammatmain(nface,ielfa,vel,gradtel,gamma,xlet,xxj,ipnei,
		       &betam,nef,flux,&num_cpus);

	mafilltmain(nef,ipnei,neifa,neiel,vfa,xxn,area,
		    au,&au[*nflnei],jq,irow,&nzs,b,vel,umel,alet,ale,gradtfa,xxi,
		    body,volume,ielfa,lakonf,ifabou,nbody,nef,
		    &dtimef,velo,veloo,cvfa,hcfa,cvel,gradvel,xload,gamma,
		    xrlfa,xxj,nactdohinv,&a1,&a2,&a3,flux,iau6,xxni,xxnj,
		    iturbulent,of2,sc);
		
	for(i=0;i<*nef;i++){rwork[i]=1./au[*nflnei+i];}
		
	/* reordering the lhs (getting rid of zeros) */
		
	reorderlhsmain(au,am,iamorig,&nz_num,&num_cpus);
		
	/* modifying the rhs (taking common contributions
	   between the cpu-blocks into account) */
		
	FORTRAN(reorderrhs,(au,&b[0],&vel[0],neighblock,&nneighblock));
		
	/* calculation of the energy residual */
		
	calcrestfluidmain(nef,am,&b[0],&au[*nflnei],iam,jam,
			  &vel[0],&relnormt,nestart,&num_cpus);

	dgmresmain(nef,&b[0],&vel[0],&nz_num,iam,jam,am,
		   &isym,&itol,&tol,&itmax,&iter,
		   &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
		   &ligw,rwork,iwork,nestart,&num_cpus,&dgmrestol);
		
	if(ierr>1){
	  printf("*WARNING in compfluid: error message from predgmres (T)=%d\n",ierr);
	  if(ierr>ierrmax)ierrmax=ierr;
	}
		
	/* extrapolation of the temperature at the element centers
	   to the face centers */
		
	extrapol_telmain(nface,ielfa,xrlfa,vel,vfa,
			 ifabou,xbounact,nef,gradtel,gradtfa,neifa,rf,area,volume,
			 xle,xxi,&icyclic,xxn,ipnei,ifatie,xload,xlet,xxj,
			 coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
			 &iitg,c,&num_cpus,&compressible,xxna,&ncfd);
		
	/* recalculating the density  */
		
	/* calculating the density at the element centers */
		
	FORTRAN(calcrhoel,(nef,vel,rhcon,nrhcon,ielmatf,ntmat_,
			   ithermal,mi));
		
	/* calculating the density at the face centers */
		
	FORTRAN(calcrhofa,(nface,vfa,rhcon,nrhcon,ielmatf,ntmat_,
			   ithermal,mi,ielfa));

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
	              au,&au[*nflnei],jq,irow,&nzs,b,vel,umfa,alet,ale,gradkfa,xxi,
	              body,volume,ielfa,lakonf,ifabou,nbody,nef,
	              &dtimef,velo,veloo,cvfa,hcfa,cvel,gradvel,xload,gamma,
	              xrlfa,xxj,nactdohinv,&a1,&a2,&a3,flux,iau6,xxni,xxnj,
	              iturbulent,f1,of2,yy,umel,gradkel,gradoel,sc);
	}else{

	  /* underrelaxation of the temperature only for the
	     compressible simple scheme */
		    
	  if(isimplec==0){
	    memcpy(&temp[6**nef],&vel[6**nef],sizeof(double)**nef);
	  }
		    
	  /* convective scheme */
		    
	  //	      if(ischeme==1){
	  hrk_udmain(nface,ielfa,vel,ipnei,nef,flux,vfa,&num_cpus);
	  /*	      }else{
		      FORTRAN(hrk_mod_smart,(nface,ielfa,vel,gradtel,gamma,
		      xlet,xxn,xxj,ipnei,&betam,nef,flux,vfa));
		      }*/
		    
	  mafillkcompmain(nef,ipnei,neifa,neiel,vfa,xxn,area,
			  au,&au[*nflnei],jq,irow,&nzs,b,vel,umfa,alet,ale,gradkfa,xxi,
			  body,volume,ielfa,lakonf,ifabou,nbody,nef,
			  &dtimef,velo,veloo,cvfa,hcfa,cvel,gradvel,xload,
			  xrlfa,xxj,nactdohinv,&a1,&a2,&a3,flux,iau6,xxni,xxnj,
			  iturbulent,f1,of2,yy,umel,gradkel,gradoel,inlet,sc);
	}
		
	for(i=0;i<*nef;i++){rwork[i]=1./au[*nflnei+i];}
		
	/* reordering the lhs (getting rid of zeros) */
		
	reorderlhsmain(au,am,iamorig,&nz_num,&num_cpus);
		
	/* modifying the rhs (taking common contributions
	   between the cpu-blocks into account) */
		
	FORTRAN(reorderrhs,(au,&b[0],&vel[6**nef],neighblock,&nneighblock));
		
	dgmresmain(nef,&b[0],&vel[6**nef],&nz_num,iam,jam,am,
		   &isym,&itol,&tol,&itmax,&iter,
		   &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
		   &ligw,rwork,iwork,nestart,&num_cpus,&dgmrestol);
	if(ierr>1){
	  printf("*WARNING in compfluid: error message from predgmres (k)=%d\n",ierr);
	  if(ierr>ierrmax)ierrmax=ierr;
	}
		
	/* underrelaxation of the temperature only for compressible
	   simple scheme */
		
	if((compressible==1)&&(isimplec==0)){
	  for(i=6**nef;i<7**nef;i++){vel[i]=0.8*vel[i]+0.2*temp[i];}
	}
		
	/* calculating the lhs and rhs of the omega-equation */
		
	DMEMSET(au,0,*nflnei+*nef,0.);
	DMEMSET(b,0,*nef,0.);
		
	if(compressible==0){
		    
	  /* calculate gamma (Ph.D. Thesis Jasak) */
		    
	  FORTRAN(calcgammao,(nface,ielfa,vel,gradoel,gamma,xlet,xxn,xxj,
			      ipnei,&betam,nef,flux));
		    
	  mafillomain(nef,ipnei,neifa,neiel,vfa,xxn,area,
		      au,&au[*nflnei],jq,irow,&nzs,b,vel,umfa,alet,ale,gradofa,xxi,
		      body,volume,ielfa,lakonf,ifabou,nbody,nef,
		      &dtimef,velo,veloo,cvfa,hcfa,cvel,gradvel,xload,gamma,
		      xrlfa,xxj,nactdohinv,&a1,&a2,&a3,flux,iau6,xxni,xxnj,
		      iturbulent,f1,of2,gradkel,gradoel,sc);
	}else{
		    
	  /* underrelaxation of the temperature only for the
	     compressible simple scheme */
		    
	  if(isimplec==0){
	    memcpy(&temp[7**nef],&vel[7**nef],sizeof(double)**nef);
	  }
		    
	  /* convective scheme */
		    
	  //	      if(ischeme==1){
	  hro_udmain(nface,ielfa,vel,ipnei,nef,flux,vfa,&num_cpus);
	  /*	      }else{
		      FORTRAN(hro_mod_smart,(nface,ielfa,vel,gradtel,gamma,
		      xlet,xxn,xxj,ipnei,&betam,nef,flux,vfa));
		      }*/
		    
	  mafillocompmain(nef,ipnei,neifa,neiel,vfa,xxn,area,
			  au,&au[*nflnei],jq,irow,&nzs,b,vel,umfa,alet,ale,gradofa,xxi,
			  body,volume,ielfa,lakonf,ifabou,nbody,nef,
			  &dtimef,velo,veloo,cvfa,hcfa,cvel,gradvel,xload,
			  xrlfa,xxj,nactdohinv,&a1,&a2,&a3,flux,iau6,xxni,xxnj,
			  iturbulent,f1,of2,gradkel,gradoel,inlet,sc);
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
		   &ligw,rwork,iwork,nestart,&num_cpus,&dgmrestol);
	if(ierr>1){
	  printf("*WARNING in compfluid: error message from predgmres (om)=%d\n",ierr);
	  if(ierr>ierrmax)ierrmax=ierr;
	}
		
	/* underrelaxation of the temperature only for compressible
	   simple scheme */
		
	if((compressible==1)&&(isimplec==0)){
	  for(i=7**nef;i<8**nef;i++){vel[i]=0.8*vel[i]+0.2*temp[i];}
	}
		
	/* extrapolation of the turbulence variables at the element centers
	   to the face centers */
	/*	      sum=0;
		      for(i=6**nef;i<7**nef;i++){
		      printf("i=%d,vel=%e\n",i+1,vel[i]);
		      if(fabs(vel[i])>sum) sum=fabs(vel[i]);
		      }
		      printf("sum=%e\n\n",sum);*/
		
	extrapol_kelmain(nface,ielfa,xrlfa,vel,vfa,
			 ifabou,xbounact,nef,gradkel,gradkfa,neifa,rf,area,volume,
			 xle,xxi,&icyclic,xxn,ipnei,ifatie,xlet,xxj,
			 coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
			 umfa,physcon,&iitg,c,&num_cpus,&compressible,xxna,&ncfd,
			 inlet);
	  
	extrapol_oelmain(nface,ielfa,xrlfa,vel,vfa,
			 ifabou,xbounact,nef,gradoel,gradofa,neifa,rf,area,volume,
			 xle,xxi,&icyclic,xxn,ipnei,ifatie,xlet,xxj,
			 coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
			 umfa,physcon,&iitg,c,dy,&num_cpus,&compressible,xxna,
			 &ncfd,inlet);
		
      }

      /* extrapolating the velocity from the elements centers to the face
	 centers,thereby taking the boundary conditions into account */
	    
      /*   extrapol_velmain(nface,ielfa,xrlfa,vel,vfa,
	   ifabou,xbounact,ipnei,nef,&icyclic,c,ifatie,xxn,gradvel,
	   gradvfa,neifa,rf,area,volume,xle,xxi,xxj,xlet,
	   coefmpc,nmpc,labmpc,ipompc,nodempc,ifaext,nfaext,nactdoh,
	   &iitg,&num_cpus,&compressible,xxna,&ncfd);*/
	    
      //      FORTRAN(writevfa,(vfa,nface,nactdohinv,ielfa));

      /* end subiterations */

      relnormmax=relnormt;
      if(relnormv>relnormmax){relnormmax=relnormv;}
      if(relnormp>relnormmax){relnormmax=relnormp;}
	    
      //      fprintf(f3,"%d %d %11.4e %11.4e %11.4e %11.4e\n",iincf,iit,dtimef,relnormt,relnormv,relnormp);
	    
      if(*nmethod==1){
		
	/* steady state flow:
	   calculate the velocity only once in each increment */
		
	if(relnormmax<1.e-10) iconvergence=1;
	if(compressible==0){
		    
	  /* check whether dtimef was changed in the last
	     increment */
		    
	  if(beta>0.){
	    a1=1.5/dtimef;
	    a2=-2./dtimef;
	    a3=0.5/dtimef;
	    beta=0.;
	  }
	  fprintf(f3,"%d %d %11.4e %11.4e %11.4e %11.4e\n",iincf,iit,dtimef,relnormt,relnormv,relnormp);
	  break;
	}else if((relnormmax<1.e-3)||(iit==iitf)){
	  fprintf(f3,"%d %d %11.4e %11.4e %11.4e %11.4e\n",iincf,iit,dtimef,relnormt,relnormv,relnormp);
	  break;
	}
      }else{
		
	/* transient flow:
	   calculate the velocity repeatedly in each increment */
		
	//	  if((relnormmax<1.e-8)||(iit==iitf)){
	if((relnormmax<1.e-5)||(iit==iitf)){
		    
	  /* dynamic change of time increment for transient
	     compressible flow */
		    
	  if(compressible!=0){
			
	    /* compressible flow */
			
	    if((iito<3)&&(iit<3)){
			    
	      /* fast convergence */
			    
	      dtimef*=1.05;
	      printf("increased time increment to %e\n",dtimef);
	      a1=1./dtimef;
	      a2=-a1;
	    }else if(iit>iitf-1){
			    
	      /* divergence */
			    
	      timef-=dtimef;
	      memcpy(&vel[0],&velo[0],sizeof(double)*8**nef);
	      memcpy(&velo[0],&veloo[0],sizeof(double)*8**nef);
	      dtimef*=0.25;
	      printf("divergence: recalculated increment with reduced time increment  %e\n",dtimef);
	      a1=1./dtimef;
	      a2=-a1;
	    }else if((iito>10)&&(iit>10)){
			    
	      /* slow convergence */
			    
	      dtimef*=0.95;
	      printf("decreased time increment to %e\n",dtimef);
	      a1=1./dtimef;
	      a2=-a1;
	    }
	    iito=iit;
	  }
	  fprintf(f3,"%d %d %11.4e %11.4e %11.4e %11.4e\n",iincf,iit,dtimef,relnormt,relnormv,relnormp);
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

      FORTRAN(extrapolatefluid,(nk,ipofano,ifano,inum,vfa,vold,ielfa,
				ithermal,&imach,&ikappa,xmach,xkappa,shcon,
				nshcon,ntmat_,ielmatf,physcon,mi,&iturb,xturb,
				gradtfa,gradvfa,gradpfa,gradkfa,gradofa,co,
				cofa,ifabou));
	    
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
  SFREE(xxna);SFREE(rf);SFREE(ale);SFREE(alet);SFREE(inlet);SFREE(h);
  SFREE(sc);
  if(*iturbulent>0){
    SFREE(dy);
    if(*iturbulent>2) SFREE(yy);
  }
    
  SFREE(ifabou);SFREE(umfa);SFREE(umel);
    
  SFREE(gradvel);SFREE(gradvfa);SFREE(au);SFREE(ad);SFREE(b);SFREE(advfa);
  SFREE(ap);SFREE(bp);SFREE(gradpel);SFREE(rwork);SFREE(gradpfa);
  SFREE(gradpcel);SFREE(gradpcfa);SFREE(bvcp);SFREE(bcp);
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
    
  SFREE(inum);
    
  SFREE(ipofano);SFREE(ifano);
    
  if(*nbody>0) SFREE(body);
    
  *ithermal=ithermalref;
    
  SFREE(temp);SFREE(gamma);
    
  SFREE(iam);SFREE(jam);SFREE(iamorig);SFREE(am);
  SFREE(nestart);SFREE(ineighblock);SFREE(neighblock);

  SFREE(iponofa);SFREE(inofa);
  
  *cop=co;*ipkonfp=ipkonf;*lakonfp=lakonf;*voldp=vold;
    
  return;
    
} 
