/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2019 Guido Dhondt                          */

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

#ifdef __WIN32
  _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

#ifdef CALCULIX_MPI
#include <spoolesMPI.h>
#endif

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"

#ifdef CALCULIX_MPI
ITG myid = 0, nproc = 0;
#endif

int main(int argc,char *argv[])
{
  
FILE *f1;
    
char *sideload=NULL, *set=NULL, *matname=NULL, *orname=NULL, *amname=NULL,
     *filab=NULL, *lakon=NULL, *labmpc=NULL, *prlab=NULL, *prset=NULL, 
     jobnamec[792]="",jobnamef[132]="",output[5]="asc ", *typeboun=NULL,
     *inpc=NULL,*tieset=NULL,*cbody=NULL,fneig[132]="",*sideloadtemp=NULL,
     kind1[2]="T",kind2[2]="T",*heading=NULL,*objectset=NULL;
  
ITG *kon=NULL, *nodeboun=NULL, *ndirboun=NULL, *ipompc=NULL,
    *nodempc=NULL, *nodeforc=NULL, *ndirforc=NULL,
    *nelemload=NULL,im,*inodesd=NULL,nload1,*idefforc=NULL,
    *nactdof=NULL, *icol=NULL,*ics=NULL,itempuser[3]={0,0,-2},
    *jq=NULL, *mast1=NULL, *irow=NULL, *rig=NULL,*idefbody=NULL,
    *ikmpc=NULL, *ilmpc=NULL, *ikboun=NULL, *ilboun=NULL,
    *nreorder=NULL,*ipointer=NULL,*idefload=NULL,
    *istartset=NULL, *iendset=NULL, *ialset=NULL, *ielmat=NULL,
    *ielorien=NULL, *nrhcon=NULL, *nodebounold=NULL, *ndirbounold=NULL,
    *nelcon=NULL, *nalcon=NULL, *iamforc=NULL,  *iamload=NULL,
    *iamt1=NULL, *namta=NULL, *ipkon=NULL, *iamboun=NULL,
    *nplicon=NULL, *nplkcon=NULL, *inotr=NULL, *iponor=NULL, *knor=NULL,
    *ikforc=NULL, *ilforc=NULL, *iponoel=NULL, *inoel=NULL, *nshcon=NULL,
    *ncocon=NULL,*ibody=NULL,*ielprop=NULL,*islavsurf=NULL,
    *ipoinpc=NULL,mt,nxstate,nload0,iload,*iuel=NULL,*ne2boun=NULL;
     
ITG nk,ne,nboun,nmpc,nforc,nload,nprint=0,nset,nalset,nentries=17,
  nmethod,neq[3]={0,0,0},i,mpcfree=1,mei[4],j,nzl,nam,nbounold=0,
  nforcold=0,nloadold=0,nbody,nbody_=0,nbodyold=0,network=0,nheading_=0,
  k,nzs[3],nmpc_=0,nload_=0,nforc_=0,istep,istat,nboun_=0,nintpoint=0,
  iperturb[2]={0,0},nmat,ntmat_=0,norien,ithermal[2]={0,0},nmpcold,
  iprestr,kode,isolver=0,nslavs=0,nkon_=0,ne0,nkon0,mortar=0,
  jout[2]={1,1},nlabel,nkon=0,idrct,jmax[2],iexpl,nevtot=0,ifacecount=0,
  iplas=0,npmat_=0,mi[3]={0,3,1},ntrans,mpcend=-1,namtot_=0,iumat=0,
  icascade=0,maxlenmpc,mpcinfo[4],ne1d=0,ne2d=0,infree[4]={0,0,0,0},
  callfrommain,nflow=0,jin=0,irstrt[2]={0,0},nener=0,jrstrt=0,nenerold,
  nline,*ipoinp=NULL,*inp=NULL,ntie,ntie_=0,mcs=0,nprop_=0,
  nprop=0,itpamp=0,iviewfile,nkold,nevdamp_=0,npt_=0,cyclicsymmetry,
  nmethodl,iaxial=1,inext=0,icontact=0,nobject=0,nobject_=0,iit=-1,
  nzsprevstep[3],memmpcref_,mpcfreeref=-1,maxlenmpcref,*nodempcref=NULL,
  *ikmpcref=NULL,isens=0,namtot=0,nstam=0,ndamp=0,nef=0;

ITG *meminset=NULL,*rmeminset=NULL;

ITG nzs_,nk_=0,ne_=0,nset_=0,nalset_=0,nmat_=0,norien_=0,nam_=0,
    ntrans_=0,ncs_=0,nstate_=0,ncmat_=0,memmpc_=0,nprint_=0,nuel_=0;
    
double *co=NULL, *xboun=NULL, *coefmpc=NULL, *xforc=NULL,*clearini=NULL,
	*xload=NULL, *xbounold=NULL, *xforcold=NULL,
	*vold=NULL, *sti=NULL, *xloadold=NULL, *xnor=NULL,
	*reorder=NULL,*dcs=NULL, *thickn=NULL, *thicke=NULL, *offset=NULL,
	*elcon=NULL, *rhcon=NULL, *alcon=NULL, *alzero=NULL, *t0=NULL, *t1=NULL,
	*prestr=NULL, *orab=NULL, *amta=NULL, *veold=NULL, *accold=NULL,
        *t1old=NULL, *eme=NULL, *plicon=NULL, *pslavsurf=NULL, *plkcon=NULL,
	*xstate=NULL, *trab=NULL, *ener=NULL, *shcon=NULL, *cocon=NULL,
        *cs=NULL,*tietol=NULL,*fmpc=NULL,*prop=NULL,*t0g=NULL,*t1g=NULL,
        *xbody=NULL,*xbodyold=NULL,*coefmpcref=NULL,*dacon=NULL,*vel=NULL,
        *velo=NULL,*veloo=NULL,energy[5]={0.,0.,0.,0.,0.};
    
 double ctrl[57]={4.5,8.5,9.5,16.5,10.5,4.5,0.,5.5,0.,0.,0.25,0.5,0.75,0.85,0.,0.,1.5,0.,0.005,0.01,0.,0.,0.02,1.e-5,1.e-3,1.e-8,1.e30,1.5,0.25,1.01,1.,1.,5.e-7,5.e-7,1.e-4,5.e-7,5.e-7,5.e-7,5.e-7,-1.,1.e20,1.e20,1.e20,1.e20,1.e20,1.e20,1.e20,1.5,0.5,20.5,0.5,1.5,1.5,0.001,0.1,100.5,60.5};

double fei[3],*xmodal=NULL,timepar[5],
  alpha[2]={0.,0.5},ttime=0.,qaold[2]={0.,0.},physcon[14]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

#ifdef CALCULIX_MPI
MPI_Init(&argc, &argv) ;
MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;
MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
#endif

if(argc==1){printf("Usage: CalculiX.exe -i jobname\n");FORTRAN(stop,());}
else{
  for(i=1;i<argc;i++){
    if(strcmp1(argv[i],"-i")==0) {
    strcpy(jobnamec,argv[i+1]);strcpy1(jobnamef,argv[i+1],132);jin++;break;}
    if(strcmp1(argv[i],"-v")==0) {
	printf("\nThis is Version 2.16\n\n");
	FORTRAN(stop,());
    }
  }
  if(jin==0){strcpy(jobnamec,argv[1]);strcpy1(jobnamef,argv[1],132);}

  for(i=1;i<argc;i++){
    if(strcmp1(argv[i],"-o")==0) {
    strcpy(output,argv[i+1]);break;}
  }
}
//setenv("CCX_JOBNAME_GETJOBNAME",jobnamec,1);
putenv("CCX_JOBNAME_GETJOBNAME=jobnamec");

#ifdef BAM
ITG lop=0,lrestart=0,kstep=1,kinc=1;
double time[2],dtime;
FORTRAN(uexternaldb,(&lop,&lrestart,time,&dtime,&kstep,&kinc));
#endif

FORTRAN(openfile,(jobnamef,output));

printf("\n************************************************************\n\n");
printf("CalculiX Version 2.16, Copyright(C) 1998-2019 Guido Dhondt\n");
printf("CalculiX comes with ABSOLUTELY NO WARRANTY. This is free\n");
printf("software, and you are welcome to redistribute it under\n");
printf("certain conditions, see gpl.htm\n\n");
printf("************************************************************\n\n");
printf("You are using an executable made on Mo 25. Nov 18:56:47 CET 2019\n");
fflush(stdout);

istep=0;
istat=0;
iprestr=0;
kode=0;

/* default solver */

#if defined(SGI)
 isolver=4;
#elif defined(PARDISO)
 isolver=7;
#elif defined(SPOOLES)
 isolver=0;
#elif defined(TAUCS)
 isolver=5;
#else
 isolver=3;
#endif

 NNEW(ipoinp,ITG,2*nentries);

/* conservative estimate of the fields to be allocated */

 readinput(jobnamec,&inpc,&nline,&nset_,ipoinp,&inp,&ipoinpc,ithermal,&nuel_); 

NNEW(set,char,81*nset_);
NNEW(meminset,ITG,nset_);
NNEW(rmeminset,ITG,nset_);
NNEW(iuel,ITG,4*nuel_);

FORTRAN(allocation,(&nload_,&nforc_,&nboun_,&nk_,&ne_,&nmpc_,&nset_,&nalset_,
   &nmat_,&ntmat_,&npmat_,&norien_,&nam_,&nprint_,mi,&ntrans_,
   set,meminset,rmeminset,&ncs_,&namtot_,&ncmat_,&memmpc_,&ne1d,
   &ne2d,&nflow,jobnamec,irstrt,ithermal,&nener,&nstate_,&istep,
   inpc,ipoinp,inp,&ntie_,&nbody_,&nprop_,ipoinpc,&nevdamp_,&npt_,&nslavs,
   &nkon_,&mcs,&mortar,&ifacecount,&nintpoint,infree,&nheading_,&nobject_,
   iuel,&iprestr,&nstam,&ndamp,&nef,&nbounold,&nforcold,&nloadold,&nbodyold));

SFREE(set);SFREE(meminset);SFREE(rmeminset);mt=mi[1]+1;
NNEW(heading,char,66*nheading_);

nzs_=20000000;

nload=0;nbody=0;nforc=0;nboun=0;nk=0;nmpc=0;nam=0;

/* caveat: if changing next line:
   - change noelfiles appropriately
   - change nlabel in geomview.f, expand.c, storeresidual.f
     and createmddof.f
   - change the dimension of label in geomview.f
   - change the documentation (tex-file)  */

nlabel=48;

while(istat>=0) {

  fflush(stdout);

  /* in order to reduce the number of variables to be transferred to
     the subroutines, the max. field sizes are (for most fields) copied
     into the real sizes */

  nzs[1]=nzs_;
//  nprint=nprint_;

  if((istep==0)||(irstrt[0]<0)) {
    ne=ne_;
    nset=nset_;
    nalset=nalset_;
    nmat=nmat_;
    norien=norien_;
    ntrans=ntrans_;
    ntie=ntie_;
//    nobject=nobject_;

    /* allocating space before the first step */

    /* coordinates and topology */

    NNEW(co,double,3*nk_);
    NNEW(kon,ITG,nkon_);
    NNEW(ipkon,ITG,ne_);
    NNEW(lakon,char,8*ne_);

    /* property cards */

    if(nprop_>0){
	NNEW(ielprop,ITG,ne_);
	for(i=0;i<ne_;i++) ielprop[i]=-1;
	NNEW(prop,double,nprop_);
    }

    /* fields for 1-D and 2-D elements */

    if((ne1d!=0)||(ne2d!=0)){
	NNEW(iponor,ITG,2*nkon_);
	for(i=0;i<2*nkon_;i++) iponor[i]=-1;
	NNEW(xnor,double,36*ne1d+24*ne2d);
	NNEW(knor,ITG,24*(ne1d+ne2d)*(mi[2]+1));
	NNEW(thickn,double,2*nk_);
	NNEW(thicke,double,mi[2]*nkon_);
	NNEW(offset,double,2*ne_);
	NNEW(iponoel,ITG,nk_);
	NNEW(inoel,ITG,9*ne1d+24*ne2d);
	NNEW(rig,ITG,nk_);
	NNEW(ne2boun,ITG,2*nk_);
	if(infree[2]==0)infree[2]=1;
    }

    /* SPC's */

    NNEW(nodeboun,ITG,nboun_);
    NNEW(ndirboun,ITG,nboun_);
    NNEW(typeboun,char,nboun_+1);
    if((istep == 0)||((irstrt[0]<0)&&(nam_>0)))NNEW(iamboun,ITG,nboun_);
    NNEW(xboun,double,nboun_);
    NNEW(ikboun,ITG,nboun_);
    NNEW(ilboun,ITG,nboun_);

    /* MPC's */

    NNEW(ipompc,ITG,nmpc_);
    NNEW(nodempc,ITG,3*memmpc_);
    for(i=0;i<3*memmpc_;i+=3){nodempc[i+2]=i/3+2;}
    nodempc[3*memmpc_-1]=0;
    NNEW(coefmpc,double,memmpc_);
    NNEW(labmpc,char,20*nmpc_+1);
    NNEW(ikmpc,ITG,nmpc_);
    NNEW(ilmpc,ITG,nmpc_);
    NNEW(fmpc,double,nmpc_);

    /* nodal loads */

    NNEW(nodeforc,ITG,2*nforc_);
    NNEW(ndirforc,ITG,nforc_);
    if((istep == 0)||((irstrt[0]<0)&&(nam_>0)))NNEW(iamforc,ITG,nforc_);
    NNEW(idefforc,ITG,nforc_);
    NNEW(xforc,double,nforc_);
    NNEW(ikforc,ITG,nforc_);
    NNEW(ilforc,ITG,nforc_);

    /* distributed facial loads */

    NNEW(nelemload,ITG,2*nload_);
    if((istep == 0)||((irstrt[0]<0)&&(nam_>0)))NNEW(iamload,ITG,2*nload_);
    NNEW(idefload,ITG,nload_);
    NNEW(sideload,char,20*nload_);
    NNEW(xload,double,2*nload_);

    /* distributed volumetric loads */

    NNEW(cbody,char,81*nbody_);
    NNEW(idefbody,ITG,nbody_);
    NNEW(ibody,ITG,3*nbody_);
    NNEW(xbody,double,7*nbody_);
    NNEW(xbodyold,double,7*nbody_);

    /* printing output */

    NNEW(prlab,char,6*nprint_);
    NNEW(prset,char,81*nprint_);

    /* set definitions */

    NNEW(set,char,81*nset);
    NNEW(istartset,ITG,nset);
    NNEW(iendset,ITG,nset);
    NNEW(ialset,ITG,nalset);

    /* (hyper)elastic constants */

    NNEW(elcon,double,(ncmat_+1)*ntmat_*nmat);
    NNEW(nelcon,ITG,2*nmat);

    /* density */

    NNEW(rhcon,double,2*ntmat_*nmat);
    NNEW(nrhcon,ITG,nmat);

    /* damping */

    if(ndamp>0){NNEW(dacon,double,nmat);}

    /* specific heat */

    NNEW(shcon,double,4*ntmat_*nmat);
    NNEW(nshcon,ITG,nmat);

    /* thermal expansion coefficients */

    NNEW(alcon,double,7*ntmat_*nmat);
    NNEW(nalcon,ITG,2*nmat);
    NNEW(alzero,double,nmat);

    /* conductivity */

    NNEW(cocon,double,7*ntmat_*nmat);
    NNEW(ncocon,ITG,2*nmat);

    /* isotropic and kinematic hardening coefficients*/

    if(npmat_>0){
	NNEW(plicon,double,(2*npmat_+1)*ntmat_*nmat);
	NNEW(nplicon,ITG,(ntmat_+1)*nmat);
	NNEW(plkcon,double,(2*npmat_+1)*ntmat_*nmat);
	NNEW(nplkcon,ITG,(ntmat_+1)*nmat);
    }

    /* linear dynamic properties */
    
    NNEW(xmodal,double,11+nevdamp_);
    xmodal[10]=nevdamp_+0.5;

    /* internal state variables (nslavs is needed for restart
       calculations) */

    if(mortar!=1){
	NNEW(xstate,double,nstate_*mi[0]*(ne+nslavs));
	nxstate=nstate_*mi[0]*(ne+nslavs);
    }else if(mortar==1){
	NNEW(xstate,double,nstate_*mi[0]*(ne+nintpoint));
	nxstate=nstate_*mi[0]*(ne+nintpoint);
    }

    /* material orientation */

    if((istep == 0)||((irstrt[0]<0)&&(norien>0))) {
	NNEW(orname,char,80*norien);
	NNEW(orab,double,7*norien);
	NNEW(ielorien,ITG,mi[2]*ne_);
    }

    /* transformations */

    if((istep == 0)||((irstrt[0]<0)&&(ntrans>0))) {
	NNEW(trab,double,7*ntrans);
	NNEW(inotr,ITG,2*nk_);
    }

    /* amplitude definitions */

    if((istep == 0)||((irstrt[0]<0)&&(nam_>0))) {
	NNEW(amname,char,80*nam_);
	NNEW(amta,double,2*namtot_);
	NNEW(namta,ITG,3*nam_);
    }

    if((istep == 0)||((irstrt[0]<0)&&(ithermal[0]>0))) {
	NNEW(t0,double,nk_);
	NNEW(t1,double,nk_);
	if((ne1d!=0)||(ne2d!=0)){
	    NNEW(t0g,double,2*nk_);
	    NNEW(t1g,double,2*nk_);
	}
    }

    /* the number in next line is NOT 1.2357111317 -> points
       to user input; instead it is a generic nonzero
       initialization */

    if(istep==0){
	DMEMSET(t0,0,nk_,1.2357111319);
	DMEMSET(t1,0,nk_,1.2357111319);
    }
    
    if((istep == 0)||((irstrt[0]<0)&&(ithermal[0]>0)&&(nam_>0)))NNEW(iamt1,ITG,nk_);

    if((istep==0)||((irstrt[0]<0)&&(iprestr>0)))NNEW(prestr,double,6*mi[0]*ne_);

    NNEW(vold,double,mt*nk_);
    NNEW(veold,double,mt*nk_);

    /* CFD-results */

    NNEW(vel,double,8*nef);
    NNEW(velo,double,8*nef);
    NNEW(veloo,double,8*nef);

    NNEW(ielmat,ITG,mi[2]*ne_);

    NNEW(matname,char,80*nmat);

    NNEW(filab,char,87*nlabel);

    /* tied constraints */

    if(ntie_>0){
      NNEW(tieset,char,243*ntie_);
      NNEW(tietol,double,3*ntie_);
      NNEW(cs,double,17*ntie_);
    }

    /* objectives for sensitivity analysis */

    if(nobject_>0){
      NNEW(objectset,char,324*nobject_);
      for(i=0;i<324*nobject_;i++){objectset[i]=' ';}
    }
    
    /* temporary fields for cyclic symmetry calculations */

    if((ncs_>0)||(npt_>0)){
      if(2*npt_>24*ncs_){
	NNEW(ics,ITG,2*npt_);
      }else{
	NNEW(ics,ITG,24*ncs_);
      }
      if(npt_>30*ncs_){
	  NNEW(dcs,double,npt_);
      }else{
	  NNEW(dcs,double,30*ncs_);
      }
    }

    /* slave faces */

    NNEW(islavsurf,ITG,2*ifacecount+2);

  }
  else {

      /* allocating and reallocating space for subsequent steps */
      
      if((nmethod!=4)&&(nmethod!=5)&&(nmethod!=8)&&(nmethod!=9)&& 
       ((abs(nmethod)!=1)||(iperturb[0]<2))){
        NNEW(veold,double,mt*nk_);
    }
    else{
      RENEW(veold,double,mt*nk_);
      DMEMSET(veold,mt*nk,mt*nk_,0.);
    }
    RENEW(vold,double,mt*nk_);
    DMEMSET(vold,mt*nk,mt*nk_,0.);

    RENEW(nodeboun,ITG,nboun_);
    RENEW(ndirboun,ITG,nboun_);
    RENEW(typeboun,char,nboun_+1);
    RENEW(xboun,double,nboun_);
    RENEW(ikboun,ITG,nboun_);
    RENEW(ilboun,ITG,nboun_);

    RENEW(nodeforc,ITG,2*nforc_);
    RENEW(ndirforc,ITG,nforc_);
    NNEW(idefforc,ITG,nforc_);
    RENEW(xforc,double,nforc_);
    RENEW(ikforc,ITG,nforc_);
    RENEW(ilforc,ITG,nforc_);

    RENEW(nelemload,ITG,2*nload_);
    NNEW(idefload,ITG,nload_);
    RENEW(sideload,char,20*nload_);
    RENEW(xload,double,2*nload_);

    RENEW(cbody,char,81*nbody_);
    NNEW(idefbody,ITG,nbody_);
    RENEW(ibody,ITG,3*nbody_);
    RENEW(xbody,double,7*nbody_);
    RENEW(xbodyold,double,7*nbody_);
    for(i=7*nbodyold;i<7*nbody_;i++) xbodyold[i]=0;

    if(nam > 0) {
      RENEW(iamforc,ITG,nforc_);
      RENEW(iamload,ITG,2*nload_);
      RENEW(iamboun,ITG,nboun_);
      RENEW(amname,char,80*nam_);
      RENEW(amta,double,2*namtot_);
      RENEW(namta,ITG,3*nam_);
    }

    RENEW(ipompc,ITG,nmpc_);

    RENEW(labmpc,char,20*nmpc_+1);
    RENEW(ikmpc,ITG,nmpc_);
    RENEW(ilmpc,ITG,nmpc_);
    RENEW(fmpc,double,nmpc_);

    if(ntrans > 0){
      RENEW(inotr,ITG,2*nk_);DMEMSET(inotr,2*nk,2*nk_,0);
    }

    RENEW(co,double,3*nk_);DMEMSET(co,3*nk,3*nk_,0.);

    if(ithermal[0] != 0){
	RENEW(t0,double,nk_);DMEMSET(t0,nk,nk_,0.);
	RENEW(t1,double,nk_);DMEMSET(t1,nk,nk_,0.);
	if((ne1d!=0)||(ne2d!=0)){
	    RENEW(t0g,double,2*nk_);DMEMSET(t0g,2*nk,2*nk_,0.);
	    RENEW(t1g,double,2*nk_);DMEMSET(t1g,2*nk,2*nk_,0.);
	}
	if(nam > 0) {RENEW(iamt1,ITG,nk_);}
    }

  }

  /* allocation of fields in the restart file */

  if(irstrt[0]<0){
    NNEW(nodebounold,ITG,nboun_);
    NNEW(ndirbounold,ITG,nboun_);
    NNEW(xbounold,double,nboun_);
    NNEW(xforcold,double,nforc_);
    NNEW(xloadold,double,2*nload_);
    if(ithermal[0]!=0) NNEW(t1old,double,nk_); 
    NNEW(sti,double,6*mi[0]*ne);
    NNEW(eme,double,6*mi[0]*ne);
    if(nener==1)NNEW(ener,double,mi[0]*ne*2);
    if(mcs>ntie_) RENEW(cs,double,17*mcs);
    if(mortar==1){
	NNEW(pslavsurf,double,3*nintpoint);
	NNEW(clearini,double,3*9*ifacecount);
    }
			
  }

  nenerold=nener;
  nkold=nk;

  /* opening the eigenvalue file and checking for cyclic symmetry */

  strcpy(fneig,jobnamec);
  strcat(fneig,".eig");
  cyclicsymmetry=0;
  if((f1=fopen(fneig,"rb"))!=NULL){
      if(fread(&cyclicsymmetry,sizeof(ITG),1,f1)!=1){
	  printf("*ERROR reading the information whether cyclic symmetry is involved in the eigenvalue file");
	  exit(0);
      }
      fclose(f1);
  }

  nmpcold=nmpc;

  /* reading the input file */

  if(istep==0)mortar=-1;
  FORTRAN(calinput,(co,&nk,kon,ipkon,lakon,&nkon,&ne,
            nodeboun,ndirboun,xboun,&nboun,
	    ipompc,nodempc,coefmpc,&nmpc,&nmpc_,nodeforc,ndirforc,xforc,&nforc,
	    &nforc_,nelemload,sideload,xload,&nload,&nload_,
	    &nprint,prlab,prset,&mpcfree,&nboun_,mei,set,istartset,iendset,
	    ialset,&nset,&nalset,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
	    alzero,t0,t1,matname,ielmat,orname,orab,ielorien,amname,
            amta,namta,&nam,&nmethod,iamforc,iamload,iamt1,
	    ithermal,iperturb,&istat,&istep,&nmat,&ntmat_,&norien,prestr,
	    &iprestr,&isolver,fei,veold,timepar,
	    xmodal,filab,jout,&nlabel,&idrct,
	    jmax,&iexpl,alpha,iamboun,plicon,nplicon,
	    plkcon,nplkcon,&iplas,&npmat_,mi,&nk_,trab,inotr,&ntrans,
	    ikboun,ilboun,ikmpc,ilmpc,ics,dcs,&ncs_,&namtot_,cs,&nstate_,
	    &ncmat_,&iumat,&mcs,labmpc,iponor,xnor,knor,thickn,thicke,
	    ikforc,ilforc,offset,iponoel,inoel,rig,infree,nshcon,shcon,
            cocon,ncocon,physcon,&nflow,
            ctrl,&maxlenmpc,&ne1d,&ne2d,&nener,vold,nodebounold,
            ndirbounold,xbounold,xforcold,xloadold,t1old,eme,
            sti,ener,xstate,jobnamec,irstrt,&ttime,
            qaold,output,typeboun,inpc,ipoinp,inp,tieset,tietol,
            &ntie,fmpc,cbody,ibody,xbody,&nbody,&nbody_,xbodyold,&nam_,
	    ielprop,&nprop,&nprop_,prop,&itpamp,&iviewfile,ipoinpc,
	    &nslavs,t0g,t1g,&network,&cyclicsymmetry,idefforc,idefload,
	    idefbody,&mortar,&ifacecount,islavsurf,pslavsurf,clearini,
	    heading,&iaxial,&nobject,objectset,&nprint_,iuel,&nuel_,
	    nodempcref,coefmpcref,ikmpcref,&memmpcref_,&mpcfreeref,
	    &maxlenmpcref,&memmpc_,&isens,&namtot,&nstam,dacon,vel,&nef,
	    velo,veloo,ne2boun,itempuser));

#ifdef CALCULIX_EXTERNAL_BEHAVIOURS_SUPPORT
  for(i=0;i!=nmat;++i){
    calculix_registerExternalBehaviour(matname+80*i);
  }
#endif /* CALCULIX_EXTERNAL_BEHAVIOURS_SUPPORT */
  
  if((istep==1)&&(mortar==-1)){mortar=0;}else{icontact=1;}

  nload0=nload;SFREE(idefforc);SFREE(idefload);SFREE(idefbody);

  if(nheading_>=0){
      writeheading(jobnamec,heading,&nheading_);
      SFREE(heading);
      nheading_=-1;
  }

  if((abs(nmethod)!=1)||(iperturb[0]<2))icascade=0;

//	FORTRAN(writeboun,(nodeboun,ndirboun,xboun,typeboun,&nboun));

  if(istat<0) break;

  if(istep == 1) {

    SFREE(iuel);

  /* tied contact constraints: generate appropriate MPC's */

    tiedcontact(&ntie, tieset, &nset, set,istartset, iendset, ialset,
       lakon, ipkon, kon,tietol,&nmpc, &mpcfree, &memmpc_,
       &ipompc, &labmpc, &ikmpc, &ilmpc,&fmpc, &nodempc, &coefmpc,
       ithermal, co, vold,&nef,&nmpc_,mi,&nk,&istep,ikboun,&nboun,
       kind1,kind2);

    /* reallocating space in the first step */

    /* allocating and initializing fields pointing to the previous step */

    RENEW(vold,double,mt*nk);
    NNEW(sti,double,6*mi[0]*ne);

    /* strains */

    NNEW(eme,double,6*mi[0]*ne);

    /* residual stresses/strains */

    if(iprestr==1) {
	RENEW(prestr,double,6*mi[0]*ne);
	for(i=0;i<ne;i++){
	    for(j=0;j<mi[0];j++){
		for(k=0;k<6;k++){
		    sti[6*mi[0]*i+6*j+k]=prestr[6*mi[0]*i+6*j+k];
		}
	    }
	}
    }
    else if(iprestr==2){
	RENEW(prestr,double,6*mi[0]*ne);
	for(i=0;i<ne;i++){
	    for(j=0;j<mi[0];j++){
		for(k=0;k<6;k++){
		    eme[6*mi[0]*i+6*j+k]=prestr[6*mi[0]*i+6*j+k];
		}
	    }
	}
    }
    else {
	SFREE(prestr);
    }

    NNEW(nodebounold,ITG,nboun);
    NNEW(ndirbounold,ITG,nboun);
    NNEW(xbounold,double,nboun);
    NNEW(xforcold,double,nforc);
    NNEW(xloadold,double,2*nload);

    /* initial temperatures: store in the "old" boundary conditions */

    if(ithermal[0]>1){
	for(i=0;i<nboun;i++){
	    if(strcmp1(&typeboun[i],"F")==0) continue;
	    if(ndirboun[i]==0){
		xbounold[i]=vold[mt*(nodeboun[i]-1)];
	    }
	}
    }

    /* initial temperatures: store in the "old" temperature field */

    if(ithermal[0]!=0){
      NNEW(t1old,double,nk);
      for(i=0;i<nk;i++) t1old[i]=t0[i];
    }

    /* element definition */

    RENEW(kon,ITG,nkon);
    RENEW(ipkon,ITG,ne);
    RENEW(lakon,char,8*ne);

    /* property cards */

    if(nprop>0){
	RENEW(ielprop,ITG,ne);
	RENEW(prop,double,nprop);
    }else{
	SFREE(ielprop);SFREE(prop);
    }

    /* fields for 1-D and 2-D elements */

    if((ne1d!=0)||(ne2d!=0)){
	RENEW(iponor,ITG,2*nkon);
	RENEW(xnor,double,infree[0]);
	RENEW(knor,ITG,infree[1]);
	SFREE(thickn);
	RENEW(thicke,double,mi[2]*nkon);
	RENEW(offset,double,2*ne);
	RENEW(inoel,ITG,3*(infree[2]-1));
	RENEW(iponoel,ITG,infree[3]);
	RENEW(rig,ITG,infree[3]);
	RENEW(ne2boun,ITG,2*infree[3]);
    }

    /* set definitions */ 

    RENEW(set,char,81*nset);
    RENEW(istartset,ITG,nset);
    RENEW(iendset,ITG,nset);
    RENEW(ialset,ITG,nalset);

    /* material properties */

    RENEW(elcon,double,(ncmat_+1)*ntmat_*nmat);
    RENEW(nelcon,ITG,2*nmat);

    RENEW(rhcon,double,2*ntmat_*nmat);
    RENEW(nrhcon,ITG,nmat);

    if(ndamp>0){RENEW(dacon,double,nmat);}

    RENEW(shcon,double,4*ntmat_*nmat);
    RENEW(nshcon,ITG,nmat);

    RENEW(cocon,double,7*ntmat_*nmat);
    RENEW(ncocon,ITG,2*nmat);

    RENEW(alcon,double,7*ntmat_*nmat);
    RENEW(nalcon,ITG,2*nmat);
    RENEW(alzero,double,nmat);

    RENEW(matname,char,80*nmat);
    RENEW(ielmat,ITG,mi[2]*ne);

    /* allocating space for the state variables */

    if(mortar!=1){
	RENEW(xstate,double,nstate_*mi[0]*(ne+nslavs));
	for(i=nxstate;i<nstate_*mi[0]*(ne+nslavs);i++){xstate[i]=0.;}
    }else if(mortar==1){
	RENEW(xstate,double,nstate_*mi[0]*(ne+nintpoint));
	for(i=nxstate;i<nstate_*mi[0]*(ne+nintpoint);i++){xstate[i]=0.;}
    }

    /* next statements for plastic materials and nonlinear springs */

    if(npmat_>0){
	RENEW(plicon,double,(2*npmat_+1)*ntmat_*nmat);
	RENEW(nplicon,ITG,(ntmat_+1)*nmat);
	RENEW(plkcon,double,(2*npmat_+1)*ntmat_*nmat);
	RENEW(nplkcon,ITG,(ntmat_+1)*nmat);
    }

    /* material orientation */

    if(norien > 0) {
      RENEW(orname,char,80*norien);
      RENEW(ielorien,ITG,mi[2]*ne);
      RENEW(orab,double,7*norien);
    }
    else {
	SFREE(orname);
	SFREE(ielorien);
	SFREE(orab);
    }

    /* amplitude definitions */

    if(nam > 0) {
      RENEW(amname,char,80*nam);
      RENEW(namta,ITG,3*nam);
      RENEW(amta,double,2*namta[3*nam-2]);
    }
    else {
      SFREE(amname);
      SFREE(amta);
      SFREE(namta);
      SFREE(iamforc);
      SFREE(iamload);
      SFREE(iamboun);
    }

    if(ntrans > 0){
      RENEW(trab,double,7*ntrans);
    }
    else{SFREE(trab);SFREE(inotr);}

    if(ithermal[0] == 0){SFREE(t0);SFREE(t1);SFREE(t0g);SFREE(t1g);}
    if((ithermal[0] == 0)||(nam<=0)){SFREE(iamt1);}

    if(ncs_>0){
      RENEW(ics,ITG,ncs_);
//      SFREE(dcs);
    }else if(npt_>0){SFREE(ics);}
    SFREE(dcs);

    if(mcs>0){
	RENEW(cs,double,17*mcs);
    }else{
	SFREE(cs);
    }

 }else{

    /* reallocating space in all but the first step (>1) */

    RENEW(vold,double,mt*nk);

    /* if the SPC boundary conditions were changed in the present step,
       they have to be rematched with those in the last step. Removed SPC 
       boundary conditions do not appear any more (this is different from
       forces and loads, where removed forces or loads are reset to zero;
       a removed SPC constraint does not have a numerical value any more) */
       
    NNEW(reorder,double,nboun);
    NNEW(nreorder,ITG,nboun);
    if(nbounold<nboun){
      RENEW(xbounold,double,nboun);
      RENEW(nodebounold,ITG,nboun);
      RENEW(ndirbounold,ITG,nboun);
    }
    FORTRAN(spcmatch,(xboun,nodeboun,ndirboun,&nboun,xbounold,nodebounold,
		      ndirbounold,&nbounold,ikboun,ilboun,vold,reorder,nreorder,
                      mi,typeboun));
    RENEW(xbounold,double,nboun);
    RENEW(nodebounold,ITG,nboun);
    RENEW(ndirbounold,ITG,nboun);
    SFREE(reorder); SFREE(nreorder);

    /* for additional forces or loads in the present step, the
       corresponding slots in the force and load fields of the
       previous steps are initialized */

    RENEW(xforcold,double,nforc);
    for(i=nforcold;i<nforc;i++) xforcold[i]=0;

    RENEW(xloadold,double,2*nload);
    for(i=2*nloadold;i<2*nload;i++) xloadold[i]=0;

    if(ithermal[0]!=0){
      RENEW(t1old,double,nk);
    }

    if(nam > 0) {
      RENEW(amname,char,80*nam);
      RENEW(namta,ITG,3*nam);
      RENEW(amta,double,2*namta[3*nam-2]);
    }
    
  }

  /* reallocating fields for all steps (>=1) */

  RENEW(co,double,3*nk);

  RENEW(nodeboun,ITG,nboun);
  RENEW(ndirboun,ITG,nboun);
  RENEW(typeboun,char,nboun+1);
  RENEW(xboun,double,nboun);
  RENEW(ikboun,ITG,nboun);
  RENEW(ilboun,ITG,nboun);
    
  RENEW(nodeforc,ITG,2*nforc);
  RENEW(ndirforc,ITG,nforc);
  RENEW(xforc,double,nforc);
  RENEW(ikforc,ITG,nforc);
  RENEW(ilforc,ITG,nforc);

  /* temperature loading */
  
  if(ithermal[0] != 0){
      RENEW(t0,double,nk);
      RENEW(t1,double,nk);
      if((ne1d!=0)||(ne2d!=0)){
	  RENEW(t0g,double,2*nk);
	  RENEW(t1g,double,2*nk);
      }
      if(nam > 0) {RENEW(iamt1,ITG,nk);}
  }

  RENEW(nelemload,ITG,2*nload);
  RENEW(sideload,char,20*nload);
  RENEW(xload,double,2*nload);

  RENEW(cbody,char,81*nbody);
  RENEW(ibody,ITG,3*nbody);
  RENEW(xbody,double,7*nbody);
  RENEW(xbodyold,double,7*nbody);

  RENEW(ipompc,ITG,nmpc);
  RENEW(labmpc,char,20*nmpc+1);
  RENEW(ikmpc,ITG,nmpc);
  RENEW(ilmpc,ITG,nmpc);
  RENEW(fmpc,double,nmpc);

  /* energy */

  if((nener==1)&&(nenerold==0)){
    NNEW(ener,double,mi[0]*ne*2);
    if((istep>1)&&(iperturb[0]>1)){
      printf(" *ERROR in CalculiX: in nonlinear calculations\n");
      printf("        energy output must be selected in the first step\n\n");
      FORTRAN(stop,());
    }
  }

  /* initial velocities and accelerations */

  if((nmethod==4)||(nmethod==5)||(nmethod==8)||(nmethod==9)||
     ((abs(nmethod)==1)&&(iperturb[0]>=2))){
    RENEW(veold,double,mt*nk);
  }
  else {SFREE(veold);}

  if((nmethod == 4)&&(iperturb[0]>1)) {
    NNEW(accold,double,mt*nk);
    }

  if(nam > 0) {
    RENEW(iamforc,ITG,nforc);
    RENEW(iamload,ITG,2*nload);
    RENEW(iamboun,ITG,nboun);
  }

  /* generate force convection elements */

//  if(network==1){
  if(network>0){
      ne0=ne;nkon0=nkon;nload1=nload;
      RENEW(ipkon,ITG,ne+nload);
      RENEW(lakon,char,8*(ne+nload));
      RENEW(kon,ITG,nkon+9*nload);
      NNEW(inodesd,ITG,nk);
      RENEW(nelemload,ITG,4*nload);
      RENEW(sideload,char,40*nload);
      
      FORTRAN(genadvecelem,(inodesd,ipkon,&ne,lakon,kon,&nload,
			    sideload,nelemload,&nkon,&network));
      
      SFREE(inodesd);
      RENEW(ipkon,ITG,ne);
      RENEW(lakon,char,8*ne);
      RENEW(kon,ITG,nkon);
      RENEW(sti,double,6*mi[0]*ne);
      RENEW(eme,double,6*mi[0]*ne);
      if(iprestr>0) RENEW(prestr,double,6*mi[0]*ne);
      if(nprop>0) RENEW(ielprop,ITG,ne);
      if((ne1d!=0)||(ne2d!=0)) RENEW(offset,double,2*ne);
      RENEW(nelemload,ITG,2*nload);
      RENEW(sideload,char,20*nload);
      RENEW(xload,double,2*nload);
      RENEW(xloadold,double,2*nload);
      if(nam>0){
	  RENEW(iamload,ITG,2*nload);
	  for(i=2*nload1;i<2*nload;i++)iamload[i]=0;
      }
      if(nener==1)RENEW(ener,double,mi[0]*ne*2);
      if(norien>0)RENEW(ielorien,ITG,mi[2]*ne);
      RENEW(ielmat,ITG,mi[2]*ne);
      for(i=mi[2]*ne0;i<mi[2]*ne;i++)ielmat[i]=1;
  }

  if(ntrans > 0){
    RENEW(inotr,ITG,2*nk);
  }
  
  /*   calling the user routine ufaceload (can be empty) */

  if(ithermal[1]>=2){
      NNEW(sideloadtemp,char,20*nload);
      for(i=0;i<nload;i++){
	  strcpy1(&sideloadtemp[20*i],&sideload[20*i],20);
	  if((strcmp1(&sideload[20*i]," ")==0)&&
	     (strcmp1(&sideload[20*i+1]," ")!=0)){
	      strcpy1(&sideloadtemp[20*i],"F",1);
	  }
      }
      FORTRAN(ufaceload,(co,ipkon,kon,lakon,&nboun,nodeboun,
              nelemload,sideloadtemp,&nload,&ne,&nk));
      SFREE(sideloadtemp);
  }

      /* storing the undecascaded MPC's */

  if(mpcfreeref==-1){
      memmpcref_=memmpc_;mpcfreeref=mpcfree;maxlenmpcref=maxlenmpc;
      NNEW(nodempcref,ITG,3*memmpc_);memcpy(nodempcref,nodempc,sizeof(ITG)*3*memmpc_);
      NNEW(coefmpcref,double,memmpc_);memcpy(coefmpcref,coefmpc,sizeof(double)*memmpc_);
      NNEW(ikmpcref,ITG,nmpc);memcpy(ikmpcref,ikmpc,sizeof(ITG)*nmpc);
  }
 
  /* decascading MPC's only necessary if MPC's changed */

  if(((istep == 1)||(ntrans>0)||(mpcend<0)||(nk!=nkold)||(nmpc!=nmpcold))&&(icascade==0)) {
//  if(icascade==0) {

    /* decascading the MPC's */

    printf(" Decascading the MPC's\n\n");

    callfrommain=1;
    cascade(ipompc,&coefmpc,&nodempc,&nmpc,
	    &mpcfree,nodeboun,ndirboun,&nboun,ikmpc,
	    ilmpc,ikboun,ilboun,&mpcend,
	    labmpc,&nk,&memmpc_,&icascade,&maxlenmpc,
            &callfrommain,iperturb,ithermal);
  }

  /* determining the matrix structure: changes if SPC's have changed */
  
  if((icascade==0)&&(nmethod<8)) printf(" Determining the structure of the matrix:\n");
  
  NNEW(nactdof,ITG,mt*nk);  
  NNEW(mast1,ITG,nzs[1]);
  NNEW(irow,ITG,1);
  
  if((mcs==0)||(cs[1]<0)){
      
      NNEW(icol,ITG,mt*nk);
      NNEW(jq,ITG,mt*nk+1);
      NNEW(ipointer,ITG,mt*nk);
      
      if((icascade==0)&&((nmethod<8)||(nmethod>10))){
	if((nmethod==11)||(nmethod==13)){nmethodl=2;}else{nmethodl=nmethod;}
	  mastruct(&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,&nboun,ipompc,
		   nodempc,&nmpc,nactdof,icol,jq,&mast1,&irow,&isolver,neq,
		   ikmpc,ilmpc,ipointer,nzs,&nmethodl,ithermal,
                   ikboun,ilboun,iperturb,mi,&mortar,typeboun,labmpc,
		   &iit,&icascade,&network,&iexpl);
      }
      else{neq[0]=1;neq[1]=1;neq[2]=1;}
  }
  else{
      
      NNEW(icol,ITG,8*nk);
      NNEW(jq,ITG,8*nk+1);
      NNEW(ipointer,ITG,8*nk);
      
      if(nmethod==13){nmethodl=2;}else{nmethodl=nmethod;}
      mastructcs(&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,&nboun,
		 ipompc,nodempc,&nmpc,nactdof,icol,jq,&mast1,&irow,&isolver,
		 neq,ikmpc,ilmpc,ipointer,nzs,&nmethodl,
		 ics,cs,labmpc,&mcs,mi,&mortar);
  }
  
  SFREE(ipointer);SFREE(mast1);
  if((icascade==0)&&(nmethod<8))RENEW(irow,ITG,nzs[2]);

  /* nmethod=1: static analysis   */
  /* nmethod=2: frequency analysis  */
  /* nmethod=3: buckling analysis */
  /* nmethod=4: (linear or nonlinear) dynamic analysis */
  /* nmethod=5: steady state dynamics analysis */
  /* nmethod=6: Coriolis frequency calculation */
  /* nmethod=7: flutter frequency calculation */
  /* nmethod=8:  magnetostatics */
  /* nmethod=9:  magnetodynamics */
  /* nmethod=10: electromagnetic eigenvalue problems */
  /* nmethod=11: superelement creation */
  /* nmethod=12: sensitivity analysis  */
  /* nmethod=13: Green function calculation */
     

  if((nmethod<=1)||(nmethod==11)||((iperturb[0]>1)&&(nmethod<8)))
    {
	if(iperturb[0]<2){
	
	mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
	mpcinfo[3]=maxlenmpc;

	if(icascade!=0){
	    printf(" *ERROR in CalculiX: the matrix structure may");
	    printf("        change due to nonlinear equations;");
	    printf("        a purely linear calculation is not");
	    printf("        feasible; use NLGEOM on the *STEP card.");
	    FORTRAN(stop,());
	}

	linstatic(co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,ndirboun,xboun,&nboun, 
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, nelemload,sideload,xload,&nload, 
	     nactdof,&icol,jq,&irow,neq,&nzl,&nmethod,ikmpc, 
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, vold,iperturb,sti,nzs,
	     &kode,filab,eme,&iexpl,plicon,
             nplicon,plkcon,nplkcon,&xstate,&npmat_,matname,
	     &isolver,mi,&ncmat_,&nstate_,cs,&mcs,&nkon,&ener,
             xbounold,xforcold,xloadold,amname,amta,namta,
             &nam,iamforc,iamload,iamt1,iamboun,&ttime,
             output,set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,trab,inotr,&ntrans,fmpc,cbody,ibody,xbody,&nbody,
	     xbodyold,timepar,thicke,jobnamec,tieset,&ntie,&istep,&nmat,
	     ielprop,prop,typeboun,&mortar,mpcinfo,tietol,ics,&icontact,
	     orname,itempuser);

	for(i=0;i<3;i++){nzsprevstep[i]=nzs[i];}

	memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
        maxlenmpc=mpcinfo[3];

      }

      else{

	mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
	mpcinfo[3]=maxlenmpc;

	nonlingeo(&co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,ndirboun,xboun,&nboun, 
	     &ipompc,&nodempc,&coefmpc,&labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc,&nelemload,&sideload,xload,&nload, 
	     nactdof,&icol,jq,&irow,neq,&nzl,&nmethod,&ikmpc, 
	     &ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, 
	     &vold,iperturb,sti,nzs,&kode,filab,&idrct,jmax,
	     jout,timepar,eme,xbounold,xforcold,xloadold,
	     veold,accold,amname,amta,namta,
	     &nam,iamforc,&iamload,iamt1,alpha,
             &iexpl,iamboun,plicon,nplicon,plkcon,nplkcon,
	     &xstate,&npmat_,&istep,&ttime,matname,qaold,mi,
	     &isolver,&ncmat_,&nstate_,&iumat,cs,&mcs,&nkon,&ener,
	     mpcinfo,output,
             shcon,nshcon,cocon,ncocon,physcon,&nflow,ctrl,
             set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,ikforc,ilforc,trab,inotr,&ntrans,&fmpc,
             cbody,ibody,xbody,&nbody,xbodyold,ielprop,prop,
	     &ntie,tieset,&itpamp,&iviewfile,jobnamec,tietol,&nslavs,thicke,
	     ics,&nintpoint,&mortar,
	     &ifacecount,typeboun,&islavsurf,&pslavsurf,&clearini,&nmat,
	     xmodal,&iaxial,&inext,&nprop,&network,orname,vel,&nef,
	     velo,veloo,energy,itempuser);

	memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
        maxlenmpc=mpcinfo[3];

	for(i=0;i<3;i++){nzsprevstep[i]=nzs[i];}

      }
    }else if((nmethod==2)||(nmethod==13)){
      
      /* FREQUENCY ANALYSIS */
      
      if((mcs==0)||(cs[1]<0)){
#ifdef ARPACK
	  
	  mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
	  mpcinfo[3]=maxlenmpc;
	  
	  arpack(co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,ndirboun,xboun,&nboun, 
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, nelemload,sideload,xload,&nload, 
	     nactdof,icol,jq,&irow,neq,&nzl,&nmethod,ikmpc, 
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     shcon,nshcon,cocon,ncocon,
             alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr,vold,iperturb,sti,nzs,
	     &kode,mei,fei,filab,
	     &iexpl,plicon,nplicon,plkcon,nplkcon,
	     &xstate,&npmat_,matname,mi,&ncmat_,&nstate_,&ener,jobnamec,
             output,set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,&isolver,trab,inotr,&ntrans,&ttime,fmpc,cbody,
	     ibody,xbody,&nbody,thicke,&nslavs,tietol,&nkon,mpcinfo,
	     &ntie,&istep,&mcs,ics,tieset,cs,&nintpoint,&mortar,&ifacecount,
	     &islavsurf,&pslavsurf,&clearini,&nmat,typeboun,ielprop,prop,
             orname);

	  memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
	  maxlenmpc=mpcinfo[3];

	  for(i=0;i<3;i++){nzsprevstep[i]=nzs[i];}

#else
	  printf("*ERROR in CalculiX: the ARPACK library is not linked\n\n");
	  FORTRAN(stop,());
#endif

      }else{

#ifdef ARPACK

	  mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
	  mpcinfo[3]=maxlenmpc;
	  
	  arpackcs(co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,ndirboun,
             xboun,&nboun, 
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, nelemload,sideload,xload,&nload, 
	     nactdof,icol,jq,&irow,neq,&nzl,&nmethod,ikmpc, 
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
             alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, 
	     vold,iperturb,sti,nzs,&kode,mei,fei,filab,
	     &iexpl,plicon,nplicon,plkcon,nplkcon,
	     &xstate,&npmat_,matname,mi,ics,cs,&mpcend,&ncmat_,
             &nstate_,&mcs,&nkon,jobnamec,output,set,&nset,istartset,
             iendset,ialset,&nprint,prlab,
             prset,&nener,&isolver,trab,inotr,&ntrans,&ttime,fmpc,cbody,
             ibody,xbody,&nbody,&nevtot,thicke,&nslavs,tietol,mpcinfo,
	     &ntie,&istep,tieset,&nintpoint,&mortar,&ifacecount,&islavsurf,
	     &pslavsurf,&clearini,&nmat,typeboun,ielprop,prop,orname);

	  memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
	  maxlenmpc=mpcinfo[3];

	  for(i=0;i<3;i++){nzsprevstep[i]=nzs[i];}

#else
	  printf("*ERROR in CalculiX: the ARPACK library is not linked\n\n");
	  FORTRAN(stop,());
#endif

      }
  }else if(nmethod==3){
    
#ifdef ARPACK
	arpackbu(co,&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,xboun,&nboun, 
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, 
	     nelemload,sideload,xload,&nload, 
	     nactdof,icol,jq,irow,neq,&nzl,&nmethod,ikmpc, 
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
             alcon,nalcon,alzero,ielmat,ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, 
	     vold,iperturb,sti,nzs,&kode,mei,fei,filab,
	     eme,&iexpl,plicon,nplicon,plkcon,nplkcon,
	     xstate,&npmat_,matname,mi,&ncmat_,&nstate_,ener,output,
             set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,&isolver,trab,inotr,&ntrans,&ttime,fmpc,cbody,
	     ibody,xbody,&nbody,thicke,jobnamec,&nmat,ielprop,prop,
	     orname,typeboun);
#else
            printf("*ERROR in CalculiX: the ARPACK library is not linked\n\n");
            FORTRAN(stop,());
#endif
    }
  else if(nmethod==4)
    {
	if((ne1d!=0)||(ne2d!=0)){
	    printf(" *WARNING: 1-D or 2-D elements may cause problems in modal dynamic calculations\n");
	    printf("           ensure that point loads defined in a *MODAL DYNAMIC step\n");
	    printf("           and applied to nodes belonging to 1-D or 2-D elements have been\n");
	    printf("           applied to the same nodes in the preceding FREQUENCY step with\n");
	    printf("           magnitude zero; look at example shellf.inp for a guideline.\n\n");}

      printf(" Composing the dynamic response from the eigenmodes\n\n");

      dyna(&co,&nk,&kon,&ipkon,&lakon,&ne,&nodeboun,&ndirboun,&xboun,&nboun,
	    &ipompc,&nodempc,&coefmpc,&labmpc,&nmpc,nodeforc,ndirforc,xforc,&nforc,
	    nelemload,sideload,xload,&nload,
	    &nactdof,neq,&nzl,icol,irow,&nmethod,&ikmpc,&ilmpc,&ikboun,&ilboun,
            elcon,nelcon,rhcon,nrhcon,cocon,ncocon,
            alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,&t0,
	    &t1,ithermal,prestr,&iprestr,&vold,iperturb,&sti,nzs,
	    timepar,xmodal,&veold,amname,amta,
	    namta,&nam,iamforc,iamload,&iamt1,
	    jout,&kode,filab,&eme,xforcold,xloadold,
            &t1old,&iamboun,&xbounold,&iexpl,plicon,
            nplicon,plkcon,nplkcon,&xstate,&npmat_,matname,
            mi,&ncmat_,&nstate_,&ener,jobnamec,&ttime,set,&nset,
            istartset,iendset,&ialset,&nprint,prlab,
            prset,&nener,trab,&inotr,&ntrans,&fmpc,cbody,ibody,xbody,&nbody,
            xbodyold,&istep,&isolver,jq,output,&mcs,&nkon,&mpcend,ics,cs,
	    &ntie,tieset,&idrct,jmax,ctrl,&itpamp,tietol,&nalset,
	    ikforc,ilforc,thicke,&nslavs,&nmat,typeboun,ielprop,prop,orname);
    }
  else if(nmethod==5)
    {
	  if((ne1d!=0)||(ne2d!=0)){
	      printf(" *WARNING: 1-D or 2-D elements may cause problems in steady state calculations\n");
	      printf("           ensure that point loads defined in a *STEADY STATE DYNAMICS step\n");
	      printf("           and applied to nodes belonging to 1-D or 2-D elements have been\n");
	      printf("           applied to the same nodes in the preceding FREQUENCY step with\n");
	      printf("           magnitude zero; look at example shellf.inp for a guideline.\n\n");}

      printf(" Composing the steady state response from the eigenmodes\n\n");

      steadystate(&co,&nk,&kon,&ipkon,&lakon,&ne,&nodeboun,&ndirboun,&xboun,&nboun,
	    &ipompc,&nodempc,&coefmpc,&labmpc,&nmpc,nodeforc,ndirforc,xforc,&nforc,
	    nelemload,sideload,xload,&nload,
	    &nactdof,neq,&nzl,icol,irow,&nmethod,&ikmpc,&ilmpc,&ikboun,&ilboun,
            elcon,nelcon,rhcon,nrhcon,cocon,ncocon,
            alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,&t0,
	    &t1,ithermal,prestr,&iprestr,&vold,iperturb,sti,nzs,
	    timepar,xmodal,&veold,amname,amta,
	    namta,&nam,iamforc,iamload,&iamt1,
	    jout,&kode,filab,&eme,xforcold,xloadold,
            &t1old,&iamboun,&xbounold,&iexpl,plicon,
            nplicon,plkcon,nplkcon,xstate,&npmat_,matname,
            mi,&ncmat_,&nstate_,&ener,jobnamec,&ttime,set,&nset,
            istartset,iendset,ialset,&nprint,prlab,
            prset,&nener,trab,&inotr,&ntrans,&fmpc,cbody,ibody,xbody,&nbody,
	    xbodyold,&istep,&isolver,jq,output,&mcs,&nkon,ics,cs,&mpcend,
	    ctrl,ikforc,ilforc,thicke,&nmat,typeboun,ielprop,prop,orname,
	    &ndamp,dacon);
    }
  else if((nmethod==6)||(nmethod==7))
    {

      printf(" Composing the complex eigenmodes from the real eigenmodes\n\n");

      complexfreq(&co,&nk,&kon,&ipkon,&lakon,&ne,&nodeboun,&ndirboun,&xboun,&nboun,
	    &ipompc,&nodempc,&coefmpc,&labmpc,&nmpc,nodeforc,ndirforc,xforc,&nforc,
	    nelemload,sideload,xload,&nload,
	    &nactdof,neq,&nzl,icol,irow,&nmethod,&ikmpc,&ilmpc,&ikboun,&ilboun,
            elcon,nelcon,rhcon,nrhcon,cocon,ncocon,
            alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,&t0,
	    &t1,ithermal,prestr,&iprestr,&vold,iperturb,&sti,nzs,
	    timepar,xmodal,&veold,amname,amta,
	    namta,&nam,iamforc,iamload,&iamt1,
	    jout,&kode,filab,&eme,xforcold,xloadold,
            &t1old,&iamboun,&xbounold,&iexpl,plicon,
            nplicon,plkcon,nplkcon,xstate,&npmat_,matname,
            mi,&ncmat_,&nstate_,&ener,jobnamec,&ttime,set,&nset,
            istartset,iendset,&ialset,&nprint,prlab,
            prset,&nener,trab,&inotr,&ntrans,&fmpc,cbody,ibody,xbody,&nbody,
            xbodyold,&istep,&isolver,jq,output,&mcs,&nkon,&mpcend,ics,cs,
	    &ntie,tieset,&idrct,jmax,ctrl,&itpamp,tietol,&nalset,
	    ikforc,ilforc,thicke,jobnamef,mei,&nmat,ielprop,prop,orname,
            typeboun);
    }
  else if((nmethod>7)&&(nmethod<12)){

	mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
	mpcinfo[3]=maxlenmpc;

	electromagnetics(&co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,
             ndirboun,xboun,&nboun, 
	     &ipompc,&nodempc,&coefmpc,&labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc,&nelemload,&sideload,xload,&nload, 
	     nactdof,&icol,&jq,&irow,neq,&nzl,&nmethod,&ikmpc, 
	     &ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, 
	     &vold,iperturb,sti,nzs,&kode,filab,&idrct,jmax,
	     jout,timepar,eme,xbounold,xforcold,xloadold,
	     veold,accold,amname,amta,namta,
	     &nam,iamforc,&iamload,iamt1,alpha,
             &iexpl,iamboun,plicon,nplicon,plkcon,nplkcon,
	     &xstate,&npmat_,&istep,&ttime,matname,qaold,mi,
	     &isolver,&ncmat_,&nstate_,&iumat,cs,&mcs,&nkon,&ener,
	     mpcinfo,output,
             shcon,nshcon,cocon,ncocon,physcon,&nflow,ctrl,
             &set,&nset,&istartset,&iendset,&ialset,&nprint,prlab,
             prset,&nener,ikforc,ilforc,trab,inotr,&ntrans,&fmpc,
             cbody,ibody,xbody,&nbody,xbodyold,ielprop,prop,
	     &ntie,&tieset,&itpamp,&iviewfile,jobnamec,&tietol,&nslavs,thicke,
	     ics,&nalset,&nmpc_,&nmat,typeboun,&iaxial,&nload_,&nprop,
	     &network,orname);

	memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
        maxlenmpc=mpcinfo[3];
  }

  else if(nmethod==12){
  
	sensitivity(co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,ndirboun,
	     xboun,&nboun, ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,
             ndirforc,xforc,&nforc, nelemload,sideload,xload,&nload, 
	     nactdof,icol,jq,&irow,neq,&nzl,&nmethod,ikmpc, 
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, vold,iperturb,sti,nzs,
	     &kode,filab,eme,&iexpl,plicon,
             nplicon,plkcon,nplkcon,&xstate,&npmat_,matname,
	     &isolver,mi,&ncmat_,&nstate_,cs,&mcs,&nkon,&ener,
             xbounold,xforcold,xloadold,amname,amta,namta,
             &nam,iamforc,iamload,iamt1,iamboun,&ttime,
             output,set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,trab,inotr,&ntrans,fmpc,cbody,ibody,xbody,&nbody,
	     xbodyold,timepar,thicke,jobnamec,tieset,&ntie,&istep,&nmat,
	     ielprop,prop,typeboun,&mortar,mpcinfo,tietol,ics,&icontact,
	     &nobject,&objectset,&istat,orname,nzsprevstep,&nlabel,physcon,
	     jobnamef,iponor,knor,&ne2d,iponoel,inoel,&mpcend);
  }

  SFREE(nactdof);SFREE(icol);SFREE(jq);SFREE(irow);

  /* reset tempuserflag */
  
  itempuser[0]=0;

  /* deleting the perturbation loads and temperatures */

  if((iperturb[0] == 1)&&(nmethod==3)) {
      nforc=0;
      nload=0;
      nbody=0;
      if(ithermal[0] == 1) {
	  for(k=0;k<nk;++k){
	      t1[k]=t0[k];
	  }
      }
  }else{
      nbounold=nboun;
      for (i=0;i<nboun;i++) {
	  nodebounold[i]=nodeboun[i];
	  ndirbounold[i]=ndirboun[i];
      }
      nforcold=nforc;
      nloadold=nload;
      nbodyold=nbody;
      
    /* resetting the amplitude to none except for time=total time amplitudes */

      if(nam > 0) {
	  for (i=0;i<nboun;i++) {
	      if(iamboun[i]>0){
		  if(namta[3*iamboun[i]-1]>0){
		      iamboun[i]=0;
		      xboun[i]=xbounold[i];}
	      }
	  }
	  for (i=0;i<nforc;i++){
	      if(iamforc[i]>0){
		  if(namta[3*iamforc[i]-1]>0){
		      iamforc[i]=0;
		      xforc[i]=xforcold[i];}
	      }
	  }
	  for (i=0;i<2*nload;i++){
	      if(iamload[i]>0){
		  if(namta[3*iamload[i]-1]>0){
		      iamload[i]=0;
		      xload[i]=xloadold[i];}
	      }
	  }
	  for (i=1;i<3*nbody;i=i+3){
	      if(ibody[i]>0){
		  if(namta[3*ibody[i]-1]>0){
		      ibody[i]=0;
		      xbody[7*(i-1)/3]=xbodyold[7*(i-1)/3];}
	      }
	  }
	  if(ithermal[0]==1) {
	      if(iamt1[i]>0){
		  if(namta[3*iamt1[i]-1]>0){
		      iamt1[i]=0;
		      t1[i]=t1old[i];}
	      }
	  }
      }
  }

  /* removing the advective elements, if any */

  if(network>0){
      ne=ne0;nkon=nkon0;
      RENEW(ipkon,ITG,ne);
      RENEW(lakon,char,8*ne);
      RENEW(kon,ITG,nkon);
      RENEW(sti,double,6*mi[0]*ne);
      RENEW(eme,double,6*mi[0]*ne);
      if(iprestr>0) RENEW(prestr,double,6*mi[0]*ne);
      if(nprop>0) RENEW(ielprop,ITG,ne);
      if((ne1d!=0)||(ne2d!=0)) RENEW(offset,double,2*ne);
      if(nener==1)RENEW(ener,double,mi[0]*ne*2);
      if(norien>0)RENEW(ielorien,ITG,mi[2]*ne);
      RENEW(ielmat,ITG,mi[2]*ne);

      /* reactivating the original load labels */

      for(i=nload-1;i>=nload0;i--){
	  if(strcmp2(&sideload[20*i],"                    ",20)==0){
	      iload=nelemload[2*i+1];
	      strcpy1(&sideload[20*(iload-1)],"F",1);
	  }
      }

  }

  nload=nload0;

  if((nmethod == 4)&&(iperturb[0]>1)) SFREE(accold);

  if(irstrt[0]>0){
    jrstrt++;
    if(jrstrt>=irstrt[0]){
      jrstrt=0;
      FORTRAN(restartwrite,(&istep,&nset,&nload,&nforc,&nboun,&nk,&ne,
        &nmpc,&nalset,&nmat,&ntmat_,&npmat_,&norien,&nam,&nprint, 
        mi,&ntrans,&ncs_,&namtot,&ncmat_,&mpcend,&maxlenmpc,&ne1d,
        &ne2d,&nflow,&nlabel,&iplas,&nkon,ithermal,&nmethod,iperturb,
        &nstate_,&nener,set,istartset,iendset,ialset,co,kon,ipkon,
        lakon,nodeboun,ndirboun,iamboun,xboun,ikboun,ilboun,ipompc,
        nodempc,coefmpc,labmpc,ikmpc,ilmpc,nodeforc,ndirforc,iamforc,
        xforc,ikforc,ilforc,nelemload,iamload,sideload,xload,
        elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
        alzero,plicon,nplicon,plkcon,nplkcon,orname,orab,ielorien,
        trab,inotr,amname,amta,namta,t0,t1,iamt1,veold,
        ielmat,matname,prlab,prset,filab,vold,nodebounold,
        ndirbounold,xbounold,xforcold,xloadold,t1old,eme,
        iponor,xnor,knor,thicke,offset,iponoel,inoel,rig,
        shcon,nshcon,cocon,ncocon,ics,
	sti,ener,xstate,jobnamec,infree,prestr,&iprestr,cbody,
	ibody,xbody,&nbody,xbodyold,&ttime,qaold,cs,&mcs,output,
	physcon,ctrl,typeboun,fmpc,tieset,&ntie,tietol,&nslavs,t0g,t1g,
	&nprop,ielprop,prop,&mortar,&nintpoint,&ifacecount,islavsurf,
	pslavsurf,clearini,irstrt,vel,&nef,velo,veloo,ne2boun));
    }
  } 
	  
}

FORTRAN(closefile,());

strcpy(fneig,jobnamec);
strcat(fneig,".frd");
if((f1=fopen(fneig,"ab"))==NULL){
    printf("*ERROR in frd: cannot open frd file for writing...");
    exit(0);
}
fprintf(f1," 9999\n");
fclose(f1);

/* deallocating the fields
   this section is addressed immediately after leaving calinput */

SFREE(ipoinpc);SFREE(inpc);SFREE(inp);SFREE(ipoinp);

if(ncs_>0) SFREE(ics);
if(mcs>0) SFREE(cs);
SFREE(tieset);SFREE(tietol);

SFREE(co);SFREE(kon);SFREE(ipkon);SFREE(lakon);

SFREE(nodeboun);SFREE(ndirboun);SFREE(typeboun);SFREE(xboun);SFREE(ikboun);
SFREE(ilboun);SFREE(nodebounold);SFREE(ndirbounold);SFREE(xbounold);

SFREE(ipompc);SFREE(labmpc);SFREE(ikmpc);SFREE(ilmpc);SFREE(fmpc);
SFREE(nodempc);SFREE(coefmpc);

SFREE(nodempcref);SFREE(coefmpcref);SFREE(ikmpcref);

SFREE(nodeforc);SFREE(ndirforc);SFREE(xforc);SFREE(ikforc);SFREE(ilforc);
SFREE(xforcold);

SFREE(nelemload);SFREE(sideload);SFREE(xload);SFREE(xloadold);

SFREE(cbody);SFREE(ibody);SFREE(xbody);SFREE(xbodyold);

if(nam>0){SFREE(iamboun);SFREE(iamforc);SFREE(iamload);SFREE(amname);
    SFREE(amta);SFREE(namta);}

SFREE(set);SFREE(istartset);SFREE(iendset);SFREE(ialset);

SFREE(elcon);SFREE(nelcon);SFREE(rhcon);SFREE(nrhcon);SFREE(shcon);SFREE(nshcon);
SFREE(cocon);SFREE(ncocon);SFREE(alcon);SFREE(nalcon);SFREE(alzero);
if(nprop>0){SFREE(ielprop);SFREE(prop);}
if(npmat_>0){SFREE(plicon);SFREE(nplicon);SFREE(plkcon);SFREE(nplkcon);}
if(ndamp>0){SFREE(dacon);}

if(norien>0){SFREE(orname);SFREE(orab);SFREE(ielorien);}
if(ntrans>0){SFREE(trab);SFREE(inotr);}
if(iprestr>0){SFREE(prestr);}

if(ithermal[0]!=0){
    SFREE(t0);SFREE(t1);SFREE(t1old);
    if(nam>0) SFREE(iamt1);
    if((ne1d!=0)||(ne2d!=0)){SFREE(t0g);SFREE(t1g);}
}

SFREE(prlab);SFREE(prset);SFREE(filab);SFREE(xmodal);

SFREE(ielmat);SFREE(matname);

SFREE(sti);SFREE(eme);SFREE(ener);SFREE(xstate);

SFREE(vold);SFREE(veold);SFREE(vel);SFREE(velo);SFREE(veloo);

if((ne1d!=0)||(ne2d!=0)){
    SFREE(iponor);SFREE(xnor);SFREE(knor);SFREE(thicke);SFREE(offset);
    SFREE(iponoel);SFREE(inoel);SFREE(rig);SFREE(ne2boun);
}

SFREE(islavsurf);
if(mortar==1){SFREE(pslavsurf);SFREE(clearini);}

if(nobject_>0){SFREE(objectset);}

#ifdef CALCULIX_MPI
MPI_Finalize();
#endif

#ifdef CALCULIX_EXTERNAL_BEHAVIOURS_SUPPORT
 calculix_freeExternalBehaviours();
#endif /* CALCULIX_EXTERNAL_BEHAVIOURS_SUPPORT */

 return 0;
      
}




