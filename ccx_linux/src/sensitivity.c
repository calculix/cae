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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"

void sensitivity(double *co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
		 ITG *ne,
		 ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
		 ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
		 ITG *nmpc,
		 ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
		 ITG *nelemload,char *sideload,double *xload,
		 ITG *nload,ITG *nactdof,
		 ITG *icol,ITG *jq,ITG **irowp,ITG *neq,ITG *nzl,
		 ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
		 ITG *ilboun,
		 double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
		 double *alcon,ITG *nalcon,double *alzero,ITG **ielmatp,
		 ITG **ielorienp,ITG *norien,double *orab,ITG *ntmat_,
		 double *t0,double *t1,double *t1old,
		 ITG *ithermal,double *prestr,ITG *iprestr,
		 double *vold,ITG *iperturb,double *sti,ITG *nzs, 
		 ITG *kode,char *filab,double *eme,
		 ITG *iexpl,double *plicon,ITG *nplicon,double *plkcon,
		 ITG *nplkcon,
		 double **xstatep,ITG *npmat_,char *matname,ITG *isolver,
		 ITG *mi,ITG *ncmat_,ITG *nstate_,double *cs,ITG *mcs,
		 ITG *nkon,double **enerp,double *xbounold,
		 double *xforcold,double *xloadold,
		 char *amname,double *amta,ITG *namta,
		 ITG *nam,ITG *iamforc,ITG *iamload,
		 ITG *iamt1,ITG *iamboun,double *ttime,char *output,
		 char *set,ITG *nset,ITG *istartset,
		 ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
		 char *prset,ITG *nener,double *trab,
		 ITG *inotr,ITG *ntrans,double *fmpc,ITG *ipobody,ITG *ibody,
		 double *xbody,ITG *nbody,double *xbodyold,double *timepar,
		 double *thicke,char *jobnamec,char *tieset,ITG *ntie,
		 ITG *istep,ITG *nmat,ITG *ielprop,double *prop,char *typeboun,
		 ITG *mortar,ITG *mpcinfo,double *tietol,ITG *ics,ITG *icontact,
		 ITG *nobject,char **objectsetp,ITG *istat,char *orname,
		 ITG *nzsprevstep,ITG *nlabel,double *physcon,char *jobnamef,
		 ITG *iponor2d,ITG *knor2d,ITG *ne2d,ITG *iponoel2d,ITG *inoel2d,
		 ITG *mpcend){
	     
  char description[13]="            ",*lakon=NULL,cflag[1]=" ",fneig[132]="",
    stiffmatrix[132]="",*lakonfa=NULL,*objectset=NULL;
       
  ITG *inum=NULL,k,*irow=NULL,ielas=0,icmd=0,iinc=1,nasym=0,
    mass[2]={0,0},stiffness=1,buckling=0,rhsi=1,intscheme=0,*ncocon=NULL,
    *nshcon=NULL,mode=-1,noddiam=-1,coriolis=0,iout,
    ifreebody,*itg=NULL,ntg=0,ngraph=1,mt=mi[1]+1,ne0,*integerglob=NULL,     
    icfd=0,*inomat=NULL,*islavact=NULL,*islavnode=NULL,*nslavnode=NULL,
    *islavsurf=NULL,nmethodl,*kon=NULL,*ipkon=NULL,*ielmat=NULL,nzss,
    *mast1=NULL,*irows=NULL,*jqs=NULL,*ipointer=NULL,i,iread,
    *nactdofinv=NULL,*nodorig=NULL,ndesi,iobject,*iponoel=NULL,node,
    *nodedesi=NULL,*ipoface=NULL,*nodface=NULL,*inoel=NULL,*ipoorel=NULL,
    icoordinate=0,iorientation=0,ishapeenergy=0,imass=0,idisplacement=0,
    *istartdesi=NULL,*ialdesi=NULL,*iorel=NULL,*ipoeldi=NULL,*ieldi=NULL,
    *istartelem=NULL,*ialelem=NULL,ieigenfrequency=0,cyclicsymmetry=0,
    nherm,nev,iev,inoelsize,*itmp=NULL,nmd,nevd,*nm=NULL,*ielorien=NULL,
    igreen=0,iglob=0,idesvar=0,inorm=0,irand=0,*nodedesiinv=NULL,
    *nnodes=NULL,iregion=0,*konfa=NULL,*ipkonfa=NULL,nsurfs,
    *iponor=NULL,*iponoelfa=NULL,*inoelfa=NULL,
    ifreemax,nconstraint,
    *iponexp=NULL,*ipretinfo=NULL,nfield,iforce,*iponod2dto3d=NULL,
    *iponk2dto3d=NULL,ishape=0;
      
  double *stn=NULL,*v=NULL,*een=NULL,cam[5],*xstiff=NULL,*stiini=NULL,*tper,
    *f=NULL,*fn=NULL,qa[4],*epn=NULL,*xstateini=NULL,*xdesi=NULL,
    *vini=NULL,*stx=NULL,*enern=NULL,*xbounact=NULL,*xforcact=NULL,
    *xloadact=NULL,*t1act=NULL,*ampli=NULL,*xstaten=NULL,*eei=NULL,
    *enerini=NULL,*cocon=NULL,*shcon=NULL,*qfx=NULL,*df2=NULL,
    *qfn=NULL,*cgr=NULL,*xbodyact=NULL,*springarea=NULL,*emn=NULL,        
    *clearini=NULL,ptime=0.,*emeini=NULL,*doubleglob=NULL,*au=NULL,
    *ad=NULL,*b=NULL,*aub=NULL,*adb=NULL,*pslavsurf=NULL,
    *pmastsurf=NULL,*cdn=NULL,*xstate=NULL,*fnext=NULL,*energyini=NULL,
    *energy=NULL,*ener=NULL,*dxstiff=NULL,*d=NULL,*z=NULL,
    distmin,*df=NULL,*g0=NULL,*dgdx=NULL,sigma=0,*xinterpol=NULL,
    *dgdxglob=NULL,*extnor=NULL,*veold=NULL,*accold=NULL,bet,gam,
    dtime,time,reltime=1.,*weightformgrad=NULL,*fint=NULL,*xnor=NULL,
    *dgdxdy=NULL,*senvector=NULL;

  FILE *f1;
  
#ifdef SGI
  ITG token;
#endif

  irow=*irowp;ener=*enerp;xstate=*xstatep;ipkon=*ipkonp;lakon=*lakonp;
  kon=*konp;ielmat=*ielmatp;ielorien=*ielorienp;objectset=*objectsetp;

  tper=&timepar[1];

  time=*tper;
  dtime=*tper;

  ne0=*ne;

  /* determining the global values to be used as boundary conditions
     for a submodel */

  ITG irefine=0;
  getglobalresults(&jobnamec[396],&integerglob,&doubleglob,nboun,iamboun,xboun,
		   nload,sideload,iamload,&iglob,nforc,iamforc,xforc,
                   ithermal,nk,t1,iamt1,&sigma,&irefine);
  
  /* check which design variables are active */
  
  for(i=0;i<*ntie;i++){
    if(strcmp1(&tieset[i*243+80],"D")==0){
      if(strcmp1(&tieset[i*243],"COORDINATE")==0){
	icoordinate=1;
	if(*nobject>0){
	  if(strcmp1(&objectset[90],"EDG")==0){iregion=1;}
	}
	break;
      }else if(strcmp1(&tieset[i*243],"ORIENTATION")==0){
	if(*norien==0){
	  printf(" *ERROR in sensitivity: the ORIENTATION sensitivity was requested,\n");
	  printf("        yet no orientations were defined.\n");
	  FORTRAN(stop,());
	}
	iorientation=1;
	break;
      }
    }
  }
  
  /* check which targets are active */
  
  for(i=0;i<*nobject;i++){
    if((strcmp1(&objectset[i*324],"DISPLACEMENT")==0)||
       (strcmp1(&objectset[i*324],"X-DISP")==0)||
       (strcmp1(&objectset[i*324],"Y-DISP")==0)||
       (strcmp1(&objectset[i*324],"Z-DISP")==0)){
      idisplacement=1;
    }else if(strcmp1(&objectset[i*324],"EIGENFREQUENCY")==0){
      ieigenfrequency=1;
    }else if(strcmp1(&objectset[i*324],"GREEN")==0){
      ieigenfrequency=1;
      igreen=1;
    }else if(strcmp1(&objectset[i*324],"MASS")==0){
      imass=1;
    }else if(strcmp1(&objectset[i*324],"STRAINENERGY")==0){
      ishapeenergy=1;
    }else if(strcmp1(&objectset[i*324],"STRESS")==0){
      idisplacement=1;
      //      }else if(strcmp1(&objectset[i*324],"THICKNESS")==0){
      //	  ithickness=1;
    }
  }

  /* check number of constraints */

  nconstraint=0;
  for(i=0;i<*nobject;i++){
    if((strcmp1(&objectset[i*324+18],"LE")==0)||
       (strcmp1(&objectset[i*324+18],"GE")==0)){
      nconstraint++;
    }
  }

  /* EIGENFREQUENCY as objective should not be used with any
     other objective in the same step */

  if(((idisplacement==1)||(imass==1)||(ishapeenergy==1))&&
     (ieigenfrequency==1)){
    printf(" *ERROR in sensitivity: the objective EIGENFREQUENCY\n");
    printf("        cannot be used with any other objective within\n");
    printf("        the same step\n");
    exit(0);
  }

  if(ishapeenergy==1){
    NNEW(enerini,double,mi[0]**ne);
    NNEW(emeini,double,6*mi[0]**ne); 
    NNEW(stiini,double,6*mi[0]**ne); 
      
    memcpy(&enerini[0],&ener[0],sizeof(double)*mi[0]**ne);
    memcpy(&emeini[0],&eme[0],sizeof(double)*6*mi[0]**ne);
    memcpy(&stiini[0],&sti[0],sizeof(double)*6*mi[0]**ne);
  }

  if(ieigenfrequency==1){

    /* opening the eigenvalue file and checking for cyclic symmetry */
      
    strcpy(fneig,jobnamec);
    strcat(fneig,".eig");
      
    if((f1=fopen(fneig,"rb"))==NULL){
      printf(" *ERROR in sensitivity: cannot open eigenvalue file for reading");
      printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
      
    if(fread(&cyclicsymmetry,sizeof(ITG),1,f1)!=1){
      printf(" *ERROR in sensitivity reading the cyclic symmetry flag in the eigenvalue file");
      printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
      
    if(fread(&nherm,sizeof(ITG),1,f1)!=1){
      printf(" *ERROR in sensitivity reading the Hermitian flag in the eigenvalue file");
      printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
      
    if(nherm!=1){
      printf(" *ERROR in sensitivity: the eigenvectors in the .eig-file result\n");
      printf("       from a non-Hermitian eigenvalue problem. The \n");
      printf("       sensitivity procedure cannot handle that yet\n\n");
      FORTRAN(stop,());
    }
    //      iperturbsav=iperturb[0];
    if(fread(&iperturb[0],sizeof(ITG),1,f1)!=1){
      printf(" *ERROR in sensitivity reading the perturbation flag in the eigenvalue file");
      printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
    /*      if(iperturbsav!=iperturb[0]){
	    if(iperturb[0]==0){
	    printf(" *ERROR in sensitivity: the .eig-file was created without perturbation info\n");
	    printf("        the present *SENSITIVITY step, however, contains the\n");
	    printf("        perturbation parameter on the *STEP card.\n");
	    }else if(iperturb[0]==1){
	    printf(" *ERROR in sensitivity: the .eig-file was created with perturbation info\n");
	    printf("        the present *SENSITIVITY step, however, does not contain the\n");
	    printf("        perturbation parameter on the *STEP card.\n");
	    }
	    FORTRAN(stop,());
	    }*/

    if(iperturb[0]==1){
      if(fread(vold,sizeof(double),mt**nk,f1)!=mt**nk){
	printf("*ERROR in sensitivity reading the reference displacements in the eigenvalue file...");
	exit(0);
      }
    }
  }

  /* determining the elements belonging to a given node */
  
  NNEW(iponoel,ITG,*nk);
  NNEW(inoel,ITG,2**nkon);
  FORTRAN(elementpernode,(iponoel,inoel,lakon,ipkon,kon,ne));

  if(icoordinate==1){
     
    /* find the 
       - external faces belonging to a given node
       - nodes belonging to a given external surface */

    NNEW(ipoface,ITG,*nk);
    NNEW(nodface,ITG,5*6**ne);
    NNEW(konfa,ITG,8*6**ne);
    NNEW(ipkonfa,ITG,6**ne);
    NNEW(lakonfa,char,8*6**ne);
    FORTRAN(findextsurface,(nodface,ipoface,ne,ipkon,lakon,kon,
			    konfa,ipkonfa,nk,lakonfa,&nsurfs,
			    &ifreemax));
    RENEW(nodface,ITG,5*ifreemax);
    RENEW(konfa,ITG,8*nsurfs);
    RENEW(ipkonfa,ITG,nsurfs);
    RENEW(lakonfa,char,8*nsurfs);

    /* find the external faces belonging to a given node */

    NNEW(iponoelfa,ITG,*nk);
    NNEW(inoelfa,ITG,3*8*6**ne);
    FORTRAN(extfacepernode,(iponoelfa,inoelfa,lakonfa,ipkonfa,konfa,
			    &nsurfs,&inoelsize));
    RENEW(inoelfa,ITG,3*inoelsize);
      
    /* determining the design variables */
      
    NNEW(nodedesi,ITG,*nk);
    NNEW(itmp,ITG,*nk);
    NNEW(nodedesiinv,ITG,*nk);
      
    if(*ne2d!=0){

      NNEW(iponod2dto3d,ITG,3**nk);
      NNEW(iponk2dto3d,ITG,*nk);
	 
      FORTRAN(getdesiinfo2d,(set,istartset,iendset,ialset,nset,
			     mi,nactdof,&ndesi,nodedesi,ntie,tieset,
			     nodedesiinv,lakon,ipkon,kon,iponoelfa,
			     iponod2dto3d,iponor2d,knor2d,iponoel2d,
			     inoel2d,nobject,objectset,iponk2dto3d,ne));
			   			     
      
    }else{
      
      FORTRAN(getdesiinfo3d,(set,istartset,iendset,ialset,nset,
			     mi,nactdof,&ndesi,nodedesi,ntie,tieset,
			     itmp,nmpc,nodempc,ipompc,nodedesiinv,
			     iponoel,inoel,lakon,ipkon,
			     kon,&iregion,ipoface,nodface,nk));  
    }
      
    SFREE(itmp);
    RENEW(nodedesi,ITG,ndesi);

    /* storing the elements to which each design variable belongs
       in field ialdesi */

    NNEW(istartdesi,ITG,ndesi+1);
    NNEW(ialdesi,ITG,*nkon);
    FORTRAN(elemperdesi,(&ndesi,nodedesi,iponoel,inoel,
			 istartdesi,ialdesi,lakon,ipkon,kon,
			 nodedesiinv,&icoordinate,&iregion));
    RENEW(ialdesi,ITG,istartdesi[ndesi]-1);
            
    /* calculating the normal direction for every designvariable */
      
    NNEW(extnor,double,3**nk);
      
    FORTRAN(normalsonsurface_se,(ipkon,kon,lakon,extnor,co,nk,ipoface,
				 nodface,nactdof,mi,nodedesiinv,&iregion,
				 iponoelfa,&ndesi,nodedesi,iponod2dto3d,
				 ikboun,nboun,ne2d)); 
      
    /* if the sensitivity calculation is used in a optimization script
       this script usually contains a loop consisting of:
       1. a call to CalculiX to define the sensitivities
       2. a small modification of the surface geometry in a direction which
       decrease the objective function (only the design variables)
       3. a modification of the internal mesh in order to preserve
       mesh quality
       The latter point can be done by performing a linear elastic
       calculation in which the small modification in 2. is applied
       a *boundary condition and all other nodes (on the external 
       surface but no design variables) are fixed by *equation's
       in a direction normal to the surface. At corners and edges
       there my be more than one normal. The necessary equations are
       calculated in normalsforeq_se.f and stored in jobname.equ */

    NNEW(iponor,ITG,8*nsurfs);
    for(i=0;i<8*nsurfs;i++) iponor[i]=-1;
    NNEW(xnor,double,24*nsurfs);
    NNEW(iponexp,ITG,2**nk);
    NNEW(ipretinfo,ITG,*nk);
      
    FORTRAN(normalsforequ_se,(nk,co,iponoelfa,inoelfa,konfa,ipkonfa,lakonfa,
			      &nsurfs,iponor,xnor,nodedesiinv,jobnamef,
			      iponexp,nmpc,labmpc,ipompc,nodempc,ipretinfo,
			      kon,ipkon,lakon,iponoel,inoel,iponor2d,knor2d,
			      iponod2dto3d,ipoface,nodface));
                  
    SFREE(konfa);SFREE(ipkonfa);SFREE(lakonfa);SFREE(iponor);SFREE(xnor);
    SFREE(iponoelfa);SFREE(inoelfa);SFREE(iponexp);SFREE(ipretinfo);
	  
    /* createinum is called in order to determine the nodes belonging
       to elements; this information is needed in frd_se */
      
    NNEW(inum,ITG,*nk);
    FORTRAN(createinum,(ipkon,inum,kon,lakon,nk,ne,&cflag[0],nelemload,
			nload,nodeboun,nboun,ndirboun,ithermal,co,vold,mi,ielmat,
			ielprop,prop));
	  				   
    /* storing the normal information in the frd-file for the optimizer */
          
    ++*kode;

    inorm=1;
    nfield=3;
    iforce=0;
      
    if(strcmp1(&filab[4],"I")==0){
      
      FORTRAN(map3dto1d2d,(extnor,ipkon,inum,kon,lakon,&nfield,nk,
			   ne,cflag,co,vold,&iforce,mi,ielprop,prop));
    }
      
    frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,
	    istep,
	    &iinc,&mode,&noddiam,description,mi,&ngraph,ne,cs,set,nset,
	    istartset,iendset,ialset,jobnamec,output,
	    extnor,&iobject,objectset,ntrans,inotr,trab,&idesvar,orname,
	    &icoordinate,&inorm,&irand,&ishape); 
    inorm=0;

    /* storing the normal direction for every design variable */

    NNEW(xdesi,double,3*ndesi);
    for(k=0;k<ndesi;k++){
      node=nodedesi[k]-1;
      memcpy(&xdesi[3*k],&extnor[3*node],sizeof(double)*3);
    }
      
    /* calculation of the smallest distance between nodes */
      
    FORTRAN(smalldist,(co,&distmin,lakon,ipkon,kon,ne));

    /* resizing xdesi to a length of distmin */

    for(k=0;k<3*ndesi;k++){
      xdesi[k]*=distmin;
    }

    /* calculation of gaussian random fields for robust optimization */
      
    /*    if(physcon[10]>0){      
	 
      randomfieldmain(kon,ipkon,lakon,ne,nmpc,nactdof,mi,nodedesi,&ndesi,
		      istartdesi,ialdesi,co,physcon,isolver,ntrans,nk,inotr,trab,jobnamec,
		      nboun,cs,mcs,inum,nmethod,kode,filab,nstate_,istep,description,set,
		      nset,iendset,output,istartset,ialset,extnor); 
	    	           
		      }*/

    SFREE(inum);SFREE(extnor);      

  }else if(iorientation==1){
    ndesi=3**norien;
    distmin=1.e-3;

    /* writing the design variables into the dat-file */

    FORTRAN(writedesi,(norien,orname));

    /* storing the elements with a given orientation in
       ipoorel and iorel */

    NNEW(ipoorel,ITG,*norien);
    NNEW(iorel,ITG,2**ne);
    FORTRAN(elemperorien,(ipoorel,iorel,ielorien,ne,mi));

    /* storing the orientation of the design variables
       in nodedesi (per orientation there are three design
       variables - the Euler angles) */

    NNEW(nodedesi,ITG,ndesi);
    for(i=0;i<ndesi;i++){
      nodedesi[i]=i/3+1;
    }

    /* storing the elements corresponding with a given
       design variable in istartdesi and ialdesi */

    NNEW(istartdesi,ITG,ndesi+1);
    NNEW(ialdesi,ITG,3**ne);
    FORTRAN(elemperdesi,(&ndesi,nodedesi,ipoorel,iorel,
			 istartdesi,ialdesi,lakon,ipkon,kon,
			 nodedesiinv,&icoordinate,&iregion));
    SFREE(ipoorel);SFREE(iorel);SFREE(nodedesi);
    RENEW(ialdesi,ITG,istartdesi[ndesi]-1);

  }
      
  /* storing the design variables per element
     (including 0 for the unperturbed state) */
      
  NNEW(ipoeldi,ITG,*ne);
  NNEW(ieldi,ITG,2*(istartdesi[ndesi]-1+*ne));
  NNEW(istartelem,ITG,*ne+1);
  NNEW(ialelem,ITG,istartdesi[ndesi]-1+*ne);
  
  FORTRAN(desiperelem,(&ndesi,istartdesi,ialdesi,ipoeldi,ieldi,ne,
		       istartelem,ialelem));
  
  SFREE(ipoeldi);SFREE(ieldi);
  RENEW(ialelem,ITG,istartelem[*ne]-1);

  /* allocating fields for the actual external loading */

  NNEW(xbounact,double,*nboun);
  for(k=0;k<*nboun;++k){xbounact[k]=xbounold[k];}
  NNEW(xforcact,double,*nforc);
  NNEW(xloadact,double,2**nload);
  NNEW(xbodyact,double,7**nbody);
  /* copying the rotation axis and/or acceleration vector */
  for(k=0;k<7**nbody;k++){xbodyact[k]=xbody[k];}
  if(*ithermal==1){
    NNEW(t1act,double,*nk);
    for(k=0;k<*nk;++k){t1act[k]=t1old[k];}
  }
  
  /* assigning the body forces to the elements */ 

  /* if(*nbody>0){
    ifreebody=*ne+1;
    NNEW(ipobody,ITG,2*ifreebody**nbody);
    for(k=1;k<=*nbody;k++){
      FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			 iendset,ialset,&inewton,nset,&ifreebody,&k));
      RENEW(ipobody,ITG,2*(*ne+ifreebody));
    }
    RENEW(ipobody,ITG,2*(ifreebody-1));
    }*/

  /* allocating a field for the instantaneous amplitude */

  NNEW(ampli,double,*nam);

  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
		    xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,xbodyact,
		    t1old,t1,t1act,iamt1,nk,amta,
		    namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
		    xbounold,xboun,xbounact,iamboun,nboun,
		    nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
		    co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
		    ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,istartset,
		    iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,
		    ipobody,iponoel,inoel,ipkon,kon,ielprop,prop,ielmat,
		    shcon,nshcon,rhcon,nrhcon,cocon,ncocon,ntmat_,lakon));

  /* determining the structure of the df matrix */

  nzss=20000000;
  NNEW(mast1,ITG,nzss);
  NNEW(irows,ITG,1);
  NNEW(jqs,ITG,ndesi+1);
  NNEW(ipointer,ITG,ndesi);

  mastructse(kon,ipkon,lakon,ne,ipompc,nodempc,nmpc,nactdof,jqs,
             &mast1,&irows,ipointer,&nzss,mi,mortar,nodedesi,&ndesi,
             &icoordinate,ielorien,istartdesi,ialdesi);

  SFREE(ipointer);SFREE(mast1);
  RENEW(irows,ITG,nzss);
  
  /* invert nactdof */
  
  NNEW(nactdofinv,ITG,mt**nk);
  NNEW(nodorig,ITG,*nk);
  FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			 ipkon,lakon,kon,ne));
  SFREE(nodorig);

  /* reading the stiffness matrix, mass matrix, eigenfrequencies
     and eigenmodes */

  if(ieigenfrequency==1){

    /* reading the eigenvalues / eigenmodes */
      
    if(!cyclicsymmetry){
	  
      if(fread(&nev,sizeof(ITG),1,f1)!=1){
	printf(" *ERROR in sensitivity reading the number of eigenvalues in the eigenvalue file");
	printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }
	  
      NNEW(d,double,nev);
	  
      if(fread(d,sizeof(double),nev,f1)!=nev){
	printf(" *ERROR in sensitivity reading the eigenvalues in the eigenvalue file");
	printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }
	  
      /*	  for(i=0;i<nev;i++){
		  if(d[i]>0){d[i]=sqrt(d[i]);}else{d[i]=0.;}
		  }*/
	  
      NNEW(ad,double,neq[1]);
      NNEW(adb,double,neq[1]);
      NNEW(au,double,nzsprevstep[2]);
      NNEW(aub,double,nzs[1]);
	  
      if(fread(ad,sizeof(double),neq[1],f1)!=neq[1]){
	printf(" *ERROR in sensitivity reading the diagonal of the stiffness matrix in the eigenvalue file");
	printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }
	  
      if(fread(au,sizeof(double),nzsprevstep[2],f1)!=nzsprevstep[2]){
	printf(" *ERROR in sensitivity reading the off-diagonals of the stiffness matrix in the eigenvalue file");
	printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }
	  
      if(fread(adb,sizeof(double),neq[1],f1)!=neq[1]){
	printf(" *ERROR in sensitivity reading the diagonal of the mass matrix in the eigenvalue file");
	printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }
	  
      if(fread(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
	printf(" *ERROR in sensitivity reading the off-diagonals of the mass matrix in the  eigenvalue file");
	printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }
	  
      NNEW(z,double,neq[1]*nev);
	  
      if(fread(z,sizeof(double),neq[1]*nev,f1)!=neq[1]*nev){
	printf(" *ERROR in sensitivity reading the eigenvectors in the eigenvalue file");
	printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }
    }
    else{
      nev=0;
      do{
	if(fread(&nmd,sizeof(ITG),1,f1)!=1){
	  break;
	}
	if(fread(&nevd,sizeof(ITG),1,f1)!=1){
	  printf(" *ERROR in sensitivity reading the number of eigenvalues for nodal diameter %" ITGFORMAT " in the eigenvalue file",nmd);
	  printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
	  printf("        1) the nonexistence of the .eig file\n");
	  printf("        2) other boundary conditions than in the input deck\n");
	  printf("           which created the .eig file\n\n");
	  exit(0);
	}
	if(nev==0){
	  NNEW(d,double,nevd);
	  NNEW(nm,ITG,nevd);
	}else{
	  RENEW(d,double,nev+nevd);
	  RENEW(nm,ITG,nev+nevd);
	}
	      
	if(fread(&d[nev],sizeof(double),nevd,f1)!=nevd){
	  printf(" *ERROR in sensitivity reading the eigenvalues for nodal diameter %" ITGFORMAT " in the eigenvalue file",nmd);
	  printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
	  printf("        1) the nonexistence of the .eig file\n");
	  printf("        2) other boundary conditions than in the input deck\n");
	  printf("           which created the .eig file\n\n");
	  exit(0);
	}
	      
	for(i=nev;i<nev+nevd;i++){nm[i]=nmd;}
	      
	if(nev==0){

	  /* reading stiffness and mass matrix; */

	  NNEW(ad,double,neq[1]);
	  NNEW(au,double,nzs[1]);
		  
	  if(fread(ad,sizeof(double),neq[1],f1)!=neq[1]){
	    printf(" *ERROR in sensitivity reading the diagonal of the stiffness matrix in the eigenvalue file");
	    printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
	    printf("        1) the nonexistence of the .eig file\n");
	    printf("        2) other boundary conditions than in the input deck\n");
	    printf("           which created the .eig file\n\n");
	    exit(0);
	  }
		  
	  if(fread(au,sizeof(double),nzs[1],f1)!=nzs[1]){
	    printf(" *ERROR in sensitivity reading the off-diagonals of the stiffness matrix in the eigenvalue file");
	    printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
	    printf("        1) the nonexistence of the .eig file\n");
	    printf("        2) other boundary conditions than in the input deck\n");
	    printf("           which created the .eig file\n\n");
	    exit(0);
	  }

	  NNEW(adb,double,neq[1]);
	  NNEW(aub,double,nzs[1]);
		  
	  if(fread(adb,sizeof(double),neq[1],f1)!=neq[1]){
	    printf(" *ERROR in sensitivity reading the diagonal of the mass matrix in the eigenvalue file");
	    printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
	    printf("        1) the nonexistence of the .eig file\n");
	    printf("        2) other boundary conditions than in the input deck\n");
	    printf("           which created the .eig file\n\n");
	    exit(0);
	  }
		  
	  if(fread(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
	    printf(" *ERROR in sensitivity reading the off-diagonals of the mass matrix in the eigenvalue file");
	    printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
	    printf("        1) the nonexistence of the .eig file\n");
	    printf("        2) other boundary conditions than in the input deck\n");
	    printf("           which created the .eig file\n\n");
	    exit(0);
	  }

	  if(igreen!=1){SFREE(ad);SFREE(au);SFREE(adb);SFREE(aub);}
	}
	      
	if(nev==0){
	  NNEW(z,double,neq[1]*nevd);
	}else{
	  RENEW(z,double,(long long)neq[1]*(nev+nevd));
	}
	      
	if(fread(&z[(long long)neq[1]*nev],sizeof(double),neq[1]*nevd,f1)!=neq[1]*nevd){
	  printf(" *ERROR in sensitivity reading the eigenvectors for nodal diameter %" ITGFORMAT " in the eigenvalue file",nmd);
	  printf(" *INFO  in sensitivity: if there are problems reading the .eig file this may be due to:\n");
	  printf("        1) the nonexistence of the .eig file\n");
	  printf("        2) other boundary conditions than in the input deck\n");
	  printf("           which created the .eig file\n\n");
	  exit(0);
	}
	nev+=nevd;
      }while(1);
    }
    fclose(f1);
  }else{
    nev=1;
    if((idisplacement==1)||((ishapeenergy==1)&&(iperturb[1]==1))){
	
      /* reading the stiffness matrix from previous step for sensitivity analysis */
      /* matrix stored in <jobname>.stm file */

      /* nzs,irow,jq and icol are stored too, since the static analysis
	 can involve contact, whereas in the sensitivity analysis contact is not
	 taken into account while determining the structure of the stiffness
	 matrix (in mastruct.c)
      */

      /* for mass and strain energy the stiffness matrix is not needed */
	
      strcpy(stiffmatrix,jobnamec);
      strcat(stiffmatrix,".stm");
	
      if((f1=fopen(stiffmatrix,"rb"))==NULL){
	printf("*ERROR in sensitivity: cannot open stiffness-matrix file for reading");
	exit(0);
      }
	
      if(fread(&nasym,sizeof(ITG),1,f1)!=1){
	printf("*ERROR in sensitivity reading the symmetry flag of the stiffness matrix file...");
	exit(0);
      }
	
      if(fread(nzs,sizeof(ITG),3,f1)!=3){
	printf("*ERROR in sensitivity reading the number of subdiagonal nonzeros in the stiffness matrix file...");
	exit(0);
      }
      RENEW(irow,ITG,nzs[2]);

      if(fread(irow,sizeof(ITG),nzs[2],f1)!=nzs[2]){
	printf("*ERROR in sensitivity reading irow in the stiffness matrix file...");
	exit(0);
      }

      if(fread(jq,sizeof(ITG),neq[1]+1,f1)!=neq[1]+1){
	printf("*ERROR in sensitivity reading jq in the stiffness matrix file...");
	exit(0);
      }

      if(fread(icol,sizeof(ITG),neq[1],f1)!=neq[1]){
	printf("*ERROR in sensitivity reading icol in the stiffness matrix file...");
	exit(0);
      }
	
      NNEW(ad,double,neq[1]);
      NNEW(au,double,(nasym+1)*nzs[2]);

      if(fread(ad,sizeof(double),neq[1],f1)!=neq[1]){
	printf("*ERROR in sensitivity reading the diagonal of the stiffness matrix in the .stm-file");
	exit(0);
      }
	
      if(fread(au,sizeof(double),(nasym+1)*nzs[2],f1)!=(nasym+1)*nzs[2]){
	printf("*ERROR in sensitivity reading the off-diagonals of the stiffness matrix in the .stm-file");
	exit(0);
      }  
	
      fclose(f1);
    }
  }

  /* loop over all eigenvalues, or, if it is not an eigenvalue sensitivity study,
     loop over just one value */

  for(iev=0;iev<nev;iev++){
      
    NNEW(g0,double,*nobject);
    NNEW(dgdx,double,ndesi**nobject);
    if(icoordinate==1){NNEW(dgdxglob,double,2**nk**nobject);}

    /* Reading the "raw" sensititities */

    iread=0;
    if(*nobject>0){
      if(strcmp1(&objectset[80],"R")==0){
	FORTRAN(readsen,(g0,dgdx,&ndesi,nobject,nodedesi,jobnamef));
	iread=1;
      }
    }

    if(iread==0){

      /* for cyclic symmetry calculations only the odd modes are calculated
         (modes occur in phase-shifted pairs) */

      if(cyclicsymmetry){
	if((iev/2)*2!=iev){
	  SFREE(g0);SFREE(dgdx);
	  if(icoordinate==1){SFREE(dgdxglob);}
	  continue;
	}
	mode=iev;
	noddiam=nm[iev];
      }

      /* determining the internal forces and the stiffness coefficients */
      
      NNEW(f,double,*neq); /* FAKE */
      
      /* needed for nonlinear strain energy */

      if((iperturb[1]==1)&&(ishapeenergy==1)){
	NNEW(fint,double,*neq);
      }
      
      /* allocating a field for the stiffness matrix 
	 (calculated in results_se and needed in mafillsmse */
      
      NNEW(xstiff,double,(long long)27*mi[0]**ne);
      if(iorientation==1) NNEW(dxstiff,double,(long long)3*27*mi[0]**ne);

      NNEW(fn,double,mt**nk);  /* FAKE */
      if(!cyclicsymmetry){
	NNEW(df,double,nzss);
      }else{
	NNEW(df,double,2*nzss);
      }
      NNEW(stx,double,6*mi[0]**ne);
      
      iout=-1;
      NNEW(v,double,mt**nk);
      if(iperturb[1]==1) memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
      
      results_se(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		 elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		 ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		 prestr,iprestr,filab,eme,emn,een,iperturb,
		 f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		 ndirboun,xbounact,nboun,ipompc,
		 nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,veold,accold,
		 &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		 xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		 &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
		 emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
		 iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
		 fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		 &reltime,&ne0,xforc,nforc,thicke,shcon,nshcon,
		 sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		 mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
		 islavsurf,ielprop,prop,energyini,energy,df,&distmin,
		 &ndesi,nodedesi,sti,nkon,jqs,irows,nactdofinv,
		 &icoordinate,dxstiff,istartdesi,ialdesi,xdesi,
		 &ieigenfrequency,fint,&ishapeenergy,typeboun);
	  
      iout=1;SFREE(v);
      
      nmethodl=*nmethod;
      
      /* v contains the values dK/dx has to be multiplied with */
      
      if(ieigenfrequency==0){
	NNEW(v,double,mt**nk);
	memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
      }else{
	if(!cyclicsymmetry){
	  NNEW(v,double,mt**nk);
	  FORTRAN(resultsnoddir,(nk,v,nactdof,&z[iev*neq[1]],ipompc,
				 nodempc,coefmpc,nmpc,mi));
	}else{
	  NNEW(v,double,2*mt**nk);
	  FORTRAN(resultsnoddir,(nk,v,nactdof,&z[iev*neq[1]],ipompc,
				 nodempc,coefmpc,nmpc,mi));
	  FORTRAN(resultsnoddir,(nk,&v[mt**nk],nactdof,&z[iev*neq[1]+neq[1]/2],ipompc,
				 nodempc,coefmpc,nmpc,mi));
	}
	  
	ptime=d[iev];
	if(ptime>0){ptime=sqrt(ptime)/6.283185308;}else{ptime=0.;}

	/* for an eigenfrequency objective (K-eigenvalue*M) is
	   taken instead of K (not for ORIENTATION as design 
	   variable since the mass does not change with orientation)*/

	if(icoordinate==1){
	  sigma=d[iev];
	  mass[0]=1;
	}
      }

      /* determining the system matrix and the external forces */
      
      mafillsmmain_se(co,nk,kon,ipkon,lakon,ne,nodeboun,
		      ndirboun,xbounact,nboun,
		      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		      nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
		      nbody,cgr,nactdof,neq,&nmethodl,
		      ikmpc,ilmpc,ikboun,ilboun,
		      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		      ielorien,norien,orab,ntmat_,
		      t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
		      stx,iexpl,plicon,nplicon,plkcon,nplkcon,
		      xstiff,npmat_,&dtime,matname,mi,
		      ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,physcon,
		      shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,&coriolis,
		      ibody,xloadold,&reltime,veold,springarea,nstate_,
		      xstateini,xstate,thicke,integerglob,doubleglob,
		      tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
		      pmastsurf,mortar,clearini,ielprop,prop,&ne0,fnext,
		      &distmin,&ndesi,nodedesi,df,&nzss,jqs,irows,
		      &icoordinate,dxstiff,xdesi,istartelem,ialelem,v,&sigma,
		      &cyclicsymmetry,labmpc,ics,cs,mcs,&ieigenfrequency);
      
      if(iorientation==1) SFREE(dxstiff);
      
      /* determining the values and the derivatives of the objective functions */
            
      iout=-1; 
      objectivemain_se(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		       ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		       prestr,iprestr,filab,eme,emn,een,iperturb,
		       f,fn,nactdof,&iout,qa,vold,nodeboun,
		       ndirboun,xbounact,nboun,ipompc,
		       nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,veold,accold,
		       &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		       xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		       &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
		       emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
		       iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
		       fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		       &reltime,&ne0,xforc,nforc,thicke,shcon,nshcon,
		       sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		       mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
		       islavsurf,ielprop,prop,energyini,energy,&distmin,
		       &ndesi,nodedesi,nobject,objectset,g0,dgdx,sti,
		       df,nactdofinv,jqs,irows,&idisplacement,nzs,jobnamec,
		       isolver,icol,irow,jq,kode,cs,output,istartdesi,ialdesi,
		       xdesi,orname,&icoordinate,&iev,d,z,au,ad,aub,adb,&cyclicsymmetry,
		       &nzss,&nev,&ishapeenergy,fint,nlabel,&igreen,&nasym,
		       iponoel,inoel,nodedesiinv,dgdxglob,df2,dgdxdy,nkon,iponod2dto3d,
		       iponk2dto3d,ics,mcs,mpcend,&noddiam,ipobody,ibody,xbody,nbody); 
      iout=1;

      SFREE(v);SFREE(f);SFREE(xstiff);SFREE(fn);SFREE(df);SFREE(stx);
      if((iperturb[1]==1)&&(ishapeenergy==1)){SFREE(fint);}

    }
      
    if(icoordinate==1){

      /* Storing the "raw" sensititities */
      if(*nobject>0){
	if(strcmp1(&objectset[80],"W")==0){
	  FORTRAN(writesen,(g0,dgdx,&ndesi,nobject,nodedesi,jobnamef));
	}
      }

      /* Elimination of mesh-dependency and filling of dgdxglob  */
  
      NNEW(weightformgrad,double,ndesi);
          
      FORTRAN(distributesens,(istartdesi,ialdesi,ipkon,lakon,ipoface,
			      &ndesi,nodedesi,nodface,kon,co,dgdx,nobject,
			      weightformgrad,nodedesiinv,&iregion,objectset,
			      dgdxglob,nk,physcon));
		  
      SFREE(weightformgrad);

      /* for quadratic elements: interpolation of midnode-sensitivity 
	 to the corner nodes */
	  
      NNEW(xinterpol,double,*nk**nobject);
      NNEW(nnodes,ITG,*nk);
	  
      FORTRAN(quadraticsens,(ipkon,lakon,kon,nobject,dgdxglob,
			     xinterpol,nnodes,ne,nk,nodedesiinv,objectset));
		            
      SFREE(nnodes);SFREE(xinterpol);

      /* Filtering of sensitivities */
  
      filtermain(co,dgdxglob,nobject,nk,nodedesi,&ndesi,objectset,
		 xdesi,&distmin);

      /* scaling the designnodes being in the transition between 
	 the designspace and the non-designspace */
	 
      transitionmain(co,dgdxglob,nobject,nk,nodedesi,&ndesi,objectset,
		     ipkon,kon,lakon,ipoface,nodface,nodedesiinv);
	  
      /* createinum is called in order to determine the nodes belonging
	 to elements; this information is needed in frd_se */

      NNEW(inum,ITG,*nk);
      FORTRAN(createinum,(ipkon,inum,kon,lakon,nk,ne,&cflag[0],nelemload,
			  nload,nodeboun,nboun,ndirboun,ithermal,co,vold,mi,ielmat,
			  ielprop,prop));

      /* calculate projected gradient if:
	 - any constraint functions are defined AND
	 - at least one constraint function is active */

      projectgradmain(nobject,&objectset,&dgdxglob,g0,&ndesi,nodedesi,
		      nk,isolver,co,xdesi,&distmin,&nconstraint);	     

      /* scaling the designnodes being in the transition between 
	 the designspace and the non-designspace */
	 
      transitionmain(co,dgdxglob,nobject,nk,nodedesi,&ndesi,objectset,
		     ipkon,kon,lakon,ipoface,nodface,nodedesiinv);

      /* for 2D models extrapolate the results of the midnodes 
	 to the 2 symmetry planes */

      if(*ne2d!=0){
	FORTRAN(extrapol2dto3d,(dgdxglob,iponod2dto3d,&ndesi,
				nodedesi,nobject,nk,xinterpol,nnodes,
				ipkon,lakon,kon,ne,iponoel,inoel));
      }
	    
      //++*kode;

      NNEW(senvector,double,3**nk);
	  
      for(iobject=0;iobject<*nobject;iobject++){

	/* storing the sensitivities in the frd-file for visualization 
	   and for optimization */

	++*kode;
	nfield=2;
	iforce=0;

	if(strcmp1(&filab[4],"I")==0){
           
	  FORTRAN(map3dto1d2d,(&dgdxglob[2**nk*iobject],ipkon,inum,kon,lakon,&nfield,nk,
			       ne,cflag,co,vold,&iforce,mi,ielprop,prop));
      
	}
          
	frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,
		istep,
		&iinc,&mode,&noddiam,description,mi,&ngraph,ne,cs,set,nset,
		istartset,iendset,ialset,jobnamec,output,
		dgdxglob,&iobject,objectset,ntrans,inotr,trab,&idesvar,orname,
		&icoordinate,&inorm,&irand,&ishape);

	/*for(k=0;k<ndesi;k++){ 
	  node=nodedesi[k]-1;	*/	 
	/*memcpy(&senvector[3*node],&xdesi[3*k],sizeof(double)*3);*/
	/*for(i=0;i<3;i++){
	  senvector[3*node+i]=xdesi[3*k+i];
	  senvector[3*node+i]=senvector[3*node+i]/distmin
	  *dgdxglob[2*node+1+2*iobject**nk];
	  }
	  }
               
	  ishape=1;
	  frd_sen(co,nk,stn,inum,nmethod,kode,filab,&ptime,nstate_,istep,	      
	  &iinc,&mode,&noddiam,description,mi,&ngraph,ne,cs,set,nset,
	  istartset,iendset,ialset,jobnamec,output,
	  senvector,&iobject,objectset,ntrans,inotr,trab,&idesvar,orname,
	  &icoordinate,&inorm,&irand,&ishape);
	  ishape=0;*/

	/* writing the objectives in the dat-file for the optimizer */
	  	      
	FORTRAN(writeobj,(objectset,&iobject,g0));
      }  
	  
      SFREE(inum);
	  
    }
      
    SFREE(g0);SFREE(dgdx);
    if(icoordinate==1){SFREE(dgdxglob);SFREE(senvector);}
            
  } // end loop over nev

  SFREE(iponoel);SFREE(inoel);SFREE(nodedesiinv);
  
  if(*ne2d!=0){
    SFREE(iponod2dto3d);SFREE(iponk2dto3d);
  }
  
  if(ieigenfrequency==1){
    if(!cyclicsymmetry){SFREE(d);SFREE(ad);SFREE(adb);SFREE(au);
      SFREE(aub);SFREE(z);
    }else if(igreen!=1){
      SFREE(d);SFREE(z);SFREE(nm);
    }else{
      SFREE(d);SFREE(z);SFREE(nm);SFREE(ad);SFREE(au);SFREE(adb);SFREE(aub);}
  }else if(idisplacement==1){
    SFREE(ad);SFREE(au);
  }
  
  SFREE(xbounact);SFREE(xforcact);SFREE(xloadact);SFREE(t1act);SFREE(ampli);
  SFREE(xbodyact);SFREE(nactdofinv);
  // if(*nbody>0) SFREE(ipobody);

  if(ishapeenergy==1){
    SFREE(enerini);SFREE(emeini);SFREE(stiini);
  }

  SFREE(istartdesi);SFREE(ialdesi);SFREE(istartelem);SFREE(ialelem);

  if(icoordinate==1){
    SFREE(nodedesi);SFREE(xdesi);SFREE(ipoface);SFREE(nodface);
  }

  SFREE(irows);SFREE(jqs);

  *irowp=irow;*enerp=ener;*xstatep=xstate;*ipkonp=ipkon;*lakonp=lakon;
  *konp=kon;*ielmatp=ielmat;*ielorienp=ielorien;*objectsetp=objectset;

  (*ttime)+=(*tper);
 
  return;
}
