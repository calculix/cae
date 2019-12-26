/*     CALCULIX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2019 Guido Dhondt                     */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation; either version 2 of    */
/*     the License, or (at your option) any later version.               */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <pthread.h>

#define Linux 1
#define IRIX 2
#define IRIX64 3
#define HP 4

#if ARCH == Linux
#define FORTRAN(A,B) A##_  B
#elif ARCH == IRIX || ARCH == IRIX64
#define FORTRAN(A,B) A##_##B
#elif ARCH == HP
#define FORTRAN(A,B) A##B
#endif

#if ARCH == Linux
#define CEE(A,B) A##_  B
#elif ARCH == IRIX || ARCH == IRIX64
#define CEE(A,B) A##_##B
#elif ARCH == HP
#define CEE(A,B) A##B
#endif

/* setting arrays to constant values */

/* serial */
#define DMEMSET(a,b,c,d) for(im=b;im<c;im++)a[im]=d
/* parallel */
#define DOUMEMSET(a,b,c,d) setpardou(&a[b],d,c-b,num_cpus)
#define ITGMEMSET(a,b,c,d) setparitg(&a[b],d,c-b,num_cpus)

/* memory allocation, reallocation, freeing */

/* allocating memory for double reals and initializing it to zero (parallell) */
#define DNEW(a,b,c) {a=(b *)u_malloc((c)*sizeof(b),__FILE__,__LINE__,#a); \
        DOUMEMSET(a,0,c,0.);}
/* allocating memory for ITG and initializing it to zero (parallell) */
#define INEW(a,b,c) {a=(b *)u_malloc((c)*sizeof(b),__FILE__,__LINE__,#a); \
        ITGMEMSET(a,0,c,0);}
/* allocating memory without initialization */
#define MNEW(a,b,c) a=(b *)u_malloc((c)*sizeof(b),__FILE__,__LINE__,#a)
/* allocating memory and initializing it to zero (serial) */
#define NNEW(a,b,c) a=(b *)u_calloc((c),sizeof(b),__FILE__,__LINE__,#a)
/* reallocating memory; no initialization */
#define RENEW(a,b,c) a=(b *)u_realloc((b *)(a),(c)*sizeof(b),__FILE__,__LINE__,#a)
/* freeing memory */
#define SFREE(a) u_free(a,__FILE__,__LINE__,#a)

#ifdef LONGLONG
#define ITG long long
#define ITGFORMAT "lld"
#else
#define ITG int
#define ITGFORMAT "d"
#endif

void FORTRAN(actideacti,(char *set,ITG *nset,ITG *istartset,ITG *iendset,
                         ITG *ialset,char *objectset,ITG *ipkon,ITG *ibject,
                         ITG *ne));

void FORTRAN(actideactistr,(char *set,ITG *nset,ITG *istartset,ITG *iendset,
             ITG *ialset,char *objectset,ITG *ipkon,ITG *iobject,ITG *ne,
             ITG *neinset,ITG *iponoel,ITG *inoel,ITG *nepar,
             ITG *nkinsetinv,ITG *nk));

void FORTRAN(adaptfields,(ITG *ipkonf,ITG *ipnei,ITG *ielmatf,ITG *ielorienf,
                          ITG *neiel,ITG *neifa,ITG *neij,ITG *ipkonfcp,
                          ITG *ipneicp,ITG *ielmatfcp,ITG *ielorienfcp,
                          ITG *neielcp,ITG *neifacp,ITG *neijcp,ITG *mmm,
                          ITG *nnn,ITG *nactdoh,ITG *nactdohinv,ITG *mi,
                          ITG *nef,ITG *nface,ITG *ielfa,ITG *ielfacp,
                          ITG *norien,char *lakonf,char *lakonfcp));

void FORTRAN(addimdnodecload,(ITG *nodeforc,ITG *i,ITG *imdnode,
             ITG *nmdnode,double *xforc,ITG *ikmpc,ITG *ilmpc,
             ITG *ipompc,ITG *nodempc,ITG *nmpc,ITG *imddof,ITG *nmddof,
             ITG *nactdof,ITG *mi,ITG *imdmpc,ITG *nmdmpc,ITG *imdboun,
             ITG *nmdboun,ITG *ikboun,ITG *nboun,ITG *ilboun,ITG *ithermal));

void FORTRAN(addimdnodedload,(ITG *nelemload,char *sideload,ITG *ipkon,
             ITG *kon,char *lakon,ITG *i,ITG *imdnode,ITG *nmdnode,
             ITG *ikmpc,ITG *ilmpc,
             ITG *ipompc,ITG *nodempc,ITG *nmpc,ITG *imddof,ITG *nmddof,
             ITG *nactdof,ITG *mi,ITG *imdmpc,ITG *nmdmpc,ITG *imdboun,
             ITG *nmdboun,ITG *ikboun,ITG *nboun,ITG *ilboun,ITG *ithermal));

void FORTRAN(addizdofcload,(ITG *nodeforc,ITG *ndirforc,ITG *nactdof,
             ITG *mi,ITG *izdof,ITG *nzdof,ITG *i,ITG *iznode,ITG *nznode,
             ITG *nk,ITG *imdnode,ITG *nmdnode,double *xforc,
             ITG *ntrans,ITG *inotr));

void FORTRAN(addizdofdload,(ITG *nelemload,char *sideload,ITG *ipkon,
             ITG *kon,char *lakon,ITG *nactdof,ITG *izdof,ITG *nzdof,
             ITG *mi,ITG *i,ITG *iznode,ITG *nznode,ITG *nk,
             ITG *imdnode,ITG *nmdnode));

void FORTRAN(adjustcontactnodes,(char *tieset,ITG *ntie,ITG *itietri,double *cg,
             double *straight,double *co,double *vold,double *xo,double *yo,
             double *zo,double *x,double *y,double *z,ITG *nx,ITG *ny,
             ITG *nz,ITG *istep,ITG *iinc,ITG *iit,ITG *mi,ITG *imastop,
             ITG *nslavnode,ITG *islavnode,char *set,ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,double *tietol,double *clearini,
             double *clearslavnode,ITG *itiefac,ITG *ipkon,ITG *kon,
             char *lakon,ITG *islavsurf));

void FORTRAN(allocation,(ITG *nload_,ITG *nforc_,ITG *nboun_,
             ITG *nk_,ITG *ne_,ITG *nmpc_,ITG *nset_,ITG *nalset_,
             ITG *nmat_,ITG *ntmat_,ITG *npmat_,ITG *norien_,ITG *nam_,
             ITG *nprint_,ITG *mi,ITG *ntrans_,
             char *set,ITG *meminset,
             ITG *rmeminset,ITG *ncs_,ITG *namtot_,ITG *ncmat_,
             ITG *memmpc_,ITG *ne1d,ITG *ne2d,ITG *nflow,
             char *jobnamec,ITG *irstrt,ITG *ithermal,ITG *nener,
             ITG *nstate_,ITG *istep,char *inpc,
             ITG *ipoinp,ITG *inp,ITG *ntie_,ITG *nbody_,
             ITG *nprop_,ITG *ipoinpc,ITG *nevdamp,ITG *npt_,
             ITG *nslavsm,ITG *nkon_,ITG *mcs,ITG *mortar,ITG *ifacecount,
             ITG *nintpoint,ITG *infree,ITG *nheading_,ITG *nobject_,
	     ITG *iuel,ITG *iprestr,ITG *nstam,ITG *ndamp,ITG *nef,
	     ITG *nbounold,ITG *nforcold,ITG *nloadold,ITG *nbodyold));

void FORTRAN(allocont,(ITG *ncont,ITG *ntie,char *tieset,ITG *nset,
             char *set,ITG *istartset,ITG *iendset,ITG *ialset,
             char *lakon,ITG *ncone,double *tietol,ITG *ismallsliding,
             char *kind1,char *kind2,ITG *mortar,ITG *istep));

void FORTRAN(applyboun,(ITG *ifaext,ITG *nfaext,ITG *ielfa,ITG *ikboun,
             ITG *ilboun,ITG *nboun,char *typeboun,ITG *nelemload,
             ITG *nload,char *sideload,ITG *isolidsurf,ITG *nsolidsurf,
             ITG *ifabou,ITG *nfabou,ITG *nface,ITG *nodeboun,ITG *ndirboun,
             ITG *ikmpc,ITG *ilmpc,char *labmpc,ITG *nmpc,ITG *nactdohinv,
             ITG *compressible,ITG *iatleastonepressurebc,ITG *ipkonf,
             ITG *kon,ITG *konf,ITG *inlet));

void FORTRAN(applybounfem,(int *nodeboun,int *ndirboun,int *nboun,
       double *xbounact,int *ithermal,int *nk,int *iponoel,int *inoel,
       double *vold,double *vcontu,double *t1act,int *isolidsurf,
       int *nsolidsurf,double *xsolidsurf,int *nfreestream,int *ifreestream,
       int *turbulent,double *voldaux,double *shcon,int *nshcon,
       double *rhcon,int *nrhcon,int *ielmat,int *ntmat_,double *physcon,
       double *v,int *compressible,int *ismooth,int *nmpc,int *nodempc,
       int *ipompc,double *coefmpc,int *inomat,int *mi,int *ikboun,
       int *ilboun,int *ilmpc,char *labmpc));

void FORTRAN(applybounp,(int *nodeboun,int *ndirboun,int *nboun,
       double *xbounact,int *ithermal,int *nk,int *iponoel,int *inoel,
       double *vold,double *vcontu,double *t1act,int *isolidsurf,
       int *nsolidsurf,double *xsolidsurf,int *nfreestream,int *ifreestream,
       int *turbulent,double *voldaux,double *shcon,int *nshcon,
       double *rhcon,int *nrhcon,int *ielmat,int *ntmat_,double *physcon,
       double *v,int *ipompc,int *nodempc,double *coefmpc,int *nmpc,
       int *inomat,int *mi));

void FORTRAN(applympc,(ITG *nface,ITG *ielfa,ITG *is,ITG *ie,ITG *ifabou,
                       ITG *ipompc,double *vfa,double *coefmpc,ITG *nodempc,
                       ITG *ipnei,ITG *neifa,char *labmpc,double *xbounact,
                       ITG *nactdoh,ITG *ifaext,ITG *nfaext));

void FORTRAN(applympc_dpel,(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
             double *vfa,ITG *ifabou,double *xbounact,ITG *nef,
             double *gradpcel,double *gradpcfa,ITG *neifa,double *rf,
             double *area,double *volume,double *xle,double *xxi,
             ITG *icyclic,double *xxn,ITG *ipnei,ITG *ifatie,
             double *coefmpc,ITG *nmpc,char *labmpc,ITG *ipompc,
             ITG *nodempc,ITG *ifaext,ITG *nfaext,ITG *nactdoh,
             ITG *iflag,double *xxj,double *xlet));

void FORTRAN(applympc_hfa,(ITG *nface,ITG *ielfa,ITG *is,ITG *ie,ITG *ifabou,
                       ITG *ipompc,double *hfa,double *coefmpc,ITG *nodempc,
                       ITG *ipnei,ITG *neifa,char *labmpc,double *xbounact,
                       ITG *nactdoh));

void arpack(double *co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
             ITG *ne,
             ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
             ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
             ITG *nmpc,
             ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
             ITG *nelemload,char *sideload,double *xload,
             ITG *nload,
             ITG *nactdof,
             ITG *icol,ITG *jq,ITG **irowp,ITG *neq,ITG *nzl,
             ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
             ITG *ilboun,
             double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
             double *shcon,ITG *nshcon,double *cocon,ITG *ncocon,
             double *alcon,ITG *nalcon,double *alzero,ITG **ielmatp,
             ITG **ielorienp,ITG *norien,double *orab,ITG *ntmat_,
             double *t0,double *t1,double *t1old,
             ITG *ithermal,double *prestr,ITG *iprestr,
             double *vold,ITG *iperturb,double *sti,ITG *nzs,
             ITG *kode,ITG *mei,double *fei,char *filab,
             ITG *iexpl,double *plicon,ITG *nplicon,double *plkcon,
             ITG *nplkcon,
             double **xstatep,ITG *npmat_,char *matname,ITG *mi,
             ITG *ncmat_,ITG *nstate_,double **enerp,char *jobnamec,
             char *output,char *set,ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,ITG *nener,ITG *isolver,double *trab,
             ITG *inotr,ITG *ntrans,double *ttime,double *fmpc,
             char *cbody,ITG *ibody,double *xbody,ITG *nbody,double *thicke,
             ITG *nslavs,double *tietol,ITG *nkon,ITG *mpcinfo,ITG *ntie,
             ITG *istep,ITG *mcs,ITG *ics,char *tieset,
             double *cs,ITG *nintpoint,ITG *mortar,ITG *ifacecount,
             ITG **islavsurfp,double **pslavsurfp,double **clearinip,
             ITG *nmat,char *typeboun,ITG *ielprop,double *prop,
             char *orname);

void arpackbu(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
             ITG *ne,
             ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
             ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
             ITG *nmpc,
             ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
             ITG *nelemload,char *sideload,double *xload,
             ITG *nload,
             ITG *nactdof,
             ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
             ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
             ITG *ilboun,
             double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
             double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
             ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
             double *t0,double *t1,double *t1old,
             ITG *ithermal,double *prestr,ITG *iprestr,
             double *vold,ITG *iperturb,double *sti,ITG *nzs,
             ITG *kode,ITG *mei,double *fei,
             char *filab,double *eme,
             ITG *iexpl,double *plicon,ITG *nplicon,double *plkcon,
             ITG *nplkcon,
             double *xstate,ITG *npmat_,char *matname,ITG *mi,
             ITG *ncmat_,ITG *nstate_,double *ener,char *output,
             char *set,ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,ITG *nener,ITG *isolver,double *trab,
             ITG *inotr,ITG *ntrans,double *ttime,double *fmpc,
             char *cbody,ITG *ibody,double *xbody,ITG *nbody,
             double *thicke,char *jobnamec,ITG *nmat,ITG *ielprop,
             double *prop,char *orname,char *typeboun);

void arpackcs(double *co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
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
             ITG *kode,ITG *mei,double *fei,
             char *filab,
             ITG *iexpl,double *plicon,ITG *nplicon,double *plkcon,
             ITG *nplkcon,
             double **xstatep,ITG *npmat_,char *matname,ITG *mi,
             ITG *ics,double *cs,ITG *mpcend,ITG *ncmat_,ITG *nstate_,
             ITG *mcs,ITG *nkon,char *jobnamec,
             char *output,char *set,ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,ITG *nener,ITG *isolver,double *trab,
             ITG *inotr,ITG *ntrans,double *ttime,double *fmpc,
             char *cbody,ITG *ibody,double *xbody,ITG *nbody,
             ITG *nevtot,double *thicke,ITG *nslavs,double *tietol,
             ITG *mpcinfo,ITG *ntie,ITG *istep,
             char *tieset,ITG *nintpoint,ITG *mortar,ITG *ifacecount,
             ITG **islavsurfp,double **pslavsurfp,double **clearinip,
             ITG *nmat,char *typeboun,ITG *ielprop,double *prop,
             char *orname);

void FORTRAN(assigndomtonodes,(ITG *ne,char *lakon,ITG *ipkon,ITG *kon,
             ITG *ielmat,ITG *inomat,double *elcon,ITG *ncmat_,ITG *ntmat_,
             ITG *mi,ITG *ne2));

void FORTRAN(autocovmatrix,(double *co,double *ad,double *au,ITG *jqs,
             ITG *irows,ITG *ndesi,ITG *nodedesi,double *physcon));

void FORTRAN(basis,(double *x,double *y,double *z,double *xo,double *yo,
                    double *zo,ITG *nx,ITG *ny,ITG *nz,double *planfa,
                    ITG *ifatet,ITG *nktet,ITG *netet,double *field,
                    ITG *nfield,double *cotet,ITG *kontyp,ITG *ipkon,
                    ITG *kon,ITG *iparent,double *xp,double *yp,double *zp,
                    double *value,double *ratio,ITG *iselect,ITG *nselect,
                    ITG *istartset,ITG *iendset,ITG *ialset,ITG *imastset,
                    ITG *ielemnr,ITG *nterms,ITG *konl));

void biosav(ITG *ipkon,ITG *kon,char *lakon,ITG *ne,double *co,
            double *qfx,double *h0,ITG *mi,ITG *inomat,ITG *nk);

void FORTRAN(biotsavart,(ITG *ipkon,ITG *kon,char *lakon,ITG *ne,double *co,
                         double *qfx,double *h0,ITG *mi,ITG *nka,ITG *nkb));

void *biotsavartmt(ITG *i);

void FORTRAN(blockanalysis,(char *set,ITG *nset,ITG *istartset,ITG *iendset,
             ITG *ialset,ITG *nblk,ITG *ipkon,ITG *kon,ITG *ielfa,
             ITG *nodface,ITG *neiel,ITG *neij,ITG *neifa,ITG *ipoface,
             ITG *ipnei,ITG *konf,ITG *istartblk,ITG *iendblk,ITG *nactdoh,
             ITG *nblket,ITG *nblkze,ITG *nef,ITG *ielblk,ITG *nk,
             ITG *nactdohinv));

void FORTRAN(bodyforce,(char *cbody,ITG *ibody,ITG *ipobody,ITG *nbody,
             char *set,ITG *istartset,ITG *iendset,ITG *ialset,
             ITG *inewton,ITG *nset,ITG *ifreebody,ITG *k));

void FORTRAN(calcbody,(ITG *nef,double *body,ITG *ipobody,ITG *ibody,
                       double *xbody,double *coel,double *vel,char *lakon,
                       ITG *nactdohinv));

void FORTRAN(calcenergy,(ITG *ipkon,char *lakon,ITG *kon,double *co,
                         double *ener,ITG *mi,ITG *ne,double *thicke,
                         ITG *ielmat,
                         double *energy,ITG *ielprop,double *prop,ITG *nea,
                         ITG *neb));

void *calcenergymt(ITG *i);

void FORTRAN(calcflux,(double *area,double *vfa,double *xxna,ITG *ipnei,
                        ITG *nef,ITG *neifa,double *flux,double *xxj,
                        double *gradpfa,double *xlet,double *xle,double *vel,
                        double *advfa,ITG *ielfa,ITG *neiel,ITG *ifabou,
                        double *hfa,ITG *nefa,ITG *nefb));

void *calcflux1mt(ITG *i);

void calcfluxmain(double *area,double *vfa,double *xxna,ITG *ipnei,
                  ITG *nef,ITG *neifa,double *flux,double *xxj,
                  double *gradpfa,double *xlet,double *xle,double *vel,
                  double *advfa,ITG *ielfa,ITG *neiel,ITG *ifabou,
                  double *hfa,ITG *num_cpus);

void FORTRAN(calcinitialflux,(double *area,double *vfa,double *xxna,
                              ITG *ipnei,ITG *nef,ITG *neifa,char *lakonf,
                              double *flux));

void FORTRAN(calccvel,(ITG *nef,double *vel,double *shcon,
                       ITG *nshcon,ITG *ielmatf,ITG *ntmat_,
                       ITG *mi,double *cvel,double *physcon));

void FORTRAN(calccvelcomp,(ITG *nef,double *vel,double *shcon,
                       ITG *nshcon,ITG *ielmatf,ITG *ntmat_,
                       ITG *mi,double *cvel,double *physcon));

void FORTRAN(calccvfa,(ITG *nface,double *vfa,double *shcon,
                        ITG *nshcon,ITG *ielmatf,ITG *ntmat_,
                       ITG *mi,ITG *ielfa,double *cvfa,double *physcon));

void FORTRAN(calccvfacomp,(ITG *nface,double *vfa,double *shcon,
                        ITG *nshcon,ITG *ielmatf,ITG *ntmat_,
                       ITG *mi,ITG *ielfa,double *cvfa,double *physcon));

void FORTRAN(calcgammak,(ITG *nface,ITG *ielfa,double *vel,double *gradkel,
                        double *gamma,double *xlet,double *xxn,double *xxj,
                         ITG *ipnei,double *betam,ITG *nef,double *flux));

void FORTRAN(calcgammao,(ITG *nface,ITG *ielfa,double *vel,double *gradoel,
                        double *gamma,double *xlet,double *xxn,double *xxj,
                         ITG *ipnei,double *betam,ITG *nef,double *flux));

void FORTRAN(calcgammap,(ITG *nface,ITG *ielfa,double *vel,double *gradpel,
                        double *gamma,double *xlet,double *xxn,double *xxj,
                         ITG *ipnei,double *betam,ITG *nef,double *flux));

void FORTRAN(calcgammat,(ITG *ielfa,double *vel,double *gradtel,
                        double *gamma,double *xlet,double *xxj,
                        ITG *ipnei,double *betam,ITG *nef,double *flux,
                        ITG *nfacea,ITG *nfaceb));

void *calcgammat1mt(ITG *i);

void calcgammatmain(ITG *nface,ITG *ielfa,double *vel,double *gradtel,
                    double *gamma,double *xlet,double *xxj,
                    ITG *ipnei,double *betam,ITG *nef,double *flux,
                    ITG *num_cpus);

void FORTRAN(calcgammav,(ITG *ielfa,double *vel,double *gradvel,
                        double *gamma,double *xlet,double *xxj,
                        ITG *ipnei,double *betam,ITG *nef,double *flux,
                        ITG *nfacea,ITG *nfaceb));

void *calcgammav1mt(ITG *i);

void calcgammavmain(ITG *nface,ITG *ielfa,double *vel,double *gradvel,
                    double *gamma,double *xlet,double *xxj,
                    ITG *ipnei,double *betam,ITG *nef,double *flux,
                    ITG *num_cpus);

void FORTRAN(calchcel,(double *vel,double *cocon,ITG *ncocon,ITG *ielmatf,
                       ITG *ntmat_,ITG *mi,double *hcel,ITG *nef));

void FORTRAN(calchcfa,(ITG *nface,double *vfa,double *cocon,
                        ITG *ncocon,ITG *ielmatf,ITG *ntmat_,
                        ITG *mi,ITG *ielfa,double *hcfa));

void FORTRAN(calch0interface,(ITG *nmpc,ITG *ipompc,ITG *nodempc,
                              double *coefmpc,double *h0));
                      
void FORTRAN(calcmac,(ITG *neq,double *z,double *zz,ITG *nev,double *mac,
                      double* maccpx,ITG *istartnmd,ITG *iendnmd,ITG *nmd,
                      ITG *cyclicsymmetry,ITG *neqact,double *bett,
                      double *betm,ITG *nevcomplex));

void FORTRAN(calcmach,(double *vold,double *voldaux,double *v,int *nk,
           int *ielmat,int *ntmat_,double *shcon,int *nshcon,double *rhcon,
           int *nrhcon,int *iout,int *nmethod,int *convergence,
           double *physcon,int *iponoel,int *inoel,int *ithermal,
           int *nactdoh,int *iit,int *compressible,int *ismooth,
           double *vcontu,double *vtu,int *turbulent,int *inomat,
           int *nodeboun,int *ndirboun,int *nboun,int *mi,double *co,
           double *factor));

void FORTRAN(calcmass,(ITG *ipkon,char *lakon,ITG *kon,double *co,ITG *mi,
             ITG *nelem,ITG *ne,double *thicke,ITG *ielmat,
             ITG *nope,double *t0,double *rhcon,
             ITG *nrhcon,ITG *ntmat_,ITG *ithermal,double *csmass,
             ITG *ielprop,double *prop));

void FORTRAN(calcpel,(ITG *ne,ITG *nactdoh,double *vel,double *b,ITG *nef));

void FORTRAN(calcresfluid,(ITG *n,double *a,double *b,double *au,ITG *ia,
                           ITG *ja,double *x,double *res));

void FORTRAN(calcrespfluid,(ITG *n,double *b,double *au,
                            double *x,double *res));

void FORTRAN(calcrestfluid,(ITG *n,double *a,double *b,double *au,ITG *ia,
                            ITG *ja,double *x,double *res,double *xmin,
                            double *xmax));

void calcrestfluidmain(ITG *n,double *a,double *b,double *au,ITG *ia,
                       ITG *ja,double *x,double *res,ITG *nestart,
                       ITG *num_cpus);

void *calcrestfluidmt(ITG *i);

void FORTRAN(calcresvfluid,(ITG *n,double *a,double *b,double *au,ITG *ia,
                            ITG *ja,double *x,double *res));

void FORTRAN(calcresvfluid1,(ITG *n,double *a,double *b,double *au,ITG *ia,
                             ITG *ja,double *x,double *res));

void *calcresvfluid1mt(ITG *i);

void calcresvfluidmain(ITG *n,double *a,double *b,double *au,ITG *ia,
                       ITG *ja,double *x,double *res,ITG *nestart,
                       ITG *num_cpus,ITG *ncfd);

void FORTRAN(calcresvfluid2,(ITG *n,double *x,double *xmax,ITG *nef));

void *calcresvfluid2mt(ITG *i);

void calcresidual(ITG *nmethod,ITG *neq,double *b,double *fext,double *f,
        ITG *iexpl,ITG *nactdof,double *aux2,double *vold,
        double *vini,double *dtime,double *accold,ITG *nk,double *adb,
        double *aub,ITG *icol,ITG *irow,ITG *nzl,double *alpha,
        double *fextini,double *fini,ITG *islavnode,ITG *nslavnode,
        ITG *mortar,ITG *ntie,
        double *f_cm,double *f_cs,ITG *mi,ITG *nzs,ITG *nasym,
        ITG *idamping,double *veold,double *adc,double *auc,double *cvini,
        double *cv,double *alpham,ITG *num_cpus);

void calcresidual_em(ITG *nmethod,ITG *neq,double *b,double *fext,double *f,
        ITG *iexpl,ITG *nactdof,double *aux1,double *aux2,double *vold,
        double *vini,double *dtime,double *accold,ITG *nk,double *adb,
        double *aub,ITG *icol,ITG *irow,ITG *nzl,double *alpha,
        double *fextini,double *fini,ITG *islavnode,ITG *nslavnode,
        ITG *mortar,ITG *ntie,
        double *f_cm,double *f_cs,ITG *mi,ITG *nzs,ITG *nasym,ITG *ithermal);

void FORTRAN(calcrhoel,(ITG *nef,double *vel,double *rhcon,
                        ITG *nrhcon,ITG *ielmatf,ITG *ntmat_,ITG *ithermal,
                        ITG *mi));

void FORTRAN(calcrhoel1,(ITG *nef,double *vel,double *rhcon,
                         ITG *nrhcon,ITG *ielmatf,ITG *ntmat_,ITG *ithermal,
                         ITG *mi,ITG *nefa,ITG *nefb));

void *calcrhoel1mt(ITG *i);

void calcrhoelmain(ITG *nef,double *vel,double *rhcon,
                   ITG *nrhcon,ITG *ielmatf,ITG *ntmat_,ITG *ithermal,
                   ITG *mi,ITG *num_cpus);

void FORTRAN(calcrhoelcomp,(ITG *nef,double *vel,double *shcon,
                             ITG *ielmatf,ITG *ntmat_,ITG *mi,
                             ITG *nefa,ITG *nefb));

void *calcrhoelcomp1mt(ITG *i);

void calcrhoelcompmain(ITG *nef,double *vel,double *shcon,
                                ITG *ielmatf,ITG *ntmat_,ITG *mi,
                                ITG *num_cpus);

void FORTRAN(calcrhofa,(ITG *nface,double *vfa,double *rhcon,
                        ITG *nrhcon,ITG *ielmatf,ITG *ntmat_,ITG *ithermal,
                        ITG *mi,ITG *ielfa));

void *calcrhofa1mt(ITG *i);

void calcrhofamain(ITG *nface,double *vfa,double *rhcon,
                   ITG *nrhcon,ITG *ielmatf,ITG *ntmat_,ITG *ithermal,
                   ITG *mi,ITG *ielfa,ITG *num_cpus);

void FORTRAN(hrr,(double *vfa,double *shcon,
                            ITG *ielmatf,ITG *ntmat_,ITG *mi,
                            ITG *ielfa,ITG *nfacea,ITG *nfaceb));

void *hrr1mt(ITG *i);

void hrrmain(ITG *nface,double *vfa,double *shcon,
                       ITG *ielmatf,ITG *ntmat_,ITG *mi,
                       ITG *ielfa,ITG *num_cpus);

void FORTRAN(hrr_mod_smart,(double *vfa,double *shcon,
                            ITG *ielmatf,ITG *ntmat_,ITG *mi,ITG *ielfa,
                            ITG *ipnei,double *vel,ITG *nef,double *flux,
                            double *gradpel,double *gradtel,double *xxj,
                            double *xlet,ITG *nfacea,ITG *nfaceb));

void *hrr_mod_smart1mt(ITG *i);

void hrr_mod_smartmain(ITG *nface,double *vfa,double *shcon,
                                 ITG *ielmatf,ITG *ntmat_,ITG *mi,ITG *ielfa,
                                 ITG *ipnei,double *vel,ITG *nef,double *flux,
                                 double *gradpel,double *gradtel,double *xxj,
                                 double *xlet,ITG *num_cpus);

void FORTRAN(hrr_ud,(double *vfa,double *shcon,
                                ITG *ielmatf,ITG *ntmat_,ITG *mi,ITG *ielfa,
                                ITG *ipnei,double *vel,ITG *nef,double *flux,
                                ITG *nfacea,ITG *nfaceb,double *xxi1,
                                double *xle1,double *gradpel1,double *gradtel1,
                                ITG *neij1));

void *hrr_ud1mt(ITG *i);

void hrr_udmain(ITG *nface,double *vfa,double *shcon,
                          ITG *ielmatf,ITG *ntmat_,ITG *mi,ITG *ielfa,
                          ITG *ipnei,double *vel,ITG *nef,double *flux,
                          ITG *num_cpus,double *xxi,double *xle,
                          double *gradpel,double *gradtel,ITG *neij);

void calcshapef(int *nvar_,int *ipvar,double **var,int *ne,
             char *lakon,double *co,int *ipkon,int *kon,
             int *nelemface,char *sideface,int *nface,
             int *nvarf_,int *ipvarf,double **varfp);

void FORTRAN(calcstabletimeinccont,(ITG *ne,char *lakon,ITG *kon,ITG *ipkon,
             ITG *mi,ITG *ielmat,double *elcon,ITG *mortar,double *adb,
             double *alpha,ITG *nactdof,double *springarea,ITG *ne0,
             ITG *ntmat_,ITG *ncmat_,double *dtcont,double *smscale,
             double *dtset,ITG *mscalmethod));

void FORTRAN(calcstabletimeincvol,(ITG *ne0,double *elcon,ITG *nelcon,
            double *rhcon,ITG *nrhcon,double *alcon,ITG *nalcon,double *orab,
            ITG *ntmat_,ITG *ithermal,double *alzero,double *plicon,
            ITG *nplicon,double *plkcon,ITG *nplkcon,ITG *npmat_,ITG *mi,
            double *dtime,double *xstiff,ITG *ncmat_,double *vold,ITG *ielmat,
            double *t0,double *t1,char *matname,char *lakon,
            double *xmatwavespeed,ITG *nmat,ITG *ipkon,double *co,ITG *kon,
            double *dtvol,double *alpha,double *smscale,double *dtset,
            ITG *mscalmethod));

void FORTRAN(calcstressheatflux,(double *sti,double *umel,double *gradvel,
             double *qfx,double *hcel,double *gradtel,ITG *nef,ITG *isti,
             ITG *iqfx,ITG *mi));

void FORTRAN(calcstressheatfluxfem,(int *kon,char *lakon,int *ipkon,int *ielmat,
             int *ntmat_,double *vold,char *matname,int *mi,double *shcon,
             int *nshcon,int *turbulent,int *compressible,int *ipvar,
             double *var,double *sti,double *qfx,double *cocon,int *ncocon,
             int *ne,int *isti,int *iqfx));

void FORTRAN(calcttel,(ITG *nef,double *vel,double *shcon,ITG *nshcon,
                       ITG *ielmatf,ITG *ntmat_,ITG *mi,double *physcon,
                       double *ttel));

void FORTRAN(calcttfaext,(ITG *nfaext,double *vfa,double *shcon,ITG *nshcon,
                          ITG *ielmatf,ITG *ntmat_,ITG *mi,double *physcon,
                          double *ttfa,ITG *ifaext,ITG *ielfa));

void FORTRAN(calculated,(ITG *nktet,double *d,double *dmin,
                             ITG *ipoed,ITG *iedg,double *cotet));

void FORTRAN(calculateh,(ITG *nk,double *v,double *veold,double *stn,
                         double *een,double *emn,double *epn,double *enern,
                         double *qfn,double *errn,double *h,char *filab,
                         ITG *mi,double *d,ITG *nh,double *dmin,ITG *ipoed,
                         ITG *iedg,double *cotet));

void FORTRAN(calcumel,(ITG *nef,double *vel,double *shcon,
                       ITG *nshcon,ITG *ielmatf,ITG *ntmat_,ITG *ithermal,
                       ITG *mi,double *umel));

void FORTRAN(calcumfa,(ITG *nface,double *vfa,double *shcon,
                        ITG *nshcon,ITG *ielmatf,ITG *ntmat_,ITG *ithermal,
                        ITG *mi,ITG *ielfa,double *umfa));

void FORTRAN(calcvel,(ITG *ne,ITG *nactdoh,double *vel,double *b,
                      ITG *neq,ITG *nef));

void FORTRAN(calcview,(char *sideload,double *vold,double *co,
             double *pmid,double *e1,double *e2,double *e3,
             ITG *kontri,ITG *nloadtr,double *adview,double *auview,
             double *dist,ITG *idist,double *area,ITG *ntrit,ITG *mi,ITG *jqrad,
             ITG *irowrad,ITG *nzsrad,double *sidemean,ITG *ntria,
             ITG *ntrib,char *covered,ITG *ng));

void *calcviewmt(ITG *i);

void FORTRAN(calinput,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
               ITG *nkon,ITG *ne,
               ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
               ITG *ipompc,ITG *nodempc,double *coefmpc,ITG *nmpc,
               ITG *nmpc_,ITG *nodeforc,ITG *ndirforc,double *xforc,
               ITG *nforc,ITG *nforc_,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,ITG *nload_,
               ITG *nprint,char *prlab,char *prset,ITG *mpcfree,ITG *nboun_,
               ITG *mei,char *set,ITG *istartset,
               ITG *iendset,ITG *ialset,ITG *nset,ITG *nalset,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,double *t0,
               double *t1,char *matname,
               ITG *ielmat,char *orname,double *orab,ITG *ielorien,
               char *amname,double *amta,ITG *namta,ITG *nam,
               ITG *nmethod,ITG *iamforc,ITG *iamload,
               ITG *iamt1,ITG *ithermal,ITG *iperturb,
               ITG *istat,ITG *istep,ITG *nmat,
               ITG *ntmat_,ITG *norien,double *prestr,ITG *iprestr,
               ITG *isolver,double *fei,double *veold,double *timepar,
               double *xmodal,char *filab,ITG *jout,ITG *nlabel,
               ITG *idrct,ITG *jmax,ITG *iexpl,double *alpha,ITG *iamboun,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               ITG *iplas,ITG *npmat_,ITG *mi,ITG *nk_,
               double *trab,ITG *inotr,ITG *ntrans,ITG *ikboun,
               ITG *ilboun,ITG *ikmpc,ITG *ilmpc,ITG *ics,
               double *dcs,ITG *ncs_,ITG *namtot_,double *cs,
               ITG *nstate_,ITG *ncmat_,ITG *iumat,ITG *mcs,
               char *labmpc,ITG *iponor,double *xnor,ITG *knor,
               double *thickn,double *thicke,ITG *ikforc,ITG *ilforc,
               double *offset,ITG *iponoel,ITG *inoel,ITG *rig,
               ITG *infree,ITG *nshcon,double *shcon,double *cocon,
               ITG *ncocon,double *physcon,ITG *nflow,double *ctrl,
               ITG *maxlenmpc,ITG *ne1d,ITG *ne2d,ITG *nener,
               double *vold,ITG *nodebounold,
               ITG *ndirbounold,double *xbounold,double *xforcold,
               double *xloadold,double *t1old,double *eme,
               double *sti,double *ener,
               double *xstate,char *jobnamec,ITG *irstrt,
               double *ttime,double *qaold,
               char *output,char *typeboun,char *inpc,
               ITG *ipoinp,ITG *inp,char *tieset,double *tietol,
               ITG *ntie,double *fmpc,char *cbody,ITG *ibody,double *xbody,
               ITG *nbody,ITG *nbody_,double *xbodyold,ITG *nam_,
               ITG *ielprop,ITG *nprop,ITG *nprop_,double *prop,
               ITG *itpamp,ITG *iviewfile,ITG *ipoinpc,
               ITG *nslavs,double *t0g,double *t1g,ITG *network,
               ITG *cyclicsymmetry,ITG *idefforc,ITG *idefload,
               ITG *idefbody,ITG *mortar,ITG *ifacecount,ITG *islavsurf,
               double *pslavsurf,double *clearini,char *heading,
               ITG *iaxial,ITG *nobject,char *objectset,ITG *nprint_,
               ITG *iuel,ITG *nuel_,ITG *nodempcref,double *coefmpcref,
               ITG *ikmpcref,ITG *memmpcref_,ITG *mpcfreeref,
               ITG *maxlenmpcref,ITG *memmpc_,ITG *isens,ITG *namtot,
               ITG *stam,double *dacon,double *vel,ITG *nef,
               double *velo,double *veloo,ITG *ne2boun,ITG *itempuser));    

void cascade(ITG *ipompc,double **coefmpcp,ITG **nodempcp,ITG *nmpc,
   ITG *mpcfree,ITG *nodeboun,ITG *ndirboun,ITG*nboun,ITG*ikmpc,
   ITG *ilmpc,ITG *ikboun,ITG *ilboun,ITG *mpcend,
   char *labmpc,ITG *nk,ITG *memmpc_,ITG *icascade,ITG *maxlenmpc,
   ITG *callfrommain,ITG *iperturb,ITG *ithermal);

void cascadefem(int *ipompc,double **coefmpcp,int **nodempcp,int *nmpc,
   int *mpcfree,int *nodeboun,int *ndirboun,int*nboun,int*ikmpc,
   int *ilmpc,int *ikboun,int *ilboun,int *mpcend,int *mpcmult,
   char *labmpc,int *nk,int *memmpc_,int *icascade,int *maxlenmpc,
   int *callfrommain,int *iperturb,int *ithermal);

void FORTRAN(cataloguenodes,(ITG *iponofa,ITG *inofa,ITG *ifreefa,ITG *ielfa,
             ITG *ifaboun,ITG *ipkon,ITG *kon,char *lakon,ITG *nface,ITG *ne));
    
void FORTRAN(catedges,(ITG *netet_,ITG *iedg,ITG *kontet,ITG *ipoed,
                       ITG *ifreeed,ITG *iedtet,ITG *ipoeled,ITG *ieled,
                       ITG *ifreele));

void FORTRAN(cattet,(ITG *kontet,ITG *netet_,ITG *ifac,ITG *ne,ITG *ipkon,
                     ITG *kon,ITG *ifatet,ITG *ifreetet,double *bc,
                     ITG *itetfa,ITG *ifreefa,double *planfa,ITG *ipofa,
                     double *cotet,double *cg,ITG *ipoeln,ITG *ieln,
                     ITG *ifreeln,char *lakon,ITG *kontetor,ITG *iquad));

void FORTRAN(cavity,(ITG *kontet,ITG *ifatet,ITG *ifreetet,double *bc,
                     ITG *ifac,ITG *itetfa,ITG *ifreefa,double *planfa,
                     ITG *ipofa,double *cotet,ITG *ibase,ITG *node,
                     ITG *iexternfa,ITG *ipoed,ITG *iedg,ITG *ifreeed,
                     ITG *ipoeled,ITG *ieled,ITG *ifreele,ITG *nktet,
                     ITG *netet_,ITG *ibasenewnodes,double *conewnodes,
                     ITG *ipoeln,ITG *ieln,ITG *ifreeln,ITG *nnewnodes,
                     ITG *iexternedg,ITG *iedtet,double *cg,
                     double *height));

void FORTRAN(cavityext,(ITG *kontet,ITG *ifatet,ITG *ifreetet,double *bc,
                     ITG *ifac,ITG *itetfa,ITG *ifreefa,double *planfa,
                     ITG *ipofa,double *cotet,ITG *ibase,ITG *node,
                     ITG *iexternfa,ITG *ipoed,ITG *iedg,ITG *ifreeed,
                     ITG *ipoeled,ITG *ieled,ITG *ifreele,ITG *nktet,
                     ITG *netet_,ITG *ibasenewnodes,double *conewnodes,
                     ITG *ipoeln,ITG *ieln,ITG *ifreeln,ITG *nnewnodes,
                     ITG *iexternedg,ITG *iedtet,double *cotetorig,
                     ITG *iedge,ITG *iexternnode,double *cg,
                     double *hnewnodes));

void FORTRAN(cfdconv,(double *vold,double *voldaux,double *v,int *nk,
           int *ielmat,int *ntmat_,double *shcon,int *nshcon,double *rhcon,
           int *nrhcon,int *iout,int *nmethod,int *convergence,
           double *physcon,int *iponoel,int *inoel,int *ithermal,
           int *nactdoh,int *iit,int *compressible,int *ismooth,
           double *vcontu,double *vtu,int *turbulent,int *inomat,
           int *nodeboun,int *ndirboun,int *nboun,int *mi,double *co,
           double *factor,double *vconini,double *dtimef,double *del,
           double *sum,double *sumx,double *sumxx,double *sumy,
           double *sumxy,int *nstart,double *shockcoef));

ITG cgsolver(double *A,double *x,double *b,ITG neq,ITG len,ITG *ia,ITG *iz,
                                double *eps,ITG *niter,ITG precFlg);

void checkconvergence(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
          ITG *ne,double *stn,ITG *nmethod,
          ITG *kode,char *filab,double *een,double *t1act,
          double *time,double *epn,ITG *ielmat,char *matname,
          double *enern,double *xstaten,ITG *nstate_,ITG *istep,
          ITG *iinc,ITG *iperturb,double *ener,ITG *mi,char *output,
          ITG *ithermal,double *qfn,ITG *mode,ITG *noddiam,double *trab,
          ITG *inotr,ITG *ntrans,double *orab,ITG *ielorien,ITG *norien,
          char *description,double *sti,
          ITG *icutb,ITG *iit,double *dtime,double *qa,double *vold,
          double *qam,double *ram1,double *ram2,double *ram,
          double *cam,double *uam,ITG *ntg,double *ttime,
          ITG *icntrl,double *theta,double *dtheta,double *veold,
          double *vini,ITG *idrct,double *tper,ITG *istab,double *tmax,
          ITG *nactdof,double *b,double *tmin,double *ctrl,double *amta,
          ITG *namta,ITG *itpamp,ITG *inext,double *dthetaref,ITG *itp,
          ITG *jprint,ITG *jout,ITG *uncoupled,double *t1,ITG *iitterm,
          ITG *nelemload,ITG *nload,ITG *nodeboun,ITG *nboun,ITG *itg,
          ITG *ndirboun,double *deltmx,ITG *iflagact,char *set,ITG *nset,
          ITG *istartset,ITG *iendset,ITG *ialset,double *emn,double *thicke,
          char *jobnamec,ITG *mortar,ITG *nmat,ITG *ielprop,double *prop,
          ITG *ialeatoric,ITG *kscale,
          double *energy,double *allwk,double *energyref,
          double *emax,double *enres,double *enetoll,double *energyini,
          double *allwkini,double *temax,double *reswk,ITG *ne0,
          ITG *neini,double *dampwk,double *dampwkini,double *energystartstep);

void checkconvnet(ITG *icutb,ITG *iin,
                  double *cam1t,double *cam1f,double *cam1p,
                  double *cam2t,double *cam2f,double *cam2p,
                  double *camt,double *camf,double *camp,
                  ITG *icntrl,double *dtheta,double *ctrl,
                  double *cam1a,double *cam2a,double *cama,
                  double *vamt,double *vamf,double *vamp,double *vama,
                  double *qa,double *qamt,double *qamf,double *ramt,
                  double *ramf,double *ramp,ITG *iplausi,ITG *ichannel);

void FORTRAN(checkconstraint,(ITG *nobject,char *objectset,double *g0,
             ITG *nactive,ITG *nnlconst,ITG *ipoacti,ITG *ndesi,
             double *dgdxglob,ITG *nk,ITG *nodedesi,ITG *iconstacti,
             double *objnorm,ITG *inameacti));

void checkdivergence(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
          ITG *ne,double *stn,ITG *nmethod,
          ITG *kode,char *filab,double *een,double *t1act,
          double *time,double *epn,ITG *ielmat,char *matname,
          double *enern,double *xstaten,ITG *nstate_,ITG *istep,
          ITG *iinc,ITG *iperturb,double *ener,ITG *mi,char *output,
          ITG *ithermal,double *qfn,ITG *mode,ITG *noddiam,double *trab,
          ITG *inotr,ITG *ntrans,double *orab,ITG *ielorien,ITG *norien,
          char *description,double *sti,
          ITG *icutb,ITG *iit,double *dtime,double *qa,double *vold,
          double *qam,double *ram1,double *ram2,double *ram,
          double *cam,double *uam,ITG *ntg,double *ttime,
          ITG *icntrl,double *theta,double *dtheta,double *veold,
          double *vini,ITG *idrct,double *tper,ITG *istab,double *tmax,
          ITG *nactdof,double *b,double *tmin,double *ctrl,double *amta,
          ITG *namta,ITG *itpamp,ITG *inext,double *dthetaref,ITG *itp,
          ITG *jprint,ITG *jout,ITG *uncoupled,double *t1,ITG *iitterm,
          ITG *nelemload,ITG *nload,ITG *nodeboun,ITG *nboun,ITG *itg,
          ITG *ndirboun,double *deltmx,ITG *iflagact,char *set,ITG *nset,
          ITG *istartset,ITG *iendset,ITG *ialset,double *emn,double *thicke,
          char *jobnamec,ITG *mortar,ITG *nmat,ITG *ielprop,double *prop,
          ITG *ialeatoric,ITG *kscale,
          double *energy,double *allwk,double *energyref,
          double *emax,double *enres,double *enetoll,double *energyini,
          double *allwkini,double *temax,double *reswk,ITG *ne0,
          ITG *neini,double *dampwk,double *dampwkini,double *energystartstep);

void checkinclength(double *time,double *ttime,double *theta,double *dtheta,
          ITG *idrct,double *tper,double *tmax,double *tmin,double *ctrl,
          double *amta,ITG *namta,ITG *itpamp,ITG *inext,double *dthetaref,
          ITG *itp,ITG *jprint,ITG *jout);
         
void FORTRAN(checkimpacts,(ITG *ne,ITG *neini,double *temax,
                   double *sizemaxinc,double *energyref,double *tmin,
                   double *tper,
                   ITG *idivergence,ITG *idirinctime,ITG *istab,
                   double *dtheta,double *enres,double *energy,
                   double *energyini,double *allwk,double *allwkini,
                   double *dampwk,double *dampwkini,double *emax,
                   ITG *mortar,double *maxdecay,double *enetoll));

void FORTRAN(checkinputvaluesnet,(ITG *ieg,ITG *nflow,double *prop,
                                  ITG *ielprop,char *lakon));

void FORTRAN(checkprojectgrad,(ITG *nactiveold,ITG *nactive,ITG *ipoacti,
                               ITG *ipoactiold,char *objectset,double *lambda,
                               ITG *nnlconst,ITG *iconstacti,ITG *iconstactiold,
                               ITG *inameacti,ITG *inameactiold,double *g0,
		               ITG *nobject));

void FORTRAN(checktime,(ITG *itpamp,ITG *namta,double *tinc,double *ttime,
             double *amta,double *tmin,ITG *inext,ITG *itp,ITG *istep,
             double *tper));

void FORTRAN(checktruecontact,(ITG *ntie,char *tieset,double *tietol,
             double *elcon,ITG *itruecontact,ITG *ncmat_,ITG *ntmat_));

void FORTRAN(closefile,());

void FORTRAN(closefilefluid,());
  
void FORTRAN(compdt,(int *nk,double *dt,int *nshcon,
       double *shcon,int *nrhcon,double *rhcon,double *vold,
       int *ntmat_,int *iponoel,int *inoel,double *dtimef,int *iexplicit,
       int *ielmat,double *physcon,double *dh,double *cocon,int *ncocon,
       int *ithermal,int *mi,int *ipkon,int *kon,char *lakon,
       int *ne,double *v,double *co,int *turbulent,
       double *vcontu,double *vcon));

void compfluid(double **cop,ITG *nk,ITG **ipkonp,ITG *konf,char **lakonp,
    char **sideface,ITG *ifreestream,
    ITG *nfreestream,ITG *isolidsurf,ITG *neighsolidsurf,
    ITG *nsolidsurf,ITG *nshcon,double *shcon,
    ITG *nrhcon,double *rhcon,double **voldp,ITG *ntmat_,ITG *nodeboun,
    ITG *ndirboun,ITG *nboun,ITG *ipompc,ITG *nodempc,ITG *nmpc,
    ITG *ikmpc,ITG *ilmpc,ITG *ithermal,ITG *ikboun,ITG *ilboun,
    ITG *turbulent,ITG *isolver,ITG *iexpl,double *ttime,
    double *time,double *dtime,ITG *nodeforc,ITG *ndirforc,double *xforc,
    ITG *nforc,ITG *nelemload,char *sideload,double *xload,ITG *nload,
    double *xbody,ITG *ipobody,ITG *nbody,ITG *ielmatf,char *matname,
    ITG *mi,ITG *ncmat_,double *physcon,ITG *istep,ITG *iinc,
    ITG *ibody,double *xloadold,double *xboun,
    double *coefmpc,ITG *nmethod,double *xforcold,double *xforcact,
    ITG *iamforc,ITG *iamload,double *xbodyold,double *xbodyact,
    double *t1old,double *t1,double *t1act,ITG *iamt1,double *amta,
    ITG *namta,ITG *nam,double *ampli,double *xbounold,double *xbounact,
    ITG *iamboun,ITG *itg,ITG *ntg,char *amname,double *t0,ITG **nelemface,
    ITG *nface,double *cocon,ITG *ncocon,double *xloadact,double *tper,
    ITG *jmax,ITG *jout,char *set,ITG *nset,ITG *istartset,
    ITG *iendset,ITG *ialset,char *prset,char *prlab,ITG *nprint,
    double *trab,ITG *inotr,ITG *ntrans,char *filab,char *labmpc,
    double *sti,ITG *norien,double *orab,char *jobnamef,char *tieset,
    ITG *ntie,ITG *mcs,ITG *ics,double *cs,ITG *nkon,ITG *mpcfree,
    ITG *memmpc_,double *fmpc,ITG *nef,ITG **inomat,double *qfx,
    ITG *neifa,ITG *neiel,ITG *ielfa,ITG *ifaext,double *vfa,double *vel,
    ITG *ipnei,ITG *nflnei,ITG *nfaext,char *typeboun,ITG *neij,
    double *tincf,ITG *nactdoh,ITG *nactdohinv,ITG *ielorien,
    char *jobnamec,ITG *ifatie,ITG *nstate_,double *xstate,char *orname,
    ITG *kon,double *ctrl,ITG *kode,double *velo,double *veloo,
    ITG *initial);

void compfluidfem(double **cop,int *nk,int **ipkonp,int **konp,char **lakonp,
    int *ne,char **sideface,int *ifreestream,
    int *nfreestream,int *isolidsurf,int *neighsolidsurf,
    int *nsolidsurf,int **iponoel,int **inoel,int *nshcon,double *shcon,
    int *nrhcon,double *rhcon,double **voldp,int *ntmat_,int *nodeboun,
    int *ndirboun,int *nboun,int **ipompcp,int **nodempcp,int *nmpc,
    int **ikmpcp,int **ilmpcp,int *ithermal,int *ikboun,int *ilboun,
    int *turbulent,int *isolver,int *iexpl,double *vcontu,double *ttime,
    double *time,double *dtime,int *nodeforc,int *ndirforc,double *xforc,
    int *nforc,int *nelemload,char *sideload,double *xload,int *nload,
    double *xbody,int *ipobody,int *nbody,int **ielmatp,char *matname,
    int *mi,int *ncmat_,double *physcon,int *istep,int *iinc,
    int *ibody,double *xloadold,double *xboun,
    double **coefmpcp,int *nmethod,double *xforcold,double *xforcact,
    int *iamforc,int *iamload,double *xbodyold,double *xbodyact,
    double *t1old,double *t1,double *t1act,int *iamt1,double *amta,
    int *namta,int *nam,double *ampli,double *xbounold,double *xbounact,
    int *iamboun,int *itg,int *ntg,char *amname,double *t0,int **nelemface,
    int *nface,double *cocon,int *ncocon,double *xloadact,double *tper,
    int *jmax,int *jout,char *set,int *nset,int *istartset,
    int *iendset,int *ialset,char *prset,char *prlab,int *nprint,
    double *trab,int *inotr,int *ntrans,char *filab,char **labmpcp,
    double *sti,int *norien,double *orab,char *jobnamef,char *tieset,
    int *ntie,int *mcs,int *ics,double *cs,int *nkon,int *mpcfree,
    int *memmpc_,double **fmpcp,int *nef,int **inomat,double *qfx);

void FORTRAN(complete_hel,(ITG *neq,double *b,double *hel,double *ad,
             double *au,ITG *jq,ITG *irow,ITG *nzs));

void FORTRAN(complete_hel_blk,(double *vel,double *hel,double *auv6,
             ITG *ipnei,ITG *neiel,ITG *nef,ITG *nactdohinv));

void FORTRAN(complete_hel_cyclic,(ITG *neq,double *b,double *hel,double *ad,
             double *au,ITG *jq,ITG *irow,ITG *ipnei,ITG *neiel,
             ITG *ifatie,double *c,char *lakonf,ITG *neifa,ITG *nzs));

void FORTRAN(complete_hel_cyclic_blk,(double *vel,double *hel,double *auv6,
             double *c,ITG *ipnei,ITG *neiel,ITG *neifa,ITG *ifatie,
             ITG *nef));

void complete_hel_blk_main(double *vel,double *hel,double *auv6,double *c,
      ITG *ipnei,ITG *neiel,ITG *neifa,ITG *ifatie,ITG *nef);

void complexfreq(double **cop,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,ITG *ne,
               ITG **nodebounp,ITG **ndirbounp,double **xbounp,ITG *nboun,
               ITG **ipompcp,ITG **nodempcp,double **coefmpcp,char **labmpcp,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,double *xforc,
               ITG *nforc,ITG *nelemload,char *sideload,double *xload,
               ITG *nload,
               ITG **nactdofp,ITG *neq,ITG *nzl,ITG *icol,ITG *irow,
               ITG *nmethod,ITG **ikmpcp,ITG **ilmpcp,ITG **ikbounp,
               ITG **ilbounp,double *elcon,ITG *nelcon,double *rhcon,
               ITG *nrhcon,double *cocon,ITG *ncocon,
               double *alcon,ITG *nalcon,double *alzero,
               ITG **ielmatp,ITG **ielorienp,ITG *norien,double *orab,
               ITG *ntmat_,double **t0p,
               double **t1p,ITG *ithermal,double *prestr,ITG *iprestr,
               double **voldp,ITG *iperturb,double **stip,ITG *nzs,
               double *timepar,double *xmodal,
               double **veoldp,char *amname,double *amta,
               ITG *namta,ITG *nam,ITG *iamforc,ITG *iamload,
               ITG **iamt1p,ITG *jout,
               ITG *kode,char *filab,double **emep,double *xforcold,
               double *xloadold,
               double **t1oldp,ITG **iambounp,double **xbounoldp,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstate,ITG *npmat_,char *matname,ITG *mi,
               ITG *ncmat_,ITG *nstate_,double **enerp,char *jobnamec,
               double *ttime,char *set,ITG *nset,ITG *istartset,
               ITG *iendset,ITG **ialsetp,ITG *nprint,char *prlab,
               char *prset,ITG *nener,double *trab,
               ITG **inotrp,ITG *ntrans,double **fmpcp,char *cbody,ITG *ibody,
               double *xbody,ITG *nbody,double *xbodyold,ITG *istep,
               ITG *isolver,ITG *jq,char *output,ITG *mcs,ITG *nkon,
               ITG *mpcend,ITG *ics,double *cs,ITG *ntie,char *tieset,
               ITG *idrct,ITG *jmax,
               double *ctrl,ITG *itpamp,double *tietol,ITG *nalset,
               ITG *ikforc,ITG *ilforc,double *thicke,
               char *jobnamef,ITG *mei,ITG *nmat,ITG *ielprop,double *prop,
               char *orname,char *typeboun);

void FORTRAN(con2phys,(double *vold,double *voldaux,double *v,int *nk,
           int *ielmat,int *ntmat_,double *shcon,int *nshcon,double *rhcon,
           int *nrhcon,int *iout,int *nmethod,int *convergence,
           double *physcon,int *iponoel,int *inoel,int *ithermal,
           int *nactdoh,int *iit,int *compressible,int *ismooth,
           double *vcontu,double *vtu,int *turbulent,int *inomat,
           int *nodeboun,int *ndirboun,int *nboun,int *mi,double *co,
           double *factor));

void contact(ITG *ncont,ITG *ntie,char *tieset,ITG *nset,char *set,
             ITG *istartset,ITG *iendset,ITG *ialset,ITG *itietri,
             char *lakon,ITG *ipkon,ITG *kon,ITG *koncont,ITG *ne,
             double *cg,double *straight,ITG *ifree,double *co,
             double *vold,ITG *ielmat,double *cs,double *elcon,
             ITG *istep,ITG *iinc,ITG *iit,ITG *ncmat_,ITG *ntmat_,
             ITG *ne0,double *vini,
             ITG *nmethod,
             ITG *iperturb,ITG *ikboun,ITG *nboun,ITG *mi,ITG *imastop,
             ITG *nslavnode,ITG *islavnode,ITG *islavsurf,ITG *itiefac,
             double *areaslav,ITG *iponoels,ITG *inoels,double *springarea,
             double *tietol,double *reltime,ITG *imastnode,ITG *nmastnode,
             double *xmastnor,char *filab,ITG *mcs,
             ITG *ics,ITG *nasym,double *xnoels,ITG *mortar,
             double *pslavsurf,double *pmastsurf,double *clearini,
             double *theta,double *xstateini,double *xstate,ITG *nstate_,
             ITG *icutb,ITG *ialeatoric,char *jobnamef,double *alea);

void FORTRAN(convert2slapcol,(double *au,double *ad,
                              ITG *jq,ITG *nzs,ITG *nef,double *aua));
    
void FORTRAN(coriolissolve,(double *cc,ITG *nev,double *aa,double *bb,
             double *xx,double *eiga,double *eigb,double *eigxx,
             ITG *iter,double *d,double *temp,ITG *nevcomplex));

void FORTRAN(correctflux,(ITG *nef,ITG *ipnei,ITG *neifa,ITG *neiel,
                          double *flux,
                          double *vfa,double *advfa,double *area,double *vel,
                          double *alet,ITG *ielfa,double *ale,ITG *ifabou,
                          ITG *nefa,ITG *nefb,double *xxnj,double *gradpcfa));

void *correctflux1mt(ITG *i);

void correctfluxmain(ITG *nef,ITG *ipnei,ITG *neifa,ITG *neiel,
                     double *flux,
                     double *vfa,double *advfa,double *area,double *vel,
                     double *alet,ITG *ielfa,double *ale,ITG *ifabou,
                     ITG *num_cpus,double *xxnj,double *gradpcfa);

void FORTRAN(correctfluxcomp,(ITG *nef,ITG *ipnei,ITG *neifa,ITG *neiel,
                              double *flux,
                              double *vfa,double *advfa,double *area,
                              double *vel,
                              double *alet,ITG *ielfa,double *ale,ITG *ifabou,
                              ITG *ielmatf,ITG *mi,double *shcon,ITG *ntmat_,
                              ITG *nefa,ITG *nefb,double *xxnj,
                              double *gradpcfa,ITG *inlet));

void *correctfluxcomp1mt(ITG *i);

void correctfluxcompmain(ITG *nef,ITG *ipnei,ITG *neifa,ITG *neiel,
                         double *flux,
                         double *vfa,double *advfa,double *area,double *vel,
                         double *alet,ITG *ielfa,double *ale,ITG *ifabou,
                         ITG *ielmatf,ITG *mi,double *shcon,ITG *ntmat_,
                         ITG *num_cpus,double *xxnj,double *gradpcfa,
                         ITG *inlet);

void FORTRAN(correctrhoelcomp,(ITG *nef,double *vel,double *shcon,
                        ITG *ielmatf,ITG *ntmat_,ITG *mi));

void FORTRAN(correctvel,(double *adv,ITG *nef,double *volume,double *gradpcel,
                          double *vel,ITG *nefa,ITG *nefb));

void *correctvel1mt(ITG *i);

void correctvelmain(double *adv,ITG *nef,double *volume,double *gradpcel,
                    double *vel,ITG *num_cpus);

void FORTRAN(correctvel_simplec,(double *adv,ITG *nef,double *volume,double *gradpel,
                                 double *vel,ITG *ipnei,double *auv));

void cpypardou(double *var1,double *var2,ITG *size,ITG *num_cpus);

void *cpypardoumt(ITG *i);

void cpyparitg(ITG *iva1,ITG *iva2,ITG *size,ITG *num_cpus);

void *cpyparitgmt(ITG *i);

void FORTRAN(createblock,(ITG *nef,ITG *ipnei,ITG *neiel,ITG *iam,
                          ITG *jam,ITG *iamorig,ITG *nflnei,ITG *nz_num,
                          ITG *num_cpus,ITG *nestart,ITG *ineighblock,
                          ITG *neighblock,ITG *icyclic,ITG *neifa,
                          ITG *ifatie,ITG *nneighblock));

void FORTRAN(createblock_struct,(ITG *neq,ITG *ipointers,ITG *icolpardiso,
                                double *aupardiso,ITG *nestart,ITG *num_cpus,
                                ITG *ja,ITG *nz_num));

void FORTRAN(createelemneigh,(ITG *nk,ITG *iponoel,ITG *inoel,
                          ITG *istartnneigh,ITG *ialnneigh,
                          ITG *icheckelems,ITG *istarteneigh,
                          ITG *ialeneigh));

void FORTRAN(createfint,(ITG *ne,ITG *ipkon,char *lakon,ITG *kon,
                         ITG *nactdof,ITG *mi,double *fn0,double *fint));

void FORTRAN(createialnk,(ITG *nk,ITG *iponoel,ITG *inoel,ITG *istartnk,
                          ITG *ialnk,ITG *ipkon));

void FORTRAN(createinterfacempcs,(ITG *imastnode,double *xmastnor,
             ITG *nmastnode,ITG *ikmpc,ITG *ilmpc,ITG *nmpc,ITG *ipompc,
             ITG *nodempc,double *coefmpc,char *labmpc,ITG *mpcfree,
             ITG *ikboun,ITG *nboun));

void FORTRAN(createinum,(ITG *ipkon,ITG *inum,ITG *kon,char *lakon,ITG *nk,
             ITG *ne,char *cflag,ITG *nelemload,ITG *nload,ITG *nodeboun,
             ITG *nboun,ITG *ndirboun,ITG *ithermal,double *co,
             double *vold,ITG *mi,ITG *ielmat,ITG *ielprop,double *prop));

void FORTRAN(createmddof,(ITG *imddof,ITG *nmddof,ITG *istartset,
       ITG *iendset,ITG *ialset,ITG *nactdof,ITG *ithermal,ITG *mi,
       ITG *imdnode,ITG *nmdnode,ITG *ikmpc,ITG *ilmpc,ITG *ipompc,
       ITG *nodempc,ITG *nmpc,ITG *imdmpc,
       ITG *nmdmpc,ITG *imdboun,ITG *nmdboun,ITG *ikboun,ITG *nboun,
       ITG *nset,ITG *ntie,char *tieset,char *set,char *lakon,ITG *kon,
       ITG *ipkon,char *labmpc,ITG *ilboun,char *filab,char *prlab,
       char *prset,ITG *nprint,ITG *ne,ITG *cyclicsymmetry));

void FORTRAN(createmdelem,(ITG *imdnode,ITG *nmdnode,
             ITG *ikmpc,ITG *ilmpc,ITG *ipompc,ITG *nodempc,ITG *nmpc,
             ITG *imddof,ITG *nmddof,ITG *nactdof,ITG *mi,ITG *imdmpc,
             ITG *nmdmpc,ITG *imdboun,ITG *nmdboun,ITG *ikboun,ITG *nboun,
             ITG *ilboun,ITG *ithermal,ITG *imdelem,ITG *nmdelem,
             ITG *iponoel,ITG *inoel,char *prlab,char *prset,ITG *nprint,
             char *lakon,char *set,ITG *nset,ITG *ialset,ITG *ipkon,
             ITG *kon,ITG *istartset,ITG *iendset,ITG *nforc,
             ITG *ikforc,ITG *ilforc));

void FORTRAN(createnodeneigh,(ITG *nk,ITG *istartnk,
                          ITG *ialnk,ITG *istartnneigh,ITG *ialnneigh,
                          ITG *ichecknodes,char *lakon,ITG *ipkon,ITG *kon,
                          ITG *nkinsetinv,ITG *neielemtot));

void FORTRAN(createtet,(ITG *kontet,ITG *ifatet,ITG *netet,
             ITG *inodfa,ITG *ifreefa,double *planfa,ITG *ipofa,
             ITG *nodes,double *cotet,ITG *iparentelement));

void FORTRAN(createtiedsurfs,(ITG *nodface,ITG *ipoface,char *set,
             ITG *istartset,ITG *iendset,ITG *ialset,char *tieset,
             ITG *inomat,ITG *ne,ITG *ipkon,char *lakon,ITG *kon,
             ITG *ntie,double *tietol,ITG *nalset,ITG *nk,ITG *nset,
             ITG *iactive));

void FORTRAN(create_iau6,(ITG *nef,ITG *ipnei,ITG *neiel,ITG *jq,
                          ITG *irow,ITG *nzs,ITG *iau6,char *lakonf));

void FORTRAN(dattime,(char *date,char *clock));

void CEE(ddotc,(ITG *n,double *dx,ITG *incx,double *dy,ITG *incy,
                double *funcddot));

void dam1parll(ITG *mt,ITG *nactdof,double *aux2,double *v,
                    double *vini,ITG *nk,ITG *num_cpus);

void *dam1parllmt(ITG *i);

void dam2parll(double *dampwk,double *cv,double *cvini,
                    double *aux2,ITG *neq0,ITG *num_cpus);

void *dam2parllmt(ITG *i);

void *ddotc1mt(ITG *i);

void FORTRAN(desiperelem,(ITG *ndesi,ITG *istartdesi,ITG *ialdesi,
                          ITG *ipoeldi,ITG *ieldi,ITG *ne,
                          ITG *istartelem,ITG *ialelem));

void FORTRAN(determineextern,(ITG *ifac,ITG *itetfa,ITG *iedg,ITG *ipoed,
                              ITG *iexternedg,ITG *iexternfa,ITG *iexternnode,
                              ITG *nktet_,ITG *ipofa));

void dfdbj(double *bcont,double **dbcontp,ITG *neq,ITG *nope,
           ITG *konl,ITG *nactdof,double *s,double *z,ITG *ikmpc,
           ITG *ilmpc,ITG *ipompc,ITG *nodempc,ITG *nmpc,
           double *coefmpc,double *fnl,ITG *nev,
           ITG **ikactcontp,ITG **ilactcontp,ITG *nactcont,ITG *nactcont_,
           ITG *mi,ITG *cyclicsymmetry,ITG *izdof,ITG *nzdof);
      
void FORTRAN(dgesv,(ITG *nteq,ITG *nhrs,double *ac,ITG *lda,ITG *ipiv,
                     double *bc,ITG *ldb,ITG *info)); 

void FORTRAN(dgetrs,(char *trans,ITG *nteq,ITG *nrhs,double *ac,ITG *lda,
                      ITG *ipiv,double *bc,ITG *ldb,ITG *info));

void FORTRAN(dgmres1,(ITG *n,double *b,double *x,ITG *nelt,ITG *ia,ITG *ja,
             double *a,ITG *isym,ITG *itol,double *tol,ITG *itmax,ITG *iter,
             double *err,ITG *ierr,ITG *iunit,double *sb,double *sx,
             double *rgwk,ITG *lrgw,ITG *igwk,ITG *ligw,double *rwork,
             ITG *iwork));

void dgmresmain(ITG *n,double *b,double *x,ITG *nelt,ITG *ia,ITG *ja,
             double *a,ITG *isym,ITG *itol,double *tol,ITG *itmax,ITG *iter,
             double *err,ITG *ierr,ITG *iunit,double *sb,double *sx,
             double *rgwk,ITG *lrgw,ITG *igwk,ITG *ligw,double *rwork,
             ITG *iwork,ITG *nestart,ITG *num_cpus);

void *dgmres1mt(ITG *i);

void FORTRAN(disp_sen_dv,(ITG *nodeset,ITG *istartset,ITG *iendset,ITG *ialset,
                          ITG *iobject,ITG *mi,ITG *nactdof,double *dgdu,
                          double *vold,char *objectset,ITG *nactdofinv,
                          ITG *neq,double *g0));

void FORTRAN(distributesens,(ITG *istartdesi,ITG *ialdesi,ITG *ipkon,
            char *lakon,ITG *ipoface,ITG *ndesi,ITG *nodedesi,ITG *nodface,
            ITG *kon,double *co,double *dgdx,ITG *nobject,
            double *weightformgrad,ITG *nodedesiinv,ITG *noregion,
            char *objectset,double *dgdxglob,ITG *nk,double *physcon));

void divparll(double *var1,double *var2,ITG *size,ITG *num_cpus);

void *divparllmt(ITG *i);

void FORTRAN(drfftf,(ITG *ndata,double *r,double *wsave,ITG *isave));

void FORTRAN(drffti,(ITG *ndata,double *wsave,ITG *isave));

void FORTRAN(dnaupd,(ITG *ido,char *bmat,ITG *n,char *which,ITG *nev,
             double *tol,double *resid,ITG *ncv,double *z,ITG *ldz,
             ITG *iparam,ITG *ipntr,double *workd,double *workl,
             ITG *lworkl,ITG *info));

void FORTRAN(dsaupd,(ITG *ido,char *bmat,ITG *n,char *which,ITG *nev,
             double *tol,double *resid,ITG *ncv,double *z,ITG *ldz,
             ITG *iparam,ITG *ipntr,double *workd,double *workl,
             ITG *lworkl,ITG *info));

void FORTRAN(dneupd,(ITG *rvec,char *howmny,ITG *select,double *d,
             double *di,double *z,ITG *ldz,double *sigma,double *sigmai,
             double *workev,char *bmat,ITG *neq,char *which,
             ITG *nev,double *tol,double *resid,ITG *ncv,double *v,
             ITG *ldv,ITG *iparam,ITG *ipntr,double *workd,
             double *workl,ITG *lworkl,ITG *info));

void FORTRAN(dseupd,(ITG *rvec,char *howmny,ITG *select,double *d,double *z,
             ITG *ldz,double *sigma,char *bmat,ITG *neq,char *which,
             ITG *nev,double *tol,double *resid,ITG *ncv,double *v,
             ITG *ldv,ITG *iparam,ITG *ipntr,double *workd,
             double *workl,ITG *lworkl,ITG *info));

void FORTRAN(dslugm,(ITG *nef,double *b,double *vel,ITG *nelt,ITG *ia,
		     ITG *ja,double *aua,ITG *isym,ITG *nsave,ITG *itol,
		     double *tol,ITG *itmax,ITG *iter,double *err,
		     ITG *ierr,ITG *iunit,double *rgwk,ITG *lenw,
		     ITG *igwk,ITG *leniw));

void FORTRAN(dsort,(double *dx,ITG *iy,ITG *n,ITG *kflag));

void dyna(double **cop,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,ITG *ne,
               ITG **nodebounp,ITG **ndirbounp,double **xbounp,ITG *nboun,
               ITG **ipompcp,ITG **nodempcp,double **coefmpcp,char **labmpcp,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,double *xforc,
               ITG *nforc,ITG *nelemload,char *sideload,double *xload,
               ITG *nload,
               ITG **nactdofp,ITG *neq,ITG *nzl,ITG *icol,ITG *irow,
               ITG *nmethod,ITG **ikmpcp,ITG **ilmpcp,ITG **ikbounp,
               ITG **ilbounp,double *elcon,ITG *nelcon,double *rhcon,
               ITG *nrhcon,double *cocon,ITG *ncocon,
               double *alcon,ITG *nalcon,double *alzero,
               ITG **ielmatp,ITG **ielorienp,ITG *norien,double *orab,
               ITG *ntmat_,double **t0p,
               double **t1p,ITG *ithermal,double *prestr,ITG *iprestr,
               double **voldp,ITG *iperturb,double **stip,ITG *nzs,
               double *timepar,double *xmodal,
               double **veoldp,char *amname,double *amta,
               ITG *namta,ITG *nam,ITG *iamforc,ITG *iamload,
               ITG **iamt1p,ITG *jout,
               ITG *kode,char *filab,double **emep,double *xforcold,
               double *xloadold,
               double **t1oldp,ITG **iambounp,double **xbounoldp,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double **xstatep,ITG *npmat_,char *matname,ITG *mi,
               ITG *ncmat_,ITG *nstate_,double **enerp,char *jobnamec,
               double *ttime,char *set,ITG *nset,ITG *istartset,
               ITG *iendset,ITG **ialsetp,ITG *nprint,char *prlab,
               char *prset,ITG *nener,double *trab,
               ITG **inotrp,ITG *ntrans,double **fmpcp,char *cbody,ITG *ibody,
               double *xbody,ITG *nbody,double *xbodyold,ITG *istep,
               ITG *isolver,ITG *jq,char *output,ITG *mcs,ITG *nkon,
               ITG *mpcend,ITG *ics,double *cs,ITG *ntie,char *tieset,
               ITG *idrct,ITG *jmax,
               double *ctrl,ITG *itpamp,double *tietol,ITG *nalset,
               ITG *ikforc,ITG *ilforc,double *thicke,
               ITG *nslavs,ITG *nmat,char *typeboun,ITG *ielprop,double *prop,
               char *orname);

void dynacont(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
              ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
              ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
              ITG *nmpc,ITG *nodeforc,ITG *ndirforc,double *xforc,
              ITG *nforc,ITG *nelemload,char *sideload,double *xload,
              ITG *nload,
              ITG *nactdof,ITG *neq,ITG *nzl,ITG *icol,ITG *irow,
              ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
              ITG *ilboun,double *elcon,ITG *nelcon,double *rhcon,
              ITG *nrhcon,double *cocon,ITG *ncocon,
              double *alcon,ITG *nalcon,double *alzero,
              ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,
              ITG *ntmat_,double *t0,
              double *t1,ITG *ithermal,double *prestr,ITG *iprestr,
              double *vold,ITG *iperturb,double *sti,ITG *nzs,
              double *tinc,double *tper,double *xmodal,
              double *veold,char *amname,double *amta,
              ITG *namta,ITG *nam,ITG *iamforc,ITG *iamload,
              ITG *iamt1,ITG *jout,char *filab,double *eme,double *xforcold,
              double *xloadold,
              double *t1old,ITG *iamboun,double *xbounold,ITG *iexpl,
              double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
              double *xstate,ITG *npmat_,char *matname,ITG *mi,
              ITG *ncmat_,ITG *nstate_,double *ener,char *jobnamec,
              double *ttime,char *set,ITG *nset,ITG *istartset,
              ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
              char *prset,ITG *nener,double *trab,
              ITG *inotr,ITG *ntrans,double *fmpc,char *cbody,ITG *ibody,
              double *xbody,ITG *nbody,double *xbodyold,ITG *istep,
              ITG *isolver,ITG *jq,char *output,ITG *mcs,ITG *nkon,
              ITG *mpcend,ITG *ics,double *cs,ITG *ntie,char *tieset,
              ITG *idrct,ITG *jmax,double *tmin,double *tmax,
              double *ctrl,ITG *itpamp,double *tietol,ITG *iit,
              ITG *ncont,ITG *ne0,double *reltime,double *dtime,
              double *bcontini,double *bj,double *aux,ITG *iaux,
              double *bcont,ITG *nev,double *v,
              ITG *nkon0,double *deltmx,double *dtheta,double *theta,
              ITG *iprescribedboundary,ITG *mpcfree,ITG *memmpc_,
              ITG *itietri,ITG *koncont,double *cg,double *straight,
              ITG *iinc,double *vini,
              double *aa,double *bb,double *aanew,double *d,double *z,
              double *zeta,double *b,double *time0,double *time1,
              ITG *ipobody,
              double *xforcact,double *xloadact,double *t1act,
              double *xbounact,double *xbodyact,double *cd,double *cv,
              double *ampli,double *dthetaref,double *bjp,double *bp,
              double *cstr,ITG *imddof,
              ITG *nmddof,ITG **ikactcontp,ITG *nactcont,ITG *nactcont_,
              double *aamech,double *bprev,ITG *iprev,ITG *inonlinmpc,
              ITG **ikactmechp,ITG *nactmech,ITG *imdnode,ITG *nmdnode,
              ITG *imdboun,ITG *nmdboun,ITG *imdmpc,ITG *nmdmpc,
              ITG *itp,ITG *inext,
              ITG *imastop,ITG *nslavnode,ITG *islavnode,
              ITG *islavsurf,
              ITG *itiefac,double *areaslav,ITG *iponoels,ITG *inoels,
              double *springarea,ITG *izdof,ITG *nzdof,double *fn,
              ITG *imastnode,ITG *nmastnode,double *xmastnor,
              double *xstateini,ITG *nslavs,
              ITG *cyclicsymmetry,double *xnoels,ITG *ielas,ITG *ielprop,
              double *prop);
 
void dynboun(double *amta,ITG *namta,ITG *nam,double *ampli,double *time,
             double *ttime,double *dtime,double *xbounold,double *xboun,
             double *xbounact,ITG *iamboun,ITG *nboun,ITG *nodeboun,
             ITG *ndirboun,double *ad,double *au,double *adb,
             double *aub,ITG *icol,ITG *irow,ITG *neq,ITG *nzs,
             double *sigma,double *b,ITG *isolver,
             double *alpham,double *betam,ITG *nzl,
             ITG *init,double *bact,double *bmin,ITG *jq,char *amname,
             double *bv,double *bprev,double *bdiff,
             ITG *nactmech,ITG *icorrect,ITG *iprev,double *reltime);

void FORTRAN(dynresults,(ITG *nk,double *v,ITG *ithermal,ITG *nactdof,
             double *vold,ITG *nodeboun,ITG *ndirboun,double *xboun,
             ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
             char *labmpc,ITG *nmpc,double *b,double *bp,double *veold,
             double *dtime,ITG *mi,ITG *imdnode,ITG *nmdnode,ITG *imdboun,
             ITG *nmdboun,ITG *imdmpc,ITG *nmdmpc,ITG *nmethod,double *time));

void FORTRAN(edgedivide,(ITG *nnewnodes,ITG *nktet_,ITG *ipoed,
                         ITG *iexternedg,ITG *iedg,double *d,
                         double *h,ITG *n,double *r,ITG *iext));

void FORTRAN(effectivemodalmass,(ITG *neq,ITG *nactdof,ITG *mi,double *adb,
                        double *aub,ITG *jq,ITG *irow,ITG *nev,double *z,
                        double *co,ITG *nk));

void electromagnetics(double **co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
             ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
             ITG **ipompcp,ITG **nodempcp,double **coefmpcp,char **labmpcp,
             ITG *nmpc,ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
             ITG **nelemloadp,char **sideloadp,double *xload,
             ITG *nload,ITG *nactdof,ITG **icolp,ITG **jqp,ITG **irowp,
             ITG *neq,ITG *nzl,ITG *nmethod,ITG **ikmpcp,ITG **ilmpcp,
             ITG *ikboun,ITG *ilboun,double *elcon,ITG *nelcon,
             double *rhcon,ITG *nrhcon,double *alcon,ITG *nalcon,
             double *alzero,ITG **ielmatp,ITG **ielorienp,ITG *norien,
             double *orab,ITG *ntmat_,double *t0,double *t1,double *t1old,
             ITG *ithermal,double *prestr,ITG *iprestr,double **vold,
             ITG *iperturb,double *sti,ITG *nzs,ITG *kode,char *filab,
             ITG *idrct,ITG *jmax,ITG *jout,double *timepar,double *eme,
             double *xbounold,double *xforcold,double *xloadold,
             double *veold,double *accold,char *amname,double *amta,
             ITG *namta,ITG *nam,ITG *iamforc,ITG **iamloadp,ITG *iamt1,
             double *alpha,ITG *iexpl,ITG *iamboun,double *plicon,
             ITG *nplicon,double *plkcon,ITG *nplkcon,
             double **xstatep,ITG *npmat_,ITG *istep,double *ttime,
             char *matname,double *qaold,ITG *mi,
             ITG *isolver,ITG *ncmat_,ITG *nstate_,ITG *iumat,
             double *cs,ITG *mcs,ITG *nkon,double **ener,ITG *mpcinfo,
             char *output,double *shcon,ITG *nshcon,double *cocon,ITG *ncocon,
             double *physcon,ITG *nflow,double *ctrl,
             char **setp,ITG *nset,ITG **istartsetp,
             ITG **iendsetp,ITG **ialsetp,ITG *nprint,char *prlab,
             char *prset,ITG *nener,ITG *ikforc,ITG *ilforc,double *trab,
             ITG *inotr,ITG *ntrans,double **fmpcp,char *cbody,
             ITG *ibody,double *xbody,ITG *nbody,double *xbodyold,
             ITG *ielprop,double *prop,ITG *ntie,char **tiesetp,
             ITG *itpamp,ITG *iviewfile,char *jobnamec,double **tietolp,
             ITG *nslavs,double *thicke,ITG *ics,ITG *nalset,ITG *nmpc_,
             ITG *nmat,char *typeboun,ITG *iaxial,ITG *nload_,ITG *nprop,
             ITG *network,char *orname);

void elementcpuload(ITG *neapar,ITG *nebpar,ITG *ne,ITG *ipkon,
                             ITG *num_cpus);

void FORTRAN(elementpernode,(ITG *iponoel,ITG *inoel,char *lakon,ITG *ipkon,
                             ITG *kon,ITG *ne));

void FORTRAN(elemperorien,(ITG *ipoorel,ITG *iorel,ITG *ielorien,
                              ITG *ne,ITG *mi));

void FORTRAN(elemperdesi,(ITG *ndesi,ITG *nodedesi,ITG *iponoel,
                            ITG *inoel,ITG *istartdesi,ITG *ialdesi,
                            char *lakon,ITG *ipkon,ITG *kon,ITG *nodedesiinv,
                            ITG *icoordinate,ITG *noregion));

void FORTRAN(envtemp,(ITG *itg,ITG *ieg,ITG *ntg,ITG *ntr,char *sideload,
                      ITG *nelemload,ITG *ipkon,ITG *kon,char *lakon,
                      ITG *ielmat,ITG *ne,ITG *nload,
                      ITG *kontri,ITG *ntri,ITG *nloadtr,
                      ITG *nflow,ITG *ndirboun,ITG *nactdog,
                      ITG *nodeboun,ITG *nacteq,
                      ITG *nboun,ITG *ielprop,double *prop,ITG *nteq,
                      double *v,ITG *network,double *physcon,
                      double *shcon,ITG *ntmat_,double *co,
                      double *vold,char *set,ITG *nshcon,
                      double *rhcon,ITG *nrhcon,ITG *mi,ITG *nmpc,
                      ITG *nodempc,ITG *ipompc,char *labmpc,ITG *ikboun,
                      ITG *nasym,double *ttime,double *time,ITG *iaxial));

void FORTRAN(equationcheck,(double *ac,ITG *nteq,ITG *nactdog,
                            ITG *itg,ITG *ntg,ITG *nacteq,ITG *network));

void FORTRAN(errorestimator,(double *yi,double *yn,ITG *ipkon,
             ITG *kon,char *lakon,ITG *nk,ITG *ne,ITG *mi,ITG *ielmat,
             ITG *nterms,ITG *inum,double *co,double *vold,char *cflag,
             ITG *ielprop,double *prop));

void expand(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
             ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
             ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
             ITG *nmpc,ITG *nodeforc,ITG *ndirforc,double *xforc,
             ITG *nforc,ITG *nelemload,char *sideload,double *xload,
             ITG *nload,ITG *nactdof,ITG *neq,
             ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,ITG *ilboun,
             double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
             double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
             ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
             double *t0,ITG *ithermal,double *prestr,ITG *iprestr,
             double *vold,ITG *iperturb,double *sti,ITG *nzs,
             double *adb,double *aub,char *filab,double *eme,
             double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
             double *xstate,ITG *npmat_,char *matname,ITG *mi,
             ITG *ics,double *cs,ITG *mpcend,ITG *ncmat_,
             ITG *nstate_,ITG *mcs,ITG *nkon,double *ener,
             char *jobnamec,char *output,char *set,ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,ITG *nener,double *trab,
             ITG *inotr,ITG *ntrans,double *ttime,double *fmpc,
             ITG *nev,double **z,ITG *iamboun,double *xbounold,
             ITG *nsectors,ITG *nm,ITG *icol,ITG *irow,ITG *nzl,ITG *nam,
             ITG *ipompcold,ITG *nodempcold,double *coefmpcold,
             char *labmpcold,ITG *nmpcold,double *xloadold,ITG *iamload,
             double *t1old,double *t1,ITG *iamt1,double *xstiff,ITG **icolep,
             ITG **jqep,ITG **irowep,ITG *isolver,
             ITG *nzse,double **adbep,double **aubep,ITG *iexpl,ITG *ibody,
             double *xbody,ITG *nbody,double *cocon,ITG *ncocon,
             char* tieset,ITG* ntie,ITG *imddof,ITG *nmddof,
             ITG *imdnode,ITG *nmdnode,ITG *imdboun,ITG *nmdboun,
             ITG *imdmpc,ITG *nmdmpc,ITG **izdofp,ITG *nzdof,ITG *nherm,
             double *xmr,double *xmi,char *typeboun,ITG *ielprop,double *prop,
             char *orname,ITG *itiefac);

void FORTRAN(extfacepernode,(ITG *iponoelfa,ITG *inoelfa,char *lakonfa,
                          ITG *ipkonfa,ITG *konfa,ITG *nsurfs,ITG *inoelsize));

void FORTRAN(extrapol2dto3d,(double *dgdxglob,ITG *iponod2dto3d,ITG *ndesi,
                             ITG *nodedesi,ITG *nobject,ITG *nk,
                             double *xinterpol,ITG *nnodes,ITG *ipkon,
                             char *lakon,ITG *kon,ITG *ne,ITG *iponoel,
                             ITG *inoel));

void FORTRAN(extrapolate,(double *yi,double *yn,ITG *ipkon,ITG *inum,
             ITG *kon,char *lakon,ITG *nfield,ITG *nk,ITG *ne,ITG *mi,
             ITG *ndim,double *orab,ITG *ielorien,double *co,ITG *iorienloc,
             char *cflag,double *vold,ITG *force,
             ITG *ielmat,double *thicke,ITG *ielprop,double *prop));

void FORTRAN(extrapolatefem,(double *sti,double *stn,int *ipkon,int *inum,
             int *kon,char *lakon,int *nfield,int *nk,int *ne,int *mi,
             int *ndim,double *orab,int *ielorien,double *co,int *iorienglob,
             char *cflag,int *nelemload,int *nload,int *nodeboun,int *nboun,
             int *ndirboun,double *vold,int *ithermal,int *force,
             int *cfd,int *ielmat,double *thicke,char *filab));

void FORTRAN(extrapolate_d_v_simple,(ITG *ielfa,double *xrlfa,
                              double *adv,double *advfa,double *hfa,
                              ITG *icyclic,double *c,ITG *ifatie,
                              double *vel,ITG *nef,double *volume,
                              ITG *nfacea,ITG *nfaceb,ITG *ncfd));

void *extrapolate_d_v_simple1mt(ITG *i);

void extrapolate_d_v_simplemain(ITG *nface,ITG *ielfa,double *xrlfa,
                                double *adv,double *advfa,double *hfa,
                                ITG *icyclic,double *c,ITG *ifatie,
                                double *vel,ITG *nef,double *volume,
                                ITG *num_cpus,ITG *ncfd);

void FORTRAN(extrapolate_d_v_simplec,(ITG *ielfa,double *xrlfa,
                              double *adv,double *advfa,double *hfa,
                              ITG *icyclic,double *c,ITG *ifatie,
                              double *vel,ITG *nef,double *volume,
                              double *auv,ITG *ipnei,ITG *nfacea,
                              ITG *nfaceb,ITG *ncfd));

void *extrapolate_d_v_simplec1mt(ITG *i);

void extrapolate_d_v_simplecmain(ITG *nface,ITG *ielfa,double *xrlfa,
                              double *adv,double *advfa,double *hfa,
                              ITG *icyclic,double *c,ITG *ifatie,
                              double *vel,ITG *nef,double *volume,
                              double *auv,ITG *ipnei,ITG *num_cpus,
                              ITG *ncfd);

void FORTRAN(extrapolate_se,(double *yi,double *yn,ITG *ipkon,ITG *inum,
             ITG *kon,char *lakon,ITG *nfield,ITG *nk,ITG *ne,ITG *mi,
             ITG *ndim,double *orab,ITG *ielorien,double *co,ITG *iorienloc,
             char *cflag,double *vold,ITG *force,
             ITG *ielmat,double *thicke,ITG *ielprop,double *prop,
             ITG *ialeneigh,ITG *neaneigh,ITG *nebneigh,ITG *ialnneigh,
             ITG *naneigh,ITG *nbneigh));

void FORTRAN(extrapolatefluid,(ITG *nk,ITG *iponofa,ITG *inofa,ITG *inum,
             double *vfa,double *v,ITG *ielfa,ITG *ithermal,
             ITG *imach,ITG *ikappa,double *xmach,double *xkappa,
             double *shcon,ITG *nshcon,ITG *ntmat_,ITG *ielmatf,
             double *physcon,ITG *mi,ITG *iturb,double *xturb,
             double *gradtfa,double *gradvfa,double *gradpfa,
             double *gradkfa,double *gradofa,double *co,double *cofa,
             ITG *ifabou));

void FORTRAN(extrapol_dpel,(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
             double *vfa,ITG *ifabou,double *xbounact,ITG *nef,
             double *gradpcel,double *gradpcfa,ITG *neifa,double *rf,
             double *area,double *volume,double *xle,double *xxi,
             ITG *icyclic,double *xxn,ITG *ipnei,ITG *ifatie,
             double *coefmpc,ITG *nmpc,char *labmpc,ITG *ipompc,
             ITG *nodempc,ITG *ifaext,ITG *nfaext,ITG *nactdoh,
             ITG *iflag,double *xxj,double *xlet,double *c));

void extrapol_dpelmain(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
             double *vfa,ITG *ifabou,double *xbounact,ITG *nef,
             double *gradpcel,double *gradpcfa,ITG *neifa,double *rf,
             double *area,double *volume,double *xle,double *xxi,
             ITG *icyclic,double *xxn,ITG *ipnei,ITG *ifatie,
             double *coefmpc,ITG *nmpc,char *labmpc,ITG *ipompc,
             ITG *nodempc,ITG *ifaext,ITG *nfaext,ITG *nactdoh,
             ITG *iflag,double *xxj,double *xlet,double *c,
             ITG *num_cpus,double *xxna,ITG *ncfd,ITG *ip);

void FORTRAN(extrapol_dpel1,(ITG *ielfa1,double *xrlfa1,double *vfap1,
                            double *vel1,ITG *ifabou1,
                             ITG *nef1,ITG *nfacea,ITG *nfaceb));

void *extrapol_dpel1mt(ITG *i);

void FORTRAN(extrapol_dpel2,(ITG *ipnei1,ITG *neifa1,double *vfap1,
                             double *area1,double *xxna1,double *volume1,
                             double *gradpcel1,ITG *nefa,ITG *nefb,
                             ITG *ncfd1));

void *extrapol_dpel2mt(ITG *i);

void FORTRAN(extrapol_dpel3,(ITG *ielfa1,double *xrlfa1,ITG *icyclic1,
                             ITG *ifatie1,double *gradpfa1,double *gradpel1,
                             double *c1,ITG *ipnei1,double *xxi1,ITG *nfacea,
                             ITG *nfaceb,ITG *ncfd1));

void *extrapol_dpel3mt(ITG *i);

void FORTRAN(extrapol_dpel5,(ITG *ielfa1,ITG *ipnei1,double *vel1,
                             double *xlet1,double *gradpcfa1,double *xxj1,
                             ITG *nef1,ITG *nfacea,ITG *nfaceb,ITG *ncfd1));

void *extrapol_dpel5mt(ITG *i);

void FORTRAN(extrapol_kel1,(ITG *ielfa1,double *xrlfa1,double *vfap1,
                            double *vel1,ITG *ifabou1,
                            ITG *nef1,double *umfa1,double *constant1,
                            ITG *inlet,ITG *nfacea,ITG *nfaceb));

void *extrapol_kel1mt(ITG *i);

void FORTRAN(extrapol_kel2,(ITG *ipnei1,ITG *neifa1,double *vfa1,
                            double *area1,double *xxna1,double *volume1,
                            double *gradkel1,ITG *nefa,ITG *nefb,
                            ITG *ncfd1));

void *extrapol_kel2mt(ITG *i);

void FORTRAN(extrapol_kel3,(ITG *ielfa1,double *xrlfa1,ITG *icyclic1,
                            ITG *ifatie1,double *gradkfa1,double *gradkel1,
                            double *c1,ITG *ipnei1,double *xxi1,
                            ITG *nfacea,ITG *nfaceb,ITG *ncfd1));

void *extrapol_kel3mt(ITG *i);

void FORTRAN(extrapol_kel4,(ITG *ielfa1,double *vfa1,double *vfap1,
                            double *gradkfa1,double *rf1,ITG *ifabou1,
                            ITG *ipnei1,double *vel1,
                            double *xxi1,double *xle1,ITG *nef1,
                            ITG *inlet1,ITG *nfacea,ITG *nfaceb));

void *extrapol_kel4mt(ITG *i);

void FORTRAN(extrapol_kel5,(ITG *ielfa1,ITG *ipnei1,double *vel1,
                            double *xlet1,double *gradkfa1,double *xxj1,
                            ITG *nef1,ITG *nfacea,ITG *nfaceb,ITG *ncfd1));

void *extrapol_kel5mt(ITG *i);

void extrapol_kelmain(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
             double *vfa,ITG *ifabou,double *xbounact,ITG *nef,
             double *gradkel,double *gradkfa,ITG *neifa,double *rf,
             double *area,double *volume,double *xle,double *xxi,
             ITG *icyclic,double *xxn,ITG *ipnei,ITG *ifatie,
             double *xlet,double *xxj,double *coefmpc,
             ITG *nmpc,char *labmpc,ITG *ipompc,ITG *nodempc,ITG *ifaext,
             ITG *nfaext,ITG *nactdoh,double *umfa,double *physcon,
             ITG *iitg,double *c,ITG *num_cpus,ITG *compressible,
             double *xxna,ITG *ncfd,ITG *inlet);

void FORTRAN(extrapol_oel1,(ITG *ielfa1,double *xrlfa1,double *vfap1,
                            double *vel1,ITG *ifabou1,
                            ITG *nef1,double *umfa1,double *constant1,
                            double *dy1,double *vfa1,ITG *inlet1,
                            ITG *nfacea,ITG *nfaceb));

void *extrapol_oel1mt(ITG *i);

void FORTRAN(extrapol_oel2,(ITG *ipnei1,ITG *neifa1,double *vfa1,
                            double *area1,double *xxna1,double *volume1,
                            double *gradoel1,ITG *nefa,ITG *nefb,
                            ITG *ncfd1));

void *extrapol_oel2mt(ITG *i);

void FORTRAN(extrapol_oel3,(ITG *ielfa1,double *xrlfa1,ITG *icyclic1,
                            ITG *ifatie1,double *gradofa1,double *gradoel1,
                            double *c1,ITG *ipnei1,double *xxi1,
                            ITG *nfacea,ITG *nfaceb,ITG *ncfd1));

void *extrapol_oel3mt(ITG *i);

void FORTRAN(extrapol_oel4,(ITG *ielfa1,double *vfa1,double *vfap1,
                            double *gradofa1,double *rf1,ITG *ifabou1,
                            ITG *ipnei1,double *vel1,
                            double *xxi1,double *xle1,ITG *nef1,
                            ITG *inlet1,ITG *nfacea,ITG *nfaceb));

void *extrapol_oel4mt(ITG *i);

void FORTRAN(extrapol_oel5,(ITG *ielfa1,ITG *ipnei1,double *vel1,
                            double *xlet1,double *gradofa1,double *xxj1,
                            ITG *nef1,ITG *nfacea,ITG *nfaceb,ITG *ncfd1));

void *extrapol_oel5mt(ITG *i);

void extrapol_oelmain(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
             double *vfa,ITG *ifabou,double *xbounact,ITG *nef,
             double *gradoel,double *gradofa,ITG *neifa,double *rf,
             double *area,double *volume,double *xle,double *xxi,
             ITG *icyclic,double *xxn,ITG *ipnei,ITG *ifatie,
             double *xlet,double *xxj,double *coefmpc,
             ITG *nmpc,char *labmpc,ITG *ipompc,ITG *nodempc,ITG *ifaext,
             ITG *nfaext,ITG *nactdoh,double *umfa,double *physcon,
             ITG *iitg,double *c,double *dy,ITG *num_cpus,
             ITG *compressible,double *xxna,ITG *ncfd,ITG *inlet);

void FORTRAN(extrapol_pel,(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
             double *vfa,ITG *ifabou,double *xbounact,ITG *nef,
             double *gradpel,double *gradpfa,ITG *neifa,double *rf,
             double *area,double *volume,double *xle,double *xxi,
             ITG *icyclic,double *xxn,ITG *ipnei,ITG *ifatie,
             double *coefmpc,ITG *nmpc,char *labmpc,ITG *ipompc,
             ITG *nodempc,ITG *ifaext,ITG *nfaext,ITG *nactdoh,
             ITG *iitg,double *c));

void FORTRAN(extrapol_pel1,(ITG *ielfa1,double *xrlfa1,double *vfap1,
                            double *vel1,ITG *ifabou1,double *xbounact1,
                            ITG *nef1,ITG *nfacea,ITG *nfaceb));

void *extrapol_pel1mt(ITG *i);

void FORTRAN(extrapol_pel2,(ITG *ipnei1,ITG *neifa1,double *vfa1,
                            double *area1,double *xxna1,double *volume1,
                            double *gradpel1,ITG *nefa,ITG *nefb,
                            ITG *ncfd1));

void *extrapol_pel2mt(ITG *i);

void FORTRAN(extrapol_pel3,(ITG *ielfa1,double *xrlfa1,ITG *icyclic1,
                            ITG *ifatie1,double *gradpfa1,double *gradpel1,
                            double *c1,ITG *ipnei1,double *xxi1,ITG *nfacea,
                            ITG *nfaceb,ITG *ncfd1));

void *extrapol_pel3mt(ITG *i);

void FORTRAN(extrapol_pel4,(ITG *ielfa1,double *vfa1,double *vfap1,
                            double *gradpfa1,double *rf1,ITG *ifabou1,
                            ITG *ipnei1,double *vel1,
                            double *xxi1,double *xle1,ITG *nef1,
                            ITG *nfacea,ITG *nfaceb));

void *extrapol_pel4mt(ITG *i);

void extrapol_pelmain(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
             double *vfa,ITG *ifabou,double *xbounact,ITG *nef,
             double *gradpel,double *gradpfa,ITG *neifa,double *rf,
             double *area,double *volume,double *xle,double *xxi,
             ITG *icyclic,double *xxn,ITG *ipnei,ITG *ifatie,
             double *coefmpc,ITG *nmpc,char *labmpc,ITG *ipompc,
             ITG *nodempc,ITG *ifaext,ITG *nfaext,ITG *nactdoh,
             ITG *iitg,double *c,ITG *num_cpus,ITG *compressible,
             double *xxna,ITG *ncfd);

void FORTRAN(extrapol_tel,(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
             double *vfa,ITG *ifabou,double *xbounact,ITG *nef,
             double *gradtel,double *gradtfa,ITG *neifa,double *rf,
             double *area,double *volume,double *xle,double *xxi,
             ITG *icyclic,double *xxn,ITG *ipnei,ITG *ifatie,
             double *xload,double *xlet,double *xxj,double *coefmpc,
             ITG *nmpc,char *labmpc,ITG *ipompc,ITG *nodempc,ITG *ifaext,
             ITG *nfaext,ITG *nactdoh,ITG *iitg,double *c));

void FORTRAN(extrapol_tel1,(ITG *ielfa1,double *xrlfa1,double *vfap1,
                            double *vel1,ITG *ifabou1,double *xbounact1,
                            ITG *nef1,ITG *nfacea,ITG *nfaceb));

void *extrapol_tel1mt(ITG *i);

void FORTRAN(extrapol_tel2,(ITG *ipnei1,ITG *neifa1,double *vfa1,
                            double *area1,double *xxna1,double *volume1,
                            double *gradtel1,ITG *nefa,ITG *nefb,
                            ITG *ncfd1));

void *extrapol_tel2mt(ITG *i);

void FORTRAN(extrapol_tel3,(ITG *ielfa1,double *xrlfa1,ITG *icyclic1,
                            ITG *ifatie1,double *gradtfa1,double *gradtel1,
                            double *c1,ITG *ipnei1,double *xxi1,ITG *ifabou1,
                            double *xxn1,double *xload1,ITG *nfacea,
                            ITG *nfaceb,ITG *ncfd1));

void *extrapol_tel3mt(ITG *i);

void FORTRAN(extrapol_tel4,(ITG *ielfa1,double *vfa1,double *vfap1,
                            double *gradtfa1,double *rf1,ITG *ifabou1,
                            ITG *ipnei1,double *vel1,
                            double *xxi1,double *xle1,ITG *nef1,
                            ITG *nfacea,ITG *nfaceb));

void *extrapol_tel4mt(ITG *i);

void FORTRAN(extrapol_tel5,(ITG *ielfa1,ITG *ipnei1,double *vel1,
                            double *xlet1,double *gradtfa1,double *xxj1,
                            ITG *nef1,ITG *nfacea,ITG *nfaceb,ITG *ncfd1));

void *extrapol_tel5mt(ITG *i);

void extrapol_telmain(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
             double *vfa,ITG *ifabou,double *xbounact,ITG *nef,
             double *gradtel,double *gradtfa,ITG *neifa,double *rf,
             double *area,double *volume,double *xle,double *xxi,
             ITG *icyclic,double *xxn,ITG *ipnei,ITG *ifatie,
             double *xload,double *xlet,double *xxj,double *coefmpc,
             ITG *nmpc,char *labmpc,ITG *ipompc,ITG *nodempc,ITG *ifaext,
             ITG *nfaext,ITG *nactdoh,ITG *iitg,double *c,ITG *num_cpus,
             ITG *compressible,double *xxna,ITG *ncfd);

void FORTRAN(extrapol_vel,(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
             double *vfa,ITG *ifabou,double *xbounact,ITG *ipnei,
             ITG *nef,ITG *icyclic,double *c,ITG *ifatie,double *xxn,
             double *gradvel,double *gradvfa,ITG *neifa,double *rf,
             double *area,double *volume,double *xle,double *xxi,
             double *xxj,double *xlet,double *coefmpc,ITG *nmpc,
             char *labmpc,ITG *ipompc,ITG *nodempc,ITG *ifaext,
             ITG *nfaext,ITG *nactdoh,ITG *iitg));

void FORTRAN(extrapol_vel1,(ITG *ielfa1,double *xrlfa1,ITG *icyclic1,
                            ITG *ifatie1,double *vfap1,double *vel1,
                            double *c1,ITG *ifabou1,double *xbounact1,
                            ITG *ipnei1,double *xxn1,ITG *nef1,
                            ITG *nfacea,ITG *nfaceb,ITG *ncfd1));

void *extrapol_vel1mt(ITG *i);

void FORTRAN(extrapol_vel2,(ITG *ipnei1,ITG *neifa1,double *vfa1,
                            double *area1,double *xxna1,double *volume1,
                            double *gradvel1,ITG *nefa,ITG *nefb,ITG *ncfd1));

void *extrapol_vel2mt(ITG *i);

void FORTRAN(extrapol_vel3,(ITG *ielfa1,double *xrlfa1,ITG *icyclic1,
                            ITG *ifatie1,double *gradvfa1,double *gradvel1,
                            double *c1,ITG *ipnei1,double *xxi1,ITG *nfacea,
                            ITG *nfaceb,ITG *ncfd1));

void *extrapol_vel3mt(ITG *i);

void FORTRAN(extrapol_vel4,(ITG *ielfa1,double *vfa1,double *vfap1,
                            double *gradvfa1,double *rf1,ITG *ifabou1,
                            ITG *ipnei1,double *xxn1,double *vel1,
                            double *xxi1,double *xle1,ITG *nef1,
                            ITG *nfacea,ITG *nfaceb,ITG *ncfd1));

void *extrapol_vel4mt(ITG *i);

void FORTRAN(extrapol_vel5,(ITG *ielfa1,ITG *ipnei1,double *vel1,
                            double *xlet1,double *gradvfa1,double *xxj1,
                            ITG *nef1,ITG *nfacea,ITG *nfaceb,ITG *ncfd1));

void *extrapol_vel5mt(ITG *i);

void extrapol_velmain(ITG *nface,ITG *ielfa,double *xrlfa,double *vel,
                      double *vfa,ITG *ifabou,double *xbounact,ITG *ipnei,
                      ITG *nef,ITG *icyclic,double *c,ITG *ifatie,double *xxn,
                      double *gradvel,double *gradvfa,ITG *neifa,double *rf,
                      double *area,double *volume,double *xle,double *xxi,
                      double *xxj,double *xlet,double *coefmpc,ITG *nmpc,
                      char *labmpc,ITG *ipompc,ITG *nodempc,ITG *ifaext,
                      ITG *nfaext,ITG *nactdoh,ITG *iitg,ITG *num_cpus,
                      ITG *compressible,double *xxna,ITG *ncfd);

void FORTRAN(extrapol_xel,(ITG *nface,ITG *ielfa,double *xrlfa,double *xel,
             double *xfa,ITG *ifabou,double *xbounact,ITG *nef,
             double *gradxel,double *gradxfa,ITG *neifa,double *rf,
             double *area,double *volume,double *xle,double *xxi,
             ITG *icyclic,double *xxn,ITG *ipnei,ITG *ifatie,
             double *xload,double *xlet,double *xxj,double *coefmpc,
             ITG *nmpc,char *labmpc,ITG *ipompc,ITG *nodempc,ITG *ifaext,
             ITG *nfaext,ITG *nactdoh,ITG *iitg));

void FORTRAN(fcrit,(double *time,double *tend,double *aai,double *bbi,
                      double *zetaj,double *dj,double *ddj,
                      double *h1,double *h2,double *h3,double *h4,
                      double *func,double *funcp));

void FORTRAN(fill_neiel,(ITG *nef,ITG *ipnei,ITG *neiel,ITG *neielcp));

void FORTRAN(filter,(double *dgdxglob,ITG *nobject,ITG *nk,ITG *nodedesi,
                     ITG *ndesi,char *objectset,double *xo,double *yo,
                     double *zo,double *x,double *y,double *z,ITG *nx,
                     ITG *ny,ITG *nz,ITG *neighbor,double *r,ITG *ndesia,
                     ITG *ndesib,double *xdesi,double *distmin));

void filtermain(double *co,double *dgdxglob,ITG *nobject,ITG *nk,
                ITG *nodedesi,ITG *ndesi,char *objectset,double *xdesi,
                double *distmin);

void *filtermt(ITG *i);

void FORTRAN(findextsurface,(ITG *nodface,ITG *ipoface,ITG *ne,ITG *ipkon,
                             char *lakon,ITG *kon,ITG *konfa,ITG *ipkonfa,
                             ITG *nk,char *lakonfa,ITG *nsurfs,ITG *ifreemax));

void FORTRAN(findslavcfd,(int *nmpc,char *labmpc,int *ipompc,int *nodempc,
                          int *islav,int *nslav,int *inoslav,int *inomast,
                          int *ics,double *cs,int *imast,int *nmast,double *co,
                          int *inomat,int *nr,int *nz,double *rcs,double *zcs,
                          double *rcs0,double *zcs0,int *ncs));

void FORTRAN(findsurface,(ITG *ipoface,ITG *nodface,ITG *ne,ITG *ipkon,ITG *kon,
                     char *lakon,ITG *ntie,char *tieset));

void FORTRAN(fixnode,(ITG *nobject,ITG *nk,char *set,ITG *nset,ITG *istartset,
                      ITG *iendset,ITG *ialset,ITG *iobject,ITG *nodedesiinv,
                      double *dgdxglob,char *objectset));                       

void FORTRAN (flowoutput,(ITG *itg,ITG *ieg,ITG *ntg,ITG *nteq,
                          double *bc,char *lakon,
                          ITG *ntmat_,double *v,double *shcon,ITG *nshcon,
                          ITG *ipkon,ITG *kon,double *co,ITG *nflow,
                          double *dtime,double *ttime,double *time,
                          ITG *ielmat,double *prop,
                          ITG *ielprop,ITG *nactdog,ITG *nacteq,ITG *iin,
                          double *physcon,double *camt,double *camf,double *camp,
                          double *rhcon,ITG *nrhcon,
                          double *vold,char *jobnamef,char *set,ITG *istartset,
                          ITG *iendset,ITG *ialset,ITG *nset,ITG *mi,
                          ITG *iaxial,ITG *istep,ITG *iit));

void FORTRAN(flowresult,(ITG *ntg,ITG *itg,double *cam,double *vold,
              double *v,
              ITG *nload,char *sideload,ITG *nelemload,
              double *xloadact,ITG *nactdog,ITG *network,ITG *mi,
              ITG *ne,ITG *ipkon,char *lakon,ITG *kon));

void FORTRAN(forcesolve,(double *zc,ITG *nev,double *aa,double *bb,
             double *xx,double *eiga,double *eigb,double *eigxx,
             ITG *iter,double *d,ITG *neq,double *z,ITG *istartnmd,
             ITG *iendnmd,ITG *nmd,ITG *cyclicsymmetry,ITG *neqact,
             ITG *igeneralizedforce));

void forparll(ITG *mt,ITG *nactdof,double *aux2,double *veold,
                    ITG *nk,ITG *num_cpus);

void *forparllmt(ITG *i);

void frd(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne0,
         double *v,double *stn,ITG *inum,ITG *nmethod,ITG *kode,
         char *filab,double *een,double *t1,double *fn,double *time,
         double *epn,ITG *ielmat,char *matname,double *enern,
         double *xstaten,ITG *nstate_,ITG *istep,ITG *iinc,
         ITG *ithermal,double *qfn,ITG *mode,ITG *noddiam,
         double *trab,ITG *inotr,ITG *ntrans,double *orab,
         ITG *ielorien,ITG *norien,char *description,ITG *ipneigh,
         ITG *neigh,ITG *mi,double *stx,double *vr,double *vi,
         double *stnr,double *stni,double *vmax,double *stnmax,
         ITG *ngraph,double *veold,double *ener,ITG *ne,double *cs,
         char *set,ITG *nset,ITG *istartset,ITG *iendset,ITG *ialset,
         double *eenmax,double *fnr,double *fni,double *emn,
         double *thicke,char *jobnamec,char *output,double *qfx,
         double *cdn,ITG *mortar,double *cdnr,double *cdni,ITG *nmat,
         ITG *ielprop,double *prop);

void frdcyc(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,double *v,
            double *stn,ITG *inum,ITG *nmethod,ITG *kode,char *filab,
            double *een,double *t1,double *fn,double *time,double *epn,
            ITG *ielmat,char *matname,double *cs,ITG *mcs,ITG *nkon,
            double *enern,double *xstaten,ITG *nstate_,ITG *istep,
            ITG *iinc,ITG *iperturb,double *ener,ITG *mi,char *output,
            ITG *ithermal,double *qfn,ITG *ialset,ITG *istartset,
            ITG *iendset,double *trab,ITG *inotr,ITG *ntrans,double *orab,
            ITG *ielorien,ITG *norien,double *sti,double *veold,ITG *noddiam,
            char *set,ITG *nset,double *emn,double *thicke,char *jobnamec,
            ITG *ne0,double *cdn,ITG *mortar,ITG *nmat,double *qfx,
            ITG *ielprop,double *prop);

void frd_norm_se(double *co,ITG *nk,double *stn,ITG *inum,ITG *nmethod,
         ITG *kode,char *filab,double *fn,double *time,ITG *nstate_,
         ITG *istep,ITG *iinc,ITG *mode,ITG *noddiam,char *description,
         ITG *mi,ITG *ngraph,ITG *ne,double *cs,char *set,ITG *nset,
         ITG *istartset,ITG *iendset,ITG *ialset,double *thicke,
         char *jobnamec,char *output,double *dgdxtotglob,ITG *numobject,
         char *objectset,double *extnor,ITG *ntrans,double *trab,
         ITG *inotr);

void frd_sen(double *co,ITG *nk,double *dstn,ITG *inum,ITG *nmethod,
         ITG *kode,char *filab,double *time,ITG *nstate_,
         ITG *istep,ITG *iinc,ITG *mode,ITG *noddiam,char *description,
         ITG *mi,ITG *ngraph,ITG *ne,double *cs,char *set,ITG *nset,
         ITG *istartset,ITG *iendset,ITG *ialset,
         char *jobnamec,char *output,double *v,ITG *iobject,
         char *objectset,ITG *ntrans,ITG *inotr,double *trab,
         ITG *idesvar,char *orname,ITG *icoordinate,ITG *inorm,
         ITG *irand,ITG *ishape);

void frd_sen_se(double *co,ITG *nk,double *stn,ITG *inum,ITG *nmethod,
         ITG *kode,char *filab,double *fn,double *time,ITG *nstate_,
         ITG *istep,ITG *iinc,ITG *mode,ITG *noddiam,char *description,
         ITG *mi,ITG *ngraph,ITG *ne,double *cs,char *set,ITG *nset,
         ITG *istartset,ITG *iendset,ITG *ialset,double *thicke,
         char *jobnamec,char *output,double *dgdxglob,ITG *iobject,
         char *objectset);

void frd_orien_se(double *co,ITG *nk,double *stn,ITG *inum,ITG *nmethod,
         ITG *kode,char *filab,double *fn,double *time,ITG *nstate_,
         ITG *istep,ITG *iinc,ITG *mode,ITG *noddiam,char *description,
         ITG *mi,ITG *ngraph,ITG *ne,double *cs,char *set,ITG *nset,
         ITG *istartset,ITG *iendset,ITG *ialset,double *thicke,
         char *jobnamec,char *output,double *dgdxtotglob,ITG *numobject,
         char *objectset,ITG *ntrans,ITG *inotr,double *trab,
         ITG *idesvar,char *orname);

void FORTRAN(frdfluid,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
             ITG *ne,double *v,ITG *kode,double *time,
             ITG *ielmatf,char *matname,char *filab,
             ITG *inum,ITG *ntrans,ITG *inotr,double *trab,ITG *mi,
             ITG *istep,double *stn,double *qfn,ITG *nactdofinv,
             double *xmach,double *xkappa,double *physcon,double *xturb,
             double *coel,double *vel,double *cofa,double *vfa,ITG *nface));

void FORTRAN(frdfluidfem,(double *co,int *nk,int *kon,int *ipkon,char *lakon,
             int *ne,double *v,double *vold,int *kode,double *time,
             int *ielmat,char *matname,int *nnstep,double *vtu,
             double *vcontu,double *voldaux,double *physcon,char *filab,
             int *inomat,int *ntrans,int *inotr,double *trab,int *mi,
             double *stn,double *qfn,int *istep));

void frdgeneralvector(double *v,ITG *iset,ITG *ntrans,char * filabl,
               ITG *nkcoords,
               ITG *inum,char *m1,ITG *inotr,double *trab,double *co,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *mi,ITG *ngraph,
               FILE *f1,char *output,char *m3);

void frdheader(ITG *icounter,double *oner,double *time,double *pi,
               ITG *noddiam,double *cs,ITG *null,ITG *mode,
               ITG *noutloc,char *description,ITG *kode,ITG *nmethod,
               FILE *f1,char *output,ITG *istep,ITG *iinc);

void FORTRAN(frditeration,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
             ITG *ne,double *v,double *time,ITG *ielmat,char *matname,
             ITG *mi,ITG *istep,ITG *iinc,ITG *ithermal));

void frdselect(double *field1,double *field2,ITG *iset,ITG *nkcoords,ITG *inum,
     char *m1,ITG *istartset,ITG *iendset,ITG *ialset,ITG *ngraph,ITG *ncomp,
     ITG *ifield,ITG *icomp,ITG *nfield,ITG *iselect,char *m2,FILE *f1,
     char *output,char *m3);

void frdset(char *filabl,char *set,ITG *iset,ITG *istartset,ITG *iendset,
            ITG *ialset,ITG *inum,ITG *noutloc,ITG *nout,ITG *nset,
            ITG *noutmin,ITG *noutplus,ITG *iselect,ITG *ngraph);

void frdvector(double *v,ITG *iset,ITG *ntrans,char * filabl,ITG *nkcoords,
               ITG *inum,char *m1,ITG *inotr,double *trab,double *co,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *mi,ITG *ngraph,
               FILE *f1,char *output,char *m3);

void FORTRAN(frictionheating,(ITG *ne0,ITG *ne,ITG *ipkon,char *lakon,ITG *ielmat,
                     ITG *mi,double *elcon,ITG *ncmat_,ITG *ntmat_,
                     ITG *kon,ITG *islavsurf,double *pmastsurf,
                     double *springarea,double *co,double *vold,
                     double *veold,double *pslavsurf,double *xload,
                     ITG *nload,ITG *nload_,ITG *nelemload,ITG *iamload,
                     ITG *idefload,char *sideload,double *stx,ITG *nam,
                     double *time,double *ttime,char *matname,ITG *istep,
                     ITG *iinc));

void FORTRAN(fsub,(double *time,double *tend,double *aai,double *bbi,
                   double *ddj,double *h1,double *h2,double *h3,double *h4,
                   double *func,double *funcp));

void FORTRAN(fsuper,(double *time,double *tend,double *aai,double *bbi,
                       double *h1,double *h2,double *h3,double *h4,
                       double *h5,double *h6,double *func,double *funcp));

void FORTRAN(gasmechbc,(double *vold,ITG *nload,char *sideload,
                        ITG *nelemload,double *xload,ITG *mi));

void FORTRAN(genadvecelem,(ITG *inodesd,ITG *ipkon,ITG *ne,char *lakon,
             ITG *kon,ITG *nload,char *sideload,ITG *nelemload,ITG *nkon,
             ITG *network));

void FORTRAN(gencontelem_f2f,(char *tieset,ITG *ntie,ITG *itietri,ITG *ne,
             ITG *ipkon,ITG *kon,char *lakon,double *cg,double *straight,
             ITG *ifree,ITG *koncont,double *co,double *vold,double *xo,
             double *yo,double *zo,double *x,double *y,double *z,ITG *nx,
             ITG *ny,ITG *nz,ITG *ielmat,double *elcon,ITG *istep,ITG *iinc,
             ITG *iit,ITG *ncmat_,ITG *ntmat_,ITG *mi,ITG *imastop,
             ITG *islavsurf,ITG *itiefac,double *springarea,double *tietol,
             double *reltime,char *filab,ITG *nasym,
             double *pslavsurf,double *pmastsurf,double *clearini,
             double *theta,double *xstateini,double *xstate,ITG *nstate_,
             ITG *ne0,ITG *icutb,ITG *ialeatoric,ITG *nmethod,
             char *jobnamef,double *alea));

void FORTRAN(gencontelem_n2f,(char *tieset,ITG *ntie,ITG *itietri,ITG *ne,
     ITG *ipkon,ITG *kon,char *lakon,
     double *cg,double *straight,ITG *ifree,ITG *koncont,
     double *co,double *vold,double *xo,double *yo,double *zo,
     double *x,double *y,double *z,ITG *nx,ITG *ny,ITG *nz,
     ITG *ielmat,double *elcon,ITG *istep,ITG *iinc,ITG *iit,
     ITG *ncmat_,ITG *ntmat_,
     ITG *nmethod,ITG *mi,ITG *imastop,ITG *nslavnode,
     ITG *islavnode,ITG *islavsurf,ITG *itiefac,double *areaslav,
     ITG *iponoels,ITG *inoels,double *springarea,
     char *set,ITG *nset,ITG *istartset,ITG *iendset,ITG *ialset,
     double *tietol,double *reltime,
     char* filab,ITG *nasym,double *xnoels,ITG *icutb,ITG *ne0,
     char *jobnamef));

void FORTRAN(gencycsymelemcfd,(double *cs,int *islav,int *nslav,
         int *imast,int *nmast,int *inomat,int *nk,double *co,int *ne,
         int *ipkon,char *lakon,int *kon,int *nkon,int *mi,int *ielmat,
         double *vold,int *ielslav,int *ielmast,int *inoslav,int *inomast,
         int *iponoel,int *inoel));

void FORTRAN(generateeminterfaces,(ITG *istartset,ITG *iendset,
             ITG *ialset,ITG *iactive,ITG *ipkon,char *lakon,ITG *kon,
             ITG *ikmpc,ITG *nmpc,ITG *nafaces));

void  FORTRAN(gennactdofinv,(ITG *nactdof,ITG *nactdofinv,ITG *nk,
       ITG *mi,ITG *nodorig,ITG *ipkon,char *lakon,ITG *kon,ITG *ne));

unsigned long genrand();

void FORTRAN(gentiedmpc,(char *tieset,ITG *ntie,ITG *itietri,
          ITG *ipkon,ITG *kon,char *lakon,char *set,ITG *istartset,
          ITG *iendset,ITG *ialset,double *cg,double *straight,
          ITG *koncont,double *co,double *xo,double *yo,double *zo,
          double *x,double *y,double *z,ITG *nx,ITG *ny,ITG *nz,ITG *nset,
          ITG *ifaceslave,ITG *istartfield,ITG *iendfield,ITG *ifield,
          ITG *ipompc,ITG *nodempc,double *coefmpc,ITG *nmpc,ITG *nmpc_,
          ITG *mpcfree,ITG *ikmpc,ITG *ilmpc,char *labmpc,ITG *ithermal,
          double *tietol,ITG *icfd,ITG *ncont,ITG *imastop,ITG *ikboun,
          ITG *nboun,char *kind));

void FORTRAN(geomview,(double *vold,double *co,double *pmid,double *e1,
             double *e2,double *e3,ITG *kontri,double *area,double *cs,
             ITG *mcs,ITG *inocs,ITG *ntrit,ITG *nk,ITG *mi,double *sidemean));

void FORTRAN(getdesiinfo2d,(char *set,ITG *istartset,ITG *iendset,ITG *ialset,
            ITG *nset,ITG *mi,ITG *nactdof,ITG *ndesi,
            ITG *nodedesi,ITG *ntie,char *tieset,ITG *nodedesiinv,
            char *lakon,ITG *ipkon,ITG *kon,ITG *iponoelfa,ITG *iponod2dto3d,
            ITG *iponor2d,ITG *knor2d,ITG *iponoel2d,ITG *inoel2d,ITG *nobject,
            char *objectset,ITG *iponk2dto3d,ITG *ne));  

void FORTRAN(getdesiinfo3d,(char *set,ITG *istartset,ITG *iendset,ITG *ialset,
            ITG *nset,ITG *mi,ITG *nactdof,ITG *ndesi,
            ITG *nodedesi,ITG *ntie,char *tieset,ITG *itmp,ITG *nmpc,
            ITG *nodempc,ITG *ipompc,ITG *nodedesiinv,ITG *iponoel,
            ITG *inoel,char *lakon,ITG *ipkon,ITG *kon,ITG *noregion,
            ITG *ipoface,ITG *nodface,ITG *nk));  

void getglobalresults (char *jobnamec,ITG **integerglobp,double **doubleglobp,
                       ITG *nboun,ITG *iamboun,double *xboun,ITG *nload,
                       char *sideload,ITG *iamload,ITG *iglob,ITG *nforc,
                       ITG *iamforc,double *xforc,ITG *ithermal,ITG *nk,
                       double *t1,ITG *iamt1,double *sigma);

void getlocalresults(ITG **integerglobp,double **doubleglobp,ITG *nktet,
                     double *cotet,double *h,ITG *netet_,ITG *kontet,
                     ITG *ifatet,double *planfa,ITG *kontetor);

ITG getSystemCPUs();;

void FORTRAN(hrk_ud,(ITG *ielfa,double *vel,
                      ITG *ipnei,ITG *nef,double *flux,
                      double *vfa,ITG *nfacea,ITG *nfaceb));

void *hrk_ud1mt(ITG *i);

void hrk_udmain(ITG *nface,ITG *ielfa,double *vel,
                ITG *ipnei,ITG *nef,double *flux,
                double *vfa,ITG *num_cpus);

void FORTRAN(hro_ud,(ITG *ielfa,double *vel,
                      ITG *ipnei,ITG *nef,double *flux,
                      double *vfa,ITG *nfacea,ITG *nfaceb));

void *hro_ud1mt(ITG *i);

void hro_udmain(ITG *nface,ITG *ielfa,double *vel,
                ITG *ipnei,ITG *nef,double *flux,
                double *vfa,ITG *num_cpus);

void FORTRAN(hrp_mod_smart,(ITG *ielfa,double *vel,double *gradvel,
                             double *xlet,double *xxj,
                             ITG *ipnei,ITG *nef,double *flux,
                             double *vfa,ITG *nfacea,ITG *nfaceb));

void *hrp_mod_smart1mt(ITG *i);

void hrp_mod_smartmain(ITG *nface,ITG *ielfa,double *vel,double *gradvel,
                        double *xlet,double *xxj,
                        ITG *ipnei,ITG *nef,double *flux,
                       double *vfa,ITG *num_cpus);

void FORTRAN(hrp_ud,(ITG *ielfa,double *vel,
                      ITG *ipnei,ITG *nef,double *flux,
                      double *vfa,ITG *nfacea,ITG *nfaceb,
                      double *xxi1,double *xle1,double *gradpel1,
                      ITG *neij1));

void *hrp_ud1mt(ITG *i);

void hrp_udmain(ITG *nface,ITG *ielfa,double *vel,
                ITG *ipnei,ITG *nef,double *flux,
                double *vfa,ITG *num_cpus,double *xxi,
                double *xle,double *gradpel,ITG *neij);

void FORTRAN(hrt_ud,(ITG *ielfa,double *vel,
                      ITG *ipnei,ITG *nef,double *flux,
                      double *vfa,ITG *nfacea,ITG *nfaceb,
                      double *xxi1,double *xle1,double *gradtel1,
                      ITG *neij1));

void *hrt_ud1mt(ITG *i);

void hrt_udmain(ITG *nface,ITG *ielfa,double *vel,
                ITG *ipnei,ITG *nef,double *flux,
                double *vfa,ITG *num_cpus,double *xxi,
                double *xle,double *gradtel,ITG *neij);

void FORTRAN(hrv_mod_smart,(ITG *ielfa,double *vel,double *gradvel,
                             double *xlet,double *xxj,
                             ITG *ipnei,ITG *nef,double *flux,
                             double *vfa,ITG *nfacea,ITG *nfaceb));

void *hrv_mod_smart1mt(ITG *i);

void hrv_mod_smartmain(ITG *nface,ITG *ielfa,double *vel,double *gradvel,
                        double *xlet,double *xxj,
                        ITG *ipnei,ITG *nef,double *flux,
                       double *vfa,ITG *num_cpus);

void FORTRAN(hrv_ud,(ITG *ielfa1,double *vel1,ITG *ipnei1,ITG *nef1,
                      double *flux1,double *vfa1,ITG *nfacea,ITG *nfaceb,
                      double *xxi1,double *xle1,double *gradvel1,ITG *neij1));

void *hrv_ud1mt(ITG *i);

void hrv_udmain(ITG *nface,ITG *ielfa,double *vel,ITG *ipnei,ITG *nef,
                double *flux,double *vfa,ITG *num_cpus,double *xxi,
                double *xle,double *gradvel,ITG *neij);

void FORTRAN(identamta,(double *amta,double *reftime,ITG *istart,ITG *iend,
               ITG *id));

void FORTRAN(identifytiedface,(char *tieset,ITG *ntie,char *set,ITG *nset,
                               ITG *faceslave,char *kind));

void FORTRAN(includefilename,(char *buff,char *includefn,ITG *lincludefn));

void FORTRAN(inicalcbody,(ITG *nef,double *body,ITG *ipobody,ITG *ibody,
                          double *xbody,double *coel,double *vel,char *lakon,
                          ITG *nactdohinv,ITG *icent));

void inicont(ITG* nk,ITG *ncont,ITG *ntie,char *tieset,ITG *nset,char *set,
             ITG *istartset,ITG *iendset,ITG *ialset,ITG **itietrip,
             char *lakon,ITG *ipkon,ITG *kon,ITG **koncontp,
             ITG *ncone,double *tietol,ITG *ismallsliding,ITG **itiefacp,
             ITG **islavsurfp,ITG **islavnodep,ITG **imastnodep,
             ITG **nslavnodep,ITG **nmastnodep,ITG *mortar,
             ITG **imastopp,ITG *nkon,ITG **iponoels,ITG **inoelsp,
             ITG **ipep,ITG **imep,ITG *ne,ITG *ifacecount,
             ITG *iperturb,ITG *ikboun,ITG *nboun,double *co,ITG *istep,
             double **xnoelsp);

void iniparll(ITG *mt,ITG *nactdof,double *b,double *v,
              double *veold,double *accold,double *bet,
              double *gam,double *dtime,double *cam,
              ITG *nk,ITG *num_cpus);

void *iniparllmt(ITG *i);

void FORTRAN(init,(ITG *nktet,ITG *inodfa,ITG *ipofa,ITG *netet_));

void FORTRAN(initialcfd,(ITG *ne,ITG *ipkon,ITG *kon,char *lakon,
             double *co,double *coel,double *cofa,ITG *nface,ITG *ielfa,
             double *area,ITG *ipnei,ITG *neiel,double *xxn,double *xxi,
             double *xle,double *xlen,double *xlet,double *xrlfa,double *cosa,
             double *volume,ITG *neifa,double *xxj,double *cosb,
             double *dmin,ITG *ifatie,double *cs,char *tieset,ITG *icyclic,
             double *c,ITG *neij,double *physcon,ITG *isolidsurf,
             ITG *nsolidsurf,double *dy,double *xxni,double *xxnj,
             double *xxicn,ITG *nflnei,ITG *iturbulent,double *rf,
             double *yy,double *vel,double *velo,double *veloo,
             double *xxna,double *ale,double *alet,double *h));
  
void FORTRAN(initialcfdfem,(double *yy,int *nk,double *co,int *ne,int *ipkon,
       int *kon,char *lakon,double *x,double *y,double *z,double *x0,
       double *y0,double *z0,int *nx,int *ny,int *nz,int *isolidsurf,
       int *neighsolidsurf,double *xsolidsurf,double *dt,int *nshcon,
       double *shcon,int *nrhcon,double *rhcon,double *vold,double *voldaux,
       int *ntmat_,int *iponoel,int *inoel,int *iexplicit,
       int *ielmat,int *nsolidsurf,int *turbulent,double *physcon,
       int *compressible,char *matname,int *inomat,double *vcontu,
       int *mi,int *euler,int *ithermal));

void FORTRAN(initialchannel,(ITG *itg,ITG *ieg,ITG *ntg,double *ac,double *bc,
                         char *lakon,double *v,ITG * ipkon,ITG *kon,
                         ITG *nflow,ITG *ikboun,ITG *nboun,double *prop,
                         ITG *ielprop,ITG *nactdog,ITG *ndirboun,
                         ITG *nodeboun,double *xbounact,ITG *ielmat,
                         ITG *ntmat_,double *shcon,ITG *nshcon,
                         double *physcon,ITG *ipiv,ITG *nteq,
                         double *rhcon,ITG *nrhcon,ITG *ipobody,ITG *ibody,
                         double *xbody,double *co,ITG *nbody,ITG *network,
                         ITG *iin_abs,double *vold,char *set,ITG *istep,
                         ITG *iit,ITG *mi,ITG *ineighe,ITG *ilboun,
                         double *ttime,double *time,ITG *iaxial));

void FORTRAN(initialnet,(ITG *itg,ITG *ieg,ITG *ntg,double *ac,double *bc,
                         char *lakon,double *v,ITG * ipkon,ITG *kon,
                         ITG *nflow,ITG *ikboun,ITG *nboun,double *prop,
                         ITG *ielprop,ITG *nactdog,ITG *ndirboun,
                         ITG *nodeboun,double *xbounact,ITG *ielmat,
                         ITG *ntmat_,double *shcon,ITG *nshcon,
                         double *physcon,ITG *ipiv,ITG *nteq,
                         double *rhcon,ITG *nrhcon,ITG *ipobody,ITG *ibody,
                         double *xbody,double *co,ITG *nbody,ITG *network,
                         ITG *iin_abs,double *vold,char *set,ITG *istep,
                         ITG *iit,ITG *mi,ITG *ineighe,ITG *ilboun,
                         ITG *channel,ITG *iaxial,ITG *nmpc,char *labmpc,
                         ITG *ipompc,ITG *nodempc,double *coefmpc,
                         double *ttime,double *time,ITG *iponoel,ITG *inoel));

void FORTRAN(initincf,(ITG *nface,double *dmin,double *vfa,
                             double *umfa,double *cvfa,double *hcfa,
                             ITG *ithermal,double *tincfguess,
                             ITG *compressible));

void insert(ITG *ipointer,ITG **mast1p,ITG **mast2p,ITG *i1,
            ITG *i2,ITG *ifree,ITG *nzs_);

void insertfem(int *ipointer,int **mast1p,int **mast2p,int *i1,
            int *i2,int *ifree,int *nzs_);

void insertfreq(ITG *ipointer,ITG **mast1p,ITG **nextp,ITG *i1,
                ITG *i2,ITG *ifree,ITG *nzs_);

void insertrad(ITG *ipointer,ITG **mast1p,ITG **mast2p,ITG *i1,
            ITG *i2,ITG *ifree,ITG *nzs_);

void FORTRAN(integral_boundary,(double *sumfix,double *sumfree,ITG *ifaext,
                                ITG *nfaext,ITG *ielfa,ITG *ifabou,double *vfa,ITG *ipnei,
                                double *xxn));

void FORTRAN(interpolatestate,(ITG *ne,ITG *ipkon,ITG *kon,char *lakon,
             ITG *ne0,ITG *mi,double *xstate,
             double *pslavsurf,ITG *nstate_,double *xstateini,
             ITG *islavsurf,ITG *islavsurfold,
             double *pslavsurfold,char *tieset,ITG *ntie,ITG *itiefac));

void interpolcycsymcfd(int *nkold,double *cotet,int *neold,int *ipkon,
     int *kon,int **nodempcp,int *ipompc,int *nmpc,
     int *ikmpc,int *ilmpc,double **coefmpcp,char *labmpc,
     int *mpcfree,int *memmpc_,char *lakon,int *ncs,int *nslav,
     int *ithermal,double *cs,int *inoslav,int *inomast,int *ics,int *islav);

void FORTRAN(islavactive,(char *tieset,ITG *ntie,ITG *itietri,double *cg,
              double *straight,double *co,double *vold,double *xo,
              double *yo,double *zo,double *x,double *y,double *z,
              ITG *nx,ITG *ny,ITG *nz,ITG *mi,ITG *imastop,ITG *nslavnode,
              ITG *islavnode,ITG *islavact));

void FORTRAN(isortid,(ITG *ix,double *dy,ITG *n,ITG *kflag));

void FORTRAN(isortii,(ITG *ix,ITG *iy,ITG *n,ITG *kflag));

void FORTRAN(isortiid,(ITG *ix,ITG *iy,double *dy,ITG *n,ITG *kflag));

void FORTRAN(isortiddc,(ITG *ix,double *dy1,double *dy2,char *cy,ITG *n,
                         ITG *kflag));

void FORTRAN(isortiiddc,(ITG *ix1,ITG *ix2,double *dy1,double *dy2,
                         char *cy,ITG *n,ITG *kflag));

void FORTRAN(jouleheating,(ITG *ipkon,char *lakon,ITG *kon,double *co,
             double *elcon,ITG *nelcon,ITG *mi,ITG *ne,double *sti,
             ITG *ielmat,ITG *nelemload,char *sideload,double *xload,
             ITG *nload,ITG *nload_,ITG *iamload,ITG *nam,ITG *idefload,
             ITG *ncmat_,ITG *ntmat_,double *alcon,ITG *nalcon,ITG *ithermal,
             double *vold,double *t1));

void FORTRAN(keystart,(ITG *ifreeinp,ITG *ipoinp,ITG *inp,char *name,
           ITG *iline,ITG *ikey));
  
void linstatic(double *co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
             ITG *ne,
             ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
             ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
             ITG *nmpc,
             ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
             ITG *nelemload,char *sideload,double *xload,
             ITG *nload,ITG *nactdof,
             ITG **icolp,ITG *jq,ITG **irowp,ITG *neq,ITG *nzl,
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
             ITG *inotr,ITG *ntrans,double *fmpc,char *cbody,ITG *ibody,
             double *xbody,ITG *nbody,double *xbodyold,double *timepar,
             double *thicke,char *jobnamec,char *tieset,ITG *ntie,
             ITG *istep,ITG *nmat,ITG *ielprop,double *prop,char *typeboun,
             ITG *mortar,ITG *mpcinfo,double *tietol,ITG *ics,ITG *icontact,
             char *orname,ITG *itempuser);

void FORTRAN(localaxes,(ITG *ibody,ITG *nbody,double *xbody,double *e1,
                        double *e2,double *xn));

void FORTRAN(localaxescs,(double *cs,ITG *mcs,double *e1,
                        double *e2,double *xn));

void FORTRAN(lump,(double *adb,double *aub,double *adl,int *irow,int *jq,
                   int *neq));

void FORTRAN(mafillcorio,(double *co,ITG *nk,ITG *kon,ITG *ipkon,
               char *lakon,
               ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
               ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
               double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,
               double *ad,double *au,ITG *nactdof,
               ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
               ITG *ilboun,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,
               double *prestr,ITG *iprestr,double *vold,
               ITG *iperturb,double *sti,ITG *nzs,double *stx,
               double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
               ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *ibody,ITG *ielprop,double *prop));

void FORTRAN(mafilldm,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
               ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
               double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,
               double *ad,double *au,ITG *nactdof,
               ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
               ITG *ilboun,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,
               double *prestr,ITG *iprestr,double *vold,
               ITG *iperturb,double *sti,ITG *nzs,double *stx,
               double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
               ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *ibody,double *clearini,
               ITG *mortar,double *springarea,double *pslavsurf,
               double *pmastsurf,double *reltime,ITG *nasym));

void FORTRAN(mafilldmss,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,
               ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,
               double *ad,double *au,ITG *nactdof,
               ITG *jq,ITG *irow,ITG *neq,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,double *vold,
               ITG *iperturb,double *sti,double *stx,
               ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
               ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,double *physcon,
               double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *ibody,
               double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
               double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               ITG *nea,ITG *neb,
               double *freq,ITG *ndamp,double *dacon));

void mafilldmssmain(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,
               double *ad,double *au,ITG *nactdof,
               ITG *jq,ITG *irow,ITG *neq,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,double *vold,
               ITG *iperturb,double *sti,ITG *nzs,double *stx,
               ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
               ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,double *physcon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *ibody,
               double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
               double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               double *freq,ITG *ndamp,double *dacon);

void *mafilldmssmt(ITG *i);

void FORTRAN(mafillem,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
               ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
               double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,
               double *ad,double *au,double *bb,ITG *nactdof,
               ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
               ITG *ilboun,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,
               double *prestr,ITG *iprestr,double *vold,
               ITG *iperturb,double *sti,ITG *nzs,double *stx,
               double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
               ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
               double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
               double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,ITG *iactive,double *h0,double *pslavsurf,
               double *pmastsurf,ITG *mortar,double *clearini,
               ITG *ielprop,double *prop,ITG *iponoel,ITG *inoel,
               ITG *network));

void FORTRAN(mafillfreq_em,(double *ad,double *au,double *adb,double *aub,
             ITG *irow,ITG *jq,ITG *neq,double *adfreq,double *aubfreq,
             ITG *irowfreq,ITG *iaux,ITG *jqfreq,ITG *icolfreq,ITG *neqfreq,
             ITG *nzsfreq,double *om,ITG *symmetryflag,ITG *inputformat,
             double *b,double *bfreq));

void FORTRAN(mafillk,(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
               double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umfa,double *alet,
               double *ale,double *gradkfa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cpfa,double *hcfa,double *cvel,
               double *gradvel,double *xload,double *gammat,double *xrlfa,
               double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
               double *flux,ITG *nefa,ITG *nefb,ITG *iau6,double *xxni,
               double *xxnj,ITG *iturbulent,double *f1,double *of2,
	       double *yy,double *umel,double *gradkel,double *gradoel,
	       double *sc));

void FORTRAN(mafillkcomp,(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
               double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umfa,double *alet,
               double *ale,double *gradkfa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cvfa,double *hcfa,double *cvel,
               double *gradvel,double *xload,double *xrlfa,
               double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
               double *flux,ITG *nefa,ITG *nefb,ITG *iau6,double *xxni,
               double *xxnj,ITG *iturbulent,double *f1,double *of2,
               double *yy,double *umel,double *gradkel,double *gradoel,
	       ITG *inlet,double *sc));

void mafillkcompmain(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
               double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umfa,double *alet,
               double *ale,double *gradkfa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cpfa,double *hcfa,double *cvel,
               double *gradvel,double *xload,double *xrlfa,
               double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
               double *flux,ITG *iau6,double *xxni,double *xxnj,
               ITG *iturbulent,double *f1,double *of2,double *yy,
               double *umel,double *gradkel,double *gradoel,ITG *inlet,
	       double *sc);

void *mafillkcompmt(ITG *i);

void mafillkmain(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
               double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umfa,double *alet,
               double *ale,double *gradkfa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cpfa,double *hcfa,double *cvel,
               double *gradvel,double *xload,double *gammat,double *xrlfa,
               double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
               double *flux,ITG *iau6,double *xxni,double *xxnj,
               ITG *iturbulent,double *f1,double *of2,double *yy,
	       double *umel,double *gradkel,double *gradoel,double *sc);

void *mafillkmt(ITG *i);

void FORTRAN(mafillnet,(ITG *itg,ITG *ieg,ITG *ntg,
                        double *ac,ITG *nload,char *sideload,
                        ITG *nelemload,double *xloadact,char *lakon,
                        ITG *ntmat_,double *v,double *shcon,ITG *nshcon,
                        ITG *ipkon,ITG *kon,double *co,ITG *nflow,
                        ITG *iinc,ITG *istep,
                        double *dtime,double *ttime,double *time,
                        ITG *ielmat,ITG *nteq,double *prop,
                        ITG *ielprop,ITG *nactdog,ITG *nacteq,
                        double *physcon,double *rhcon,ITG *nrhcon,
                        ITG *ipobody,ITG *ibody,double *xbody,ITG *nbody,
                        double *vold,double *xloadold,double *reltime,
                        ITG *nmethod,char *set,ITG *mi,ITG *nmpc,
                        ITG *nodempc,ITG *ipompc,double *coefmpc,
                        char *labmpc,ITG *iaxial,double *cocon,ITG *ncocon,
                        ITG *iponoel,ITG *inoel));

void FORTRAN(mafillo,(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
               double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umfa,double *alet,
               double *ale,double *gradofa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cpfa,double *hcfa,double *cvel,
               double *gradvel,double *xload,double *gammat,double *xrlfa,
               double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
               double *flux,ITG *nefa,ITG *nefb,ITG *iau6,double *xxni,
               double *xxnj,ITG *iturbulent,double *f1,double *of2,
	       double *gradkel,double *gradoel,double *sc));

void FORTRAN(mafillocomp,(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
               double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umfa,double *alet,
               double *ale,double *gradofa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cpfa,double *hcfa,double *cvel,
               double *gradvel,double *xload,double *xrlfa,
               double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
               double *flux,ITG *nefa,ITG *nefb,ITG *iau6,double *xxni,
               double *xxnj,ITG *iturbulent,double *f1,double *of2,
               double *gradkel,double *gradoel,ITG *inlet,double *sc));

void mafillocompmain(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
               double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umel,double *alet,
               double *ale,double *gradofa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cpfa,double *hcfa,double *cvel,
               double *gradvel,double *xload,double *xrlfa,
               double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
               double *flux,ITG *iau6,double *xxni,double *xxnj,
               ITG *iturbulent,double *f1,double *of2,double *gradkel,
               double *gradoel,ITG *inlet,double *sc);

void *mafillocompmt(ITG *i);

void mafillomain(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
               double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umel,double *alet,
               double *ale,double *gradofa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cpfa,double *hcfa,double *cvel,
               double *gradvel,double *xload,double *gammat,double *xrlfa,
               double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
               double *flux,ITG *iau6,double *xxni,double *xxnj,
               ITG *iturbulent,double *f1,double *of2,double *gradkel,
	       double *gradoel,double *sc);

void *mafillomt(ITG *i);

void FORTRAN(mafillp,(ITG *ne,char *lakonf,ITG *ipnei,
             ITG *neifa,ITG *neiel,double *vfa,double *area,double *adfa,
             double *xlet,double *cosa,double *volume,double *au,double *ad,
             ITG *jq,ITG *irow,double *ap,ITG *ielfa,ITG *ifabou,
             double *xle,double *b,double *xxn,
             ITG *neq,ITG *nzs,double *hfa,double *gradpel,
             double *bp,double *xxi,ITG *neij,double *xlen,double *cosb,
             ITG *nefa,ITG *nefb,ITG *iau6,double *xxicn,double *flux));

void FORTRAN(mafillpbc,(ITG *nef,double *au,double *ad,ITG *jq,ITG *irow,
                        double *b,ITG *iatleastonepressurebc,ITG *nzs));

void FORTRAN(mafillpcomp,(ITG *ne,char *lakonf,ITG *ipnei,
             ITG *neifa,ITG *neiel,double *vfa,double *area,double *adfa,
             double *xlet,double *cosa,double *volume,double *au,double *ad,
             ITG *jq,ITG *irow,double *ap,ITG *ielfa,ITG *ifabou,
             double *xle,double *b,double *xxn,ITG *neq,
             ITG *nzs,double *hfa,double *gradpel,
             double *bp,double *xxi,ITG *neij,double *xlen,double *cosb,
             ITG *ielmatf,ITG *mi,double *a1,double *a2,double *a3,double *velo,
             double *veloo,double *dtimef,double *shcon,ITG *ntmat_,double *vel,
             ITG *nactdohinv,double *xrlfa,double *flux,ITG *nefa,ITG *nefb,
             ITG *iau6,double *xxicn,double *gamma,ITG *inlet));

void mafillpcompmain(ITG *ne,char *lakonf,ITG *ipnei,
             ITG *neifa,ITG *neiel,double *vfa,double *area,double *adfa,
             double *xlet,double *cosa,double *volume,double *au,double *ad,
             ITG *jq,ITG *irow,double *ap,ITG *ielfa,ITG *ifabou,
             double *xle,double *b,double *xxn,ITG *neq,
             ITG *nzs,double *hfa,double *gradpel,
             double *bp,double *xxi,ITG *neij,double *xlen,double *cosb,
             ITG *ielmatf,ITG *mi,double *a1,double *a2,double *a3,double *velo,
             double *veloo,double *dtimef,double *shcon,ITG *ntmat_,double *vel,
             ITG *nactdohinv,double *xrlfa,double *flux,ITG *iau6,
             double *xxicn,double *gamma,ITG *inlet);

void *mafillpcompmt(ITG *i);

void mafillpmain(ITG *ne,char *lakonf,ITG *ipnei,
             ITG *neifa,ITG *neiel,double *vfa,double *area,double *adfa,
             double *xlet,double *cosa,double *volume,double *au,double *ad,
             ITG *jq,ITG *irow,double *ap,ITG *ielfa,ITG *ifabou,
             double *xle,double *b,double *xxn,
             ITG *neq,ITG *nzs,double *hfa,double *gradpel,
             double *bp,double *xxi,ITG *neij,double *xlen,double *cosb,
             ITG *iatleastonepressurebc,ITG *iau6,double *xxicn,
             double *flux);

void *mafillpmt(ITG *i);

void FORTRAN(mafillsm,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
               ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
               double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,
               double *ad,double *au,double *bb,ITG *nactdof,
               ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
               ITG *ilboun,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,
               double *prestr,ITG *iprestr,double *vold,
               ITG *iperturb,double *sti,ITG *nzs,double *stx,
               double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
               ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
               double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
               double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               double *fnext,ITG *nea,ITG *neb,ITG *kscale,ITG *iponoel,
               ITG *inoel,ITG *network,double *smscale,ITG *mscalmethod));

void FORTRAN(mafillsmcsse,(double *co,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nelemload,char *sideload,double *xload,
               ITG *nload,double *xbody,ITG *ipobody,ITG *nbody,
               double *cgr,ITG *nactdof,ITG *neq,ITG *nmethod,ITG *ikmpc,
               ITG *ilmpc,double *elcon,ITG *nelcon,double *rhcon,
               ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,
               ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,
               ITG *ntmat_,double *t0,double *t1,ITG *ithermal,
               ITG *iprestr,double *vold,ITG *iperturb,double *sti,
               double *stx,ITG *iexpl,double *plicon,ITG *nplicon,
               double *plkcon,ITG *nplkcon,double *xstiff,ITG *npmat_,
               double *dtime,char *matname,ITG *mi,ITG *ncmat_,ITG *mass,
               ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *ttime,double *time,
               ITG *istep,ITG *iinc,ITG *coriolis,ITG *ibody,
               double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
               double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               ITG *nea,ITG *neb,double *distmin,ITG *ndesi,
               ITG *nodedesi,double *df,ITG *jqs,
               ITG *irows,double *dfminds,ITG *icoordinate,
               double *dxstiff,double *xdesi,ITG *istartelem,ITG *ialelem,
               double *v,double *sigma,char *labmpc,ITG *ics,double *cs,
               ITG *mcs,ITG *nk,ITG *nzss));

void FORTRAN(mafillsmse,(double *co,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nelemload,char *sideload,double *xload,
               ITG *nload,double *xbody,ITG *ipobody,ITG *nbody,
               double *cgr,ITG *nactdof,ITG *neq,ITG *nmethod,ITG *ikmpc,
               ITG *ilmpc,double *elcon,ITG *nelcon,double *rhcon,
               ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,
               ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,
               ITG *ntmat_,double *t0,double *t1,ITG *ithermal,
               ITG *iprestr,double *vold,ITG *iperturb,double *sti,
               double *stx,ITG *iexpl,double *plicon,ITG *nplicon,
               double *plkcon,ITG *nplkcon,double *xstiff,ITG *npmat_,
               double *dtime,char *matname,ITG *mi,ITG *ncmat_,ITG *mass,
               ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *ttime,double *time,
               ITG *istep,ITG *iinc,ITG *coriolis,ITG *ibody,
               double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
               double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               ITG *nea,ITG *neb,double *distmin,ITG *ndesi,
               ITG *nodedesi,double *df,ITG *jqs,
               ITG *irows,double *dfl,ITG *icoordinate,
               double *dxstiff,double *xdesi,ITG *istartelem,ITG *ialelem,
               double *v,double *sigma,ITG *ieigenfrequency));

void FORTRAN(mafillsmse2,(double *co,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nelemload,char *sideload,double *xload,
               ITG *nload,double *xbody,ITG *ipobody,ITG *nbody,
               double *cgr,ITG *nactdof,ITG *neq,ITG *nmethod,ITG *ikmpc,
               ITG *ilmpc,double *elcon,ITG *nelcon,double *rhcon,
               ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,
               ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,
               ITG *ntmat_,double *t0,double *t1,ITG *ithermal,
               ITG *iprestr,double *vold,ITG *iperturb,double *sti,
               double *stx,ITG *iexpl,double *plicon,ITG *nplicon,
               double *plkcon,ITG *nplkcon,double *xstiff,ITG *npmat_,
               double *dtime,char *matname,ITG *mi,ITG *ncmat_,ITG *mass,
               ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *ttime,double *time,
               ITG *istep,ITG *iinc,ITG *coriolis,ITG *ibody,
               double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
               double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               ITG *nea,ITG *neb,double *distmin,ITG *ndesi,
               ITG *nodedesi,double *df2,ITG *jqs2,
               ITG *irows2,double *dfl2,ITG *icoordinate,
               double *dxstiff,double *xdesi,ITG *istartelem,ITG *ialelem,
               double *v,double *sigma,ITG *ieigenfrequency));

void *mafillsmmt(ITG *i);

void *mafillsmse2mt(ITG *i);

void *mafillsmsemt(ITG *i);

void mafillsmmain(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
               ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
               double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,
               double *ad,double *au,double *bb,ITG *nactdof,
               ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
               ITG *ilboun,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,
               double *prestr,ITG *iprestr,double *vold,
               ITG *iperturb,double *sti,ITG *nzs,double *stx,
               double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
               ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
               double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
               double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               double *fnext,ITG *kscale,ITG *iponoel,ITG *inoel,
               ITG *network,ITG *ntrans,ITG *inotr,double *trab,
               double *smscale,ITG *mscalmethod);

void mafillsmmain_se(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
               ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
               double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,ITG *nactdof,ITG *neq,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
               ITG *ilboun,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,
               double *prestr,ITG *iprestr,double *vold,
               ITG *iperturb,double *sti,double *stx,
               ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
               ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
               double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
               double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               double *fnext,double *distmin,ITG *ndesi,ITG *nodedesi,
               double *df,ITG *nzss,ITG *jqs,ITG *irows,
               ITG *icoordinate,double *dxstiff,double *xdesi,
               ITG *istartelem,ITG *ialelem,double *v,double *sigma,
               ITG *cyclicsymmetry,char *labmpc,ITG *ics,double *cs,
               ITG *mcs,ITG *ieigenfrequency);

void mafillsmmain_se2(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
               ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
               double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,ITG *nactdof,ITG *neq,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
               ITG *ilboun,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,
               double *prestr,ITG *iprestr,double *vold,
               ITG *iperturb,double *sti,double *stx,
               ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
               ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
               double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
               double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               double *fnext,double *distmin,ITG *ndesi,ITG *nodedesi,
               double *df2,ITG *nzss2,ITG *jqs2,ITG *irows2,
               ITG *icoordinate,double *dxstiff,double *xdesi,
               ITG *istartelem,ITG *ialelem,double *v,double *sigma,
               ITG *cyclicsymmetry,char *labmpc,ITG *ics,double *cs,
               ITG *mcs,ITG *ieigenfrequency);

void FORTRAN(mafillsmas,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
               ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
               double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,
               double *ad,double *au,double *bb,ITG *nactdof,
               ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
               ITG *ilboun,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,
               double *prestr,ITG *iprestr,double *vold,
               ITG *iperturb,double *sti,ITG *nzs,double *stx,
               double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
               ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
               double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
               double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               ITG *kscale,ITG *iponoel,ITG *inoel,ITG *network));

void FORTRAN(mafillsmas1,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
               ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
               double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,
               double *ad,double *au,double *bb,ITG *nactdof,
               ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
               ITG *ilboun,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,
               double *prestr,ITG *iprestr,double *vold,
               ITG *iperturb,double *sti,ITG *nzs,double *stx,
               double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
               ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
               double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
               double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               ITG *nea,ITG *neb,ITG *kscale));

void mafillsmasmain(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
               ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
               double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,
               double *ad,double *au,double *bb,ITG *nactdof,
               ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
               ITG *ilboun,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,
               double *prestr,ITG *iprestr,double *vold,
               ITG *iperturb,double *sti,ITG *nzs,double *stx,
               double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
               ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
               double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
               double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               ITG *kscale);

void *mafillsmasmt(ITG *i);

void FORTRAN(mafillsmcs,(double *co,ITG *nk,ITG *kon,ITG *ipkon,
               char *lakon,
               ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
               ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
               double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,
               double *ad,double *au,double *bb,ITG *nactdof,
               ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
               ITG *ilboun,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,
               double *prestr,ITG *iprestr,double *vold,
               ITG *iperturb,double *sti,ITG *nzs,double *stx,
               double *adb,double *aub,ITG *iexpl,double *plicon,
               ITG *nplicon,double *plkcon,ITG *nplkcon,double *xstiff,
               ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ics,double *cs,ITG *nm,ITG *ncmat_,char *labmpc,
               ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,ITG *mcs,ITG *coriolis,ITG *ibody,
               double *xloadold,double *reltime,ITG *ielcs,double *veold,
               double *springarea,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               ITG *kscale,double *xstateini,double *xstate,ITG *nstate_));

void FORTRAN(mafillsmcsas,(double *co,ITG *nk,ITG *kon,ITG *ipkon,
               char *lakon,
               ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
               ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
               double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,
               double *ad,double *au,double *bb,ITG *nactdof,
               ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
               ITG *ilboun,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,
               double *prestr,ITG *iprestr,double *vold,
               ITG *iperturb,double *sti,ITG *nzs,double *stx,
               double *adb,double *aub,ITG *iexpl,double *plicon,
               ITG *nplicon,double *plkcon,ITG *nplkcon,double *xstiff,
               ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ics,double *cs,ITG *nm,ITG *ncmat_,char *labmpc,
               ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,ITG *mcs,ITG *coriolis,ITG *ibody,
               double *xloadold,double *reltime,ITG *ielcs,double *veold,
               double *springarea,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,ITG *nstate_,double *xstateini,double *xstate,
               double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               ITG *kscale));

void FORTRAN(mafillsmforc,(ITG *nforc,ITG *ndirforc,ITG *nodeforc,
             double *xforc,ITG *nactdof,double *fext,ITG *ipompc,
             ITG *nodempc,double *coefmpc,ITG *mi,ITG *rhsi,double *fnext,
             ITG *nmethod,ITG *ntrans,ITG *inotr,double *trab,double *co));

void FORTRAN(mafillsm_company,(double *co,ITG *nk,ITG *kon,ITG *ipkon,
               char *lakon,
               ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun,
               ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
               double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,
               double *ad,double *au,double *bb,ITG *nactdof,
               ITG *icol,ITG *jq,ITG *irow,ITG *neq,ITG *nzl,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
               ITG *ilboun,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,
               double *prestr,ITG *iprestr,double *vold,
               ITG *iperturb,double *sti,ITG *nzs,double *stx,
               double *adb,double *aub,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff,
               ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhs,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *kinc,ITG *coriolis,ITG *ibody,
               double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
               double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
               ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
               double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               double *fnext,ITG *kscale,ITG *iponoel,
               ITG *inoel,ITG *network,ITG *ntrans,ITG *inotr,double *trab));

void FORTRAN(mafillt,(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
               double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umel,double *alet,
               double *ale,double *gradtfa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cpfa,double *hcfa,double *cvel,
               double *gradvel,double *xload,double *gammat,double *xrlfa,
               double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
               double *flux,ITG *nefa,ITG *nefb,ITG *iau6,double *xxni,
               double *xxnj,ITG *iturbulent,double *of2,double *sc));

void FORTRAN(mafilltcomp,(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
               double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umel,double *alet,
               double *ale,double *gradtfa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cpfa,double *hcfa,double *cvel,
               double *gradvel,double *xload,double *gammat,double *xrlfa,
               double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
               double *flux,ITG *nefa,ITG *nefb,ITG *iau6,double *xxni,
	       double *xxnj,ITG *iturbulent,double *of2,double *sc));

void mafilltcompmain(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
               double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umel,double *alet,
               double *ale,double *gradtfa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cpfa,double *hcfa,double *cvel,
               double *gradvel,double *xload,double *gammat,double *xrlfa,
               double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
               double *flux,ITG *iau6,double *xxni,double *xxnj,
	       ITG *iturbulent,double *of2,double *sc);

void *mafilltcompmt(ITG *i);

void mafilltmain(ITG *nef,ITG *ipnei,ITG *neifa,
               ITG *neiel,double *vfa,double *xxn,double *area,
               double *au,double *ad,ITG *jq,ITG *irow,ITG *nzs,
               double *b,double *vel,double *umel,double *alet,
               double *ale,double *gradtfa,double *xxi,double *body,
               double *volume,ITG *ielfa,char *lakonf,
               ITG *ifabou,ITG *nbody,ITG *neq,double *dtimef,double *velo,
               double *veloo,double *cpfa,double *hcfa,double *cvel,
               double *gradvel,double *xload,double *gammat,double *xrlfa,
               double *xxj,ITG *nactdohinv,double *a1,double *a2,double *a3,
               double *flux,ITG *iau6,double *xxni,double *xxnj,
               ITG *iturbulent,double *of2,double *sc);

void *mafilltmt(ITG *i);

void FORTRAN(mafillv,(ITG *nef,ITG *ipnei,ITG *neifa,ITG *neiel,
             double *vfa,double *xxn,double *area,double *au,double *ad,
             ITG *jq,ITG *irow,ITG *nzs,double *b,double *vel,double *cosa,
             double *umfa,double *alet,double *ale,double *gradvfa,
             double *xxi,double *body,double *volume,
             ITG *ielfa,char *lakonf,ITG *ifabou,ITG *nbody,
             double *dtimef,double *velo,double *veloo,
             double *sel,double *xrlfa,double *gamma,double *xxj,
             ITG *nactdohinv,double *a1,double *a2,double *a3,double *flux,
             ITG *nefa,ITG *nefb,ITG *icyclic,double *c,ITG *ifatie,
             ITG *iau6,double *xxna,double *xxnj,ITG *iturbulent,
             double *gradvel,double *of2,double *yy,double *umel,
	     ITG *ncfd,ITG *inlet,double *sc));

void FORTRAN(mafillvcomp,(ITG *nef,ITG *ipnei,ITG *neifa,ITG *neiel,
             double *vfa,double *xxn,double *area,double *au,double *ad,
             ITG *jq,ITG *irow,ITG *nzs,double *b,double *vel,double *cosa,
             double *umfa,double *alet,double *ale,double *gradvfa,
             double *xxi,double *body,double *volume,
             ITG *ielfa,char *lakonf,ITG *ifabou,ITG *nbody,
             double *dtimef,double *velo,double *veloo,
             double *sel,double *xrlfa,double *gamma,double *xxj,
             ITG *nactdohinv,double *a1,double *a2,double *a3,double *flux,
             ITG *nefa,ITG *nefb,ITG *icyclic,double *c,ITG *ifatie,
             ITG *iau6,double *xxni,double *xxnj,ITG *iturbulent,
             double *gradvel,double *of2,double *yy,double *umel,
	     ITG *ncfd,ITG *inlet,double *sc));

void mafillvcompmain(ITG *nef,ITG *ipnei,ITG *neifa,ITG *neiel,
             double *vfa,double *xxn,double *area,double *au,double *ad,
             ITG *jq,ITG *irow,ITG *nzs,double *b,double *vel,double *cosa,
             double *umfa,double *alet,double *ale,double *gradvfa,
             double *xxi,double *body,double *volume,
             ITG *ielfa,char *lakonf,ITG *ifabou,ITG *nbody,
             double *dtimef,double *velo,double *veloo,
             double *sel,double *xrlfa,double *gamma,double *xxj,
             ITG *nactdohinv,double *a1,double *a2,double *a3,
             double *flux,ITG *icyclic,double *c,ITG *ifatie,ITG *iau6,
             double *xxna,double *xxnj,ITG *iturbulent,double *gradvel,
	     double *of2,double *yy,double *umel,ITG *ncfd,ITG *inlet,
	     double *sc);

void *mafillvcompmt(ITG *i);

void mafillvmain(ITG *nef,ITG *ipnei,ITG *neifa,ITG *neiel,
             double *vfa,double *xxn,double *area,double *au,double *ad,
             ITG *jq,ITG *irow,ITG *nzs,double *b,double *vel,double *cosa,
             double *umfa,double *alet,double *ale,double *gradvfa,
             double *xxi,double *body,double *volume,
             ITG *ielfa,char *lakonf,ITG *ifabou,ITG *nbody,
             double *dtimef,double *velo,double *veloo,
             double *sel,double *xrlfa,double *gamma,double *xxj,
             ITG *nactdohinv,double *a1,double *a2,double *a3,
             double *flux,ITG *icyclic,double *c,ITG *ifatie,
             ITG *iau6,double *xxni,double *xxnj,ITG *iturbulent,
             double *gradvel,double *of2,double *yy,double *umel,
             ITG *ncfd,ITG *inlet,double *sc);

void *mafillvmt(ITG *i);

void FORTRAN(map3dto1d2d,(double *extnor,ITG *ipkon,ITG *inum,ITG *kon,
                          char *lakon,ITG *nfield,ITG *nk,ITG *ne,
                          char *cflag,double *co,double *vold,ITG *iforce,
                          ITG *mi,ITG *ielprop,double *prop));

void mastruct(ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
              ITG *nodeboun,ITG *ndirboun,ITG *nboun,ITG *ipompc,
              ITG *nodempc,ITG *nmpc,ITG *nactdof,ITG *icol,
              ITG *jq,ITG **mast1p,ITG **irowp,ITG *isolver,ITG *neq,
              ITG *ikmpc,ITG *ilmpc,ITG *ipointer,ITG *nzs,ITG *nmethod,
              ITG *ithermal,ITG *ikboun,ITG *ilboun,ITG *iperturb,
              ITG *mi,ITG *mortar,char *typeboun,char *labmpc,
              ITG *iit,ITG *icascade,ITG *network,ITG *iexpl);

void mastructcs(ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,ITG *nodeboun,
               ITG *ndirboun,ITG *nboun,ITG *ipompc,ITG *nodempc,
               ITG *nmpc,ITG *nactdof,ITG *icol,ITG *jq,ITG **mast1p,
               ITG **irowp,ITG *isolver,ITG *neq,
               ITG *ikmpc,ITG *ilmpc,ITG *ipointer,
               ITG *nzs,ITG *nmethod,ITG *ics,double *cs,
               char *labmpc,ITG *mcs,ITG *mi,ITG *mortar);

void mastructem(ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
              ITG *nodeboun,ITG *ndirboun,ITG *nboun,ITG *ipompc,
              ITG *nodempc,ITG *nmpc,ITG *nactdof,ITG *icol,
              ITG *jq,ITG **mast1p,ITG **irowp,ITG *isolver,ITG *neq,
              ITG *ikmpc,ITG *ilmpc,ITG *ipointer,ITG *nzs,
              ITG *ithermal,ITG *mi,ITG *ielmat,double *elcon,ITG *ncmat_,
              ITG *ntmat_,ITG *inomat,ITG *network);

void mastructf(ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
               ITG *icol,ITG *jq,ITG **mast1p,ITG **irowp,
               ITG *isolver,ITG *ipointer,ITG *nzs,
               ITG *ipnei,ITG *ineiel,ITG *mi);

void mastructffem(int *nk,int *kon,int *ipkon,char *lakon,int *ne,
              int *nodeboun,int *ndirboun,int *nboun,int *ipompc,
              int *nodempc,int *nmpc,int *nactdoh,int *icolt,
              int *icolv,int *icolp,int *icolk,int *jqt,int *jqv,int *jqp,
              int *jqk,int **mast1p,int **irowtp,int **irowvp,int **irowpp,
              int **irowkp,int *isolver,int *neqt,int *neqv,int *neqp,
              int *neqk,int *ikmpc,int *ilmpc,int *ipointer,
              int *nzst,int *nzsv,int *nzsp,int *nzsk,
              int *ithermal,int *ikboun,int *ilboun,int *turbulent,
              int *nactdok,int *ifreestream,int *nfreestream,
              int *isolidface,int *nsolidface,int *nzs,int *iexplicit,
             int *ielmat,int *inomat,char *labmpc);

void mastructrad(ITG *ntr,ITG *nloadtr,char *sideload,ITG *ipointerrad,
              ITG **mast1radp,ITG **irowradp,ITG *nzsrad,
              ITG *jqrad,ITG *icolrad);

void mastructrand(ITG *icols,ITG *jqs,ITG **mast1p,ITG **irowsp,
                  ITG *ipointer,ITG *nzss,
                  ITG *ndesi,double *physcon,double *xo,double *yo,
                  double *zo,double *x,double *y,double *z,ITG *nx,
                  ITG *ny,ITG *nz);

void mastructse(ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
              ITG *ipompc,ITG *nodempc,ITG *nmpc,
              ITG *nactdof,ITG *jqs,ITG **mast1p,ITG **irowsp,
              ITG *ipointer,ITG *nzss,ITG *mi,ITG *mortar,
              ITG *nodedesi,ITG *ndesi,ITG *icoordinate,ITG *ielorien,
              ITG *istartdesi,ITG *ialdesi);

void mastructse2(ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
              ITG *ipompc,ITG *nodempc,ITG *nmpc,
              ITG *nactdof,ITG *jqs2,ITG **mast1p,ITG **irows2p,
              ITG *ipointer,ITG *nzss2,ITG *mi,ITG *mortar,
              ITG *nodedesi,ITG *ndesi,ITG *icoordinate,ITG *ielorien,
              ITG *istartdesi,ITG *ialdesi);

void FORTRAN(materialdata_cfd,(ITG *nef,double *vel,double *shcon,
             ITG *nshcon,ITG *ielmat,ITG *ntmat_,ITG *mi,double *cvel,
             double *vfa,double *cocon,ITG *ncocon,double *physcon,
             double *cvfa,ITG *ithermal,ITG *nface,double *umel,
             double *umfa,ITG *ielfa,double *hcfa,double *rhcon,
             ITG *nrhcon));

void FORTRAN(materialdata_cfd1,(ITG *nef,double *vel,double *shcon,
             ITG *nshcon,ITG *ielmat,ITG *ntmat_,ITG *mi,double *cvel,
             double *physcon,ITG *ithermal,double *umel,double *rhcon,
             ITG *nrhcon,ITG *nefa,ITG *nefb));

void *materialdata_cfd1mt(ITG *i);

void FORTRAN(materialdata_cfd2,(double *shcon,
             ITG *nshcon,ITG *ielmat,ITG *ntmat_,ITG *mi,
             double *vfa,double *cocon,ITG *ncocon,double *physcon,
             double *cvfa,ITG *ithermal,
             double *umfa,ITG *ielfa,double *hcfa,double *rhcon,
             ITG *nrhcon,ITG *nfacea,ITG *nfaceb));

void *materialdata_cfd2mt(ITG *i);

void FORTRAN(materialdata_cfd_comp,(ITG *nef,double *vel,double *shcon,
             ITG *nshcon,ITG *ielmat,ITG *ntmat_,ITG *mi,double *cvel,
             double *vfa,double *cocon,ITG *ncocon,double *physcon,
             double *cvfa,ITG *ithermal,ITG *nface,double *umel,
             double *umfa,ITG *ielfa,double *hcfa));

void FORTRAN(materialdata_cfd_comp1,(ITG *nef,double *vel,double *shcon,
             ITG *nshcon,ITG *ielmat,ITG *ntmat_,ITG *mi,double *cvel,
             double *physcon,ITG *ithermal,double *umel,
             ITG *nefa,ITG *nefb));

void *materialdata_cfd_comp1mt(ITG *i);

void FORTRAN(materialdata_cfd_comp2,(double *shcon,
             ITG *nshcon,ITG *ielmat,ITG *ntmat_,ITG *mi,
             double *vfa,double *cocon,ITG *ncocon,double *physcon,
             double *cvfa,ITG *ithermal,double *umfa,ITG *ielfa,
             double *hcfa,ITG *nfacea,ITG *nfaceb));

void *materialdata_cfd_comp2mt(ITG *i);

void materialdata_cfd_compmain(ITG *nef,double *vel,double *shcon,
             ITG *nshcon,ITG *ielmatf,ITG *ntmat_,ITG *mi,double *cvel,
             double *vfa,double *cocon,ITG *ncocon,double *physcon,
             double *cvfa,ITG *ithermal,ITG *nface,double *umel,
             double *umfa,ITG *ielfa,double *hcfa,ITG *num_cpus);

void materialdata_cfdmain(ITG *nef,double *vel,double *shcon,
             ITG *nshcon,ITG *ielmatf,ITG *ntmat_,ITG *mi,double *cvel,
             double *vfa,double *cocon,ITG *ncocon,double *physcon,
             double *cvfa,ITG *ithermal,ITG *nface,double *umel,
             double *umfa,ITG *ielfa,double *hcfa,double *rhcon,
             ITG *nrhcon,ITG *num_cpus);

void matrixstorage(double *ad,double **aup,double *adb,double *aub,
                double *sigma,ITG *icol,ITG **irowp,
                ITG *neq,ITG *nzs,ITG *ntrans,ITG *inotr,
                double *trab,double *co,ITG *nk,ITG *nactdof,
                char *jobnamec,ITG *mi,ITG *ipkon,char *lakon,
                ITG *kon,ITG *ne,ITG *mei,ITG *nboun,ITG *nmpc,
                double *cs,ITG *mcs,ITG *ithermal,ITG *nmethod);

void FORTRAN(meannode,(ITG *nk,ITG *inum,double *v));

void FORTRAN(midexternaledges,(ITG *iexternedg,ITG *nexternedg,ITG *iedgext,
                              ITG *ifreeed,ITG *ieled,ITG *ipoeled,
                              ITG *iedg,ITG *iedtet,ITG *kontetor));

void FORTRAN(midexternalfaces,(ITG *iexternfa,ITG *nexternfa,ITG *ifacext,
                               ITG *ifreefa,ITG *itetfa,
                               ITG *ifac,ITG *kontet,ITG *kontetor,
                               ITG *ialsetexternel,ITG *nexternel,
                               ITG *iedgextfa,ITG *ifacexted,
                               ITG *ipoed,ITG *iedg,ITG *iexternedg));

void FORTRAN(mpcrem,(ITG *i,ITG *mpcfree,ITG *nodempc,ITG *nmpc,ITG *ikmpc,
                     ITG *ilmpc,char *labmpc,double *coefmpc,ITG *ipompc));

void mtseed(ITG *iseed);

void FORTRAN(mult,(double *matrix,double *trans,ITG *n));

void FORTRAN(negativepressure,(ITG *ne0,ITG *ne,ITG *mi,double *stx,
                               double *pressureratio));

void FORTRAN(networkelementpernode,(ITG *iponoel,ITG *inoel,char *lakon,
             ITG *ipkon,ITG *kon,ITG *inoelsize,ITG *nflow,ITG *ieg,
             ITG *ne,ITG *network));

void FORTRAN(networkinum,(ITG *ipkon,ITG *inum,ITG *kon,char *lakon,
       ITG *ne,ITG *itg,ITG *ntg));

void FORTRAN(newnodes,(ITG *nktet_,ITG *ipoed,ITG *n,
                       ITG *iedg,double *h,double *d,double *r,
                       double *conewnodes,double *cotet,ITG *ibasenewnodes,
                       ITG *ipoeled,ITG *ieled,double *doubleglob,
                       ITG *integerglob,ITG *nnewnodes,ITG *iedgnewnodes,
                       double *hnewnodes,ITG *iexternedg));

void FORTRAN(newtincf,(ITG *ithermal,double *tincfguess,
                       ITG *compressible,double *vel,double *hcel,
                       double *umel,double *cvel,double *h,
                       double *sc,ITG *iturbulent,ITG *ipkonf,
                       ITG *nmethod,ITG *nef,char *lakonf,
                       double *xxn,ITG *ipnei));

void FORTRAN(nident,(ITG *x,ITG *px,ITG *n,ITG *id));

void FORTRAN(nidentll,(long long *x,long long *px,ITG *n,ITG *id));

void FORTRAN(nodestiedface,(char *tieset,ITG *ntie,ITG *ipkon,ITG *kon,
       char *lakon,char *set,ITG *istartset,ITG *iendset,ITG *ialset,
       ITG *nset,ITG *faceslave,ITG *istartfield,ITG *iendfield,
       ITG *ifield,ITG *nconf,ITG *ncone,char *kind));

void nonlingeo(double **co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
             ITG *ne,
             ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
             ITG **ipompcp,ITG **nodempcp,double **coefmpcp,char **labmpcp,
             ITG *nmpc,
             ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
             ITG **nelemloadp,char **sideloadp,double *xload,
             ITG *nload,ITG *nactdof,
             ITG **icolp,ITG *jq,ITG **irowp,ITG *neq,ITG *nzl,
             ITG *nmethod,ITG **ikmpcp,ITG **ilmpcp,ITG *ikboun,
             ITG *ilboun,
             double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
             double *alcon,ITG *nalcon,double *alzero,ITG **ielmatp,
             ITG **ielorienp,ITG *norien,double *orab,ITG *ntmat_,
             double *t0,double *t1,double *t1old,
             ITG *ithermal,double *prestr,ITG *iprestr,
             double **vold,ITG *iperturb,double *sti,ITG *nzs,
             ITG *kode,char *filab,ITG *idrct,
             ITG *jmax,ITG *jout,double *timepar,
             double *eme,double *xbounold,
             double *xforcold,double *xloadold,
             double *veold,double *accold,
             char *amname,double *amta,ITG *namta,ITG *nam,
             ITG *iamforc,ITG **iamloadp,
             ITG *iamt1,double *alpha,ITG *iexpl,
             ITG *iamboun,double *plicon,ITG *nplicon,double *plkcon,
             ITG *nplkcon,
             double **xstatep,ITG *npmat_,ITG *istep,double *ttime,
             char *matname,double *qaold,ITG *mi,
             ITG *isolver,ITG *ncmat_,ITG *nstate_,ITG *iumat,
             double *cs,ITG *mcs,ITG *nkon,double **ener,ITG *mpcinfo,
             char *output,
             double *shcon,ITG *nshcon,double *cocon,ITG *ncocon,
             double *physcon,ITG *nflow,double *ctrl,
             char *set,ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,ITG *nener,ITG *ikforc,ITG *ilforc,double *trab,
             ITG *inotr,ITG *ntrans,double **fmpcp,char *cbody,
             ITG *ibody,double *xbody,ITG *nbody,double *xbodyold,
             ITG *ielprop,double *prop,ITG *ntie,char *tieset,
             ITG *itpamp,ITG *iviewfile,char *jobnamec,double *tietol,
             ITG *nslavs,double *thicke,ITG *ics,
             ITG *nintpoint,ITG *mortar,ITG *ifacecount,char *typeboun,
             ITG **islavsurfp,double **pslavsurfp,double **clearinip,
             ITG *nmat,double *xmodal,ITG *iaxial,ITG *inext,ITG *nprop,
             ITG *network,char *orname,double *vel,ITG *nef,
             double *velo,double *veloo,double *energy,ITG *itempuser);

void FORTRAN(nonlinmpc,(double *co,double *vold,ITG *ipompc,ITG *nodempc,
                   double *coefmpc,char *labmpc,ITG *nmpc,ITG *ikboun,
                   ITG *ilboun,ITG *nboun,double *xbounact,double *aux,
                   ITG *iaux,ITG *maxlenmpc,ITG *ikmpc,ITG *ilmpc,
                   ITG *icascade,ITG *kon,ITG *ipkon,char *lakon,
                   ITG *ne,double *reltime,ITG *newstep,double *xboun,
                   double *fmpc,ITG *newinc,ITG *idiscon,ITG *ncont,
                   double *trab,ITG *ntrans,ITG *ithermal,ITG *mi));

void FORTRAN(norm,(double *vel,double *velnorm,ITG *nef));

void FORTRAN(normalsforequ_se,(ITG *nk,double *co,ITG *iponoelfa,
                               ITG *inoelfa,ITG *konfa,ITG *ipkonfa,
                               char *lakonfa,ITG *ne,ITG *ipnor,
                               double *xnor,ITG *nodedesiinv,char *jobnamef,
                               ITG *iponexp,ITG *nmpc,char *labmpc,
                               ITG *ipompc,ITG *nodempc,ITG *ipretinfo,
                               ITG *kon,ITG *ipkon,char *lakon,ITG *iponoel,
                               ITG *inoel,ITG *iponor2d,ITG *knor2d,
                               ITG *iponod2dto3d,ITG *ipoface,ITG *nodeface));

void FORTRAN(normalsoninterface,(ITG *istartset,ITG *iendset,
             ITG *ialset,ITG *imast,ITG *ipkon,ITG *kon,char *lakon,
             ITG *imastnode,ITG *nmastnode,double *xmastnor,double *co));

void FORTRAN(normalsonsurface_se,(ITG *ipkon,ITG *kon,char*lakon,
             double *extnor,double *co,ITG *nk,ITG *ipoface,
             ITG *nodface,ITG *nactdof,ITG *mi,ITG *nodedesiinv,
             ITG *iregion,ITG *iponoelfa,ITG *ndesi,ITG *nodedesi,
             ITG *iponod2dto3d,ITG *ikboun,ITG *nboun,ITG *ne2d));

void FORTRAN(objective_disp_tot,(double *dgdx,double *df,ITG *ndesi,
                                 ITG *iobject,ITG *jqs,ITG *irows,
                                 double *dgdu));

void objectivemain_se(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
             ITG *ne,double *v,double *stn,ITG *inum,
             double *stx,
             double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
             double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
             ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
             double *t0,double *t1,ITG *ithermal,double *prestr,
             ITG *iprestr,char *filab,double *eme,double *emn,
             double *een,ITG *iperturb,double *f,double *fn,ITG *nactdof,
             ITG *iout,double *qa,
             double *vold,ITG *nodeboun,ITG *ndirboun,
             double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,
             double *coefmpc,char *labmpc,ITG *nmpc,ITG *nmethod,
             double *cam,ITG *neq,double *veold,double *accold,
             double *bet,double *gam,double *dtime,double *time,
             double *ttime,double *plicon,
             ITG *nplicon,double *plkcon,ITG *nplkcon,
             double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
             double *epn,char *matname,ITG *mi,ITG *ielas,
             ITG *icmd,ITG *ncmat_,ITG *nstate_,double *stiini,
             double *vini,ITG *ikboun,ITG *ilboun,double *ener,
             double *enern,double *emeini,double *xstaten,double *eei,
             double *enerini,double *cocon,ITG *ncocon,char *set,
             ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,double *qfx,double *qfn,double *trab,
             ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,
             ITG *nload,ITG *ikmpc,ITG *ilmpc,ITG *istep,ITG *iinc,
             double *springarea,double *reltime,ITG *ne0,double *xforc,
             ITG *nforc,double *thicke,
             double *shcon,ITG *nshcon,char *sideload,double *xload,
             double *xloadold,ITG *icfd,ITG *inomat,double *pslavsurf,
             double *pmastsurf,ITG *mortar,ITG *islavact,double *cdn,
             ITG *islavnode,ITG *nslavnode,ITG *ntie,double *clearini,
             ITG *islavsurf,ITG *ielprop,double *prop,double *energyini,
             double *energy,double *distmin,
             ITG *ndesi,ITG *nodedesi,ITG *nobject,
             char *objectset,double *g0,double *dgdx,double *sti,
             double *df,ITG *nactdofinv,ITG *jqs,ITG *irows,
             ITG *idisplacement,ITG *nzs,char *jobnamec,ITG *isolver,
             ITG *icol,ITG *irow,ITG *jq,ITG *kode,double *cs,char *output,
             ITG *istartdesi,ITG *ialdesi,double *xdesi,char *orname,
             ITG *icoordinate,ITG *iev,double *d,double *z,double *au,
             double *ad,double *aub,double *adb,ITG *cyclicsymmetry,
             ITG *nzss,ITG *nev,ITG *ishapeenergy,double *fint,
             ITG *nlabel,ITG *igreen,ITG *nasym,ITG *iponoel,ITG *inoel,
             ITG *nodedesiinv,double *dgdxglob,double *df2,double *dgdxdy,
             ITG *nkon,ITG *iponod2dto3d,ITG *iponk2dto3d,ITG *ics,
             ITG *mcs,ITG *mpcend,ITG *noddiam);

void *objectivemt_shapeener_dx(ITG *i);

void *objectivemt_mass_dx(ITG *i);

void FORTRAN(objective_disp,(ITG *nodeset,ITG *istartset,ITG *iendset,
             ITG *ialset,ITG *nk,ITG *idesvar,ITG *iobject,ITG *mi,
             double *g0,ITG *nobject,double *vold,char *objectset));

void FORTRAN(objective_disp_dx,(ITG *nodeset,ITG *istartset,ITG *iendset,
             ITG *ialset,ITG *nk,ITG *idesvar,ITG *iobject,ITG *mi,
             ITG *nactdof,double *dgdx,ITG *ndesi,ITG *nobject,
             double *vold,double *b,char *objectset));

void FORTRAN(objective_freq,(double *dgdx,double *df,double *vold,
             ITG *ndesi,ITG *iobject,ITG *mi,ITG *nactdofinv,ITG *jqs,
             ITG *irows));

void FORTRAN(objective_freq_cs,(double *dgdx,double *df,double *vold,
             ITG *ndesi,ITG *iobject,ITG *mi,ITG *nactdofinv,ITG *jqs,
             ITG *irows,ITG *nk,ITG *nzss));

void FORTRAN(objective_mass_dx,(double *co1,ITG *kon1,ITG *ipkon1,char *lakon1,
             ITG *nelcon1,double *rhcon1,ITG *ielmat1,
             ITG *ielorien1,ITG *norien1,ITG *ntmat1_,
             char *matname1,ITG *mi1,double *thicke1,ITG *mortar1,ITG *nea,
             ITG *neb,ITG *ielprop1,double *prop1,double *distmin1,
             ITG *ndesi1,ITG *nodedesi1,ITG *nobject1,
             double *g01,double *dgdx1,ITG *iobject1,double *xmass1,
             ITG *istartdesi1,ITG *ialdesi1,double *xdesi1,ITG *idesvar));
    
void FORTRAN(objective_shapeener_dx,(double *co1,ITG *kon1,ITG *ipkon1,
             char *lakon1,ITG *ne1,double *stx1,double *elcon1,
             ITG *nelcon1,double *rhcon1,ITG *nrhcon1,double *alcon1,
             ITG *nalcon1,double *alzero1,ITG *ielmat1,ITG *ielorien1,
             ITG *norien1,double *orab1,ITG *ntmat1_,double *t01,
             double *t11,ITG *ithermal1,double *prestr1,ITG *iprestr1,
             ITG *iperturb1,ITG *iout1,double *vold1,ITG *nmethod1,
             double *veold1,double *dtime1,double *time1,double *ttime1,
             double *plicon1,ITG *nplicon1,double *plkcon1,ITG *nplkcon1,
             double *xstateini1,double *xstiff1,double *xstate1,
             ITG *npmat1_,char *matname1,ITG *mi1,ITG *ielas1,ITG *icmd1,
             ITG *ncmat1_,ITG *nstate1_,double *stiini1,double *vini1,
             double *ener1,double *enerini1,ITG *istep1,ITG *iinc1,
             double *springarea1,double *reltime1,ITG *calcul_qa1,
             ITG *iener1,ITG *ikin1,ITG *ne01,double *thicke1,
             double *emeini1,double *pslavsurf1,double *pmastsurf1,
             ITG *mortar1,double *clearini1,ITG *nea,ITG *neb,
             ITG *ielprop1,double *prop1,double *distmin1,ITG *ndesi1,
             ITG *nodedesi1,ITG *nobject1,double *g01,
             double *dgdx1,ITG *iobject1,double *sti1,double *xener1,
             ITG *istartdesi1,ITG *ialdesi1,double *xdesi1,ITG *idesvar));

void FORTRAN(objective_shapeener_tot,(ITG *ne,ITG *kon,ITG *ipkon,char *lakon,
             double *fint,double *vold,ITG *iperturb,ITG *mi,ITG *nactdof,
             double *dgdx,double *df,ITG *ndesi,ITG *iobject,ITG *jqs,
             ITG *irows,double *vec,ITG *iponk2dto3d));

void FORTRAN(objective_stress,(ITG *nodeset,ITG *istartset,ITG *iendset,
             ITG *ialset,ITG *nk,ITG *idesvar,ITG *iobject,ITG *mi,
             double *g0,ITG *nobject,double *stn,char *objectset,
             double *expks));

void FORTRAN(objective_stress_dx_dy,(ITG *nodeset,ITG *istartset,ITG *iendset,
             ITG *ialset,ITG *nk,ITG *idesvar1,ITG *idesvar2,ITG *iobject,
             double *dgdx,ITG *ndesi,ITG *nobject,double *stn,double *dstn1,
             double *dstn2,double *dstn12,char *objectset,double *g0,
             double *dgdxdy));

void FORTRAN(objective_stress_se,(ITG *nk,ITG *iobject,ITG *mi,double *dstn,
             char *objectset,ITG *ialnneigh,ITG *naneigh,ITG *nbneigh,
             double *stn,double *dksper));
              
void FORTRAN(objective_stress_tot,(double *dgdx,double *df,ITG *ndesi,
                                   ITG *iobject,ITG *jqs,ITG *irows,
                                   double *dgdu));

void FORTRAN(op,(ITG *n,double *x,double *y,double *ad,double *au,ITG *jq,ITG *irow));

void FORTRAN(opas,(ITG *n,double *x,double *y,double *ad,double *au,ITG *jq,
                   ITG *irow,ITG *nzs));

void FORTRAN(op_corio,(ITG *n,double *x,double *y,double *ad,double *au,
                       ITG *jq,ITG *irow));

void FORTRAN(openfile,(char *jobname,char *output));

void FORTRAN(openfilefluid,(char *jobname));

void FORTRAN(postprojectgrad,(ITG *ndesi,ITG *nodedesi,double *dgdxglob,
                              ITG *nactive,ITG *nobject,ITG *nnlconst,
                              ITG *ipoacti,ITG *nk,ITG *iconst,char *objectset,
                              ITG *iconstacti,ITG *inameacti));

void FORTRAN(posttransition,(double *dgdxglob,ITG *nobject,ITG *nk,
             ITG *nodedesi,ITG *ndesi,char *objectset));

void FORTRAN(postview,(ITG *ntr,char *sideload,ITG *nelemload,ITG *kontri,
             ITG *ntri,ITG *nloadtr,double *tenv,double *adview,double *auview,
             double *area,double *fenv,ITG *jqrad,ITG *irowrad,ITG *nzsrad));

void preparll(ITG *mt,double *dtime,double *veold,double *scal1,
                   double *accold,double *uam,ITG *nactdof,double *v,
                   double *vold,double *scal2,ITG *nk,ITG *num_cpus);

void *preparllmt(ITG *i);

void FORTRAN(precfd,(ITG *ne,ITG *ipkon,ITG *kon,char *lakon,ITG *ipnei,
                     ITG *neifa,ITG *neiel,ITG *ipoface,ITG *nodface,
                     ITG *ielfa,ITG *nkonnei,ITG *nface,ITG *ifaext,
                     ITG *nfaext,ITG *isolidsurf,ITG *nsolidsurf,char *set,
                     ITG *nset,ITG *istartset,ITG *iendset,ITG *ialset,
                     double *vel,double *vold,ITG *mi,ITG *neij,ITG *nef,
                     ITG *nactdoh,ITG *ipkonf,char *lakonf,ITG *ielmatf,
                     ITG *ielmat,ITG *ielorienf,ITG *ielorien,ITG *norien,
                     double *cs,ITG *mcs,char *tieset,double *x,double *y,
                     double *z,double *xo,double *yo,double *zo,ITG *nx,
                     ITG *ny,ITG *nz,double *co,ITG *ifatie,double *velo,
                     double *veloo,ITG *initial));

void FORTRAN(precfdcyc,(int *nelemface,char *sideface,int *nface,int *ipoface,
        int *nodface,int *ne,int *ipkon,int *kon,char *lakon,
        int *ikboun,int *ilboun,double *xboun,int *nboun,int *nk,
        int *isolidsurf,int *nsolidsurf,int *ifreestream,int *nfreestream,
        int *neighsolidsurf,int *iponoel,int *inoel,int *inoelfree,
        int *nef,double *co,int *ipompc,int *nodempc,int *ikmpc,
        int *ilmpc,int *nmpc,char *set,int *istartset,int *iendset,
        int *ialset,int *nset,int *iturbulent));

void precontact(ITG *ncont,ITG *ntie,char *tieset,ITG *nset,char *set,
        ITG *istartset,ITG *iendset,ITG *ialset,ITG *itietri,
        char *lakon,ITG *ipkon,ITG *kon,ITG *koncont,ITG *ne,
        double *cg,double *straight,double *co,double *vold,
        ITG *istep,ITG *iinc,ITG *iit,ITG *itiefac,
        ITG *islavsurf,ITG *islavnode,ITG *imastnode,
        ITG *nslavnode,ITG *nmastnode,ITG *imastop,ITG *mi,
        ITG *ipe,ITG *ime,double *tietol,ITG *iflagact,
        ITG *nintpoint,double **pslavsurfp,double *xmastnor,double *cs,
        ITG *mcs,ITG *ics,double *clearini,ITG *nslavs);

void FORTRAN(preconvert2slapcol,(ITG *irow,ITG *ia,
                              ITG *jq,ITG *ja,ITG *nzs,ITG *nef));

void FORTRAN(predgmres,(ITG *n,double *b,double *x,ITG *nelt,ITG *ia,ITG *ja,
             double *a,ITG *isym,ITG *itol,double *tol,ITG *itmax,ITG *iter,
             double *err,ITG *ierr,ITG *iunit,double *sb,double *sx,
             double *rgwk,ITG *lrgw,ITG *igwk,ITG *ligw,double *rwork,
             ITG *iwork));

void FORTRAN(predgmres_struct,(ITG *neq,double *b,double *sol,ITG *nelt,
                               ITG *irow,ITG *jq,double *au,ITG *isym,
                               ITG *itol,double *tol,ITG *itmax,ITG *iter,double *err,
                               ITG *ierr,ITG *iunit,double *sb,double *sx,
                               double *rgwk,ITG *lrgw,ITG *igwk,ITG *ligw,
                               double *rwork,ITG *iwork));

void predgmres_struct_mt(double *ad,double **au,double *adb,double *aub,
         double *sigma,double *b,ITG *icol,ITG *irow,
         ITG *neq,ITG *nzs,ITG *symmetryflag,ITG *inputformat,ITG *jq,
         ITG *nzs3,ITG *nrhs);

void prediction(double *uam,ITG *nmethod,double *bet,double *gam,double *dtime,
               ITG *ithermal,ITG *nk,double *veold,double *accold,double *v,
               ITG *iinc,ITG *idiscon,double *vold,ITG *nactdof,ITG *mi,
               ITG *num_cpus);

void prediction_em(double *uam,ITG *nmethod,double *bet,double *gam,double *dtime,
               ITG *ithermal,ITG *nk,double *veold,double *v,
               ITG *iinc,ITG *idiscon,double *vold,ITG *nactdof,ITG *mi);

void FORTRAN(prefilter,(double *co,ITG *nodedesi,ITG *ndesi,double *xo,
                        double *yo,double *zo,double *x,double *y,
                        double *z,ITG *nx,ITG *ny,ITG *nz));

void preiter(double *ad,double **aup,double *b,ITG **icolp,ITG **irowp,
             ITG *neq,ITG *nzs,ITG *isolver,ITG *iperturb);

void FORTRAN(preprojectgrad,(double *vector,ITG *ndesi,ITG *nodedesi,
                             double *dgdxglob,
                             ITG *nactive,ITG *nobject,ITG *nnlconst,
                             ITG *ipoacti,ITG *nk,double *rhs,ITG *iconst,
                             char *objectset,double *xtf));

void FORTRAN(presgradient,(int *iponoel,int *inoel,double *sa,double *sav,
            int *nk,double *dt,double *shockcoef,double *dtimef,int *ipkon,
            int *kon,char *lakon,double *vold,int *mi,int *compressible,
            int *nmethod,double *dtl,int *isolidsurf,int *nsolidsurf,
            double *co,int *euler));

void FORTRAN(prethickness,(double *co,double *xo,double *yo,
             double *zo,double *x,double *y,double *z,ITG *nx,ITG *ny,
             ITG *nz,ITG *ifree,ITG *nodedesiinv,ITG *ndesiboun,
             ITG *nodedesiboun,char *set,ITG *nset,char *objectset,
             ITG *iobject,ITG *istartset,ITG *iendset,ITG *ialset));       

void FORTRAN(pretransition,(ITG *ipkon,ITG *kon,char *lakon,double *co,
             ITG *nk,ITG *ipoface,ITG *nodface,ITG *nodedesiinv,double *xo,
             double *yo,double *zo,double *x,double *y,double *z,ITG *nx,
             ITG *ny,ITG *nz,ITG *ifree));

void FORTRAN(printoutfluid,(char *set,ITG *nset,ITG *istartset,ITG *iendset,
             ITG *ialset,ITG *nprint,char *prlab,char *prset,
             ITG *ipkon,
             char *lakon,double *stx,double *eei,double *xstate,
             double *ener,ITG *mi,ITG *nstate_,
             double *co,ITG *kon,double *qfx,double *ttime,double *trab,
             ITG *inotr,ITG *ntrans,double *orab,ITG *ielorien,
             ITG *norien,double *vold,
             ITG *ielmatf,double *thicke,double *eme,double *vcontu,
             double *physcon,ITG *nactdoh,ITG *ielpropf,double *prop,
             double *xkappa,double *xmach,ITG *ithermal,char *orname));

void FORTRAN(printoutfluidfem,(char *set,int *nset,int *istartset,int *iendset,
             int *ialset,int *nprint,char *prlab,char *prset,
             double *v,double *t1,double *fn,int *ipkon,
             char *lakon,double *stx,double *eei,double *xstate,
             double *ener,int *mi,int *nstate_,int *ithermal,
             double *co,int *kon,double *qfx,double *ttime,double *trab,
             int *inotr,int *ntrans,double *orab,int *ielorien,
             int *norien,int *nk,int *ne,int *inum,char *filab,double *vold,
             int *ikin,int *ielmat,double *thicke,double *eme,
             double *vcontu,double *physcon));

void FORTRAN(printoutface,(double *co,double *rhcon,ITG *nrhcon,ITG *ntmat_,
            double *vold,double *shcon,ITG *nshcon,double *cocon,
            ITG *ncocon,ITG *icompressible,ITG *istartset,ITG *iendset,
            ITG *ipkon,char *lakon,ITG *kon,ITG *ialset,char *prset,
            double *timef,ITG *nset,char *set,ITG *nprint,char *prlab,
            ITG *ielmat,ITG *mi,ITG *ithermal,ITG *nactdoh,ITG *icfd,
            double *time,double *stn));

void FORTRAN(projectgrad,(double *vector,ITG *ndesi,ITG *nodedesi,
                          double *dgdxglob,
                          ITG *nactive,ITG *nobject,ITG *nnlconst,
                          ITG *ipoacti,ITG *nk,double *rhs,ITG *iconst,
                          char *objectset,double *lambda,double *xtf,
                          double *objnorm));

void projectgradmain(ITG *nobject,char **objectsetp,double **dgdxglobp,
                     double *g0,ITG *ndesi,ITG *nodedesi,ITG *nk,ITG *isolver,
                     double *co,double *xdesi,double *distmin,ITG *nconstraint);

void FORTRAN(projectnodes,(ITG *nktet_,ITG *ipoed,ITG *iedgmid,
             ITG *iexternal,ITG *iedgext,double *cotet,ITG *nktet,
             ITG *iedg,ITG *iquad,ITG *integerglob,double *doubleglob,
             ITG *iexternfa,ITG *ifacext,ITG *itreated,
             ITG *ialsetexternel,ITG *nexternel,ITG *ilist,ITG *isharp,
             ITG *ipofa,ITG *ifac,ITG *iedgextfa,ITG *ifacexted,
             ITG *nexternedg));
                       
void FORTRAN(propertynet,(ITG *ieg,ITG *nflow,double *prop,ITG *ielprop,
                          char *lakon,ITG *iin,double *prop_store,
                          double *ttime,double *time,ITG *nam,char *amname,
                          ITG *namta,double *amta));

int pthread_create (pthread_t *thread_id,const pthread_attr_t *attributes,
                    void *(*thread_function)(void *),void *arguments);

int pthread_join (pthread_t thread,void **status_ptr);
          
void FORTRAN(quadraticsens,(ITG *ipkon,
             char *lakon,ITG *kon,ITG *nobject,double *dgdxglob,
             double *xinterpol,ITG *nnodes,ITG *ne,ITG *nk,
             ITG *nodedesiinv,char *objectset));

void radcyc(ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
            double *cs,ITG *mcs,ITG *nkon,ITG *ialset,ITG *istartset,
            ITG *iendset,ITG **kontrip,ITG *ntri,
            double **cop,double **voldp,ITG *ntrit,ITG *inocs,ITG *mi);

void radflowload(ITG *itg,ITG *ieg,ITG *ntg,ITG *ntr,double *adrad,
       double *aurad,
       double *bcr,ITG *ipivr,double *ac,double *bc,ITG *nload,
       char *sideload,ITG *nelemload,double *xloadact,char *lakon,ITG *ipiv,
       ITG *ntmat_,double *vold,double *shcon,ITG *nshcon,ITG *ipkon,
       ITG *kon,double *co,ITG *kontri,ITG *ntri,
       ITG *nloadtr,double *tarea,double *tenv,double *physcon,double *erad,
       double **adviewp,double **auviewp,
       ITG *nflow,ITG *ikboun,double *xboun,ITG *nboun,ITG *ithermal,
       ITG *iinc,ITG *iit,double *cs,ITG *mcs,ITG *inocs,ITG *ntrit,
       ITG *nk,double *fenv,ITG *istep,double *dtime,double *ttime,
       double *time,ITG *ilboun,ITG *ikforc,ITG *ilforc,double *xforc,
       ITG *nforc,double *cam,ITG *ielmat,ITG *nteq,double *prop,
       ITG *ielprop,ITG *nactdog,ITG *nacteq,ITG *nodeboun,ITG *ndirboun,
       ITG *network,double *rhcon,ITG *nrhcon,
       ITG *ipobody,ITG *ibody,double *xbody,ITG *nbody,ITG *iviewfile,
       char *jobnamef,double *ctrl,double *xloadold,double *reltime,
       ITG *nmethod,char *set,ITG *mi,ITG * istartset,ITG* iendset,
       ITG *ialset,ITG *nset,ITG *ineighe,ITG *nmpc,ITG *nodempc,
       ITG *ipompc,double *coefmpc,char *labmpc,ITG *iemchange,ITG *nam,
       ITG *iamload,ITG *jqrad,ITG *irowrad,ITG *nzsrad,ITG *icolrad,
       ITG *ne,ITG *iaxial,double *qa,double *cocon,ITG *ncocon,
       ITG *iponoel,ITG *inoel,ITG *nprop,char *amname,ITG *namta,
       double *amta);

void FORTRAN (radmatrix,(ITG *ntr,double *adrad,double *aurad,double *bcr,
       char *sideload,ITG *nelemload,double *xloadact,char *lakon,
       double *vold,ITG *ipkon,ITG *kon,double *co,ITG *nloadtr,
       double *tarea,double *tenv,double *physcon,double *erad,
       double *adview,double *auview,ITG *ithermal,ITG *iinc,
       ITG *iit,double *fenv,ITG *istep,
       double *dtime,double *ttime,double *time,ITG *iviewfile,
       double *xloadold,double *reltime,ITG *nmethod,
       ITG *mi,ITG *iemchange,ITG *nam,ITG *iamload,ITG *jqrad,
       ITG *irowrad,ITG *nzsrad));

void FORTRAN(radresult,(ITG *ntr,double *xloadact,double *bcr,
       ITG *nloadtr,double *tarea,double * tenv,double *physcon,double *erad,
       double *auview,double *fenv,ITG *irowrad,ITG *jqrad,
       ITG *nzsrad,double *q));

unsigned ITG randmt();

void randomfieldmain(ITG *kon,ITG *ipkon,char *lakon,ITG *ne,ITG *nmpc,
         ITG *nactdof,ITG *mi,ITG *nodedesi,ITG *ndesi,ITG *istartdesi,
         ITG *ialdesi,double *co,double *physcon,ITG *isolver,ITG *ntrans,
         ITG *nk,ITG *inotr,double *trab,char *jobnamec,ITG *nboun,
         double *cs,         
         ITG *mcs,ITG *inum,ITG *nmethod,ITG *kode,char *filab,ITG *nstate_,
         ITG *istep,char *description,char *set,ITG *nset,ITG *iendset,
         char *output,ITG *istartset,ITG *ialset,double *extnor);

void FORTRAN(randomval,(double *randval,ITG *nev));

void FORTRAN(readforce,(double *zc,ITG *neq,ITG *nev,ITG *nactdof,
             ITG *ikmpc,ITG *nmpc,ITG *ipompc,ITG *nodempc,ITG *mi,
             double *coefmpc,char *jobnamec,double *aa,
             ITG *igeneralizedforce));

void readinput(char *jobnamec,char **inpcp,ITG *nline,ITG *nset,ITG *ipoinp,
               ITG **inpp,ITG **ipoinpcp,ITG *ithermal,ITG *nuel_); 

void FORTRAN(readsen,(double *g0,double *dgdx,ITG *ndesi,ITG *nobject,
                       ITG *nodedesi,char *jobnamef));

void FORTRAN(readview,(ITG *ntr,double *adview,double *auview,double *fenv,
             ITG *nzsrad,ITG *ithermal,char *jobnamef));

void FORTRAN(rearrange,(double *au,ITG *irow,ITG *icol,ITG *ndim,ITG *neq));

void FORTRAN(rectcyl,(double *co,double *v,double *fn,double *stn,
                      double *qfn,double *een,double *cs,ITG *nk,
                      ITG *icntrl,double *t,char *filab,ITG *imag,
                      ITG *mi,double *emn));

void FORTRAN(rectcylexp,(double *co,double *v,double *fn,double *stn,
                      double *qfn,double *een,double *cs,ITG *nkt,
                      ITG *icntrl,double *t,char *filab,ITG *imag,ITG *mi,
                      ITG *iznode,ITG *nznode,ITG *nsectors,ITG *nk,
                      double *emn));

void FORTRAN(rectcyltrfm,(ITG *node,double *co,double *cs,ITG *cntrl,
             double *fin,double *fout));

void FORTRAN(rectcylvi,(double *co,double *v,double *fn,double *stn,
                      double *qfn,double *een,double *cs,ITG *nk,
                      ITG *icntrl,double *t,char *filab,ITG *imag,ITG *mi,
                      double *emn));

void FORTRAN(rectcylvold,(double *co,double *vold,double *cs,
                      ITG *icntrl,ITG *mi,
                      ITG *iznode,ITG *nznode,ITG *nsectors,ITG *nk));

void refinemesh(ITG *nk,ITG *ne,double *co,ITG *ipkon,ITG *kon,
                double *v,double *veold,double *stn,double *een,
                double *emn,double *epn,double *enern,double *qfn,
                double *errn,char *filab,ITG *mi,char *lakon,
                char *jobnamec);

void FORTRAN(reinit2,(ITG *kontet,ITG *ifac,ITG *ieln,ITG *netet_,
                      ITG *newsize,ITG *ifatet,ITG *itetfa,ITG *iedg,
                      ITG *ieled));

void remastruct(ITG *ipompc,double **coefmpcp,ITG **nodempcp,ITG *nmpc,
              ITG *mpcfree,ITG *nodeboun,ITG *ndirboun,ITG *nboun,
              ITG *ikmpc,ITG *ilmpc,ITG *ikboun,ITG *ilboun,
              char *labmpc,ITG *nk,
              ITG *memmpc_,ITG *icascade,ITG *maxlenmpc,
              ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
              ITG *nactdof,ITG *icol,ITG *jq,ITG **irowp,ITG *isolver,
              ITG *neq,ITG *nzs,ITG *nmethod,double **fp,
              double **fextp,double **bp,double **aux2p,double **finip,
              double **fextinip,double **adbp,double **aubp,ITG *ithermal,
              ITG *iperturb,ITG *mass,ITG *mi,ITG *iexpl,ITG *mortar,
              char *typeboun,double **cvp,double **cvinip,ITG *iit,
              ITG *network);

void remastructar(ITG *ipompc,double **coefmpcp,ITG **nodempcp,ITG *nmpc,
              ITG *mpcfree,ITG *nodeboun,ITG *ndirboun,ITG *nboun,
              ITG *ikmpc,ITG *ilmpc,ITG *ikboun,ITG *ilboun,
              char *labmpc,ITG *nk,
              ITG *memmpc_,ITG *icascade,ITG *maxlenmpc,
              ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
              ITG *nactdof,ITG *icol,ITG *jq,ITG **irowp,ITG *isolver,
              ITG *neq,ITG *nzs,ITG *nmethod,ITG *ithermal,
              ITG *iperturb,ITG *mass,ITG *mi,ITG *ics,double *cs,
              ITG *mcs,ITG *mortar,char *typeboun,ITG *iit,ITG *network,
              ITG *iexpl);

void remastructem(ITG *ipompc,double **coefmpcp,ITG **nodempcp,ITG *nmpc,
              ITG *mpcfree,ITG *nodeboun,ITG *ndirboun,ITG *nboun,
              ITG *ikmpc,ITG *ilmpc,ITG *ikboun,ITG *ilboun,
              char *labmpc,ITG *nk,
              ITG *memmpc_,ITG *icascade,ITG *maxlenmpc,
              ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
              ITG *nactdof,ITG *icol,ITG *jq,ITG **irowp,ITG *isolver,
              ITG *neq,ITG *nzs,ITG *nmethod,double **fp,
              double **fextp,double **bp,double **aux2p,double **finip,
              double **fextinip,double **adbp,double **aubp,ITG *ithermal,
              ITG *iperturb,ITG *mass,ITG *mi,ITG *ielmat,double *elcon,
              ITG *ncmat_,ITG *ntmat_,ITG *inomat,ITG *network);

void FORTRAN(renumber,(ITG *nnn,ITG *npn,ITG *adj,ITG *xadj,ITG *iw,
                       ITG *mmm,ITG *xnpn,ITG *inum1,ITG *ipnei,ITG *nef,
                       ITG *neiel,ITG *neifa,ITG *ifatie));

void renumbermain(ITG *nef,ITG *ipnei,ITG *neiel,ITG *ipkonf,ITG *ielmatf,
                  ITG *ielorienf,ITG *neifa,ITG *neij,ITG *nflnei,
                  ITG *nactdoh,ITG *nactdohinv,ITG *nface,ITG *ielfa,
                  ITG *mi,ITG *ifatie,ITG *norien,char *lakonf);

void FORTRAN(reorderlhs,(double *au,double *am,ITG *iamorig,ITG *nz_numa,
                 ITG *nz_numb));

void *reorderlhs1mt(ITG *i);

void reorderlhsmain(double *au,double *am,ITG *iamorig,ITG *nz_num,
    ITG *num_cpus);

void FORTRAN(reorderrhs,(double *am,double *b,double *vold,ITG *neighblock,
                    ITG *nneighblock));

void res1parll(ITG *mt,ITG *nactdof,double *aux2,double *vold,
                    double *vini,double *dtime,double *accold,
                    ITG *nk,ITG *num_cpus);

void *res1parllmt(ITG *i);

void res2parll(double *b,double *scal1,double *fext,double *f,
                    double *alpha,double *fextini,double *fini,
                    double *adb,double *aux2,ITG *neq0,ITG *num_cpus);

void *res2parllmt(ITG *i);

void res3parll(ITG *mt,ITG *nactdof,double *f,double *fn,
                    ITG *nk,ITG *num_cpus);

void *res3parllmt(ITG *i);

void res4parll(double *cv,double *alpham,double *adb,double *aux2,
                    double *b,double *scal1,double *alpha,double *cvini,
                    ITG *neq0,ITG *num_cpus);

void *res4parllmt(ITG *i);

void FORTRAN(restartshort,(ITG *nset,ITG *nload,ITG *nbody,ITG *nforc,
    ITG *nboun,
    ITG *nk,ITG *ne,ITG *nmpc,ITG *nalset,ITG *nmat,ITG *ntmat,ITG *npmat,
    ITG *norien,ITG *nam,ITG *nprint,ITG *mint,ITG *ntrans,ITG *ncs,
    ITG *namtot,ITG *ncmat,ITG *memmpc,ITG *ne1d,ITG *ne2d,ITG *nflow,
    char *set,ITG *meminset,ITG *rmeminset,char *jobnamec,ITG *irestartstep,
    ITG *icntrl,ITG *ithermal,ITG *nener,ITG *nstate_,ITG *ntie,ITG *nslavs,
    ITG *nkon,ITG *mcs,ITG *nprop,ITG *mortar,ITG *ifacecount,ITG *nintpoint,
    ITG *infree,ITG *nef));

void FORTRAN(restartwrite,(ITG *istep,ITG *nset,ITG*nload,ITG *nforc,
  ITG * nboun,ITG *nk,ITG *ne,ITG *nmpc,ITG *nalset,ITG *nmat,ITG *ntmat_,
  ITG *npmat_,ITG *norien,ITG *nam,ITG *nprint,ITG *mi,
  ITG *ntrans,ITG *ncs_,ITG *namtot,ITG *ncmat_,ITG *mpcend,
  ITG *maxlenmpc,ITG *ne1d,
  ITG *ne2d,ITG *nflow,ITG *nlabel,ITG *iplas,ITG *nkon,ITG *ithermal,
  ITG *nmethod,ITG *iperturb,ITG *nstate_,ITG *nener,char *set,
  ITG *istartset,ITG *iendset,ITG *ialset,double *co,ITG *kon,ITG *ipkon,
  char *lakon,ITG *nodeboun,ITG *ndirboun,ITG *iamboun,double *xboun,
  ITG *ikboun,ITG *ilboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
  char *labmpc,ITG *ikmpc,ITG *ilmpc,ITG *nodeforc,ITG *ndirforc,
  ITG *iamforc,double *xforc,ITG *ikforc,ITG *ilforc,ITG *nelemload,
  ITG *iamload,char *sideload,double *xload,
  double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,double *alcon,
  ITG *nalcon,double *alzero,double *plicon,ITG *nplicon,double *plkcon,
  ITG *nplkcon,char *orname,double *orab,ITG *ielorien,double *trab,
  ITG *inotr,char *amname,double *amta,ITG *namta,double *t0,double *t1,
  ITG *iamt1,double *veold,ITG *ielmat,char *matname,
  char *prlab,char *prset,char *filab,double *vold,
  ITG *nodebounold,ITG *ndirbounold,double *xbounold,double *xforcold,
  double *xloadold,double *t1old,double *eme,ITG *iponor,
  double *xnor,ITG *knor,double *thicke,double *offset,
  ITG *iponoel,ITG *inoel,ITG *rig,
  double *shcon,ITG *nshcon,double *cocon,ITG *ncocon,
  ITG *ics,double *sti,double *ener,double *xstate,
  char *jobnamec,ITG *infree,double *prestr,ITG *iprestr,
  char *cbody,ITG *ibody,double *xbody,ITG *nbody,double *xbodyold,
  double *ttime,double *qaold,double *cs,
  ITG *mcs,char *output,double *physcon,double *ctrl,char *typeboun,
  double *fmpc,char *tieset,ITG *ntie,double *tietol,ITG *nslavs,
  double *t0g,double *t1g,ITG *nprop,ITG *ielprop,double *prop,ITG *mortar,
  ITG *nintpoint,ITG *ifacecount,ITG *islavsurf,double *pslavsurf,
  double *clearini,ITG *irstrt,double *vel,ITG *nef,double *velo,
  double *veloo,ITG *ne2boun));

void FORTRAN(resultnet,(ITG *itg,ITG *ieg,ITG *ntg,
                        double *bc,ITG *nload,char *sideload,
                        ITG *nelemload,double *xloadact,char *lakon,
                        ITG *ntmat_,double *v,double *shcon,ITG *nshcon,
                        ITG *ipkon,ITG *kon,double *co,ITG *nflow,
                        ITG *iinc,ITG *istep,
                        double *dtime,double *ttime,double *time,
                        ITG *ikforc,ITG *ilforc,
                        double *xforcact,ITG *nforc,
                        ITG *ielmat,ITG *nteq,double *prop,
                        ITG *ielprop,ITG *nactdog,ITG *nacteq,ITG *iin,
                        double *physcon,double *camt,double *camf,
                        double *camp,double *rhcon,ITG *nrhcon,
                        ITG *ipobody,ITG *ibody,double *xbody,ITG *nbody,
                        double *dtheta,double *vold,double *xloadold,
                        double *reltime,ITG *nmethod,char *set,ITG *mi,
                        ITG *ineighe,double *cama,double *vamt,
                        double *vamf,double *vamp,double *vama,
                        ITG *nmpc,ITG *nodempc,ITG *ipompc,double *coefmpc,
                        char *labmpc,ITG *iaxial,double *qat,double *qaf,
                        double *ramt,double *ramf,double *ramp,
                        double *cocon,ITG *ncocon,ITG *iponoel,ITG *inoel,
                        ITG *iplausi));

void results(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
             ITG *ne,double *v,double *stn,ITG *inum,
             double *stx,
             double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
             double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
             ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
             double *t0,double *t1,ITG *ithermal,double *prestr,
             ITG *iprestr,char *filab,double *eme,double *emn,
             double *een,ITG *iperturb,double *f,double *fn,ITG *nactdof,
             ITG *iout,double *qa,
             double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
             double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,
             double *coefmpc,char *labmpc,ITG *nmpc,ITG *nmethod,
             double *vmax,ITG *neq,double *veold,double *accold,
             double *beta,double *gamma,double *dtime,double *time,
             double *ttime,double *plicon,
             ITG *nplicon,double *plkcon,ITG *nplkcon,
             double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
             double *epl,char *matname,ITG *mi,ITG *ielas,
             ITG *icmd,ITG *ncmat_,ITG *nstate_,double *stiini,
             double *vini,ITG *ikboun,ITG *ilboun,double *ener,
             double *enern,double *emeini,double *xstaten,double *eei,
             double *enerini,double *cocon,ITG *ncocon,char *set,
             ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,double *qfx,double *qfn,double *trab,
             ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,
             ITG *nload,ITG *ikmpc,ITG *ilmpc,ITG *istep,ITG *iinc,
             double *springarea,double *reltime,ITG *ne0,double *thicke,
             double *shcon,ITG *nshcon,char *sideload,double *xload,
             double *xloadold,ITG *icfd,ITG *inomat,double *pslavsurf,
             double *pmastsurf,ITG *mortar,ITG *islavact,double *cdn,
             ITG *islavnode,ITG *nslavnode,ITG *ntie,double *clearini,
             ITG *islavsurf,ITG *ielprop,double *prop,double *energyini,
             double *energy,ITG *kscale,ITG *iponoel,ITG *inoel,ITG *nener,
             char *orname,ITG *network,ITG *ipobody,double *xbodyact,
             ITG *ibody,char *typeboun,ITG *itiefac,char *tieset,
             double *smscale,ITG *mscalmethod);

void FORTRAN(resultsem,(double *co,ITG *kon,ITG *ipkon,char *lakon,
             double *v,double *elcon,ITG *nelcon,ITG *ielmat,ITG *ntmat_,
             double *vold,double *dtime,char *matname,ITG *mi,ITG *ncmat_,
             ITG *nea,ITG *neb,double *sti,double *alcon,
             ITG *nalcon,double *h0,ITG *istartset,ITG *iendset,ITG *ialset,
             ITG *iactive,double *fn));

void *resultsemmt(ITG *i);

void resultsforc(ITG *nk,double *f,double *fn,ITG *nactdof,
       ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,ITG *nmpc,
       ITG *mi,double *fmpc,ITG *calcul_fn,ITG *calcul_f,ITG *num_cpus);

void  FORTRAN(resultsforc_em,(ITG *nk,double *f,double *fn,ITG *nactdof,
       ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,ITG *nmpc,
       ITG *mi,double *fmpc,ITG *calcul_fn,ITG *calcul_f,ITG *inomat));
       
void  FORTRAN(resultsforc_se,(ITG *nk,double *dfn,ITG *nactdof,
       ITG *ipompc,ITG *nodempc,double *coefmpc,ITG *nmpc,
       ITG *mi,double *fmpc,ITG *calcul_fn,ITG *calcul_f,ITG *idesvar,
       double *df,ITG *jqs,ITG *irows,double *distmin));

void resultsini(ITG *nk,double *v,ITG *ithermal,char *filab,
       ITG *iperturb,double *f,double *fn,ITG *nactdof,ITG *iout,
       double *qa,double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
       double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
       char *labmpc,ITG *nmpc,ITG *nmethod,double *cam,ITG *neq,
       double *veold,double *accold,double *bet,double *gam,double *dtime,
       ITG *mi,double *vini,ITG *nprint,char *prlab,ITG *intpointvar,
       ITG *calcul_fn,ITG *calcul_f,ITG *calcul_qa,ITG *calcul_cauchy,
       ITG *ikin,ITG *intpointvart,char *typeboun,ITG *num_cpus);

void FORTRAN(resultsini_em,(ITG *nk,double *v,ITG *ithermal,char *filab,
       ITG *iperturb,double *f,double *fn,ITG *nactdof,ITG *iout,
       double *qa,double *b,ITG *nodeboun,ITG *ndirboun,
       double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
       char *labmpc,ITG *nmpc,ITG *nmethod,double *cam,ITG *neq,
       double *veold,double *dtime,
       ITG *mi,double *vini,ITG *nprint,char *prlab,ITG *intpointvar,
       ITG *calcul_fn,ITG *calcul_f,ITG *calcul_qa,ITG *calcul_cauchy,
       ITG *iener,ITG *ikin,ITG *intpointvart,double *xforc,ITG *nforc));

void FORTRAN(resultsk,(int *nk,int *nactdoh,double *vtu,double *solk,
      double *solt,int *ipompc,int *nodempc,double *coefmpc,int *nmpc));

void FORTRAN(resultsmech,(double *co,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
          double *v,double *stx,double *elcon,ITG *nelcon,double *rhcon,
          ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
          ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,double *t0,
          double *t1,ITG *ithermal,double *prestr,ITG *iprestr,double *eme,
          ITG *iperturb,double *fn,ITG *iout,double *qa,double *vold,
          ITG *nmethod,double *veold,double *dtime,double *time,
          double *ttime,double *plicon,ITG *nplicon,double *plkcon,
          ITG *nplkcon,double *xstateini,double *xstiff,double *xstate,
          ITG *npmat_,char *matname,ITG *mi,ITG *ielas,ITG *icmd,ITG *ncmat_,
          ITG *nstate_,double *stiini,double *vini,double *ener,double *eei,
          double *enerini,ITG *istep,ITG *iinc,double *springarea,
          double *reltime,ITG *calcul_fn,ITG *calcul_qa,ITG *calcul_cauchy,
          ITG *iener,ITG *ikin,ITG *nal,ITG *ne0,double *thicke,
          double *emeini,double *pslavsurf,double *pmastsurf,ITG *mortar,
          double *clearini,ITG *nea,ITG *neb,ITG *ielprop,double *prop,
          ITG *kscale,ITG *list,ITG *ilist,double *smscale,ITG *mscalmethod,
          double *enerscal));

void *resultsmechmt(ITG *i);

void *resultsmechmtstr(ITG *i);

void *resultsmechmt_se(ITG *i);

void FORTRAN(resultsmech_se,(double *co,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
          double *v,double *stx,double *elcon,ITG *nelcon,double *rhcon,
          ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
          ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,double *t0,
          double *t1,ITG *ithermal,double *prestr,ITG *iprestr,double *eme,
          ITG *iperturb,double *fn,ITG *iout,double *vold,
          ITG *nmethod,double *veold,double *dtime,double *time,
          double *ttime,double *plicon,ITG *nplicon,double *plkcon,
          ITG *nplkcon,double *xstateini,double *xstiff,double *xstate,
          ITG *npmat_,char *matname,ITG *mi,ITG *ielas,ITG *icmd,ITG *ncmat_,
          ITG *nstate_,double *stiini,double *vini,double *ener,double *eei,
          double *enerini,ITG *istep,ITG *iinc,double *springarea,
          double *reltime,ITG *calcul_fn,ITG *calcul_cauchy,
          ITG *iener,ITG *ikin,ITG *ne0,double *thicke,
          double *emeini,double *pslavsurf,double *pmastsurf,ITG *mortar,
          double *clearini,ITG *nea,ITG *neb,ITG *ielprop,double *prop,
          double *dfn,ITG *idesvar,ITG *nodedesi,
          double *fn0,double *sti,ITG *icoordinate,
          double *dxstiff,ITG *ialdesi,double *xdesi));

void FORTRAN(resultsnoddir,(ITG *nk,double *v,ITG *nactdof,double *b,
       ITG *ipompc,ITG *nodempc,double *coefmpc,ITG *nmpc,ITG *mi));

void FORTRAN(resultsp,(int *nk,int *nactdoh,double *v,double *sol,
             int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *mi));

void  FORTRAN(resultsprint,(double *co,ITG *nk,ITG *kon,ITG *ipkon,
       char *lakon,ITG *ne,double *v,double *stn,ITG *inum,double *stx,
       ITG *ielorien,ITG *norien,double *orab,double *t1,ITG *ithermal,
       char *filab,double *een,ITG *iperturb,double *fn,ITG *nactdof,
       ITG *iout,double *vold,ITG *nodeboun,ITG *ndirboun,ITG *nboun,
       ITG *nmethod,double *ttime,double *xstate,double *epn,ITG *mi,
       ITG *nstate_,double *ener,double *enern,double *xstaten,double *eei,
       char *set,ITG *nset,ITG *istartset,ITG *iendset,ITG *ialset,ITG *nprint,
       char *prlab,char *prset,double *qfx,double *qfn,double *trab,ITG *inotr,
       ITG *ntrans,ITG *nelemload,ITG *nload,ITG *ikin,ITG *ielmat,
       double *thicke,double *eme,double *emn,double *rhcon,ITG *nrhcon,
       double *shcon,ITG *nshcon,double *cocon,ITG *ncocon,ITG *ntmat_,
       char *sideload,ITG *icfd,ITG *inomat,double *pslavsurf,
       ITG *islavact,double *cdn,ITG *mortar,ITG *islavnode,ITG *nslavnode,
       ITG *ntie,ITG *islavsurf,double *time,ITG *ielprop,double *prop,
       double *veold,ITG *ne0,ITG *nmpc,ITG *ipompc,ITG *nodempc,
       char *labmpc,double *energyini,double *energy,char *orname,
       double *xload,ITG *itiefac,double *pmastsurf,double *springarea,
       char *tieset));

void resultsstr(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
             ITG *ne,double *v,double *stn,ITG *inum,
             double *stx,
             double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
             double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
             ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
             double *t0,double *t1,ITG *ithermal,double *prestr,
             ITG *iprestr,char *filab,double *eme,double *emn,
             double *een,ITG *iperturb,double *f,double *fn,ITG *nactdof,
             ITG *iout,double *qa,
             double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
             double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,
             double *coefmpc,char *labmpc,ITG *nmpc,ITG *nmethod,
             double *vmax,ITG *neq,double *veold,double *accold,
             double *beta,double *gamma,double *dtime,double *time,
             double *ttime,double *plicon,
             ITG *nplicon,double *plkcon,ITG *nplkcon,
             double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
             double *epl,char *matname,ITG *mi,ITG *ielas,
             ITG *icmd,ITG *ncmat_,ITG *nstate_,double *stiini,
             double *vini,ITG *ikboun,ITG *ilboun,double *ener,
             double *enern,double *emeini,double *xstaten,double *eei,
             double *enerini,double *cocon,ITG *ncocon,char *set,
             ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,double *qfx,double *qfn,double *trab,
             ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,
             ITG *nload,ITG *ikmpc,ITG *ilmpc,ITG *istep,ITG *iinc,
             double *springarea,double *reltime,ITG *ne0,double *xforc,
             ITG *nforc,double *thicke,
             double *shcon,ITG *nshcon,char *sideload,double *xload,
             double *xloadold,ITG *icfd,ITG *inomat,double *pslavsurf,
             double *pmastsurf,ITG *mortar,ITG *islavact,double *cdn,
             ITG *islavnode,ITG *nslavnode,ITG *ntie,double *clearini,
             ITG *islavsurf,ITG *ielprop,double *prop,double *energyini,
             double *energy,ITG *kscale,ITG *nener,
             char *orname,ITG *network,ITG *neapar,ITG *nebpar);
             
void results_se(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
             ITG *ne,double *v,double *stn,ITG *inum,
             double *stx,
             double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
             double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
             ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
             double *t0,double *t1,ITG *ithermal,double *prestr,
             ITG *iprestr,char *filab,double *eme,double *emn,
             double *een,ITG *iperturb,double *f,double *fn,ITG *nactdof,
             ITG *iout,double *qa,
             double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
             double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,
             double *coefmpc,char *labmpc,ITG *nmpc,ITG *nmethod,
             double *vmax,ITG *neq,double *veold,double *accold,
             double *beta,double *gamma,double *dtime,double *time,
             double *ttime,double *plicon,
             ITG *nplicon,double *plkcon,ITG *nplkcon,
             double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
             double *epl,char *matname,ITG *mi,ITG *ielas,
             ITG *icmd,ITG *ncmat_,ITG *nstate_,double *stiini,
             double *vini,ITG *ikboun,ITG *ilboun,double *ener,
             double *enern,double *emeini,double *xstaten,double *eei,
             double *enerini,double *cocon,ITG *ncocon,char *set,
             ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,double *qfx,double *qfn,double *trab,
             ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,
             ITG *nload,ITG *ikmpc,ITG *ilmpc,ITG *istep,ITG *iinc,
             double *springarea,double *reltime,ITG *ne0,double *xforc,
             ITG *nforc,double *thicke,
             double *shcon,ITG *nshcon,char *sideload,double *xload,
             double *xloadold,ITG *icfd,ITG *inomat,double *pslavsurf,
             double *pmastsurf,ITG *mortar,ITG *islavact,double *cdn,
             ITG *islavnode,ITG *nslavnode,ITG *ntie,double *clearini,
             ITG *islavsurf,ITG *ielprop,double *prop,double *energyini,
             double *energy,double *df,double *distmin,
             ITG *ndesi,ITG *nodedesi,double *sti,
             ITG *nkon,ITG *jqs,ITG *irows,ITG *nactdofinv,
             ITG *icoordinate,double *dxstiff,ITG *istartdesi,
             ITG *ialdesi,double *xdesi,ITG *ieigenfrequency,
             double *fint,ITG *ishapeenergy,char *typeboun);

void FORTRAN(resultst,(int *nk,int *nactdoh,double *v,double *sol,
             int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *mi));

void FORTRAN(resultstherm,(double *co,ITG *kon,ITG *ipkon,
       char *lakon,double *v,
       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,ITG *ielmat,
       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,double *t0,
       ITG *iperturb,double *fn,double *shcon,ITG *nshcon,ITG *iout,
       double *qa,double *vold,ITG *ipompc,ITG *nodempc,
       double *coefmpc,ITG *nmpc,double *dtime,
       double *time,double *ttime,double *plkcon,ITG *nplkcon,double *xstateini,
       double *xstiff,double *xstate,ITG *npmat_,char *matname,
       ITG *mi,ITG *ncmat_,ITG *nstate_,double *cocon,ITG *ncocon,
       double *qfx,ITG *ikmpc,ITG *ilmpc,ITG *istep,
       ITG *iinc,double *springarea,ITG *calcul_fn,ITG *calcul_qa,ITG *nal,
       ITG *nea,ITG *neb,ITG *ithermal,ITG *nelemload,ITG *nload,
       ITG *nmethod,double *reltime,char *sideload,double *xload,
       double *xloadold,double *pslavsurf,double *pmastsurf,ITG *mortar,
       double *clearini,double *plicon,ITG *nplicon,ITG *ielprop,
       double *prop,ITG *iponoel,ITG *inoel,ITG *network,ITG *ipobody,
       double *xbodyact,ITG *ibody));

void *resultsthermemmt(ITG *i);

void *resultsthermmt(ITG *i);

void *resultsthermmt_se(ITG *i);

void resultsinduction(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
             ITG *ne,double *v,double *stn,ITG *inum,
             double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
             double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
             ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
             double *t0,double *t1,ITG *ithermal,double *prestr,
             ITG *iprestr,char *filab,double *eme,double *emn,
             double *een,ITG *iperturb,double *f,double *fn,ITG *nactdof,
             ITG *iout,double *qa,
             double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
             double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,
             double *coefmpc,char *labmpc,ITG *nmpc,ITG *nmethod,
             double *vmax,ITG *neq,double *veold,double *accold,
             double *beta,double *gamma,double *dtime,double *time,
             double *ttime,double *plicon,
             ITG *nplicon,double *plkcon,ITG *nplkcon,
             double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
             double *epl,char *matname,ITG *mi,ITG *ielas,
             ITG *icmd,ITG *ncmat_,ITG *nstate_,double *sti,
             double *vini,ITG *ikboun,ITG *ilboun,double *ener,
             double *enern,double *emeini,double *xstaten,double *eei,
             double *enerini,double *cocon,ITG *ncocon,char *set,
             ITG *nset,ITG *istartset,
             ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
             char *prset,double *qfx,double *qfn,double *trab,
             ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,
             ITG *nload,ITG *ikmpc,ITG *ilmpc,ITG *istep,ITG *iinc,
             double *springarea,double *reltime,ITG *ne0,double *xforc,
             ITG *nforc,double *thicke,
             double *shcon,ITG *nshcon,char *sideload,double *xload,
             double *xloadold,ITG *icfd,ITG *inomat,double *h0,
             ITG *islavnode,ITG *nslavnode,ITG *ntie,ITG *ielprop,
             double *prop,ITG *iactive,double *energyini,double *energy,
             ITG *iponoel,ITG *inoel,char *orname,ITG *network,
             ITG *ipobody,double *xbody,ITG *ibody);

void FORTRAN(resultsv1,(int *nk,int *nactdoh,double *v,double *sol,
            int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *mi));

void FORTRAN(resultsv2,(int *nk,int *nactdoh,double *v,double *sol,
             int *ipompc,int *nodempc,double *coefmpc,int *nmpc,int *mi));

void FORTRAN(rhs,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
               ITG *ne,ITG *ipompc,ITG *nodempc,double *coefmpc,
               ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
               double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
               double *xload,ITG *nload,double *xbody,ITG *ipobody,
               ITG *nbody,double *cgr,double *bb,ITG *nactdof,ITG *neq,
               ITG *nmethod,ITG *ikmpc,ITG *ilmpc,
               double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
               double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
               ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
               double *t0,double *t1,ITG *ithermal,
               ITG *iprestr,double *vold,ITG *iperturb,ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               ITG *npmat_,double *ttime,double *time,ITG *istep,
               ITG *iinc,double *dtime,double *physcon,ITG *ibody,
               double *xbodyold,double *reltime,double *veold,
               char *matname,ITG *mi,ITG *ikactmech,ITG *nactmech,
               ITG *ielprop,double *prop,double *sti,double *xstateini,
               double *xstate,ITG *nstate_,ITG *ntrans,ITG *inotr,
               double *trab,ITG *nea,ITG *neb));

void FORTRAN(rhsnodef,(double *co,ITG *kon,ITG *ne,ITG *ipompc,ITG *nodempc,
             double *coefmpc,ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
             double *xforc,
             ITG *nforc,double *fext,ITG *nactdof,ITG *nmethod,ITG *ikmpc,
             ITG *ntmat_,ITG *iperturb,ITG *mi,ITG *ikactmech,ITG *nactmech,
             ITG *ntrans,ITG *inotr,double *trab));

void rhsmain(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
             ITG *ipompc,ITG *nodempc,double *coefmpc,ITG *nmpc,
             ITG *nodeforc,ITG *ndirforc,double *xforcact,ITG *nforc,
             ITG *nelemload,char *sideload,double *xloadact,ITG *nload,
             double *xbodyact,ITG *ipobody,ITG *nbody,double *cgr,
             double *fext,
             ITG *nactdof,ITG *neq,ITG *nmethod,ITG *ikmpc,ITG *ilmpc,
             double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
             double *alcon,
             ITG *nalcon,double *alzero,ITG *ielmat,ITG *ielorien,
             ITG *norien,
             double *orab,ITG *ntmat_,double *t0,double *t1act,
             ITG *ithermal,
             ITG *iprestr,double *vold,ITG *iperturb,ITG *iexpl,
             double *plicon,
             ITG *nplicon,double *plkcon,ITG *nplkcon,ITG *npmat_,
             double *ttime,
             double *time,ITG *istep,ITG *iinc,double *dtime,
             double *physcon,
             ITG *ibody,double *xbodyold,double *reltime,double *veold,
             char *matname,
             ITG *mi,ITG *ikactmech,ITG *nactmech,ITG *ielprop,
             double *prop,
             double *sti,double *xstateini,double *xstate,ITG *nstate_,
             ITG *ntrans,ITG *inotr,double *trab);

void *rhsmt(ITG *i);

void FORTRAN(rhsp,(ITG *ne,char *lakon,ITG *ipnei,
             ITG *neifa,ITG *neiel,double *vfa,double *area,double *adfa,
             double *xlet,double *cosa,double *volume,double *au,double *ad,
             ITG *jq,ITG *irow,double *ap,ITG *ielfa,ITG *ifabou,
             double *xle,double *b,double *xxn,
             ITG *neq,ITG *nzs,double *hfa,double *gradpel,double *bp,
             double *xxi,ITG *neij,double *xlen,ITG *nefa,ITG *nefb,
             double *xxicn,double *flux,double *xxnj,double *gradpcfa,
             double *cosb));

void FORTRAN(rhspcomp,(ITG *ne,char *lakonf,ITG *ipnei,
             ITG *neifa,ITG *neiel,double *vfa,double *area,double *adfa,
             double *xlet,double *cosa,double *volume,double *au,double *ad,
             ITG *jq,ITG *irow,double *ap,ITG *ielfa,ITG *ifabou,
             double *xle,double *b,double *xxn,ITG *neq,
             ITG *nzs,double *hfa,double *gradpel,
             double *bp,double *xxi,ITG *neij,double *xlen,double *cosb,
             ITG *ielmatf,ITG *mi,double *a1,double *a2,double *a3,double *velo,
             double *veloo,double *dtimef,double *shcon,ITG *ntmat_,double *vel,
             ITG *nactdohinv,double *xrlfa,double *flux,ITG *nefa,ITG *nefb,
             ITG *iau6,double *xxicn,double *gamma,double *xxnj,
             double *gradpcfa));

void rhspcompmain(ITG *ne,char *lakonf,ITG *ipnei,
             ITG *neifa,ITG *neiel,double *vfa,double *area,double *adfa,
             double *xlet,double *cosa,double *volume,double *au,double *ad,
             ITG *jq,ITG *irow,double *ap,ITG *ielfa,ITG *ifabou,
             double *xle,double *b,double *xxn,ITG *neq,
             ITG *nzs,double *hfa,double *gradpel,
             double *bp,double *xxi,ITG *neij,double *xlen,double *cosb,
             ITG *ielmatf,ITG *mi,double *a1,double *a2,double *a3,double *velo,
             double *veloo,double *dtimef,double *shcon,ITG *ntmat_,double *vel,
             ITG *nactdohinv,double *xrlfa,double *flux,ITG *iau6,
                  double *xxicn,double *gamma,double *xxnj,double *gradpcfa);

void *rhspcompmt(ITG *i);

void rhspmain(ITG *ne,char *lakon,ITG *ipnei,
             ITG *neifa,ITG *neiel,double *vfa,double *area,double *adfa,
             double *xlet,double *cosa,double *volume,double *au,double *ad,
             ITG *jq,ITG *irow,double *ap,ITG *ielfa,ITG *ifabou,
             double *xle,double *b,double *xxn,
             ITG *neq,ITG *nzs,double *hfa,double *gradpel,double *bp,
             double *xxi,ITG *neij,double *xlen,ITG *iatleastonepressurebc,
             double *xxicn,double *flux,double *xxnj,double *gradpcfa,
             double *cosb);

void *rhspmt(ITG *i);

void FORTRAN(rotationvector,(double *a,double *euler));

void FORTRAN(rotationvectorinv,(double *a,double* euler));
       
void FORTRAN(scalesen,(double *dgdxglob,ITG *nobject,ITG *nk,ITG *nodedesi,
                       ITG *ndesi,char *objectset,ITG *iscaleflag,
                       ITG *ipoacti,ITG *nnlconst,ITG *nactive));

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
             ITG *inotr,ITG *ntrans,double *fmpc,char *cbody,ITG *ibody,
             double *xbody,ITG *nbody,double *xbodyold,double *timepar,
             double *thicke,char *jobnamec,char *tieset,ITG *ntie,
             ITG *istep,ITG *nmat,ITG *ielprop,double *prop,char *typeboun,
             ITG *mortar,ITG *mpcinfo,double *tietol,ITG *ics,
             ITG *icontact,ITG *nobject,char **objectsetp,ITG *istat,
             char *orname,ITG *nzsfreq,ITG *nlabel,double *physcon,
             char *jobnamef,ITG *iponor2d,ITG *knor2d,ITG *ne2d,
             ITG *iponoel2d,ITG *inoel2d,ITG *mpcend);

void FORTRAN(sensitivity_glob,(double *dgdxtot,double *dgdxtotglob,ITG *nobject,
             ITG *ndesi,ITG *nodedesi,ITG *nk));

void sensitivity_out(char *jobnamec,double *dgdxtotglob,ITG *neq,
                     ITG *nobject,double *g0);

void setpardou(double *var1,double var2,ITG size,ITG num_cpus);

void *setpardoumt(ITG *i);

void setparitg(ITG *iva1,ITG iva2,ITG size,ITG num_cpus);

void *setparitgmt(ITG *i);
                  
void sgenrand(unsigned long seed);

void FORTRAN(shape3tri,(double *xi,double *et,double *xl,double *xsj,
                      double *xs,double *shp,ITG *iflag));

void FORTRAN(shape4q,(double *xi,double *et,double *xl,double *xsj,
                      double *xs,double *shp,ITG *iflag));

void FORTRAN(shape4tet,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,ITG *iflag));

void FORTRAN(shape6tri,(double *xi,double *et,double *xl,double *xsj,
                      double *xs,double *shp,ITG *iflag));

void FORTRAN(shape6w,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,ITG *iflag));

void FORTRAN(shape8h,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,ITG *iflag));

void FORTRAN(shape8q,(double *xi,double *et,double *xl,double *xsj,
                      double *xs,double *shp,ITG *iflag));

void FORTRAN(shape10tet,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,ITG *iflag));

void FORTRAN(shape15w,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,ITG *iflag));

void FORTRAN(shape20h,(double *xi,double *et,double *ze,double *xl,
             double *xsj,double *shp,ITG *iflag));
       
void FORTRAN(slavintpoints,(ITG *ntie,ITG *itietri,ITG *ipkon,
        ITG *kon,char *lakon,double *straight,
        ITG *nintpoint,ITG *koncont,double *co,double *vold,double *xo,
        double *yo,double *zo,double *x,double *y,double *z,ITG *nx,
        ITG *ny,ITG *nz,ITG *islavsurf,
        ITG *islavnode,ITG *nslavnode,ITG *imastop,
        ITG *mi,ITG *ncont,ITG *ipe,ITG *ime,double *pslavsurf,
        ITG *i,ITG *l,ITG *ntri));

void FORTRAN(smalldist,(double *co,double *distmin,char *lakon,
             ITG *ipkon,ITG *kon,ITG *ne));

void FORTRAN(smoothshock,(double *adbv,double *aubv,double *adl,
            double *addiv,double *sol,double *aux,int *icolv,int *irowv,
            int *jqv,int *neqv,int *nzlv,double *sa));

void FORTRAN(solveeq,(double *adbv,double *aubv,double *adl,double *addiv,
            double *b,double *sol,double *aux,int *icolv,int *irowv,
            int *jqv,int *neqv,int *nzsv,int *nzlv));

void FORTRAN(solveexplicitly,(ITG *nef,double *vel,double *bv,double *auv,
             ITG *ipnei,ITG *neiel,ITG *nflnei));

void FORTRAN(sortev,(ITG *nev,ITG *nmd,double *eigxx,ITG *cyclicsymmetry,
                     double *xx,double *eigxr,ITG *pev,
                     ITG *istartnmd,ITG *iendnmd,double *aa,double *bb,
                     ITG *nevcomplex));

void FORTRAN(spcmatch,(double *xboun,ITG *nodeboun,ITG *ndirboun,ITG *nboun,
               double *xbounold,ITG *nodebounold,ITG *ndirbounold,
               ITG *nbounold,ITG *ikboun,ITG *ilboun,double *vold,
	       double *reorder,ITG *nreorder,ITG *mi,char *typeboun));

void FORTRAN(splitline,(char *text,char *textpart,ITG *n));

void spooles(double *ad,double *au,double *adb,double *aub,
             double *sigma,double *b,
             ITG *icol,ITG *irow,ITG *neq,ITG *nzs,ITG *symmtryflag,
             ITG *inputformat,ITG *nzs3);

void FORTRAN(springforc_n2f,(double *xl,ITG *konl,double *vl,ITG *imat,
             double *elcon,ITG *nelcon,double *elas,double *fnl,ITG *ncmat_,
             ITG *ntmat_,ITG *nope,char *lakonl,double *t1l,ITG *kode,
             double *elconloc,double *plicon,ITG *nplicon,ITG *npmat_,
             double *senergy,ITG *iener,double *cstr,ITG *mi,
             double *springarea,ITG *nmethod,ITG *ne0,ITG *nstate_,
             double *xstateini,double *xstate,double *reltime,
             ITG *ielas,double *venergy,ITG *ielorien,double *orab,
             ITG *norien,ITG *nelem));

void FORTRAN(springstiff_n2f,(double *xl,double *elas,ITG *konl,double *voldl,
             double *s,ITG *imat,double *elcon,ITG *nelcon,ITG *ncmat_,
             ITG *ntmat_,ITG *nope,char *lakonl,double *t1l,ITG *kode,
             double *elconloc,double *plicon,ITG *nplicon,ITG *npmat_,
             ITG *iperturb,double *springarea,ITG *nmethod,ITG *mi,ITG *ne0,
             ITG *nstate_,double *xstateini,double *xstate,double *reltime,
             ITG *nasym,ITG *ielorien,double *orab,ITG *norien,ITG *nelem));

void steadystate(double **co,ITG *nk,ITG **kon,ITG **ipkon,char **lakon,ITG *ne,
          ITG **nodeboun,ITG **ndirboun,double **xboun,ITG *nboun,
          ITG **ipompcp,ITG **nodempcp,double **coefmpcp,char **labmpcp,ITG *nmpc,
          ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
          ITG *nelemload,char *sideload,double *xload,
          ITG *nload,
          ITG **nactdof,ITG *neq,ITG *nzl,ITG *icol,ITG *irow,
          ITG *nmethod,ITG **ikmpcp,ITG **ilmpcp,ITG **ikboun,
          ITG **ilboun,
          double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
          double *cocon,ITG *ncocon,
          double *alcon,ITG *nalcon,double *alzero,ITG **ielmat,
          ITG **ielorien,ITG *norien,double *orab,ITG *ntmat_,
          double **t0,
          double **t1,ITG *ithermal,double *prestr,ITG *iprestr,
          double **voldp,ITG *iperturb,double *sti,ITG *nzs,
          double *timepar,double *xmodal,
          double **veoldp,char *amname,double *amta,
          ITG *namta,ITG *nam,ITG *iamforc,ITG *iamload,
          ITG **iamt1,ITG *jout,ITG *kode,char *filab,
          double **emep,double *xforcold,double *xloadold,
          double **t1old,ITG **iamboun,
          double **xbounold,ITG *iexpl,double *plicon,ITG *nplicon,
          double *plkcon,ITG *nplkcon,
          double *xstate,ITG *npmat_,char *matname,ITG *mi,
          ITG *ncmat_,ITG *nstate_,double **enerp,char *jobnamec,
          double *ttime,char *set,ITG *nset,ITG *istartset,
          ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
          char *prset,ITG *nener,double *trab,
          ITG **inotr,ITG *ntrans,double **fmpcp,char *cbody,ITG *ibody,
          double *xbody,ITG *nbody,double *xbodyold,ITG *istep,
          ITG *isolver,ITG *jq,char *output,ITG *mcs,ITG *nkon,
          ITG *ics,double *cs,ITG *mpcend,double *ctrl,
          ITG *ikforc,ITG *ilforc,double *thicke,ITG *nmat,
          char *typeboun,ITG *ielprop,double *prop,char *orname,
          ITG *ndamp,double *dacon);

void FORTRAN(stop,());

void FORTRAN(stopwithout201,());

void storecontactdof(ITG *nope,ITG *nactdof,ITG *mt,ITG *konl,
          ITG **ikactcontp,
          ITG *nactcont,ITG *nactcont_,double *bcont,double *fnl,
          ITG *ikmpc,ITG *nmpc,ITG *ilmpc,ITG *ipompc,ITG *nodempc,
          double *coefmpc);

void FORTRAN(storeresidual,(ITG *nactdof,double *b,double *fn,char *filab,
             ITG *ithermal,ITG *nk,double *sti,double *stn,
             ITG *ipkon,ITG *inum,ITG *kon,char *lakon,
             ITG *ne,ITG *mi,double *orab,ITG *ielorien,
             double *co,ITG *itg,ITG *ntg,double *vold,
             ITG *ielmat,double *thicke,ITG *ielprop,double *prop));
     
void FORTRAN(storecontactprop,(ITG *ne,ITG *ne0,char *lakon,ITG *kon,
         ITG *ipkon,ITG *mi,ITG *ielmat,double *elcon,ITG *mortar,
         double *adb,ITG *nactdof,double *springarea,ITG *ncmat_,
         ITG *ntmat_,double *stx,double *temax));

ITG strcmp1(const char *s1,const char *s2);

ITG strcmp2(const char *s1,const char *s2,ITG length);

ITG strcpy1(char *s1,const char *s2,ITG length);

void stress_sen(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
       double *stn,double *elcon,ITG *nelcon,
       double *rhcon,ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,
       ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
       double *t0,
       double *t1,ITG *ithermal,double *prestr,ITG *iprestr,char *filab,
       double *emn,
       double *een,ITG *iperturb,double *f,ITG *nactdof,
       double *vold,ITG *nodeboun,ITG *ndirboun,
       double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
       char *labmpc,ITG *nmpc,ITG *nmethod,double *cam,ITG *neq,double *veold,
       double *accold,double *bet,double *gam,double *dtime,double *time,
       double *ttime,double *plicon,ITG *nplicon,double *plkcon,
       ITG *nplkcon,double *xstateini,double *xstate,ITG *npmat_,
       double *epn,char *matname,ITG *mi,ITG *ielas,ITG *ncmat_,
       ITG *nstate_,
       double *stiini,double *vini,ITG *ikboun,ITG *ilboun,
       double *enern,double *emeini,double *xstaten,double *enerini,
       double *cocon,ITG *ncocon,char *set,ITG *nset,ITG *istartset,
       ITG *iendset,
       ITG *ialset,ITG *nprint,char *prlab,char *prset,double *qfx,double *qfn,
       double *trab,
       ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,ITG *nload,
       ITG *ikmpc,ITG *ilmpc,
       ITG *istep,ITG *iinc,double *springarea,double *reltime,ITG *ne0,
       double *xforc,ITG *nforc,double *thicke,
       double *shcon,ITG *nshcon,char *sideload,double *xload,
       double *xloadold,ITG *icfd,ITG *inomat,double *pslavsurf,
       double *pmastsurf,ITG *mortar,ITG *islavact,double *cdn,
       ITG *islavnode,ITG *nslavnode,ITG *ntie,double *clearini,
       ITG *islavsurf,ITG *ielprop,double *prop,double *energyini,
       double *energy,ITG *kscale,char *orname,ITG *network,
       ITG *nestart,ITG *neend,ITG *jqs,ITG *irows,ITG *nodedesi,
       double *xdesi,ITG *ndesi,ITG *iobject,ITG *nobject,char *objectset,
       double *g0,double *dgdx,ITG *idesvara,ITG *idesvarb,ITG *nasym,
       ITG *isolver,double *distmin,ITG *nodeset,double *b);

void *stress_senmt(ITG *i);

void stress_sen_2ndorder(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
       ITG *ne,
       double *stn,double *elcon,ITG *nelcon,
       double *rhcon,ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,
       ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
       double *t0,
       double *t1,ITG *ithermal,double *prestr,ITG *iprestr,char *filab,
       double *emn,
       double *een,ITG *iperturb,double *f,ITG *nactdof,
       double *vold,ITG *nodeboun,ITG *ndirboun,
       double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
       char *labmpc,ITG *nmpc,ITG *nmethod,double *cam,ITG *neq,double *veold,
       double *accold,double *bet,double *gam,double *dtime,double *time,
       double *ttime,double *plicon,ITG *nplicon,double *plkcon,
       ITG *nplkcon,double *xstateini,double *xstate,ITG *npmat_,
       double *epn,char *matname,ITG *mi,ITG *ielas,ITG *ncmat_,
       ITG *nstate_,
       double *stiini,double *vini,ITG *ikboun,ITG *ilboun,
       double *enern,double *emeini,double *xstaten,double *enerini,
       double *cocon,ITG *ncocon,char *set,ITG *nset,ITG *istartset,
       ITG *iendset,
       ITG *ialset,ITG *nprint,char *prlab,char *prset,double *qfx,double *qfn,
       double *trab,
       ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,ITG *nload,
       ITG *ikmpc,ITG *ilmpc,
       ITG *istep,ITG *iinc,double *springarea,double *reltime,ITG *ne0,
       double *xforc,ITG *nforc,double *thicke,
       double *shcon,ITG *nshcon,char *sideload,double *xload,
       double *xloadold,ITG *icfd,ITG *inomat,double *pslavsurf,
       double *pmastsurf,ITG *mortar,ITG *islavact,double *cdn,
       ITG *islavnode,ITG *nslavnode,ITG *ntie,double *clearini,
       ITG *islavsurf,ITG *ielprop,double *prop,double *energyini,
       double *energy,ITG *kscale,char *orname,ITG *network,
       ITG *nestart,ITG *neend,ITG *jqs,ITG *irows,ITG *nodedesi,
       double *xdesi,ITG *ndesi,ITG *iobject,ITG *nobject,char *objectset,
       double *g0,double *dgdx,ITG *idesvar1,ITG *idesvar2,ITG *nasym,
       ITG *isolver,double *distmin,ITG *nodeset,double *b1,double *b2,
       double *dgdxdy);

void *stress_sen_2ndordermt(ITG *i);

void stress_sen_dv(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
       double *dstn,double *elcon,ITG *nelcon,
       double *rhcon,ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,
       ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
       double *t0,double *t1,ITG *ithermal,double *prestr,ITG *iprestr,
       char *filab,ITG *iperturb,double *dv,ITG *nmethod,double *dtime,
       double *time,double *ttime,double *plicon,ITG *nplicon,double *plkcon, 
       ITG *nplkcon,double *xstateini,double *xstate,ITG *npmat_,char *matname,
       ITG *mi,ITG *ielas,ITG *ncmat_,ITG *nstate_,double *stiini,double *vini,  
       double *emeini,double *enerini,ITG *istep,ITG *iinc,double *springarea,
       double *reltime,ITG *ne0,double *thicke,double *pslavsurf,
       double *pmastsurf,ITG *mortar,double *clearini,ITG *ielprop,double *prop,   
       ITG *kscale,ITG *iobject,double *g0,ITG *nea,ITG *neb,ITG *nasym, 
       double *distmin,double *dstx,ITG *ialnk,double *dgdu,ITG *ialeneigh,   
       ITG *neaneigh,ITG *nebneigh,ITG *ialnneigh,ITG *naneigh,ITG *nbneigh,
       double *stn,double *expks,char *objectset,ITG *idof,ITG *node,
       ITG *idir,double *vold,double *dispmin);       

void *stress_sen_dvmt(ITG *i);

void stress_sen_dx(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
       double *dstn,double *elcon,ITG *nelcon,
       double *rhcon,ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,
       ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
       double *t0,double *t1,ITG *ithermal,double *prestr,ITG *iprestr,
       char *filab,ITG *iperturb,double *vold,ITG *nmethod,double *dtime,   
       double *time,double *ttime,double *plicon,ITG *nplicon,double *plkcon,  
       ITG *nplkcon,double *xstateini,double *xstate,ITG *npmat_,char *matname,
       ITG *mi,ITG *ielas,ITG *ncmat_,ITG *nstate_,double *stiini,double *vini,
       double *emeini,double *enerini,ITG *istep,ITG *iinc,double *springarea,
       double *reltime,ITG *ne0,double *thicke,double *pslavsurf,
       double *pmastsurf,ITG *mortar,double *clearini,ITG *ielprop,double *prop,
       ITG *kscale,ITG *iobject,char *objectset,double *g0,double *dgdx,
       ITG *nea,ITG *neb,ITG *nasym,double *distmin,ITG*idesvar,double *dstx,   
       ITG *ialdesi,ITG *ialeneigh,ITG *neaneigh,ITG *nebneigh,ITG *ialnneigh,
       ITG *naneigh,ITG *nbneigh,double *stn,double *expks,ITG *ndesi);  

void *stress_sen_dxmt(ITG *i);

void FORTRAN(subspace,(double *d,double *aa,double *bb,double *cc,
             double *alpham,double *betam,ITG *nev,
             double *xini,double *cd,double *cv,double *time,
             double *rwork,ITG *lrw,ITG *k,ITG *jout,double *rpar,
             double *bj,ITG *iwork,ITG *liw,ITG *iddebdf,double *bjp));

void FORTRAN(tempload,(double *xforcold,double *xforc,double *xforcact,
               ITG *iamforc,ITG *nforc,double *xloadold,double *xload,
               double *xloadact,ITG *iamload,ITG *nload,ITG *ibody,
               double *xbody,ITG *nbody,double *xbodyold,double *xbodyact,
               double *t1old,double *t1,double *t1act,ITG *iamt1,
               ITG *nk,double *amta,ITG *namta,ITG *nam,double *ampli,
               double *time,double *reltime,double *ttime,double *dtime,
               ITG *ithermal,ITG *nmethod,
               double *xbounold,double *xboun,double *xbounact,
               ITG *iamboun,ITG *nboun,ITG *nodeboun,
               ITG *ndirboun,ITG *nodeforc,ITG *ndirforc,ITG *istep,
               ITG *iint,double *co,double *vold,ITG *itg,ITG *ntg,
               char *amname,ITG *ikboun,ITG *ilboun,ITG *nelemload,
               char *sideload,ITG *mi,ITG *ntrans,double *trab,
               ITG *inotr,double *veold,ITG *integerglob,
               double *doubleglob,char *tieset,ITG *istartset,
               ITG *iendset,ITG *ialset,ITG *ntie,ITG *nmpc,ITG *ipompc,
               ITG *ikmpc,ITG *ilmpc,ITG *nodempc,double *coefmpc,
               ITG *ipobody,ITG *iponoel,ITG *inoel,ITG *ipkon,ITG *kon,
               ITG *ielprop,double *prop,ITG *ielmat,double *shcon,
               ITG *nshcon,double *rhcon,ITG *nrhcon,double *cocon,
               ITG *ncocon,ITG *ntmat_,char *lakon));

void FORTRAN(tempload_em,(double *xforcold,double *xforc,double *xforcact,
               ITG *iamforc,ITG *nforc,double *xloadold,double *xload,
               double *xloadact,ITG *iamload,ITG *nload,ITG *ibody,
               double *xbody,ITG *nbody,double *xbodyold,double *xbodyact,
               double *t1old,double *t1,double *t1act,ITG *iamt1,
               ITG *nk,double *amta,ITG *namta,ITG *nam,double *ampli,
               double *time,double *reltime,double *ttime,double *dtime,
               ITG *ithermal,ITG *nmethod,
               double *xbounold,double *xboun,double *xbounact,
               ITG *iamboun,ITG *nboun,ITG *nodeboun,
               ITG *ndirboun,ITG *nodeforc,ITG *ndirforc,ITG *istep,
               ITG *iint,double *co,double *vold,ITG *itg,ITG *ntg,
               char *amname,ITG *ikboun,ITG *ilboun,ITG *nelemload,
               char *sideload,ITG *mi,ITG *ntrans,double *trab,
               ITG *inotr,double *veold,ITG *integerglob,
               double *doubleglob,char *tieset,ITG *istartset,
               ITG *iendset,ITG *ialset,ITG *ntie,ITG *nmpc,ITG *ipompc,
               ITG *ikmpc,ITG *ilmpc,ITG *nodempc,double *coefmpc,
               double *h0scale,ITG *inomat,ITG *ipobody,ITG *iponoel,
               ITG *inoel,ITG *ipkon,ITG *kon,char *lakon,ITG *ielprop,
               double *prop,ITG *ielmat,double *shcon,ITG *nshcon,
               double *rhcon,ITG *nrhcon,ITG *ntmat_,double *cocon,
               ITG *ncocon));

void FORTRAN(temploaddiff,(double *xforcold,double *xforc,double *xforcact,
               ITG *iamforc,ITG *nforc,double *xloadold,double *xload,
               double *xloadact,ITG *iamload,ITG *nload,ITG *ibody,
               double *xbody,ITG *nbody,double *xbodyold,double *xbodyact,
               double *t1old,double *t1,double *t1act,ITG *iamt1,
               ITG *nk,double *amta,ITG *namta,ITG *nam,double *ampli,
               double *time,double *reltime,double *ttime,double *dtime,
               ITG *ithermal,ITG *nmethod,
               double *xbounold,double *xboun,double *xbounact,
               ITG *iamboun,ITG *nboun,ITG *nodeboun,
               ITG *ndirboun,ITG *nodeforc,ITG *ndirforc,ITG *istep,
               ITG *iint,double *co,double *vold,ITG *itg,ITG *ntg,
               char *amname,ITG *ikboun,ITG *ilboun,ITG *nelemload,
               char *sideload,ITG *mi,double *xforcdiff,double *xloaddiff,
               double *xbodydiff,double *t1diff,double *xboundiff,
               ITG *icorrect,ITG *iprescribedboundary,ITG *ntrans,
               double *trab,ITG *inotr,double *veold,ITG *nactdof,
               double *bcont,double *fn,ITG *ipobody,ITG *iponoel,
               ITG *inoel,ITG *ipkon,ITG *kon,char *lakon,ITG *ielprop,
               double *prop,ITG *ielmat,double *shcon,ITG *nshcon,
               double *rhcon,ITG *nrhcon,ITG *ntmat_,double *cocon,
               ITG *ncocon));

void FORTRAN(temploadfem,(double *xforcold,double *xforc,double *xforcact,
               int *iamforc,int *nforc,double *xloadold,double *xload,
               double *xloadact,int *iamload,int *nload,int *ibody,
               double *xbody,int *nbody,double *xbodyold,double *xbodyact,
               double *t1old,double *t1,double *t1act,int *iamt1,
               int *nk,double *amta,int *namta,int *nam,double *ampli,
               double *time,double *reltime,double *ttime,double *dtime,
               int *ithermal,int *nmethod,
               double *xbounold,double *xboun,double *xbounact,
               int *iamboun,int *nboun,int *nodeboun,
               int *ndirboun,int *nodeforc,int *ndirforc,int *istep,
               int *iint,double *co,double *vold,int *itg,int *ntg,
               char *amname,int *ikboun,int *ilboun,int *nelemload,
               char *sideload,int *mi,int *ntrans,double *trab,
               int *inotr,double *veold,int *integerglob,
               double *doubleglob,char *tieset,int *istartset,
               int *iendset,int *ialset,int *ntie,int *nmpc,int *ipompc,
               int *ikmpc,int *ilmpc,int *nodempc,double *coefmpc));

void FORTRAN(temploadmodal,(double *amta,ITG *namta,ITG *nam,double *ampli,
         double *timemin,double *ttimemin,double *dtime,double *xbounold,
         double *xboun,double *xbounmin,ITG *iamboun,ITG *nboun,
         ITG *nodeboun,ITG *ndirboun,char *amname,double *reltime));

void thicknessmain(double *co,double *dgdx,ITG *nobject,ITG *nk,
                ITG *nodedesi,ITG *ndesi,char *objectset,ITG *ipkon,
                ITG *kon,char *lakon,char *set,ITG *nset,ITG *istartset,
                ITG *iendset,ITG *ialset,ITG *iobject,ITG *nodedesiinv,
                double *dgdxglob);

void *thicknessmt(ITG *i);
    
void FORTRAN(tiefaccont,(char *lakon,ITG *ipkon,ITG *kon,ITG *ntie,
       char *tieset,ITG *nset,char *set,ITG *istartset,ITG *iendset,
       ITG *ialset,ITG *itiefac,ITG *islavsurf,ITG *islavnode,
       ITG *imastnode,ITG *nslavnode,ITG *nmastnode,ITG *nslavs,
       ITG *nmasts,ITG *ifacecount,ITG *iponoels,ITG *inoels,ITG *ifreenoels,
       ITG *mortar,ITG *ipoface,ITG *nodface,ITG *nk,double *xnoels));   

void tiedcontact(ITG *ntie,char *tieset,ITG *nset,char *set,
               ITG *istartset,ITG *iendset,ITG *ialset,
               char *lakon,ITG *ipkon,ITG *kon,double *tietol,
               ITG *nmpc,ITG *mpcfree,ITG *memmpc_,
               ITG **ipompcp,char **labmpcp,ITG **ikmpcp,ITG **ilmpcp,
               double **fmpcp,ITG **nodempcp,double **coefmpcp,
               ITG *ithermal,double *co,double *vold,ITG *nef,
               ITG *nmpc_,ITG *mi,ITG *nk,ITG *istep,ITG *ikboun,
               ITG *nboun,char *kind1,char *kind2);
        
void FORTRAN(transformatrix,(double *xab,double *p,double *a));

void FORTRAN(transition,(double *dgdxglob,ITG *nobject,ITG *nk,ITG *nodedesi,
                ITG *ndesi,char *objectset,double *xo,double *yo,double *zo,
                double *x,double *y,double *z,ITG *nx,ITG *ny,ITG *nz,
                double *co,ITG *ifree,ITG *ndesia,ITG *ndesib));                         

void transitionmain(double *co,double *dgdxglob,ITG *nobject,ITG *nk,
                ITG *nodedesi,ITG *ndesi,char *objectset,ITG *ipkon,
                ITG *kon,char *lakon,ITG *ipoface,ITG *nodface,
                ITG *nodedesiinv);

void *transitionmt(ITG *i);

void FORTRAN(trianeighbor,(ITG *ipe,ITG *ime,ITG *imastop,ITG *ncont,
               ITG *koncont,ITG *ifreeme));

void FORTRAN(triangucont,(ITG *ncont,ITG *ntie,char *tieset,ITG *nset,
          char *set,ITG *istartset,ITG *iendset,ITG *ialset,ITG *itietri,
          char *lakon,ITG *ipkon,ITG *kon,ITG *koncont,char *kind1,
          char *kind2,double *co,ITG *nk));

void FORTRAN(tridiagonal_nrhs,(double *a,double *b,ITG *n,ITG *m,
             ITG *nrhs));

void FORTRAN(turningdirection,(double *v,double *e1,double *e2,double *xn,
                               ITG *mi,
                               ITG *nk,char *turdir,char *lakon,ITG *ipkon,
                               ITG *kon,ITG *ne,double *co));

#ifdef BAM
void FORTRAN(uexternaldb,(ITG *lop,ITG *lrestart,double *time,double *dtime,
                          ITG *kstep,ITG *kinc));
#endif

void FORTRAN(ufaceload,(double *co,ITG *ipkon,ITG *kon,char *lakon,
                        ITG *nboun,ITG *nodeboun,
                        ITG *nelemload,char *sideload,ITG *nload,
                        ITG *ne,ITG *nk));

void FORTRAN(uinit,());

void FORTRAN(uiter,(ITG *iit));

void FORTRAN(uout,(double *v,ITG *mi,ITG *ithermal,char *filab,
                 ITG *kode));

void FORTRAN(updatecon,(double *vold,double *voldaux,double *v,int *nk,
           int *ielmat,int *ntmat_,double *shcon,int *nshcon,double *rhcon,
           int *nrhcon,int *iout,int *nmethod,int *convergence,
           double *physcon,int *iponoel,int *inoel,int *ithermal,
           int *nactdoh,int *iit,int *compressible,int *ismooth,
           double *vcontu,double *vtu,int *turbulent,int *inomat,
           int *nodeboun,int *ndirboun,int *nboun,int *mi,double *co,
           double *factor));

void FORTRAN(updatecont,(ITG *koncont,ITG *ncont,double *co,double *vold,
                         double *cg,double *straight,ITG *mi));

void FORTRAN(updatecontpen,(ITG *koncont,ITG *ncont,double *co,double *vold,
                         double *cg,double *straight,ITG *mi,ITG *imastnode,
                         ITG *nmastnode,double *xmastnor,ITG *ntie,
                         char *tieset,ITG *nset,char *set,ITG *istartset,
                         ITG *iendset,ITG *ialset,ITG *ipkon,char *lakon,
                         ITG *kon,double *cs,ITG *mcs,ITG *ics));

void *u_calloc(size_t num,size_t size,const char *file,const int line,const char *ptr_name);

void *u_free(void* num,const char *file,const int line,const char *ptr_name);

void *u_malloc(size_t size,const char *file,const int line,const char *ptr_name);

void *u_realloc(void* num,size_t size,const char *file,const int line,const char *ptr_name);

void utempread(double *t1,ITG *istep,char *jobnamec);

void FORTRAN(velinireltoabs,(ITG *ibody,double *xbody,char *cbody,ITG *nbody,
                             char *set,ITG *istartset,
                             ITG *iendset,ITG *ialset,ITG *nset,
                             double *veold,ITG *mi,ITG *ipkon,ITG *kon,
                             char *lakon,double *co,ITG *itreated));

void FORTRAN(velsolve,(ITG *nef,ITG *ipnei,double *bv,double *auv,double *adv,
                       double *vel,double *temp,ITG *neiel));

void worparll(double *allwk,double *fnext,ITG *mt,double *fnextini,
                   double *v,double *vini,ITG *nk,ITG *num_cpus);

void *worparllmt(ITG *i);

void writeBasisParameter(FILE *f,ITG *istep,ITG *iinc);

void FORTRAN(writeboun,(ITG *nodeboun,ITG *ndirboun,double *xboun,
      char *typeboun,ITG *nboun));

void FORTRAN(writebv,(double *,ITG *));

void FORTRAN(writecvg,(ITG *itep,ITG *iinc,ITG *icutb,ITG *iit,ITG *ne,ITG *ne0,
                       double *ram,double *qam,double *cam,double *uam,
                       ITG *ithermal));

void FORTRAN(writedeigdx,(ITG *iev,double *d,ITG *ndesi,char*orname,
                          double *dgdx));

void FORTRAN(writedesi,(ITG *norien,char *orname));

void FORTRAN(writeelem,(ITG *i,char *lakon));

void FORTRAN(writeev,(double *,ITG *,double *,double *));

void FORTRAN(writeevcomplex,(double *eigxx,ITG *nev,double *fmin,double *fmax));

void FORTRAN(writeevcs,(double *,ITG *,ITG *,double *,double *));

void FORTRAN(writeevcscomplex,(double *eigxx,ITG *nev,ITG *nm,double *fmin,
            double *fmax));

void FORTRAN(writehe,(ITG *));

void writeheading(char *jobnamec,char *heading,ITG *nheading);

void FORTRAN(writeim,());

void FORTRAN(writeinput,(char *inpc,ITG *ipoinp,ITG *inp,ITG *nline,ITG *ninp,
                         ITG *ipoinpc));

void FORTRAN(writelm,(ITG *iter,double *lambda,ITG *nactive,ITG *nnlconst,
                      char *objectset,ITG *nobject,ITG *ipacti,
                      ITG *iconstacti,ITG *inameacti));

void FORTRAN(writemac,(double *mac,ITG *nev,ITG *nevcomplex));

void FORTRAN(writemaccs,(double *mac,ITG *nev,ITG* nm));

void FORTRAN(writemeshinp,(ITG *kontet,ITG *netet_,double *cotet,ITG *nktet,
                           ITG *ij));

void FORTRAN(writempc,(ITG *,ITG *,double *,char *,ITG *));

void FORTRAN(writeobj,(char *objectset,ITG *iobject,double *g0));

void FORTRAN(writepf,(double *d,double *bjr,double *bji,double *freq ,
                      ITG *nev,ITG *mode,ITG *nherm));

void FORTRAN(writerandomfield,(double *d,ITG *nev,double *abserr,
             double *relerr));

void FORTRAN(writere,());

void FORTRAN(writerefinemesh,(ITG *kontet,ITG *netet_,double *cotet,ITG *nktet,
                              char *jobnamec,ITG *ipkon,ITG *kon,char *lakon,
                              ITG *iquad,ITG *iedtet,ITG *iedgmid,
                              ITG *itreated,ITG *ne));

void FORTRAN(writesen,(double *g0,double *dgdx,ITG *ndesi,ITG *nobject,
                       ITG *nodedesi,char *jobnamef));

void FORTRAN(writesta,(ITG *istep,ITG *j,ITG *icutb,ITG *l,double *ttime,
                       double *time,double *dtime));

void FORTRAN(writestadiv,(ITG *istep,ITG *j,ITG *icutb,ITG *l,double *ttime,
                       double *time,double *dtime));

void FORTRAN(writesubmatrix,(double *submatrix,ITG *noderetain,
             ITG *ndirretain,ITG *nretain,char *jobnamec));

void FORTRAN(writesummary,(ITG *istep,ITG *j,ITG *icutb,ITG *l,double *ttime,
                   double *time,double *dtime));

void FORTRAN(writesummarydiv,(ITG *istep,ITG *j,ITG *icutb,ITG *l,double *ttime,
                   double *time,double *dtime));

void FORTRAN(writetetmesh,(ITG *kontet,ITG *netet_,double *cotet,
     ITG *nktet,double *field,ITG *nfield));
                        
void FORTRAN(writeturdir,(double *xn,char *turdir,ITG *nev));
                        
void FORTRAN(writeturdircs,(double *xn,char *turdir,ITG *nev,ITG *nm));
                        
void FORTRAN(writeview,(ITG *ntr,double *adview,double *auview,double *fenv,
            ITG *nzsrad,char *jobnamef));

void FORTRAN(zienzhu,(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
                      ITG *ne,double *stn,ITG *ipneigh,ITG *neigh,
                      double *sti,ITG *mi));

void FORTRAN(znaupd,(ITG *ido,char *bmat,ITG *n,char *which,ITG *nev,
             double *tol,double *resid,ITG *ncv,double *z,ITG *ldz,
             ITG *iparam,ITG *ipntr,double *workd,double *workl,
             ITG *lworkl,double *rwork,ITG *info));

void FORTRAN(zneupd,(ITG *rvec,char *howmny,ITG *select,double *d,
             double *z,ITG *ldz,double *sigma,
             double *workev,char *bmat,ITG *neq,char *which,
             ITG *nev,double *tol,double *resid,ITG *ncv,double *v,
             ITG *ldv,ITG *iparam,ITG *ipntr,double *workd,
             double *workl,ITG *lworkl,double *rwork,ITG *info));

#ifdef CALCULIX_EXTERNAL_BEHAVIOURS_SUPPORT

/*!
 * list of supported external behaviours
 */
typedef enum {
  //! interface to standard calculix behaviours
  CALCULIX_STANDARD_INTERFACE,
  //! interface to abaqus behaviours in small strain
  CALCULIX_ABAQUS_INTERFACE,
  //! interface to abaqus behaviours in finite strain
  CALCULIX_ABAQUSNL_INTERFACE
} CalculixInterface;
/*!
 * a structure describing an external behaviour
 */
typedef struct
{
  // name of the material
  const char *n;
  // interface
  CalculixInterface i;
  // function pointer
  const void* ptr;
}  CalculixExternalBehaviour;
/*!
 * \return the description of an external beahviour
 * \param[in] n : external behaviour name
 */
const CalculixExternalBehaviour*
calculix_searchExternalBehaviour(const char*);
/*!
 * \param[in] n: material name
 */
void calculix_registerExternalBehaviour(const char *);
/*!
 * \brief free the memory associated with external behaviours
 * treatment.
 */
void calculix_freeExternalBehaviours();

#endif /* CALCULIX_EXTERNAL_BEHAVIOURS_SUPPORT */
