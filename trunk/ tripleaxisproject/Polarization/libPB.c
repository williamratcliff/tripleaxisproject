/*
  C library for correcting polarized beam data using He-3 CELLS.
  Main entry is PBcorrectData.
*/

#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <math.h>

#include "PB.h"

/*
  local declarations do not need to be in PB.h
  PB.h only contains stuff required by user of libPB
  i.e. decl of PBcorrectData  PBsim PBreadflags and their args
 */

/*
  Set up a data structure to handle up to 4 linear equations including
  error propagation for pol beam corrections, to solve for S(cross-sections).
  Y(countrates) = Cmatrix * S(cross-sections)
  Will need flags to show which equations are active, since linear constraints
  can be added. Also need to be able to copy the equ structure.
  Equations can be combined as well.
  Original 4 equations indexed 0-3 are for C++ C-- C+- C--   +- is front flip
  When equations are added these indices are no longer meaningful.
  Keep indices of active equations and free variables.
*/

typedef struct {
  unsigned long sec[4] ;
  double Y[4], Yesq[4], lambI, lambF ;
  double C[4][4], Cesq[4][4] ;
  double S[4], Sesq[4];
  int Nactive, activeEq[4];
  int Nfree, freeS[4];
} PBdatapt;

typedef struct {
  int freeToConstrain; /* index 0-3 of free variable to constrain */
  int Nfree, freeS[4];
  double Coef[4];
} ConstraintEq;

/*
  constraint eq
  S[freeToConstrain] = SUM(j to Nfree) Coef[j] * S[freeVar[j]]
*/  

typedef struct {
  char   name[64];
  double tEmpty; /* at 1.77 A, may eventually have to put in wavelength dep */
  double tEmptySlope; /* wavelength dependence */
  double L;      /* max gas length in cm */
  double D;      /* diameter in cm for id purposes only */
  double R;      /* radius of curvature of windows, may want 2 values */
  double nsL0;   /* uncorrected nsL at 1.77 A */
  double nsL0err;
  double nsL;    /* corrected for curvature */
  double nsLerr;
} He3CELL;       /* base He3 CELL constants */

typedef struct {
  double PHe;    /* init He3 polarization */
  double PHeErr; 
  double hrsBeam; /* beamtime hours logged this exp */
  double T;      /* cell decay time in hours in holding fld */
  double Terr;
  unsigned long startSecs; /* conv of startDate to seconds since Jan1 1971 */
} He3CELLpol;

typedef struct {
  double waveRelWidth;   /* std-dev-Lambda/Lambda */
  double angleVwidth;    /* vertical angular resolution stnd dev */
  double angleHwidth;    /* horizontal angular resolution stnd dev */
  double usedRadius;     /* eff beam radius in cm for cell */
} expResol;

typedef struct {
  He3CELL cell;
  He3CELLpol pol;
  expResol res;
} He3CELLexp;


typedef struct {
  double teff;  /* transport efficiency */
  double terr;
  double feff;  /* flipper efficiency */
  double ferr;
} efficiency;

typedef struct {
  He3CELLexp P;
  He3CELLexp A;
  efficiency eP;
  efficiency eA;
} PBsetup;

double Dn = 2.072141789 ;
double TWOPI = 6.283185307 ;

static double correctionCoef(He3CELLexp *ex, double tau) ;
static int transfac(He3CELLexp *ex, unsigned long secs, double lambda,
	     double *tp, double *tm, double *tpesq, double *tmesq) ;
static int PBcoef(PBdatapt *d, PBsetup *s) ;
static int PBcorrect(PBdatapt *d) ;
static int PBcorrectDatapt(PBdatapt *d) ;
static int applyConstraint(PBdatapt *d, ConstraintEq *eq) ;
static int combineEqs(PBdatapt *d, int eq1, int eq2) ;
static int deleteEq(PBdatapt *d, int eq) ;
static int constrainResult(PBdatapt *d) ;
static void constraintTOeqs() ;
static double poisson(int idum, double d1, double d2, double *distval,
		      double min, double max) ;
static double uniform(int idum, double d1, double d2, double *distval,
		      double min, double max) ;
static double gaussian(int idum, double ave, double sig, double *distval,
		       double d1, double d2) ;

static int monitorCorrectPG(PBdatapt *d) ;
static int polmonitorCorrectPG(PBdatapt *d, He3CELLexp *pol) ;
static int SIMmonitorCorrectPG(PBdatapt *d) ;
static int SIMpolmonitorCorrectPG(PBdatapt *d, He3CELLexp *pol) ;


static int PBdefinePolarizer(char *filename) ;
static int PBdefineAnalyzer(char *filename) ;
static int PBsetflags(PBflags *flgs) ;

/* create one instance of a PB setup to store He3CELL info etc */
static PBsetup s ;
/* and a default flags structure */
static PBflags flags = {
  0, 0, 0, 100000, 0,
  {1, 1, 1, 1},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0}
} ;
static PBflags defltflags = {
  0, 0, 0, 100000, 0,
  {1, 1, 1, 1},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0}
} ;

static ConstraintEq eqs[4] ;
static int Nconstraint = 0 ;

static int *DBG = &(flags.Debug) ;






/*
  correct PB data passed in structure
  Any of the passed structure pointers can be NULL
  except that if PBindata is not NULL then
  PBoutdata must be also non-NULL or error return (!=0).
  It is of course also an error if any of the double pointers
  in the passed PBindata or PBoutdata are NULL.
  Caller is responsible for storage allocation in
  PBindata and PBoutdata which must have double arrays
  that can hold up to npts.
*/

int PBcorrectData(char *PCellFile, char *ACellFile, PBflags *flgs,
		  int npts, PBindata *in, PBoutdata *out)
{
  /* storage for PB data */
  static PBdatapt d ;

  int i, ierr, j ;
  double temp ;
  FILE *fp ;
  printf("reached\n");
  printf("%s\n",PCellFile);
  printf("%s\n",ACellFile);
  printf("npts %d\n",npts);
  printf("in\n");
  //for (i=0; i <npts;i++){
  //    printf("%d %f %f %f %f %u %f %f %u\n",i,in->Ei[i],in->Ef[i],in->Cpm[i],in->Epm[i],in->tpm[i],in->Cmp[i],in->Emp[i],in->tmp[i]);
  //    }
  printf("MonitorCorrect=%d\n",flgs->MonitorCorrect);
  printf("PolMonitorCorrect=%d\n",flgs->PolMonitorCorrect);
  printf("Debug=%d\n",flgs->Debug);
  printf("SimFlux=%d\n",flgs->SimFlux);
  printf("SimDeviate=%d\n",flgs->SimDeviate);
  printf("CountsEnable=%d %d %d %d\n",flgs->CountsEnable[0],flgs->CountsEnable[1],flgs->CountsEnable[2],flgs->CountsEnable[3]);
  printf("SConstrain=%d %d %d %d\n",flgs->Sconstrain[0],flgs->Sconstrain[1],flgs->Sconstrain[2],flgs->Sconstrain[3]);
  printf("Smp=%f %f %f %f\n",flgs->Smp[0],flgs->Smp[1],flgs->Smp[2],flgs->Smp[3]);
  printf("Spm=%f %f %f %f\n",flgs->Spm[0],flgs->Spm[1],flgs->Spm[2],flgs->Spm[3]);
  printf("Spp=%f %f %f %f\n",flgs->Spp[0],flgs->Spp[1],flgs->Spp[2],flgs->Spp[3]);
  printf("Smm=%f %f %f %f\n",flgs->Smm[0],flgs->Smm[1],flgs->Smm[2],flgs->Smm[3]);
  
  if( PCellFile != NULL ) if( PBdefinePolarizer(PCellFile) ) return 1 ;
  if( ACellFile != NULL ) if( PBdefineAnalyzer(ACellFile) ) return 2 ;
  
  printf("Cells defined\n");
  if( flgs != NULL ) PBsetflags(flgs) ;
  if( in == NULL || npts < 1 ) return 0 ;
  if( out == NULL ) return 3 ;

  constraintTOeqs() ;

  /*
    process each data point
  */
  
  for( i=0 ; i<npts ; i++ ) {
    /* first just copy the input to the PBdatapt struct d */
    d.Y[0] = in->Cpp[i] ;
    d.Y[1] = in->Cmm[i] ;
    d.Y[2] = in->Cpm[i] ;
    d.Y[3] = in->Cmp[i] ;
    d.Yesq[0] = in->Epp[i]*in->Epp[i] ;
    d.Yesq[1] = in->Emm[i]*in->Emm[i] ;
    d.Yesq[2] = in->Epm[i]*in->Epm[i] ;
    d.Yesq[3] = in->Emp[i]*in->Emp[i] ;

    d.sec[0] = in->tpp[i] ;
    d.sec[1] = in->tmm[i] ;
    d.sec[2] = in->tpm[i] ;
    d.sec[3] = in->tmp[i] ;
    temp = in->Ei[i] ;
    if( temp <= 0. ) return 4 ;
    temp /= Dn ;
    d.lambI = TWOPI/sqrt(temp) ;
    temp = in->Ef[i] ;
    if( temp <= 0. ) return 4 ;
    temp /= Dn ;
    d.lambF = TWOPI/sqrt(temp) ;

    if( (ierr = PBcorrectDatapt(&d)) ) {
      printf("ERROR in PBcorrectDatapt = %d\n",ierr) ;
    }
    /* transfer the corrected S to the PBoutdata struct */
    for( j=0 ; j<4 ; j++ ) {
      if( j == 0 ) {
	out->Spp[i] = d.S[0] ;
	if(d.Sesq[0]>0.) out->Epp[i] = sqrt(d.Sesq[0]) ; else out->Epp[i]=0. ;
      } else if( j == 1 ) {
	out->Smm[i] = d.S[1] ;
	if(d.Sesq[1]>0.) out->Emm[i] = sqrt(d.Sesq[1]) ; else out->Emm[i]=0. ;
      } else if( j == 2 ) {
	out->Spm[i] = d.S[2] ;
	if(d.Sesq[2]>0.) out->Epm[i] = sqrt(d.Sesq[2]) ; else out->Epm[i]=0. ;
      } else if( j == 3 ) {
	out->Smp[i] = d.S[3] ;
	if(d.Sesq[3]>0.) out->Emp[i] = sqrt(d.Sesq[3]) ; else out->Emp[i]=0. ;
      }
    }

  }
  return 0 ;
}

static int PBcorrectDatapt(PBdatapt *d)
{
  int j, k, ierr ;
  if( flags.MonitorCorrect ) {
    monitorCorrectPG(d) ;
  } else if ( flags.PolMonitorCorrect ) {
    polmonitorCorrectPG(d, &(s.P));
  }
  
  /* intial all 4 equations active and all S(cross-sections) free unknowns */
  d->Nactive = 4 ;
  d->Nfree = 4 ;
  for( j=0 ; j<4 ; j++ ) {
    d->activeEq[j] = j ;
    d->freeS[j] = j ;
  }
  
  /* calc the correction coeficients */
  if( (ierr = PBcoef(d, &s)) > 0 ) {
    printf("error in PBcoef = %d\n", ierr) ;
    return 6 ;
  }
  
  /* apply any constraints on the S */
  for( j=0 ; j<Nconstraint ; j++ ) {
    if( applyConstraint(d, eqs+j) > 0 )
      printf("constraint %d failed\n",j) ;
  }

  /*
    add equations if requested, only allowed to add 2 equations twice
    The second equation becomes inactive with new effective coefs etc
    occupying the first equation's position.
  */
  for( j=0 ; j<4 ; j++ ) {
    if( flags.CountsAdd1[j] ) {
      for( k=j+2 ; k<4 ; k++ ) {
	if( flags.CountsAdd1[k] ) {
	  if( combineEqs(d, j, k) ) printf("combineEqs failed\n") ;
	}
      }
    }
    if( flags.CountsAdd2[j] ) {
      for( k=j+2 ; k<4 ; k++ ) {
	if( flags.CountsAdd2[k] ) {
	  if( combineEqs(d, j, k) ) printf("combineEqs failed\n") ;
	}
      }
    }
  }

  /* delete any equations where counts are disabled */
  for( j=0 ; j<4 ; j++ ) { if( ! flags.CountsEnable[j] ) deleteEq(d, j) ; }
    
  /* check that Nactive equations == Nfree Unknowns */
  if( d->Nactive != d->Nfree ) {
    printf("Nactive != Nfree\n") ;
    return 7 ;
  }
  
  if( *DBG ) {
    printf("after CTRLS activeEqindices and freeSindices are\n") ;
    for( j=0 ; j<d->Nactive ; j++ ) printf("%1d ",d->activeEq[j]) ;
    printf("\n") ;
    for( j=0 ; j<d->Nfree ; j++ ) printf("%1d ",d->freeS[j]) ;
    printf("\n\n") ;
  }
  
  /* call PBcorrect to solve for the cross-sections by Gaussian elim */
  if( (ierr = PBcorrect(d)) > 0 ) {
    if( ierr == 2 ) {
      printf("invalid Nactive equations or Nactive != Nfree\n") ;
      return 8 ;
    }
    printf("error in PBcorrect = %d\n", ierr) ;
    return 9 ;
  }
  
  /* apply any constraints to the final result */
  constrainResult(d) ;
  return 0 ;
}

static void constraintTOeqs()
{
  /*
    process the S(cross-section) contraints in global flags
  */
  int i, j, k, nc ;
  nc = 0 ;
  for( i=0 ; i<4 ; i++ ) {
    if( flags.Sconstrain[i] ) {
      eqs[nc].freeToConstrain = i ;
      k = 0 ;
      for( j=0 ; j<4 ; j++ ) {
	if( j == i ) continue ;
	if( i == 0 && fabs(flags.Spp[j]) > 0. ) {
	  eqs[nc].freeS[k] = j ;
	  eqs[nc].Coef[k] = flags.Spp[j] ;
	  k++ ;
	} else if( i == 1 && fabs(flags.Smm[j]) > 0. ) {
	  eqs[nc].freeS[k] = j ;
	  eqs[nc].Coef[k] = flags.Smm[j] ;
	  k++ ;
	} else if( i == 2 && fabs(flags.Spm[j]) > 0. ) {
	  eqs[nc].freeS[k] = j ;
	  eqs[nc].Coef[k] = flags.Spm[j] ;
	  k++ ;
	} else if( i == 3 && fabs(flags.Smp[j]) > 0. ) {
	  eqs[nc].freeS[k] = j ;
	  eqs[nc].Coef[k] = flags.Smp[j] ;
	  k++ ;
	}
      }
      eqs[nc].Nfree = k ;
      nc++ ;
    }
  }
  Nconstraint = nc ;
}

/*
  I want to change transfac so that it returns the 1+PHe3 and 1-PHe3
  cell transmissions since it looks like this routine will be needed
  to calculate transmissions for other wavelengths e.g. when counting
  against a monitor after the polarizer.
  Why not go ahead and calc the transmissions for lambda/2 and lambda/3?
  The overhead is not very high.
  But then the function would return 6 transmissions along with
  6 error values.
  In any case it is probably a good idea to just return tp and tm
  and their separate errors instead of ts and ta.
*/

static int transfac(He3CELLexp *ex, unsigned long secs, double lambda,
		     double *tp, double *tm, double *tpesq, double *tmesq)
{
  /*
    given cell pointer, data-time-stamp and the wavelength for that cell
    compute the symmetric and antisymmetric transmission factors and the
    error sum.
  */
  He3CELL *cell ;
  He3CELLpol *pol ;
  expResol *res ;


  double lamfac ;
  double nsL, nsLerr ;
  double tEmpty ;

  unsigned long tdiff ;
  double tr, pre, tre;
  double PHe3, PHe3esq, tau, tauesq;
  double C, etermA, etermP ;


  if( ex == NULL | tp == NULL || tm == NULL || tpesq == NULL || tmesq == NULL )
    return 1 ;

  cell = &(ex->cell) ;
  pol = &(ex->pol) ;
  res = &(ex->res) ;

  if( fabs(pol->T) <= 0. ) return 2 ;

  lamfac = lambda/1.77 ;
  nsL = cell->nsL * lamfac ;
  nsLerr = cell->nsLerr * lamfac ;
  tEmpty = cell->tEmpty + lamfac*cell->tEmptySlope ;

  /*
    Here error in He3 polarization increases with time.
    Another procedure should be used if flipping ratios are
    used to help determine the dacay.
   */
  tdiff = secs - pol->startSecs ;
  tr = (tdiff/3600.)/pol->T ;
  PHe3 = exp(-tr) ;
  pre = pol->PHeErr ;
  tre = pol->Terr/pol->T ;
  PHe3esq = PHe3*PHe3*pre*pre ;
  PHe3 *= pol->PHe ;
  PHe3esq += PHe3*PHe3*tr*tr*tre*tre ;
  tauesq = nsL*nsL*PHe3esq + nsLerr*nsLerr ;

  tau = nsL*(1-PHe3) ;
  C = correctionCoef(ex, tau) ;
  *tp = C*tEmpty*exp(-tau) ;
  *tpesq = (*tp)*(*tp)*tauesq ;

  tau = nsL*(1+PHe3) ;
  C = correctionCoef(ex, tau) ;
  *tm = C*tEmpty*exp(-tau) ;
  *tmesq = (*tm)*(*tm)*tauesq ;

  /* *ts = tp + tm ; */
  /* symmetric transmission coef */
  /* *ta = tp - tm ; */
  /* antisymmetric transmission coef */
  /* *tesq = tpesq + tmesq ; */
  return 0 ;
}

static int PBcoef(PBdatapt *d, PBsetup *s)
{
  int i, j, ierr ;

  /* index order is ++ -- +- -+ */
  static double alpmu[4] = { 1., -1.,  1., -1. } ;
  static double betnu[4] = { 1., -1., -1.,  1. } ;

  double alp1sq, bet1sq ;
  double eAalp, eAFalp, ePbet, ePFbet ;
  double eAalpsq, eAFalpsq, ePbetsq, ePFbetsq ;

  double cA, cAsq, cP, cPsq ;
  double etermA, etermP ;

  double etAsq, etAesq, efA, efAesq ;
  double etPsq, etPesq, efP, efPesq ;

  double alp, bet, mu, nu ;

  double tpA, tmA, tpAesq, tmAesq, tpP, tmP, tpPesq, tmPesq ;
  double tsA, taA, tAesq, tsP, taP, tPesq ;

  if( d == NULL || s == NULL ) return 1;
  if( fabs(s->P.pol.T) <= 0. || fabs(s->A.pol.T) <= 0. ) return 2;

  /* cell decay constants must be > 0 */

  /*
    first calc the cell He3 polarizations at datapt time with errs
    use s->(P,A).pol.(PHe,hr0,T) along with d->hr
    to calc current He3 polarizations and errs

    calc nsL*(1+-PHe3) = tau+- and its err which is +- independent
    use tau+- to calc correction coefs along with He3exp (expA, expP)

    finally can calculate transmission factors and errors
    and then coefs and errors
  */


  /* efficiency stuff */
  etPsq = s->eP.teff ;
  etPsq *= etPsq ;
  etPesq = s->eP.terr ;
  etPesq *= etPesq ;
  efP = s->eP.feff ;
  efPesq = s->eP.ferr ;
  efPesq *= efPesq ;

  etAsq = s->eA.teff ;
  etAsq *= etAsq ;
  etAesq = s->eA.terr ;
  etAesq *= etAesq ;
  efA = s->eA.feff ;
  efAesq = s->eA.ferr ;
  efAesq *= efAesq ;

  /*
    Now set up loop to calculate the coefs and errsq
                    Pol  Ana
    index 0  + +    OFF  OFF
    index 1  - -    ON   ON
    index 2  + -    ON   OFF
    index 3  - +    OFF  ON
    but first compute the transport factors, eAalp, ePbet
  */

  for( i=0 ; i<4 ; i++ ) {
    /* i is row index so compute alp and bet for i */
    alp = alpmu[i] ;
    alp1sq = (1.- alp)*(1.- alp) ;
    eAFalp = (1.- alp)*s->eA.feff - 1. ;
    eAFalpsq = eAFalp*eAFalp ;
    eAalp = s->eA.teff*eAFalp ;
    eAalpsq = eAalp*eAalp ;

    bet = betnu[i] ;
    bet1sq = (1.- bet)*(1.- bet) ;
    ePFbet = (1.- bet)*s->eP.feff - 1. ;
    ePFbetsq = ePFbet*ePFbet ;
    ePbet = s->eP.teff*ePFbet ;
    ePbetsq = ePbet*ePbet ;

    /* use the dataPt time for each measurement (row) */

    /*
      given cell pointer, data-time-stamp and the wavelength for that cell
      compute the symmetric and antissymmetric transmission factors and the
      error sum.
      Now transfrac returns the + - transmission factors and errors.
    */

    if( (ierr = transfac(&(s->A), d->sec[i], d->lambF,
			 &tpA, &tmA, &tpAesq, &tmAesq)) ) return ierr ;
    if( (ierr = transfac(&(s->P), d->sec[i], d->lambI,
			 &tpP, &tmP, &tpPesq, &tmPesq)) ) return ierr ;
    taA = tpA - tmA ;
    tsA = tpA + tmA ;
    tAesq = tpAesq + tmAesq ;

    taP = tpP - tmP ;
    tsP = tpP + tmP ;
    tPesq = tpPesq + tmPesq ;

    /*
      the only error factors that depend on mu nu
      are the same factors as in C so we can compute them before the
      column j loop
    */

    etermA = (ePFbetsq*etPesq + bet1sq*etPsq*efPesq)*taP*taP +
      (1.+ePbetsq)*tPesq ;
    etermP = (eAFalpsq*etAesq + alp1sq*etAsq*efAesq)*taA*taA +
      (1.+eAalpsq)*tAesq ;

    for( j=0 ; j<4 ; j++ ) {
      /* j is column index so compute mu nu */
      mu = alpmu[j] ;
      nu = betnu[j] ;

      cA = tsA - mu*eAalp*taA ;
      cAsq = cA*cA ;
      cP = tsP - nu*ePbet*taP ;
      cPsq = cP*cP ;

      d->C[i][j] = cA*cP/8. ; /* include 1/2 from N/2 */
      d->Cesq[i][j] = (cAsq*etermA + cPsq*etermP)/64. ;
    }
  }

  return 0;
}

static double correctionCoef(He3CELLexp *ex, double tau)
{
  
  He3CELL *cell ;
  expResol *res ;
  double L, R, r ;
  double sw, sV, sH ;

  if( ex == NULL ) return 1.;
  cell = &(ex->cell) ;
  res = &(ex->res) ;
  L = cell->L ;
  if( L <= 0. ) return 1.;
  R = cell->R ;
  if( R <= 0. ) R = 1.e12 ;
  r = res->usedRadius ;
  sw = res->waveRelWidth ;
  sV = res->angleVwidth ;
  sH = res->angleHwidth ;
  return
    1. + tau*tau*sw*sw/2. - tau*(1 - L/R/2. + r*r/(R*R)/2.)*(sV*sV + sH*sH)/2.;
}

/*
  Correct a PolarizedBeamDatapoint.
  CountRates have been measured for some set of the standard polarize beam
  cross-sections, ++ -- +- -+, and the correction coeficients and their errors
  have been determined (time dependent).
  order here is ++=PoffAoff --=PonAon +-=PonAoff -+=PoffAon
  Now solve the system of eqs
  Y(CountRates)+-sigY = CorrectionMatrix+-sigC * S(cross-sections)
  for the S cross-sections including error propagation given the
  errors in CountRates and in the CorrectionMatrix elements.
  Do this by Gaussian elimination
  S and Serr are solved for.
  N.B. all errs are sigsq
  for this version meant for pol beam CS extraction from up to 4 measurements
  so limit N <= 4
  err return codes
  0 success
  1 NULL passed storage
  2 Nfree != Nactive
  3 got zero diagonal CorrectionMatrix value
 */


static int PBcorrect(PBdatapt *d)
{
  PBdatapt D;

  int N, m, row, col, i, im, j, jm, itemp;
  int indx[4];
  double Amax, temp, tempsq, temp1, temp2, sigsq;

  if( d == NULL ) return 1;

  if( d->Nactive < 1 || d->Nactive > 4 || d->Nactive != d->Nfree ) return 2;
  /* also check the activeEq and freeS indices for invalid range */
  for( i=0 ; i<d->Nactive ; i++ ) {
    if( d->activeEq[i] < 0 || d->activeEq[i] > 3 ||
	d->freeS[i] < 0 || d->freeS[i] > 3 ) return 2;
  }

  /* pack the passed datapt into the local datastruct D */
  for( i=0 ; i<d->Nactive ; i++ ) {
    im = d->activeEq[i];
    D.Y[i] = d->Y[im];
    D.Yesq[i] = d->Yesq[im];
    for( j=0 ; j<d->Nfree ; j++ ) {
      jm = d->freeS[j];
      D.C[i][j] = d->C[im][jm];
      D.Cesq[i][j] = d->Cesq[im][jm];
    }
  }
  N = d->Nactive;

  for( col=0 ; col<N ; col++ ) indx[col] = col ;

  for( col=0 ; col<N-1 ; col++ ) {
    Amax = fabs(D.C[col][col]);
    m = col;
    for( row=col+1 ; row<N ; row++ ) {
      if( Amax < fabs(D.C[row][col]) ) {
	Amax = D.C[row][col];
	m = row;
      }
    }
    /* swap rows if m != col */

    if( m != col ) {
      for( i=col ; i<N ; i++ ) {
	temp = D.C[col][i] ;
	D.C[col][i] = D.C[m][i] ;
	D.C[m][i] = temp ;
	temp = D.Cesq[col][i] ;
	D.Cesq[col][i] = D.Cesq[m][i] ;
	D.Cesq[m][i] = temp ;
      }
      temp = D.Y[col] ;
      D.Y[col] = D.Y[m] ;
      D.Y[m] = temp ;
      temp = D.Yesq[col] ;
      D.Yesq[col] = D.Yesq[m] ;
      D.Yesq[m] = temp ;

      itemp = indx[col] ;
      indx[col] = indx[m] ;
      indx[m] = itemp ;
    }

    if( D.C[col][col] == 0. ) { return 3; }

    for( row=col+1 ; row<N ; row++ ) {

      temp = -D.C[row][col]/D.C[col][col] ;
      tempsq = temp*temp ;
      temp1 = D.C[row][col]*D.C[row][col] ;
      temp2 = D.C[col][col]*D.C[col][col] ;
      sigsq = (D.Cesq[row][col] + D.Cesq[col][col]*temp1/temp2)/temp2 ;

      for( i=col ; i<N ; i++ ) {
	D.C[row][i] += temp*D.C[col][i] ;
	D.Cesq[row][i] += tempsq*D.Cesq[col][i] +
	  D.C[col][i]*D.C[col][i]*sigsq ;
      }
      D.Y[row] += temp*D.Y[col];
      D.Yesq[row] += tempsq*D.Yesq[col] + D.Y[col]*D.Y[col]*sigsq ;
    }
  }

  /* now obtain the soln */
  for( row=N-1 ; row>=0 ; row-- ) {
    D.S[row] = D.Y[row] ;
    D.Sesq[row] = D.Yesq[row] ;
    for( i=row+1 ; i<N ; i++ ) {
      tempsq = D.C[row][i]*D.C[row][i] ;
      D.S[row] -= D.C[row][i]*D.S[i] ;
      D.Sesq[row] += tempsq*D.Sesq[i] + D.S[i]*D.S[i]*D.Cesq[row][i] ;
    }
    D.S[row] /= D.C[row][row] ;
    tempsq = D.C[row][row]*D.C[row][row] ;
    D.Sesq[row] = 
      (D.Sesq[row] + D.Cesq[row][row]*D.S[row]*D.S[row]/tempsq)/tempsq ;
  }

  /* move soln back to original indx positions using Y as temp storage */
  for( row=0 ; row<N ; row++ ) {
    D.Y[row] = D.S[row] ; D.Yesq[row] = D.Sesq[row] ;
  }
  for( row=0 ; row<N ; row++ ) {
    D.S[row] = D.Y[indx[row]] ; D.Sesq[row] = D.Yesq[indx[row]] ;
  }

  /* move solutions into freeS indices */
  for( row=0 ; row<N ; row++ ) {
    d->S[d->freeS[row]] = D.S[row];
    d->Sesq[d->freeS[row]] = D.Sesq[row];
  }
  return 0;
}

/*
  applying constraint modifies PBdatapt
  reduces Nfree by one
*/

static int applyConstraint(PBdatapt *d, ConstraintEq *eq)
{
  int i, j, ok ;
  double Coef[4] ;
  double Cf ;

  if( d == NULL || eq == NULL ) return 1;
  if( eq->Nfree + 1 > d->Nfree ) return 2;

  /*
    also check that freeToConstrain and each eq->freeS index is
    in the list d->freeS
  */
  ok = 0 ;
  for( j=0 ; j<d->Nfree ; j++ ) {
    if( eq->freeToConstrain == d->freeS[j] ) { ok = 1 ; break ; }
  }
  if( ! ok ) return 3 ;

  for( i=0 ; i<eq->Nfree ; i++ ) {
    ok = 0 ;
    for( j=0 ; j<d->Nfree ; j++ ) {
      if( eq->freeS[i] == d->freeS[j] ) { ok = 1 ; break ; }
    }
    if( ! ok ) return 3 ;
  }

  /*
    Now move the eq->Coef into ordered locations in local Coef
    with others == 0
  */
  for( i=0 ; i<4 ; i++ ) Coef[i] = 0.;
  for( i=0 ; i<eq->Nfree ; i++ ) Coef[eq->freeS[i]] = eq->Coef[i] ;

  /*
    remove freeToConstrain index from d->freeS
  */
  for( i=0, j=0 ; i<d->Nfree ; i++ ) {
    if( d->freeS[i] == eq->freeToConstrain ) continue ;
    d->freeS[j] = d->freeS[i] ;
    j++ ;
  }
  d->Nfree-- ;

  /*
    Now compute the new coefs for each index in the new d->freeS list
    i.e. except for freeToConstrain (set that coef to zero for now)
    for each active equ in d->activeEq
  */
  for( i=0 ; i<d->Nactive ; i++ ) {
    for( j=0 ; j<d->Nfree ; j++ ) {
      Cf = Coef[d->freeS[j]] ;
      d->C[d->activeEq[i]][d->freeS[j]] += 
	Cf*d->C[d->activeEq[i]][eq->freeToConstrain] ;
      d->Cesq[d->activeEq[i]][d->freeS[j]] += 
	Cf*Cf*d->Cesq[d->activeEq[i]][eq->freeToConstrain] ;
    }
    d->C[d->activeEq[i]][eq->freeToConstrain] = 0. ;
    d->Cesq[d->activeEq[i]][eq->freeToConstrain] = 0. ;
  }
  return 0;
}

static int combineEqs(PBdatapt *d, int eq1, int eq2)
{
  int i, j, ok ;

  if( d == NULL ) return 1 ;
  /*
    First make sure eq1 and eq2 are in the d->activeEq list
  */
  if( eq1 == eq2 ) return 2 ;
  ok = 0 ;
  for( i=0 ; i<d->Nactive ; i++ ) {
    if( eq1 == d->activeEq[i] ) { ok = 1 ; break ; }
  }
  if( ! ok ) return 2 ;
  ok = 0 ;
  for( i=0 ; i<d->Nactive ; i++ ) {
    if( eq2 == d->activeEq[i] ) { ok = 1 ; break ; }
  }
  if( ! ok ) return 2 ;

  /*
    Now add equations eq1 and eq2 putting result into eq1
    and remove eq2 from activeEq list and reduce Nactive
  */

  /* add Y values */
  d->Y[eq1] += d->Y[eq2] ;
  d->Yesq[eq1] += d->Yesq[eq2] ;

  /* add coefs */
  for( i=0 ; i<4 ; i++ ) {
    d->C[eq1][i] += d->C[eq2][i] ;
    d->Cesq[eq1][i] += d->Cesq[eq2][i] ;
  }

  for( i=0, j=0 ; i<d->Nactive ; i++ ) {
    if( eq2 == d->activeEq[i] ) { continue ; }
    d->activeEq[j] = d->activeEq[i] ;
    j++ ;
  }
  d->Nactive-- ;
  return 0 ;
}

static int deleteEq(PBdatapt *d, int eq)
{
  int i, j, ok ;

  if( d == NULL ) return 1 ;
  /*
    First make sure eq is in the d->activeEq list
  */

  ok = 0 ;
  for( i=0 ; i<d->Nactive ; i++ ) {
    if( eq == d->activeEq[i] ) { ok = 1 ; break ; }
  }
  if( ! ok ) return 2 ;

  /*
    Now delete equation eq
  */

  for( i=0, j=0 ; i<d->Nactive ; i++ ) {
    if( eq == d->activeEq[i] ) { continue ; }
    d->activeEq[j] = d->activeEq[i] ;
    j++ ;
  }
  d->Nactive-- ;
  return 0 ;
}

static int constrainResult(PBdatapt *d)
{
  int i, j ;
  double Cf ;
  double *S, *Sesq ;
  ConstraintEq *eq ;

  if( d == NULL || eqs == NULL ) return 1;
  if( Nconstraint < 1 ) return 0;


  for( i=0 ; i<Nconstraint ; i++ ) {
    eq = eqs + i ;

    S = d->S+(eq->freeToConstrain) ;
    *S = 0. ;
    Sesq = d->Sesq+(eq->freeToConstrain) ;
    *Sesq = 0. ;
    for( j=0 ; j<eq->Nfree ; j++ ) {
      Cf = eq->Coef[j] ;
      *S += Cf*d->S[eq->freeS[j]] ;
      *Sesq += Cf*Cf*d->Sesq[eq->freeS[j]] ;
    }
  }
  return 0;
}

static double mygauss(double x, double H, double P, double FW)
{
  double z, zsq ;
  z = 2.*(x - P)/FW ;
  zsq = z*z ;
  return H*pow(0.5,zsq) ;
}

static double monFuncBT4PG(double E)
{
  /* from fitting two Gaussians to the BT4 correction data */
  return
    1.0025 + 
    mygauss(E, 0.7087,  8.6305, 5.2114) + 
    mygauss(E, 0.40319,12.488, 22.381) ;
}

static int monitorCorrect(PBdatapt *d, double (*monFunc)())
{
  /*
    The monFunc returns the correction factor
    (SUMn a_n/n)/a_1
    as a function of the monochromator energy in meV as its argument
    where a_n are the wavelength order fractions in the incident beam
    from the monochromator.
  */

  int i ;
  double EI, kI, cor, corsq ;

  if( d == NULL || monFunc == NULL ) return 1 ;
  kI = (TWOPI/d->lambI) ;
  EI = Dn*kI*kI ;
  cor = monFunc(EI) ;
  corsq = cor*cor ;
  for( i=0 ; i<4 ; i++ ) { d->Y[i] *= cor ; d->Yesq[i] *= corsq ; }
  return 0 ;
}
static int SIMmonitorCorrect(PBdatapt *d, double (*monFunc)())
{
  /*
    The monFunc returns the correction factor
    (SUMn a_n/n)/a_1
    as a function of the monochromator energy in meV as its argument
    where a_n are the wavelength order fractions in the incident beam
    from the monochromator.
    Here correct simulated counts lower because counting against
    fixed monitor total which includes lambdahalf
  */

  int i ;
  double EI, kI, cor, corsq ;

  if( d == NULL || monFunc == NULL ) return 1 ;
  kI = (TWOPI/d->lambI) ;
  EI = Dn*kI*kI ;
  cor = monFunc(EI) ;
  corsq = cor*cor ;
  for( i=0 ; i<4 ; i++ ) { d->Y[i] /= cor ; d->Yesq[i] /= corsq ; }
  return 0 ;
}


static int monitorCorrectPG(PBdatapt *d)
{
  return monitorCorrect(d, monFuncBT4PG) ;
}
static int SIMmonitorCorrectPG(PBdatapt *d)
{
  return SIMmonitorCorrect(d, monFuncBT4PG) ;
}

static int polmonitorCorrect(PBdatapt *d, He3CELLexp *pol, double (*monFunc)())
{
  /*
    This calcs the counts correction factor when using a beam monitor
    placed after the polarizer and assuming only second order contamination.
    The monFunc returns the correction factor
    (SUMn a_n/n)/a_1
    as a function of the monochromator energy in meV as its argument
    where a_n are the wavelength order fractions in the incident beam
    from the monochromator.
  */

  int i ;
  double EI, kI, lamb, lamb2, cor ;
  double a1, a2 ;
  double tp, tm, tpesq, tmesq, t1, t2 ;

  if( d == NULL || monFunc == NULL ) return 1 ;
  lamb = d->lambI ;
  lamb2 = lamb/2. ;
  kI = (TWOPI/lamb) ;
  EI = Dn*kI*kI ;
  cor = monFunc(EI) ;
  /* extract the order fractions from the standard corfac */
  a1 = 1./(2*cor - 1) ;
  a2 = 1 - a1 ;
  /* correct each of the four cross-section counts */
  for( i=0 ; i<4 ; i++ ) {
    if( ! transfac(pol, d->sec[i], lamb, &tp, &tm, &tpesq, &tmesq) ) return 1 ;
    t1 = (tp + tm)/2. ;
    if( ! transfac(pol, d->sec[i], lamb2, &tp, &tm, &tpesq, &tmesq) ) return 1 ;
    t2 = (tp + tm)/2. ;
    cor = a1*t1 + a2*t2/2. ;
    d->Y[i] *= cor ;
    d->Yesq[i] *= cor*cor ;
  }
  return 0 ;
}
static int SIMpolmonitorCorrect(PBdatapt *d, He3CELLexp *pol, double (*monFunc)())
{
  /*
    This calcs the counts correction factor when using a beam monitor
    placed after the polarizer and assuming only second order contamination.
    The monFunc returns the correction factor
    (SUMn a_n/n)/a_1
    as a function of the monochromator energy in meV as its argument
    where a_n are the wavelength order fractions in the incident beam
    from the monochromator.
    Correct SIM data to lower value when counting against monitor
  */

  int i ;
  double EI, kI, lamb, lamb2, cor ;
  double a1, a2 ;
  double tp, tm, tpesq, tmesq, t1, t2 ;

  if( d == NULL || monFunc == NULL ) return 1 ;
  lamb = d->lambI ;
  lamb2 = lamb/2. ;
  kI = (TWOPI/lamb) ;
  EI = Dn*kI*kI ;
  cor = monFunc(EI) ;
  /* extract the order fractions from the standard corfac */
  a1 = 1./(2*cor - 1) ;
  a2 = 1 - a1 ;
  /* correct each of the four cross-section counts */
  for( i=0 ; i<4 ; i++ ) {
    if( ! transfac(pol, d->sec[i], lamb, &tp, &tm, &tpesq, &tmesq) ) return 1 ;
    t1 = (tp + tm)/2. ;
    if( ! transfac(pol, d->sec[i], lamb2, &tp, &tm, &tpesq, &tmesq) ) return 1 ;
    t2 = (tp + tm)/2. ;
    cor = a1*t1 + a2*t2/2. ;
    d->Y[i] /= cor ;
    d->Yesq[i] /= cor*cor ;
  }
  return 0 ;
}

static int polmonitorCorrectPG(PBdatapt *d, He3CELLexp *pol)
{
  return polmonitorCorrect(d, pol, monFuncBT4PG) ;
}
static int SIMpolmonitorCorrectPG(PBdatapt *d, He3CELLexp *pol)
{
  return SIMpolmonitorCorrect(d, pol, monFuncBT4PG) ;
}

/*

  Given a PBsetup and a PBdatapt array that contains
  simulated PB data, compute coefs and underlying CS with errs.
  data set is in the PBdatapt array.
  The data are input from a columnar data file where
  and the PBcorrection requires columns for
  EI EF CountsOFF timestampOFFOFF CountsONON timestampONON CountsONOFF ...
  We will assume the time stamps are in UTC seconds.
  the first column is the time-stamp value in hours,
  and the next 4 columns are ++ +- -+ -- data per datapoint.
  This program then writes the calculated CS and err to stdout.

  Read in PBsetup info from 2 files one for polarizer and one for analyzer
  The files contain 8 lines each which are 4 header lines each followed
  by a data line
  line1 is    He3CELL: name tEmpty tEmpSlope Lcm Dcm Rcrv-cm nsL0-1.77A nsL0err
  line2 data for these entires
  line3 is He3CELLpol: PHe  PHeErr hrsBeam T  Terr startTimeSecs  startTimeString
  line4 is data for these
  line5 is     He3exp: waveRelWidth angleVwidth angleHwidth usedRadius
  line6 is data for these
  line7 is efficiency: transportEff terr flipperEff ferr
  line8 is data for these

  NB startTimeString should be converted to startTimeSecs under UNIX with cmnd
  date --date="startTimeString" +%s

  reader calculates nsL corrected and nsLerr putting them into PBsetup
  so program input is 3 filenames  on cmnd line:
  Pcellfile Acellfile DATAfile [flag1 [flag2]]
  If there is a lambdaFlag(L) the init and final Energies are wavelengths in
  Angstroms instead of the default meV
  There may also be a flag indicating that counts were obtained using
  a normalizing monitor place after the polarizer cell.
  This is not very good since such a monitor will also count lambda/2.
  If no monitor is available before the polarizing cell, it is probably best
  to count against time, but then how to account for flux changes with
  changing incident energy. We would need to know the fraction of lambda/2
  as a function of incident energy
*/




static int PBdefineCell(char *filename, He3CELLexp *cellexp, efficiency *eff)
{
  int i, nscan ;
  double cellfac ;
  FILE *fp ;
  char buf[512] ;

  He3CELL *cell ;
  He3CELLpol *pol ;
  expResol *exper ;

  if( (fp = fopen(filename, "r")) == NULL ) {
    printf("failed to open %s\n", filename) ;
    return 1 ;
  }
  for( i=0 ; i<4 ; i++ ) {
    if( (fgets(buf, 511, fp)) == NULL ) {
      printf("failed to read from %s\n", filename) ;
      fclose(fp) ;
      return 2 ;
    }
    if( (fgets(buf, 511, fp)) == NULL ) {
      printf("failed to read from %s\n", filename) ;
      fclose(fp) ;
      return 2 ;
    }
    if( i == 0 ) {
      cell = &(cellexp->cell) ;
      if( sscanf(buf,"%s %lf %lf %lf %lf %lf %lf %lf",
		 cell->name,&(cell->tEmpty),&(cell->tEmptySlope),
		 &(cell->L),&(cell->D),&(cell->R),
		 &(cell->nsL0),&(cell->nsL0err)) < 8 ) {
	printf("failed to read cell data from %s\n", filename) ;
	fclose(fp) ;
	return 3 ;
      }
    } else if( i == 1 ) {
      pol = &(cellexp->pol) ;
      nscan = sscanf(buf,"%lf %lf %lf %lf %lf %ul",
		     &(pol->PHe),&(pol->PHeErr),&(pol->hrsBeam),
		     &(pol->T),&(pol->Terr),&(pol->startSecs)) ;
      if( nscan < 6 ) {
	if( nscan < 5 ) printf("failed to read cell data from %s\n", filename) ;
	else printf("failed to read CELL startSecs from %s\n", filename) ;
	fclose(fp) ;
	return 3 ;
      }
      /*
	convert the startDate to seconds since Jan 1 1971
	replace - with space and pass the string to system call
	date --date="startDate" +%s
	to return the time in seconds
	This date command syntax is probably Unix specific
	NB unsigned long max value is about 4294967295
	Jan 1 2008 is                       1199163600
	difference is                       3095803695
	so run out of data space in about 100 years
	ICE may add fraction of second to the time stamp but
	the ul format will strip that fraction.
      */
    } else if( i == 2 ) {
      exper = &(cellexp->res) ;
      if( sscanf(buf,"%lf %lf %lf %lf",
		 &(exper->waveRelWidth),
		 &(exper->angleVwidth),&(exper->angleHwidth),
		 &(exper->usedRadius)) < 4 ) {
	printf("failed to read cell data from %s\n", filename) ;
	fclose(fp) ;
	return 3 ;
      }
    } else if( i == 3 ) {
      if( sscanf(buf,"%lf %lf %lf %lf",
		 &(eff->teff),&(eff->terr),
		 &(eff->feff),&(eff->ferr)) < 4 ) {
	printf("failed to read cell data from %s\n", filename) ;
	fclose(fp) ;
	return 3 ;
      }
    }
  }
  fclose(fp) ;
  /* compute curvature corrected nsL and nsLerr */
  cellfac = 1. ;
  if( cell->R > 0. && cell->L > 0. )
    cellfac *= (1. - exper->usedRadius*exper->usedRadius/cell->R/cell->L) ;
  cell->nsL = cellfac*cell->nsL0 ;
  cell->nsLerr = cellfac*cell->nsL0err ;
  return 0 ;
}
static int PBdefinePolarizer(char *filename)
{
  return PBdefineCell(filename, &(s.P), &(s.eP)) ;
}
static int PBdefineAnalyzer(char *filename)
{
  return PBdefineCell(filename, &(s.A), &(s.eA)) ;
}



static int PBsetflags(PBflags *flgs)
{
  /* set the global flag structure from the called one */
  flags = *flgs ;
  return 0 ;
}


/*
  PBsim assumes Pcell and Acell data have already been read into
  global by using PBcorrectData(Pcell, Acell, ...)

  PBsim first reads EI EF and cross-section values from filename
  then calculates SIM count-rates.
  then corrects the count rates back to cross-section values.
  The SIM data goes into filename.SIM
  and the corrected CS to filname.CS
  any ctrls should be in filename.CTRL
*/

int PBsim(char *filename)
{

  PBdatapt d ;

  int i, ierr, j, k, n, dataOK, Npt ;
  char buf[512], labl[64] ;
  char outfile[512] ;
  double S[4] ;
  double dval, d1, d2 ;
  double iout, fout, err[4] ;
  double hrsecs, HrStep ;
  double EI, EF, EIstart, EIstep, EFstart, EFstep ;

  unsigned long ulhrs ;
  FILE *fp, *fpout, *fpcorrect ;

  if( filename == NULL ) return 1 ;

  strcpy(outfile, filename) ;
  strcat(outfile, ".CTRL") ;
  ierr = PBreadflags(outfile) ;
  if( ierr == 1 ) {
    printf("NO .CTRL file\n") ;
  } else if( ierr > 1 ) {
    printf("ERROR in PBreadflags = %d\n", ierr) ;
    return 2 ;
  }

  /* make sure flags are applied to constraints */
  constraintTOeqs() ;


  if( (fp = fopen(filename, "r")) == NULL ) {
    printf("failed to open %s\n", filename) ;
    return 2 ;
  }
  n = 0 ;
  /*
    simulation data files
    contain:
    Npts: Npt
    TimestepHrs: HRstep
    EIseq: EIstart EIstep
    EFseq: EFstart EFstep
    SPoffAoff SPonAon SPonAoff SPoffAon
    SPoffAoff SPonAon SPonAoff SPoffAon
    ...

    This simulation assumes that a set of cross-sections is 
    measured at the same time.
    hr gets converted to seconds and added to the POLARIZER
    utc seconds to generate the necessary utc for each CS msrment
   */

  dataOK = 0 ;
  while( fgets(buf, 511, fp) ) {
    if( buf[0] == '#' ) continue ;
    if( n == 0 ) {
      if( sscanf(buf,"%s %d", labl, &Npt) < 2 ) {
	printf("failed to read Npts from %s\n", filename);
	fclose(fp) ;
	return 3 ;
      }
    } else if( n == 1 ) {
      if( sscanf(buf,"%s %lf", labl, &HrStep) < 2 ) {
	printf("failed to read HrStep from %s\n", filename);
	fclose(fp) ;
	return 3 ;
      }
    } else if( n == 2 ) {
      if( sscanf(buf,"%s %lf %lf", labl, &EIstart, &EIstep) < 3 ) {
	printf("failed to read EIstart EIstep from %s\n", filename);
	fclose(fp) ;
	return 3 ;
      }
    } else if( n == 3 ) {
      if( sscanf(buf,"%s %lf %lf", labl, &EFstart, &EFstep) < 3 ) {
	printf("failed to read EFstart EFstep from %s\n", filename);
	fclose(fp) ;
	return 3 ;
      }
      dataOK = 1 ;
      break ;
    }
    n++ ;
  }
  if( ! dataOK ) {
    printf("sim datafile requires 3 start lines followed by data lines\n") ;
    fclose(fp) ;
    return 3 ;
  }
  dataOK = 0 ;
  n = 0 ;

  strcpy(outfile, filename) ;
  strcat(outfile, ".SIM") ;
  if( (fpout = fopen(outfile, "w")) == NULL ) {
    printf("failed to open simulation data file %s\n", outfile) ;
    return 2 ;
  }
  strcpy(outfile, filename) ;
  strcat(outfile, ".CS") ;
  if( (fpcorrect = fopen(outfile, "w")) == NULL ) {
    printf("failed to open simulation correctCS %s\n", outfile) ;
    return 2 ;
  }


  while( fgets(buf, 511, fp) ) {
    if( buf[0] == '#' ) continue ;
    if( sscanf(buf,"%lf %lf %lf %lf",S,S+1,S+2,S+3) < 4 ) {
      printf("failed to read 4 CS values from %s\n", filename);
      fclose(fp) ;
      fclose(fpout) ;
      return 3 ;
    }
    for( i=0 ; i<4 ; i++ ) d.S[i] = S[i] ;
    dataOK = 1 ;
    EI = EIstart + n*EIstep ;
    d.lambI = TWOPI/sqrt(EI/Dn) ;
    EF = EFstart + n*EFstep ;
    d.lambF = TWOPI/sqrt(EF/Dn) ;

    /* convert the sim hrs into utc offset from POLARIZER utc */
    hrsecs = n*HrStep*3600. ;
    ulhrs = hrsecs ;
    for( i=0 ; i< 4 ; i++ ) {
      d.sec[i] = s.P.pol.startSecs + ulhrs ;
    }

    if( (ierr = PBcoef(&d, &s)) > 0 ) {
      printf("error in PBcoef = %d for point index= %d\n", ierr, i) ;
      fclose(fp) ;
      fclose(fpout) ;
      return 4 ;
    }

    /* calc the count rates and then deviate them if requested */
    for( j=0 ; j<4 ; j++ ) {
      d.Y[j] = 0. ;
      for( k=0 ; k<4 ; k++ ) {
	d.Y[j] += d.C[j][k]*d.S[k] ;
      }
      d.Y[j] *= fabs(flags.SimFlux) ; /* put the factor of 1/2 into coefs */
      d.Yesq[j] = d.Y[j] ;
      err[j] = 1. ;
      if( d.Yesq[j] > 0 ) err[j] = sqrt(d.Yesq[j]) ;
      if( flags.SimDeviate ) d.Y[j] = poisson(0, d.Y[j], err[j], &dval, d1, d2) ;
      if( flags.MonitorCorrect ) SIMmonitorCorrectPG(&d) ;
      else if( flags.PolMonitorCorrect ) SIMpolmonitorCorrectPG(&d, &(s.P)) ;
    }

    fprintf(fpout, "%9g %9g  ", EI,EF) ;
    for( j=0 ; j<4 ; j++ ) fprintf(fpout, "  %11g", d.Y[j]) ;
    for( j=0 ; j<4 ; j++ ) fprintf(fpout, "  %11g", err[j]) ;
    for( j=0 ; j<4 ; j++ ) fprintf(fpout, "  %11d", d.sec[j]) ;
    fprintf(fpout, "\n") ;


    if( (ierr = PBcorrectDatapt(&d)) ) {
      printf("ERROR in PBcorrectDatapt = %d\n", ierr) ;
    }
    for( j=0 ; j<4 ; j++ ) {
      err[j] = 1. ;
      if( d.Sesq[j] > 0 ) err[j] = sqrt(d.Sesq[j]) ;
    }
    if( flags.SimFlux < 1 ) flags.SimFlux = 1 ;
    for( j=0 ; j<4 ; j++ ) fprintf(fpcorrect, "%11g  ",d.S[j]/flags.SimFlux);
    for( j=0 ; j<4 ; j++ ) fprintf(fpcorrect, "%11g  ",err[j]/flags.SimFlux);
    fprintf(fpcorrect, "\n") ;

    n++ ;
    if( n > Npt ) {
      fclose(fp) ;
      fclose(fpout) ;
      fclose(fpcorrect) ;
      return 0 ;
    }
  }
  /* finished reading S values but can continue to Npt with last S values */
  fclose(fp) ;

  for( i=n ; i<Npt ; i++ ) {
    /* put the last read S values into d */
    for( j=0 ; j<4 ; j++ ) d.S[j] = S[j] ;

    EI = EIstart + i*EIstep ;
    d.lambI = TWOPI/sqrt(EI/Dn) ;
    EF = EFstart + i*EFstep ;
    d.lambF = TWOPI/sqrt(EF/Dn) ;

    /* convert the sim hrs into utc offset from POLARIZER utc */
    hrsecs = i*HrStep*3600. ;
    ulhrs = hrsecs ;
    for( j=0 ; j< 4 ; j++ ) {
      d.sec[j] = s.P.pol.startSecs + ulhrs ;
    }

    if( (ierr = PBcoef(&d, &s)) > 0 ) {
      printf("error in PBcoef = %d for point index= %d\n", ierr, i) ;
      fclose(fpout) ;
      return 4 ;
    }

    /* calc the count rates and then deviate them if requested */
    for( j=0 ; j<4 ; j++ ) {
      d.Y[j] = 0. ;
      for( k=0 ; k<4 ; k++ ) {
	d.Y[j] += d.C[j][k]*d.S[k] ;
      }
      d.Y[j] *= fabs(flags.SimFlux) ; /* put the factor of 1/2 into coefs */
      d.Yesq[j] = d.Y[j] ;
      err[j] = 1. ;
      if( d.Yesq[j] > 0 ) err[j] = sqrt(d.Yesq[j]) ;
      if( flags.SimDeviate ) d.Y[j] = poisson(0, d.Y[j], err[j], &dval, d1, d2) ;
      if( flags.MonitorCorrect ) SIMmonitorCorrectPG(&d) ;
      else if( flags.PolMonitorCorrect ) SIMpolmonitorCorrectPG(&d, &(s.P)) ;
    }

    fprintf(fpout, "%9g %9g  ", EI,EF) ;
    for( j=0 ; j<4 ; j++ ) fprintf(fpout, "  %11g", d.Y[j]) ;
    for( j=0 ; j<4 ; j++ ) fprintf(fpout, "  %11g", err[j]) ;
    for( j=0 ; j<4 ; j++ ) fprintf(fpout, "  %11d", d.sec[j]) ;
    fprintf(fpout, "\n") ;


    if( (ierr = PBcorrectDatapt(&d)) ) {
      printf("ERROR in PBcorrectDatapt = %d\n", ierr) ;
    }
    for( j=0 ; j<4 ; j++ ) {
      err[j] = 1. ;
      if( d.Sesq[j] > 0 ) err[j] = sqrt(d.Sesq[j]) ;
    }
    if( flags.SimFlux < 1 ) flags.SimFlux = 1 ;
    for( j=0 ; j<4 ; j++ ) fprintf(fpcorrect, "%11g  ",d.S[j]/flags.SimFlux);
    for( j=0 ; j<4 ; j++ ) fprintf(fpcorrect, "%11g  ",err[j]/flags.SimFlux);
    fprintf(fpcorrect, "\n") ;

  }
  fclose(fpout) ;
  fclose(fpcorrect) ;
  return 0 ;
}

static double poisson(int idum, double ave, double sig, double *distval,
		double d1, double d2)
{
  /*
    p(x) = e(-lam) lam^x/x!
    e(-lam) SUM(x=0,inf) lam^x/x! = 1
    so pick a uniform deviate and find the x where SUM crosses e(lam)
  */ 
  double scale, mean, dval ;
  double udev ;
  double sum, fac, lam, target ;
  int i ;

  if ( sig <= 0. || ave <= 0. ) { return 0. ; }
  scale = ave/sig/sig ;
  /*
    mean is the lam in the distribution function
    we use ave = lam
    and    sig^2 = lam
    so sig is the first moment of the distribution
    and mean = lam = (ave/sig)^2
  */
  mean = scale*ave ;
  if ( mean >= 10. )
    {
      dval = gaussian(0, ave, sig, distval, d1, d2) ;
      return dval ;
    }
  udev = uniform(0, d1, d2, distval, 0., 1.) ;
  sum = 0. ;
  fac = 1. ;
  lam = 1. ;
  target = udev*exp(mean) ;
  for ( i=1 ; i<100 ; i++ )
    {
      sum += lam/fac ;
      if ( sum >= target ) { return (i-1)/scale ; }
      lam *= mean ;
      fac *= i ;
    }
  return 0. ;
}

static double gaussian(int idum, double ave, double sig, double *distval,
		double d1, double d2)
{
  double v1, v2, r, fac ;
  static int already = 0 ;
  static double deviate[2], dev, dv ;
  static double norm = 0.39894228 ;

  if( sig <= 0. ) { *distval = 0. ;  return 0. ; }
  if( already ) 
    {
      already = 0 ;
      dev = deviate[1] ;
      *distval = (norm/sig) * exp(-0.5*dev*dev) ;
      return ave + sig*dev ;
    }
  do
    {
      v1 = uniform(0, d1, d2, &dv, -1., 1.) ;
      v2 = uniform(0, d1, d2, &dv, -1., 1.) ;
      r = v1*v1 + v2*v2 ;
    }
  while( r >= 1. ) ;
  fac = sqrt( -2.*log(r)/r) ;
  deviate[0] = v1*fac ;
  deviate[1] = v2*fac ;
  already = 1 ;
  dev = deviate[0] ;
  *distval = (norm/sig) * exp(-0.5*dev*dev) ;
  return ave + sig*dev ;
}

static double uniform(int idum, double d1, double d2, double *distval,
	       double min, double max)
{
  int i, j ;
  double ran, range ;

  static double r[97] ;
  static int m1  = 259200 ;
  static int ia1 =   7141 ;
  static int ic1 =  54773 ;
  static double rm1 = 3.8580247e-6 ;
  static int m2  = 134456 ;
  static int ia2 =   8121 ;
  static int ic2 =  28411 ;
  static double rm2 = 7.4373773e-6 ;
  static int m3  = 243000 ;
  static int ia3 =   4561 ;
  static int ic3 =  51349 ;

  static int iff = 0 ;

  static int ix1, ix2, ix3 ;


  range = max - min ;
  if( range <= 0. ) { *distval = 0. ; return 0. ; }

  if( idum < 0 || iff == 0 )
    {
      iff = 1 ;
      ix1 = (ic1-idum) % m1 ;
      ix1 = (ia1*ix1+ic1) % m1 ;
      ix2 = ix1 % m2 ;
      ix1 = (ia1*ix1+ic1) % m1 ;
      ix3 = ix1 % m3 ;
      for( i=0 ; i<97 ; i++ )
	{
	  ix1 = (ia1*ix1+ic1) % m1 ;
	  ix2 = (ia2*ix2+ic2) % m2 ;
	  r[i] = (ix1 + ix2*rm2)*rm1 ;
	}
    }

  ix1 = (ia1*ix1+ic1) % m1 ;
  ix2 = (ia2*ix2+ic2) % m2 ;
  ix3 = (ia3*ix3+ic3) % m3 ;
  j = (97*ix3)/m3 ;
  ran = r[j] ;
  r[j] = (ix1 + ix2*rm2)*rm1 ;
  *distval = 1/range ;
  return min + range*ran ;
}

int PBreadflags(char *filename)
{
  /*
    process cmnds from filename to update the global flags
    R = reset to default flags
    C index  dpp dmm dpm dmp    Sindex = 0 + dpp*Spp + dmm*Smm + dpm*Spm + dmp*Smp
       NB the Sxx corresponding to Sindex is assumed to have dxx=0
       constrain a cross-section
       to disable the constraint on Sindex use
       default is no constraints
       NB indices are 1-4 corresponding to pp mm pm mp
    C -index
    A Cindex1 Cindex2    add two Counts equations
    B Cindex1 Cindex2    add two other Counts equations
    D Cindex1 Cindex2 ...  delete Count equations         deflt all equations enable
    D -Cindex1 ...  to undelete Count equations
    M intflag   flag for monitorCorrection                deflt OFF
    P intflag   flag for monitorafterPolarizerCorrection  deflt OFF
    F int       set simFlux        deflt 100000
    E int       flag sim deviates  deflt OFF
  */

  int k, nr, si, ic[4] ;
  char buf[512] ;
  double cf[4] ;
  FILE *fp ;

  if( (fp = fopen(filename, "r")) == NULL ) return 1 ;

  while( fgets(buf, 511, fp) ) {
    if( buf[0] == '#' ) continue ;
    if( buf[0] == 'C' || buf[0] == 'c' ) {
      /* constraint */
      if( (nr = sscanf(buf+1,"%d %lf %lf %lf %lf",
		       &si, cf, cf+1, cf+2, cf+3)) < 1 ) continue;
      if( abs(si) < 1 || abs(si) > 4 ) {
	printf("S index out of range, must be 1-4\n") ;
	fclose(fp) ;
	return 2 ;
      }
      if( si < 0 ) {
	flags.Sconstrain[abs(si)-1] = 0 ;
      } else {
	if( nr < 5 ) {
	  printf("set a constraint requires Sindex dpp dmm dpm dmp\n") ;
	  fclose(fp) ;
	  return 2 ;
	}
	flags.Sconstrain[si-1] = 1 ;
	if( si == 1 ) {
	  for( k=0 ; k<4 ; k++ ) flags.Spp[k] = cf[k] ;
	  flags.Spp[0] = 0. ;
	} else if( si == 2 ) {
	  for( k=0 ; k<4 ; k++ ) flags.Smm[k] = cf[k] ;
	  flags.Smm[1] = 0. ;
	} else if( si == 3 ) {
	  for( k=0 ; k<4 ; k++ ) flags.Spm[k] = cf[k] ;
	  flags.Spm[2] = 0. ;
	} else if( si == 4 ) {
	  for( k=0 ; k<4 ; k++ ) flags.Smp[k] = cf[k] ;
	  flags.Smp[3] = 0. ;
	}
      }

    } else if( buf[0] == 'A' || buf[0] == 'a' ) {
      /* add1 equations */
      if( (nr = sscanf(buf+1,"%d %d", ic, ic+1)) < 2 ) continue;
      if( ic[0] < 1 || ic[0] > 4 || ic[1] < 1 || ic[1] > 4 ) {
	printf("C index out of range, must be 1-4\n") ;
	fclose(fp) ;
	return 2 ;
      }
      for( k=0 ; k<4 ; k++ ) flags.CountsAdd1[k] = 0 ;
      flags.CountsAdd1[ic[0]-1] = 1 ;
      flags.CountsAdd1[ic[1]-1] = 1 ;
      
    } else if( buf[0] == 'B' || buf[0] == 'b' ) {
      /* add2 equations */
      if( (nr = sscanf(buf+1,"%d %d", ic, ic+1)) < 2 ) continue;
      if( ic[0] < 1 || ic[0] > 4 || ic[1] < 1 || ic[1] > 4 ) {
	printf("C index out of range, must be 1-4\n") ;
	fclose(fp) ;
	return 2 ;
      }
      for( k=0 ; k<4 ; k++ ) flags.CountsAdd2[k] = 0 ;
      flags.CountsAdd2[ic[0]-1] = 1 ;
      flags.CountsAdd2[ic[1]-1] = 1 ;
      
    } else if( buf[0] == 'D' || buf[0] == 'd' ) {
      if( (nr = sscanf(buf+1,"%d %d %d", ic, ic+1, ic+2)) < 1 ) continue;
      for( k=0 ; k<nr ; k++ ) {
	if( abs(ic[k]) < 1 || abs(ic[k]) > 4 ) {
	  printf("C index out of range, must be 1-4\n") ;
	  fclose(fp) ;
	  return 2 ;
	}
	if( ic[k] > 0 ) flags.CountsEnable[ic[k]-1] = 0 ;
	else flags.CountsEnable[abs(ic[k])-1] = 1 ;
      }
      
    } else if( buf[0] == 'M' || buf[0] == 'm' ) {
      if( (nr = sscanf(buf+1,"%d", ic)) < 1 ) continue;	
      if( ic[0] ) flags.MonitorCorrect = 1 ;
      else flags.MonitorCorrect = 0 ;
    } else if( buf[0] == 'P' || buf[0] == 'p' ) {
      if( (nr = sscanf(buf+1,"%d", ic)) < 1 ) continue;	
      if( ic[0] ) flags.PolMonitorCorrect = 1 ;
      else flags.PolMonitorCorrect = 0 ;
    } else if( buf[0] == 'D' || buf[0] == 'd' ) {
      if( (nr = sscanf(buf+1,"%d", ic)) < 1 ) continue;	
      if( ic[0] ) flags.Debug = 1 ;
      else flags.Debug = 0 ;
    } else if( buf[0] == 'F' || buf[0] == 'f' ) {
      if( (nr = sscanf(buf+1,"%d", ic)) < 1 ) continue;	
      if( ic[0] > 0 ) flags.SimFlux = ic[0] ;
    } else if( buf[0] == 'E' || buf[0] == 'e' ) {
      if( (nr = sscanf(buf+1,"%d", ic)) < 1 ) continue;	
      if( ic[0] ) flags.SimDeviate = 1 ;
      else flags.SimDeviate = 0 ;
    }
  }
  
  fclose(fp) ;
  return 0 ;
}
