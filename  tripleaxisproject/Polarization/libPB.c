/*
  C library for correcting polarized beam data using He-3 CELLS.
  Main entry is PBcorrectData.
  Modified PB structures to hold multiple setups to accomodate
  new philosophy of putting all He3 cell info for an experiment
  into a single file with one CELL per line as grabbed from the
  He3logger spreadsheet.
*/

#include <stdlib.h>
#include <stdio.h>
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
  double lambI, lambF ;
  double Y[4], Yesq[4], Yr[4], Yresq[4] ;
  double C[4][4], Cesq[4][4] ;
  double S[4], Sesq[4];
  double R ;
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
  double xsigsq;         /* std-dev-Lambda/Lambda squared */
  double vsigsq;         /* vertical angular resolution factor squared */
  double hsigsq;         /* horizontal angular resolution factor squared */
  double angcor;         /* part of linear tau correction coef */
  double t1, t2;         /* linear and quadratic tau correction coefs */
  double beamArea;
  double usedRadius;     /* eff beam radius in cm for cell */
} expResol;

typedef struct {
  double teff;  /* transport efficiency */
  double terr;
  double feff;  /* flipper efficiency */
  double ferr;
} efficiency;

typedef struct {
  int PorA ;
  He3CELL cell;
  He3CELLpol pol;
  expResol res;
  efficiency eff;
} He3CELLexp;

typedef struct {
  He3CELLexp P;
  He3CELLexp A;
} PBsetup;

typedef struct {
  int cells[2][4] ;
} PBcells ;

double Dn = 2.072141789 ;
double TWOPI = 6.283185307 ;

static double correctionCoef(He3CELLexp *ex, double tau) ;
static int transfac(He3CELLexp *ex, unsigned long secs, double lambda,
	     double *tp, double *tm, double *tpesq, double *tmesq) ;
static int PBcoef(PBdatapt *d) ;
static int PBcorrect(PBdatapt *d) ;
static int PBcorrectDatapt(PBdatapt *d) ;
static findcellsFORdatapt(PBdatapt *d, PBcells *c) ;
static int applyConstraint(PBdatapt *d, ConstraintEq *eq) ;
static int combineEqs(PBdatapt *d, int eq1, int eq2) ;
static int deleteEq(PBdatapt *d, int eq) ;
static int constrainResult(PBdatapt *d) ;
static void constraintTOeqs() ;
static int checkNegativeCS(PBdatapt *d) ;
static int fixNegativeCS(int k, PBdatapt *d) ;
static double poisson(int idum, double d1, double d2, double *distval,
		      double min, double max) ;
static double uniform(int idum, double d1, double d2, double *distval,
		      double min, double max) ;
static double gaussian(int idum, double ave, double sig, double *distval,
		       double d1, double d2) ;

static int monitorCorrect(PBdatapt *d) ;
static int polmonitorCorrect(PBdatapt *d) ;
static int SIMmonitorCorrect(PBdatapt *d) ;
static int SIMpolmonitorCorrect(PBdatapt *d) ;

/* monitor correction functions */
static double monFuncBT4PG(double E) ;
static double monFuncPG2cmfilter(double E) ;
static double (*monCorFunc[])()= {
  monFuncBT4PG, monFuncPG2cmfilter
} ;
static int NmonoCorFunc = 2 ;

static int PBdefineCells(char *filename) ;
static int PBsetflags(PBflags *flgs) ;

/* create a storage poll for CELL info read from file. 24 cells better suff */
static int MAXCELLS = 24 ;
static He3CELLexp expcells[24] ;
static int Ncells = 0 ;  /* actual number of cells read */

/* and a default flags structure */
static PBflags flags = {
  0, 0, 0, 0, 100000, 0, 0, 0,
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
  0, 0, 0, 0, 100000, 0, 0, 0,
  {1, 1, 1, 1},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0},
  {0, 0, 0, 0}
} ;
static PBflags flagsSave ;

static ConstraintEq eqs[4] ;
static int Nconstraint = 0 ;
static ConstraintEq eqsS[4] ;
static int NconstraintS = 0 ;

static int *DBG = &(flags.Debug) ;

typedef struct {
  int n ;
  int i[4] ;
} ilst ;
typedef struct {
  double A[4][4] ;
} matrix ;
typedef struct {
  matrix m[4][4] ;
} matrixmatrix ;


double deter(int n, matrix *m) ;
double invert(int n, matrix *m, matrix *inv) ;
double invertP(int n, matrix *m, matrix *inv, matrixmatrix *p) ;

static int NN = 4 ;

static matrix zeromatrix = {
  {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}}
} ;

static int CalcMethod = 0 ;

/*
  main entrypoint: PBcorrectData
  correct PB data passed in structure PBindata.
  Any of the passed structure pointers can be NULL
  except that if PBindata is not NULL then
  PBoutdata must be also non-NULL or error return (!=0).
  It is of course also an error if any of the double pointers
  in the passed PBindata or PBoutdata are NULL.
  Caller is responsible for storage allocation in
  PBindata and PBoutdata which must have double arrays
  that can hold up to npts.
*/

int PBcorrectData(char *CellFile, PBflags *flgs,
		  int npts, PBindata *in, PBoutdata *out)
{
  /* storage for PB data */
  static PBdatapt d ;

  static char Cstrs[4][4] = {"Cpp","Cmm","Cpm","Cmp"} ;
  static char Sstrs[4][4] = {"Spp","Smm","Spm","Smp"} ;

  int i, ierr, j, k, ic, is ;
  unsigned long secs ;
  double temp, err ;
  FILE *fp ;

  if( CellFile != NULL ) if( PBdefineCells(CellFile) || Ncells < 1 ) return 1 ;

  if( flgs != NULL ) PBsetflags(flgs) ;

  if( in == NULL || npts < 1 ) return 0 ;
  if( out == NULL ) return 3 ;

  constraintTOeqs() ;

  /*
    process each data point
  */
  for( i=0 ; i<npts ; i++ ) {
    /* first just copy the input to the PBdatapt struct d */
    if( in->Cpp == NULL ) {
      if( flags.CountsEnable[0] ) {
	if( *DBG ) printf("NULL Cpp disables that equation\n") ;
	flags.CountsEnable[0] = 0 ;
      }
      d.Y[0] = d.Yr[0] = 0. ; d.Yesq[0] = d.Yresq[0] = 0. ;
    } else {
      if( in->tpp == NULL ) return 3 ;
      secs = in->tpp[i] ;
      d.Y[0] = d.Yr[0] = in->Cpp[i] ;
      if( in->Epp == NULL )
	if( d.Y[0] > 0. ) d.Yesq[0] = d.Yresq[0] = d.Y[0] ;
	else d.Yesq[0] = d.Yresq[0] = 1. ;
      else d.Yesq[0] = d.Yresq[0] = in->Epp[i]*in->Epp[i] ;
    }
    if( in->Cmm == NULL ) {
      if( flags.CountsEnable[1] ) {
	if( *DBG ) printf("NULL Cmm disables that equation\n") ;
	flags.CountsEnable[1] = 0 ;
      }
      d.Y[1] = d.Yr[1] = 0. ; d.Yesq[1] = d.Yresq[1] = 0. ;
    } else {
      if( in->tmm == NULL ) return 3 ;
      secs = in->tmm[i] ;
      d.Y[1] = d.Yr[1] = in->Cpp[i] ;
      if( in->Emm == NULL )
	if( d.Y[1] > 0. ) d.Yesq[1] = d.Yresq[1] = d.Y[1] ;
	else d.Yesq[1] = d.Yresq[1] = 1. ;
      else d.Yesq[1] = d.Yresq[1] = in->Emm[i]*in->Emm[i] ;
    }
    if( in->Cpm == NULL ) {
      if( flags.CountsEnable[2] ) {
	if( *DBG ) printf("NULL Cpm disables that equation\n") ;
	flags.CountsEnable[2] = 0 ;
      }
      d.Y[2] = d.Yr[2] = 0. ; d.Yesq[2] = d.Yresq[2] = 0. ;
    } else {
      if( in->tpm == NULL ) return 3 ;
      secs = in->tpm[i] ;
      d.Y[2] = d.Yr[2] = in->Cpm[i] ;
      if( in->Epm == NULL )
	if( d.Y[2] > 0. ) d.Yesq[2] = d.Yresq[2] = d.Y[2] ;
	else d.Yesq[2] = d.Yresq[2] = 1. ;
      else d.Yesq[2] = d.Yresq[2] = in->Epm[i]*in->Epm[i] ;
    }
    if( in->Cmp == NULL ) {
      if( flags.CountsEnable[3] ) {
	if( *DBG ) printf("NULL Cmp disables that equation\n") ;
	flags.CountsEnable[3] = 0 ;
      }
      d.Y[3] = d.Yr[3] = 0. ; d.Yesq[3] = d.Yresq[3] = 0. ;
    } else {
      if( in->tmp == NULL ) return 3 ;
      secs = in->tmp[i] ;
      d.Y[3] = d.Yr[3] = in->Cmp[i] ;
      if( in->Emp == NULL )
	if( d.Y[3] > 0. ) d.Yesq[3] = d.Yresq[3] = d.Y[3] ;
	else d.Yesq[3] = d.Yresq[3] = 1. ;
      else d.Yesq[3] = d.Yresq[3] = in->Emp[i]*in->Emp[i] ;
    }

    if( in->Cpp == NULL ) d.sec[0] = secs ; else d.sec[0] = in->tpp[i] ;
    if( in->Cmm == NULL ) d.sec[1] = secs ; else d.sec[1] = in->tmm[i] ;
    if( in->Cpm == NULL ) d.sec[2] = secs ; else d.sec[2] = in->tpm[i] ;
    if( in->Cmp == NULL ) d.sec[3] = secs ; else d.sec[3] = in->tmp[i] ;


    temp = in->Ei[i] ;
    if( temp <= 0. ) return 4 ;
    temp /= Dn ;
    d.lambI = TWOPI/sqrt(temp) ;
    temp = in->Ef[i] ;
    if( temp <= 0. ) return 4 ;
    temp /= Dn ;
    d.lambF = TWOPI/sqrt(temp) ;


    if( (ierr = PBcorrectDatapt(&d)) ) {
      if( *DBG ) printf("ERROR in PBcorrectDatapt = %d\n",ierr) ;
    }

    if( *DBG > 1 ) {
      /* print the equations and soln for this datapt */
      if( d.Nfree == d.Nactive ) {
	for( j=0 ; j<d.Nactive ; j++ ) {
	  ic = d.activeEq[j] ;
	  err = 0. ;
	  if( d.Yesq[ic] > 0. ) err = sqrt(d.Yesq[ic]) ;
	  printf("%3s(%10g %10g) = ",Cstrs[ic],d.Y[ic],err);
	  for( k=0 ; k<d.Nfree ; k++ ) {
	    is = d.freeS[k] ;
	    err = 0. ;
	    if( d.Cesq[ic][is] > 0. ) err = sqrt(d.Cesq[ic][is]) ;
	    if( k > 0 ) printf(" + ") ;
	    printf("(%10g %10g)*%3s",d.C[ic][is],err,Sstrs[is]) ;
	  }
	  is = d.freeS[j] ;
	  err = 0. ;
	  if( d.Sesq[is] > 0. ) err = sqrt(d.Sesq[is]) ;
	  printf("   %s(%10g %10g)\n",Sstrs[is],d.S[is],err) ;
	}
      } else {
	printf("Nfree NOT equal Nactive\n") ;
      }
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
    out->R[i] = d.R ;
  }
  return 0 ;
}

static int PBcorrectDatapt(PBdatapt *d)
{
  int j, k, ierr ;

  
  /* intial all 4 equations active and all S(cross-sections) free unknowns */
  d->Nactive = 4 ;
  d->Nfree = 4 ;
  for( j=0 ; j<4 ; j++ ) {
    d->activeEq[j] = j ;
    d->freeS[j] = j ;
  }
  
  /* calc the correction coeficients. this is also where monitor correct */
  if( (ierr = PBcoef(d)) > 0 ) {
    if( *DBG ) printf("error in PBcoef = %d\n", ierr) ;
    return 6 ;
  }
  
  /* apply any constraints on the S */
  for( j=0 ; j<Nconstraint ; j++ ) {
    if( applyConstraint(d, eqs+j) > 0 )
      if( *DBG ) printf("constraint %d failed\n",j) ;
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
	  if( combineEqs(d, j, k) ) if( *DBG ) printf("combineEqs failed\n") ;
	}
      }
    }
    if( flags.CountsAdd2[j] ) {
      for( k=j+2 ; k<4 ; k++ ) {
	if( flags.CountsAdd2[k] ) {
	  if( combineEqs(d, j, k) ) if( *DBG ) printf("combineEqs failed\n") ;
	}
      }
    }
  }

  /* delete any equations where counts are disabled */
  for( j=0 ; j<4 ; j++ ) { if( ! flags.CountsEnable[j] ) deleteEq(d, j) ; }
    
  /* check that Nactive equations == Nfree Unknowns */
  if( d->Nactive != d->Nfree ) {
    if( *DBG ) printf("Nactive != Nfree\n") ;
    return 7 ;
  }
  if( d->Nactive < 1 ) {
    if( *DBG ) printf("NO remaining degrees of freedom\n") ;
    return 7 ;
  }
  
  if( *DBG ) {
    printf("after CTRLS activeEqindices and freeSindices are ( ") ;
    for( j=0 ; j<d->Nactive ; j++ ) printf("%1d ",d->activeEq[j]+1) ;
    printf(")  ( ") ;
    for( j=0 ; j<d->Nfree ; j++ ) printf("%1d ",d->freeS[j]+1) ;
    printf(")\n") ;
  }
  
  /* call PBcorrect to solve for the cross-sections by Gaussian elim */
  if( (ierr = PBcorrect(d)) > 0 ) {
    if( ierr == 2 ) {
      if( *DBG ) printf("invalid Nactive equations or Nactive != Nfree\n") ;
      return 8 ;
    }
    if( *DBG ) printf("error in PBcorrect = %d\n", ierr) ;
    return 9 ;
  }
  
  /* apply any constraints to the final result */
  constrainResult(d) ;

  /*
    Suppose some of the results come back with negative CS solutions
    Use a flag to require constraining negative values to zero?
    When we do this we have to throw out an equation. The logical choice
    is to throw out the equation for the counts corresponding to the
    CS value that we are now going to constrain.
    Note that this may only leave one equation in one unknown.
  */

  if( d->Nactive > 1 && flags.NoNegativeCS && checkNegativeCS(d) > 0 ) {
    /* save the global constraint etc info */
    flagsSave = flags ;
    NconstraintS = Nconstraint ;
    for( j=0 ; j<Nconstraint ; j++ ) eqsS[j] = eqs[j] ;

    while(d->Nactive > 1 && flags.NoNegativeCS && (k=checkNegativeCS(d)) > 0 )
      fixNegativeCS(k,d) ;

    flags = flagsSave ;
    Nconstraint = NconstraintS ;
    for( j=0 ; j<Nconstraint ; j++ ) eqs[j] = eqsS[j] ;
  }

  return 0 ;
}

static int checkNegativeCS(PBdatapt *d)
{
  /*
    called only if d->Nactive > 1 so we can add a constraint
    and throw away an equation
  */

  int i ;
  for( i=0 ; i<4 ; i++ ) {
    if( ! flags.Sconstrain[i] && flags.CountsEnable[i] && d->S[i] < 0 )
      return i+1 ;
  }
  return 0 ;
}
static int fixNegativeCS(int k, PBdatapt *d)
{
  int i, j, ierr ;

  i = k - 1 ;

  flags.Sconstrain[i] = 1 ;
  if( i==0 ) {
    for( j=0 ; j<4 ; j++ ) flags.Spp[j] = 0. ;
  } else if( i==1 ) {
    for( j=0 ; j<4 ; j++ ) flags.Smm[j] = 0. ;
  } else if( i==2 ) {
    for( j=0 ; j<4 ; j++ ) flags.Spm[j] = 0. ;
  } else if( i==3 ) {
    for( j=0 ; j<4 ; j++ ) flags.Smp[j] = 0. ;
  }
  flags.CountsEnable[i] = 0 ;

  constraintTOeqs() ;

  if( (ierr = PBcorrectDatapt(d)) ) {
    if( *DBG ) printf("ERROR in PBcorrectDatapt = %d\n",ierr) ;
    return ierr ;
  }

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
    return 11 ;

  cell = &(ex->cell) ;
  pol = &(ex->pol) ;
  res = &(ex->res) ;

  if( fabs(pol->T) <= 0. ) return 12 ;

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

static findcellsFORdatapt(PBdatapt *d, PBcells *c)
{
  int i, j, k, ic ;

  if( d == NULL || c == NULL ) return 1;

  for( i=0 ; i<4 ; i++ ) {

    /*
      find the cells for the 4 CS of datapt based on the datapt timestamp
    */
    c->cells[0][i] = c->cells[1][i] = -1 ;
    if( ! flags.CountsEnable[i] ) continue ;
    for( j=Ncells-1 ; j>=0 ; j-- ) {
      if( expcells[j].pol.startSecs <= d->sec[i] ) {
	/* ic=PorA==0 means found a Polarizer */
	ic = expcells[j].PorA ;
	c->cells[ic][i] = j ;
	/* look for a second cell even if half-polarized */

	for( k=j-1 ; k>= 0 ; k-- ) {
	  if(expcells[k].pol.startSecs <= d->sec[i] && expcells[k].PorA != ic) {
	    c->cells[1-ic][i] = k ;
	    break ;  /* break the k cells loop */
	  }
	}
	break ; /* done so break the j cells loop */
      }
    } /* ends j cells loop */
    if( c->cells[0][i] < 0 || c->cells[1][i] < 0 ) {
      if( *DBG ) printf("Couldnt find cells for CS %d at datapt.\n", i) ;
      return 2 ;
    }
    /* cell decay constants must be > 0 */
    if(expcells[c->cells[0][i]].pol.T <= 0. ||
       expcells[c->cells[1][i]].pol.T <= 0.) {
      if( *DBG ) printf("zero time constant for a cell\n") ;
      return 2 ;
    }

  }   /* ends i CS loop */
  return 0 ;
}

static int PBcoef(PBdatapt *d)
{

  /*
    Feb 2008 RWE
    PBcoef has been changed to go thru the expcells to find the
    cells that were active for the UNIXtime of each datapt CS
    this could be up to four different pairs of cells
    Must find both Pol and Anal although one of these may be NOCELL
    handled by zero nsL for half-polarized beam.

    Also, any data correction for counting against monitor is done here
  */

  /*
    indices in expcells for each datapt CS
    row 1 for P
    row 2 for A
    indicate missing cell with index < 0
  */

  static PBsetup setup ;
  static PBsetup *s = &setup ;
  static PBcells c ;

  int i, j, k, ierr, ic, NRsum ;

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

  double tp, tm, teff, feff, Rsum ;
  double R[4] ;

  if( d == NULL ) return 1;

  if( ierr = findcellsFORdatapt(d, &c) ) return ierr ;

  /*
    first calc the cell He3 polarizations at datapt time with errs
    use s->(P,A).pol.(PHe,hr0,T) along with d->hr
    to calc current He3 polarizations and errs

    calc nsL*(1+-PHe3) = tau+- and its err which is +- independent
    use tau+- to calc correction coefs along with He3exp (expA, expP)

    finally can calculate transmission factors and errors
    and then coefs and errors

    Jan 29 2008 added calculation of average flipping ratio R
    for the times of the 4 CS msrmnts
  */

  /* apply monitor correction if requested */
  if( flags.MonitorCorrect ) monitorCorrect(d) ;
  else if( flags.PolMonitorCorrect ) polmonitorCorrect(d) ;

  /* Now redo loop over CS loading the correct cells into s */

  NRsum = 0 ;
  Rsum = 0. ;
  for( i=0 ; i<4 ; i++ ) {
    if( c.cells[0][i] < 0 || c.cells[1][i] < 0 ) continue ;
    s->P = expcells[c.cells[0][i]] ;
    s->A = expcells[c.cells[1][i]] ;

    /* efficiency stuff */
    etPsq = s->P.eff.teff ;
    etPsq *= etPsq ;
    etPesq = s->P.eff.terr ;
    etPesq *= etPesq ;
    efP = s->P.eff.feff ;
    efPesq = s->P.eff.ferr ;
    efPesq *= efPesq ;

    etAsq = s->A.eff.teff ;
    etAsq *= etAsq ;
    etAesq = s->A.eff.terr ;
    etAesq *= etAesq ;
    efA = s->A.eff.feff ;
    efAesq = s->A.eff.ferr ;
    efAesq *= efAesq ;

    /* for calc R need product of P and A transport efficiencies and avg feff */
    teff = s->P.eff.teff * s->A.eff.teff ;
    feff = (s->P.eff.feff + s->A.eff.feff)/2. ;

    /*
      Now  calculate the coefs and errsq
      Pol  Ana
      index 0  + +    OFF  OFF
      index 1  - -    ON   ON
      index 2  + -    ON   OFF
      index 3  - +    OFF  ON
      but first compute the transport factors, eAalp, ePbet
    */


    /* i is row index so compute alp and bet for i */
    alp = alpmu[i] ;
    alp1sq = (1.- alp)*(1.- alp) ;
    eAFalp = (1.- alp)*s->A.eff.feff - 1. ;
    eAFalpsq = eAFalp*eAFalp ;
    eAalp = s->A.eff.teff*eAFalp ;
    eAalpsq = eAalp*eAalp ;

    bet = betnu[i] ;
    bet1sq = (1.- bet)*(1.- bet) ;
    ePFbet = (1.- bet)*s->P.eff.feff - 1. ;
    ePFbetsq = ePFbet*ePFbet ;
    ePbet = s->P.eff.teff*ePFbet ;
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

    tp = tsA*tsP ;
    tm = taA*taP ;
    Rsum += (tp + teff*tm)/(tp - teff*(2*feff-1)*tm) ;
    NRsum++ ;

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
  d->R = 0. ;
  if( NRsum > 0 ) d->R = Rsum/NRsum ;

  return 0;
}

static double correctionCoef(He3CELLexp *ex, double tau)
{
  expResol *res ;
  if( ex == NULL ) return 1.;
  res = &(ex->res) ;
  return (1. + res->t2*tau*tau - res->t1*tau) ;
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
  /*
    RWE Feb 14 2008
    NOT happy with error results using the Gaussian elimination and propagating
    errors for every mathematical operation.
    This is likely overestimating errors
    and in fact the results compared to exact solution confirm this.
    I now have the exact algebraic solution for 2 equations along with the err
    orig equ for counts: Cu +-sigu = (Auv +- siguv) Sv
    D = A11*A22 - A12*A21
    solution: Su = (Buv/D) Cv
    where Buv = (-1)^(u+v) Avbar,ubar    cofactors of transpose of A
    then errors are
    sigSw^2 = SUMuv (dSw/dAuv)^2 siguv^2 + SUMu (dSw/dCu)^2 sigu^2
    or
    D^2 sigSw^2 = Sw^2 SUMuv Aubarvbar^2 siguv^2
                  + Cb(Cb - 2*Sw*Abw)sigawbar^2
                  + Ca(Ca - 2*Sw*Aaw)sigbwbar^2
		  + Abwbar^2 siga^2 + Aawbar^2 sigb^2
    where indices are a=1 b=2

    We also have the exact solution for 4 equations when the spin-flipper
    efficiency is perfect and at equal times.
    In that case Aab,uv = Aabarb,ubarv = Aabbar,ubarv = Aabarbbar,ubarvbar
    In terms of 2x2 matrices, cn and cf,
    and 2-vectors Cn, Cf, Sn and Sf the orig equations are then
    Cn = (N/2) (cn cf) Sn
    Cf         (cf cn) Sf
    where cn = (c++ c--)  cf = (c+- c-+)
               (c-- c++)       (c-+ c+-)
    
    Cn = (C++)    Cf = (C+-)  and similar for Sn and Sf
         (C--)         (C-+)

    Because the matrices cn and cf are symmetric the inversion solution is
    ( cn -cf) (Cn) = N/2 ( cn^2-cf^2     0     ) (Sn)
    (-cf  cn) (Cf)       (     0     cn^2-cf^2 ) (Sf)
    where cn^2 - cf^2 = cd is also symmetric.
    The final soln is then:
    ( cd^-1 cn   -cd^-1 cf ) (Cn) = N/2 (Sn)
    (-cd^-1 cf    cd^-1 cn ) (Cf)       (Sf)

    cd = (e f)    e = c++^2 + c--^2 - c+-^2 - c-+^2
         (f e)    f = 2c++c-- - 2c+-c-+
    
    cd^-1 = ( e -f)/(e^2 - f^2)
            (-f  e)
    This soln will not work because typically the CS are at different times.

    RWE Feb 18 2008
    get algebraic inverted soln in terms of determinants via cofactors
    and along with this we can get the partials of each inverse matrix element
    wrt orig matrix elems.

    solving A S = C as
    S = INV C   Su = SUMv INVuv Cv
    and sigSu^2 = SUMv SUMmn(dINVuv/dAmn sigAmn)^2 Cv^2
                + SUMv INVuv^2 sigCv^2

  */

  static matrix M, INV ;
  static matrixmatrix P ;

  PBdatapt D;

  int N, m, row, col, i, im, j, jm, itemp, k, l;
  int indx[4];
  double Amax, temp, tempsq, temp1, temp2, sigsq;
  double det, deter, sumsiguv, Pkl, invij ;

  if( d == NULL ) return 1;

  if( d->Nactive < 1 || d->Nactive > 4 || d->Nactive != d->Nfree ) return 2;
  /* also check the activeEq and freeS indices for invalid range */
  for( i=0 ; i<d->Nactive ; i++ ) {
    if( d->activeEq[i] < 0 || d->activeEq[i] > 3 ||
	d->freeS[i] < 0 || d->freeS[i] > 3 ) return 2;
  }

  /* pack the passed datapt into the local datapt D and matix M */
  for( i=0 ; i<d->Nactive ; i++ ) {
    im = d->activeEq[i];
    D.Y[i] = d->Y[im];
    D.Yesq[i] = d->Yesq[im];
    for( j=0 ; j<d->Nfree ; j++ ) {
      jm = d->freeS[j];
      D.C[i][j] = d->C[im][jm];
      D.Cesq[i][j] = d->Cesq[im][jm];
      M.A[i][j] = d->C[im][jm] ;
    }
  }
  N = d->Nactive;

  if( CalcMethod == 0 ) {
    /* invert M returning INV and partials P */
    det = invertP(N, &M, &INV, &P) ;
    if( fabs(det) <= 0. ) return 3 ;
    /*
      Su = SUMv INVuv Cv
      sigSu^2 = SUMv SUMmn(dINVuv/dAmn sigAmn)^2 Cv^2
      + SUMv INVuv^2 sigCv^2
    */
    for( i=0 ; i<N ; i++ ) {
      D.S[i] = 0. ;
      D.Sesq[i] = 0. ;
      for( j=0 ; j<N ; j++ ) {
	invij = INV.A[i][j] ;
	D.S[i] += invij*D.Y[j] ;
	sumsiguv = 0. ;
	for( k=0 ; k<N ; k++ ) {
	  for( l=0 ; l<N ; l++ ) {
	    Pkl = P.m[i][j].A[k][l] ;
	    sumsiguv += Pkl*Pkl*D.Cesq[k][l] ;
	  }
	}
	D.Sesq[i] += sumsiguv*D.Y[j]*D.Y[j] + invij*invij*D.Yesq[j] ;
      }
    }
  } else if( N == 2 && CalcMethod == 2 ) {

    /*
      Apply the exact soln if there are 2 equations
    */

    deter = D.C[0][0]*D.C[1][1] - D.C[0][1]*D.C[1][0] ;
    if( fabs(deter) <= 0. ) return 3 ;
    sumsiguv = 0. ;
    for( i=0 ; i<2 ; i++ ) {
      for( j=0 ; j<2 ; j++ )
	sumsiguv += D.C[1-i][1-j]*D.C[1-i][1-j]*D.Cesq[i][j] ;
    }
    for( i=0 ; i<2 ; i++ ) {
      D.S[i] = D.C[1-i][1-i]*D.Y[i] - D.C[i][1-i]*D.Y[1-i] ;
      D.S[i] /= deter ;
      D.Sesq[i] = D.S[i]*D.S[i]*sumsiguv + 
	D.Y[0]*(D.Y[0]-2*D.S[i]*D.C[0][i])*D.Cesq[1][1-i] +
	D.Y[1]*(D.Y[1]-2*D.S[i]*D.C[1][i])*D.Cesq[0][1-i] +
	D.C[0][1-i]*D.C[0][1-i]*D.Yesq[1] +
	D.C[1][1-i]*D.C[1][1-i]*D.Yesq[0] ;
      D.Sesq[i] /= (deter*deter) ;
    }
    /* move solutions into freeS indices */
    for( row=0 ; row<N ; row++ ) {
      d->S[d->freeS[row]] = D.S[row];
      d->Sesq[d->freeS[row]] = D.Sesq[row];
    }
    return 0 ;
  } else {
    /* Gaussian elimination soln */
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
  }

  /* move solutions into freeS indices */
  /* this works because d->freeS[row] >= row */
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
static double monFuncPG2cmfilter(double E)
{
  /* meant for elastic condition with 2percent lambda/2 */
  return (0.98 + 0.02/2)/0.98 ;
}


static int monitorCorrect(PBdatapt *d)
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

  if( d == NULL ) return 1 ;
  kI = (TWOPI/d->lambI) ;
  EI = Dn*kI*kI ;
  cor = monCorFunc[flags.MonoSelect](EI) ;
  corsq = cor*cor ;
  for( i=0 ; i<4 ; i++ ) {
    d->Y[i] = cor*d->Yr[i] ; d->Yesq[i] = corsq*d->Yresq[i] ;
  }
  return 0 ;
}
static int SIMmonitorCorrect(PBdatapt *d)
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

  if( d == NULL ) return 1 ;
  kI = (TWOPI/d->lambI) ;
  EI = Dn*kI*kI ;
  cor = monCorFunc[flags.MonoSelect](EI) ;
  corsq = cor*cor ;
  for( i=0 ; i<4 ; i++ ) {
    d->Y[i] = d->Yr[i]/cor ; d->Yesq[i] = d->Yresq[i]/corsq ;
  }
  return 0 ;
}


static int polmonitorCorrect(PBdatapt *d)
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

  int i, ierr ;
  static PBcells c ;

  double EI, kI, lamb, lamb2, cor ;
  double a1, a2 ;
  double tp, tm, tpesq, tmesq, t1, t2 ;

  if( d == NULL ) return 1 ;
  if( ierr = findcellsFORdatapt(d, &c) ) return ierr ;

  lamb = d->lambI ;
  lamb2 = lamb/2. ;
  kI = (TWOPI/lamb) ;
  EI = Dn*kI*kI ;
  cor = monCorFunc[flags.MonoSelect](EI) ;
  /* extract the order fractions from the standard corfac */
  a1 = 1./(2*cor - 1) ;
  a2 = 1 - a1 ;
  if( a2 < 0. ) a2 = 0. ;
  /* correct cross-section counts */
  for( i=0 ; i<4 ; i++ ) {
    if( c.cells[0][i] < 0 || c.cells[1][i] < 0 ) continue ;
    if( ierr = transfac(&expcells[c.cells[0][i]],
		   d->sec[i], lamb, &tp, &tm, &tpesq, &tmesq) ) return ierr ;
    t1 = (tp + tm)/2. ;
    if( ierr = transfac(&expcells[c.cells[0][i]],
		   d->sec[i], lamb2, &tp, &tm, &tpesq, &tmesq) ) return ierr ;
    t2 = (tp + tm)/2. ;
    cor = 1. + 0.5*(a2/a1)*(t2/t1) ;
    d->Y[i] = cor*d->Yr[i] ;
    d->Yesq[i] = cor*cor*d->Yresq[i] ;
  }
  return 0 ;
}
static int SIMpolmonitorCorrect(PBdatapt *d)
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

  int i, ierr ;
  double EI, kI, lamb, lamb2, cor ;
  double a1, a2 ;
  double tp, tm, tpesq, tmesq, t1, t2 ;
  static PBcells c ;

  if( d == NULL ) return 1 ;
  if( ierr = findcellsFORdatapt(d, &c) ) return ierr ;

  lamb = d->lambI ;
  lamb2 = lamb/2. ;
  kI = (TWOPI/lamb) ;
  EI = Dn*kI*kI ;
  cor = monCorFunc[flags.MonoSelect](EI) ;
  /* extract the order fractions from the standard corfac */
  a1 = 1./(2*cor - 1) ;
  a2 = 1 - a1 ;
  if( a2 < 0. ) a2 = 0. ;
  /* correct each of the four cross-section counts */
  for( i=0 ; i<4 ; i++ ) {
    if( c.cells[0][i] < 0 || c.cells[1][i] < 0 ) continue ;
    if( ierr = transfac(&expcells[c.cells[0][i]],
		   d->sec[i], lamb, &tp, &tm, &tpesq, &tmesq) ) return ierr ;
    t1 = (tp + tm)/2. ;
    if( ierr =  transfac(&expcells[c.cells[0][i]],
		   d->sec[i], lamb2, &tp, &tm, &tpesq, &tmesq) ) return ierr ;
    t2 = (tp + tm)/2. ;
    cor = 1 + 0.5*(a2/a1)*(t2/t1) ;
    d->Y[i] = d->Yr[i]/cor ;
    d->Yesq[i] = d->Yresq[i]/cor/cor ;
  }
  return 0 ;
}


/*

  Given a PBsetup and a PBdatapt array that contains
  simulated PB data, compute coefs and underlying CS with errs.
  data set is in the PBdatapt array.
  The data are input from a columnar data file where
  and the PBcorrection requires columns for
  EI EF CountsOFF timestampOFFOFF CountsONON timestampONON CountsONOFF ...
  We will assume the time stamps are in UTC (UNIX time) seconds.
  This program then writes the calculated CS and err to stdout.

  Read in PBsetup info from cell def file
  The files contain 1 line per cell used with the following data:
  cellName	PorA	iDate	iTime	iUNIXtime	E(meV)	lambda	iPol
  iPolErr	T(hr)	Terr(hr)	A(cm2)	nsL	nsLerr	trans	tErr
  flip	fErr	tEmpty	tESlope	Pbar	Lcm	Diacm	rCRVcm	volcc	nsL0
  nsL0err	nsLE	nsLEerr	resolName	Hmos'	Vmos'	dspA	Hcols'
  Hcol2'	VcolsDeg	Vcol2Deg	omRad	Hsig^2	Vsig^2	Xsig^2
  curvCor	angCor

  NB that some of these values are recalc in this software. This is just a line
  from the He3logger spreadsheet which keeps track of cells used during an
  experiment.
  Required fields are;
  cellName PorA iDate iTime iUNIXtime E(meV)ORlambda iPol iPolErr T Terr
  nsL nsLerr trans tErr flip fErr tEmpty



  NB startTimeString should be converted to startTimeSecs UTC UNIX time
  under UNIX with cmnd
  date --date="startTimeString" +%s
  or otherwise calculated

  reader calculates nsL corrected and nsLerr putting them into PBsetup

  There may also be a flag indicating that counts were obtained using
  a normalizing monitor place after the polarizer cell.
  This is not very good since such a monitor will also count lambda/2.
  If no monitor is available before the polarizing cell, it is probably best
  to count against time, but then how to account for flux changes with
  changing incident energy. We would need to know the fraction of lambda/2
  as a function of incident energy
*/


/*
  PBdefineCells reads file for the cell definitions used during an experiment
  These are stored in global expcells ordered by install time so that for
  any dataPt the correct cells to use can be found.
*/

static int PBdefineCells(char *filename)
{
  int i, j, nscan ;
  double cellfac ;
  FILE *fp ;
  char *cp ;
  char buf[2048] ;
  char PorA[4] ;
  char date[32] ;
  char time[32] ;
  char resName[64] ;

  double curvCor, Energy, lambda, beamArea, Pbar, volcc, nsL0, nsL0err ;
  double nsLE, nsLEerr, Hmos, Vmos, dspA, Hcols, Hcol2, VcolsDeg, Vcol2Deg ;
  double omRad, hsigsq, vsigsq ;

  He3CELL *cell ;
  He3CELLpol *pol ;
  expResol *exper ;
  efficiency *eff ;

  He3CELLexp swap ;

  if( (fp = fopen(filename, "r")) == NULL ) {
    if( *DBG ) printf("failed to open %s\n", filename) ;
    return 1 ;
  }

  /* init the number of cells defined to zero */
  Ncells = 0 ;


   while( (fgets(buf, 2047, fp)) != NULL ) {
     if( buf[0] == 'c' && buf[1] == 'e' && buf[2] == 'l' ) continue ;

     cell = &(expcells[Ncells].cell) ;
     pol = &(expcells[Ncells].pol) ;
     exper = &(expcells[Ncells].res) ;
     eff = &(expcells[Ncells].eff) ;

     curvCor = 1. ;
     exper->angcor = 1. ;
     exper->t1 = exper->t2 = 0. ;
     exper->hsigsq = exper->vsigsq = exper->xsigsq = 0. ;
     cell->tEmptySlope = 0. ;

     /*
       sscanf is choking on a tab delimiter after the %ul read of startSecs
       I think, so change all tabs to space
       NOT SO maybe its the ul format should be lu
     */
     /* while( (cp = strchr(buf, '\t')) ) *cp = ' ' ; */

     nscan = sscanf(buf,"%s %s %s %s %lu %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		    cell->name,PorA,date,time,
		    &(pol->startSecs),
		    &Energy,&lambda,&(pol->PHe),&(pol->PHeErr),
		    &(pol->T),&(pol->Terr),
		    &beamArea,
		    &(cell->nsL),&(cell->nsLerr),
		    &(eff->teff),&(eff->terr),&(eff->feff),&(eff->ferr),
		    &(cell->tEmpty),&(cell->tEmptySlope),
		    &Pbar,&(cell->L),&(cell->D),&(cell->R),
		    &volcc,&nsL0,&nsL0err,&nsLE,&nsLEerr,
		    resName,&Hmos,&Vmos,&dspA,&Hcols,&Hcol2,&VcolsDeg,&Vcol2Deg,
		    &omRad,
		    &(exper->hsigsq),&(exper->vsigsq),&(exper->xsigsq),
		    &curvCor,&(exper->angcor)) ;


     if( nscan < 19 && *DBG )
       printf("failed read complete cell data %d from %s\n", Ncells,filename) ;


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
       the lu format will strip that fraction.
     */

     /*
       nsL already has been corrected for cell curvature in spreadsheet
       but need to compute t1 and t2 coefs for correction coef
     */

     if( PorA[0] == 'P' ) expcells[Ncells].PorA = 0 ;
     else expcells[Ncells].PorA = 1 ;

     exper->t1 = -0.5*exper->angcor*(hsigsq + vsigsq) ;
     exper->t2 =  0.5*exper->xsigsq ;

     Ncells++ ;
     if( Ncells >= MAXCELLS && *DBG )
       printf("read maximum number of experiment cells = %d\n", MAXCELLS) ;

   }
   fclose(fp) ;

   /*
     order the cells by install UNIXtime
   */

   for( i=0 ; i<Ncells-1 ; i++ ) {
     for( j=i ; j<Ncells ; j++ ) {
       if( expcells[j].pol.startSecs < expcells[i].pol.startSecs ) {
	 swap = expcells[j] ;
	 expcells[j] = expcells[i] ;
	 expcells[i] = swap ;
       }
     }
   }

   return 0 ;
 }




 static int PBsetflags(PBflags *flgs)
 {
   /* set the global flag structure from the called one */
   if( flgs->MonoSelect >= NmonoCorFunc ) flgs->MonoSelect = 0 ;
   flags = *flgs ;
   return 0 ;
 }


 /*
   PBsim assumes all cells data have already been read into
   global by using PBcorrectData(cellFile, ...)

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

   unsigned long ulhrs, ulsecs ;
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
       if( sscanf(buf,"%s %lu %lf", labl, &ulsecs, &HrStep) < 3 ) {
	 printf("failed to read startSecs HrStep from %s\n", filename);
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
       d.sec[i] = ulsecs + ulhrs ;
     }

     /* here just use PBcoef to get the equation coefs, any mon cor irrelevant */
     if( (ierr = PBcoef(&d)) > 0 ) {
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
       d.Yr[j] = d.Y[j] ;
       d.Yresq[j] = d.Yesq[j] ;

       /* do the reverse monitor correction */
       if( flags.MonitorCorrect ) SIMmonitorCorrect(&d) ;
       else if( flags.PolMonitorCorrect ) SIMpolmonitorCorrect(&d) ;
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
       d.sec[j] = ulsecs + ulhrs ;
     }

     if( (ierr = PBcoef(&d)) > 0 ) {
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
       d.Yr[j] = d.Y[j] ;
       d.Yresq[j] = d.Yesq[j] ;

       if( flags.MonitorCorrect ) SIMmonitorCorrect(&d) ;
       else if( flags.PolMonitorCorrect ) SIMpolmonitorCorrect(&d) ;
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
     I = reset to default flags
     C index  dpp dmm dpm dmp  Sindex = 0 + dpp*Spp + dmm*Smm + dpm*Spm + dmp*Smp
	NB the Sxx corresponding to Sindex is assumed to have dxx=0
	constrain a cross-section
	to disable the constraint on Sindex use
	default is no constraints
	NB indices are 1-4 corresponding to pp mm pm mp
     C -index
     A Cindex1 Cindex2    add two Counts equations
     B Cindex1 Cindex2    add two other Counts equations
     E Cindex1 Cindex2 ...  delete Count equations     deflt all equations enable
     E -Cindex1 ...  to undelete Count equations
     M intflag   flag for monitorCorrection                deflt OFF
     P intflag   flag for monitorafterPolarizerCorrection  deflt OFF
     S intflag   flag for selecting a monochromatorFunction used in monCorrect
     D debug flag
     F int       set simFlux        deflt 100000
     R int       flag sim deviates  deflt OFF
     X int       calcMehtod         deflt 0 determ   1 Gauss 2=2D exact
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

     } else if( buf[0] == 'E' || buf[0] == 'e' ) {
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
     } else if( buf[0] == 'S' || buf[0] == 's' ) {
       if( (nr = sscanf(buf+1,"%d", ic)) < 1 ) continue;	
       if( ic[0] < NmonoCorFunc && ic[0] >= 0 ) flags.MonoSelect = ic[0] ;
     } else if( buf[0] == 'D' || buf[0] == 'd' ) {
       if( (nr = sscanf(buf+1,"%d", ic)) < 1 ) continue;	
       if( ic[0] ) flags.Debug = abs(ic[0]) ;
       else flags.Debug = 0 ;
     } else if( buf[0] == 'I' || buf[0] == 'i' ) {
       if( (nr = sscanf(buf+1,"%d", ic)) < 1 ) continue;	
       if( ic[0] ) flags = defltflags ;
     } else if( buf[0] == 'F' || buf[0] == 'f' ) {
       if( (nr = sscanf(buf+1,"%d", ic)) < 1 ) continue;	
       if( ic[0] > 0 ) flags.SimFlux = ic[0] ;
     } else if( buf[0] == 'R' || buf[0] == 'r' ) {
       if( (nr = sscanf(buf+1,"%d", ic)) < 1 ) continue;	
       if( ic[0] ) flags.SimDeviate = 1 ;
       else flags.SimDeviate = 0 ;
     } else if( buf[0] == 'X' || buf[0] == 'x' ) {
       if( (nr = sscanf(buf+1,"%d", ic)) < 1 ) continue;
       CalcMethod = ic[0] ;
     }
   }

   fclose(fp) ;
   return 0 ;
 }

 int PBreaddata(char *filename)
 {

   /*
     functionally similar to PBcorrectData but reads from file
     instead of passed data structures
     and reads control flags etc from file.CTRL

     read lines
     Ei Ef Cpp tpp Cmm tmm Cpm tmm Cpm tpm Cmp tmp CppErr CmmErr CpmErr CmpErr
     from file and
     automatically disable any missing counts and constrain the
     corresponding CS to zero.
     Then call PBcorrectDatapt
   */

   PBdatapt d ;

   static char Cstrs[4][4] = {"Cpp","Cmm","Cpm","Cmp"} ;
   static char Sstrs[4][4] = {"Spp","Smm","Spm","Smp"} ;

   int i, ierr, j, k, n, is, ic, nw ;
   char *cp ;
   char buf[512], labl[64] ;
   char outfile[512] ;
   double S[4] ;
   double dbl ;
   double iout, fout, err ;
   double hrsecs, HrStep ;
   double EI, EF, EIstart, EIstep, EFstart, EFstep ;

   unsigned long ulhrs ;
   FILE *fp, *fpout, *fpcorrect ;


   if( filename == NULL ) return 1 ;

   flags = defltflags ;

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
   strcpy(outfile, filename) ;
   strcat(outfile, ".out") ;
   if( (fpout = fopen(outfile, "w")) == NULL ) {
     printf("failed to open %s\n", outfile) ;
     return 2 ;
   }
   fprintf(fpout,"# ") ;
   for( j=0 ; j<4 ; j++ ) {
     if( flags.Sconstrain[j] ) continue ;
     fprintf(fpout, "%6s       %3sERR    ", Sstrs[j],Sstrs[j]) ;
   }
   fprintf(fpout,"     NSFflipRatio\n") ;

   n = 0 ;

   while( fgets(buf, 511, fp) ) {
     /* find up to 14 words in buf and try to read each as double */
     cp = buf ;
     if( *cp == '#' ) continue ;
     nw = 0 ;

     while(nw < 14 && *cp != '\n' && *cp != '\0') {
       /* skip white space */
       while( *cp == ' ' || *cp == '\t' ) cp++ ;
       /* break if end of line */
       if( *cp == '\n' || *cp == '\0' ) break ;
       /* try to scan this word as a double */


       if( sscanf(cp, "%lf", &dbl)<1 ) {
	 if( nw < 2 ) {
	   printf("failed to read lambdaI or lambdaF at data pt # %d\n",n+1) ;
	   fclose(fp) ; fclose(fpout) ;
	   return 3 ;
	 } else if( nw < 10 ) {
	   if( flags.CountsEnable[nw/2 - 1] ) {
	     printf("failed to read required counts %d\n",nw/2-1) ;
	     fclose(fp) ; fclose(fpout) ;
	     return 4 ;
	   }
	 } else {
	   d.Yresq[nw-10] = d.Yesq[nw-10] = d.Y[nw-10] ;
	 }
       } else {
	 if( dbl <= 0. ) {
	   printf("read negative value\n") ;
	   fclose(fp) ; fclose(fpout) ;
	   return 4 ;
	 }
	 if( nw < 2 ) {
	   dbl /= Dn ;
	   dbl = TWOPI/sqrt(dbl) ;
	 }
	 if( nw == 0 ) {
	   EI = dbl ;
	   d.lambI = dbl ;
	 } else if( nw == 1 ) {
	   EF = dbl ;
	   d.lambF = dbl ;
	 } else if( nw < 10 && nw%2 == 0 ) {
	   d.Y[nw/2-1] = d.Yr[nw/2-1] = dbl ;
	   d.Yesq[nw/2-1] = d.Yresq[nw/2-1] = dbl ;
	 } else if( nw < 10 && sscanf(cp, "%lu",&(d.sec[nw/2-1])) < 1 ) {
	   printf("failed to read time stamp as unsigned long pt # %d\n",n+1);
	   fclose(fp) ; fclose(fpout) ;
	   return 3 ;
	 } else if( nw >= 10 && nw < 14 ) {
	   d.Yesq[nw-10] = d.Yresq[nw-10] = dbl ;
	 }
       }
       /* go to just past this word */
       while( *cp != ' ' && *cp != '\t' && *cp != '\n' && *cp != '\0' ) cp++ ;
       nw++ ;
       if( *cp == '\n' || *cp == '\0' ) break ;
     }


     /*
       if( flags.MonitorCorrect ) {
       monitorCorrect(&d) ;
       } else if ( flags.PolMonitorCorrect ) {
       polmonitorCorrect(&d);
       }

       monitor corrections now in PBcoef
     */

     if( (ierr = PBcorrectDatapt(&d)) ) {
       if( *DBG ) printf("ERROR in PBcorrectDatapt = %d\n",ierr) ;
     }

     /* write the output */

     for( j=0 ; j<d.Nfree ; j++ ) {
       is = d.freeS[j] ;
       err = 0. ;
       if( d.Sesq[is] > 0. ) err = sqrt(d.Sesq[is]) ;
       fprintf(fpout,"%10g %10g  ",d.S[is],err) ;
     }
     fprintf(fpout,"    %10g\n",d.R) ;

     if( *DBG > 1 ) {
       /* print the equations and soln for this datapt */
       if( d.Nfree == d.Nactive ) {
	 for( j=0 ; j<d.Nactive ; j++ ) {
	   ic = d.activeEq[j] ;
	   err = 0. ;
	   if( d.Yesq[ic] > 0. ) err = sqrt(d.Yesq[ic]) ;
	   printf("%3s(%10g %10g) = ",Cstrs[ic],d.Y[ic],err);
	   for( k=0 ; k<d.Nfree ; k++ ) {
	     is = d.freeS[k] ;
	     err = 0. ;
	     if( d.Cesq[ic][is] > 0. ) err = sqrt(d.Cesq[ic][is]) ;
	     if( k > 0 ) printf(" + ") ;
	     printf("(%10g %10g)*%3s",d.C[ic][is],err,Sstrs[is]) ;
	   }
	   is = d.freeS[j] ;
	   err = 0. ;
	   if( d.Sesq[is] > 0. ) err = sqrt(d.Sesq[is]) ;
	   printf("   %s(%10g %10g)\n",Sstrs[is],d.S[is],err) ;
	 }
       } else {
	 printf("Nfree NOT equal Nactive\n") ;
       }
     }

     n++ ;
   }

   if( *DBG ) printf("%d data points read and corrected\n",n) ;
   fclose(fp) ;


   return 0 ;
 }

/*
  determinant module
  for algebraic linear equ solutions
  allowing error calculation
  start with up to NN = 4 matrices
*/


static void notinLst( ilst *l, ilst *n )
{
  int i, j, found ;
  n->n = 0 ;
  for( i=0 ; i<NN ; i++ ) {
    found = 0 ;
    for( j=0 ; j<l->n ; j++ ) {
      if( i == l->i[j] ) { found = 1 ; break ; }
    }
    if( ! found ) { n->i[n->n] = i ; ++n->n ; }
  }
}

static void appendLst( ilst *l, ilst *n, int i )
{
  *n = *l ;
  if( n->n == NN ) return ;
  n->i[n->n] = i ;
  ++n->n ;
}

/* l not modified */

static double cofactor(matrix *m, ilst *r, ilst *c)
{
  /*
    recursive cofactor calculator
    each instance of cofactor on call stack will have its own
    local versions of r c = rl cl and complements rn cn
    cofactor calcs the determinant of m omitting rows and cols r, c.
  */
  ilst rl, cl ;
  ilst rn, cn ; /* indices not in rl and cl */
  int i ;
  double cof, ss ;

  rl = *r ; cl = *c ;
  notinLst(&rl, &rn) ;
  notinLst(&cl, &cn) ;

  if( rn.n == 1 ) return (m->A[rn.i[0]][cn.i[0]]) ;
  if( rn.n == 2 ) return (m->A[rn.i[0]][cn.i[0]]*m->A[rn.i[1]][cn.i[1]] -
		    m->A[rn.i[0]][cn.i[1]]*m->A[rn.i[1]][cn.i[0]]) ;
  
  /* always expand on the lowest remaining col so that sign starts at 1 */
  appendLst(c, &cl, cn.i[0]) ;
  cof = 0. ;
  ss = 1. ;
  for( i=0 ; i<rn.n ; i++ ) {
    appendLst(r, &rl, rn.i[i]) ;
    /* call cofactor recursively */
    cof += ss*m->A[rn.i[i]][cn.i[0]]*cofactor(m, &rl, &cl) ;
    ss *= -1.;
  }
  return cof ;
}

static double cofactorP(matrix *m, ilst *r, ilst *c, matrix *p)
{
  /*
    recursive cofactor calculator
    each instance of cofactor on call stack will have its own
    local versions of r c = rl cl and complements rn cn
    cofactor calcs the determinant of m omitting rows and cols r, c.
    cofactorP also calculates the partial derivative of the cofactor
    wrt each matrix coeficient and returns those partials in p.
    Now each instance of cofactor P on the call stack will have its own
    local version of a partial derivative matric pl.
    It uses this local version to get a return partials matrix from
    a recursive call to cofactorP.
  */
  ilst rl, cl ;
  ilst rn, cn ; /* indices not in rl and cl */
  int i, j, k ;
  double cof, lcof, ss, Ai0 ;
  matrix pl ;

  rl = *r ; cl = *c ;
  notinLst(&rl, &rn) ;
  notinLst(&cl, &cn) ;

  /* init the return partials matrix to zeros */
  *p = zeromatrix ;

  if( rn.n == 1 ) {
    p->A[rn.i[0]][cn.i[0]] = 1. ;
    return (m->A[rn.i[0]][cn.i[0]]) ;
  }
  if( rn.n == 2 ) {
    p->A[rn.i[0]][cn.i[0]] = m->A[rn.i[1]][cn.i[1]] ;
    p->A[rn.i[1]][cn.i[1]] = m->A[rn.i[0]][cn.i[0]] ;
    p->A[rn.i[1]][cn.i[0]] = -m->A[rn.i[0]][cn.i[1]] ;
    p->A[rn.i[0]][cn.i[1]] = -m->A[rn.i[1]][cn.i[0]] ;
    return (m->A[rn.i[0]][cn.i[0]]*m->A[rn.i[1]][cn.i[1]] -
	    m->A[rn.i[0]][cn.i[1]]*m->A[rn.i[1]][cn.i[0]]) ;
  }
  /* always expand on the lowest remaining col so that sign starts at 1 */
  appendLst(c, &cl, cn.i[0]) ;
  cof = 0. ;
  ss = 1. ;
  for( i=0 ; i<rn.n ; i++ ) {
    appendLst(r, &rl, rn.i[i]) ;
    /* call cofactor recursively */
    lcof = cofactorP(m, &rl, &cl, &pl) ;
    Ai0 = m->A[rn.i[i]][cn.i[0]] ;
    cof += ss*Ai0*lcof ;
    /*
      dcof/dA[rni][cn0] += ss*(lcof + A[rni][cn0]*pl[rni][cn0]) ;
      else dcof/dAij += ss*A[rni][cn0]*pl[i][j] ;
    */
    for( j=0 ; j<4 ; j++ ) {
      for( k=0 ; k<4 ; k++ ) {
	p->A[j][k] += ss*Ai0*pl.A[j][k] ;
      }
    }
    p->A[rn.i[i]][cn.i[0]] += ss*lcof ;

    ss *= -1.;
  }
  return cof ;
}


double deter(int n, matrix *m)
{
  /*
    compute determinant of nxn
  */
  int i, j ;
  ilst r, c ;

  for( i=n, j=0 ; i<NN ; i++, j++ ) { r.i[j] = i ; c.i[j] = i ; }
  r.n = NN - n ;
  c.n = NN - n ;
  return cofactor(m, &r, &c) ;
}

double invert(int n, matrix *m, matrix *inv)
{
  int i, j ;
  double det, ss ;
  ilst r, c ;
  ilst rl, cl ;

  for( i=n, j=0 ; i<NN ; i++, j++ ) { r.i[j] = i ; c.i[j] = i ; }
  r.n = NN - n ;
  c.n = NN - n ;
  det = cofactor(m, &r, &c) ;

  /* now we need the adj of m to get its inverse */
  for( i=0 ; i<n ; i++ ) {
    appendLst(&c, &cl, i) ;
    for( j=0 ; j<n ; j++ ) {
      appendLst(&r, &rl, j) ;
      ss = 1. ;
      if( (i+j)%2 ) ss = -1. ;
      inv->A[i][j] = ss*cofactor(m, &rl, &cl)/det ;
    }
  }
  return det ;
}


double invertP(int n, matrix *m, matrix *inv, matrixmatrix *p)
{
  /*
    invert matrix m returned in inv
    and also calc partials of the inverse matrix elems wrt each m coef
    and return in p
  */

  int i, j, k, l ;
  double det, ss, lcof ;
  ilst r, c ;
  ilst rl, cl ;
  matrix pl, dd ;

  for( i=n, j=0 ; i<NN ; i++, j++ ) { r.i[j] = i ; c.i[j] = i ; }
  r.n = NN - n ;
  c.n = NN - n ;
  det = cofactorP(m, &r, &c, &dd) ;

  for( i=0 ; i<NN ; i++ ) {
    for( j=0 ; j<NN ; j++ ) {
      p->m[i][j] = zeromatrix ;
    }
  }


  /* now we need the adj of m to get its inverse */
  for( i=0 ; i<n ; i++ ) {
    appendLst(&c, &cl, i) ;
    for( j=0 ; j<n ; j++ ) {
      appendLst(&r, &rl, j) ;
      ss = 1. ;
      if( (i+j)%2 ) ss = -1. ;
      lcof = cofactorP(m, &rl, &cl, &pl) ;
      inv->A[i][j] = ss*lcof/det ;
      /*
	dinv-ij/dAmn = ss(det*pl[m][n] - lcof*dd[m][n])/det/det
      */
      for( k=0 ; k<4 ; k++ ) {
	for( l=0 ; l<4 ; l++ ) {
	  p->m[i][j].A[k][l] = ss*(pl.A[k][l]/det - lcof*dd.A[k][l]/det/det) ;
	}
      }

    }
  }
  return det ;
}
