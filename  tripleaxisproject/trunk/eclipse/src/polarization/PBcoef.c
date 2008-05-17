/*
  Calc transmission coefs and errs for a PBdatapt using the datapt time and
  all the info about the cells and polarized beam transport.
  At this level the inputs are:

  C,+-M    correction coefs+- for each cell transmission
  tE,M     empty cell transmission each cell
  nsL,M    corrected numberdensity-sigma-gasLength for each cell with errors

  PHe3,M   initial polarization each cell with errors
  t0,M     initial polarization time hrs
  T,M      decay time for each cell with errors   this is holding fld dep

  et,M     transport efficiency front and back with errors
  eF,M     flipper efficiency front and back with errors
*/

#include <stdlib.h>
#include <math.h>
#include "PB.h"

static double correctionCoef(He3CELLexp *ex, double tau) ;

static int transfac(He3CELLexp *ex, double secs, double lambda,
		     double *ts, double *ta, double *tesq)
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
  double tp, tpesq, tm, tmesq ;


  if( ex == NULL | ts == NULL || ta == NULL || tesq == NULL ) return 1 ;

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
  tp = C*tEmpty*exp(-tau) ;
  tpesq = tp*tp*tauesq ;

  tau = nsL*(1+PHe3) ;
  C = correctionCoef(ex, tau) ;
  tm = C*tEmpty*exp(-tau) ;
  tmesq = tm*tm*tauesq ;

  *ts = tp + tm ;  /* symmetric transmission coef */
  *ta = tp - tm ;  /* antisymmetric transmission coef */
  *tesq = tpesq + tmesq ;
  return 0 ;
}

int PBcoef(PBdatapt *d, PBsetup *s)
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
    index 0  + +
    index 1  - -
    index 2  + -
    index 3  - +
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
    */

    if( (ierr = transfac(&(s->A), d->sec[i], d->lambF, &tsA, &taA, &tAesq)) )
      return ierr ;
    if( (ierr = transfac(&(s->P), d->sec[i], d->lambI, &tsP, &taP, &tPesq)) )
      return ierr ;
 
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
    1. + tau*tau*sw*sw/2. - tau*(1 - L/R/2. + r*r/(R*R)/2.)*(sV*sV + sH*sH)/2. ;
}
