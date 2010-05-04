/*
  Correct a PolarizedBeamDatapoint.
  CountRates have been measured for some set of the standard polarize beam
  cross-sections, ++ -- +- -+, and the correction coeficients and their errors
  have been determined (time dependent).
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

#include <stdlib.h>

#include "PB.h"

int PBcorrect(PBdatapt *d)
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

int applyConstraint(PBdatapt *d, ConstraintEq *eq)
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

int combineEqs(PBdatapt *d, int eq1, int eq2)
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

int deleteEq(PBdatapt *d, int eq)
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

int constrainResult(PBdatapt *d, int nc, ConstraintEq *eqs)
{
  int i, j ;
  double Cf ;
  double *S, *Sesq ;
  ConstraintEq *eq ;

  if( d == NULL || eqs == NULL ) return 1;
  if( nc < 1 ) return 0;


  for( i=0 ; i<nc ; i++ ) {
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

