#define DEBUG 0
#include <stdio.h>

/*------------------------------------------------------------------

  lsq.c

  Marquardt-Levenberg style non-linear least-squares fit
  to a non-linear function by linearization of the
  fitting function. A gradient search proportion is included.
  See Bevington 11-4 p.232-240.
  Note that SVD matrix solution does not make much sense for non-linear
  least squares since the Hessian can only be approximated successively
  and therefore any singularities are ill defined.
  This version includes constraint conditioning,
  and checks against irrelevant parameters via zero diagonal curvature.
  All input and returned values are passed via the Fit_Lsq data structure.
  Storage is only restricted by system.
  For max portability this module does no IO, uses no POSIX or other
  system calls, except memory allocation

  normal opeartion is fitting until converged

  1. initial enforcement of constraints

  2. calc yc (check yc passed flag) chi and a, b matrices.  
     If no fitting goto finish

  3. check worst param and convergence.  If converged goto finish

  4. check chi behavior and adjust damping.  If too many iterations goto finish

  5. if irrelevant parameter goto 2.

  6. calc param steps, check sign change and enforce param constraints. 
     Set xsqrprev goto 2.

  7. finish

  singlestep:
  If we want to return after each iteration, at 6. 
  just finish and return instead of goto 2.


  Determine pointer check from passed arg same value

  basic matrix equation stuff:
  
  a = sum( w(i)*gradj(yc(i))*gradk(yc(i)) )
  e = a inverse is the error matrix
  b = sum( w(i)*(y(i)-yc(i))*gradk(yc(i)) )
  dp = b*e solution gives parameter increments dp
  chisqr reduced = sum( w(i)*(y(i)-yc(i))**2 ) / (npt - npf)
  sigmasqr for p(j) = e(j,j) gives delta(chisqr) = 1 under other param minim

-----------------------------------------------------------------------------*/

/*
  forward FUNCTION DECLARATIONS
*/

static double doubleprecision() ;  /* function to get machine precision */

double	**dmatrix( int nr, int nc ) ;  /* matrix memory allocation */
void	free_dmatrix( double **m, int nr ) ;

/* lin equ */
static int luinv( int lu, double **a, int n, int *indx, double **y) ; 
int luleq( int new, double **a, double **alu, int n, int *indx, 
		  double *b, int nit, double eps ) ;

/*---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lsq.h"
/* lsq.h defines the data structure argument for the lsq function */

/*---------------*
 | lsq  FUNCTION |
 *---------------*/

int fitLsq( Fit_Lsq *lsqarg )
{ 
  static Fit_Lsq       *lsq = NULL ;
  static Fit_LsqData   *lsqd = NULL ;
  static Fit_LsqParams *lsqp = NULL ;
  static Fit_LsqCtrl   *lsqc = NULL ;

  /* volatile pointers assigned every time */

  double  (*ff)() ;

  double    *p ;
  double *part ;
  double **ppart ;
  double *pmin ;
  double *pmax ;

  Fit_Equation	*ce = NULL ;
  Fit_Equation	*ces = NULL ;
  Fit_PCoefs	*pcs = NULL ;

  int	*iceq = NULL ;
  int	*ip = NULL ;
  int   *ulim ;
  int   *llim ;

  double	*y ;
  double	*e ;

  double	*pe = NULL ;
  double	*yc ;
  double	*ec ;


  int	npt, npf, npar, nceq ;

  int	i, j, ji, k, l, jj, kk, ijx ;
  int	ierr = 0 ;
  int   ferr = 0 ;
  int	nzer, nnum ;

  int    fitting = 0 ;


  /* stuff to save in case we don't reset */

  static int	npv, npvs, ips, swap ;
  static int    igroup ;
  static int	itot, it, itg, itgtot, itntot, resetgrad ;
  static int	conv ;
  
  static double	wsum, eval, f, xcut, xsqr, dy, ss, chi, xsqrprev ;
  static double	h, efac, cmden, pes, cj, ck, dv, fdp ;
  static double convfave, xsqrat ;
  

  static  double heps = 1.e-4 ;	/* limit on sign change interval subdivision */
  static  double smallestdouble = -1. ;

  /* static local storage pool */

  static double	**a = NULL ;
  static double	*w = NULL ;
  static double	*dp = NULL ;
  static double	*b = NULL ;
  static double	*bm = NULL ;
  static double	**am = NULL ;
  static double	**alu = NULL ;
  static int	*indx = NULL ;
  static int    *rip = NULL ;
  static double	*ps = NULL ;
  static Fit_PCoefs  *pc = NULL ;
  
  static  int npalloc = 0 ;
  static  int nvalloc = 0 ;
  static  int ndalloc = 0 ;


  static char *errtxt[] = {
	"NULL data structure pointer",
	"ZERO number of points",
	"UNDEFINED function pointer AND NOFIT requested",
	"NULL x-array pointer",
	"ZERO parameters OR not enough memory allocated",
	"NULL parameter array pointer",
	"NULL equations structure pointer",
	"NULL parameter control pointers",
	"FAILED correlation matrix allocation",
	"FAILED free param local storage allocation",
	"ZERO degrees of freedom for fitting",
	"NULL parameter error array pointer",
	"NULL y or error array pointers",
	"FAILED local weights array allocation",
	"FAILED ycalc array allocation",
	"FAILED local parameter arrays allocation",
	"FAILED local equation allocation",
	"FAILED mesh arrays allocation"
	} ;

  static char *fatalerr[] = {
	"maximum allowed number of iterations exceded",
	"maximum number of gradient iterations exceded",
	"parameters insist on changing sign",
	"matrix singular !!  Try different parameters or check partial defs",
	"matrix inversion failed",
	"all free parameters became irrelevant"
	} ;


/******************************************************************************/

  /* every call error checking and assignments */

  if( lsqarg == NULL ) { ierr = 1 ; goto reterr ; }

  if( lsqarg != lsq ) {
    lsq = lsqarg ;
    ierr = 0 ;
    ferr = 0 ;
    lsqc = lsq->controls ;    /* CONTROLS */
    lsqp = lsq->params ;  /* PARAMS   */
    lsqd = lsq->data ;    /* DATA     */
    if( lsqc == NULL || lsqp == NULL || lsqd == NULL ) {
      ierr = 1 ; goto reterr ;
    }
  }

  if( (npt = lsqd->npt) <= 0 || lsqd->nalloc < npt ) {
    ierr = 2 ; goto reterr ;
  }
  ff = lsqp->f ;
  igroup = lsqc->igroup ;

  if( ! lsqc->ycpassed ) {
      if( ff == NULL )  { ierr = 3 ; goto reterr ; }
      if( lsqd->x == NULL ) { ierr = 4 ; goto reterr ; }
  }

  if( (npar = lsqp->npar) <= 0 || lsqp->npalloc < npar ) {
    //printf("lsq npar=%d npalloc=%d\n", npar, lsqp->npalloc) ;
    ierr = 5 ; goto reterr ;
  }
  if( (p = lsqp->p) == NULL ) {
    ierr = 6 ; goto reterr ;
  }
  if( (part = lsqp->part) == NULL ) {
    ierr = 6 ; goto reterr ;
  }
  if( lsqc->ycpassed && (ppart = lsqd->part) == NULL ) {
    ierr = 6 ; goto reterr ;
  }

  if( (pmin = lsqp->min) == NULL ) lsqc->plimits = 0 ;
  if( (pmax = lsqp->max) == NULL ) lsqc->plimits = 0 ;

  nceq = 0 ;
  if( lsqc->equations ) {
    if( (ce = lsqp->ceqs) != NULL ) {
      if( (iceq = lsqp->iceq) == NULL ) { ierr = 7 ; goto reterr ; }
      for( i=0 ; i<npar ; i++ ) {
	if( (ce+i)->on ) {
	  *(iceq+nceq) = i ;
	  nceq++ ;
	}
      }
    }
    lsqp->nceq = nceq ;
  }

  ulim = lsqp->ulim ;
  llim = lsqp->llim ;

  /* find out if we are fitting */
  fitting = 1 ;
  if( lsqc->nofit ) fitting = 0 ;
  npf = 0 ;

  y = lsqd->y ;
  e = lsqd->e ;
  if( y == NULL || e == NULL ) { ierr = 13 ; goto reterr ; }

  if( lsqc->special && !lsqc->ycpassed ) {
      ff( -1, igroup, lsqd->user, lsqd->x[0], npar, p, part ) ;
      lsqc->special = 0 ;
      return (0) ;
  }

  if( npt > ndalloc ) {
    if( (w = (double*) realloc(w,   npt*sizeof(double))) == NULL ) {
      ierr = 14 ; goto reterr ;
    }
    else ndalloc = npt ;
  }

  nzer = 0 ;	wsum = 0. ;
  
  for( i=0 ; i<npt ; i++ ) {
    eval = *(e+i) ;	*(w+i) = 0. ;
    if( eval <= 0. ) {
      if( *(y+i) <= 0. )	nzer++ ;
      else			*(w+i) = 1./(*(y+i)) ;
    } else {
      *(w+i) = 1./eval/eval ;
    }
    wsum += *(w+i) ;
  }

  if( nzer > 0 ) {
    nnum = npt - nzer ;
    if( nnum > 0 ) wsum /= nnum ;
    else if( nnum <= 0 || wsum <= 0. ) wsum = 1. ;
    for( i=0 ; i<npt ; i++ ) { if( *(w+i) <= 0. ) *(w+i) = wsum ; }
  }
  
  /* make sure y-calc arrays are non null */
  
  if( lsqd->yc == NULL ) {
    if( lsqc->ycpassed ) { ierr = 15 ; goto reterr ; }
    lsqd->yc = (double *) realloc(lsqd->yc, npt*sizeof(double)) ;
  }
  if( lsqd->ec == NULL )
    lsqd->ec = (double *) realloc(lsqd->ec, npt*sizeof(double)) ;
  yc = lsqd->yc ;
  ec = lsqd->ec ;
  if( yc == NULL || ec == NULL ) { ierr = 15 ; goto reterr ; }
 
  if( lsqp->ptyp == NULL || lsqp->ipfree == NULL ) {
    fitting = 0 ;
  } else {
    ip = lsqp->ipfree ;
    for( i=0 ; i<npar ; i++ ) {
      if( lsqp->ptyp[i] == 1 ) {
	ip[npf] = i ;
	npf++ ;
	for( j=0 ; j<nceq ; j++ ) { if( iceq[j] == i ) { npf-- ; break ; } }
      }
    }
    if( npf < 1 ) fitting = 0 ;
  }

  if( fitting ) {
    if( lsqp->badpar == NULL ) {
      ierr = 8 ; goto reterr ;
    }
    lsqp->npfree = npf ;
 
    if( npf > nvalloc ) {
      free_dmatrix(am, nvalloc) ;
      free_dmatrix(alu, nvalloc) ;
      am   = dmatrix( npf, npf ) ;
      alu  = dmatrix( npf, npf ) ;
      dp   = (double *) realloc( dp, npf*sizeof(double) ) ;
      b    = (double *) realloc( b,  npf*sizeof(double) ) ;
      bm   = (double *) realloc( bm, npf*sizeof(double) ) ;
      indx = (int*)     realloc(indx,npf*sizeof(int) ) ;
      nvalloc = npf ;
    }

    if(am==NULL || alu==NULL || dp==NULL || b==NULL || bm==NULL || indx==NULL)
      { ierr = 10 ; goto reterr ; }

    if( npt - npf < 1 ) { ierr = 11 ; goto reterr ; }

    if( (pe = lsqp->pe) == NULL ) { ierr = 12 ; goto reterr ; }

    if( npf > lsqp->cmalloc ) {
      if( lsqp->cm != NULL ) free_dmatrix( lsqp->cm, lsqp->cmalloc ) ;
      lsqp->cmalloc = 0 ;
      if( (lsqp->cm = dmatrix(npf, npf)) != NULL ) lsqp->cmalloc = npf ;
    }

    if( (a = lsqp->cm) == NULL ) { ierr = 9 ; goto reterr ; }

    if( npar > npalloc ) {
      if( (pc = (Fit_PCoefs*)realloc(pc, npar*sizeof(Fit_PCoefs))) == NULL ) 
	{ ierr = 16 ; goto reterr ; }
      for( i=npalloc ; i<npar ; i++ ) {
	(pc+i)->ics = NULL ; (pc+i)->cns = NULL ;
      }
      ps   = (double *) realloc(ps,   npar*sizeof(double)) ;
      rip  = (int *)    realloc(rip,  npar*sizeof(int)) ;
      if( ps == NULL || rip == NULL ) {
	ierr = 16 ; goto reterr ;
      }
      npalloc = npar ;
    }

    if( nceq > 0 ) {
      for( i=0 ; i<npar ; i++ ) (pc+i)->nc = 0 ;
      for( i=0 ; i<nceq ; i++ ) {
	ces = ce + *(iceq+i) ;
	for( j=0 ; j < ces->nt ; j++ ) {
	  ji = ces->ifr[j] ;  /* ji is the free param index */
	  ++((pc+ji)->nc) ;
	}
      }
      for( i=0 ; i<npar ; i++ ) {
	pcs = pc + i ;
	if( pcs->nc > 0 ) {
	  pcs->ics = (int *)    realloc(pcs->ics, pcs->nc*sizeof(int)) ;
	  pcs->cns = (double *) realloc(pcs->cns, pcs->nc*sizeof(double)) ;
	  if(pcs->ics==NULL || pcs->cns==NULL) { ierr = 17 ; goto reterr ; }
	  pcs->nc = 0 ;
	}
      }
      for( i=0 ; i<nceq ; i++ ) {
	ces = ce + *(iceq+i) ;
	for( j=0 ; j < ces->nt ; j++ ) {
	  ji = ces->ifr[j] ;/* ji is the free param index */
	  pcs = pc + ji ;
	  *((pcs->ics)+(pcs->nc)) = ces->inc ;
	  *((pcs->cns)+(pcs->nc)) = *((ces->cn)+j) ;
	  ++(pcs->nc) ;
	}
      }
    }
  }  /* finishes if fitting section */

  npv = npf ;
  
  if( lsqc->reset ) {
    f   = 1.e-3 ;	/* f is the gradient search proportion */
    convfave = f ;
    
    itot=0 ;
    it=1 ;
    itg=0 ;	/* counts successive gradient iterations */
    itgtot=0 ;	/* counts total gradient iterations */
    itntot=0 ;	/* counts total normal iterations */
    swap = 0 ;
    conv = 1 ;
    
    xcut = (lsqc->chiworst)*(lsqc->chiworst) ;
    
    resetgrad = 1 ;
    lsqp->nbad = 0 ;
    lsqc->converged = 0 ;
    if( smallestdouble < 0. ) smallestdouble = doubleprecision() ;
  } /* finishes reset section */


  /* enforce constraints on param values */      

  for( i=0 ; i<nceq ; i++ ) {
    ces = ce + *(iceq+i) ;
    p[ces->inc] = ces->c0 ;
    for( j=0 ; j < ces->nt ; j++ ) {
      p[ces->inc] += *((ces->cn)+j) * p[*((ces->ifr)+j)] ;
    }
  }
  
  /* enforce param limits if requested */

  if( lsqc->plimits ) {
    if( llim != NULL ) {
      for( i=0 ; i<npar ; i++ )
	if( llim[i] && p[i] < pmin[i] ) p[i] = pmin[i] ;
    }
    if( ulim != NULL ) {
      for( i=0 ; i<npar ; i++ )
	if( ulim[i] && p[i] > pmax[i] ) p[i] = pmax[i] ;
    }
  }

  if( lsqc->enforceonly ) return (0) ;
  
/******************************************************************************
    program structure is 4 labelled sections
	1. (setup)   matrices setup and convergence test
	2. (failed)  convergence failure
	3. (calc)    parameter step solution
	4. (finish)  cleanup and exit
******************************************************************************/


/*	setup and convergence check section	*/

setup:

  if( resetgrad == 1 ) itg = 0 ;

  if( fitting ) {
    for( j=0 ; j<npv ; j++ ) {
      b[j] = 0. ;
      for( k=0 ; k<npv ; k++ ) { a[j][k] = 0. ; }
    }
  }

  xsqr = 0. ;
  for( i=0 ; i<npt ; i++ ) {
    if( ! lsqc->ycpassed ) {
      yc[i] = ff( i+1, igroup, lsqd->user, lsqd->x[i], npar, p, part ) ;
    } else {
      for( j=0 ; j<npar ; j++ ) part[j] = ppart[i][j] ;
    }
    dy = y[i] - yc[i] ;
    xsqr += w[i]*dy*dy ;
    ec[i] = 0. ; /* RWE added Feb 2006 */

    if( fitting ) {
      /** COLLECT DERIVATIVES FOR CONSTRAINED PARAMETERS **/
      if( nceq > 0 ) {
	for( j=0 ; j<npv ; j++ ) {
	  for( k=0 ; k < (pc+ip[j])->nc ; k++ ) 
	    part[ip[j]] +=
	      *(((pc+ip[j])->cns)+k) * part[*(((pc+ip[j])->ics)+k)] ;
	}
      }
      
      /** CALCULATE a and b MATRICES **/

      for( j=0 ; j<npv ; j++ ) {
	b[j] += w[i]*dy*part[ip[j]] ;
	for( k=0 ; k<npv ; k++ ) {
	  a[j][k] += w[i]*part[ip[j]]*part[ip[k]] ;
	}
      }
      ss = 0. ;	ec[i] = 0. ;
      for( j=0 ; j<npv ; j++ ) {
	for( k=0 ; k<npv ; k++ ) {
	  ss += part[ip[j]]*part[ip[k]]*am[k][j] ;
	}
      }
      if( ss > 0. ) ec[i] = sqrt(ss) ;
      else ec[i] = 0. ;
      /* calc errors needs am? a? from prev iter? */
    }   /* end if fitting */
  }	/* end of NPT loop */

  /*
    if( DEBUG )
    {
    fprintf(fpt,"a matrix b in setup before convergence check\n") ;
    for( j=0 ; j<npv ; j++ )
    {
    for( k=0 ; k<npv ; k++ ) fprintf(fpt,"%10g ", a[j][k]) ;
    fprintf(fpt,"%10g\n", b[j]) ;
    }
    }
  */
  
  xsqr /= (npt-npv) ;
  chi = 0. ;
  if( xsqr > 0. ) { chi = sqrt(xsqr) ; }
  lsqp->chi = chi ;
  lsqp->niter = itot ;
  lsqp->flam = f ;

  /* incr itot before exit on nofit in case user is fitting with novary */
  ++itot ;

  if( npv < 1 ) fitting = 0 ;
  if( ! fitting ) {
    if( itot > lsqc->maxit && itot > 2 ) ferr = 1 ;
    goto finish ;  /* NO FIT */
  }  

  if( itot <= 1 ) {
    for(j=0 ; j<npar ; j++ )	ps[j] = p[j] ;
    if( lsqc->iworst>0 && lsqc->chiworst>0. && xsqr>xcut ) {
      /** FAILED xsqr CUTOFF TEST SO VARY ONLY THE WORST PARAMETER **/
      npvs=npv ;	npv=1 ;
      ips=ip[0] ;	ip[0] = lsqc->iworst ;
      swap = 1 ;
    }
  } else if( swap==1 && (xsqr/xsqrprev <= 0.9 || xsqr <= xcut)) {
    /** xsqr NOW ok SO RELOAD FREE PARAMETERS AND CONTINUE **/
    swap = 0 ;
    ip[0] = ips ;
    npv=npvs ;
  }
  
  /* CONVERGENCE TEST if after first iteration */
  
  if(itot > 1) {
    if( lsqc->cmode==0 ) {
      /* p% convergence test */
      for( j=0 ; j<npv ; j++ ) {
	fdp = fabs(dp[j]) ;
	if( fdp <= smallestdouble ) continue ;
	if( fabs(p[ip[j]]) <= lsqc->test ) {
	  if( fdp <= lsqc->test ) continue ;
	  goto failed;
	}
	if( fabs(dp[j]/p[ip[j]]) > lsqc->test) goto failed;
      }
    } else {
      /* xsqr% convergence test */
      if( xsqr > 0. && fabs(xsqrprev-xsqr)/xsqr > lsqc->test ) goto failed ;
    }

    lsqc->converged = 1 ;

    /*
      for singlestep  we can have itot > 1 and not have gone thru
      calc section at all to transfer a,b to am,bm
    */
    if( lsqc->singlestep ) {
      for( j=0 ; j<npv ; j++ ) {
	bm[j] = b[j] ;
	for( k=0 ; k<npv ; k++ ) am[j][k] = a[j][k] ;
      }
    }
    goto finish ;
  }	/* end itot > 1 convergence test */

/*************************** end setup section *******************************/



/*********************** failed convergence section **************************/

failed:

  if( itot > 1 && xsqrprev < xsqr ) {
    /* chi increase so damp */
    ++itgtot ;
    for( j=0 ; j<npar ; j++ ) p[j] = ps[j] ;
    if( swap ) {
      /* unswap since not working properly */
      ++itg ;	f *= 10. ; swap = 0 ; conv = 0 ;
      npv = npvs ; ip[0] = ips ;
      resetgrad = 0 ;
      goto setup ;
    }

    if( itg >= lsqc->maxit ) { ferr = 2 ; goto finish ; }
      
    /* attempt to accelerate convergence by watching xsqr changes */
    ++itg ;
    if( xsqrprev > 0. ) xsqrat = xsqr/xsqrprev ;
    else xsqrat = 1. ;
    if( conv )  f = convfave ;
    else        f *= 5.*xsqrat ;
    conv = 0 ;

    for( j=0 ; j<npv ; j++ ) {
      b[j] = bm[j] ;
      for( k=0 ; k<npv ; k++ ) a[j][k] = am[j][k] ;
    }

    /*
      if( DEBUG ) {
      fprintf(fpt,"a matrix b in chi increased\n") ;
      for( j=0 ; j<npv ; j++ ) {
      for( k=0 ; k<npv ; k++ ) fprintf(fpt,"%10g ", a[j][k]) ;
      fprintf(fpt,"%10g\n", b[j]) ;
      }
      }
    */ 
 
  } else if( itot > 1 && xsqrprev >= xsqr ) {
    /*
      undamp as chi decreases
      often there is a critical value of f which allows continued convergence
      so carefully undamp using weighted average f found during
      converging steps.
    */
    for( j=0 ; j<npar ; j++ ) ps[j] = p[j] ;
    ++itntot ;
    if( xsqr > 0. ) xsqrat = xsqrprev/xsqr ;
    else xsqrat = 1. ;
    convfave = (convfave + xsqrat*f)/(1. + xsqrat) ;
    if( ! conv ) f = (xsqrat*convfave + f)/(1. + xsqrat) ;
    else f = (convfave + xsqrat*(f/10.))/(1. + xsqrat) ;
    conv = 1 ;	++it ;
    if( it > lsqc->maxit ) { ferr = 1 ; goto finish ; }
  }

  /*
    if( DEBUG ) {
    fprintf(fpt,"a matrix b in chi decreased\n") ;
    for( j=0 ; j<npv ; j++ ) {
    for( k=0 ; k<npv ; k++ ) fprintf(fpt,"%10g ", a[j][k]) ;
    fprintf(fpt,"%10g\n", b[j]) ;
    }
    }
  */  
  

/****************** end failed convergence section ***************************/



/******** calc section SOLVE FOR THE PARAMETER STEPS RETURNED IN b ***********/

  //calc:

  /* ready to calculate parameter steps */
  for( j=0 ; j<npv ; j++ ) {
    /* check for singularity */
    if( a[j][j] <= 0. )	{
      /* fix jth free param and continue if possible */
      lsqp->badpar[lsqp->nbad] = *(ip+j) ;
      lsqp->nbad++ ;
      for( k=j ; k<npv-1 ; k++ ) *(ip+k) = *(ip+k+1) ;
      --npv ;
      resetgrad = 1 ; --itot ;
      if( npv <= 0 ) { ferr = 6 ; fitting = 0 ; goto finish ; }
      goto setup ;
    }
  }
  
  for( j=0 ; j<npv ; j++ ) {
    /* save a,b to am,bm and scale diag */
    bm[j] = b[j] ;
    for( k=0 ; k<npv ; k++ ) am[j][k] = a[j][k] ;
    a[j][j] *= (1.+f) ;
  }

  /*
    if( DEBUG ) {
    fprintf(fpt,"a matrix b in calc before luleq\n") ;
    for( j=0 ; j<npv ; j++ ) {
    for( k=0 ; k<npv ; k++ ) fprintf(fpt,"%10g ", a[j][k]) ;
    fprintf(fpt,"%10g\n", b[j]) ;
    }
    }
  */

  /* b will be replaced with param steps */

  if(luleq(0, a, alu, npv, indx, b, lsqc->nrefine, lsqc->epsref) != 0)
    { ferr = 4 ; goto finish ; }	/* singularity */

  /*
    if( DEBUG ) {
    fprintf(fpt,"a matrix b in calc after luleq\n") ;
    for( j=0 ; j<npv ; j++ ) {
    for( k=0 ; k<npv ; k++ ) fprintf(fpt,"%10g ", a[j][k]) ;
    fprintf(fpt,"%10g\n", b[j]) ;
    } 
    }
  */  
  
  for( j=0 ; j<npv ; j++ ) { dp[j] = b[j] ; p[*(ip+j)] += dp[j] ; }
  
  if( lsqc->fixsign ) {
    /* sign change check */
    for( j=0 ; j<npv ; j++ ) {
      k = *(ip+j) ;
      if( lsqp->sign+k ) {
	h = 1. ;
	while( p[k]*ps[k] < 0. ) {
	  h /= 2. ;
	  if( h < heps ) { ferr = 3 ; goto finish ; }
	  /** PARAMETERS INSIST ON CHANGING SIGN SO EXIT **/
	  p[k] = ps[k] + h*dp[j] ;
	}
	dp[j] *= h ;
      }
    }
  }
  
  /* update the constrained parameters */
  for( k=0 ; k<nceq ; k++ ) {
    ces = ce + *(iceq+k) ;
    p[ces->inc] = ces->c0 ;
    for( l=0 ; l < ces->nt ; l++ ) {
      p[ces->inc] += *((ces->cn)+l) * p[*((ces->ifr)+l)] ;
    }
  }

  /* enforce param limits if requested */

  if( lsqc->plimits ) {
      if( llim != NULL ) {
	for( i=0 ; i<npar ; i++ )
	  if( llim[i] && p[i] < pmin[i] ) p[i] = pmin[i] ;
      }
      if( ulim != NULL ) {
	for( i=0 ; i<npar ; i++ )
	  if( ulim[i] && p[i] > pmax[i] ) p[i] = pmax[i] ;
      }
  }

  /* update xsqrprev and return for next iteration */
  if( conv==1 || itot <= 1 ) 	{ xsqrprev = xsqr ; resetgrad = 1 ; }
  else				{ resetgrad = 0 ; }
  if( ! lsqc->singlestep ) goto setup ;

/**************************** end calc section *******************************/


/**************************** exit code section ******************************/

finish:

  if( xsqr > 0. ) lsqp->chi = sqrt(xsqr) ;
  else		lsqp->chi = 0. ;

  /*
    if( DEBUG ) {
    fprintf(fpt,"a matrix b at start of finish\n") ;
    for( j=0 ; j<npv ; j++ ) {
    for( k=0 ; k<npv ; k++ ) fprintf(fpt,"%10g ", a[j][k]) ;
    fprintf(fpt,"%10g\n", b[j]) ;
    }
    fprintf(fpt,"am matrix bm at start of finish\n") ;
    for( j=0 ; j<npv ; j++ ) {
    for( k=0 ; k<npv ; k++ ) fprintf(fpt,"%10g ", am[j][k]) ;
    fprintf(fpt,"%10g\n", bm[j]) ;
    }
    }
  */
  
  /* recalc error matrix from am which has lambda factor set to zero */
  if( fitting ) {
    if( luleq(0, am, alu, npv, indx, bm, lsqc->nrefine, lsqc->epsref) != 0 )
      { ferr = 4 ; goto reterr ; }	/* singularity */

    /*
      if( DEBUG ) {
      fprintf(fpt,"am matrix bm after luleq in finish\n") ;
      for( j=0 ; j<npv ; j++ ) {
      for( k=0 ; k<npv ; k++ ) fprintf(fpt,"%10g ", am[j][k]) ;
      fprintf(fpt,"%10g\n", bm[j]) ;
      }
      fprintf(fpt,"alu matrix after luleq in finish\n") ;
      for( j=0 ; j<npv ; j++ ) {
      for( k=0 ; k<npv ; k++ ) fprintf(fpt,"%10g ", alu[j][k]) ;
      fprintf(fpt,"\n") ;
      }
      }
    */
  
    if(luinv( 1, alu, npv, indx, am ) != 0) {
      ferr = 5 ; goto reterr ;
    }  /* put inverse in am */

    /*
      if( DEBUG ) {
      fprintf(fpt,"am inverse matrix after luinv in finish\n") ;
      for( j=0 ; j<npv ; j++ ) {
      for( k=0 ; k<npv ; k++ ) fprintf(fpt,"%10g ", am[j][k]) ;
      fprintf(fpt,"\n") ;
      }
      }
    */
  
    efac = 1. ;
    if( lsqc->ierrwgt ) efac = chi ; /* weight param errors by chi */

    for( k=0 ; k<npar ; k++ ) pe[k] = 0. ;
    for( j=0 ; j<npv ; j++ ) {
      if( am[j][j] <= 0. )	pe[ip[j]] = 0. ;
      else			pe[ip[j]] = efac*sqrt(am[j][j]) ;
    }
  
    for( i=0 ; i<npv ; i++ ) {
      /* correlation matrix */
      for( j=0 ; j<npv ; j++ ) {
	cmden = sqrt(fabs(am[i][i]))*sqrt(fabs(am[j][j])) ;
	a[i][j] = 0. ;
	if( cmden > 0. ) a[i][j] = am[i][j]/cmden ;
      }
    }
  
    /* also compute final and errors for constrained parameters */
    /* use rip to store ip inverted */
  
    for( i=0 ; i<npar ; i++ ) rip[i] = 0 ;
    for( i=0 ; i<npv  ; i++ ) rip[ip[i]] = i ;
      
    for( i=0 ; i<nceq ; i++ ) {
      ces = ce + *(iceq+i) ;
      p[ces->inc] = ces->c0 ;
      pes = 0. ;
      for( j=0 ; j < ces->nt ; j++ ) {
	cj = *((ces->cn)+j) ;
	jj = *((ces->ifr)+j) ;
	p[ces->inc] += cj * p[jj] ;
	for( k=0 ; k < ces->nt ; k++ ) {
	  ck = *((ces->cn)+k) ;
	  kk = *((ces->ifr)+k) ;
	  pes += cj*ck*pe[jj]*pe[kk]*a[rip[jj]][rip[kk]] ;
	}
      }
      pe[ces->inc] = 0. ;
      if( pes > 0. ) pe[ces->inc] = sqrt(pes) ;
    }  

    /* RWE Feb 2006 compute final yc errors */
    if ( ff != NULL ) {
      for( i=0 ; i<npt ; i++ ) {
	ss = 0. ;
	ec[i] = 0. ;
	yc[i] = ff( i+1, igroup, lsqd->user, lsqd->x[i], npar, p, part ) ;
	for( j=0 ; j<npv ; j++ )
	  for( k=0 ; k<npv ; k++ )
	    ec[i] += part[ip[j]]*part[ip[k]]*pe[ip[j]]*pe[ip[k]] ;
	if( ec[i] > 0. ) ec[i] = sqrt(ec[i]) ;
      }
    }
  }

  // for passed function calculations cant do mesh here
  if( ! lsqc->meshcalc || lsqd->m < 0 || lsqc->ycpassed ) goto reterr ;
  if( ff == NULL ) { ierr = 18 ; goto reterr ; }

  lsqd->mpt = (npt-1)*(lsqd->m+1) + 1 ;
  if( lsqd->mpt > lsqd->meshalloc || lsqd->nx != lsqd->nxalloc ) {
    lsqd->ym = (double*)realloc(lsqd->ym, lsqd->mpt*sizeof(double)) ;
    if( lsqd->xm != NULL ) free_dmatrix( lsqd->xm, lsqd->meshalloc ) ;
    lsqd->xm = dmatrix( lsqd->mpt, lsqd->nx ) ;

    if( lsqd->ym == NULL || lsqd->xm == NULL ) {
      ierr = 18 ; goto reterr ;
    }
    lsqd->meshalloc = lsqd->mpt ;
    lsqd->nxalloc = lsqd->nx ;
  }
  ijx = 0 ;
  dv = lsqd->m + 1. ;

  for( i=0 ; i<npt-1 ; i++ ) {
    for( j=0 ; j<=lsqd->m ; j++ ) {
      for( k=0 ; k < lsqd->nx ; k++ ) {
	lsqd->xm[ijx][k] = (1-j/dv)*lsqd->x[i][k] + (j/dv)*lsqd->x[i+1][k] ;
      }
      ijx++ ;
    }
  }
  for( k=0 ; k < lsqd->nx ; k++ ) lsqd->xm[ijx][k] = lsqd->x[npt-1][k] ;

  ijx = 0 ;
  for( i=0 ; i<npt-1 ; i++ ) {
    for( j=0 ; j<=lsqd->m ; j++ ) {
      lsqd->ym[ijx] = ff(ijx+1, igroup, lsqd->user, lsqd->xm[ijx],
			 npar, p, part);
      ijx++ ;
    }
  }
  lsqd->ym[ijx] = ff(ijx+1, igroup, lsqd->user, lsqd->xm[ijx], npar, p, part) ;

reterr:

  //if( lsqc->returntxt == NULL )
  //lsqc->returntxt = (char*)realloc(lsqc->returntxt, 80*sizeof(char)) ;
  //if( lsqc->returntxt != NULL ) {
  if( ierr > 0 ) { strcpy(lsqc->returntxt,errtxt[ierr-1]) ; }
  else if( ferr > 0 ) { strcpy(lsqc->returntxt,fatalerr[ferr-1]) ; }
  else if( fitting ) { strcpy(lsqc->returntxt, "converged") ; }
  else { strcpy(lsqc->returntxt,"nofit") ; }
  //}
  lsqc->status = 0 ;
  if( ierr > 0 ) lsqc->status = 100+ierr ;
  if( ferr > 0 ) lsqc->status = ferr ;
  return (lsqc->status) ;
}


/*----------------------------------------------------------------------------*/




/*---------------------------------------------------------------------------*/

/* linear equations stuff mostly from Numerical Recipes */


static double doubleprecision()
{
  /* get machines DOUBLE precision */
  double s, s1 ;
  s = 0.00001 ;
  s1 = 1. + s ;
  while( s1 != 1. ) { s /= 2 ; s1 = 1. + s ; }
  return (2*s) ;
}

#define	TINY	1.0e-20 ;

/******************************************************************************
	given nxn matrix a replace it by the LU decomposition of a rowwise perm
	of itself. a and n are input. a is output. indx is output vector
	which records the row permutation effected by the partial pivoting
	d is output as +- 1 depending on even or odd row interchanges.
	Used with lubksb to solve linear equations
******************************************************************************/

static int ludcmp( double **a, int n, int *indx, double *d )
{
  
  int i, j, k, nm1 ;
  int imax = 0 ;
  double big, dum, sum, temp ;
  double *vv ;	/* will store implicit scaling of each row */
  
  vv = (double *) malloc( (unsigned) n*sizeof(double) ) ;
  if( vv == NULL ) return (1) ;
  *d = 1.0 ;
  nm1 = n - 1 ;
  for( i=0 ; i<n ; i++ )
    {
      big = 0.0 ;
      for( j=0 ; j<n ; j++ )
	{ if( (temp=fabs(a[i][j])) > big ) big = temp ; }
      if( big == 0.0 ) return (1) ;
      vv[i] = 1.0/big ;
    }
  
  for( j=0 ; j<n ; j++ )
    {
      for( i=0 ; i<j ; i++ )
	{
	  sum = a[i][j] ;
	  for( k=0 ; k<i ; k++ )	sum -= a[i][k]*a[k][j] ;
	  a[i][j] = sum ;
	}
      big = 0.0 ;
      for( i=j ; i<n ; i++ )
	{
	  sum = a[i][j] ;
	  for( k=0 ; k<j ; k++ )	sum -= a[i][k]*a[k][j] ;
	  a[i][j] = sum ;
	  if( (dum=vv[i]*fabs(sum)) >= big )
	    { big = dum ;	imax = i ; }
	}
      if( j != imax )	/* do we need to interchange rows */
	{	/* yes */
	  for( k=0 ; k<n ; k++ )
	    {
	      dum = a[imax][k] ;
	      a[imax][k] = a[j][k] ;
	      a[j][k] = dum ;
	    }
	  *d = -(*d) ;
	  vv[imax] = vv[j] ;	/* also interchange scl fact */
	}
      indx[j] = imax ;
      if( a[j][j] == 0.0 ) { a[j][j] = TINY ; return (1) ; }
      if( j != nm1 )
	{
	  dum = 1.0/(a[j][j]) ;
	  for( i=j+1 ; i<n ; i++ ) a[i][j] *= dum ;
	}
    }
  free(vv) ;
  return (0) ;
}



/******************************************************************************
   solves the set of n linear equations Ax = B. Here a is input as LU
   decomp of A. b is input as the rhs vector B and return with the
   solution vector X. a, n and indx are not modified.
******************************************************************************/

static void lubksb( double **a, int n, int *indx, double *b )
{
  int i, ii = -1 ;
  int ip, j, nm1 ;
  double sum ;
  
  nm1 = n - 1 ;
  for( i=0 ; i<n ; i++ )
    {
      ip = indx[i] ;
      sum = b[ip] ;	b[ip] = b[i] ;
      if( ii >= 0 )	for( j=ii ; j<=i-1 ; j++ ) sum -= a[i][j]*b[j] ;
      else if( sum != 0. )	ii = i ;
      b[i] = sum ;
    }
  for( i=nm1 ; i>=0 ; i-- )
    {
      sum = b[i] ;
      for( j=i+1 ; j<=n-1 ; j++ ) sum -= a[i][j]*b[j] ;
      b[i] = sum/a[i][i] ;
    }
}



/***********************************************************
 inverse of matrix A with dimension N
 note A is returned as LU decompostion
 call with lu = 1 if a is already LU decomp
***********************************************************/

static int luinv( int lu, double **a, int n, int *indx, double **y)
{
  int	i, j ;
  double	d ;
  double	*col ;
  
  if( (col = (double *) calloc( n, sizeof(double) )) == NULL ) return (1);
  if( ! lu )
    if( ludcmp( a, n, indx, &d ) != 0 ) return (1) ;
  
  for( j=0 ; j<n ; j++ )
    {
      for( i=0 ; i<n ; i++ ) col[i] = 0.0 ;
      col[j] = 1.0 ;
      lubksb( a, n, indx, col ) ;
      for( i=0 ; i<n ; i++ ) y[i][j] = col[i] ;
    }
  free( col ) ;
  return (0) ;
}


/************************************
  iterate on solution x
************************************/

static void mprove( double **a, double **alu, int n, int *indx,
		   double *b, double *x, int nit, double eps )
{
  int i, j, k ;
  double *r ;
  double rmax, rmax2 ;
  
  r = (double *) calloc( n, sizeof(double) ) ;
  
  for( k=0 ; k<nit ; k++ )
    {
      for( i=0 ; i<n ; i++ )
	{
	  r[i] = -b[i] ;
	  for( j=0 ; j<n ; j++ ) r[i] += a[i][j]*x[j] ;
	}
      lubksb( alu, n, indx, r ) ;
      rmax = fabs(r[0]) ;
      for( i=0 ; i<n ; i++ )
	{
	  x[i] -= r[i] ;
	  rmax2 = fabs(r[i]) ;
	  if( rmax2 > rmax ) rmax = rmax2 ;
	}
      if( rmax <= eps ) break ;
    }
  
  free(r) ;
}


/***********************************************************
 linear equ solution of Ax=B with dimension N
 if ! lu must do LU decomposition and put into alu
 a is left unchanged
***********************************************************/

int luleq( int lu, double **a, double **alu,
		 int n, int *indx, double *b, int nit, double eps )
{
  int i, j ;
  double d ;
  double *bb = NULL ;
  
  if( ! lu )
    {
      for( i=0 ; i<n ; i++ )
	for( j=0 ; j<n ; j++ ) alu[i][j] = a[i][j] ;
      if( ludcmp( alu, n, indx, &d ) != 0 ) return (1) ;
    }
  
  if( nit > 0 )
    {
      bb = (double *) calloc( n, sizeof(double) ) ;
      for( i=0 ; i<n ; i++ ) bb[i] = b[i] ;
    }
  
  lubksb( alu, n, indx, b ) ;	/* get solution into b */
  
  if( nit > 0 )
    {
      mprove( a, alu, n, indx, bb, b, nit, eps ) ;
      free( bb ) ;
    }
  return (0) ;
}

double **dmatrix( int nr, int nc )
{
  int i ;
  double **m ;
  
  m = (double **) malloc( (unsigned) nr*sizeof(double*) ) ;
  if( m == NULL ) return (NULL) ;
  for( i=0 ; i<nr ; i++ )
    {
      *(m+i) = (double *) malloc( (unsigned) nc*sizeof(double) ) ;
      if( *(m+i) == NULL ) return (NULL) ;
    }
  return (m) ;
}

void free_dmatrix( double **m, int nr )
{
  int i ;
  if( m == NULL ) return ;
  for( i=nr-1 ; i>=0 ; i-- )	free(*(m+i)) ;
  free(m) ;
}





























