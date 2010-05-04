/*
  lsq.h

  data structures used by nonlinear least squares fitting routine, lsq
  marked with @ if lsq can mem alloc it
  only  cm and mesh arrays can be reallocated by lsq
  since params and data must be OK
*/


typedef struct		/* constraint equation information */
{
  int        on ;       /* flags active constraint equation */
  int 	    inc ;	/* array index of constrained parameter */
  double     c0 ;	/* constant term in cnstrnt equ */
  int 	     nt ;	/* number of free params on rhs of cnstrnt equ */
  double  cn[4] ;	/* coeficients for each term */
  int	 ifr[4] ;	/* indices of free params on rhs */
} Fit_Equation ;

typedef struct
{
  int	     nc ;	/* number of cnstrnt equs using this param  */
  int	   *ics ;	/* param indices using this param */
  double   *cns ;	/* coefs for each param using this param */
} Fit_PCoefs ;

typedef struct		/* parameters data structure */
{
  double     *p ;	/* array of params */
  double  *part ;       /* array of partials */
  double    *pe ;	/* array of param errors */
  double   *min ;       /* array of param minima allowed */
  double   *max ;       /* array of param maxima allowed */
  double   **cm ;	/* @ correlation matrix of free params */
  double    chi ;	/* chi */
  double outlie ;       /* prob of sig3 * prob sig4 outliers */
  double   sequ ;       /* prob of max positive variation sequence */
  double   flam ;       /* current flambda for Marquardt-Levenberg */
  int   npalloc ;       /* params allocated */
  int   cmalloc ;       /* correlation matrix allocated */
  int	   npar ;	/* number of parameters */
  int    npfree ;	/* number of free params, calculated in lsq */
  int      nbad ;
  int	  niter ;	/* number of iterations */
  int	   nceq ;	/* number of constraint equations, calculated in lsq */
  int      nsec ;       /* lsq calc time in secs */
  int	  *ptyp ;	/* 0=fix 1=vary 2=constrain */
  int   *ipfree ;       /* indices of free params */
  int   *badpar ;       /* irrelevant params */
  int     *sign ;       /* array of param sign change flags */
  int     *ulim ;       /* upper limit on/off flags */
  int     *llim ;       /* lower limit on/off flags */
  int	  *iceq ;	/* indices of params constrained */
  Fit_Equation *ceqs ;	/* pointer to cnstrnt equations */
  Fit_PCoefs *cons ;
  int    imodel ;
  char   *model ;       /* model name */
  double (*f)() ;       /* function pointer */
} Fit_LsqParams ;

typedef double *Fit_Vec ;

typedef struct		/* xy data structure */
{
  int	  npt ;		/* number of data points */
  int  nalloc ;         /* points allocated */
  int meshalloc ;       /* mesh points allocated */
  int nxalloc ;         /* x components allocated */
  int    npar ;         /* copy from params for part alloc */
  int	   *i ;		/* integer x index */
  int	   nx ;		/* number of x-components */
  int      ix ;         /* x index used for plotting, peakfit and autofit */
  int       m ;         /* mesh for calc curve = npts between two data */
  int      mf ;         /* mesh flag set to 1 when yc updated */
  int     mpt ;         /* number of points total for mesh curve */
  Fit_Vec  *x ;		/* x vector for each data point */
  double   *y ;		/* y value */
  double   *e ;		/* y errors */
  void    *user ;       /* per point user info address of C-struc array data */
  double  *yc ;		/* @ calculated y values */
  double  *ec ;		/* @ fit values calc errors */
  Fit_Vec *part ;       /* @ per point partials */
  Fit_Vec *xm ;         /* @ x mesh */
  double  *ym ;         /* @ calc y mesh */
} Fit_LsqData ;


typedef struct			/* fitting controls */
{
  int      status ;     /* 0 normal -1 in background  -2 bgerror >0 errcode */
  int   converged ;     /* converged flag */
  int       reset ;     /* flags reset iteration counts and damping */
  int       nofit ;     /* flags no fit just calc */
  int  singlestep ;     /* flags single iteration */
  int    ycpassed ;     /* flags yc passed instead of calculated */
  int    meshcalc ;     /* flags do mesh calc if f != NULL */
  int enforceonly ;
  int     special ;     /* special flag to request call ff(-1) only */
  int   closepipe ;
  int      igroup ;     /* group index when calculating for a param group */
  int      ngroup ;     /* store ngroups */
  int     waitmax ;     /* max secs to wait for pipe ycalc, 0=wait forever */
  int          bg ;     /* flags background job */
  int   timeoutbg ;     /* timeoutbg wait for bg lsq to finish */
  int      exitbg ;     /* flags closing pipe to background shells */
  int       maxit ;     /* max iterations allowed. deflt 20 */
  int     fixsign ;     /* flags sign change checking */
  int     plimits ;     /* 1 to enforce param limits */
  int   equations ;     /* 1 to use equations */
  int     nrefine ;     /* number of lin eq root refine iters def 0 */
  int       cmode ;     /* converge mode 0=param%  1=chi%  deflt 0 */
  int     ierrwgt ;     /* 1 for wgt param errors by chi deflt 0 */
  int      iworst ;     /* index of worst param to vary */
  double     test ;     /* converge frac change criterion def 0.001 */
  double chiworst ;     /* chi cutoff for vary only worst param def 7.*/
  double   epsref ;     /* eps criterion for lin eq refinement def 0. */
  int    retalloc ;
  char returntxt[80] ;	/* return message */
} Fit_LsqCtrl ;

typedef struct		     /* fitting controls params and data */
{
  Fit_LsqCtrl   *controls ;  /* structure containing management and contrl */
  Fit_LsqParams *params ;    /* structure containing param info */
  Fit_LsqData   *data ;	     /* structure containing data x,y,e,yc,npt etc */
} Fit_Lsq ;





