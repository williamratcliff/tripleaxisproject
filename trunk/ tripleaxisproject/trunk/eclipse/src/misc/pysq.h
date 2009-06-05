/*
 * pysq.h --
 *
 *     build dependencies should be in this file.
 *     for example, the build process may prepare a pysqConfig.h file
 *     that does #def for the preprocessor to tell the compiler
 *     where and which include files are available.
 */

#ifndef _PYSQ_H
#define _PYSQ_H




#define	 PI		3.14159265
#define	 PIOVER2	1.5707963
#define  TWOPI		6.2831853
#define  DEGTORAD	0.017453293
#define	 RADTODEG	57.295779513

#define  DNEUT          2.072141789

#define  smallestSQ     0.0001

#define  MAXFIL         128

#define  ZERO           { 0., 0., 0. }
#define  XVEC           { 1., 0., 0. }
#define  YVEC           { 0., 1., 0. }
#define  ZVEC           { 0., 0., 1. }
#define  UNITMATRIX     { XVEC, YVEC, ZVEC }





#ifdef WIN32
#define STRICT
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef STRICT
#undef WIN32_LEAN_AND_MEAN
#include <windowsx.h>
#endif /* WIN32 */

/*
  build may define the directory where pysq auxilliary files are.

  although the plan is to have all auxilliary files for this extension
  in the same location so that if pysq can be imported any of the
  auxilliary files can also be imported.
*/
#define PYSQ_LIBRARY "/usr/local/lib/pysq"

/*
  build should define where the python header file is.
  This can be done via -I(directory) during make
*/
#include <Python.h>

#ifndef EXPORT
#define EXPORT
#endif

#undef EXTERN

#ifdef __cplusplus
#   define EXTERN extern "C" EXPORT
#else
#   define EXTERN extern EXPORT
#endif

#ifndef _ANSI_ARGS_
#   define _ANSI_ARGS_(x)       ()
#endif



/*
  essential includes
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
/* pysq uses strcat strcpy strncmp */



/*
  N.B. the assert library may be replaced by Python memory exception call
  so it is not essential, but used for debugging mem problems
  see below
  For now we are not using assert
*/

//#include <assert.h>

/*
 * The includes are setup so that there should be no
 * build dependencies.
 * If an include is not found then we try to make adjustments below
 * assuming that that include will be commented out
 */

/* do we need errno */
//#include <errno.h>

/* character types, I dont think this is used in pyfit */
//#include <ctype.h>

/* need memory features.h */
//#include <memory.h>

/* do we need unistd */
//#include <unistd.h>

/* may need malloc.h on some systems, normally in stdlib */
//#include <malloc.h>

/* float defines things like DBL_MIN DBL_MAX  but not essential */
#include <float.h>
#include <limits.h>


#ifndef M_PI
#define M_PI    	3.14159265358979323846
#endif /* M_PI */

#ifndef M_PI_2
#define M_PI_2		1.57079632679489661923
#endif

#ifndef M_SQRT2
#define M_SQRT2		1.41421356237309504880
#endif /* M_SQRT2 */

#ifndef M_SQRT1_2
#define M_SQRT1_2	0.70710678118654752440
#endif /* M_SQRT1_2 */

#ifndef SHRT_MAX
#define SHRT_MAX	0x7FFF
#endif /* SHRT_MAX */

#ifndef SHRT_MIN
#define SHRT_MIN	-(SHRT_MAX)
#endif /* SHRT_MAX */

#ifndef USHRT_MAX
#define	USHRT_MAX	0xFFFF
#endif /* USHRT_MAX */

#ifndef INT_MAX
#define INT_MAX		2147483647
#endif /* INT_MAX */

/* suppose float.h is not available, then _FLOAT_H___ will be undef */
#ifndef _FLOAT_H___

/*
 * ----------------------------------------------------------------------
 *
 * DBL_MIN, DBL_MAX --
 *
 * 	DBL_MAX and DBL_MIN are the largest and smaller double
 * 	precision numbers that can be represented by the floating
 * 	point hardware. If the compiler is ANSI, they can be found in
 * 	float.h.  Otherwise, we use HUGE_VAL or HUGE to determine
 * 	them.
 *
 * ----------------------------------------------------------------------
 */
/*
 * Don't want to include __infinity (definition of HUGE_VAL (SC1.x))
 */

#ifdef HUGE_VAL
#define DBL_MAX		HUGE_VAL
#define DBL_MIN		(1/HUGE_VAL)
#else
#ifdef HUGE
#define DBL_MAX		HUGE
#define DBL_MIN		(1/HUGE)
#else
/*
 * Punt: Assume values simple and relatively small
 */
#define DBL_MAX		3.40282347E+38
#define DBL_MIN		1.17549435E-38
#endif /*HUGE*/
#endif /*HUGE_VAL*/
#endif /*! _FLOAT_H___ */

#undef INLINE
#ifdef __GNUC__
#define INLINE inline
#else
#define INLINE
#endif
#undef EXPORT
#define EXPORT

#undef MIN
#define MIN(a,b)	(((a)<(b))?(a):(b))

#undef MAX
#define MAX(a,b)	(((a)>(b))?(a):(b))

/*
 * ----------------------------------------------------------------------
 *
 *  	The following are macros replacing math library functions:
 *  	"fabs", "fmod", "abs", "rint", and "exp10".
 *
 *  	Although many of these routines may be in your math library,
 *  	they may not be in Python and so may not get loaded when
 *      this extension library is loaded. This makes it
 *  	difficult to dynamically load this library as a shared
 *  	object unless the math library is also shared (which isn't
 *  	true on several systems).  We can avoid the problem by
 *  	replacing the "exotic" math routines with macros.
 *
 * ----------------------------------------------------------------------
 */
#undef ABS
#define ABS(x)		(((x)<0)?(-(x)):(x))

#undef EXP10
#define EXP10(x)	(pow(10.0,(x)))

#undef FABS
#define FABS(x) 	(((x)<0.0)?(-(x)):(x))

#undef SIGN
#define SIGN(x)		(((x) < 0.0) ? -1 : 1)

/*
 * Be careful when using the next two macros.  They both assume the floating
 * point number is less than the size of an int.  That means, for example, you
 * can't use these macros with numbers bigger than than 2^31-1.
 */
#undef FMOD
#define FMOD(x,y) 	((x)-(((int)((x)/(y)))*y))

#undef ROUND
#define ROUND(x) 	((int)((x) + (((x)<0.0) ? -0.5 : 0.5)))

#define TRUE 	1
#define FALSE 	0
/*
 * The macro below is used to modify a "char" value (e.g. by casting
 * it to an unsigned character) so that it can be used safely with
 * macros such as isspace.
 */
#define UCHAR(c) ((unsigned char) (c))

#undef VARARGS

#ifdef __cplusplus
#define ANYARGS (...)
#define VARARGS(first)  (first, ...)
#define VARARGS2(first, second)  (first, second, ...)
#else
#define ANYARGS ()
#define VARARGS(first) ()
#define VARARGS2(first, second) ()
#endif /* __cplusplus */


/*
 * Would like to set the Py memory exception and return
 */


#define isNULL(EX)  ( !(EX) ? 1 : 0 )



/* Im not sure we need stdarg or varargs */

#if defined(__STDC__)
#include <stdarg.h>
#else
#include <varargs.h>
#endif

/* COMPLEX declarations */
typedef struct { double r,i ; } DCMPLX ;

static DCMPLX Cadd( DCMPLX a, DCMPLX b ) ;
static DCMPLX Csub( DCMPLX a, DCMPLX b ) ;
static DCMPLX Cmul( DCMPLX a, DCMPLX b ) ;
static DCMPLX Cdiv( DCMPLX a, DCMPLX b ) ;
static DCMPLX Cmplx( double r, double i ) ;
static double Cabs( DCMPLX c ) ;
static double Cabsq( DCMPLX c ) ;
static DCMPLX Cconj( DCMPLX c ) ;
static DCMPLX Csqrt( DCMPLX c ) ;
static DCMPLX RCmul( double f, DCMPLX c ) ;
static DCMPLX Cexp( DCMPLX c ) ;
static DCMPLX RCexp( double f ) ;
static DCMPLX CRpow( DCMPLX c, double f ) ;


/*
  data structures for ATOMdef ATOMSlist and Qlists
*/


typedef struct
{
  double m[3][3] ;
} M33 ;
typedef struct
{
  int n ;
  double j2frac ;
  double coefs[14] ;
} MagFormFactor ;
typedef struct
{
  int iso ;
  double u[3][3] ;
} DebyeWallerFactor ;

/* structure for an atom definition to go in the atomDef hashtable */
typedef struct
{
  int i ;       /* index to keep track of creation order */
  int atomic ;  /* atomic number */
  char *name ;
  char *descript ;
  DCMPLX b ;    /* complex scattering amplitude */
  double Eres, HW ; /* Eres > 0 when there is a resonance */
  int nbcoef ;
  double bcoef[9] ; /* GSAS bRealNonRes and coefs for E dep of b */
  double m ;    /* magnetic moment in Bohr magnetons */
  MagFormFactor *ff ;     /* magnetic form factor */
  DebyeWallerFactor *dw ;
} ATOMdef ;

/* data structure for an atom placed in a material structure */
typedef struct
{
  int omit ;
  DCMPLX b ;
  double r[3] ; /* position */
  double s[3] ; /* moment direction in reduced coords */
  double ss[3] ; /* moment direction in sample Cartesian */
  double m ;    /* override magnetic moment */
  double c ;    /* site occupation */
} myATOM ;
typedef struct
{
  int i ;       /* creation index */
  int omit ;
  int done ;    /* calculation done flag */
  char *name ;  /* from the name field of ATOMdef or some other descriptor */
  ATOMdef *atom ;
  int spcgrp ;
  int setting ;
  int isub ;
  char Wsym ;   /* possible Wyckoff symmetry symbol */
  char *ptsym ; /* point symmetry from Wyckoff tables */
  int mult ;    /* site multiplicity */
  int natoms ;  /* actual number of atoms generated */
  ATOM *atoms ;
} ATOMgroup ;
typedef struct
{
  int n ;
  int nalloc ;
  ATOMgroup *atomslist ;
} ATOMSlist ;
typedef struct
{
  double hkl[3] ;
  double qun[3] ;
  double q ;
  double tt, om ;
  double lor ;
  DCMPLX S0, Sx, Sy, Sz ;
  DCMPLX Sp0, Spx, Spy, Spz ;
  double csp, csm, pp, mm, pm, mp ;
  double pbx[3], pby[3], pbz[3] ;
  double outofplane ;
} Q ;

typedef struct
{
  int n ;
  double hkl[3] ;
  double step[3] ;
  Q *q ;
} Qs ;
typedef struct
{
  int n ;
  int nalloc ;
  Qs *qlist ;
} Qlist ;



typedef struct
{
  int n ;
  int *i ;
} Int_List ;

typedef struct
{
  int n ;
  char **list ;
} String_List ;

typedef struct
{
  int i ;
  char *name ;
} Flag ;

/*
 * spacegroup data structures
 */

typedef struct {
  char *ptstr ;
  double m[3][3] ;
  double t[3] ;
} Subpt ;
typedef struct {
  int index ;
  int npts ;
  char Wyck ;
  char sitesym[12] ;
  Subpt *pts ;
} Subgrp ;
typedef struct {
  char symbol[16] ;
  char origin[16] ;
  int spcnum ;
  int setting[6] ;
  int ntrans ;
  Subpt trans[4] ;
  int nsub ;
  Subgrp *sub ;
} Spcgrp ;

#define NSPACEGROUPS 1024
//Spcgrp spcgrps[NSPACEGROUPS] ;
//Tcl_HashTable spcgrpnTable ;
static int NspaceGrps = 230 ;

/* utility function declarations*/

static int memExc(int isnull, char *loc) ;
static int setString(char **s, const char *src) ;


typedef struct {
  char *word ;
  int  size ;
  int  nalloc ;
} Word ;
typedef struct {
  Word *words ;
  double *nums ;
  int *isnum ;
  int nwords ;
  int nnum ;
  int nalloc ;
} WordList ;


/*--------- function declarations for vector stuff --------------------------*/

static void vecpro( double vec1[], double vec2[], double prod[] ) ;
static double dotpro( double vec1[], double vec2[] ) ;
static void sclpro( double scl, double vec[], double prod[] ) ;
static double recvec( double hkl[], double uni[] ) ;
static double vecmag( double vec[] ) ;
static double unitvec( double vec[] ) ;
static void matvec( double mm[3][3], double vec[], double prod[] ) ;
/*---------------------------------------------------------------------------*/

#define SQ_CHAR(c)	((isalnum(UCHAR(c))) || \
	(c == '_') || (c == ':') || (c == '@') || (c == '.'))

static int readWyck(char *fil, Spcgrp *grps, int Ngrps) ;

#endif /*PYSQ_H*/





