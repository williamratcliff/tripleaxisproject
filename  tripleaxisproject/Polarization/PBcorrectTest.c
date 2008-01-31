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
  line3 is He3CELLpol: PHe  PHeErr hrsBeam T  Terr
  line4 is data for these
  line5 is     He3exp: waveRelWidth angleVwidth angleHwidth usedRadius
  line6 is data for these
  line7 is efficiency: transportEff terr flipperEff ferr
  line8 is data for these

  reader calculates nsL corrected and nsLerr putting them into PBsetup
  so program input is 3 filenames  on cmnd line:
  Pcellfile Acellfile DATAfile [flag1 [flag2]]
  If there is a lambdaFlag(L) the init and final Energies are wavelengths in
  Angstroms instead of the default meV
  There may also be a flag indicating that counts were obtained using
  a normalizing monitor place after the polarizer cell.
  This is not very good since such a monitor will also count lambda/2.
  If no monitor is available before the polarizing cell, it is probably best
  to count against time.
*/

#define DEBUG 0

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "PB.h"

static double uniform(int idum, double d1, double d2, double *distval,
		      double min, double max) ;
static double gaussian(int idum, double ave, double sig, double *distval,
		       double d1, double d2) ;
static double poisson(int idum, double d1, double d2, double *distval,
		      double min, double max) ;

int PBcoef(PBdatapt *d, PBsetup *s) ;
int PBcorrect(PBdatapt *d) ;
int applyConstraint(PBdatapt *d, ConstraintEq *eq) ;
int combineEqs(PBdatapt *d, int eq1, int eq2) ;
int deleteEq(PBdatapt *d, int eq) ;
int constrainResult(PBdatapt *d, int nc, ConstraintEq *eqs) ;

int main(int argc, char *argv[])
{
  static int MAXDAT = 1024 ;
  PBdatapt d[MAXDAT] ;
  ConstraintEq eqs[4] ;

  static double Dn = 2.072141789 ;
  static double TWOPI = 6.283185307 ;
  int wavelengthData ;
  int i, ierr, j, k, n, nr, nc, nw ;
  int nscan ;
  int activeEq[4], nactive ;
  int freeS[4], nfree ;
  int si, ic[4] ;
  double cf[4] ;
  char *cp ;
  char buf[512] ;
  char outfile[512] ;
  char startDate[256];
  char syscmnd[512];
  double N, err[4], cellfac ;
  double dval, d1, d2, dbl ;
  double iout, fout, kout;
  FILE *fp ;

  PBsetup s ;
  He3CELL *cell ;
  He3CELLpol *pol ;
  expResol *exper ;
  efficiency *eff ;

  if( argc < 4 ) {
    printf("PBcorrectTest requires 3 filenames (Pcell Acell & DATfile\n") ;
    exit(0) ;
  }
  wavelengthData = 0 ;
  if( argc > 4 && (*argv[4] == 'L' || *argv[4] == 'l') ) wavelengthData = 1 ;

  N = 1000. ;

  if( (fp = fopen(argv[1], "r")) == NULL ) {
    printf("failed to open %s\n", argv[1]) ;
    exit(0) ;
  }
  for( i=0 ; i<4 ; i++ ) {
    if( (fgets(buf, 511, fp)) == NULL ) {
      printf("failed to read from %s\n", argv[1]) ;
      fclose(fp) ;
      exit(0) ;
    }
    if( (fgets(buf, 511, fp)) == NULL ) {
      printf("failed to read from %s\n", argv[1]) ;
      fclose(fp) ;
      exit(0) ;
    }
    if( i == 0 ) {
      cell = &(s.P.cell) ;
      if( sscanf(buf,"%s %lf %lf %lf %lf %lf %lf %lf",
		 cell->name,&(cell->tEmpty),&(cell->tEmptySlope),
		 &(cell->L),&(cell->D),&(cell->R),
		 &(cell->nsL0),&(cell->nsL0err)) < 8 ) {
	printf("failed to read cell data from %s\n", argv[1]) ;
	fclose(fp) ;
	exit(0) ;
      }
    } else if( i == 1 ) {
      pol = &(s.P.pol) ;
      nscan = sscanf(buf,"%lf %lf %lf %lf %lf %ul",
		     &(pol->PHe),&(pol->PHeErr),&(pol->hrsBeam),
		     &(pol->T),&(pol->Terr),&(pol->startSecs)) ;
      if( nscan < 6 ) {
	if( nscan < 5 ) printf("failed to read cell data from %s\n", argv[1]) ;
	else printf("failed to read CELL startSecs from %s\n", argv[1]);
	fclose(fp) ;
	exit(0) ;
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
      exper = &(s.P.res) ;
      if( sscanf(buf,"%lf %lf %lf %lf",
		 &(exper->waveRelWidth),
		 &(exper->angleVwidth),&(exper->angleHwidth),
		 &(exper->usedRadius)) < 4 ) {
	printf("failed to read cell data from %s\n", argv[1]) ;
	fclose(fp) ;
	exit(0) ;
      }
    } else if( i == 3 ) {
      eff = &(s.eP) ;
      if( sscanf(buf,"%lf %lf %lf %lf",
		 &(eff->teff),&(eff->terr),
		 &(eff->feff),&(eff->ferr)) < 4 ) {
	printf("failed to read cell data from %s\n", argv[1]) ;
	fclose(fp) ;
	exit(0) ;
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

  if( (fp = fopen(argv[2], "r")) == NULL ) {
    printf("failed to open %s\n", argv[2]) ;
    exit(0) ;
  }
  for( i=0 ; i<4 ; i++ ) {
    if( (fgets(buf, 511, fp)) == NULL ) {
      printf("failed to read from %s\n", argv[2]) ;
      fclose(fp) ;
      exit(0) ;
    }
    if( (fgets(buf, 511, fp)) == NULL ) {
      printf("failed to read from %s\n", argv[2]) ;
      fclose(fp) ;
      exit(0) ;
    }
    if( i == 0 ) {
      cell = &(s.A.cell) ;
      if( sscanf(buf,"%s %lf %lf %lf %lf %lf %lf %lf",
		 cell->name,&(cell->tEmpty),&(cell->tEmptySlope),
		 &(cell->L),&(cell->D),&(cell->R),
		 &(cell->nsL0),&(cell->nsL0err)) < 8 ) {
	printf("failed to read cell data from %s\n", argv[2]) ;
	fclose(fp) ;
	exit(0) ;
      }
    } else if( i == 1 ) {
      pol = &(s.A.pol) ;
      nscan = sscanf(buf,"%lf %lf %lf %lf %lf %ul",
		     &(pol->PHe),&(pol->PHeErr),&(pol->hrsBeam),
		     &(pol->T),&(pol->Terr),&(pol->startSecs)) ;
      if( nscan < 6 ) {
	if( nscan < 5 ) printf("failed to read cell data from %s\n", argv[2]) ;
	else printf("failed to read CELL startSecs from %s\n", argv[2]);
	fclose(fp) ;
	exit(0) ;
      }
    } else if( i == 2 ) {
      exper = &(s.A.res) ;
      if( sscanf(buf,"%lf %lf %lf %lf",
		 &(exper->waveRelWidth),
		 &(exper->angleVwidth),&(exper->angleHwidth),
		 &(exper->usedRadius)) < 4 ) {
	printf("failed to read cell data from %s\n", argv[2]) ;
	fclose(fp) ;
	exit(0) ;
      }
    } else if( i == 3 ) {
      eff = &(s.eA) ;
      if( sscanf(buf,"%lf %lf %lf %lf",
		 &(eff->teff),&(eff->terr),
		 &(eff->feff),&(eff->ferr)) < 4 ) {
	printf("failed to read cell data from %s\n", argv[2]) ;
	fclose(fp) ;
	exit(0) ;
      }
    }
  }
  fclose(fp) ;

  cellfac = 1. ;
  if( cell->R > 0. && cell->L > 0. )
    cellfac *= (1. - exper->usedRadius*exper->usedRadius/cell->R/cell->L) ;
  cell->nsL = cellfac*cell->nsL0 ;
  cell->nsLerr = cellfac*cell->nsL0err ;


  /*
    ready to read the data values file
    for these tests limit to MAXDAT data lines
  */

  if( (fp = fopen(argv[3], "r")) == NULL ) {
    printf("failed to open %s\n", argv[3]) ;
    exit(0) ;
  }

  n = 0 ;
  for( i=0 ; i<4 ; i++ ) activeEq[i] = 1 ; /* initial all Eqs as active */

  while( fgets(buf, 511, fp) ) {
    /* find up to 10 words in buf and try to read each as double */
    cp = buf ;
    nw = 0 ;

    d[n].Nactive = 4 ;
    d[n].Nfree = 4 ;
    for( i=0 ; i<4 ; i++ ) {
      d[n].activeEq[i] = i ;
      d[n].freeS[i] = i ;
    }

    while(nw < 10 && *cp != '\n' && *cp != '\0') {
      /* skip white space */
      while( *cp == ' ' || *cp == '\t' ) cp++ ;
      /* break if end of line */
      if( *cp == '\n' || *cp == '\0' ) break ;
      /* try to scan this word as a double */


      if( sscanf(cp, "%lf", &dbl)<1 ) {
	if( nw < 2 ) {
	  printf("failed to read lambdaI or lambdaF at data pt # %d\n",n+1) ;
	  fclose(fp) ;
	  exit(0) ;
	} else {
	  /* mark this CS as unmeasured i.e. inactive Eq */
	  activeEq[nw/2-1] = 0 ;
	}
      } else {
	if( nw < 2 && ! wavelengthData ) {
	  dbl /= Dn ;
	  dbl = TWOPI/sqrt(dbl) ;
	}
	if( nw == 0 ) d[n].lambI = dbl ;
	else if( nw == 1 ) d[n].lambF = dbl ;
	else if( nw%2 == 0 ) d[n].Y[nw/2-1] = dbl ;
	else if( sscanf(cp, "%ul",&(d[n].sec[nw/2-1])) < 1 ) {
	  printf("failed to read time stamp as unsigned long pt # %d\n",n+1);
	  fclose(fp) ;
	  exit(0) ;
	}
      }
      /* go to just past this word */
      while( *cp != ' ' && *cp != '\t' && *cp != '\n' && *cp != '\0' ) cp++ ;
      nw++ ;
      if( *cp == '\n' || *cp == '\0' ) break ;
    }
    /* mark any missing words as inactiveEq */
    for( i=nw ; i<10 ; i++ ) activeEq[i/2-1] = 0 ;

    n++ ;
    if( n >= MAXDAT ) {
      printf("%d point limit reached\n",MAXDAT) ;
      break ;
    }
  }
  fclose(fp) ;

  nactive = 0 ;
  nfree = 4 ;
  for( j=0 ; j<4 ; j++ ) {
    if( activeEq[j] ) {
      activeEq[nactive] = j ;
      ++nactive ;
    }
    freeS[j] = j ;
  }

  for( i=0 ; i<n ; i++ ) {
    d[i].Nactive = nactive ;
    for( j=0 ; j<nactive ; j++ ) {
      d[i].activeEq[j] = activeEq[j] ;
      d[i].Yesq[activeEq[j]] = d[i].Y[activeEq[j]] ;
    }
  }

  /* Now go ahead and compute coefs */
  for( i=0 ; i<n ; i++ ) {
    if( (ierr = PBcoef(d+i, &s)) > 0 ) {
      printf("error in PBcoef = %d for point index= %d\n", ierr, i) ;
      //fclose(fp) ;
      exit(0) ;
    }
  }

  /*
    Now process any controls file = DATAfile.CTRLS
    start by printing nactive-activeEq and nfree-freeS
    read the controls file which allows following commands
    add a constraint on the underlying cross-sections
    C sindexToConstrain sindexFr coef ...
    add two equations
    A eqindex1 eqindex2

    print out the new nactive-activeEq and nfree-freeS after each cmnd
  */

  nc = 0 ;

  if( DEBUG ) {
    printf("starting activeEqindices and freeSindices are\n") ;
    for( i=0 ; i<nactive ; i++ ) printf("%1d ",activeEq[i]) ;
    printf("\n") ;
    for( i=0 ; i<nfree ; i++ ) printf("%1d ",freeS[i]) ;
    printf("\n\n") ;
  }

  strcpy(outfile, argv[3]) ;
  strcat(outfile, ".CTRLS") ;
  if( (fp = fopen(outfile, "r")) == NULL ) {
    if ( DEBUG ) printf("NO CTRLS file\n") ;
  } else {
    if ( DEBUG ) printf("processing CTRLS file %s\n",outfile) ;

    while( fgets(buf, 511, fp) ) {
      if( buf[0] == '#' ) continue ;
      if( buf[0] == 'C' || buf[0] == 'c' ) {
	/* constraint */
	if( (nr = sscanf(buf+1,"%d %d %lf %d %lf %d %lf",
			 &si, ic, cf, ic+1, cf+1, ic+2, cf+2)) < 3 ) continue;
	eqs[nc].freeToConstrain = si ;
	eqs[nc].Nfree = (nr - 1)/2;
	for( j=0 ; j<eqs[nc].Nfree ; j++ ) {
	  eqs[nc].freeS[j] = ic[j] ;
	  eqs[nc].Coef[j] = cf[j] ;
	}

	for( i=0 ; i<n ; i++ ) {
	  if( applyConstraint(d+i, eqs+nc) > 0 ) {
	    printf("constraint failed at point %d\n",i) ;
	  }
	}
	nc++ ;

	if ( DEBUG ) {
	  printf("after applying constraint activeEqindices and freeSindices are\n") ;
	  for( i=0 ; i<d[0].Nactive ; i++ ) printf("%1d ",d[0].activeEq[i]) ;
	  printf("\n") ;
	  for( i=0 ; i<d[0].Nfree ; i++ ) printf("%1d ",d[0].freeS[i]) ;
	  printf("\n\n") ;
	}

      } else if( buf[0] == 'A' || buf[0] == 'a' ) {
	/* add equations */
	if( (nr = sscanf(buf+1,"%d %d", ic, ic+1)) < 2 ) continue;
	if ( DEBUG ) printf("about to add equations %d and %d\n", ic[0], ic[1]) ;
	for( i=0 ; i<n ; i++ ) {
	  combineEqs(d+i, ic[0], ic[1]) ;
	}
	if ( DEBUG ) {
	  printf("after adding equations activeEqindices and freeSindices are\n") ;
	  for( i=0 ; i<d[0].Nactive ; i++ ) printf("%1d ",d[0].activeEq[i]) ;
	  printf("\n") ;
	  for( i=0 ; i<d[0].Nfree ; i++ ) printf("%1d ",d[0].freeS[i]) ;
	  printf("\n\n") ;
	  printf("first data point equations are:\n") ;
	  for( i=0 ; i<4 ; i++ ) {
	    printf("%11g   ",d[0].Y[i]) ;
	    for( j=0 ; j<4 ; j++ ) printf("%11g ",d[0].C[i][j]) ;
	    printf("\n") ;
	  }
	  printf("\n") ;
	}
      } else if( buf[0] == 'D' || buf[0] == 'd' ) {
	if( (nr = sscanf(buf+1,"%d %d %d", ic, ic+1, ic+2)) < 1 ) continue;
	if ( DEBUG ) printf("about to delete %d equations\n", nr) ;
	for( i=0 ; i<n ; i++ ) {
	  for( j=0 ; j<nr ; j++ ) deleteEq(d+i, ic[j]) ;
	}
	if ( DEBUG ) {
	  printf("after deleting equations activeEqindices and freeSindices are\n") ;
	  for( i=0 ; i<d[0].Nactive ; i++ ) printf("%1d ",d[0].activeEq[i]) ;
	  printf("\n") ;
	  for( i=0 ; i<d[0].Nfree ; i++ ) printf("%1d ",d[0].freeS[i]) ;
	  printf("\n\n") ;	
	}
      }
    }
    fclose(fp) ;
  }

  if( DEBUG ) {
    printf("after CTRLS activeEqindices and freeSindices are\n") ;
    for( i=0 ; i<d[0].Nactive ; i++ ) printf("%1d ",d[0].activeEq[i]) ;
    printf("\n") ;
    for( i=0 ; i<d[0].Nfree ; i++ ) printf("%1d ",d[0].freeS[i]) ;
    printf("\n\n") ;
  }

  /* Now go ahead and do the correction */
  for( i=0 ; i<n ; i++ ) {

    /* call PBcorrect to get the cross-sections */
    if( (ierr = PBcorrect(d+i)) > 0 ) {
      if( ierr == 2 ) {
	printf("invalid Nactive equations or Nactive != Nfree\n") ;
	exit(0) ;
      }
      printf("error in PBcorrect = %d for point index= %d\n", ierr, i) ;
      exit(0) ;
    }

    /* apply any constraints to result */
    constrainResult(d+i, nc, eqs) ;

    for( j=0 ; j<4 ; j++ ) err[j] = sqrt(d[i].Sesq[j]) ;

    iout = d[i].lambI ;
    fout = d[i].lambF ;
    if( ! wavelengthData ) {
      kout = TWOPI/iout ;
      iout = Dn*kout*kout ;
      kout = TWOPI/fout ;
      fout = Dn*kout*kout ;
    }

    printf("%8g %8g   ",iout,fout) ;
    for( j=0 ; j<4 ; j++ ) printf("%10g %10g %10d  ",d[i].S[j],err[j],d[i].sec[j]) ;
    printf("\n") ;
  }
  //fclose(fp) ;
  exit(0) ;
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
