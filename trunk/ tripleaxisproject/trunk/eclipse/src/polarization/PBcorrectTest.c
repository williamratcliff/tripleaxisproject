/*
  PBcorrectTest cellsFile Datafile
  uses PBcorrectData to read cellsFile
  then uses PBreaddata to read the Datafile and its .CTRL flags
  putting the results in Datafile.out
*/

#define DEBUG 0

#include <stdlib.h>
#include <stdio.h>

#include "PB.h"


int main(int argc, char *argv[])
{

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


  if( argc < 3 ) {
    printf("PBcorrectTest requires 2 filenames (cellsFile & DATfile\n") ;
    exit(0) ;
  }

  if( PBcorrectData(argv[1], NULL, 0, NULL, NULL) ) {
    printf("PBcorrectTest failed to read cell files.\n") ;
    exit(0) ;
  }

  if( (ierr = PBreaddata(argv[2])) ) {
    printf("PBcorrectTest failed PBreaddata= %d\n",ierr) ;
  }
  exit(0) ;
}
