/*
  simTest Pcellfilename Acellfilename CSfile
  link this with the libPB.so shared lib
*/

#define DEBUG 0

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "PB.h"


/*
  simTest Pcellfile Acellfile CSfile
  call PBcorrectData(Pcellfile, Acellfile, NULL, 0, NULL, NULL)
*/

int main(int argc, char *argv[])
{
  int ierr ;

  if( argc < 4 ) {
    printf("simTest requires 3 filenames (Pcell Acell & DATfile\n") ;
    exit(0) ;
  }

  if( (ierr=PBcorrectData(argv[1], argv[2], NULL, 0, NULL, NULL)) ) {
    printf("PBcorrectData err= %d\n",ierr);
    exit(0) ;
  }

  if( (ierr = PBsim(argv[3])) ) {
    printf("PBsim failure = %d\n",ierr) ;
  }
  exit(0) ;
}
