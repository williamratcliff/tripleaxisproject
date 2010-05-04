#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct {
  unsigned long *y;
} Mystruct;

int PBcorrectData(Mystruct *in,int npts)
{
  /* storage for PB data */

  int i, ierr, j ;

  printf("reached\n");
  printf("npts %d\n",npts);
  printf("in\n");
  for (i=0; i <npts;i++){
      printf("%d %u\n",i,in->y[i]);
      }
}
