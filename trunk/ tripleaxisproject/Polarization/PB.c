#include <stdio.h>
/* to compile
gcc -shared -opb.dll pb.c
*/
typedef struct {
  unsigned long sec[4] ;
  double Y[4];
} PBdatapt;

typedef struct {
  double Y[4][4];
} matrix;

int PBcorrect(PBdatapt* pt){
 pt->sec[0]=5;   
}

int correct(matrix* pt){
 int i=2,j=1;
 
 printf("hi %d\n",pt->Y[0][0]);
 pt->Y[i][j]=5.0; 
 printf("%d\n",pt->Y[i][j]);  
}
