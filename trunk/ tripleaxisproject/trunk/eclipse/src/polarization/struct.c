typedef struct {
   int n, k;
   double A[4][4], B[4][4];
} embed_array;
void call(embed_array *pv)
{
  pv->n = 5;
  pv->A[1][1] = 18.652;
}
