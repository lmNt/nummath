#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Eintrag (row,col) aus der Matrix a auslesen */
double get_entry(double* a, int ldim, int row, int col);

/* Eintrag (row,col) der Matrix a als value setzen */
void set_entry(double* a, int ldim, int row, int col, double value);

/* Matrix auf dem Bildschirm ausgeben */
void print_matrix(double* a, int rows, int cols);

/* Berechnet die lr-Zerlegung einer tridiagonalen Matrix mit unterer Nebendiagonalen ld, Diagonalen d und oberer
 * Nebendiagonalen ud und ueberspeichert diese mit der lr-Zerlegung */
void lr_decomp_tridiag(double* ld, double* d, double* ud, int n);

/* Loest das Gleichungssystem a*x = b mit einer tirdiagonalen Matrix a in lr-Darstellung ohne Pivotsuche */
void solve_lr_decomp_tridiag(double* ld, double *d, double *ud, int n, double* b, double* x);

/* Inverse Iteration unter Verwendung der lr-Zerlegung für tridiagonale Matrizen auf eine tridiagonale Matrix 
 *anwenden und bei Genauigkeit eps abbrechen */
void invit_adaptive_tridiag(double* ld, double* d, double* u, int n, double* x, double eps);

/* Berechnet die euklidische Norm */
double norm_2(double* x, int n);

/* Berechnet das euklidische Skalarprodukt */
double scalar_2(double* x, double* y, int n);

void div_vec_by_const(double* res, double* vec, const double c, int n);


/*#######################################################################################################*/

int main(void){
  
  int n, i;
  double *ld, *d, *ud, *x, *b;
  double eps, h, ih;
  FILE *F; 
  /* Matrixdimension festlegen */
  n = 100;
  
  /* Genauigkeit der Vektoriteration festlegen */
  eps = 1e-10;
  
  /* Speicher anfordern */
  ld = (double*) malloc((n-1)*sizeof(double)); 
  d  = (double*) malloc (n   *sizeof(double));
  ud = (double*) malloc((n-1)*sizeof(double));
  x  = (double*) malloc (n   *sizeof(double));

  /* Testmatrix A erzeugen  */

  for(i=0;i<=n-2;i++){
    d[i] = 2;
    ld[i] = -1;
    ud[i] = -1;
  }
   d[n-1] = 2;

  /* Zufaelligen Startvektor erzeugen */
  for(i=0;i<n;i++)
    x[i] = (rand() % 90)/10.0 + 1.0;

  /* Matrix ausgeben */
  printf("\nld:\n");
  print_matrix(ld,n-1,1);
  printf("\nd:\n");
  print_matrix(d,n,1);
  printf("\nud:\n");
  print_matrix(ud,n-1,1);
  

  /* lr-Zerlegung berechnen */
  //lr_decomp_tridiag(ld, d, ud,n);


  /* Matrix ausgeben */
 
  printf("lr-zerlegung\n");
  printf("\nld:\n");
  print_matrix(ld,n-1,1);
  printf("\nd:\n");
  print_matrix(d,n,1);
  printf("\nud:\n");
  print_matrix(ud,n-1,1);

  /*Inverse Iteration anwenden*/
  invit_adaptive_tridiag(ld, d, ud, n, x, eps);
  
    /* Loesung in Datei schreiben */
  F = fopen("./solution.txt", "w");
  h = 1.0/(n+1);
  fprintf(F,"%f %f\n",0.0,0.0);
  for(i=1;i<=n;i++){
    ih = i*h;
    fprintf(F,"%f %f\n",ih,x[i-1]);
  }
  fprintf(F,"%f %f\n",1.0,0.0);
  fclose(F);
  
  /* Speicher freigeben */
  free(d);
  free(ld);
  free(ud);
  free(x);
  return 0;
}

/*#######################################################################################################*/

double get_entry(double* a, int ldim, int row, int col){
  return a[col*ldim + row];
}

void set_entry(double* a, int ldim, int row, int col, double value){
  a[col*ldim + row] = value;
}

void print_matrix(double* a, int rows, int cols){
  int i, j;
  
  for(i=0;i<rows;i++){
    for(j=0;j<cols;j++)
      printf(" %.2f\t ", get_entry(a,rows,i,j));
    printf("\n");
  }
}

void lr_decomp_tridiag(double* ld, double* d, double* ud, int n)
{
  int i;
  double *a, *b;
  a = (double*) malloc((n-1)*sizeof(double));
  b = (double*) malloc(n*sizeof(double));

  b[0] = d[0];

  for (i = 1; i < n; i++) {
    a[i-1] = ld[i-1] / b[i-1];
    b[i] = d[i-1] - a[i-1] * ud[i-1];
  }

  memcpy(ld, a, (n-1)*sizeof(double));
  memcpy(d, b, n*sizeof(double));
}

void solve_lr_decomp_tridiag(double* ld, double* d, double* ud, int n, double* b, double* x)
{
  int i;
  double *y;
  y = (double*) malloc(n*sizeof(double));

  /* Forward substitution */
  y[0] = b[0];
  for (i = 1; i < n; i++) {
    y[i] = b[i] - ld[i-1]*y[i-1];
  }

  /* Backward substitution */
  x[n-1] = y[n-1] / d[n-1];
  for (i = n-2; i >= 0; i--) {
    x[i] = (y[i] - ud[i]*x[i+1])/d[i];
  }
}

double norm_2(double* x, int n){
   return sqrt(scalar_2(x,x,n));
}

double scalar_2(double* x, double* y, int n){
   int i;
   double res = 0;

   for (i = 0; i < n; i++)
   {
      res += x[i]*y[i];
   }
   return res;
}

void div_vec_by_const(double* res, double* vec, const double c, int n)
{
   int i;
   for (i = 0; i < n; i++)
   {
      res[i] = vec[i]/c;
   }
}

void invit_adaptive_tridiag(double* ld, double* d, double* u, int n, double* x, double eps)
{
  double *y, *test, *b;
  double gamma, lambda;
  int i;

  y = (double*) malloc(n*sizeof(double));
  test = (double*) malloc(n*sizeof(double));
  b = (double*) malloc(n*sizeof(double)); 

  gamma = norm_2(x, n);
  div_vec_by_const(x, x, gamma, n);

  lr_decomp_tridiag(ld, d, u, n);
  solve_lr_decomp_tridiag(ld, d, u, n, x, y);

  lambda = scalar_2(y,x,n);

  while (1)
  {
    for (i = 0; i < n; i++)
    {
      test[i] = lambda*x[i] - y[i];
    }
    if ( norm_2(test, n) <= (eps * norm_2(y,n)) ) break;

    gamma = norm_2(y,n);
    div_vec_by_const(x, y, gamma, n);

    solve_lr_decomp_tridiag(ld, d, u, n, x, y);
    lambda = scalar_2(y,x,n);
  }

  /* Hilfsspeicher freigeben */
  free(y);
  free(test);
  free(b);
}


















