#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Eintrag (row,col) aus der Matrix a auslesen */
double get_entry(double* a, int ldim, int row, int col);

/* Eintrag (row,col) der Matrix a als value setzen */
void set_entry(double* a, int ldim, int row, int col, double value);

/* Matrix auf dem Bildschirm ausgeben */
void print_matrix(double* a, int rows, int cols);

/* Berechnet die lr-Zerlegung der Matrix a mit Pivotsuche, ueberspeichert a und merkt sich die Vertauschungen in p*/
void lr_pivot(double* a, int n, int ldim, int* p);

/* Loest das Gleichungssystem a*x = b mit einer Matrix a in lr-Darstellung mit Pivotsuche */
void solve_lr_pivot(double* a, int n, int ldim, int* p, double* b, double* x);

/*#######################################################################################################*/

int main(void){

   int n, i;
   double *a;
   int *p;

   /* Matrixdimension festlegen */
   n = 3;

   /* Speicher anfordern */
   a = (double*) malloc(n*n*sizeof(double)); /* Matrix */
   p = (int*) malloc(n*sizeof(int));         /* Permutationsvektor */

   /* Testmatrix B erzeugen  */
   set_entry(a, 3, 0, 0, 0);
   set_entry(a, 3, 0, 1, 1);
   set_entry(a, 3, 0, 2, 3);
   set_entry(a, 3, 1, 0, 4);
   set_entry(a, 3, 1, 1, 6);
   set_entry(a, 3, 1, 2, 0);
   set_entry(a, 3, 2, 0, 2);
   set_entry(a, 3, 2, 1, 6);
   set_entry(a, 3, 2, 2, 3);

   /* Matrix ausgeben */
   printf("\na:\n");
   print_matrix(a,n,n);

   /* lr-Zerlegung mit Pivotsuche berechnen */
   lr_pivot(a,n,n,p);

   /* Matrix ausgeben */
   printf("\nlr:\n");
   print_matrix(a,n,n);

   /* Permutationsvektor ausgeben */
   printf("\np:\n");
   for(i=0;i<n;i++)
      printf(" %d \n",p[i]);

   /* Speicher freigeben */
   free(a);
   free(p);

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
         printf(" %.2f\t ",get_entry(a,rows,i,j));
      printf("\n");
   }
}

void lr_pivot(double* a, int n, int ldim, int* p){
   int k, i, j, max;
   for (k=0; k<n; k++)
   {
      max = k;
      
   }
}

void solve_lr_pivot(double* a, int n, int ldim, int* p, double* b, double* x){
   int i, j;
   double gamma;
   double* y;

   y = (double*) malloc(n*sizeof(double));

   /*Vektor b umsortieren */
   for(i=0;i<n;i++){
      gamma = b[i];
      b[i] = b[p[i]];
      b[p[i]] = gamma;
   }

   /* Loese l*y = b */
   for(j=0;j<n;j++){
      y[j] = b[j];
      for(i=j+1;i<n;i++)
         b[i] = b[i] - get_entry(a,ldim,i,j) * y[j];
   }

   /* Loese r*x = y */
   for(j=n-1;j>=0;j--){
      x[j] = y[j]/get_entry(a,ldim,j,j);
      for(i=0;i<j;i++)
         y[i] = y[i] - get_entry(a,ldim,i,j) * x[j];
   }

   /* Speicher freigeben */
   free(y);
}
