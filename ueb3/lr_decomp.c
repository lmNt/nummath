#include <stdio.h>
#include <stdlib.h>

/* Eintrag (row,col) aus der Matrix a auslesen */
double get_entry(double* a, int ldim, int row, int col);

/* Eintrag (row,col) der Matrix a als value setzen */
void set_entry(double* a, int ldim, int row, int col, double value);

/* Matrix auf dem Bildschirm ausgeben */
void print_matrix(double* a, int rows, int cols);

/* Berechnet die lr-Zerlegung der Matrix a und ueberspeichert a */
void lr_decomp(double* a, int n, int ldim);

/* Loest das Gleichungssystem a*x = b mit einer Matrix a in lr-Darstellung ohne Pivotsuche */
void solve_lr_decomp(double* a, int n, int ldim, double* b, double* x);

/*#######################################################################################################*/

int main(void){

   int n;
   double *a;

   /* Matrixdimension festlegen */
   n = 4;

   /* Speicher anfordern */
   a = (double*) malloc(n*n*sizeof(double)); /* Matrix */

   /* Testmatrix A erzeugen  */
   set_entry(a, 4, 0, 0, 2);
   set_entry(a, 4, 0, 1, -1);
   set_entry(a, 4, 0, 2, -3);
   set_entry(a, 4, 0, 3, 3);
   set_entry(a, 4, 1, 0, 4);
   set_entry(a, 4, 1, 1, 0);
   set_entry(a, 4, 1, 2, -3);
   set_entry(a, 4, 1, 3, 1);
   set_entry(a, 4, 2, 0, 6);
   set_entry(a, 4, 2, 1, 1);
   set_entry(a, 4, 2, 2, -1);
   set_entry(a, 4, 2, 3, 6);
   set_entry(a, 4, 3, 0, -2);
   set_entry(a, 4, 3, 1, -5);
   set_entry(a, 4, 3, 2, 4);
   set_entry(a, 4, 3, 3, 1);
   

   /* Matrix ausgeben */
   printf("\na:\n");
   print_matrix(a,n,n);

   /* lr-Zerlegung berechnen */
   lr_decomp(a,n,n);

   /* Matrix ausgeben */
   printf("\nlr:\n");
   print_matrix(a,n,n);  

   /* Speicher freigeben */
   free(a);

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

void lr_decomp(double* a, int n, int ldim){
   int k, i, j;
   
   for (k=0; k<n; k++)
   {
      for (i=k+1; i<n; i++)
      {
         set_entry(a, n, i, k, (get_entry(a, n, i, k) / get_entry(a, n, k, k)));
         for (j=k+1; j<n; j++)
         {
            set_entry(a, n, i, j, (get_entry(a, n, i, j) - get_entry(a, n, i, k)*get_entry(a, n, k, j)));
         }
      }
   }
}

void solve_lr_decomp(double* a, int n, int ldim, double* b, double* x){
   int i, j;
   double* y;

   y = (double*) malloc(n*sizeof(double));

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
