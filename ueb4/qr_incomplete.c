#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* Eintrag (row,col) aus der Matrix a auslesen */
double get_entry(double* a, int ldim, int row, int col);

/* Eintrag (row,col) der Matrix a als value setzen */
void set_entry(double* a, int ldim, int row, int col, double value);

/* Matrix auf dem Bildschirm ausgeben */
void print_matrix(double* a, int rows, int cols);

/* Berechnet das Minimum der Zahlen n und m */
int min_int(int n, int m);

/* Berechnet die qr-Zerlegung der (rows x cols)-Matrix a */
void qr_decomp(double* a, int rows, int cols, int ldim);

/* Testen der qr-Zerlegung */
void test_qr(double* a, double* qr, int rows, int cols, int ldim);

/*#######################################################################################################*/

int main(void){
  
  int i, j, rows, cols;
  double *a, *qr;
  time_t t;
  time(&t);
  srand(t);
  
  /* Matrixdimension festlegen (rows = Zeilen, cols = Spalten) */
  rows = 3;
  cols = 3;
  
  /* Speicher anfordern */
  a = (double*) malloc(rows*cols*sizeof(double));
  qr = (double*) malloc(rows*cols*sizeof(double));
  
  /* Matrix a erzeugen mit (double) Eintraegen zwischen 0 und 10 */
  for(i=0;i<rows;i++)
    for(j=0;j<cols;j++){
      set_entry(a,rows,i,j,(rand() % 200)/20.0);
      set_entry(qr,rows,i,j,get_entry(a,rows,i,j));
    }
  
  /* rows und cols ausgeben */
  printf("\nrows = %d, cols = %d\n",rows,cols);
  
  /* Minimum von rows und cols ausgeben */
  printf("Minimum = %d\n",min_int(rows,cols));
  
  /* Matrix a ausgeben */
  printf("\na:\n");
  print_matrix(qr,rows,cols);
  
  /* qr-Zerlegung */
  qr_decomp(qr,rows,cols,rows);
  
  /* qr-Zerlegung ausgeben */
  printf("\nqr:\n");
  print_matrix(qr,rows,cols);
    
  /* Speicher freigeben */
  free(a);
  free(qr);
  
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

int min_int(int n, int m){
   return (n >= m ? m : n);
}

void qr_decomp(double* a, int rows, int cols, int ldim){
   double rho, s, c, tau, alpha;
   int i, k, j;
   for (k = 0; k < min_int(rows, cols); k++)
   {
      for (i = k + 1; i < rows; i++)
      {
         if (get_entry(a, ldim, i, k) == 0)
         {
            rho = 1.0;
            c = 1.0;
            s = 0.0;
         }
         else if (fabs(get_entry(a, ldim, k, k)) >= fabs(get_entry(a, ldim, i, k)))
         {
            tau = get_entry(a, ldim, i, k) / get_entry(a, ldim, k, k);
            rho = tau / (sqrt(tau*tau + 1.0));
            s = rho;
            c = sqrt(1.0 - (s*s));
         }
         else
         {
            tau = get_entry(a, ldim, k, k) / get_entry(a, ldim, i, k);
            rho = sqrt(tau*tau + 1.0) / tau;
            c = 1.0 / rho;
            s = sqrt(1.0 - (c*c));
         }
         set_entry(a, ldim, k, k, c*get_entry(a, ldim, k, k) + s*get_entry(a, ldim, i, k));
         set_entry(a, ldim, i, k, rho);

         for (j = k + 1; j < cols; j++)
         {
            alpha = get_entry(a, ldim, k, j);
            set_entry(a, ldim, k, j, c*alpha + s * get_entry(a, ldim, i, j));
            set_entry(a, ldim, i, j, -s*alpha + c* get_entry(a, ldim, i, j));
         }
      }
   }
}

