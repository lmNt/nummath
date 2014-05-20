#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* Eintrag (row,col) aus der Matrix a auslesen */
double get_entry(double* a, int ldim, int row, int col);

/* Eintrag (row,col) der Matrix a als value setzen */
void set_entry(double* a, int ldim, int row, int col, double value);

/* Matrix auf dem Bildschirm ausgeben (m = rows, n = cols)*/
void print_matrix(double* a, int m, int n);

/* Berechnet das Minimum der Zahlen n und m */
int min_int(int n, int m);

/* Berechnet die qr-Zerlegung der (rows x cols)-Matrix a */
void qr_decomp(double* a, int m, int n, int ldim);

/* Loest ein lineares Gleichungssystem durch Rueckwaertseinsetzen */
void backward_subst(double* r, /*int n,*/ int ldim, double* b, double* x);

/* Multiplikation von b mit Q* */
void qr_transform(double* qr, int m, int n, int ldim, double* b);

/* Loest das Gleichungssystem a*x = b mit einer Matrix a in qr-Darstellung */
void solve_qr_decomp(double* qr, int m, int n, int ldim, double* b, double *x);

/*#######################################################################################################*/

int main(void){

   int i, n, j;
   double *a, *qr, *b, *x;
   time_t t;
   time(&t);
   srand(t);

   /* Matrixdimension festlegen */
   n = 5;

   int rows = 5;
   int cols = 5;

   /* Speicher anfordern */
   a = (double*)malloc(n*n*sizeof(double));
   qr = (double*)malloc(n*n*sizeof(double));
   b = (double*)malloc(n*sizeof(double));
   x = (double*)malloc(n*sizeof(double));


   /* Matrix a erzeugen mit (double) Eintraegen zwischen 0 und 10 */
   for (i = 0; i<rows; i++)
   {
      for (j = 0; j<cols; j++){
         set_entry(a, rows, i, j, (rand() % 200) / 20.0);
         set_entry(qr, rows, i, j, get_entry(a, rows, i, j));
      }
   }

   printf("\na:\n");
   print_matrix(a, n, n);

   for (i = 0; i<n; i++)
      b[i] = 0.0;

   b[0] = 1.0;
   b[n - 1] = 1.0;


   printf("\nb:\n");
   print_matrix(b, n, 1);

   /* qr-Zerlegung berechnen */
   qr_decomp(qr, n, n, n);


   printf("\nqr:\n");
   print_matrix(qr, n, n);

   /* Gleichungssystem mit qr-Zerlegung loesen */
   solve_qr_decomp(qr, n, n, n, b, x);

   /* Ergebnisse ausgeben */
   printf("\nLoesung mit qr:\n");
   print_matrix(x, n, 1);

   printf("\n");

   /* Speicher freigeben */
   free(a);
   free(qr);
   free(b);
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

void print_matrix(double* a, int m, int n){
   int i, j;

   for (i = 0; i<m; i++){
      for (j = 0; j<n; j++)
         printf(" %.2f\t ", get_entry(a, m, i, j));
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


void backward_subst(double* r, /*int n,*/ int ldim, double* b, double* x){
   int i, j;
   for (j = ldim - 1; j >= 0; j--)
   {
      x[j] = b[j] / get_entry(r, ldim, j, j);
      for (i = 0; i < j; i++)
      {
         b[i] = b[i] - get_entry(r, ldim, i, j) * x[j];
      }
   }
}


void qr_transform(double* qr, int m, int n, int ldim, double* b)
{
   double rho, c, s, alpha;
   int k, i;
   for (k = 0; k < min_int(m, n); k++)
   {
      for (i = k + 1; i < m; i++)
      {
         rho = get_entry(qr, ldim, i, k);
         if (rho == 1)
         {
            c = 1;
            s = 0;
         }
         else if (fabs(rho)<1)
         {
            s = rho;
            c = sqrt(1 - (s*s));
         }
         else
         {
            c = 1 / rho;
            s = sqrt(1 - c*c);
         }

         alpha = b[k];
         b[k] = c*alpha + s*b[i];
         b[i] = -s*alpha + c*b[i];
      }
   }
}

void solve_qr_decomp(double* qr, int m, int n, int ldim, double* b, double *x)
{
   double *z;
   int j, i;

   qr_transform(qr, m, n, ldim, b);
   printf("\nz:\n");
   print_matrix(b, n, 1);

   backward_subst(qr, ldim, b, x);
}
