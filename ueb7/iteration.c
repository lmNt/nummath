#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

/* Eintrag (row,col) aus der Matrix a auslesen */
double get_entry(double* a, int ldim, int row, int col);

/* Eintrag (row,col) der Matrix a als value setzen */
void set_entry(double* a, int ldim, int row, int col, double value);

/* Matrix auf dem Bildschirm ausgeben */
void print_matrix(double* a, int rows, int cols);

/* Matrix-Vektor-Multiplikation (Ax = y) */
void mvm(double* a, int n, double* x, double* y);

/* Berechnet das Minimum der Zahlen n und m */
int min_int(int n, int m);

/* Erstellen der Matrix des 1d-Modellproblems */
void build_matrix1d(double* a, int n);

/* Loest ein lineares Gleichungssystem durch Rueckwaertseinsetzen */
void backward_subst(double* r, int n, int ldim, double* b, double* x);

/* Berechnet die qr-Zerlegung der (rows x cols)-Matrix a */
void qr_decomp(double* a, int m, int n, int ldim);

/* Multiplikation von b mit Q* */
void qr_transform(double* qr, int m, int n, int ldim, double* b);

/* Loest das Gleichungssystem a*x = b mit einer Matrix a in qr-Darstellung */
void solve_qr_decomp(double* qr, int m, int n, int ldim, double* b, double *x);

/* Berechnet die euklidische Norm */
double norm_2(double* x, int n);

/* Berechnet das euklidische Skalarprodukt */
double scalar_2(double* x, double* y, int n);

/* Vektoriteration auf a anwenden und bei Genauigkeit eps abbrechen */
void power_adaptive(double* a, int n, double* x, double eps);

/* Inverse Iteration auf a anwenden und bei Genauigkeit eps abbrechen */
void invit_adaptive(double* a, int n, double* x, double eps);

void div_vec_by_const(double* res, double* vec, const double c, int n);

/*#######################################################################################################*/

int main(int argc, char* argv[])
{
   int n, i;
   double eps, h, ih;
   double *a, *x;
   time_t t;
   
   FILE *F; 

   time(&t);
   srand(t);

   /* Problemdimension festlegen */
   n = 100;

   /* Genauigkeit der Vektoriteration festlegen */
   eps = 1e-10;

   /* Speicheranforderungen */
   a = (double*) malloc(n*n*sizeof(double));
   x = (double*) malloc(n*sizeof(double));

   /* Zufaelligen Startvektor erzeugen */
   for(i=0;i<n;i++)
      x[i] = (rand() % 90)/10.0 + 1.0;

   /* Matrix des 1d-Modellproblems erzeugen */
   build_matrix1d(a,n);

   //power_adaptive(a,n,x,eps);
   invit_adaptive(a,n,x,eps);

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
   
   /* Direktes plotten aus C in Windows falls gnuplot installiert */
#if 0
   FILE *pipe = _popen("pgnuplot -persist", "w");
   if (pipe)
   {
      fprintf(pipe , "plot \"./solution.txt\" with lines\n");
      fflush(pipe);
   }
   fclose(pipe);
#endif

   /* Speicher freigeben */
   free(a);
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
         printf(" %.2f\t ",get_entry(a,rows,i,j));
      printf("\n");
   }
}

void mvm(double* a, int n, double* x, double* y){
   int i, j;

   for(i=0;i<n;i++){
      y[i] = 0.0;
      for(j=0;j<n;j++)
         y[i] = y[i] + get_entry(a,n,i,j)*x[j];
   }
}

int min_int(int n, int m){
   if(n<m)
      return n;
   else
      return m;
}

void build_matrix1d(double* a, int n){
   int i, j;

   for (i = 0; i < n; i++)
   {
      for (j = 0; j < n; j++)
      {
         /* Hauptdiagonale */
         if (i==j)
         {
            set_entry(a, n, i, j, 2);
         }
         /* Untere und obere Nebendiagonale */
         else if (abs(i-j)==1)
         {
            set_entry(a, n, i, j, -1);
         }
         /* Rest */
         else
         {
            set_entry(a, n, i, j, 0);
         }
      }
   }
}

void backward_subst(double* r, int n, int ldim, double* b, double* x){
   int i, j;

   for(j=n-1;j>=0;j--){
      x[j] = b[j]/get_entry(r,ldim,j,j);
      for(i=0;i<j;i++)
         b[i] = b[i] - get_entry(r,ldim,i,j) * x[j];
   }
}

void qr_decomp(double* a, int m, int n, int ldim){

   /* Benoetigte Variablen */
   int k, i, j, min;
   double rho, tau, c, s, gamma;

   /* Minimum aus Zeilen und Spalten */
   min = min_int(n,m);

   /* qr-Algorithmus */
   for(k=0;k<min;k++){
      for(i=k+1;i<m;i++){

         if(get_entry(a,ldim,i,k)==0.0){
            rho = 1.0; c = 1.0; s = 0.0;
         }
         else if(fabs(get_entry(a,ldim,k,k))>=fabs(get_entry(a,ldim,i,k))){
            tau = get_entry(a,ldim,i,k)/get_entry(a,ldim,k,k);
            rho = tau/sqrt(tau*tau+1.0);
            s = rho;
            c = sqrt(1.0-s*s);
         }
         else{
            tau = get_entry(a,ldim,k,k)/get_entry(a,ldim,i,k);
            rho = sqrt(tau*tau+1.0)/tau;
            c = 1.0/rho;
            s = sqrt(1.0-c*c);
         }
         set_entry(a,ldim,k,k,c*get_entry(a,ldim,k,k)+s*get_entry(a,ldim,i,k));
         set_entry(a,ldim,i,k,rho);

         for(j=k+1;j<n;j++){
            gamma = get_entry(a,ldim,k,j);
            set_entry(a,ldim,k,j,c*gamma+s*get_entry(a,ldim,i,j));
            set_entry(a,ldim,i,j,-1.0*s*gamma+c*get_entry(a,ldim,i,j));
         } 
      }
   }
}

void qr_transform(double* qr, int m, int n, int ldim, double* b){
   int k, i, min;
   double c, s, rho, alpha;

   min = min_int(m,n);

   for(k=0;k<min;k++){
      for(i=k+1;i<m;i++){
         rho = get_entry(qr,ldim,i,k);

         if(rho==1.0){
            c = 1.0; s = 0.0;
         }
         else if(fabs(rho)<1.0){
            s = rho; c = sqrt(1.0 - s*s);
         }
         else{
            c = 1.0/rho; s = sqrt(1.0 - c*c);
         }

         alpha = b[k];
         b[k]  = c*alpha + s*b[i];
         b[i]  = -1.0*s*alpha + c*b[i];
      }
   }
}

void solve_qr_decomp(double* qr, int m, int n, int ldim, double* b, double *x){
   /* Q^*b berechnen */  
   qr_transform(qr,m,n,ldim,b);

   /* Gleichungssystem loesen */
   backward_subst(qr,m,ldim,b,x);
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

void power_adaptive(double* a, int n, double *x, double eps){
   double *y, *test;
   double gamma, lambda;
   int i;

   /* Hilfsspeicher */
   y = (double*) malloc(n*sizeof(double));
   test = (double*) malloc(n*sizeof(double)); //(lambda*x - y)

   gamma = norm_2(x, n);
   div_vec_by_const(x, x, gamma, n);
   mvm(a, n, x, y);
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
      mvm(a, n, x, y);
      lambda = scalar_2(y,x,n);
   }

   /* Hilfsspeicher freigeben */
   free(y);
   free(test);
}

void invit_adaptive(double* a, int n, double* x, double eps){
   double *y, *test, *b;
   double gamma, lambda;
   int i;

   /* Hilfsspeicher */
   y    = (double*) malloc(n*sizeof(double));
   test = (double*) malloc(n*sizeof(double)); // (lambda*x - y)
   b    = (double*) malloc(n*sizeof(double)); // x abspeichern, da rechte Seite bei solve_qr_decomp ueberschrieben wird

   gamma = norm_2(x, n);
   div_vec_by_const(x, x, gamma, n);
   
   memcpy(b, x, n*sizeof(double));
   qr_decomp(a, n, n, n);
   solve_qr_decomp(a, n, n, n, x, y);
   memcpy(x, b, n*sizeof(double));
   
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

      memcpy(b, x, n*sizeof(double));
      solve_qr_decomp(a, n, n, n, x, y);
      memcpy(x, b, n*sizeof(double));
      
      lambda = scalar_2(y,x,n);
   }

   /* Hilfsspeicher freigeben */
   free(y);
   free(test);
   free(b);
}
