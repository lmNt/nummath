#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdafx.h>

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

/* Berechnet die qr-Zerlegung der (rows x cols)-Matrix a */
void qr_decomp(double* a, int m, int n, int ldim);

/* Multiplikation von b mit Q* */
void qr_transform(double* qr, int m, int n, int ldim, double* b);

/* Multiplikation der (ma x na)-Matrix a mit der (mb x nb)-Matrix b */
/* Matrix transponieren (trans = 1), Matrix nicht transponieren (trans = 0) */
void m_mult(double* a, int ma, int na, int transa, double* b, int mb, int nb, int transb, double* c);

/* Berechnet die Frobenius-Norm der Matrix a */
double norm_frob(double* a, int ldim, int m, int n);

/* Berechnet die Matrix Q der orthogonalen Iteration aus der QR-Zerlegung qr */
void get_q(double* qr, int ldim, int m, int n, double* q_hat, int k);

/* Berechnet die Eigenvektoren zu den betragsgroessten k Eigenwerten mit Genauigkeit eps*/
void orthoit(double* a, int lda, int n, double* x, int k, double eps);

/*#######################################################################################################*/

int main(void){
  int n, k, i, j;
  double eps, h, ih;
  double *a, *x;
  
  time_t t;
  FILE *F; 
  time(&t);
  srand(t);
  
  /* Problemdimension festlegen */
  n = 50;
  
  /* Anzahl der Startvektoren festlegen */
  k = 4;
  
  /* Genauigkeit der Vektoriteration festlegen */
  eps = 1e-10;
  
  /* Speichernaforderungen */
  a = (double*) malloc(n*n*sizeof(double));
  x = (double*) malloc(n*k*sizeof(double));
  
  /* Zufaellige Startvektoren erzeugen und in Matrix x zusammenfassen*/
   for(i=0;i<n*k;i++)
     x[i] = (rand() % 90)/10.0 + 1.0;
  
  /* Matrix des 1d-Modellproblems erzeugen */
  build_matrix1d(a,n);
    
  /* Orthogonale Iteration durchfuehren */
  orthoit(a,n,n,x,k,eps); 
  
  /* Loesung in Datei schreiben */
  F = fopen("./solution.txt", "w");
  h = 1.0/(n+1);
  fprintf(F,"%f ",0.0);
  for(j=0;j<k;j++)
    fprintf(F,"%f ",0.0);
  fprintf(F,"\n");
  for(i=1;i<=n;i++){
    ih = i*h;
    fprintf(F,"%f ",ih);
    for(j=0;j<k;j++)
      fprintf(F,"%f ",get_entry(x,n,i-1,j));
    fprintf(F,"\n");
      
  }
  fprintf(F,"%f ",1.0);
  for(j=0;j<k;j++)
    fprintf(F,"%f ",0.0);
  fprintf(F,"\n");
  fclose(F);
  
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
  
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if(i==j)
        set_entry(a,n,i,j,2.0);
      else if(fabs(i-j)==1)
        set_entry(a,n,i,j,-1.0);
      else
        set_entry(a,n,i,j,0.0);
    }
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

void m_mult(double* a, int ma, int na, int transa, double* b, int mb, int nb, int transb, double* c){
  int i, j, k, mc, nc, a_cols, b_rows;
  double value;
  
  mc = transa*na + (1 - transa)*ma;
  nc = transb*mb + (1 - transb)*nb;
  
  a_cols = transa*ma + (1 - transa)*na;
  b_rows = transb*nb + (1 - transb)*mb;
  
  if(a_cols!=b_rows){
    printf("Fehler in m_mult: %d != %d\n",a_cols,b_rows);
    return;
  }
  
  for(i=0;i<mc;i++){
    for(j=0;j<nc;j++){
      value = 0.0;
      for(k=0;k<a_cols;k++)
        value = value + get_entry(a,ma,(1-transa)*i+transa*k,(1-transa)*k+transa*i)*get_entry(b,mb,(1-transb)*k+transb*j,(1-transb)*j+transb*k);
      set_entry(c,mc,i,j,value);
    }
  }
  
}

double norm_frob(double* a, int ldim, int m, int n){
  int i, j;
  double tmp;
  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
      tmp += pow(fabs(get_entry(a, ldim, i, j)), 2);
    }
  }
  return sqrt(tmp);
}

void get_q(double* qr, int ldim, int m, int n, double* q_hat, int k){
   
}

void orthoit(double* a, int lda, int n, double* x, int k, double eps){
  double *q_hat, *y, *lambda, *test;
  
  q_hat = (double*) malloc(n*k*sizeof(double));
  y = (double*) malloc(n*k*sizeof(double));
  lambda = (double*) malloc(k*k*sizeof(double));
  test = (double*) malloc(n*k*sizeof(double));  /* Matrix der Abbruchbedingung */
  
  /*########################*/
  /*# Quelltext einfuegen! #*/
  /*########################*/
  
  free(q_hat);
  free(y);
  free(lambda);
  free(test);
}
