#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#define M_PI 3.141592653589793

typedef double (*RealFunction)(double arg);

/* Funktion f1 in x auswerten */
double f1(double x);

/* Funktion f2 in x auswerten */
double f2(double x);

/* Funktion F in den Punkten x[0],...,x[m] auswerten und die Funktionswerte in f speichern */
void eval(RealFunction F, double *x, double *f, int m);

/* Newtons dividierte Differenzen aus den Interpolationspunkten x[0], ... x[m] un den Stützwerten f[0], .. f[m]
 *berechnen  und in f speichern*/
void newton_differenzen(double *x, double *f, int m);

/* Interpolant in y mit Newtons dividierten Differenzen d berechnen */
double eval_newton(double y, double *x, double *d, int m);

/* aequidistante Interpolationspunkte im Intervall [a,b] berechnen und in x speichern*/
void equidistant(double *x, int m, double a, double b);

/* Tschebyscheff-Interpolationspunkte im Intervall [a,b] berechnen und in x speichern*/
void tschebyscheff(double *x, int m, double a, double b);

/*#######################################################################################################*/

int main(void){
  int m, i, n;
  double a, b, h, y;
  double *x, *fx, *t, *ft;
  FILE *file;
  RealFunction F;
  
  /* Intervall festlegen */
  
  /*fuer f1*/
  //a = 0.0;
  //b = 1.0;
  
  /*fuer f2*/
  a = -6.0;
  b = 6.0;
  
  /* Polynomgrad */
  m = 12;
  
  /* Speicher anfordern */
  x = (double*) malloc((m+1)*sizeof(double));
  fx = (double*) malloc((m+1)*sizeof(double));
  t = (double*) malloc((m+1)*sizeof(double));
  ft = (double*) malloc((m+1)*sizeof(double));
  
  /* Interpolationspunkte setzen */
  
  /*aequidistant*/
  equidistant(x, m, a,b);
  
  /* Tchebyscheff*/
  tschebyscheff(t,m,a,b);

  /*zu interpolierende Funktion festlegen*/
  
  //F = &f1;
  F = &f2;
  
  /* Stuetzwerte berechnen */
  eval(F,x,fx,m);
  eval(F,t,ft,m);

  /* Testpunkt festlegen */
  y = 0.0;

  /* Dividierte Differenzen ausrechnen */
  newton_differenzen(x,fx,m);
  newton_differenzen(t,ft,m);
  
  /* Testausgabe */
  printf("\n-----------Newton---------\n");
  printf("Interpolant in %.2f : %f\n",y,eval_newton(y,t,ft,m));
  printf("Originalfunktion in %.2f : %f\n",y,F(y));
  printf("---------------------------\n\n");
  
  /* Loesung in Datei schreiben */
  file = fopen("./f2.txt", "w");
  n = 200;
  h = (b-a)/n;
  for(i=0;i<n+1;i++){
    y = a+i*h;
    fprintf(file,"%f %f %f %f\n",y,F(y),eval_newton(y,x,fx,m),eval_newton(y,t,ft,m));
  }  
  fclose(file);  
  printf("fertig\n");
  /* Speicher freigeben */
  free(x);
  free(fx);
  free(t);
  free(ft);
  return 0;
}

/*#######################################################################################################*/

double f1(double x)
{
   /* k<=1 */
   int k = 1; 
   return sin(k*M_PI*x);
}

double f2(double x)
{
   return 1.0/(1.0 + x*x);
}

void eval(RealFunction F, double *x, double *f, int m)
{
   int i;
   for (i = 0; i < m+1; i++)
   {
      f[i] = F(x[i]);
   }
}


void newton_differenzen(double *x, double *f, int m)
{
   int i, j, n;

   for (n = 1; n < m+1; n++)
   {
      for (j = m; j >= n; j--)
      {
         i = j - n;
         f[j] = (f[j]-f[j-1])/(x[j]-x[i]);
      }
   }
}

double eval_newton(double y, double *x, double *d, int m)
{
   double s;
   int i;

   s = d[m];
   for (i = m-1; i >= 0; i--)
   {
      s = d[i] + (y - x[i]) * s;
   }
   return s;
}

void equidistant(double *x, int m, double a, double b)
{
   int i;
   for (i = 0; i < m+1; i++)
   {
      x[i] = 0.5*(b+a) + 0.5*(b-a) * (-1 + (2.0/(double)m) * (double)i);
   }
}


void tschebyscheff(double *x, int m, double a, double b)
{
   int i;
   for (i = 0; i < m+1; i++)
   {
      x[i] = 0.5*(b+a) + 0.5*(b-a) * cos( M_PI * (2*i+1)/(2*m+2) );
   }
}