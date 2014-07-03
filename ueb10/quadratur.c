#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M_PI 3.141592653589793

typedef double (*RealFunction)(double arg);

/* Funktion f1 in x auswerten */
double f1(double x);

/* Stammfunktion von f1 in x auswerten */
double F1(double x);

/* Trapezregel */
double CompositeTrapezoidRule(int l, double a, double b, double (*f)(double));

/*#######################################################################################################*/

int main(void){
  int l;
  double a, b, qe, qt;
  
  /* Intervall festlegen */
  a = 0.0;
  b = M_PI;
  
  /* Anzahl der Unterteilungen festlegen */
  l = 400;
  
	/* Berechnung des Integrals */
	qe = F1(b) - F1(a);
	qt = CompositeTrapezoidRule(l, a, b, f1);
  
	/* Ausgabe */
	printf("Exakt: %f\n", qe);
	printf("Trapezregel     : %f (Fehler: %f)\n", qt, fabs(qt - qe));
	
  return 0;
}

/*#######################################################################################################*/

double
f1(double x)
{
  return sin(x);
}

double
F1(double x)
{
  return -cos(x);
}

double
CompositeTrapezoidRule(int l, double a, double b, RealFunction f)
{
  int i;
  double h, sum = 0;

  h = (b-a)/l;

  for (i = 0; i < l; i++) {
    sum += (f(a + i*h) + f(a + i*h + h));
  }

  return sum*h/2;
}
