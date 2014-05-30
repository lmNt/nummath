#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x)
{
   double w = 5;
   return (x*x-w);
}

double bisection(double a, double b, double eps)
{
   double fa, fb, fx, x;
   fa = f(a);
   fb = f(b);

   while (b-a > eps)
   {
      x  = 0.5*(a+b);
      fx = f(x);

      if (fa*fx < 0)
      {
         b  = x;
         fb = fx;
      }
      else
      {
         a  = x;
         fa = fx;
      }
   }
   return x;
}


int main(int argc, _TCHAR* argv[])
{
   double a, b, eps, res;
   a = 0.0;
   b = 5.0;
   eps = 10e-8;
   res = bisection(a,b,eps);
   printf("Result: %.15f\n", res);
   return 0;
}