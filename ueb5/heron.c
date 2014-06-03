#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double heron(double w, double x0, double eps)
{
   int iter = 0;
   double x = x0;

   while (fabs(x*x - w) > eps)
   {
      x = 0.5 * (x + w/x);
      iter++;
   }
   x = x > 0 ? x : -x;

   printf("Number of iterations: %d \n", iter);
   return x;
}


int main(int argc, char* argv[])
{
   double res, w, x0, eps;
   w   = 2;
   eps = 10e-8;
   x0  = -100;

   res = heron(w, x0, eps);
   printf("Result: %.15f\n", res);
   printf("Deviation: %.15f\n", (sqrt(w)-res));

   return 0;
}