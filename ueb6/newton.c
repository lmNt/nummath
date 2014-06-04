#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*Modellierung vektorwertiger Funktionen*/
typedef void VectorFunctionName(double *x, double *fx, int n);
typedef VectorFunctionName *VectorFunction;

typedef void MatrixFunctionName(double *x, double *dfx, int n);
typedef MatrixFunctionName *MatrixFunction;

/*Dies sind Funktionen aus vorherigen Aufgaben*/

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

/* Loest ein lineares Gleichungssystem durch Rueckwaertseinsetzen */
void backward_subst(double* r, int n, int ldim, double* b, double* x);

/* Multiplikation von b mit Q* */
void qr_transform(double* qr, int m, int n, int ldim, double* b);

/* Loest das Gleichungssystem a*x = b mit einer Matrix a in qr-Darstellung */
void solve_qr_decomp(double* qr, int m, int n, int ldim, double* b, double *x);


/*Hier beginnen die neuen Funktionen*/

/*Beispiel aus T10*/

/*Berechnet fx = f(x)*/
void F_T10(double *x, double *fx, int n);

/*Berechnet dfx = Df(x)*/
void DF_T10(double *x, double *dfx, int n);

/*Berechnet die euklidische Norm eines Vektors*/
double norm2(double *x, int n);

/* Newton-Verfahren angwendet auf die Funktion f mit Jacobi-Matrix df, Startvektor x, Genauigkeit eps
 und maximal maxSteps Schritten
 Rückgabe: Anzahl der benötigten Schritte bis zur Konvergenz*/
int newton(double *x, VectorFunction f, MatrixFunction df,double eps, int maxSteps, int n);

/*#######################################################################################################*/

int main(void)
{
    double  eps;
    int maxSteps, n, Steps,i;
    double *start;
    VectorFunction f;
    MatrixFunction df;

    /*Dimension festlegen*/
    n = 2;

    /* Genauigkeit festlegen */
    eps = 1e-12;

    /*Anzahl der maximalen Schritte festlegen*/
    maxSteps = 15;

    /*Startvektor erzeugen und belegen*/
    start = (double *) malloc(n * sizeof(double));

    /*Für T10*/
    start[0] = 1.0;
    start[1] = 0.0;

    /*Beispiel aus T10*/
    f = F_T10;
    df = DF_T10;

    /*Beispiel aus T11*/
    /*f = F_T11;
    df = DF_T11;*/

    /* Newton-Verfahren durchführen */
    Steps = newton(start, f, df, eps, maxSteps, n);

    printf("Nullstelle:\n");
    for(i=0; i<n; i++)
        printf("%12.6f ",start[i]);
    printf("\n");
    printf("Anzahl der Schritte = %u\n",Steps);

    /*Speicher freigeben*/
    free(start);

    return 0;
}

/*#######################################################################################################*/



double norm2(double *x, int n)
{
    double norm;
    
    norm = sqrt(x[0]*x[0] + x[1]*x[1]);

    return norm;
}


void F_T10(double *x, double *fx, int n)
{
    fx[0] = cos(x[1]) - 0.2 - exp(1-x[0]);
    fx[1] = sin(x[1]) + 0.2 - x[0]*x[0] - x[1] + x[0]*(1.0 + x[1]);
}

void DF_T10(double *x, double *dfx, int n)
{
    set_entry(dfx, n, 0, 0, exp(1.0 - x[0]));
    set_entry(dfx, n, 0, 1, -sin(x[1]));
    set_entry(dfx, n, 1, 0, -2*x[0] + x[1] + 1.0);
    set_entry(dfx, n, 1, 1, cos(x[1]) + x[0] - 1.0);
}

int newton(double *x, VectorFunction f, MatrixFunction df, double eps, int maxSteps, int n)
{
    int Steps;
    int i;
    double *fx;
    double *dfx;
    double *d;
    double *x_old;

    d = malloc(n*sizeof(double));
    fx = malloc(n*sizeof(double));
    dfx = malloc(n*n*sizeof(double));

    Steps = 0;
    f(x, fx, n);

    while (norm2(fx, n) > eps) {
      f(x, fx, n);
      df(x, dfx, n);

      printf("fx\n");
      print_matrix(fx, n, 1);

      printf("Dfx\n");
      print_matrix(dfx, n, n);

      qr_decomp(dfx, n, n, n);

      for (i=0; i<n; i++) {
        fx[i] = -1 * fx[i];
      }
      
      solve_qr_decomp(dfx, n, n, n, fx, d);

      for (i=0; i<n; i++) {
        x[i] = x[i] + d[i];
      }

      if (Steps == maxSteps) {
        break;
      }
      Steps++;
    }
    
    
    return Steps;
}

double get_entry(double* a, int ldim, int row, int col)
{
    return a[col*ldim + row];
}

void set_entry(double* a, int ldim, int row, int col, double value)
{
    a[col*ldim + row] = value;
}

void print_matrix(double* a, int rows, int cols)
{
    int i, j;

    for(i=0; i<rows; i++) {
        for(j=0; j<cols; j++)
            printf(" %.2f\t ",get_entry(a,rows,i,j));
        printf("\n");
    }
}

int min_int(int n, int m)
{
    return (n >= m ? m : n);
}

void qr_decomp(double* a, int rows, int cols, int ldim)
{
    double rho, s, c, tau, alpha;
    int i, k, j;
    for (k = 0; k < min_int(rows, cols); k++) {
        for (i = k + 1; i < rows; i++) {
            if (get_entry(a, ldim, i, k) == 0) {
                rho = 1.0;
                c = 1.0;
                s = 0.0;
            } else if (fabs(get_entry(a, ldim, k, k)) >= fabs(get_entry(a, ldim, i, k))) {
                tau = get_entry(a, ldim, i, k) / get_entry(a, ldim, k, k);
                rho = tau / (sqrt(tau*tau + 1.0));
                s = rho;
                c = sqrt(1.0 - (s*s));
            } else {
                tau = get_entry(a, ldim, k, k) / get_entry(a, ldim, i, k);
                rho = sqrt(tau*tau + 1.0) / tau;
                c = 1.0 / rho;
                s = sqrt(1.0 - (c*c));
            }
            set_entry(a, ldim, k, k, c*get_entry(a, ldim, k, k) + s*get_entry(a, ldim, i, k));
            set_entry(a, ldim, i, k, rho);

            for (j = k + 1; j < cols; j++) {
                alpha = get_entry(a, ldim, k, j);
                set_entry(a, ldim, k, j, c*alpha + s * get_entry(a, ldim, i, j));
                set_entry(a, ldim, i, j, -s*alpha + c* get_entry(a, ldim, i, j));
            }
        }
    }
}

void backward_subst(double* r, int n, int ldim, double* b, double* x)
{
    int i, j;
    for (j = ldim - 1; j >= 0; j--) {
        x[j] = b[j] / get_entry(r, ldim, j, j);
        for (i = 0; i < j; i++) {
            b[i] = b[i] - get_entry(r, ldim, i, j) * x[j];
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

    backward_subst(qr, n, ldim, b, x);
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

