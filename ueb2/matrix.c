#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Eintrag (row,col) aus der Matrix a auslesen */
double get_entry(double* a, int ldim, int row, int col){
	return a[ldim*col+row];
}

/* Eintrag (row,col) der Matrix a als value setzen */
void set_entry(double* a, int ldim, int row, int col, double value){
	a[ldim*col+row] = value;
}

/* Matrix auf dem Bildschirm ausgeben */
void print_matrix(double* a, int rows, int cols){
	int r,c;
	for(r = 0; r < rows; r++)
	{
		for(c = 0; c < cols; c++)
		{
			printf("%.5f ", get_entry(a, rows, r, c));
		}
		printf("\n");
	}
}
 
/* Matrix-Vektor-Multiplikation (Ax = y) */
void mvm(double* a, int rows, int cols, double* x, double* y){
	int r,c;
	double temp;
	for(r = 0; r < rows; r++)
	{
		temp = 0.0;
		for(c = 0; c < cols; c++)
		{
			temp += get_entry(a, rows, r, c) * x[c];
		}
		y[r] = temp;
	}
}


int main(void){
  
  int i, j, n, m;
  double *a, *x, *y;
  time_t t;
  time(&t);
  srand(t);
  
  /* Matrixdimension festlegen (n = Zeilen, m = Spalten) */
  n = 3;
  m = 3;
  
  /* Speicher anfordern */
  a = (double*) malloc(n*m*sizeof(double));
  x = (double*) malloc(m*sizeof(double));
  y = (double*) malloc(n*sizeof(double));
  
  /* Matrix a erzeugen mit (double) Eintraegen zwischen 0 und 10 */
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      set_entry(a,n,i,j,(rand() % 200)/20.0);
  
  /* Vektor x erzeugen mit (double) Eintraegen zwischen 0 und 10 */
  for(i=0;i<m;i++)
    set_entry(x,m,i,0,(rand() % 200)/20.0); 
  
  /* Matrix a ausgeben */
  printf("\na:\n");
  print_matrix(a,n,m);
  
  /* Vektor x ausgeben */
  printf("\nx:\n");
  print_matrix(x,m,1);
  
  /* Matrix-Vektor-Produkt a*x = y */
  mvm(a,n,m,x,y);
  
  /* Vektor y ausgeben */
  printf("\ny:\n");
  print_matrix(y,n,1);
  
  /* Speicher freigeben */
  free(a);
  free(x);
  free(y);
  
  return 0;
}
