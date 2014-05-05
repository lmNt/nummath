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



/* Loest ein lineares Gleichungssystem durch Vorwaertseinsetzen */
void forward_subst(int ldim, double* l, double* b, double *x){
	int j,i;
	for(j = 0; j < ldim; j++)
	{
		x[j] = b[j]/get_entry(l, ldim, j, j);
		for(i = j+1; i < ldim; i++)
		{
			b[i] = b[i] - (get_entry(l, ldim, i, j)*x[j]);
		}
	}
}

int main(void){
  
  int i, j, n;
  double *a, *b, *x, *y;
  time_t t;
  time(&t);
  srand(t);
  
  /* Matrixdimension festlegen */
  n = 3;
  
  /* Speicher anfordern */
  a = (double*) malloc(n*n*sizeof(double)); /* Matrix */
  b = (double*) malloc(n*sizeof(double));   /* Rechte Seite */
  x = (double*) malloc(n*sizeof(double));   /* Loesungsvektor */
  y = (double*) malloc(n*sizeof(double));   /* Pruefvektor */
  
  /* Untere Dreiecksmatrix erzeugen mit Eintraegen zwischen 1 und 10 */
  for(i=0;i<n;i++)
    for(j=0;j<i+1;j++)
      /* Keine 0 auf der Diagonalen zulassen */
      set_entry(a,n,i,j,(rand() % 90)/10.0 + 1.0);
      
  /* Rechte Seite erzeugen mit Eintraegen zwischen 0 und 10 */
  for(i=0;i<n;i++)
    set_entry(b,n,i,0,(rand() % 200)/20.0);
  
  /* Matrix ausgeben */
  printf("\na:\n");
  print_matrix(a,n,n);
  
  /* Rechte Seite ausgeben */
  printf("\nb:\n");
  print_matrix(b,n,1);
  
  /* Loesen durch Vorwaertseinsetzen */
  forward_subst(n,a,b,x);
  
  /* Loesung ausgeben */
  printf("\nx:\n");
  print_matrix(x,n,1);
  
  /* Probe */
  mvm(a,n,n,x,y);
  
  /* Pruefvektor ausgeben */
  printf("\ny:\n");
  print_matrix(y,n,1);
  
  /* Speicher freigeben */
  free(a);
  free(b);
  free(x);
  free(y);
  
  return 0;
}
