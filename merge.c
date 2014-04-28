#include <stdio.h>
#include <stdlib.h>
#include <time.h>
	
void my_mergesort(int n, int* a)
{
  int m1, m2, i, k, l;
  int *b, *c;
  
  if(n>1)
  {
    m1 = (int) n/2;
    m2 = n - m1;
    
    /*Teilen der Liste a in die Listen b und c*/
    b = (int*) malloc(m1*sizeof(int));
    c = (int*) malloc(m2*sizeof(int));
	 
	 for(i = 0; i < m1; i++)
	 {
		 b[i] = a[i];
	 }
    
	 for(i = 0; i < m2; i++)
	 {
		 c[i] = a[i+m1];
	 }
	 
	 my_mergesort(m1, b);
	 my_mergesort(m2, c);
	 
	 i = k = l = 0;
	 
	 while(i < n)
	 {
		 while((k < m1) && (l==m2 || b[k]<=c[l]))
		 {
			 a[i] = b[k];
			 k++;
			 i++;
		 }
		 while((l < m2) && (k==m1 || c[l] <= b[k]))
		 {
			 a[i] = c[l];
			 l++;
			 i++;
		 }
	 }
	 
    /* Speicher wieder freigeben */
    free(c);
    free(b);
  }

}

/********** Hauptprogramm **********/
int main(void)
{
  int n, i;
  int* a;
  time_t t;
  time(&t);
  srand(t);
  
  n = 16;
  
  /* a mit Zufallszahlen zwischen 0 und n-1 füllen. */
  a = (int*) malloc(n*sizeof(int));
  for(i=0;i<n;i++)
  {
    a[i] = rand() % n;
  }
  
  /* Ausgabe der Eingabeliste */
  printf("\nEingabe: ");
  for(i=0;i<n;i++){
    printf(" %d ",a[i]);
  }
  printf("\n\n");
  
  /* my_mergesort-Verfahren durchfuehren */
  my_mergesort(n,a);
  
  /*Ausgabe der Liste nach Sortierung*/
  printf("\nAusgabe: ");
  for(i=0;i<n;i++){
    printf(" %d ",a[i]);
  }
  printf("\n");
  
  /* Speicher freigeben */
  free(a);
  
  return 0;
}
