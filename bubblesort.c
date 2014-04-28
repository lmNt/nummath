#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/********** bubblesort-Verfahren **********/
/* - Eingabeliste im Array a mit Laenge n.*/
/* - Eingabeliste wird mit der sortierten */
/*      Liste ueberschrieben.             */
/******************************************/
void bubblesort(int n, int* a)
{
	int sigma, h, i;
	while (1)
	{
		sigma = 0;
		for(i = 0; i < n-1; i++)
		{
			if(a[i]>a[i+1])
			{
				h = a[i+1];
				a[i+1] = a[i];
				a[i] = h;
				sigma++;
			}
		}
		if(sigma==0)
		{
			break;
		}
	}
}

/********** Hauptprogramm **********/
int main(void)
{
	int* a;
	int i, n;
	time_t t;
	time(&t);
	srand(t);
  
	n = 12;
  
  /* a mit Zufallszahlen zwischen 0 und n-1 füllen. */
	a = (int*) malloc(n*sizeof(int));
	for(i=0;i<n;i++)
	{
		a[i] = rand() % n;
	}
  
  /* Ausgabe der Eingabeliste */
	printf("Eingabe: ");
	for(i=0;i<n;i++){
		printf(" %d ",a[i]);
	}
	printf("\n");
  
  /* bubblesort-Verfahren durchfuehren */
	bubblesort(n,a);
  
  /*Ausgabe der Liste nach Sortierung*/
	printf("Ausgabe: ");
	for(i=0;i<n;i++){
		printf(" %d ",a[i]);
	}
	printf("\n");
  
  /* Speicher freigeben */
	free(a);
  
	return 0;
}
