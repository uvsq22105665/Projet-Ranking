#include "aitkenprim.h"

void estimation(int n , double* piK , double* piKmoins1 , double* piKmoins2 )
{
	double g;
	double h;
	
	//On stocke dans piKmoins2 la valeur de piK car on en a plus besoin, gain de memoire
	int i ;
	for(i = 0 ; i < n ; i++)
	{
		g = (piKmoins1[i] - piKmoins2[i]) * (piKmoins1[i] - piKmoins2[i]);
		h = piK[i] - 2 * piKmoins1[i] + piKmoins2[i];
		
		if(h != 0 && piK[i] > (g/h)) 
			piK[i] = piK[i] - (g/h);
	}
}

int aitken(Matrice *M, double *piK, double *piKmoins1, double *piKmoins2) {
	
	// k = nombre d'itérations (comme dans piK)
	int k = 0;
	int i, j;
	
	printf("# Iteration - Norme\n");
	double somme, norme, sigma;
	do {
		norme = 0.0;
		sigma = calcul_sigma(M , piK);
		
		//ancien pi prend nouveau pi 
		for (i = 0; i < M->n; i++) {
			piKmoins2[i] = piKmoins1[i];
			piKmoins1[i] = piK[i];
			piK[i] = 0.0;
		}
		
		for (i = 0; i < M->n; i++) {
			
			//multiplication ligne colonne
			for (j = M->debCol[i]; j < M->debCol[i + 1]; j++) {
				piK[i] += piKmoins1[M->T[j].i] * M->T[j].p;
			}
			piK[i] = ALPHA * piK[i] + (1 - ALPHA + sigma * ALPHA) / M->n;
			norme += fabs(piK[i] - piKmoins1[i]);
		}
		//Périodiquement, on utilise une estimation de pi
		if(k == 10) estimation(M->n, piK, piKmoins1, piKmoins2);
		
		k++;		
		printf("%d %f\n", k, log(norme));		
	}while(norme > 10e-6);
 	
 	return k;
}
