#include "puissances.h"

int puissances(Matrice *M, double *piK, double *piKmoins1) {
	
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
		
		//on incrémente l'itération
		k++;
		printf("%d %f\n", k, log(norme));
	}while(norme > 10e-6);
 	
 	return k;
}
