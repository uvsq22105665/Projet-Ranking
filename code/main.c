#include "def.h"
#include "matrice.h"
#include "puissances.h"
#include "aitkenprim.h"

int main(int argc, char *argv[]) {
	
	Matrice *M;
	time_t debut, fin;
		
	FILE *web = NULL;
	web = fopen(argv[2], "r");
	
	if(web != NULL) {
		
		//initialisation de la matrice
		M = malloc(sizeof(Matrice));
		import_matrice(web, M); fclose(web);
		
		double *piK = malloc(M->n * sizeof(double));
		double *piKmoins1 = malloc(M->n * sizeof(double));
		double *piKmoins2;
		
		if(strcmp(argv[1], "puissances") == 0) init_distrib(M->n, piK, piKmoins1);
		else if(strcmp(argv[1], "aitken") == 0) {
			piKmoins2 = malloc(M->n * sizeof(double));
			init_distrib_aitken(M->n, piK, piKmoins1, piKmoins2);
		}
		
		
		//calcul des pertinences
		time(&debut);
		
		int nbIterations = 0;
		if(strcmp(argv[1], "puissances") == 0) nbIterations = puissances(M, piK, piKmoins1);
		else if(strcmp(argv[1], "aitken") == 0) nbIterations = aitken(M, piK, piKmoins1, piKmoins2);
		else printf("Erreur arguments.\n");
		
		time(&fin);
		printf("# nbIterations - secondes\n");
		printf("# %d %f\n", nbIterations, difftime(fin, debut));

		free_matrice(M); free(M);
		free(piK); free(piKmoins1);
		if(strcmp(argv[1], "aitken") == 0) free(piKmoins2);
	}
	else
		printf("Erreur arguments.\n");
	return 0;
	
	return 0;
}
