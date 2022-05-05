#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define alpha 0.85
#define EPSILON 1e-6

int * f;
//pour compiler : gcc aitken.c -o aitken -lm
// Structure représentant un arc du graphe
typedef struct Arc Arc;
struct Arc {
  int ligne;
  int col;
  double coeff;
  Arc * next;
};


// Structure représentant la liste chaînée
typedef struct Liste Liste;
struct Liste {
  Arc ** list;
  int dimension;
};


void cleanListe(Liste *);
void addToList(Arc *, Liste *);
double norme1(double * npi, double * opi, int dimension);



Liste * loadFile(char* filename) {

    FILE * file;
    file = fopen(filename, "r");
    if (file == NULL) {
        perror("Couldn't open the file.\n");
        return NULL;
    }

    int dimension;
    int nbArcs;
    fscanf(file, "%d %d", &dimension, &nbArcs);
    
    Liste * res = (Liste * ) malloc(sizeof(Liste));
    if (res == NULL) {
        perror("Erreur\n");
        return NULL;
    }
    res -> dimension = dimension;
    res -> list = (Arc **) malloc(dimension * sizeof(Arc *));
    if (res -> list == NULL) {
        perror("Couldn't initialize linked list.\n");
        return NULL;
    }

    f = (int *) malloc(dimension * sizeof(int));
    if (f == NULL) {
        perror("Couldn't create f vector.\n");
        return NULL;
    }

    int col;
    double coeff;

    for (int i = 0; i < dimension; i++) {

        int numLigne;
        int degreLigne;
        fscanf(file, "%d %d", &numLigne, &degreLigne);

        if (degreLigne == 0)
            f[i] = 1;
        else
            f[i] = 0;

        for (int j = 0; j < degreLigne; j++) {
            
            
            fscanf(file, "%d %lf", &col, &coeff);

            Arc * arc = (Arc *) malloc(sizeof(Arc));
            if (arc == NULL) {
                perror("Couldn't create the arc.\n");
                return NULL;
            }
            arc -> ligne = numLigne;
            arc -> col = col;
            arc -> coeff = coeff;
            arc -> next = NULL;

            addToList(arc, res);
        }
    }
    
    fclose(file);
    printf("Lecture finie\n");
    return res;
}


void addToList(Arc * arc, Liste * list) {
    
    int col = arc -> col - 1;
    if (list -> list[col] == NULL) {
        list -> list[col] = arc;
    } 
    else {
        Arc * tmp = list -> list[col];
        while (tmp != NULL && tmp -> next != NULL)  //ICI
            tmp = tmp -> next;
        tmp -> next = arc;
    }
}


void cleanAll(Liste * list) {
    
    Arc * ptr;
    int n = list -> dimension;
    
    for (int i = 0; i < n; i++) {
        ptr = list -> list[i];
        while (ptr != NULL) {
            Arc * next = ptr -> next;
            free(ptr);
            ptr = next;
        }
    }
    free(list -> list);
    free(list);
}

double double_abs(double value){
    if (value < 0) return -value;
    else return value;
}


double norme1(double * npi, double * opi, int dimension) {
  
  double norme = 0.;
  
  for (int i = 0; i < dimension; i++)
    norme += double_abs(npi[i] - opi[i]);
  
  return norme;
}


double * product(double * pi, Arc ** P, int dimension) {

    double * res = (double *) malloc(dimension * sizeof(double));
    if (res == NULL) {
        perror("Error during product.\n");
        return NULL;
    }
    
    for (int i = 0; i < dimension; i++) {
        
        Arc * ptr = P[i];
        double sum = 0;
        
        for (int j = 0; j < dimension; j++) {
            if (ptr == NULL || ptr -> ligne != j + 1)
                continue;
            sum += pi[j] * (ptr -> coeff);
            ptr = ptr -> next;
        }

        res[i] = sum;
    }

    return res;
}


double * puissance(Liste * P, double epsilon) {

    int dimension = P -> dimension;
    double const_term = (1. - alpha) * (1. / dimension);

    double * opi = (double *) malloc(dimension * sizeof(double));
    if (opi == NULL) {
        perror("Couldn't create opi.\n");
        return NULL;
    }

    for (int i = 0; i < dimension; i++)
        opi[i] = 1. / dimension;
    
    double * npi = product(opi, P -> list, dimension);
    if (npi == NULL) {
        free(opi);
        perror("Error. Terminating program.\n");
        return NULL;
    }
    

    double sum = 0;
    for (int j = 0; j < dimension; j++)
        sum += opi[j] * f[j];
    sum *= alpha * (1. / dimension);
    for (int i = 0; i < dimension; i++) {
        npi[i] *= alpha;
        npi[i] += const_term;
        npi[i] += sum;
    }
    printf("\n npi : \n");
         double somme = 0;
    for (int i = 0; i < dimension; ++i)
    {
        somme+= npi[i];
        //printf(" %lf ", npi[i] );
    }
    printf(" sum : %lf\n", somme);
    

    

    int i = 2;

    while (norme1(npi, opi, dimension) >= epsilon) {
        
        free(opi);
        opi = npi;

        npi = product(opi, P -> list, dimension);
        if (npi == NULL) {
            free(opi);
            perror("Error. Terminating program.\n");
            return NULL;
        }

        double sum = 0;
        for (int j = 0; j < dimension; j++)
            sum += opi[j] * f[j];
        sum *= alpha * (1. / dimension);
        for (int i = 0; i < dimension; i++) {
            npi[i] *= alpha;
            npi[i] += const_term;
            npi[i] += sum;
        }
        
         printf("npi : \n");
         double somme = 0;
    for (int i = 0; i < dimension; ++i)
    {
        somme+= npi[i];
        //printf(" %lf ", npi[i] );
    }
    printf(" \nsum : %lf\n", somme);
        printf("Itération %d\n", i);
        i++;
    }
    
    free(opi);
    return npi;
}

double *Aitken(double * pi,double * pi1,double * pi2, int dimension){
    printf("acceleration :\n");
    double nomi = 0;
    double deno = 0;
    double * res = (double *) malloc(dimension * sizeof(double));
    for (int i = 0; i < dimension; ++i)
    {
        nomi = (pi1[i] - pi[i]) * (pi1[i] - pi[i]);
        deno = (pi2[i] - (2 * pi1[i]) + pi[i]);
        if(deno!=0 && pi2[i] > (nomi/deno) )
            res[i] = pi2[i] - (nomi/deno);
        else
            return pi2;
        
        //printf(" %lf ", res[i] );
    }
    return res;
}

void aitken_puissance(Liste * P, double epsilon) {

   int dimension = P -> dimension;
    double const_term = (1. - alpha) * (1. / dimension);

    double * opi = (double *) malloc(dimension * sizeof(double));
    if (opi == NULL) {
        perror("Couldn't create opi.\n");
        exit(0);
    }
    for (int i = 0; i < dimension; i++)
        opi[i] = 1. / dimension;
    
    double * npi = product(opi, P -> list, dimension);
    if (npi == NULL) {
        free(opi);
        perror("Error. Terminating program.\n");
        exit(0);
    }
    double sum = 0;
    for (int j = 0; j < dimension; j++)
        sum += opi[j] * f[j];
    sum *= alpha * (1. / dimension);
    for (int i = 0; i < dimension; i++) {
        npi[i] *= alpha;
        npi[i] += const_term;
        npi[i] += sum;
    }
    printf("\n npi : \n");
         double somme = 0;
    for (int i = 0; i < dimension; ++i)
    {
        somme+= npi[i];
        //printf(" %lf ", npi[i] );
    }

    printf(" sum : %lf\n", somme);

    double * pi2 = product(npi, P -> list, dimension);
    if (pi2 == NULL) {
        free(npi);
        perror("Error. Terminating program.\n");
        exit(0);
    }
    sum = 0;
    for (int j = 0; j < dimension; j++)
        sum += npi[j] * f[j];
    sum *= alpha * (1. / dimension);
    for (int i = 0; i < dimension; i++) {
        pi2[i] *= alpha;
        pi2[i] += const_term;
        pi2[i] += sum;
    }
    int i = 1;
    int period = 40;
    double delta =1;
    while (delta > epsilon) {
        
        //free(opi);
        opi = npi;
        npi= pi2;
        pi2 = product(npi, P -> list, dimension);
    if (pi2 == NULL) {
        free(npi);
        perror("Error. Terminating program.\n");
        exit(0);
    }
    sum = 0;
    for (int j = 0; j < dimension; j++)
        sum += npi[j] * f[j];
    sum *= alpha * (1. / dimension);
    for (int i = 0; i < dimension; i++) {
        pi2[i] *= alpha;
        pi2[i] += const_term;
        pi2[i] += sum;
    }
        
         printf("pi2 : \n");
         double somme = 0;
    for (int i = 0; i < dimension; ++i)
    {
        somme+= pi2[i];
        //printf(" %lf ", npi[i] );
    }
    delta = norme1(npi, pi2, dimension);
    printf(" \nsum : %lf\n", somme);
        printf("Itération %d\n", i);
        
        if(i%period ==0)
            pi2 = Aitken(opi,npi,pi2,dimension);
        i++;
        
    }
    
    free(opi);
    free(npi);
    free(pi2);
}

//pour compiler : gcc aitken.c -o aitken -lm
int main() {
    Liste * list = loadFile("data/wb-cs-stanford.txt");
    //puissance(list, EPSILON);
    aitken_puissance(list,EPSILON);
    cleanAll(list);
    free(f);
    return 0;
}
