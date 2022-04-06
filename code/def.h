#ifndef DEF_H
#define DEF_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define ALPHA 0.8

typedef struct {
	//coordonees
    int i;
	int j;
    //valeur
	double p;
}Element;

typedef struct {
	int m;
	//tableau d'elements
    Element *T; 
	//taille de la matrice carre
    int n;
	int *debCol;
}Matrice;

#endif