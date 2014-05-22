/*
* lu.c
* Par Guillaume Lahaie et Guy Francoeur
*
* Librairie pour la decomposition LU
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "lu.h"

#define MAXRAND (32767.0)
#define EPSILON 0.001

void testMatrix(int n, int in, double *a){

  if (in) { // in != 0 est vrai (matrice source sinon resultat)
    switch(n) {
      case 1: {
        double a1[] = {3.00, -7.00, -2.00, 2.00, -3.00, 5.00, 1.00, 0.00, 6.00, -4.00, 0.00, -5.00, -9.00, 5.00, -5.00, 12.00};
        a = (double*)malloc(16*sizeof(double));
        a = a1;
        break;
      }
      case 2: {
        double a2[] = {-5.000, 3.000, 3.000, 4.000, 2.000, -4.000, -2.000, 3.000, 2.000, -7.000, -3.000, 9.000, 6.000, -9.000, -5.000, 8.000};
        a = (double*)malloc(16*sizeof(double));
        a = a2;
        break;
      }
    }
  } else {
    switch(n) {
      case 1: {
        double s1[] = {3.000, -7.000, -2.000, 2.000, -1.000, -2.000, -1.000, 2.000, 2.000, -5.000, -1.000, 1.000, -3.000, 8.000, 3.000, -1.000};
        a = (double*)malloc(16*sizeof(double));
        a = s1;
        break;
      }
      case 2: {
        double s2[] = {-5.000, 3.000, 3.000, 4.000,-0.400, -2.800, -0.800, 4.600,-0.400, 2.071, -0.143, 1.071,-1.200, 1.929, -1.000, 5.000};
        a = (double*)malloc(16*sizeof(double));
        a = s2;
        break;
      }
    }
  }
}

void initializeMatrix(int n, int m, double *a) {
    srand(time(NULL));
    for (int i=0; i<m; i++ ) {
        for (int j=0; j<n; j++ ) {
            a[i*n+j] = (double) rand() / MAXRAND;
        }
    }
}

void printMatrix(int n, int m, double *a) {

    for (int i=0; i<m; ++i ) {
        for (int j=0; j<n; ++j ) {
            printf("%8.2f ", a[i*n+j]);
        }
        printf("\n");
    }
    fflush(stdout);
}


void writeMatrix(int n, int m, double *a, char *s) {
  FILE *file;
  file = fopen(s,"w");

  for (int i=0; i<m; ++i ) {
    for (int j=0; j<n; ++j ) {
      fprintf(file, "%8.2f ", a[i*n+j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
}

double *readMatrix(const char *nomFichier, int matrix_size, int nbLines) {
    char *line = malloc(256 * sizeof(char));
    size_t len;
    size_t end;
    char *endNb;
    FILE *file;
    int i = 0;
    int j;
    char *ptr;

    double *matrice = (double*)malloc(matrix_size*matrix_size*sizeof(double));

    file = fopen(nomFichier, "r");

    if(file == NULL) {
        printf("Erreur fichier\n");
        exit(-1);
    }

    //On lit la première ligne:
    while((end = getline(&line, &len, file)) != -1) {
      for(j = 0, ptr = strtok(line, " "); ptr != NULL; j++, ptr = strtok(NULL, " ")) {
        matrice[i*matrix_size+j] = (double)strtod(ptr, &endNb);
      }
      i++;
    }
    fclose(file);

    return matrice;
}

void checkResult( double* resultat, double *a, int n , int nbLines, MPI_Comm comm) {
    int i, j;

    for(i = 0; i < nbLines; i++) {
        for(j = 0; j < n; j++) {
            if(fabs(resultat[i*n+j] - a[i*n+j]) >= EPSILON) {
               fprintf(stderr,"Erreur - resultat de la factorisation incorrect\n");
               MPI_Abort(comm, 1);
            }
        }
    }
}


int trouverMax(double *maxes, int nodes) {
    int i, 
        maxPos = 0;
    double max = maxes[0];
    for(i = 1; i < nodes; i++) {
        if(fabs(maxes[i]) > fabs(maxes[maxPos])) {
            maxPos = i;
        }
    }
    return maxPos;
}





