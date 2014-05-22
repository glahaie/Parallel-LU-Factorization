#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "lu.h"
#define EPSILON 0.001

typedef struct Arguments {
    int matrix_size;
    char *matrice;
    char * solution;
    int imprimer;
    int pivot;
    char * write;
} Arguments;

//Effectue la factorisation, avec ou sans pivot
void factoriserLU(double *a, int pivot, int size);

//Vérifie si on obtient bien le bon résultat pour un test
void checkResultSeq( double* resultat, double *a, int n, int m );

/*---------------------------------------------------------------------------*/
int main( int argc, char *argv[] ) {
  double *a;
  int matrix_size;
  int i, j, k, noTest;
  Arguments arguments;

  arguments.matrix_size = -1;
  arguments.matrice = arguments.solution = NULL;
  arguments.imprimer = 0;
  arguments.pivot = 0;
  arguments.write = NULL;

  for(i = 1; i < argc; i++) {
      if(argv[i][0] == '-') {
          switch(argv[i][1]) {
              case 'm': if(i+1 < argc) {
                            arguments.matrix_size = atoi(argv[++i]);
                        }
                        break;
              case 'f': if(i+1 < argc) {
                            arguments.matrice = argv[++i];
                        }
                        break;
              case 's': if(i+1 < argc) {
                            arguments.solution = argv[++i];
                        }
                        break;
              case 'p': arguments.pivot = 1;
                        break;
              case 'i': arguments.imprimer = 1;
                        break; 
              case 'w': if(i+1 < argc) {
                            arguments.write = argv[++i];
                        }
                        break;
              default: fprintf(stderr, "Usage: %s -m <tailleMatrice> [-f nomFichier] " \
                               "[-s solution]\n", argv[0]);
                       exit(EXIT_FAILURE);
          }
      }
  }

  if(arguments.matrix_size < 0) {
      fprintf(stderr,"Erreur: vous devez fournir une taille de matrice.\n");
      exit(EXIT_FAILURE);
  }

  if(arguments.matrice) {
      a = readMatrix(arguments.matrice, arguments.matrix_size, arguments.matrix_size);
  } else {  
    a = (double *) malloc( arguments.matrix_size*arguments.matrix_size*sizeof(double) );
    initializeMatrix( arguments.matrix_size, arguments.matrix_size,a);
  }

  clock_t debut = clock();

  factoriserLU(a, arguments.pivot, arguments.matrix_size);


  clock_t fin = clock();
  double  tempsEcoule = ((double) fin - debut) / CLOCKS_PER_SEC;
  printf( "Factorisation LU %s pivot effectuee en %.3f msecs\n", 
              (arguments.pivot?"avec":"sans"), 1000*tempsEcoule);

  if(arguments.solution) {
      //on vérifie le résultat
      double *resultat = readMatrix(arguments.solution, arguments.matrix_size,
              arguments.matrix_size);
      checkResultSeq(resultat, a, arguments.matrix_size, arguments.matrix_size);
      printf("Résultat de la factorisation LU correct\n");
      if(resultat) free(resultat);
  }

  if(arguments.imprimer) {
      printMatrix(arguments.matrix_size, arguments.matrix_size, a);
  }

  if(arguments.write) {
      writeMatrix(arguments.matrix_size, arguments.matrix_size, a, arguments.write);
  }

  if(a) free( a );
  return 0;
}

void factoriserLU(double *a, int pivot, int size) {
    
    int i, j, k, posPivot;
    double mult;
    double *temp = (double*)malloc(size*sizeof(double));


    for ( i=0; i< size-1; i++ ) {

      if(pivot) {
          posPivot = i;
          for(j = i+1; j < size; j++) {
              if(fabs(a[posPivot*size + i]) < fabs(a[j*size + i])) {
                  posPivot = j;
              }
          }
          //On échange les lignes si nécessaire
          if(posPivot != i) {
              memcpy(temp, &a[posPivot*size], size*sizeof(double));
              memcpy(&a[posPivot*size], &a[i*size], size*sizeof(double));
              memcpy(&a[i*size], temp, size*sizeof(double));
          }
      }

      //On fait les réductions
      for(j = i+1; j < size; j++) {
          mult = a[j*size+i] / a[i*size+i];
          a[j*size+i] = mult;
          for(k = i+1; k < size; k++) {
              a[j*size+k] = a[j*size+k] - mult * a[i*size + k];
          }
      }
    }
    if(temp) free(temp);
}

void checkResultSeq( double* resultat, double *a, int n, int m ) {
  int i, j;

  for(i = 0; i < m; i++) {
      for(j = 0; j < n; j++) {
          assert((fabs(resultat[i*n+j] - a[i*n+j]) <= EPSILON) &&
                  "Erreur - resultat de la factorisation incorrect.");
      }
  }
}

