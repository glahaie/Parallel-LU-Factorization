/*
 * lu_mpi_row_adjacent.c
 * Par Guillaume Lahaie et Guy Francoeur
 *
 * Version de lu_mpi où chaque noeud n'a que sa partie de la matrice.
 * Les valeurs sous la matrice triangulaire superieure contiennent les
 * valeurs de la matrice L (la diagonale de L étant toujours des 1, on
 * peut facilement reconstruire, et la matrice triangulaire supérieure
 * représente la matrice U.
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include "lu.h"

typedef struct Arguments {
    int matrix_size;
    char *matrice;
    char * solution;
    int imprimer;
    int pivot;
    char* write;
} Arguments;

//Lecture d'un fichier de matrices pour un test précis
double * readMatrixAdjacent(const char *nomFichier, int matrix_size, int NbLines, int rank);

//Effectue la factorisation de la matrice a, avec ou sans pivot
void factoriserLU(double *a, int size, int nbLines, MPI_Comm comm, int rank, 
        int pivot, int nodes);


int main(int argc, char* argv[]) {
    int myrank;
    int nodes;
    double *a;

    int i, j, k, noTest;
    Arguments arguments;

    arguments.matrix_size = -1;
    arguments.matrice = arguments.solution = NULL;
    arguments.imprimer = 0;
    arguments.pivot = 0;
    arguments.write = NULL;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nodes);

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
                default: fprintf(stderr, "Usage: %s -m <tailleMatrice> "\
                                 "[-f nomFichier] [-s solution]\n", argv[0]);
                         exit(EXIT_FAILURE);
            }
        }
    }

    if(arguments.matrix_size < 0) {
        fprintf(stderr,"Erreur: vous devez fournir une taille de matrice.\n");
        exit(EXIT_FAILURE);
    }

    if(arguments.matrix_size % nodes != 0) {
        fprintf(stderr, "Erreur: le nombre de noeuds doit diviser la"\
                "taille de la matrice.\n");
        exit(EXIT_FAILURE);
    }

    MPI_Barrier( MPI_COMM_WORLD );
    double tempsEcoule;
    tempsEcoule = -MPI_Wtime();

    int nbLines = arguments.matrix_size / nodes;

    if(arguments.matrice) {
        a = readMatrixAdjacent(arguments.matrice, arguments.matrix_size, nbLines, myrank);
    } else {  
        a = (double *) malloc( nbLines*arguments.matrix_size*sizeof(double) );
        initializeMatrix( arguments.matrix_size, nbLines, a);
    }

    //Factorisation
    factoriserLU(a, arguments.matrix_size, nbLines, MPI_COMM_WORLD, myrank,
            arguments.pivot, nodes);

    //On ramene le temps max au processeur 0
    double tempsMax;
    tempsEcoule += MPI_Wtime();
    MPI_Reduce(&tempsEcoule, &tempsMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(arguments.solution) {
        //on vérifie le résultat
        double *resultat = readMatrixAdjacent(arguments.solution, arguments.matrix_size, 
                nbLines, myrank);
        checkResult(a, resultat, arguments.matrix_size, nbLines, MPI_COMM_WORLD);

        if(myrank == 0) {
            printf("Resultat de la factorisation lu correcte.\n");
        }
        if(resultat) free(resultat);
    }

    if (myrank == 0) {
        printf( "Factorisation LU %s pivot effectuee en %.3f msecs\n", 
                (arguments.pivot?"avec":"sans"),1000*tempsEcoule);
    }

    if(arguments.imprimer || arguments.write) {
        //on reconstruit la matrice
        double *matrice;
        if (myrank == 0) {
            matrice = (double *) malloc(arguments.matrix_size*
                    arguments.matrix_size*sizeof(double));
        }

        MPI_Gather(a, nbLines*arguments.matrix_size, MPI_DOUBLE, 
                matrice, nbLines*arguments.matrix_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(myrank == 0) {
            if(arguments.imprimer) {
                printMatrix(arguments.matrix_size, arguments.matrix_size, matrice);
            }
            if(arguments.write) {
                writeMatrix(arguments.matrix_size, arguments.matrix_size,
                        matrice, arguments.write);
            }
            if(matrice) free(matrice);
        }
    }

#ifndef NDEBUG
    printf("rank G:%i\n",myrank);
#endif
    free(a);
    MPI_Finalize();
    return 0;
}


double* readMatrixAdjacent(const char *nomFichier, int matrix_size, int nbLines, int rank) {

    char *line = malloc(256);
    size_t len;
    ssize_t end;
    char *endNb;
    FILE *file;
    int i = 0;
    int j;
    char *ptr;

    double *matrice = readMatrix(nomFichier, matrix_size, nbLines);

    double *bloc = (double*)malloc(nbLines*matrix_size*sizeof(double));

    //On prend la partie qu'on a besoin
    for(i = 0; i < nbLines; i++) {
        for(j = 0; j < matrix_size; j++) {
            bloc[i*matrix_size + j] = matrice[(rank*nbLines+i)*matrix_size + j];
        }
    }

    free(matrice);
    return bloc;

}

void effectuerPivot(double *a, MPI_Comm comm, int rank, int currentLine, 
        int size, int nbLines, int nodes) {
      
    int i;
    double  maxLocal = 0;
    int posMaxLocal;

    double *line = (double*)malloc(size*sizeof(double));
    double *maxes = (double*)malloc(nodes*sizeof(double));

      for(i = 0; i < nbLines; i++) {
          if((rank*nbLines +i >= currentLine)) {
              if (fabs(a[i*size+currentLine]) > maxLocal) {
                  maxLocal = fabs(a[i*size+currentLine]);
                  posMaxLocal = i;
              }
          }
      }

      //Maintenant on peut trouver le max
      MPI_Allgather(&maxLocal, 1, MPI_DOUBLE, maxes, 1, MPI_DOUBLE, comm);

      int noeudMax = trouverMax(maxes, nodes);

      //Différent traitement ici
      if (noeudMax == rank && currentLine / nbLines != rank) {
          //Je dois envoyer ma ligne
          memcpy(line, &a[posMaxLocal*size], size*sizeof(double));
          MPI_Send(line, size, MPI_DOUBLE, currentLine / nbLines, 0, comm);

          MPI_Status status;
          //Maintenant je recois
          MPI_Recv(&a[posMaxLocal*size], size, MPI_DOUBLE, 
                  currentLine / nbLines, 0, comm, &status);
      } else if(noeudMax != rank && (currentLine / nbLines) == rank) {
          //On recoit dans c
          memcpy(line, &a[(currentLine%nbLines)*size], size*sizeof(double));

          MPI_Send(line, size, MPI_DOUBLE, noeudMax, 0, MPI_COMM_WORLD);
          MPI_Status status;
          MPI_Recv(&a[(currentLine%nbLines)*size], size, 
                  MPI_DOUBLE, noeudMax, 0, comm, &status);

      } else if (noeudMax == rank && (currentLine/nbLines) == rank) {
          //Pas besoin de communication -- on pourrait vérifier 
          //si un swap est vraiment nécessaire
          if(posMaxLocal != (currentLine%nbLines)) {
              memcpy(line, &a[posMaxLocal*size], size*sizeof(double));
              memcpy(&a[posMaxLocal*size], &a[(currentLine%nbLines)*size], 
                      size*sizeof(double));
              memcpy(&a[(currentLine%nbLines)*size], line, size*sizeof(double));
          }
      }
      free(line);
      free(maxes);

}


void factoriserLU(double *a, int size, int nbLines, MPI_Comm comm, int rank, 
        int pivot, int nodes) {

    double *c = (double*)malloc(size*sizeof(double));
    int i, j, k;

    for (j=0; j<size-1; j++) {

      if(pivot) {
          effectuerPivot(a, comm, rank, j, size, nbLines, nodes);
      }


      if( j / nbLines == rank) {
          //C'est moi qui a les valeurs
          memcpy(c, &a[(j % nbLines)*size + j], (size-j)*sizeof(double));
      }
      MPI_Bcast(c,size-j,MPI_DOUBLE,(j / nbLines),comm);


      //Maintenant on met à jour nos lignes
      for (i=0; i<nbLines; i++) {
          if (rank*nbLines+i > j) {
              double mult = a[i*size+j] / c[0];
              //On sauve la valeur pour L
              a[i*size + j] = mult;
              for(k = 1; k < (size-j); k++) {
                  a[i*size+j+k] = a[i*size+j+k] - mult*c[k];
              }
          }
      }
    }
    if(c) free(c);
}


