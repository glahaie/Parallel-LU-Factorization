/*
 * lu_mpi_row_cyclic.c
 * Par Guillaume Lahaie et Guy Francoeur
 *
 * Version de lu_mpi où chaque noeud n'a que sa partie de la matrice
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include "lu.h"

typedef struct Arguments {
    int matrix_size;
    char *matrice;
    char * solution;
    int imprimer;
    int pivot;
    char *write;
} Arguments;

//Effectue la factorisation de la matrice a, avec ou sans pivot
void factoriserLU(double *a, int size, int nbLines, MPI_Comm comm, int rank, 
        int pivot, int nodes);

//Lecture d'un fichier de matrices pour un test précis
double * readMatrixCyclic(char *nomFichier, int matrix_size, int nbLines, 
        int rank, int nodes);

//Pivoter les lignes, si on utilise le pivot partiel
void effectuerPivot(double *a, MPI_Comm comm, int rank, int currentLine, 
        int size, int nbLines, int nodes);

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

    int nbLines = arguments.matrix_size / nodes;

    if(arguments.matrice) {
        a = readMatrixCyclic(arguments.matrice, arguments.matrix_size, 
                nbLines, myrank, nodes);
    } else {  
        a = (double *) malloc( nbLines*arguments.matrix_size*sizeof(double) );
        initializeMatrix( arguments.matrix_size, nbLines, a);
    }
    double tempsEcoule;
    MPI_Barrier( MPI_COMM_WORLD );
    tempsEcoule = -MPI_Wtime();

    factoriserLU(a, arguments.matrix_size, nbLines, MPI_COMM_WORLD, myrank,
            arguments.pivot, nodes);

    double tempsMax;
    tempsEcoule += MPI_Wtime();

    //On ramene le temps max au processeur 0
    MPI_Reduce(&tempsEcoule, &tempsMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(arguments.solution) {
        //on vérifie le résultat
        double *resultat = readMatrixCyclic(arguments.solution, 
                arguments.matrix_size, nbLines, myrank, nodes);
        checkResult(resultat, a, arguments.matrix_size, nbLines, MPI_COMM_WORLD);
        if(myrank == 0) {
            printf("Resultat de la factorisation lu correcte.\n");
        }

        if(resultat) free(resultat);
    }

    if (myrank == 0) {
        printf( "Factorisation LU %s pivot effectuee en %.3f msecs\n", 
                (arguments.pivot?"avec":"sans"), 1000*tempsEcoule);
    }

    if(arguments.imprimer || arguments.write) {

        //On test pour récupérer la matrice
        double *matriceComplete;
        if(myrank == 0) {
            matriceComplete = (double *) malloc(arguments.matrix_size*
                    arguments.matrix_size*sizeof(double));
        }
        for( i = 0; i < nbLines; i++) {
            MPI_Gather(&a[i*arguments.matrix_size], arguments.matrix_size, MPI_DOUBLE, 
                    &matriceComplete[i*nodes*arguments.matrix_size], 
                    arguments.matrix_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        if( myrank == 0) {
            if(arguments.imprimer) {
                printMatrix(arguments.matrix_size, arguments.matrix_size, 
                        matriceComplete);
            }
            if(arguments.write) {
                writeMatrix(arguments.matrix_size, arguments.matrix_size, 
                        matriceComplete, arguments.write);
            }
            if(matriceComplete) free(matriceComplete);
        }


    }

#ifndef NDEBUG
    printf("rank G:%i\n",myrank);
#endif

    if(a)free(a);
    MPI_Finalize();
    return 0;
}

double* readMatrixCyclic(char *nomFichier, int matrix_size, int nbLines, int rank, 
        int nodes) 
{

    char *line = malloc(256);
    size_t len;
    ssize_t end;
    char *endNb;
    FILE *file;
    int i = 0;
    int j;
    char *ptr;

    double *matrice = readMatrix(nomFichier, matrix_size, matrix_size);

    double *bloc = (double*)malloc(nbLines*matrix_size*sizeof(double));

    //On prend la partie qu'on a besoin
    for(i = 0; i < nbLines; i++) {
        for(j = 0; j < matrix_size; j++) {
            bloc[i*matrix_size + j] = matrice[(i*nodes+rank)*matrix_size + j];
        }
    }

    if(matrice) free(matrice);
    return bloc;

}

void factoriserLU(double *a, int size, int nbLines, MPI_Comm comm, int rank, 
        int pivot, int nodes) {

    double *c = (double*)malloc(size * sizeof(double));

    int i, j, k;
    for (j=0; j<size-1; j++) {

        if(pivot) {
            effectuerPivot( a, comm, rank, j, size, nbLines, nodes);
        }

        if( j % nodes == rank) {
            //C'est moi qui a les valeurs
            int line = j / nodes;
            memcpy(c, &a[line*size + j], (size-j)*sizeof(double));
        }
        MPI_Bcast(c,size-j,MPI_DOUBLE,(j % nodes),MPI_COMM_WORLD);

        //Maintenant on met à jour nos lignes
        for (i= 0; i< nbLines; i++) {
            if (i*nodes+rank > j) {
                double mult = a[i*size+j] / c[0];
                a[i*size+j] = mult;
                for(k = 1; k < (size-j); k++) {
                    a[i*size+j+k] = a[i*size+j+k] - mult*c[k];
                }
            }
        }
    }
    if(c) free(c);
}

void effectuerPivot(double *a, MPI_Comm comm, int rank, int currentLine, 
        int size, int nbLines, int nodes) {
 
    int i, posMaxLocal;
    double maxLocal = 0;

    double *maxes = (double*)malloc(nodes*sizeof(double));
    double *c = (double*)malloc(size*sizeof(double));

      for(i = 0; i < nbLines; i++) {
          if((i*nodes +rank >= currentLine)) {
              if (fabs(a[i*size+currentLine]) > maxLocal) {
                  maxLocal = fabs(a[i*size+currentLine]);
                  posMaxLocal = i;
              }
          }
      }

      //Maintenant on peut trouver le max
      MPI_Allgather(&maxLocal, 1, MPI_DOUBLE, maxes, 1, MPI_DOUBLE, MPI_COMM_WORLD);

      int noeudMax = trouverMax(maxes, nodes);

      //Différent traitement ici
      if (noeudMax == rank && currentLine % nodes != rank) {
          //Je dois envoyer ma ligne
          memcpy(c, &a[posMaxLocal*size], size*sizeof(double));
          MPI_Send(c, size, MPI_DOUBLE, currentLine % nodes, 0, MPI_COMM_WORLD);

          MPI_Status status;
          //Maintenant je recois
          MPI_Recv(&a[posMaxLocal*size], size, MPI_DOUBLE, 
                  currentLine % nodes, 0, MPI_COMM_WORLD, &status);
      } else if(noeudMax != rank && (currentLine % nodes) == rank) {
          //On recoit dans c
          memcpy(c, &a[(currentLine/nodes)*size], size*sizeof(double));
          MPI_Send(c, size, MPI_DOUBLE, noeudMax, 0, MPI_COMM_WORLD);
          MPI_Status status;
          MPI_Recv(&a[(currentLine/nodes)*size], size, MPI_DOUBLE, 
                  noeudMax, 0, MPI_COMM_WORLD, &status);

      } else if (noeudMax == rank && (currentLine%nodes) == rank) {
          
          if(posMaxLocal != (currentLine/nodes)) {
              memcpy(c, &a[posMaxLocal*size], size*sizeof(double));
              memcpy(&a[posMaxLocal*size], 
                      &a[(currentLine/nodes)*size], size*sizeof(double));
              memcpy(&a[(currentLine/nodes)*size], c, size*sizeof(double));
          }
      }
      if(maxes) free(maxes);
      if(c) free(c);
}

