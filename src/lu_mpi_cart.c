/*
 * lu_mpi_cart.c
 *
 * Par Guillaume Lahaie et Guy Francoeur
 *
 * Implémentation de la décomposition LU par blocks pour
 * une matrices N * N. J'utlise le créateur de topologie
 * cartésienne de MPI pour gérér la répartition. Il faut fournir
 * à la ligne de commande la répartition de la matrice, donc le nombre
 * de blocs en x et y, et on doit avoir x*y = nbProcs
 */


#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
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
double* readMatrixCart(char *nomFichier, int matrix_size, int nbCols, int nbLines, 
        int coords[]);

void factoriserLU(double *a, int size, int nbLinesParBloc, int nbColsParBloc,
        MPI_Comm lignes, MPI_Comm colonnes, int coord[]);

int main(int argc, char* argv[]) {

    int myrank;
    int nodes;
    double *a;
    int nbDims = 2; 
    int dims[2];


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
                case 'x': if(i+1 < argc) {
                             dims[0] = atoi(argv[++i]);
                          }
                         break; 
                case 'y': if(i+1 < argc) {
                             dims[1] = atoi(argv[++i]);
                          }
                         break; 
                case 'w': if(i+1 < argc) {
                              arguments.write = argv[++i];
                          }
                          break;
                default: fprintf(stderr, "Usage: %s -m <tailleMatrice> "\
                                 "-x <lignes> -y <colonnes> [-w nomFichier]" \
                                 "[-f nomFichier] [-s solution]\n", argv[0]);
                         exit(EXIT_FAILURE);
            }
        }
    }

    if(arguments.matrix_size < 0) {
        fprintf(stderr,"Erreur: vous devez fournir une taille de matrice.\n");
        exit(EXIT_FAILURE);
    }

    if((arguments.matrix_size % nodes != 0) || (nodes % dims[0] !=0)
            || (nodes % dims[1] !=0)) {
        fprintf(stderr, "Erreur: le nombre de noeuds doit diviser la"\
                "et les colonnes et lignes doivent diviser le nombre de noeuds .\n");
        exit(EXIT_FAILURE);
    }

    //On part la minuterie
    double tempsEcoule = -MPI_Wtime(); 

    //On crée la topologie
    int period[] = {1,0};
    int reorder = 1;
    int coord[2]; //Pour ma position
    MPI_Comm comm;

    MPI_Cart_create(MPI_COMM_WORLD, nbDims, dims, period, reorder, &comm);

    //On imprime 
    MPI_Cart_coords(comm, myrank, nbDims, coord);

    //Maintenant on crée les matrices
    int nbLinesParBloc = arguments.matrix_size / dims[0];
    int nbColsParBloc = arguments.matrix_size / dims[1];

    if(arguments.matrice) {
        a = readMatrixCart(arguments.matrice, arguments.matrix_size, 
                nbColsParBloc, nbLinesParBloc, coord);
    } else {
        a = (double*)malloc(nbLinesParBloc*nbColsParBloc*sizeof(double));
        initializeMatrix(nbColsParBloc, nbLinesParBloc,a);
    }

    //On crée les communicateurs pour les colonnes et les lignes
    MPI_Comm lignes;
    MPI_Comm colonnes;
    int cols[2] = {1,0};
    int lines[2] = {0,1};

    MPI_Cart_sub(comm, lines, &lignes);
    MPI_Cart_sub(comm, cols,  &colonnes);

    factoriserLU(a, arguments.matrix_size, nbLinesParBloc, nbColsParBloc,
            lignes, colonnes, coord);

    //On fait la décomposition
    double tempsMax;
    tempsEcoule += MPI_Wtime();
    MPI_Reduce( &tempsEcoule, &tempsMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );


    if(arguments.solution) {
        //on vérifie le résultat
        double *resultat =  readMatrixCart(arguments.solution, arguments.matrix_size, 
                nbColsParBloc, nbLinesParBloc, coord);
        checkResult(resultat, a, nbColsParBloc, nbLinesParBloc, MPI_COMM_WORLD);

        if(myrank==0) {
            printf("Resultat de la factorisation lu correcte.\n");
        }
        if(resultat) free(resultat);
    }

    if (myrank == 0) {
        printf( "Factorisation LU sans pivot effectuee en %.3f msecs\n", 
                1000*tempsEcoule);
    }

    if(arguments.imprimer || arguments.write) {

        //Si mon rang est 0 pour mon communicateur Lignes
        double *lines;
        if(coord[1] == 0) {
            lines = (double *)malloc(nbLinesParBloc*arguments.matrix_size*sizeof(double));
        }

        //On rebatit les lignes
        for(i = 0; i < nbLinesParBloc; i++) {
            MPI_Gather(&a[i*nbColsParBloc], nbColsParBloc, MPI_DOUBLE, 
                    &lines[i*arguments.matrix_size], nbColsParBloc, MPI_DOUBLE, 0, lignes);
        }

        double *matriceComplete;
        if(myrank ==0) { //Je suis la racine, rank ==0 aussi
            matriceComplete = (double *)malloc(arguments.matrix_size*
                    arguments.matrix_size*sizeof(double));
        }
        //On fait la communication
        if(coord[1] ==0){ //Je suis dans la colonne 0 -- à vérifier
            MPI_Gather(lines, arguments.matrix_size*nbLinesParBloc, 
                    MPI_DOUBLE, matriceComplete, arguments.matrix_size*
                    nbLinesParBloc, MPI_DOUBLE, 0, colonnes);
            free(lines);
        }

        if(myrank == 0) {
            if(arguments.imprimer) {
                printMatrix(arguments.matrix_size, arguments.matrix_size, 
                        matriceComplete);
            }
            if(arguments.write) {
                writeMatrix(arguments.matrix_size, arguments.matrix_size,
                        matriceComplete, arguments.write);
            }
            free(matriceComplete);
        }
    }

    MPI_Finalize();
    if(a) free(a);
    return 0;
}

double* readMatrixCart(char *nomFichier, int matrix_size, int nbCols, int nbLines, 
        int coords[]) 
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

    double *bloc = (double*)malloc(nbLines*nbCols*sizeof(double));

    //On prend la partie qu'on a besoin
    for(i = 0; i < nbLines; i++) {
        for(j = 0; j < nbCols; j++) {
            int pos = (coords[0]*nbLines+i)*matrix_size+(coords[1]*nbCols) + j;
            bloc[i*nbCols + j] = matrice[pos];
        }
    }

    free(matrice);
    return bloc;

}

void factoriserLU(double *a, int size, int nbLinesParBloc, int nbColsParBloc,
        MPI_Comm lignes, MPI_Comm colonnes, int coord[]) {

    int i, j, k;
    int posRacine[2];
    int racine;

    double *c = (double*)malloc(nbColsParBloc*sizeof(double));
    double *l = (double*)malloc(nbLinesParBloc*sizeof(double));

    for(j = 0; j < size-1; j++) {
        posRacine[0] = j / nbLinesParBloc;
        posRacine[1] = j / nbColsParBloc;
        //On trouve le processeur qui doit s'occuper d'être la racine du broadcast

        if(posRacine[0] == coord[0]) {
            //C'est moi qui envoi des informations!
            memcpy(c, &a[(j%nbLinesParBloc)*nbColsParBloc], 
                    nbColsParBloc*sizeof(double));
        }

        MPI_Bcast(c, nbColsParBloc, MPI_DOUBLE, posRacine[0], colonnes);

        if(posRacine[1] == coord[1]) {
            //C'est à moi de préparer les l
            int posInterne = j % nbColsParBloc;
            for(i = 0; i < nbLinesParBloc; i++) {
                l[i] = a[i*nbColsParBloc+posInterne] / c[posInterne];
                if(coord[0]*nbLinesParBloc + i > j) {
                    a[i*nbColsParBloc+posInterne] = l[i];
                }
            }
        }

        MPI_Bcast(l, nbLinesParBloc, MPI_DOUBLE, posRacine[1], lignes);

        //On est prêt à faire la réduction
        for(i = 0; i < nbLinesParBloc; i++) {
            if(coord[0]*nbLinesParBloc + i > j) {
                //On a un traitement à faire
                for(k = 0; k < nbColsParBloc; k++) {
                    if(coord[1]*nbColsParBloc + k > j) {
                        a[i*nbColsParBloc+k] = 
                            a[i*nbColsParBloc+k] - l[i]*c[k];
                    }
                }
            }
        }

    }
    if(c) free(c);
    if(l) free(l);
}

