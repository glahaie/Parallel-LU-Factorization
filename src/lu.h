/*
* lu.h
* Par Guillaume Lahaie et Guy Francoeur
*
* Librairie pour la decomposition LU
*/
#ifndef LU_H_
#define LU_H_

#include <mpi.h>

//Lecture d'un fichier de matrices pour un test précis -- li la matrice au
//complet, c'est le rôle du noeud de prendre les bonnes lignes ensuite.
double *readMatrix(const char *nomFichier, int matrix_size, int NbLines);

//pour avoir une matrice de test ou resultat
void testMatrix(int n, int in, double *a);

//ecriture de la matrice
void writeMatrix(int n, int m, double *a, char *s);

//Impression d'une matrice
void printMatrix(int n, int m, double *a);

//Initialiser une matrice de dimension n*n aqvec des nombres aléatoires
void initializeMatrix(int n, int m, double *a);

//Vérifie si on obtient bien le bon résultat pour un test
void checkResult(double *resultat, double *a, int n, int nbLines, MPI_Comm comm);

//Trouve le noeud contenant la valeur maximale
int trouverMax(double *maxes, int nodes);
//---------------------------------------------------------------------

#endif
