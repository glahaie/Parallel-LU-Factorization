#include <string.h>
#include <stdlib.h>
#include <stdio.h>



int main(int agc, char* argv[]) {
    char *line = malloc(256);
    size_t len;
    ssize_t end;
    char *endNb;
    FILE *file;
    int i = 0;
    int j;
    char *ptr;

    double **matrice = (double**)malloc(6*sizeof(double));

    file = fopen("solution3.txt", "r");

    if(file == NULL) {
        printf("Erreur fichier\n");
        exit(-1);
    }

    while((end = getline(&line, &len, file)) != -1) {
        matrice[i] = (double *)malloc(6*sizeof(double));

        for(j = 0, ptr = strtok(line, " "); ptr != NULL; j++, ptr = strtok(NULL, " ")) {
            matrice[i][j] = (double)strtod(ptr, &endNb);
        }
        i++;
    }

    for(i = 0; i < 6; i++) {
        for(j = 0; j < 6; j++) {
            printf("%f ", matrice[i][j]);
        }
        printf("\n");
    }

    fclose(file);

    return 0;
}
