#
# Makefile pour les différentes versions de la factorisation LU
#
# Par Guillaume Lahaie et Guy Francoeur
#
# Pour rouler les tests: make tests
#
#
NP = 30
MATRIX_SIZE = 1500
ROWS = 5
COLS = 6

all: lib lu_seq lu_mpi_row_adjacent lu_mpi_row_cyclic lu_mpi_cart

lib: src/lu.c src/lu.h
		mpicc -O2 -o src/lu.o -c src/lu.c -std=gnu99

lu_seq: src/lu_seq.c lib 
		mpicc -O2 -o lu_seq -std=gnu99 src/lu_seq.c src/lu.o

lu_mpi_row_adjacent: src/lu_mpi_row_adjacent.c lib
		mpicc -O2 -DNDEBUG -o lu_mpi_row_adjacent -std=gnu99 src/lu_mpi_row_adjacent.c src/lu.o

lu_mpi_row_cyclic: src/lu_mpi_row_cyclic.c lib
		mpicc -O2 -DNDEBUG -o lu_mpi_row_cyclic -std=gnu99 src/lu_mpi_row_cyclic.c src/lu.o

lu_mpi_cart: src/lu_mpi_cart.c lib
		mpicc -O2 -DNDEBUG -o lu_mpi_cart -std=gnu99 src/lu_mpi_cart.c src/lu.o

tests: test1 test2 test3 test4

#Tests d'une matrice 4*4
test1: all
	./lu_seq -m 4 -f "test/matrice1.txt" -s "test/solution1.txt"
	./lu_seq -m 4 -f "test/matrice1.txt" -s "test/solution1_pivot.txt" -p
	mpirun -n 2 ./lu_mpi_row_adjacent -m 4 -f "test/matrice1.txt" -s "test/solution1.txt"
	mpirun -n 2 ./lu_mpi_row_adjacent -m 4 -f "test/matrice1.txt" \
		-s "test/solution1_pivot.txt" -p
	mpirun -n 2 ./lu_mpi_row_cyclic -m 4 -f "test/matrice1.txt" -s "test/solution1.txt"
	mpirun -n 2 ./lu_mpi_row_cyclic -m 4 -f "test/matrice1.txt" \
		-s "test/solution1_pivot.txt" -p
	mpirun -n 4 ./lu_mpi_cart -m 4 -x 2 -y 2 -f "test/matrice1.txt" -s "test/solution1.txt"

#Tests d'une matrice 4*4
test2: all
	./lu_seq -m 4 -f "test/matrice2.txt" -s "test/solution2.txt"
	./lu_seq -m 4 -f "test/matrice2.txt" -s "test/solution2_pivot.txt" -p
	mpirun -n 2 ./lu_mpi_row_adjacent -m 4 -f "test/matrice2.txt" -s "test/solution2.txt"
	mpirun -n 2 ./lu_mpi_row_adjacent -m 4 -f "test/matrice2.txt" -s "test/solution2_pivot.txt" -p
	mpirun -n 2 ./lu_mpi_row_cyclic -m 4 -f "test/matrice2.txt" -s "test/solution2.txt"
	mpirun -n 2 ./lu_mpi_row_cyclic -m 4 -f "test/matrice2.txt" -s "test/solution2_pivot.txt" -p
	mpirun -n 4 ./lu_mpi_cart -m 4 -x 2 -y 2 -f "test/matrice2.txt" -s "test/solution2.txt"

#Tests d'une matrice 6*6
test3: all
	./lu_seq -m 6 -f "test/matrice3.txt" -s "test/solution3.txt"
	./lu_seq -m 6 -f "test/matrice3.txt" -s "test/solution3_pivot.txt" -p
	mpirun -n 3 ./lu_mpi_row_adjacent -m 6 -f "test/matrice3.txt" -s "test/solution3.txt"
	mpirun -n 3 ./lu_mpi_row_adjacent -m 6 -f "test/matrice3.txt" \
		-s "test/solution3_pivot.txt" -p
	mpirun -n 3 ./lu_mpi_row_cyclic -m 6 -f "test/matrice3.txt" -s "test/solution3.txt"
	mpirun -n 3 ./lu_mpi_row_cyclic -m 6 -f "test/matrice3.txt" \
		-s "test/solution3_pivot.txt" -p
	mpirun -n 6 ./lu_mpi_cart -m 6 -x 2 -y 3 -f "test/matrice3.txt" -s "test/solution3.txt"

#Tests d'une matrice 30*30
test4: all
	./lu_seq -m 30 -f "test/matrice4.txt" -s "test/solution4.txt"
	./lu_seq -m 30 -f "test/matrice4.txt" -s "test/solution4_pivot.txt" -p
	mpirun -n 10 ./lu_mpi_row_adjacent -m 30 -f "test/matrice4.txt" -s "test/solution4.txt"
	mpirun -n 10 ./lu_mpi_row_adjacent -m 30 -f "test/matrice4.txt" \
		-s "test/solution4_pivot.txt" -p
	mpirun -n 10 ./lu_mpi_row_cyclic -m 30 -f "test/matrice4.txt" -s "test/solution4.txt"
	mpirun -n 10 ./lu_mpi_row_cyclic -m 30 -f "test/matrice4.txt" \
		-s "test/solution4_pivot.txt" -p
	mpirun -n 30 ./lu_mpi_cart -m 30 -x 10 -y 3 -f "test/matrice4.txt" -s "test/solution4.txt"

#Tests de performance pour une valeur de grandeur de matrice
#Il faut spécifier aussi le nombre de noeuds et la topologie
#pour l'approche par blocs
mesures: all
	./lu_seq -m ${MATRIX_SIZE}
	./lu_seq -m ${MATRIX_SIZE} -p
	mpirun -n ${NP} ./lu_mpi_row_adjacent -m ${MATRIX_SIZE}
	mpirun -n ${NP} ./lu_mpi_row_adjacent -m ${MATRIX_SIZE} -p
	mpirun -n ${NP} ./lu_mpi_row_cyclic -m ${MATRIX_SIZE}
	mpirun -n ${NP} ./lu_mpi_row_cyclic -m ${MATRIX_SIZE} -p
	mpirun -n ${NP} ./lu_mpi_cart -m ${MATRIX_SIZE} -x ${ROWS} -y ${COLS}

clean:
		rm -f lu_seq_pivot
		rm -f lu_seq
		rm -f lu_mpi_row_adjacent
		rm -f lu_mpi_row_adjacent_pivot
		rm -f lu_mpi_row_cyclic
		rm -f lu_mpi_row_cyclic_pivot
		rm -f lu_mpi_cart
		rm -f src/*.o
