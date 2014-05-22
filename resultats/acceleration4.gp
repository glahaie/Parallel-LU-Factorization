set encoding utf8
set title "Accélération selon la taille de la matrice, avec pivot partiel" 
set xlabel "Taille de la matrice"
set ylabel "Accélération"

set yrange [0:10]
set xrange [0:10000]


set terminal postscript landscape color
set output "taille_pivot.eps"

plot 'taille_acceleration.dat' u 1:3 title 'Lignes adjacentes' with linespoints, \
'taille_acceleration.dat' u 1:5 title 'Lignes cyclique' with linespoints