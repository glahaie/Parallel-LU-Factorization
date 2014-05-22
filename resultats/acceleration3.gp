set encoding utf8
set title "Accélération selon la taille de la matrice, sans pivot partiel" 
set xlabel "Taille de la matrice"
set ylabel "Accélération"

set yrange [0:10]
set xrange [0:10000]


set terminal postscript landscape color
set output "taille.eps"

plot 'taille_acceleration.dat' u 1:2 title 'Lignes adjacentes' with linespoints, \
'taille_acceleration.dat' u 1:4 title 'Lignes cyclique' with linespoints, \
'taille_acceleration.dat' u 1:6 title 'Blocs' with linespoints