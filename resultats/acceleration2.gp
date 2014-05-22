set encoding utf8
set title "Accélération selon le nombre de noeuds, avec pivot partiel" 
set xlabel "Nombre de noeuds"
set ylabel "Accélération"

set yrange [0:10]
set xrange [0:35]


set terminal postscript landscape color
set output "noeuds_pivot.eps"

plot 'noeuds_acc.dat' u 1:3 title 'Lignes adjacentes' with linespoints, \
'noeuds_acc.dat' u 1:5 title 'Lignes cyclique' with linespoints