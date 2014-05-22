set encoding utf8
set title "Accélération selon le nombre de noeuds, sans pivot partiel" 
set xlabel "Nombre de noeuds"
set ylabel "Accélération"

set yrange [0:10]
set xrange [0:35]


set terminal postscript landscape color
set output "noeuds.eps"

plot 'noeuds_acc.dat' u 1:2 title 'Lignes adjacentes' with linespoints, \
'noeuds_acc.dat' u 1:4 title 'Lignes cyclique' with linespoints, \
'noeuds_acc.dat' u 1:6 title 'Bloc' with linespoints