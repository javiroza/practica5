set term png
set output "P5-1920-fig3.png"

set xlabel "x"
set ylabel "p(x)"
set xrange[0:3]
set yrange[0:3]

set title "Sampling per acceptació i rebuig"
set key top left

# Primer es fa el plot de les barres d'error, després el de l'histograma, i després el de la funció exacta
plot "aux.dat" using 1:2:4 w yerrorbars t"Montecarlo","" u 1:2 w histeps t"PDF",(pi*exp(-pi*x)) t"exact"
#pause -1



 
