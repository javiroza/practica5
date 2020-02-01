set term png
set output "P5-1920-fig1.png"

set xlabel "x"
set ylabel "p(x)"
#set yrange[0:2]

set title "Sampling per acceptació i rebuig"
set key top left

# Primer es fa el plot de l'histograma amb barres d'error i després el de la funció exacta
plot "P5-1920-res.dat" index 0 using 1:2:4 w yerrorbars t"Montecarlo","" u 1:2 w histeps t"PDF",(12/(pi*(2*(pi**2)-3)))*(x**2)*((sin(x))**2) t"exact"
#pause -1



 
