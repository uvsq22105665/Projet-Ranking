set term png

set xlabel "itérations"
set ylabel "ln ( norme )"

set xrange [0:80]
set yrange [-12:0]

set title "wb-cs-stanford"
set output "plot-wb-cs-stanford.png"
plot "graphes/puissances-wb-cs-stanford.txt" using 1:2 title "puissances" with lines, \
		 "graphes/aitken-wb-cs-stanford.txt" using 1:2 title "aitken" with lines

set title "Stanford"
set output "plot-Stanford.png"
plot "graphes/puissances-Stanford.txt" using 1:2 title "puissances" with lines, \
		"graphes/aitken-Stanford.txt" using 1:2 title "aitken" with lines

