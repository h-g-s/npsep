set title "lower bound improvement - medium instances"
set ylabel "average gap (%)"
set xlabel "time (sec)"
set xrange [10:600]
set yrange [34:60]
set logscale x

set size 0.6,0.6

set terminal jpeg
set output "lbImprMed.jpg"
plot "avImprMed.txt" using 1:2 title "eclq" with lines, "avImprMedCGL.txt" using 1:2 title "cgl" with lines


set terminal postscript enhanced 
set output "lbImprMed.eps"
plot "avImprMed.txt" using 1:2 title "eclq" with lines, "avImprMedCGL.txt" using 1:2 title "cgl" with lines

