set title "lower bound improvement - sprint instances"
set ylabel "average gap (%)"
set xlabel "time (sec)"
set logscale x
set xrange [1:100]
set yrange [5:20]

set size 0.6,0.6

set terminal jpeg
set output "lbImprSprint.jpg"
plot "avImprSprint.txt" using 1:2 title "eclq" with lines, "avImprSprintCGL.txt" using 1:2 title "cgl" with lines


set terminal postscript enhanced 
set output "lbImprSprint.eps"
plot "avImprSprint.txt" using 1:2 title "eclq" with lines, "avImprSprintCGL.txt" using 1:2 title "cgl" with lines
