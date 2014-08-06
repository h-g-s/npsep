set title "lower bound improvement - long instances"
set ylabel "average gap (%)"
set xlabel "time (sec)"
set xrange [10:600]
set logscale x
set size 0.6,0.6

set terminal jpeg
set output "lbImprLong.jpg"
plot "avImprLong.txt" using 1:2 title "eclq" with lines, "avImprLongCGL.txt" using 1:2 title "cgl" with lines


set terminal postscript enhanced 
set output "lbImprLong.eps"
plot "avImprLong.txt" using 1:2 title "eclq" with lines, "avImprLongCGL.txt" using 1:2 title "cgl" with lines
