set title "Lower Bound Improvement: long_hidden01"
set xlabel "time (sec)"
set terminal jpeg
set output "long_hidden01Perc.jpg"
set ylabel "lb improvement (%)"
plot "long_hidden01.gpd" using 3:5 title "npsep" with linespoints,          "long_hidden01_cgl.gpd" using 3:5 title "CGL" with linespoints
