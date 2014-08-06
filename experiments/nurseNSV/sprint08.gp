set title "Lower Bound Improvement: sprint08"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint08Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint08.gpd" using 3:5 title "npsep" with linespoints,          "sprint08_cgl.gpd" using 3:5 title "CGL" with linespoints
