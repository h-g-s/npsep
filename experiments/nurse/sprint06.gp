set title "Lower Bound Improvement: sprint06"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint06Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint06_j.gpd" using 3:5 title "npsep" with linespoints,          "sprint06_j_cgl.gpd" using 3:5 title "CGL" with linespoints
