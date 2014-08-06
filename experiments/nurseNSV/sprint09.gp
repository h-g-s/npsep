set title "Lower Bound Improvement: sprint09"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint09Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint09.gpd" using 3:5 title "npsep" with linespoints,          "sprint09_cgl.gpd" using 3:5 title "CGL" with linespoints
