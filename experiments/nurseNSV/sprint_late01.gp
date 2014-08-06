set title "Lower Bound Improvement: sprint_late01"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_late01Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_late01.gpd" using 3:5 title "npsep" with linespoints,          "sprint_late01_cgl.gpd" using 3:5 title "CGL" with linespoints
