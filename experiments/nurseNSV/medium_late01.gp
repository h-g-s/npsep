set title "Lower Bound Improvement: medium_late01"
set xlabel "time (sec)"
set terminal jpeg
set output "medium_late01Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium_late01.gpd" using 3:5 title "npsep" with linespoints,          "medium_late01_cgl.gpd" using 3:5 title "CGL" with linespoints
