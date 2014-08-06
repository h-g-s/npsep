set title "Lower Bound Improvement"
set xlabel "time (sec)"
set ylabel "lb improvement (%)"
plot "sprint_hidden01_j.gpd"  using 3:5 title "npsep" with linespoints, \
     "sprint_hidden01_j_cgl.gpd" using 3:5 title "cgl" with linespoints
