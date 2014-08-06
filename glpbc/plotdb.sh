cat plotdb.scr | sed "s/air04/$1/g" > plot$1.scr
gnuplot plot$1.scr 
