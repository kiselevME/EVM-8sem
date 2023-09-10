set terminal png size 1280,941 crop
set output 'graphics.png'

set autoscale

set title "My polynom"
set grid
set xrange [-1:1]
set yrange [-1:1]

plot x*x title "Function", \
"data_1.dat" title "My polynom"
