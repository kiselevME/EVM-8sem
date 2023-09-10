set terminal png size 1280,941 crop
set output 'plot.png'

set autoscale

set title "Plot"
set grid
set xrange [0:1]
set yrange [0:0.0000001]

plot exp(x) title "Real function", \
"data.dat" title "My estimate"
