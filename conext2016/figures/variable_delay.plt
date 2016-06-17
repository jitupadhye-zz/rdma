load 'style.gnu'

set output "variable_delay.eps"

set xrange [0:0.2]
set yrange [1:140]

set xlabel "Time(s)"
set ylabel "Queue(KB)"

set key top right

plot \
        "timely.fixed.2.100.dat" using ($1):($4) ti "TIMELY" w lines ls 2 lc rgb "red",\
        "unstable.2.0.100.dat" using ($1):($4) ti "DCQCN" w lines ls 1 lc rgb "blue"
