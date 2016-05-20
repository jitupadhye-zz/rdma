load 'style.gnu'

set output "dcqcn_bufferbloat.eps"

set xrange [0:0.2]
set yrange [1:140]

set xlabel "Time(s)"
set ylabel "Queue(KB)"

set key top right

plot \
        "unstable.2.85.0.dat" using ($1):($4) ti "Mark on dequeuing" w lines ls 1 lc rgb "blue",\
        "unstable.2.85.0.bufferbloat.dat" using ($1):($4) ti "Mark on enqueuing" w lines ls 2 lc rgb "red"
