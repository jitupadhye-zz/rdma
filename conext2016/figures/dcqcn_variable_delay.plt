load 'style.gnu'

set output "dcqcn_variable_delay.eps"

set xrange [0:0.1]
set yrange [0:1000]

set xlabel "Time(s)"
set ylabel "Queue(KB)"

set key top right

plot \
        "unstable.10.85.15.dat" using ($1):($12) ti "[85us,100us] feedback delay" w lines ls 1 lc rgb "blue"
