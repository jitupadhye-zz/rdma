load 'style.gnu'

set output "dcqcn_variable_delay.eps"

set xrange [0:0.2]
set yrange [1:140]

set xlabel "Time(s)"
set ylabel "Queue(KB)"

set key top right

plot \
        "unstable.2.0.100.dat" using ($1):($4) ti "[0,100{/Symbol m}s] feedback delay" w lines ls 1 lc rgb "blue",\
        "unstable.2.0.200.dat" using ($1):($4) ti "[0,200{/Symbol m}s] feedback delay" w lines ls 2 lc rgb "red"
