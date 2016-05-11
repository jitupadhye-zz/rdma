load 'style.gnu'

set output "timely_variable_delay.eps"

set xrange [0:0.2]
set yrange [0:140]

set xlabel "Time(s)"
set ylabel "Queue (KB)"

set key top right

plot \
        "timely.fixed.2.15.dat" using ($1):($4) ti "[0,15{/Symbol m}s] additional delay" w lines ls 1 lc rgb "blue",\
        "timely.fixed.2.100.dat" using ($1):($4) ti "[0,100{/Symbol m}s] additional delay" w lines ls 2 lc rgb "red"

