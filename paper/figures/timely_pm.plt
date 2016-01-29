load 'style.gnu'

set output "timely_stability.eps"

set xrange [0:50]
set yrange [-30:90]

set ytics (-30,0,30,60,90)

set xlabel "Number of Flows"
set ylabel "Phase Margin"

set key top right

plot \
        "pm_timely_delay0_fixed.txt" using ($1):($2) ti "0us delay" w lines ls 4 lc rgb "blue",\
        "pm_timely_delay50_fixed.txt" using ($1):($2) ti "50us delay" w lines ls 2 lc rgb "dark-green",\
        "pm_timely_delay100_fixed.txt" using ($1):($2) ti "100us delay" w lines ls 1 lc rgb "red"
        
