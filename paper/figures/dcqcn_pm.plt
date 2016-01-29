load 'style.gnu'

set output "dcqcn_stability.eps"

set xrange [0:100]
set yrange [-30:90]

set ytics (-30,0,30,60,90)

set xlabel "Number of Flows"
set ylabel "Phase Margin"

set key top left

plot \
        "pm_rai40_delay0.txt" using ($1):($2) ti "0us delay" w lines ls 4 lc rgb "blue",\
        "pm_rai40_delay50.txt" using ($1):($2) ti "50us delay" w lines ls 2 lc rgb "dark-green",\
        "pm_rai40_delay100.txt" using ($1):($2) ti "100us delay" w lines ls 1 lc rgb "red"
        
