load 'style.gnu'

set output "dcqcn_convergence_time.eps"

set xrange [0:100]
set yrange [0:0.4]

set xlabel "Number of Flows"
set ylabel "Convergence Time (s)"

set key top right

plot \
        "convergence_rai10.txt" using ($1):($2) ti "Rai=10Mbps" w lines ls 4 lc rgb "blue",\
        "convergence_kmax1000.txt" using ($1):($2) ti "Kmax=1000KB" w lines ls 3 lc rgb "green",\
        "convergence_g16.txt" using ($1):($2) ti "g=1/16" w lines ls 2 lc rgb "magenta",\
        "convergence_default.txt" using ($1):($2) ti "Default" w lines ls 1 lc rgb "red"
