set terminal pdfcairo size 6,4
set output '../plots/amplitud.pdf'

set xlabel 'cm'
set xrange [-9:9]
set grid
set multiplot layout 3,1 margins 0.1,0.95,0.1,0.95 spacing 0.08,0.15

set ylabel 'Amplitude'
set yrange [0:0.1]
plot '../data/amp_data.txt' skip 1 using 1:2 with line title '50 MHz' linewidth 1.2 linecolor '#000000', '../data/bessel.txt' using 1:2 with points pt 9 title 'Bessel'


set ylabel 'Amplitude'
set yrange [0:0.5]
plot '../data/amp_data.txt' skip 1 using 1:3 with line title '100 MHz' linewidth 1.2 linecolor '#000000', '../data/bessel.txt' using 1:3 with points pt 9 title 'Bessel'

set ylabel 'Amplitud'
set yrange[0:0.7]
set ytics 0.2
plot '../data/amp_data.txt' skip 1 using 1:4 with line title '200 MHz' linewidth 1.2 linecolor '#000000', '../data/bessel.txt' using 1:4 with points pt 9 title 'Bessel'

unset multiplot
set output
