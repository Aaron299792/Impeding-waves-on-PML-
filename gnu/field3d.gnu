set terminal pdfcairo enhanced size 10,3 font 'Arial,12'
set output '../plots/magnetic_field3dY38.pdf'
set title 'Magnetic Field Y Component'
set cblabel 'Field Strength'
set view 60, 45, 1, 1

set dgrid3d 40,40
set hidden3d

#set palette defined (0 "dark-blue", 0.3 "purple", 0.6 "magenta", 0.8 "orange", 1.0 "yellow" )
#set palette defined (0 "dark-blue", 0.3 "blue", 0.6 "cyan", 0.8 "green", 1.0 "yellow")
set palette defined (0 "blue", 0.3 "cyan", 0.6 "green", 0.8 "orange", 1.0 "yellow")
set cblabel 'Field Strength'
set pm3d

set multiplot layout 1,2

#set xlabel 'X (cm)'
#set ylabel 'Y (cm)'
#set zlabel 'Ez'
#splot '../data/xy_plot26.txt' using 1:2:6 with pm3d title 'Z = 26'

set xlabel 'X (cm)'
set ylabel 'Z (cm)'
set zlabel 'Hz'
splot '../data/xzh_plot38.txt' using 1:3:6 with pm3d title 'Y = 38'

set xlabel 'Y (cm)'
set ylabel 'Z (cm)'
set zlabel 'Hz'
splot '../data/yzh_plot38.txt' using 2:3:5 with pm3d title 'X = 38'

unset multiplot
