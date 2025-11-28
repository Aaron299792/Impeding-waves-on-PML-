# Create animation through Z slices
set terminal gif animate delay 50 size 800,600
set output "./plots/field_animation.pdf"
set pm3d map
set palette rgbformulae 33,13,10
unset key

do for [i=0:39] {
    set xlabel "X (cm)"
    set ylabel "Y (cm)"
    set cblabel "|E| field" rotate parallel
    set title sprintf("Z = %d", i)
    plot "E_field.txt" every ::(i*1600+1)::((i+1)*1600) using 1:2:7 with image
}
