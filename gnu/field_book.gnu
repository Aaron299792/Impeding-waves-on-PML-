set terminal pdf size 6,4 enhanced color font 'Helvetica,8'
set output "field_flipbook.pdf"

set pm3d map
set palette rgbformulae 33,13,10
unset key
set size ratio 1

do for [i=0:39] {
    # Add page number and Z value as title
    set label 1 sprintf("Page %d/40 - Z = %d cm", i+1, i) at graph 0.5, 1.05 center
    set xlabel "X (cm)"
    set ylabel "Y (cm)"
    
    start_row = i * 1600 + 1
    end_row = (i + 1) * 1600
    
    plot "E_field.txt" every ::start_row::end_row using 1:2:7 with image
    
    # Remove label for next plot
    unset label 1
}
