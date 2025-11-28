set terminal dumb

set print 'xy_combined.dat'
do for [y=0:38] {
    do for [x=0:38] {
        line_num = 60843 + 40*y + x
        # Use system call to get the exact line
        cmd = sprintf("sed -n '%dp' ../data/H_field.txt", line_num)
        system(cmd)
    }
    print ""
}
set print

