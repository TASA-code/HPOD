set terminal png size 5000,3600 crop
set output '../../plot/3d.png'

set hidden3d
set grid
show grid


set multiplot layout 2,2 spacing 0.4

# Plot 1: ECI data
set title 'ECI Data'
splot 'output.txt' using 3:4:5 lt -1 lw 2

# Plot 2: ECEF data
set title 'ECEF Data'
splot 'output.txt' using 9:10:11 lt -1 lw 2

# Plot 3: GEO data
set title 'ECI_GEO Data (2D)'
set xrange [-180:180]
set yrange [-90:90]
plot 'output.txt' using 15:16 with points pointsize 1 pt 7 lc rgb "black"

# Plot 4: GEO data
set title 'GEO Data (2D)'
# set xrange [-180:180]
set yrange [-90:90]
plot 'output.txt' using 17:18 with points pointsize 1 pt 7 lc rgb "black"


unset multiplot

