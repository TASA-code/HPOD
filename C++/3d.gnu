set terminal png size 2000,2400 crop
set output '3d.jpg'

set hidden3d
set grid
show grid

load 'blues.pal'

set multiplot layout 3,1 spacing 0.2

# Plot 1: ECI data
set title 'ECI Data'
splot 'build/bin/ECI.txt' using 3:4:5 lt -1 lw 2

# Plot 2: ECEF data
set title 'ECEF Data'
splot 'build/bin/ECEF.txt' using 1:2:3 lt -1 lw 2

# Plot 3: GEO data
set title 'GEO Data (2D)'
plot 'build/bin/GEO.txt' using 1:2 with points pointsize 1 pt 7 lc rgb "black"

unset multiplot

