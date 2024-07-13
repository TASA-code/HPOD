set terminal png size 5000,3600 crop
set output '../../plot/3d.png'

set hidden3d
set grid
show grid
set view equal xyz 

set multiplot layout 2,2 spacing 0.4

# Plot 1: ECI data with sphere
set title 'ECI Data'
set parametric
set angles degrees
set isosamples 25
set urange [0:360]
set vrange [-90:90]
set view 75,320

# Define sphere radius
r = 6378000

splot r*cos(v)*cos(u), r*cos(v)*sin(u), r*sin(v) with lines lt 2 lc rgb '#000000' notitle, \
      'output.txt' using 3:4:5 with points lt -1 ps 100 title 'ECI Data'

unset parametric

# Plot 2: ECEF data with sphere
set title 'ECEF Data'
set parametric
set angles degrees
set isosamples 25
set urange [0:360]
set vrange [-90:90]

splot r*cos(v)*cos(u), r*cos(v)*sin(u), r*sin(v) with lines lt 2 lc rgb '#000000' notitle, \
      'output.txt' using 9:10:11 with points lt -1 ps 100 title 'ECEF Data'

unset parametric

# Plot 3: GEO data
set title 'ECI_GEO Data (2D)'
set xrange [-180:180]
set yrange [-90:90]
plot 'output.txt' using 15:16 with points pointsize 1 pt 7 lc rgb "black" title 'ECI_GEO Data'

# Plot 4: GEO data with world map
set title 'GEO Data (2D)'
set yrange [-90:90]
plot 'output.txt' using 17:18 with points pointsize 1 pt 7 lc rgb "black" title 'GEO Data'

unset multiplot