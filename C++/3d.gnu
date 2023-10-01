set terminal png size 1200,800 crop
set output '3d.jpg

set hidden3d
set grid
show grid

load 'blues.pal'
splot 'output.txt' using 1:2:5 lt -1
