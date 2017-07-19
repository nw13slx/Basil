#!/share/apps/gnuplot/bin/gnuplot  -p
set title "DOS-CeO2 bulk" font ",15"

set term png
set output "pdos.png"
set xrange[-50:10]
pl 'pdos' u 1:5 w l lw 2 lc rgb 'red' title "Ce f",\
   'pdos' u 1:(-$17) w l lw 2 lc rgb 'red' title "",\
   'pdos' u 1:4 w l lw 2 lc rgb 'blue' title "Ce d",\
   'pdos' u 1:(-$16) w l lw 2 lc rgb 'blue' title "",\
   'pdos' u 1:9 w l lw 2 lc rgb 'green' title "O p",\
   'pdos' u 1:(-$21) w l lw 2 lc rgb 'green' title ""
set term x11
set out
replot
