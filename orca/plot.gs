#!/share/apps/gnuplot/5.0.6/bin/gnuplot -p
set title "pc QQQ" font ",15"
set style line 1 dt 2 lw 1
set term png
set output "pdos.png"
set xrange[-6:10]
pl 'pdos' u 1:2 w l lw 2 lc rgb 'black' title "Ce s"  ,\
'pdos' u 1:3 w l lw 2 lc rgb 'red' title "Ce p"  ,\
'pdos' u 1:4 w l lw 2 lc rgb 'green' title "Ce d"  ,\
'pdos' u 1:5 w l lw 2 lc rgb 'blue' title "Ce f"  ,\
'pdos' u 1:8 w l lw 2 dt "." lc rgb 'black' title "O s"   ,\
'pdos' u 1:9 w l lw 2 dt "." lc rgb 'red'  title "O p"   ,\
'pdos' u 1:15 w l lw 2 lc rgb 'black'title ""  ,\
'pdos' u 1:16 w l lw 2 lc rgb 'red'title ""  ,\
'pdos' u 1:17 w l lw 2 lc rgb 'green'title ""  ,\
'pdos' u 1:18 w l lw 2 lc rgb 'blue'title ""  ,\
'pdos' u 1:21 w l lw 2 dt "." lc rgb 'black'title ""  ,\
'pdos' u 1:22 w l lw 2 dt "." lc rgb 'red'title ""  
set term x11
set out
replot
