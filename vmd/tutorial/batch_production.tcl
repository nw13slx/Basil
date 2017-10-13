set count 0
set fp [open "list" r]
set file_data [read $fp]
foreach name $file_data {
  puts $name
  set output $name
  mol new $name type POSCAR first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
  render Tachyon $output "/share/apps/visualization/vmd/lib/vmd/tachyon_LINUXAMD64 -aasamples 12 %s -format PNG -o %s.png"
  mol delete $count
  incr count
}
