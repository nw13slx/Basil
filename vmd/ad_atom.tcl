proc ladd L {expr [join $L +]+0} ;
set dist 2

proc addbond_group {string1 string2 radius1 radius2} {
  global dist
  set sel1 [atomselect top "$string1" ]
  set coord1 [ $sel1 get {x y z} ]
  set id1 [ $sel1 get index]
  set n1 [ $sel1 num ]

  set bonds 0
  set count 0
  set O_add [ open "position" "w"]
  for { set i 0 } { $i < $n1 } { incr i } {
    set i1 [ lindex $id1 $i ]
    set sel2 [atomselect top "$string2 and within $radius2 of index $i1" ]
    set n2 [ $sel2 num ]
    if { $n2 != 6 } {
      set x [ lindex $coord1 $i ]
      set x1 [ lindex $x 0 ]
      set y1 [ lindex $x 1 ]
      set z1 [ lindex $x 2 ]
      set x2 [ ladd [ $sel2 get x ] ]
      set y2 [ ladd [ $sel2 get y ] ]
      set z2 [ ladd [ $sel2 get z ] ]
      set dx [ expr $x2 / $n2 -$x1 ]
      set dy [ expr $y2 / $n2 -$y1 ]
      set dz [ expr $z2 / $n2 -$z1 ]
      set d [ expr sqrt ($dx*$dx+$dy*$dy+$dz*$dz ) ]
      set nx [ expr $x1 - ($dx)/$d*$dist ]
      set ny [ expr $y1 - ($dy)/$d*$dist ]
      set nz [ expr $z1 - ($dz)/$d*$dist ]
      puts $O_add "$nx $ny $nz"
      incr count
    }
  }
  puts [format "%g oxygen can be added" $count ]
  close $O_add
}
