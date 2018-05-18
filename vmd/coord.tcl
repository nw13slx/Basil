proc coordination {string1 string2 radius1 } {
  set sel1 [atomselect top "$string1" ]
  set coord1 [ $sel1 get {x y z} ]
  set id1 [ $sel1 get index]
  set n1 [ $sel1 num ]

  set bonds 0

  for { set i 0 } { $i < $n1 } { incr i } {
    set i1 [ lindex $id1 $i ]
    set n_ngh [ [atomselect top "($string2) and (not index $i1) and within $radius1 of index $i1" ] num ]
    puts "$i1 $n_ngh"
  }
}

