proc addbond {type1 type2 radius} {
  set sel1 [atomselect top "type $type1" ]
  set coord1 [ $sel1 get {x y z} ]
  set id1 [ $sel1 get index]
  set n1 [ $sel1 num ]

  set sel2 [atomselect top "type $type2" ]
  set coord2 [ $sel2 get {x y z} ]
  set id2 [ $sel2 get index]
  set n2 [ $sel2 num ]

  puts [ format "%s type has %g atoms, %s type has %g atoms" $type1 $n1 $type2 $n2 ]

  set bonds 0

  for { set i 0 } { $i < $n1 } { incr i } {
    for { set j 0 } { $j < $n2 } { incr j } {
      set i1 [ lindex $id1 $i ]
      set i2 [ lindex $id2 $j ]
      if { $i1 != $i2 } {
        set c1 [ lindex $coord1 $i ]
        set c2 [ lindex $coord2 $j ]
        set dx [ veclength [vecsub $c1 $c2 ]]
        if { $dx < $radius } {
          topo delbond $i1 $i2
          topo addbond $i1 $i2
          #puts $dx
          incr bonds
        }
      }
    }
  }
  puts [format "%g bonds are added" $bonds ]
}
