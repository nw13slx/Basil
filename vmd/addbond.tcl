# this code is to add or remove bonds in vmd 
# 
# addbond (type1, type2, R1, R2) function is  used to add bonds to the A-B atom 
# pairs that have is in the range of (R1,R2). A, B should be a type name and 
# (R1,R2) is the minimum/maximum cut-off radius. If the bond already exist, it 
# will delete the bond and re-add it. It shouldn't be too hard to change A and B 
# into random selection defined by "atomselect"
#
# deletebond (type1, type2, R1, R2) function is to undo the above function.
#
# no periodic boundary condition is counted. see more notes in the code.
#
# Author: Lixin Sun, nw13mifaso@gmail.com


proc addbond_selection {sel1 sel2 radius1 radius2} {
  ## for periodic boundary
  ## see more note below
  #set pbc 0
  #set cell [ pbc get ]
  #if { [ llength $cell ] ==1 } {
  #  set vcell [ lindex $cell 0 ]
  #  if { [ llength $vcell ] == 6 } {
  #    set pbc 1
  #  }
  #}
  set coord1 [ $sel1 get {x y z} ]
  set id1 [ $sel1 get index]
  set n1 [ $sel1 num ]

  set coord2 [ $sel2 get {x y z} ]
  set id2 [ $sel2 get index]
  set n2 [ $sel2 num ]

  puts [ format "group 1 has %g atoms, group 2 has %g atoms"  $n1 $n2 ]

  set bonds 0

  for { set i 0 } { $i < $n1 } { incr i } {
    for { set j 0 } { $j < $n2 } { incr j } {
      set i1 [ lindex $id1 $i ]
      set i2 [ lindex $id2 $j ]
      if { $i1 != $i2 } {
        set c1 [ lindex $coord1 $i ]
        set c2 [ lindex $coord2 $j ]
        for { set k 0 } { $k < 3 } { incr k } {
          set dk [expr [lindex $c1 $k ] - [lindex $c2 $k]]
          ## this one is used for periodic boundary. 
          ## however, vmd doesn't really have a good way to visualize 
          ## bond across periodic boundary unless you use the dynamic 
          ## bond representation
          #set dk [expr abs($dk)]
          #if { $pbc == 1 } {
          #  while { $dk > [expr 0.5*[lindex $vcell $k]] } {
          #    set dk [expr $dk-[lindex $vcell $k]]
          #  }
          #}
          set x$k $dk
        }
        set dx [ expr sqrt($x1*$x1+$x2*$x2+$x0*$x0)] 
        if { ( $dx < $radius2 ) && ( $dx > $radius1 ) } {
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
proc deletebond {type1 type2 radius1 radius2} {
  set pbc 0
  set cell [ pbc get ]
  if { [ llength $cell ] ==1 } {
    set vcell [ lindex $cell 0 ]
    if { [ llength $vcell ] == 6 } {
      set pbc 1
    }
  }

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
        for { set k 0 } { $k < 3 } { incr k } {
          set dk [expr [lindex $c1 $k ] - [lindex $c2 $k]]
          set dk [expr abs($dk)]
          if { $pbc == 1 } {
            while { $dk > [expr 0.5*[lindex $vcell $k]] } {
              set dk [expr $dk-[lindex $vcell $k]]
            }
          }
          set x$k $dk
        }
        set dx [ expr sqrt($x1*$x1+$x2*$x2+$x0*$x0)] 
        if { {$dx < $radius2} && {$dx > $radius1} } {
          topo delbond $i1 $i2
          #puts $dx
          incr bonds
        }
      }
    }
  }
  puts [format "%g bonds are deleted" $bonds ]
}

proc addbond_group {string1 string2 radius1 radius2} {
  set sel1 [atomselect top "$string1" ]
  set sel2 [atomselect top "$string2" ]
  addbond_selection $sel1 $sel2 $radius1 $radius2
}

proc addbond {type1 type2 radius1 radius2} {
  set sel1 [atomselect top "type $type1" ]
  set sel2 [atomselect top "type $type2" ]
  addbond_selection $sel1 $sel2 $radius1 $radius2
}
