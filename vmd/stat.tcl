proc dipole { def } {
  set a [ atomselect top "$def" ]
  set x [ $a get x ]
  set y [ $a get y ]
  set z [ $a get z ]
  set q [ $a get charge ]
  set n [ $a num ]
  set dx 0
  set dy 0
  set dz 0
  for { set i 0 } { $i < $n } { incr i } {
    set dx [ expr $dx + [ lindex $q $i ] * [lindex $x $i ] ]
    set dy [ expr $dy + [ lindex $q $i ] * [lindex $y $i ] ]
    set dz [ expr $dz + [ lindex $q $i ] * [lindex $z $i ] ]
  }
  #puts "dipole $dx $dy $dz"
  return [ list $dx $dy $dz ]
}

proc moment2nd { centerid a } {
  set coord0 [ lindex [ [atomselect top "index $centerid"] get {x y z} ] 0]
  set x0 [ lindex $coord0 0]
  set y0 [ lindex $coord0 1]
  set z0 [ lindex $coord0 2]
  set x [ $a get x ]
  set y [ $a get y ]
  set z [ $a get z ]
  set q [ $a get charge ]
  set n [ $a num ]
  set tally 0
  for { set i 0 } { $i < $n } { incr i } {
    set r2 [ expr ([lindex $x $i ]-$x0)**2 +  ([lindex $y $i ]-$y0)**2 +  ([lindex $z $i ]-$z0)**2 ]
    set tally [ expr $tally + $r2 ]
  }
  return $tally
}

proc total_coordination { buffer a } {
  set type [ $a get type ]
  set index [ $a get index ]
  set n [ $a num ]
  set buffer 3
  set tot 0
  for { set i 0 } { $i < $n } { incr i } {
    set id [ lindex $index $i ]
    set t [ lindex $type $i ]
    set ngh [ atomselect top "(not type $t) and (index $index) and (within $buffer of index $id)" ]
    set tot [ expr $tot+[$ngh num]]
  }
  return $tot
}

proc count_undercoord { buffer a c } {
  set type [ $a get type ]
  set index [ $a get index ]
  set n [ $a num ]
  set d 0
  for { set i 0 } { $i < $n } { incr i } {
    set id [ lindex $index $i ]
    set t [ lindex $type $i ]
    set ngh [ atomselect top "(not type $t) and (index $index) and (within $buffer of index $id)" ]
    set n_ngh [ $ngh num ]
    if { $n_ngh < $c  } {
      incr d
    }
  }
  return $d
}

proc ladd L {expr [join $L +]+0} ;

set K 167101.001981787

proc count_ele { list } {
  set counters {}
  foreach item $list {
    dict incr counters $item
  }
  return $counters 
}

proc coulomb { sel } {
  set index [ $sel get index ]
  set q [ $sel get charge ]
  set n [ $sel num ]
  set pot 0
  for { set i 0 } { $i < $n } { incr i } {
    for { set j [expr $i+1] } { $j < $n } { incr j } {
      set id1 [ lindex $index $i ] 
      set id2 [ lindex $index $j ] 
      set b [ measure bond [ list $id1 $id2 ] ]
      set pot [ expr $pot + [ lindex $q $i ]*[lindex $q $j]/$b ] 
    }
  }
  return $pot
}

proc coulomb_cross { QMi sel } {
  set index [ $sel get index ]
  set q [ $sel get charge ]
  set n [ $sel num ]
  set pot 0
  set QMindex [ $QMi get index ]
  set QMq [ $QMi get charge ]
  set QMn [ $QMi num ]
  for { set i 0 } { $i < $n } { incr i } {
    for { set j 0 } { $j < $QMn } { incr j } {
      set id1 [ lindex $index $i ] 
      set id2 [ lindex $QMindex $j ] 
      set b [ measure bond [ list $id1 $id2 ] ]
      set pot [ expr $pot + [ lindex $q $i ]*[lindex $QMq $j]/$b ] 
    }
  }
  return $pot
}

proc hash_coulomb_cross { rest sel } {
  set index [ $sel get index ]
  set q [ $sel get charge ]
  set n [ $sel num ]
  set restindex [ $rest get index ]
  set restq [ $rest get charge ]
  set restn [ $rest num ]
  set coul {}
  for { set i 0 } { $i < $n } { incr i } {
    set pot 0
    set id1 [ lindex $index $i ] 
    for { set j 0 } { $j < $restn } { incr j } {
      set id2 [ lindex $restindex $j ] 
      set b [ measure bond [ list $id1 $id2 ] ]
      set pot [ expr $pot + [ lindex $q $i ]*[lindex $restq $j]/$b ] 
    }
    puts "... $id1 $pot"
    dict append coul $id1 $pot
  }
  return $coul
}

proc sum_coulomb_cross { coul sel } {
  set pot 0
  foreach id $sel {
    if { [ dict exist $coul $id ] >0 } {
      set pot [ expr $pot + [ dict get $coul $id ] ]
    } else {
      puts "ERROR the coulombic force is not computed ahead"
      return 0
    }
  }
  return $pot
}


