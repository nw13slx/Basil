proc selectQM {s_string name } {

  set radius 3
  set fo_QM [ open [format "%s.xyz" $name] "w"]
  set fo_MM [ open [format "%s.pc" $name] "w"]
  set fo_info [ open [format "%s.info" $name] "w"]
  
  # for periodic boundary
  # see more note below
  set pbc 0
  set cell [ pbc get ]
  if { [ llength $cell ] ==1 } {
    set vcell [ lindex $cell 0 ]
    if { [ llength $vcell ] == 6 } {
      set pbc 1
      for { set k 0 } { $k < 3 } { incr k } {
        set pbc [expr $pbc*[lindex $vcell $k]] 
      }
    }
  }
  
  set QM [ atomselect top $s_string ]
  set MM [ atomselect top [ format " not %s" $s_string ] ]
  
  set x_QM [ $QM get {x y z} ]
  set type_QM [ $QM get type ]
  set n_QM [ $QM num ]
  
  set x_MM [ $MM get {x y z} ]
  set type_MM [ $MM get type ]
  set n_MM [ $MM num ]
  puts [ format "selection string: %s " $s_string ]
  puts [ format "QM: %g MM:%g" $n_QM $n_MM ]
  puts $fo_info [ format "selection string: %s " $s_string ]
  puts $fo_MM "count_MM"
  
  for { set j 0 } { $j < $n_QM } { incr j } {
    set c2 [ lindex $x_QM $j ]
    puts $fo_QM [ format "%s %.5f %.5f %.5f" [ lindex $type_QM $j ] [ lindex $c2 0 ] [lindex $c2 1] [lindex $c2 2] ]
  }

  set count_pot 0 
  set count_MM 0


  for { set i 0 } { $i < $n_MM } { incr i } {
    set is_pot 0
    for { set j 0 } { $j < $n_QM } { incr j } {
      set c1 [ lindex $x_MM $i ]
      set c2 [ lindex $x_QM $j ]
      set type1 [ lindex $type_MM $i ] 
      if { [ string match $type1 "Ce" ] == 1 } {
        for { set k 0 } { $k < 3 } { incr k } {
          set dk [expr [lindex $c1 $k ] - [lindex $c2 $k]]
          set dk [expr abs($dk)]
          if { $pbc > 0 } {
            while { $dk > [expr 0.5*[lindex $vcell $k]] } {
              set dk [expr $dk-[lindex $vcell $k]]
            }
          }
          set x$k $dk
        }
        set dx [ expr sqrt($x1*$x1+$x2*$x2+$x0*$x0)]
        if { $dx < $radius  } {
          incr is_pot
        }
      }
    }
    if { $is_pot > 0 } {
      puts $fo_QM [ format "%s> pQ %.5f %.5f %.5f" [ lindex $type1 ] [ lindex $c1 0 ] [lindex $c1 1] [lindex $c1 2] ]
      incr count_pot
    } else {
      if { [string match $type1 "Ce" ] == 1 } {
        puts $fo_MM [ format "pQ %.5f %.5f %.5f" [ lindex $c1 0 ] [lindex $c1 1] [lindex $c1 2] ]
      } else {
        puts $fo_MM [ format "nQ %.5f %.5f %.5f" [ lindex $c1 0 ] [lindex $c1 1] [lindex $c1 2] ]
      }
      incr count_MM
    }
  }

  puts [ format "corrected MM %g pot %g" $count_MM $count_pot ]
  puts $fo_info [ format "QM: %g\nMM: %g\npot: %g" $n_QM $count_MM $count_pot ]
  

  close $fo_QM
  close $fo_MM
  close $fo_info

  set string [ format "sed -i {s/count_MM/%g/} %s.pc" $count_MM $name]
  exec /bin/sh -c $string
}
