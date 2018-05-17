proc comp_coulomb { def } {
  set a [ atomselect top "$def" ]
  return [ coulomb $a ]
}
proc checklooseatom { t_ele buffer def } {
  #check whethere there are unconnected atoms
  set dangling false
  set QM1 [ atomselect top "not type $t_ele and $def" ]
  set QM2 [ atomselect top "type $t_ele and $def" ]
  set id_QM1 [ $QM1 get index ]
  set id_QM2 [ $QM2 get index ]
  foreach id2 $id_QM2 {
    set ngh [ atomselect top "(not type $t_ele) and ($def) and (within $buffer of index $id2)" ]
    set n_ngh [ $ngh num ]
    if { $n_ngh < 1  } {
      set dangling true
    }
  }
  foreach id1 $id_QM1 {
    set ngh [ atomselect top "(type $t_ele) and ($def) and (within $buffer of index $id1)" ]
    set n_ngh [ $ngh num ]
    if { $n_ngh < 1  } {
      set dangling true
    }
  }
  return $dangling

}
proc filter { input output } {

  set fd [open $input "r"]
  set filedata [read $fd]
  close $fd
  set data [split $filedata "\n"]   
  set last [ expr [ llength $data ]-1 ]
  set data [ lreplace  $data $last $last ]
  
  set t_ele O
  set buffer 3
  
  set alleng {}
  
  set fout [ open $output "w" ]
  foreach line  $data {
  
    #read the structure
    set items [split $line]
    set last [ expr [ llength $items ]-1 ]
    set lab [ lindex $items 0 ]
    set definition [ lrange $items 2 $last ]
    set def [join $definition ]
  
    set eng [ comp_coulomb $def ]
    #check whether it is computed before
    if { [ dict exists $alleng $eng ] < 1 } {
      dict incr alleng $eng
      set dangling [ checklooseatom $t_ele $buffer $def ]
  
      if { $dangling == "false" } {
        puts "$lab $eng $def"
        puts $fout "$lab $eng $def"
      } 
    }
  }
  close $fout
}


proc autogeneration { input } {

  global active_string
  global ecp_string
  global ecp_q
  global shell_Q
  global scheme 
  set fd [open $input "r"]
  set filedata [read $fd]
  close $fd
  set data [split $filedata "\n"]   
  set last [ expr [ llength $data ]-1 ]
  set data [ lreplace  $data $last $last ]
  
  set t_ele O
  set buffer 3
  
  
  foreach line  $data {
  
    #read the structure
    set items [split $line]
    set last [ expr [ llength $items ]-1 ]
    set lab [ lindex $items 0 ]
    set definition [ lrange $items 2 $last ]
    set def [join $definition ]
  
    selectQM "$def" "$active_string" "$ecp_string" $ecp_q $shell_Q "$lab" $scheme 
    puts "$lab $def"
  }
}

