proc filter { centerid input output } {

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
    set sel [ atomselect top "$def" ] 
    set eng [ coulomb $sel ]
    #check whether it is computed before
    if { [ dict exists $alleng $eng ] < 1 } {
      dict incr alleng $eng
      set dangling [ count_undercoord $buffer $sel 1 ]
      puts $def
  
      if { $dangling < 1 } {
        #set m2 [ moment2nd $centerid $sel ]
        #puts "$lab $eng $m2 $def"
        #puts $fout "$lab $eng $m2 $def"
        puts "$lab $eng $def"
        puts $fout "$lab $eng $def"
      } else {
        puts "$lab has dangling bond"
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

