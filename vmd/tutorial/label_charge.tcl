#this code is adapted from label_atom http://www.ks.uiuc.edu/Research/vmd/current/ug/node127.html
#and pdbbfactor http://www.ks.uiuc.edu/Research/vmd/script_library/scripts/pdbbfactor/
#it only label the last frame though...
#
#usage: source nbo.tcl
#       label_bfactor <pdb file name>

proc label_bfactor { fname } {
  mol new $fname waitfor all
  set all [atomselect top all]
  set frame 0
  set in [open $fname r]
  set beta {}
  while { [gets $in line] != -1 } {
    switch -- [string range $line 0 3] {
      END {
        $all frame $frame
        $all set user $beta
        set beta {}
        incr frame
      }
      ATOM -
      HETA {
        lappend beta [string range $line 61 66]
      }
    }
  }

  set sel [atomselect top "all"]
  if {[$sel num] < 1} {
    error "label_atom: '$selection_string' must select at least 1 atom"
  }
  #get the coordinates of the atom
  set x [$sel get {x y z}] 
  set n [ $sel num]
  set q [$sel get beta]
  for { set i 0 } { $i < $n } { incr i } {
  # and draw the text
    set qq [lindex $q $i]
    puts $i 
    puts $qq
    puts [lindex $x $i ]
    if  { $qq > 0 } {
      draw color red
      draw text [lindex $x $i ] [ format "%.2f" $qq ] 
    } elseif { $qq < 0 } {
      draw color blue
      #draw text [lindex $x $i ] [ format "%.2f" $qq ]
      draw text [lindex $x $i ] [ format "%.2f" $qq ] 
    } else {
      draw color black
    }
  }
}
  
