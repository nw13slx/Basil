proc ladd L {expr [join $L +]+0} ;

#requirement: type should be element symbol instead of number. point charge should be labelled as "F"
#             the configuration does not have periodic boundary

#argument 1: QM_string definition of the QM zone. can be anything that vmd recognize
#argument 2: active_string. region 3
#argument 3: buffer_string the definition of potential group. i.e. "type Ti" or "none"
#argument 5: name i.e. "CaO2"
#argument 6: partial or formal charge scheme. either "p" or "f"

#partial charge scheme:
#   region 1: QM atom, with charge 0
#   region 2: potential, use partial charge or formal charge
#
#   compensated charge computed from: region 1 original charge + region 2 current charge
#   delta charge: compensated charge divided by the number of atoms in region 3 and region 4
#
#   region 3: active MM, use partial charge, and partial charge shell
#   region 4: frozen MM, use original charge, no shell

#formal charge scheme:
#   region 1: QM atom, with charge 0
#   region 2: potential, use partial charge or formal charge
#
#   compensated charge computed from: region 1 original charge + region 2 current charge
#   delta charge: compensated charge divided by the number of atoms in region 3 and region 4
#
#   region 3: active MM, each 
#   region 4: frozen MM, use original charge

proc selectQM_link { QM_string active_string buffer_string shell_Q name scheme } {



  #repeats the argument
  puts [ format "original command: selectQM_link \"%s\" \"%s\" \"%s\" %s %s " $QM_string $active_string $buffer_string $name $scheme]
  puts [ format "selection string: %s " $QM_string ]

  #thickness for the potential embedding region
  set thickness    3.000000000000
  set formal_pQ    4.0000000000000
  set formal_nQ   -2.0000000000000
  set partial_pQ   2.23
  set partial_nQ   -1.115
  set pQ_type  Ti
  set nQ_type  O
  set element  { " " "Ti" "O" }
  set QM_size      20

  if { [ string match $scheme "p" ] == 1 } {
    set core_Q [ expr $partial_nQ - $shell_Q ]
    set pQ $partial_pQ
    set nQ $partial_nQ
  } else {
    set core_Q [ expr $formal_nQ - $shell_Q ]
    set pQ $formal_pQ
    set nQ $formal_nQ
  }

  puts [ format "pQ %g nQ %g core_Q %g shell_Q %g" $pQ $nQ $core_Q $shell_Q]
  
  puts "find all possible sites for QM termination"
  set QMp [ atomselect top [ format "%s and type %s" $QM_string $pQ_type ] ]
  set x_QMp [ $QMp get {x y z} ]
  set type_QMp [ $QMp get type ]
  set n_QMp [ $QMp num ]
  set nbuffer [ atomselect top [ format "%s and (not %s) and type %s" $buffer_string $QM_string $nQ_type ] ]
  set index_nbuffer [ $nbuffer get index ]
  set x_nbuffer [ $nbuffer get {x y z} ]
  set q_nbuffer [ $nbuffer get charge ]
  set type_nbuffer [ $nbuffer get type ]
  set n_nbuffer [ $nbuffer num ]
  set pot_def "index "
  set n_pot 0 
  #for each atom
  for { set i 0 } { $i < $n_nbuffer } { incr i } {
    set is_pot 0
    set type1 [ lindex $type_nbuffer $i ] 
    set c1 [ lindex $x_nbuffer $i ]
    set in [ lindex $index_nbuffer $i ]
    for { set j 0 } { $j < $n_QMp } { incr j } {
      set c2 [ lindex $x_QMp $j ]
      for { set k 0 } { $k < 3 } { incr k } {
        set dk [expr [lindex $c1 $k ] - [lindex $c2 $k]]
        set dk [expr abs($dk)]
        set x$k $dk
      }
      set dx [ expr sqrt($x1*$x1+$x2*$x2+$x0*$x0)]
      #if the distance is smaller than the boundary thickness
      if { $dx < $thickness  } {
        incr is_pot
      }
    }
    if { $is_pot > 0 } {
      set pot_def "$pot_def $in"
      puts "$n_pot: potential atom found $in "
      incr n_pot
    } 
  }

  if { $n_pot > 0 } {
    set region1 [  format "(%s) or (%s)" $QM_string $pot_def  ] 
    set region3 [ format "(%s) and (not (%s)) and (not (%s))" $active_string $QM_string $pot_def ] 
  } else {
    set region1 [  format "(%s) " $QM_string ] 
    set region3 [ format "(%s) and (not (%s)) " $active_string $QM_string  ] 
  }

  puts "region 1: QM"
  set QM [ atomselect top $region1 ]
  set x_QM [ $QM get {x y z} ]
  set type_QM [ $QM get type ]
  set n_QM [ $QM num ]
  puts "region 3: active"
  set MM_active [ atomselect top $region3 ]
  set n_active [ $MM_active num ]
  set x_active [ $MM_active get {x y z} ]
  set type_active [ $MM_active get type ]
  set q_active [ $MM_active get charge ]


  puts "region 4: frozen"
  set region4 [ format "(not (%s)) and (not (%s))" $QM_string $active_string ] 
  set frozen [ atomselect top $region4 ]
  set n_frozen [ $frozen num ]
  set x_frozen [ $frozen get {x y z} ]
  set type_frozen [ $frozen get type ]
  set q_frozen [ $frozen get charge ]

  set nQ_QM [ [ atomselect top [ format "(type %s ) and (%s)" $nQ_type $region1 ] ] num ]
  set pQ_QM [ [ atomselect top [ format "(type %s ) and (%s)" $pQ_type $region1 ] ] num ]

  puts "output basic statistic"
  puts [ format "QM: %g boundary: %g MM_active:%g MM_frozen:%g" $n_QM $n_pot $n_active $n_frozen ]
  puts [ format "CA: %g AN: %g" $pQ_QM $nQ_QM ]

  set name [ format "%s-%g-%g" $name $pQ_QM $nQ_QM ]

  set fo_info [ open [format "%s.info" $name] "w"]
  puts $fo_info [ format "original command: selectQM_link \"%s\" \"%s\" \"%s\" %s %s" $QM_string $active_string $buffer_string  $name $scheme]
  puts $fo_info [ format "selection string: %s " $QM_string ]
  puts $fo_info [ format "pQ %g nQ %g core_Q %g shell_Q %g" $pQ $nQ $core_Q $shell_Q]
  puts $fo_info [ format "QM: %g boundary: %g MM_active:%g MM_frozen:%g" $n_QM $n_pot $n_active $n_frozen ]
  puts $fo_info [ format "CA: %g AN: %g" $pQ_QM $nQ_QM ]

  puts "find all possible sites for link atoms"
  set QMn [ atomselect top [ format "%s and type %s" $QM_string $nQ_type ] ]
  set x_QMn [ $QMn get {x y z} ]
  set type_QMn [ $QMn get type ]
  set n_QMn [ $QMn num ]
  set pbuffer [ atomselect top [ format "%s and (not %s) and type %s" $buffer_string $QM_string $pQ_type ] ]
  set index_pbuffer [ $pbuffer get index ]
  set x_pbuffer [ $pbuffer get {x y z} ]
  set q_pbuffer [ $pbuffer get charge ]
  set type_pbuffer [ $pbuffer get type ]
  set n_pbuffer [ $pbuffer num ]
  set link_def "index "
  set n_link 0 
  #for each atom
  for { set i 0 } { $i < $n_pbuffer } { incr i } {
    set is_link 0
    set type1 [ lindex $type_pbuffer $i ] 
    set c1 [ lindex $x_pbuffer $i ]
    set in [ lindex $index_pbuffer $i ]
    for { set j 0 } { $j < $n_QMn } { incr j } {
      set c2 [ lindex $x_QMn $j ]
      for { set k 0 } { $k < 3 } { incr k } {
        set dk [expr [lindex $c1 $k ] - [lindex $c2 $k]]
        set dk [expr abs($dk)]
        set x$k $dk
      }
      set dx [ expr sqrt($x1*$x1+$x2*$x2+$x0*$x0)]
      #if the distance is smaller than the boundary thickness
      if { $dx < $thickness  } {
        incr in
        puts  "connect [expr $j+1] $in "
        puts  $fo_info " [ expr $j+1 ] $in "
        incr n_link
      }
    }
  }
  puts  $fo_info "$n_link link atoms"

  close $fo_info

  #non-periodic boundary
  set fo_noshell [ open [format "%s_noshell.chm" $name] "w"]
  #with shell predefined
  set fo_shell [ open [format "%s_shell.chm" $name] "w"]
  #for orca testing
  set fo_QM [ open [format "%s.xyz" $name] "w"]
  set fo_MM [ open [format "%s.pc" $name] "w"]
  #force field pre-relaxation
  set fo_ff [ open [format "%s_ff.chm" $name] "w"]


  puts $fo_noshell "c_create coords=${name}_np.pun {"
  puts $fo_noshell "connect ionic\ncoordinates angstrom"
  puts $fo_shell "c_create coords=${name}_s.pun {"
  puts $fo_shell "connect ionic\ncoordinates angstrom"
  puts $fo_QM "*xyz 0 1 " 
  puts $fo_MM [ expr $n_active+$n_frozen ]
  puts $fo_ff "c_create coords=${name}_ff.pun {"
  puts $fo_ff "connect ionic\ncoordinates angstrom"

  puts "output region 1 ..."
  for { set j 0 } { $j < $n_QM } { incr j } {
    set c2 [ lindex $x_QM $j ]
    set x [ lindex $c2 0 ]
    set y [ lindex $c2 1 ]
    set z [ lindex $c2 2 ]
    set ts [ lindex $type_QM $j ]
    puts $fo_noshell [ format "%s1 %g %g %g 0" $ts $x $y $z ]
    puts $fo_shell [ format "%s1 %g %g %g 0" $ts $x $y $z ]
    puts $fo_QM [ format "%s %g %g %g" $ts $x $y $z ]
    if { $ts == $pQ_type } {
      puts $fo_ff [ format "%s1 %g %g %g %g" $ts $x $y $z $pQ ]
    } else {
      puts $fo_ff [ format "%s1 %g %g %g %g" $ts $x $y $z $core_Q ]
    }
  }

  puts "output region 3..."
  for { set j 0 } { $j < $n_active } { incr j } {
    set c2 [ lindex $x_active $j ]
    set x [ lindex $c2 0 ]
    set y [ lindex $c2 1 ]
    set z [ lindex $c2 2 ]
    set q [ lindex $q_active $j ]
    set ts [ lindex $type_active $j ]
    if { [string match $ts $pQ_type ] == 1 } {
      puts $fo_noshell [ format "%s3 %g %g %g %g" $ts $x $y $z $pQ ]
      puts $fo_shell [ format "%s3 %g %g %g %g" $ts $x $y $z $pQ ]
      puts $fo_MM [ format "%g %g %g %g" $pQ $x $y $z ]
      puts $fo_ff [ format "%s3 %g %g %g %g" $ts $x $y $z $pQ ]
    } elseif { [ string match $ts "F" ] == 1 } {
      puts $fo_noshell [ format "%s3 %g %g %g %g" $ts $x $y $z $q ]
      puts $fo_shell [ format "%s3 %g %g %g %g" $ts $x $y $z $q ]
      puts $fo_MM [ format "%g %g %g %g" $q $x $y $z ]
      puts $fo_ff [ format "%s3 %g %g %g %g" $ts $x $y $z $q ]
    } else {
      puts $fo_noshell [ format "%s3 %g %g %g %g" $ts $x $y $z $nQ ]
      puts $fo_shell [ format "%s3 %g %g %g %g" $ts $x $y $z $core_Q ]
      puts $fo_MM [ format "%g %g %g %g" $nQ $x $y $z ]
      puts $fo_ff [ format "%s3 %g %g %g %g" $ts $x $y $z $core_Q ]
    }
  }

  puts "output region 4..."
  for { set j 0 } { $j < $n_frozen } { incr j } {
    set c2 [ lindex $x_frozen $j ]
    set x [ lindex $c2 0 ]
    set y [ lindex $c2 1 ]
    set z [ lindex $c2 2 ]
    set q [ lindex $q_frozen $j ]
    set ts [ lindex $type_frozen $j ]
    puts $fo_noshell [ format "%s4 %g %g %g %g" $ts $x $y $z $q ]
    puts $fo_shell [ format "%s4 %g %g %g %g" $ts $x $y $z $q ]
    puts $fo_MM [ format "%g %g %g %g" $q $x $y $z ]
    puts $fo_ff [ format "%s4 %g %g %g %g" $ts $x $y $z $q ]
  }

  puts $fo_shell "shells"
  puts $fo_ff "shells"

  puts "output region 1 ..."
  for { set j 0 } { $j < $n_QM } { incr j } {
    set c2 [ lindex $x_QM $j ]
    set x [ lindex $c2 0 ]
    set y [ lindex $c2 1 ]
    set z [ lindex $c2 2 ]
    set ts [ lindex $type_QM $j ]
    if { [string match $ts $nQ_type ] == 1 } {
        puts $fo_ff [ format "%s1 %g %g %g %g" $ts $x $y $z $shell_Q ]
    }
  }

  puts "output shell region 3..."
  for { set j 0 } { $j < $n_active } { incr j } {
    set c2 [ lindex $x_active $j ]
    set x [ lindex $c2 0 ]
    set y [ lindex $c2 1 ]
    set z [ lindex $c2 2 ]
    set q [ lindex $q_active $j ]
    set ts [ lindex $type_active $j ]
    if { [string match $ts $nQ_type ] == 1 } {
      puts $fo_shell [ format "%s3 %g %g %g %g" $ts $x $y $z $shell_Q ]
      puts $fo_ff [ format "%s3 %g %g %g %g" $ts $x $y $z $shell_Q ]
    }
  }

  puts $fo_noshell "}"
  puts $fo_noshell "write_xyz coords=${name}_np.pun file=${name}_np_chemshell.xyz"
  puts $fo_shell "}"
  puts $fo_shell "write_xyz coords=${name}_s.pun file=${name}_s_chemshell.xyz"
  puts $fo_ff "}"
  puts $fo_ff "write_xyz coords=${name}_ff.pun file=${name}_chemff.xyz"
  puts $fo_QM "*"

  close $fo_noshell
  close $fo_shell
  close $fo_QM
  close $fo_MM
  close $fo_ff


  mol delrep 0 top

  puts "visualize region 3"
  mol color name
  mol representation CPK 1.00 0.000000 32.000000 12.000000
  mol material Opaque
  mol selection $region3
  mol addrep top

  puts "visualize region 4"
  mol material Transparent
  mol representation CPK 0.500 0.000000 32.000000 12.000000
  mol selection $region4 
  mol addrep top

  puts "visualize region 1"
  mol material Opaque
  mol representation CPK 3.000000 0.000000 32.000000 12.000000
  mol selection $region1
  mol addrep top

}
