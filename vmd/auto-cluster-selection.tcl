#requirement: type should be element symbol instead of number. point charge should be labelled as "F"
#             the configuration does not have periodic boundary
#             needs to have "q" "charge" data

#argument 1: center atom
#argument 2: inner radius
#argument 3: buffer thickness
#argument 4: temperature
#argument 5: number of attempts

#example: autoselectQM 950 7.5 O 3 1000 5 10000

#autoselectQM2 offer a more close shell option. the initial inner cluster is cleaned up a bit, where dangling atoms are removed.

set K 167101.001981787

set hash_coul234 {}
set hash_coul14  0

proc objective { centerid rmin t_ele buffer QMi sel sel_id unsel rest } {

  ##case 1, use self coulombic energy
  #return [ expr [ coulomb $sel $sel ] + [ coulomb_cross $QMi $sel ] ]

  ##case 2, interaction energy
  #global hash_coul234
  #global hash_coul14
  #if { [ dict size $hash_coul ] == 0 } {
  #  set hash_coul [ hash_coulomb_cross $rest $QMo ]
  #  set hash_coul14 [ coulomb_cross $QMi $rest ]
  #}
  #return [ expr [ sum_coulomb_cross $hash_coul234 $sel_id ] + [ coulomb_cross $QMi $unsel ] + [ coulomb_cross $sel $unsel] +$hash_coul14 + 
  
  ##case 3, total coordination
  #return [ count_coordination "index $sel_id $QM1_sel $QM2_sel "] ]
  
  #case 4, 2nd order momentum
  return [ moment2nd $centerid $sel ]
}

proc autoselectQM { centerid rmin t_ele buffer T attempt} {

  global K

  #set up the inner sphere
  set QMi [ atomselect top "within $rmin of index $centerid" ]
  set n_QMi [ $QMi num ]
  set Q_QMi [ ladd [ $QMi get charge ] ]
  set type_QMi [ $QMi get type ]
  set counters [ count_ele $type_QMi ]
  set stat ""
  dict for {item count} $counters {
    set stat  "$stat${item}$count"
  }
  puts "$stat netcharge $Q_QMi"

  set fo_info [ open [format "$stat-$T.1autoinfo" ] "w"]
  set fo_config [ open [format "$stat-$T.1config" ] "w"]
  #repeats the argument
  puts "original command: autoselectQM $centerid $rmin $t_ele $buffer $T $attempt "
  puts $fo_info "original command: autoselectQM $centerid $rmin $t_ele $buffer $T $attempt "
  puts $fo_info "$stat netcharge $Q_QMi"

  #set outer sphere
  set rmax [ expr $rmin+$buffer ]
  set QMo [ atomselect top "type $t_ele and (not within $rmin of index $centerid ) and within $buffer of (within $rmin of index $centerid)" ]
  set n_QMo [ $QMo num ]
  set id_QMo [ $QMo get index ]
  set rest [ atomselect top " (not (index $id_QMo or (within $rmin of index $centerid))) and within 10 of index $centerid " ]

  #check charge and find out number of atoms to add
  set q_QMo [ $QMo get charge ]
  set Q_QMo [ ladd $q_QMo ]
  set countq [ count_ele $q_QMo ]
  if { [ dict size $countq ] != 1 } {
    puts "ERROR, there are more than one type of charges for type $t_ele "
    puts $fo_info "ERROR, there are more than one type of charges for type $t_ele "
    return 1
  }
  set shellq [ lindex $q_QMo 0 ]
  set nshell [ expr int(-$Q_QMi/$shellq) ]
  if { $nshell <0 } {
    puts "ERROR, no way to balance the charge, change the rmin or t_ele"
    puts $fo_info "ERROR, no way to balance the charge, change the rmin or t_ele"
    return 1
  } elseif { $nshell == 0 } {
    puts "Hooray, automatically charge compensate... "
    puts $fo_info "Hooray, automatically charge compensate... "
    return 0
  } else {
    puts "need to choose $nshell from $n_QMo "
    puts $fo_info  "need to choose $nshell from $n_QMo "
    if { $nshell < 2 } {
      set attempt [ expr $n_QMo * ($n_QMo-1) ]
    }
  }

  if { $nshell == 2 } {

    set minpot 0
    set minsel {}
    set id_QMo [ $QMo get index ]
    for { set i 0 } { $i < $n_QMo } {incr i} {
      for { set j [ expr $i+1] } { $j < $n_QMo } {incr j} {
        set sel_id [ list [lindex $id_QMo $i] [lindex $id_QMo $j] ]
        set select [ atomselect top "index $sel_id "]
        set unselect [ atomselect top "index $id_QMo and not index $sel_id" ] 
        set pot [ objective $centerid $rmin $t_ele $buffer $QMi $select $sel_id $unselect $rest ]
        if { $minpot >$pot } {
          set minpot $pot
          set minsel $sel_id
        }
      }
    }

    puts "minimum potential $minpot "
    puts "(index $minsel ) or (within $rmin of index $centerid)"
    puts $fo_info "final potential $minpot "
    puts $fo_info "(index $minsel ) or (within $rmin of index $centerid)"

    close $fo_info
    close $fo_config
    return
  }

  if { $nshell == 1 } {

    set minpot 0
    set minsel {}
    set id_QMo [ $QMo get index ]
    for { set i 0 } { $i < $n_QMo } {incr i} {
      set sel_id [ list [lindex $id_QMo $i] ]
      set select [ atomselect top "index $sel_id "]
      set unselect [ atomselect top "index $id_QMo and not index $sel_id" ] 
      set pot [ objective $centerid $rmin $t_ele $buffer $QMi $select $sel_id $unselect $rest ]
      if { $minpot >$pot } {
        set minpot $pot
        set minsel $sel_id
      }
    }

    puts "minimum potential $minpot "
    puts "(index $minsel ) or (within $rmin of index $centerid)"
    puts $fo_info "final potential $minpot "
    puts $fo_info "(index $minsel ) or (within $rmin of index $centerid)"

    close $fo_info
    close $fo_config
    return
  }

  #monte carlo sampling
  #output minimum and final guess

  #prepare initial guess
  set rand [ expr srand(1356) ]
  set sel_id {}
  set i 0
  set tmp 0
  set unsel_id [ $QMo get index ]
  set n_uns [ llength $unsel_id ]
  while { $i < $nshell && $tmp < [ expr ${nshell}*10] } {
    set pos [ expr int(floor($n_uns*rand())) ] 
    set rand [ lindex $unsel_id $pos ]
    lappend sel_id $rand
    set unsel_id [ lreplace $unsel_id $pos $pos ]
    set n_uns [ llength $unsel_id ]
    incr i
    incr tmp
  }


  set select [ atomselect top "index $sel_id "]
  set unselect [ atomselect top "index $unsel_id "]
  set sel_id0 $sel_id
  set pot [ objective $centerid $rmin $t_ele $buffer $QMi $select $sel_id $unselect $rest ]
  puts "initial gues: $pot "
  puts $fo_info "initial gues: $pot "
 
  set naccept 0
  set minpot $pot
  set minsel $sel_id
  set tmp 0
  while { $tmp < $attempt } {
    #swap
    set pos1 [ expr int(floor([llength $sel_id]*rand())) ] 
    set pos2 [ expr int(floor([llength $unsel_id]*rand())) ] 
    set ele1 [ lindex $sel_id $pos1 ]
    set ele2 [ lindex $unsel_id $pos2 ]
    set sel_id [ lreplace $sel_id $pos1 $pos1 ]
    set unsel_id [ lreplace $unsel_id $pos2 $pos2 ]
    lappend sel_id $ele2
    lappend unsel_id $ele1

    set select [ atomselect top "index $sel_id "]
    set unselect [ atomselect top "index $unsel_id "]
    set newpot [ objective $centerid $rmin $t_ele $buffer $QMi $select $sel_id $unselect $rest ]
    puts "initial gues: $pot "

    set accept 0

    if { $newpot < $pot } {
      set accept 1
    } else {
      set delta [ expr ($newpot - $pot)*$K/$T ]
      set exp [ expr exp(-$delta ) ]
      set p [ expr rand() ]
      if { $p < $exp } {
        set accept 1
        #puts "roll a dice $p < $exp ,$delta"
        #puts $fo_info "roll a dice $p < $exp ,$delta"
      } 
    }

    if { $accept == 1 } {
      set pot $newpot
      set sel_id0 $sel_id
      if { $minpot >$pot } {
        set minpot $pot
        set minsel $sel_id
      }
      incr naccept
    } else { 
      #unswap
      set pos1 [ expr [llength $sel_id] -1 ]
      set pos2 [ expr [llength $unsel_id] -1 ]
      set sel_id [ lreplace $sel_id $pos1 $pos1 ]
      set unsel_id [ lreplace $unsel_id $pos2 $pos2 ]
      lappend sel_id $ele1
      lappend unsel_id $ele2
    }

    if { $tmp % 100 == 0 } {
      puts "$tmp $pot $minpot "
      puts $fo_info "$tmp $pot $minpot"
    }

    incr tmp
  }

  puts "accepting ratio $naccept $tmp [ expr $naccept/double($tmp) ]"
  puts $fo_info "accepting ratio [ expr $naccept/double($tmp) ]"

  puts "final potential $pot "
  puts "index $sel_id or (within $rmin of index $centerid)" 
  puts $fo_info "final potential $pot "
  puts $fo_info "index $sel_id or (within $rmin of index $centerid)" 

  puts "minimum potential $minpot "
  puts "index $minsel or (within $rmin of index $centerid)" 
  puts $fo_info "final potential $minpot "
  puts $fo_info "index $minsel or (within $rmin of index $centerid)" 

  close $fo_info
  close $fo_config
}
  

proc autoselectQM2 { centerid rmin t_ele buffer T attempt } {

  global K
  global debug

  #first label all cation
  set QM1 [ atomselect top "not type $t_ele and within $rmin of index $centerid" ]
  set QM2 [ atomselect top "type $t_ele and within $rmin of index $centerid" ]
  set id_QM1 [ $QM1 get index ]
  set id_QM2 [ $QM2 get index ]
  set n_QM1 [ $QM1 num ]
  set n_QM2 [ $QM2 num ]

  #remove the dangling oxygen
  set QM2_sel {}
  foreach id2 $id_QM2 {
    set ngh [ atomselect top "(not type $t_ele) and (within $rmin of index $centerid) and (within $buffer of index $id2)" ]
    set n_ngh [ $ngh num ]
    if { $n_ngh >1  } {
      lappend QM2_sel $id2
    }
  }

  #remove the dangling cation
  set QM1_sel {}
  foreach id1 $id_QM1 {
    set ngh [ atomselect top "(index $QM2_sel) and (within $buffer of index $id1)" ]
    set n_ngh [ $ngh num ]
    if { $n_ngh >1  } {
      lappend QM1_sel $id1
    }
  }

  #set up the inner sphere
  set QMi [ atomselect top "index $QM2_sel $QM1_sel" ]
  set n_QMi [ $QMi num ]
  set Q_QMi [ ladd [ $QMi get charge ] ]
  set type_QMi [ $QMi get type ]
  set counters [ count_ele $type_QMi ]
  set stat ""
  dict for {item count} $counters {
    set stat  "$stat${item}$count"
  }
  puts "$stat netcharge $Q_QMi"

  set fo_info [ open [format "$stat-$T.autoinfo" ] "w"]
  set fo_config [ open [format "$stat-$T.config" ] "w"]
  #repeats the argument
  puts "original command: autoselectQM2 $centerid $rmin $t_ele $buffer $T $attempt "
  puts $fo_info "original command: autoselectQM2 $centerid $rmin $t_ele $buffer $T $attempt "
  puts $fo_info "$stat netcharge $Q_QMi"

  #set outer sphere
  set rmax [ expr $rmin+$buffer ]
  set QMo [ atomselect top "type $t_ele and (not index $QM2_sel $QM1_sel ) and (within $buffer of (index $QM2_sel $QM1_sel))" ]
  set n_QMo [ $QMo num ]
  set id_QMo [ $QMo get index ]
  set rest [ atomselect top " (not index $id_QMo $QM1_sel $QM2_sel) and within 10 of index $centerid " ]
  puts "known atoms $QM2_sel $QM1_sel "
  puts "select atoms from $id_QMo"
  puts $fo_info "known atoms $QM2_sel $QM1_sel "
  puts $fo_info "select atoms from $id_QMo"

  #check charge and find out number of atoms to add
  set q_QMo [ $QMo get charge ]
  set Q_QMo [ ladd $q_QMo ]
  set countq [ count_ele $q_QMo ]
  if { [ dict size $countq ] != 1 } {
    puts "ERROR, there are more than one type of charges for type $t_ele "
    puts $fo_info "ERROR, there are more than one type of charges for type $t_ele "
    return 1
  }
  set shellq [ lindex $q_QMo 0 ]
  set nshell [ expr int(-$Q_QMi/$shellq) ]
  if { $nshell <0 } {
    puts "ERROR, no way to balance the charge, change the rmin or t_ele"
    puts $fo_info "ERROR, no way to balance the charge, change the rmin or t_ele"
    return 1
  } elseif { $nshell == 0 } {
    puts "Hooray, automatically charge compensate... "
    puts $fo_info "Hooray, automatically charge compensate... "
    return 0
  } else {
    puts "need to choose $nshell from $n_QMo "
    puts $fo_info  "need to choose $nshell from $n_QMo "
    if { $nshell < 2 } {
      set attempt [ expr $n_QMo * ($n_QMo-1) ]
    }
  }

  
  if { $nshell == 2 } {

    set minpot 0
    set minsel {}
    set id_QMo [ $QMo get index ]
    for { set i 0 } { $i < $n_QMo } {incr i} {
      for { set j [ expr $i+1] } { $j < $n_QMo } {incr j} {
        set sel_id [ list [lindex $id_QMo $i] [lindex $id_QMo $j] ]
        set select [ atomselect top "index $sel_id" ]
        set unselect [ atomselect top "index $id_QMo and not (index $sel_id)" ]
        set pot [ objective $centerid $rmin $t_ele $buffer $QMi $select $sel_id $unselect $rest ]
        if { $minpot >$pot } {
          set minpot $pot
          set minsel $sel_id
        }
      }
    }

    puts "minimum potential $minpot "
    puts "index $minsel $QM2_sel $QM1_sel" 
    puts $fo_info "final potential $minpot "
    puts $fo_info "index $minsel $QM2_sel $QM1_sel" 

    close $fo_info
    close $fo_config
    return
  }

  if { $nshell == 1 } {

    set minpot 0
    set minsel {}
    set id_QMo [ $QMo get index ]
    for { set i 0 } { $i < $n_QMo } {incr i} {
      set sel_id [ list [lindex $id_QMo $i] ]
      set select [ atomselect top "index $sel_id" ]
      set unselect [ atomselect top "index $id_QMo and not (index $sel_id)" ]
      set pot [ objective $centerid $rmin $t_ele $buffer $QMi $select $sel_id $unselect $rest ]
      if { $minpot >$pot } {
        set minpot $pot
        set minsel $sel_id
      }
    }

    puts "minimum potential $minpot "
    puts "index $minsel $QM2_sel $QM1_sel" 
    puts $fo_info "final potential $minpot "
    puts $fo_info "index $minsel $QM2_sel $QM1_sel" 

    close $fo_info
    close $fo_config
    return
  }

  #monte carlo sampling
  #output minimum and final guess

  ##input initial guess
  #set sel_id [ [atomselect top "($initial) and (not index $QM1_sel $QM2_sel)" ] get index ]
  #puts $sel_id
  #set unsel_id [ [atomselect top "(index $id_QMo) and (not ($initial))" ] get index ]
  #puts $unsel_id

  #prepare initial guess
  set rand [ expr srand(1356) ]
  set sel_id {}
  set i 0
  set tmp 0
  set unsel_id [ $QMo get index ]
  set n_uns [ llength $unsel_id ]
  while { $i < $nshell && $tmp < [ expr ${nshell}*10] } {
    set pos [ expr int(floor($n_uns*rand())) ] 
    set rand [ lindex $unsel_id $pos ]
    lappend sel_id $rand
    set unsel_id [ lreplace $unsel_id $pos $pos ]
    set n_uns [ llength $unsel_id ]
    incr i
    incr tmp
  }


  set sel_id0 $sel_id
  set select   [ atomselect top "index $sel_id" ]
  set unselect [ atomselect top "index $unsel_id" ]
  set pot [ objective $centerid $rmin $t_ele $buffer $QMi $select $sel_id $unselect $rest ]

  puts "initial gues: $pot "
  puts $fo_info "initial gues: $pot "

  set naccept 0
  set minpot $pot
  set minsel $sel_id
  set alleng {}
  set tmp 0
  while { $tmp < $attempt } {
    #swap
    set pos1 [ expr int(floor([llength $sel_id]*rand())) ] 
    set pos2 [ expr int(floor([llength $unsel_id]*rand())) ] 
    set ele1 [ lindex $sel_id $pos1 ]
    set ele2 [ lindex $unsel_id $pos2 ]
    set sel_id [ lreplace $sel_id $pos1 $pos1 ]
    set unsel_id [ lreplace $unsel_id $pos2 $pos2 ]
    lappend sel_id $ele2
    lappend unsel_id $ele1

    set select   [ atomselect top "index $sel_id" ]
    set unselect [ atomselect top "index $unsel_id" ]
    set newpot [ objective $centerid $rmin $t_ele $buffer $QMi $select $sel_id $unselect $rest ]

    set accept 0

    if { $newpot < $pot } {
      set accept 1
    } else {
      set delta [ expr ($newpot - $pot)*$K/$T ]
      set exp [ expr exp(-$delta ) ]
      set p [ expr rand() ]
      if { $p < $exp } {
        set accept 1
        #puts "roll a dice $p < $exp ,$delta"
        #puts $fo_info "roll a dice $p < $exp ,$delta"
      } 
    }

    if { $debug == "t" } {
      if { [ dict exists $alleng $newpot ] < 1 } {
        dict incr alleng $newpot
        puts "$tmp $newpot index $QM1_sel $QM2_sel $sel_id"
        puts $fo_config "$tmp $newpot index $QM2_sel $QM1_sel $sel_id"
      }
    }

    if { $accept == 1 } {
      set pot $newpot
      set sel_id0 $sel_id
      if { $minpot >$pot } {
        set minpot $pot
        set minsel $sel_id
      }
      incr naccept
    } else { 
      #unswap
      set pos1 [ expr [llength $sel_id] -1 ]
      set pos2 [ expr [llength $unsel_id] -1 ]
      set sel_id [ lreplace $sel_id $pos1 $pos1 ]
      set unsel_id [ lreplace $unsel_id $pos2 $pos2 ]
      lappend sel_id $ele1
      lappend unsel_id $ele2
    }

    if { $tmp % 100 == 0 } {
      puts "$tmp $minpot $pot"
      puts $fo_info "$tmp $minpot $pot"
    }


    incr tmp
  }

  puts "accepting ratio $naccept $tmp [ expr $naccept/double($tmp) ]"
  puts $fo_info "accepting ratio [ expr $naccept/double($tmp) ]"

  puts "final potential $pot "
  puts "index $sel_id0 $QM2_sel $QM1_sel" 
  puts $fo_info "final potential $pot "
  puts $fo_info "index $sel_id0 $QM2_sel $QM1_sel" 

  puts "minimum potential $minpot "
  puts "index $minsel $QM2_sel $QM1_sel" 
  puts $fo_info "minimum potential $minpot "
  puts $fo_info "index $minsel $QM2_sel $QM1_sel" 

  close $fo_info
  close $fo_config
}

proc autoselectQM3 { centerid QM1_def t_ele buffer coord T attempt } {

  global K
  global debug

  set QM1 [ atomselect top "not type $t_ele and ($QM1_def)" ]
  set QM1_sel [ $QM1 get index ]
  set n_QM1 [ $QM1 num ]
  puts "label all cation $QM1_sel"

  set QM2 [ atomselect top "type $t_ele and within $buffer of (index $QM1_sel)" ]
  set id_QM2 [ $QM2 get index ]
  set n_QM2 [ $QM2 num ]

  set QM2_sel {}
  for { set i 0 } { $i < $n_QM2 } { incr i } {
    set id [ lindex $id_QM2 $i ]
    set ngh [ atomselect top "(not type $t_ele) and (index $QM1_sel) and (within $buffer of index $id)" ]
    set n_ngh [ $ngh num ]
    if { $n_ngh >= $coord  } {
      lappend QM2_sel $id
    }
  }
  puts "find all the linked anion $QM2_sel"

  #set up the inner sphere
  set QMi [ atomselect top "index $QM2_sel $QM1_sel" ]
  set n_QMi [ $QMi num ]
  set Q_QMi [ ladd [ $QMi get charge ] ]
  set type_QMi [ $QMi get type ]
  set counters [ count_ele $type_QMi ]
  set stat ""
  dict for {item count} $counters {
    set stat  "$stat${item}$count"
  }
  puts "$stat netcharge $Q_QMi"

  set fo_info [ open [format "$stat-$T.autoinfo" ] "w"]
  set fo_config [ open [format "$stat-$T.config" ] "w"]
  #repeats the argument
  puts "original command: autoselectQM2 $centerid $t_ele $buffer $T $attempt "
  puts $fo_info "original command: autoselectQM2 $centerid $t_ele $buffer $T $attempt "
  puts $fo_info "$stat netcharge $Q_QMi"

  #set outer sphere
  set QMo [ atomselect top "type $t_ele and (not index $QM2_sel $QM1_sel ) and (within $buffer of (index $QM2_sel $QM1_sel))" ]
  set n_QMo [ $QMo num ]
  set id_QMo [ $QMo get index ]
  puts "known atoms $QM2_sel $QM1_sel "
  puts "select atoms from $id_QMo"
  puts $fo_info "known atoms $QM2_sel $QM1_sel "
  puts $fo_info "select atoms from $id_QMo"

  #check charge and find out number of atoms to add
  set q_QMo [ $QMo get charge ]
  set Q_QMo [ ladd $q_QMo ]
  set countq [ count_ele $q_QMo ]
  if { [ dict size $countq ] != 1 } {
    puts "ERROR, there are more than one type of charges for type $t_ele "
    puts $fo_info "ERROR, there are more than one type of charges for type $t_ele "
    return 1
  }
  set shellq [ lindex $q_QMo 0 ]
  set nshell [ expr int(-$Q_QMi/$shellq) ]
  if { $nshell <0 } {
    puts "ERROR, no way to balance the charge, change t_ele"
    puts $fo_info "ERROR, no way to balance the charge, change t_ele"
    return 1
  } elseif { $nshell == 0 } {
    puts "Hooray, automatically charge compensate... "
    puts $fo_info "Hooray, automatically charge compensate... "
    return 0
  } else {
    puts "need to choose $nshell from $n_QMo "
    puts $fo_info  "need to choose $nshell from $n_QMo "
    if { $nshell < 2 } {
      set attempt [ expr $n_QMo * ($n_QMo-1) ]
    }
  }

  set rmin 0
  set rest [ atomselect top " (not index $id_QMo $QM1_sel $QM2_sel) and within 10 of index $centerid " ]
  
  if { $nshell == 2 } {

    set minpot 0
    set minsel {}
    set id_QMo [ $QMo get index ]
    for { set i 0 } { $i < $n_QMo } {incr i} {
      for { set j [ expr $i+1] } { $j < $n_QMo } {incr j} {
        set sel_id [ list [lindex $id_QMo $i] [lindex $id_QMo $j] ]
        set select   [ atomselect top "index $sel_id" ]
        set unselect [ atomselect top "index $id_QMo and not index $sel_id" ]
        set pot [ objective $centerid $rmin $t_ele $buffer $QMi $select $sel_id $unselect $rest ]
        if { $minpot >$pot } {
          set minpot $pot
          set minsel $sel_id
        }
      }
    }

    puts "minimum potential $minpot "
    puts "index $minsel $QM2_sel $QM1_sel" 
    puts $fo_info "final potential $minpot "
    puts $fo_info "index $minsel $QM2_sel $QM1_sel" 

    close $fo_info
    close $fo_config
    return
  }

  if { $nshell == 1 } {

    set minpot 0
    set minsel {}
    set id_QMo [ $QMo get index ]
    for { set i 0 } { $i < $n_QMo } {incr i} {
      set sel_id [ list [lindex $id_QMo $i] ]
      set select   [ atomselect top "index $sel_id" ]
      set unselect [ atomselect top "index $id_QMo and not index $sel_id" ]
      set pot [ objective $centerid $rmin $t_ele $buffer $QMi $select $sel_id $unselect $rest ]
      if { $minpot >$pot } {
        set minpot $pot
        set minsel $sel_id
      }
    }

    puts "minimum potential $minpot "
    puts "index $minsel $QM2_sel $QM1_sel" 
    puts $fo_info "final potential $minpot "
    puts $fo_info "index $minsel $QM2_sel $QM1_sel" 

    close $fo_info
    close $fo_config
    return
  }

  #monte carlo sampling
  #output minimum and final guess

  ##input initial guess
  #set sel_id [ [atomselect top "($initial) and (not index $QM1_sel $QM2_sel)" ] get index ]
  #puts $sel_id
  #set unsel_id [ [atomselect top "(index $id_QMo) and (not ($initial))" ] get index ]
  #puts $unsel_id

  #prepare initial guess
  set rand [ expr srand(1356) ]
  set sel_id {}
  set i 0
  set tmp 0
  set unsel_id [ $QMo get index ]
  set n_uns [ llength $unsel_id ]
  while { $i < $nshell && $tmp < [ expr ${nshell}*10] } {
    set pos [ expr int(floor($n_uns*rand())) ] 
    set rand [ lindex $unsel_id $pos ]
    lappend sel_id $rand
    set unsel_id [ lreplace $unsel_id $pos $pos ]
    set n_uns [ llength $unsel_id ]
    incr i
    incr tmp
  }


  set sel_id0 $sel_id
  set select   [ atomselect top "index $sel_id" ]
  set unselect [ atomselect top "index $unsel_id" ]
  set rest [ atomselect top " (not index $id_QMo $QM1_sel $QM2_sel) and within 10 of index $centerid " ]
  set rmin 0
  set pot [ objective $centerid $rmin $t_ele $buffer $QMi $select $sel_id $unselect $rest ]

  puts "initial gues: $pot "
  puts $fo_info "initial gues: $pot "

  set naccept 0
  set minpot $pot
  set minsel $sel_id
  set alleng {}
  set tmp 0
  while { $tmp < $attempt } {
    #swap
    set pos1 [ expr int(floor([llength $sel_id]*rand())) ] 
    set pos2 [ expr int(floor([llength $unsel_id]*rand())) ] 
    set ele1 [ lindex $sel_id $pos1 ]
    set ele2 [ lindex $unsel_id $pos2 ]
    set sel_id [ lreplace $sel_id $pos1 $pos1 ]
    set unsel_id [ lreplace $unsel_id $pos2 $pos2 ]
    lappend sel_id $ele2
    lappend unsel_id $ele1

    set select   [ atomselect top "index $sel_id" ]
    set unselect [ atomselect top "index $unsel_id" ]
    set newpot [ objective $centerid $rmin $t_ele $buffer $QMi $select $sel_id $unselect $rest ]

    set accept 0

    if { $newpot < $pot } {
      set accept 1
    } else {
      set delta [ expr ($newpot - $pot)*$K/$T ]
      set exp [ expr exp(-$delta ) ]
      set p [ expr rand() ]
      if { $p < $exp } {
        set accept 1
        #puts "roll a dice $p < $exp ,$delta"
        #puts $fo_info "roll a dice $p < $exp ,$delta"
      } 
    }

    if { $debug == "t" } {
      if { [ dict exists $alleng $newpot ] < 1 } {
        dict incr alleng $newpot
        puts "$tmp $newpot index $QM1_sel $QM2_sel $sel_id"
        puts $fo_config "$tmp $newpot index $QM2_sel $QM1_sel $sel_id"
      }
    }

    if { $accept == 1 } {
      set pot $newpot
      set sel_id0 $sel_id
      if { $minpot >$pot } {
        set minpot $pot
        set minsel $sel_id
      }
      incr naccept
    } else { 
      #unswap
      set pos1 [ expr [llength $sel_id] -1 ]
      set pos2 [ expr [llength $unsel_id] -1 ]
      set sel_id [ lreplace $sel_id $pos1 $pos1 ]
      set unsel_id [ lreplace $unsel_id $pos2 $pos2 ]
      lappend sel_id $ele1
      lappend unsel_id $ele2
    }

    if { $tmp % 100 == 0 } {
      puts "$tmp $minpot $pot"
      puts $fo_info "$tmp $minpot $pot"
    }


    incr tmp
  }

  puts "accepting ratio $naccept $tmp [ expr $naccept/double($tmp) ]"
  puts $fo_info "accepting ratio [ expr $naccept/double($tmp) ]"

  puts "final potential $pot "
  puts "index $sel_id0 $QM2_sel $QM1_sel" 
  puts $fo_info "final potential $pot "
  puts $fo_info "index $sel_id0 $QM2_sel $QM1_sel" 

  puts "minimum potential $minpot "
  puts "index $minsel $QM2_sel $QM1_sel" 
  puts $fo_info "minimum potential $minpot "
  puts $fo_info "index $minsel $QM2_sel $QM1_sel" 

  close $fo_info
  close $fo_config
}

#}
#
#
