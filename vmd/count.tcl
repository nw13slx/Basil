proc ladd L {expr [join $L +]+0} ;

proc count {gstring} {
  set group [ atomselect top $gstring ]
  puts [ $group num ]
}

