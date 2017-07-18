set count 0
foreach i {2 4 5} {
  foreach j {"sub" "top"} {
    set name [format "./110%s-%g.poscar" $j $i]
    set output [format "110%s-%g" $j $i]
    mol new $name type POSCAR first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
    pbc box -color black
    mol delrep 0 top
    mol representation CPK 1.000000 0.000000 32.000000 12.000000
    mol material BrushedMetal
    mol selection {z>20 and type O}
    mol addrep top
    mol selection {z>20 and type Ce}
    mol addrep top
    mol selection {z>20 and type La}
    mol addrep top
    mol material Transparent
    mol selection {z<20 and z>19 and type O}
    mol addrep top
    mol selection {z<20 and z>19 and type Ce}
    mol addrep top
    mol selection {z<20 and z>19 and type La}
    mol addrep top
    render Tachyon $output "/usr/local/lib/vmd/tachyon_LINUXAMD64 -aasamples 12 %s -format TARGA -o %s.tga"
    mol delete $count
    pbc box -off
    incr count
  }
}
