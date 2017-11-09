source addbond.tcl
mol delrep 0 top

pbc box -color black
rotate y by 90
topo clearbonds
addbond 3 4 0 3
addbond 2 3 0 3

mol representation CPK 1.000000 0.000000 32.000000 12.000000
mol material BrushedMetal
mol selection {z>4}
mol addrep top

mol representation CPK 1.000000 0.500000 32.000000 12.000000
mol material Transparent
mol selection {z<=4}
mol addrep top

mol representation CPK 0.000000 0.500000 32.000000 12.000000
mol material BrushedMetal
mol selection {type 3 4}
mol addrep top

#mol delete $count
#pbc box -off
