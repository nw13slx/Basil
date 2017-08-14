a=5.5055
xx=$(echo "scale=8;$a/sqrt(2)" | bc)

cat > 110.lammpstrj <<EOF
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
6
ITEM: BOX BOUNDS pp pp pp
0 $xx
0 $xx
0 $a
ITEM: ATOMS id type xs ys zs
1   1   0.25  0.25  0.25
2   2   0.25  0.75  0.499
3   2   0.25  0.75  0.000001
4   1   0.75  0.75  0.75
5   2   0.75  0.25  0.501
6   2   0.75  0.25  0.999999
EOF

atomsk -expand 10 10 10 110.lammpstrj a.lmp
