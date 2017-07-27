#!/usr/bin/env python

#usage: ./converlq.py <lammps data input> <output file name> element1 element2 element3
import numpy as np
import sys

f_in=sys.argv[1]
f_out=sys.argv[2]
name=sys.argv[3:]

lmp=open(f_in,"r")

line=lmp.readline()
line=lmp.readline()

natom=int(lmp.readline().strip().split()[0])
ntype=int(lmp.readline().strip().split()[0])

line=lmp.readline()

lx=map(float,lmp.readline().strip().split()[0:2])
ly=map(float,lmp.readline().strip().split()[0:2])
lz=map(float,lmp.readline().strip().split()[0:2])
alat=lx[1]-lx[0]
cell2=(ly[1]-ly[0])/alat
cell3=(lz[1]-lz[0])/alat

#mass
line=lmp.readline()
line=lmp.readline()
line=lmp.readline()
for i in range(ntype):
  line=lmp.readline()

#poten
line=lmp.readline()
line=lmp.readline()
line=lmp.readline()
for i in range(ntype):
  line=lmp.readline()

line=lmp.readline()
line=lmp.readline()
line=lmp.readline()

f = open(f_out, 'w')
print >> f, natom, "\n", ntype, "\n",alat,"\n",cell2, "\n",cell3

line_atom=[]
line_type=[]
for i in range(natom):
    line=lmp.readline().strip().split()
    print >> f, name[int(line[1])-1], (float(line[3])-lx[0])/alat, (float(line[4])-ly[0])/alat,(float(line[5])-lz[0])/alat
