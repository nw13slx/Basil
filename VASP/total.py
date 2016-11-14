#!/usr/bin/env python

#usage: split_dos_atom argument-name argument
#argument name, parameter
# partial orbital: "s, p, d, f", # of element; i.e. p 1 f 0
# input element species: ntype, [# of types, type 0 #, type 1 #, ...]
#  this setting overrides the one in poscar/contcar
# write_DOS0 yes/no
# write_PDOS yes/no
# centerEf yes/no
# zoomin emin emax

#this file is adapted from vtst tool
#the original source code is buggy and kind of slow (especially when you have >100 atoms!)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl
import sys
import os

#default setting
write_DOS0=True
write_PDOS=True
peratom=False
centerEf=True
zoomIn=False

#initialization
ntype=0
energyshift=0
par_element=[-1,-1,-1,-1] #partial orbital
atomSpe=None
position=None
symbol=None
name=""
zoomEmax=-1
zoomEmin=-1

#main function
def main():

  import io
  from parse import parseInput
  end=parseInput(sys.argv)
  if end:
    print "end condition: ",end
    return

  f = open("DOSCAR", 'r')
  #read the header of DOSCAR

  io.dosf=f
  io.read_dosfile()
  io.read_dos0()

  #pdos=io.read_pdos()
  #if (pdos==None):
  #  io.read_parorbit()
