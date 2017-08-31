#!/bin/env python
from control import *
from iodos import *
from data import *
from plot import *
from analysis import *

import sys
def main():
  #initialization
  atoms=atom_data()
  dos=dos_data()
  cont=control_knob(sys.argv,atoms,dos)
  plt=plot(cont)
  ana=analysis(atoms,dos)

  dumpclean(cont.__dict__)

  #read in overall dos
  io=iodos(cont,plt)
  io.read_poscar()
  io.read_tot_dosfile()
  plt.tot_dos()

  #write
  if cont.write_dos0:
      io.write_tot_dosfile()
  io.delete_tot_dosfile()

  #read in partial dos
  if cont.run_pdos:
    io.read_pdos()
    plt.pdos()
    if cont.write_pdos:
        io.write_pdosfile()
    io.delete_pdosfile()

  #debug
  #dumpclean(io.__dict__)
  #dumpclean(atoms.__dict__)
  #dumpclean(dos.__dict__)

if __name__ == '__main__':
  main()
