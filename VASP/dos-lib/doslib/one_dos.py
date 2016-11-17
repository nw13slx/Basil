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

  #read in overall dos
  io=iodos(cont,plt)
  io.read_poscar()
  io.read_tot_dosfile()
  ana.bandgap()

  #write
  if cont.write_dos0:
      io.write_tot_dosfile()
  io.delete_tot_dosfile()

  #read in partial dos
  io.read_pdos()
  plt.pdos()
  if cont.write_pdos:
      io.write_pdosfile()
  io.delete_pdosfile()

  #debug
  #dumpclean(cont.__dict__)
  #dumpclean(io.__dict__)
  #dumpclean(atoms.__dict__)
  #dumpclean(dos.__dict__)

if __name__ == '__main__':
  main()
