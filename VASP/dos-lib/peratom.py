from doslib.control import *
from doslib.iodos import *
from doslib.data import *
from doslib.plot import *
from doslib.analysis import *

import sys

#here, define your own per-atom plot function, as in slab-plot.py
import slab_plot as sp
plot.atom_start=sp.plot_atom_start
plot.atom=sp.plot_atom
plot.atom_end=sp.plot_atom_end

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
  ana.peak_finder(-32.5,5)

  #write
  if cont.write_DOS0:
      io.write_tot_dosfile()
  io.delete_tot_dosfile()

  #read in partial dos
  io.read_pdos()
  plt.pdos()
  if cont.write_PDOS:
      io.write_pdosfile()
  io.delete_pdosfile()

  #debug
  #dumpclean(cont.__dict__)
  #dumpclean(io.__dict__)
  #dumpclean(atoms.__dict__)
  #dumpclean(dos.__dict__)

if __name__ == '__main__':
  main()
