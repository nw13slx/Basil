from doslib.control import *
from doslib.iodos import *
from doslib.data import *
from doslib.plot import *
from doslib.analysis import *

import sys

def main():
  #initialization
  args=Namespace(peratom=True, perspecies=False, poscar='', s=-1, write_dos0=False, write_pdos=False, zoomin=None)
  atoms=atom_data()
  dos=dos_data()
  cont=control_knob(args,atoms,dos)
  plt=plot(cont)
  ana=analysis(atoms,dos)

  ##read in overall dos
  io=iodos(cont,plt)
  io.read_poscar()

#  io.read_tot_dosfile()
#  ana.bandgap()
#  print ana.peak_weight_center(-32.5,5)
#
#  #write
#  if cont.write_dos0:
#      io.write_tot_dosfile()
#  io.delete_tot_dosfile()
#
#  #read in partial dos
#  if cont.run_pdos:
#    io.read_pdos()
#    plt.pdos()
#    if cont.write_pdos:
#        io.write_pdosfile()
#    io.delete_pdosfile()

class Namespace:
  def __init__(self, **kwargs):
    self.__dict__.update(kwargs)

if __name__ == '__main__':
  main()
