from control import control_knob
from iodos import iodos

import sys

def main():

  cont=control_knob(sys.argv)

  io=iodos(cont)
  io.read_poscar()

  io.read_tot_dosfile()
  if cont.write_DOS0:
      io.write_tot_dosfile()
  io.delete_tot_dosfile()
  io.read_pdos()

  cont.dumpclean()
  io.dumpclean()

if __name__ == '__main__':
  main()
