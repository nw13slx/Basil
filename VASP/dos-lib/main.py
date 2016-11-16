from control import control_knob
from io import iodos

import sys

def main():

  cont=control_knob(sys.argv)
  if cont.end:
    print "end condition: ",cont.end
    return
  cont.dumpclean()

  io=iodos(cont)
  if io.end:
    print "end condition: ",io.end
    return
  io.read_poscar()
  io.read_tot_dosfile()
  io.dumpclean()

if __name__ == '__main__':
  main()
