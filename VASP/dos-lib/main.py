from control import control_knob
from iodos import iodos
from data import atom_data,dos_data
from plot import plot

import sys

def dumpclean(obj):
  '''this source code comes from http://stackoverflow.com/questions/15785719/how-to-print-a-dictionary-line-by-line-in-python'''
  if type(obj) == dict:
    for k, v in obj.items():
      if type(v) == list:  #I changed it
          print k," : ",v  #I changed it
      elif hasattr(v, '__iter__'):
         print k
         dumpclean(v)
      else:
         if v:
           print '%s : %s' % (k, v)
  else:
    print obj


def main():
  atoms=atom_data()
  dos=dos_data()

  cont=control_knob(sys.argv,atoms,dos)
  plt=plot(cont)

  io=iodos(cont,plt)
  io.read_poscar()

  io.read_tot_dosfile()
  if cont.write_DOS0:
      io.write_tot_dosfile()
  io.delete_tot_dosfile()

  io.read_pdos()
  plt.pdos()

  dumpclean(cont.__dict__)
  dumpclean(io.__dict__)
  dumpclean(atoms.__dict__)
  dumpclean(dos.__dict__)
#
if __name__ == '__main__':
  main()
