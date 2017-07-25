#!/usr/bin/env python

#usage: split_dos_atom argument-name argument
#argument name, parameter
# "s, p, d, f": # of element that you want to look into for certain orbital; i.e. p 1 f 0
# ntype: [# of types, type 0 #, type 1 #, ...]
#  this setting overrides the one in poscar/contcar
# write_DOS0 yes/no
# write_PDOS yes/no
# centerEf yes/no
# zoomin emin emax
# atomDOS yes read from peratomDOS instead of DOSCAR

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
atomDOS=False

#initialization
ntype=0
natoms=0
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

  end=parseInput(sys.argv)
  if end:
    print "end condition: ",end
    return
  print atomDOS,zoomIn,zoomEmin,zoomEmax,zoomEmax-zoomEmin

  if (atomDOS == False):
    f = open("DOSCAR", 'r')
    natoms, nedos, efermi = read_dosfile_header(f)
    ef = write_dos0(f, nedos, efermi)
    start = f.tell()
    line = f.readline()
    line = f.readline().strip().split()
    ncols = int(len(line))
    if ncols==7 or ncols==19 or ncols==9 or ncols==33:
      f.seek(start)
      par_orbital, d_t2g, d_eg, All = write_spin(f,position,atomSpe,nedos, natoms, ncols, efermi,ef)
      is_spin=True
      if (write_PDOS==True):
        write_pdos(nedos,ef,par_orbital,name)
      f.close()
      del par_orbital,d_t2g,d_eg,All,ef,nedos
  else:
    f = open("DOSCAR", 'r')
    natoms, nedos, efermi = read_dosfile_header(f)
    f.close()
    plt.figure(0)
    cm = pl.get_cmap('winter')
    ef=None
    for atomi in range(natoms):
      if (atomSpe[atomi]==1):
        f=open("DOS"+str(atomi)+".dat")
        data=np.loadtxt(f)
        ef=data[:,0]
        color = cm(1.*atomi/float(natoms))
        plt.plot(data[:,0],data[:,1],c=color)
        plt.plot(data[:,0],-data[:,2],c=color)
        f.close()
        del data
    if (centerEf==True):
      plt.xlim(-7,7)
    if (zoomIn==True):
      plt.xlim(zoomEmin,zoomEmax)
    pl.show()
    pl.savefig("corestate.png")
    plt.close()

def parseInput(argv):
  if (os.path.exists("DOSCAR") !=True):
    return "DOSCAR does not exist"
  global ntype, atomSpe, position,par_element,natoms,zoomEmax,zoomEmin
  global centerEf,peratom,write_PDOS,write_DOS0,zoomIn
  global atomDOS
  #if there is poscar
  pos = None
  if ( os.path.exists("POSCAR") == True ):
    pos = "POSCAR"
  if ( os.path.exists("CONTCAR") == True ):
    pos = "CONTCAR"
  if pos:
    position= read_POSCAR(pos)
    atomSpe,ntype=read_atom_Spef(pos)
    natoms=len(atomSpe)

  i=0
  par_symbol=["s","p","d","f"]
  while i<len(argv):
    i0=i
    if (argv[i] =="ntype"):
      i+=1
      ntype=int(argv[i])
      i+=1
      if len(argv)>(ntype+i):
        atomSpe = []
        for j in range(ntype):
          atomSpe += ([j]*int(argv[j+i]))
        i+=ntype
        natoms=len(atomSpe)
    elif (argv[i].lower()=="energyshift"):
      i+=1
      energyshift=float(argv[i])
      i+=1
    elif (argv[i].lower()=="peratom"):
      i+=1
      if (argv[i].lower()=="yes"):
        peratom=True
      i+=1
    elif (argv[i].lower()=="atomdos"):
      i+=1
      if (argv[i].lower()=="yes"):
        atomDOS=True
    elif (argv[i].lower()=="zoomin"):
      zoomIn=True
      i+=1
      zoomEmin=float(argv[i])
      i+=1
      zoomEmax=float(argv[i])
      i+=1
    elif (argv[i].lower()=="write_dos0"):
      i+=1
      if (argv[i].lower()=="no"):
        write_DOS0=False
      i+=1
    elif (argv[i].lower()=="centeref"):
      i+=1
      if (argv[i].lower()=="no"):
        centerEf=False
      i+=1
    elif (argv[i].lower()=="write_pdos"):
      i+=1
      if (argv[i].lower()=="no"):
        write_PDOS=False
    else:
      j=0
      while(i<len(argv) and j <4):
        if argv[i].lower() == par_symbol[j]:
          i+=1
          par_element[j]=int(argv[i])
          i+=1
        j+=1
    if (i0==i):
      i+=1
  if (atomSpe==None):
    return "END, not definition of atom types"
  for j in range(4):
    if (par_element[j]>=ntype):
      return "END, the partial orbital "+str(par_element[j])+"is higher than ntype "+str(ntype)
  if (atomDOS==True):
    for atomi in range(natoms):
      if (os.path.exists("DOS"+str(atomi)+".dat") == False ):
        return "END, the DOS for atom "+str(atomi)+"does not exist"
  return None

def read_dosfile_header(f):
  line = f.readline()
  natoms = int(line.strip().split()[0])
  for i in range(5):
    line = f.readline()
  split=line.strip().split()
  nedos = int(split[2])
  efermi = float(split[3])
  return natoms, nedos, efermi

def read_atom_Spef(pos):
  '''
  from POSCAR or CONTCAR, read the types of atoms
  '''
  a=open(pos)
  if (pos == None):
    return None
  lines=a.readlines()
  #for older version, this line number may not be valid
  global symbol
  symbol=lines[5].strip().split()
  number=map(int,lines[6].strip().split())
  ntype=len(number)
  del a, lines
  atomSpe=[]
  for i in range(len(number)):
    atomSpe += ([i]*number[i])
  return atomSpe,ntype

def read_POSCAR(pos):
  '''
  this function is to read in the atomic coordinations from VASP "POSCAR" or "CONTCAR"; if POSCAR and CONTCAR do not exist, the number of atoms will be read in from command line argument; "python split_dos.py #_of_elements #_of_La #Co #O #others..."
  '''
  print "read from",pos
  import aselite
  atoms = aselite.read_vasp(pos)
  position = atoms.get_positions()
  return position

def
if __name__ == '__main__':
  main()
