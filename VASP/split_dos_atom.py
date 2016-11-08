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

  end=parseInput(sys.argv)
  if end:
    print "end condition: ",end
    return
  print zoomIn,zoomEmin,zoomEmax,zoomEmax-zoomEmin

  f = open("DOSCAR", 'r')
  #read the header of DOSCAR
  natoms, nedos, efermi = read_dosfile(f)

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

def parseInput(argv):
  if ( os.path.exists("DOSCAR") !=True):
    return "DOSCAR does not exist"
  global ntype, atomSpe, position,par_element,centerEf,peratom,write_PDOS,write_DOS0,zoomIn,zoomEmax,zoomEmin
  #if there is poscar
  pos = None
  if ( os.path.exists("POSCAR") == True ):
    pos = "POSCAR"
  if ( os.path.exists("CONTCAR") == True ):
    pos = "CONTCAR"
  if pos:
    position= read_POSCAR(pos)
    atomSpe,ntype=read_atom_Spef(pos)

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
    elif (argv[i].lower()=="energyshift"):
      i+=1
      energyshift=float(argv[i])
      i+=1
    elif (argv[i].lower()=="peratom"):
      i+=1
      if (argv[i].lower()=="yes"):
        peratom=True
      i+=1
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
  return None


# READ DOSCAR
def read_dosfile(f):
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

def write_dos0(f,nedos, efermi):
  '''
  WRITE DOS0 CONTAINING TOTAL DOS
  '''
  #attempt to read the first line and analysis the number of column
  start = f.tell()
  line = f.readline().strip().split()
  f.seek(start)
  ncols = int(len(line))

  ef = np.zeros(nedos)
  Current=np.zeros((nedos,2))

  #write the DOS0 file
  chunck = []
  for n in xrange(nedos):
      chunck.append(f.readline())
  data = np.loadtxt(chunck)
  e = data[:,0]
  if (centerEf==True):
    ef = (e-efermi)
  else:
    ef=e
  Current[:,:] = data[:,1:2]
  matrix = np.hstack([ef.reshape([len(e),1]),data[:,1:]])
  if write_DOS0:
    np.savetxt('DOS-all.dat',matrix,fmt='%15.8f')

  #plot the DOS0 file
  plt.figure(0)
  plt.plot(ef,data[:,1],'r-',ef,-data[:,2],'b-')
  if (centerEf==True):
    plt.xlim(-7,7)
    loc_down = np.argmin(abs(ef+7))
    loc_up = np.argmin(abs(ef-7))
    plt.ylim(-np.amax(data[loc_down:loc_up,2]), np.amax(data[loc_down:loc_up,1]))
    fermiN = np.argmin(abs(ef))
  else:
    fermiN = np.argmin(abs(ef-efermi))
  if (zoomIn==True):
    plt.xlim(zoomEmin,zoomEmax)
    loc_down = np.argmin(abs(ef-zoomEmin))
    loc_up = np.argmin(abs(ef-zoomEmax))
    plt.ylim(-np.amax(data[loc_down:loc_up,2]), np.amax(data[loc_down:loc_up,1]))
  if (energyshift !=0):
    delta = efermi-energyshift
  else:
    delta = 0
  plt.fill_between((ef+delta)[:fermiN], 0, data[:fermiN,1],facecolor='red')
  plt.fill_between((ef+delta)[:fermiN], 0, -data[:fermiN,2],facecolor='blue')
  pl.show()
  pl.savefig("DOS0.png")
  plt.close()

  #try to calculate the band gap
  EGMIN_u=30
  EGMIN_d=30
  EGMAX_u=-30
  EGMAX_d=-30
  foundEg=0
  finfo = open("band-info", 'w')
  finfo.write('between 0 %15.8f %15.8f \n' % (ef[fermiN],ef[fermiN+1]))
  for n in range(fermiN, 0, -1):
     if ((Current[n,0]!=0) and (EGMIN_u >0)):
        EGMIN_u=ef[n]
     if ((Current[n,1]!=0) and (EGMIN_d >0)):
        EGMIN_d=ef[n]
  for n in range(fermiN, nedos, 1):
     if ((Current[n,0]!=0) and (EGMAX_u <0)):
        EGMAX_u=ef[n]
     if ((Current[n,1]!=0) and (EGMAX_d <0)):
        EGMAX_d=ef[n]
  if (EGMIN_u >EGMIN_d):
      EGMIN = EGMIN_u
  else:
      EGMIN = EGMIN_d
  if (EGMAX_u < EGMAX_d):
      EGMAX = EGMAX_u
  else:
      EGMAX = EGMAX_d
  if ((Current[fermiN,0]!=0) and (Current[fermiN+1,0]!=0)):
      Eg_u=0
  else:
      Eg_u = EGMAX_u-EGMIN_u
  if ((Current[fermiN,1]!=0) and (Current[fermiN+1,1]!=0)):
      Eg_d=0
  else:
      Eg_d = EGMAX_d-EGMIN_d
  finfo.write('Eg_u %15.8f %15.8f %15.8f\n' % (EGMIN_u,EGMAX_u,Eg_u))
  finfo.write('Eg_d %15.8f %15.8f %15.8f\n' % (EGMIN_d,EGMAX_d,Eg_d))
  finfo.write('Eg %15.8f %15.8f %15.8f\n' % (EGMIN,EGMAX,EGMAX-EGMIN))
  finfo.close()
  return ef

def write_spin(f, positions, atomSpe, nedos, natoms, ncols, efermi,ef):
  print "hello"
  if  ((np.sum(par_element)==-4) and (peratom==False)):
    return None,None,None,None

  if (centerEf==True):
    loc_down = np.argmin(abs(ef+7))
    loc_up = np.argmin(abs(ef-7))
  if (zoomIn==True):
    loc_down = np.argmin(abs(ef-zoomEmin))
    loc_up = np.argmin(abs(ef-zoomEmax))
  nsites = (ncols -1)/2

  if (par_element[2]>=0):
    d_t2g=np.zeros((nedos,2))
    d_eg=np.zeros((nedos,2))
  else:
    d_t2g=None
    d_eg=None
  par_orbital=[] #store the partial s, p, d, f orbital
  for i in range(4):
    if par_element[i]>=0:
      par_orbital+=[np.zeros((nedos,2))]
    else:
      par_orbital+=[None]
  All=np.zeros((nedos,ntype))

  for atomi in xrange(natoms):
    #skip the first line, and loop over nedos
    line = f.readline()
    chunck=[]
    for n in xrange(nedos):
      chunck.append(f.readline())
    element = np.loadtxt(chunck)
    #ef = element[:,0] - efermi

    Current=np.zeros((nedos,2))
    for site in range(nsites):
      Current[:,0] += element[:,site*2+1]
      Current[:,1] += element[:,site*2+2]
      #collection
    if (atomSpe[atomi] == par_element[0]):
      par_orbital[0][:,0]+=element[:,1]
      par_orbital[0][:,1]-=element[:,2]
    if (atomSpe[atomi] == par_element[1]):
      par_orbital[1][:,0]+=element[:,3] + element[:,5] + element[:,7]
      par_orbital[1][:,1]-=(element[:,4] + element[:,6] + element[:,8])
    if (atomSpe[atomi] == par_element[2]):
      par_orbital[2][:,0]+=element[:,9] + element[:,11] + element[:,13] + element[:,15] + element[:,17]
      par_orbital[2][:,1]-=(element[:,10] + element[:,12] + element[:,14] + element[:,16] + element[:,18])
      d_t2g[:,0]+=element[:,9] + element[:,11] + element[:,15]
      d_t2g[:,1]-=(element[:,10] + element[:,12] + element[:,16])
      d_eg[:,0]+=element[:,13] + element[:,17]
      d_eg[:,1]-=(element[:,14] + element[:,18])
    if (atomSpe[atomi] == par_element[3]):
      par_orbital[3][:,0]+=element[:,19] + element[:,21] + element[:,23] + element[:,25] + element[:,27]+element[:,29]+element[:,31]
      par_orbital[3][:,1]-=(element[:,20] + element[:,22] + element[:,24] + element[:,26] + element[:,28]+element[:,30]+element[:,32])
    All[:,atomSpe[atomi]] += Current[:,0] + Current[:,1]

    #print data and png per atom
    if (peratom==True):
      matrix = np.hstack([ef.reshape([len(ef),1]),element[:,1:]])
      np.savetxt('DOS'+str(atomi)+".dat",matrix,fmt='%15.8f')
      plt.figure(i)
      plt.plot(ef,Current[:,0],'r-',ef,Current[:,1],'b-')
      if (centerEf==True):
        plt.xlim(-7,7)
      if (zoomIn==True):
        plt.xlim(zoomEmin,zoomEmax)
        #plt.ylim(-np.amax(data[loc_down:loc_up,2]), np.amax(data[loc_down:loc_up,1]))
      pl.show()
      pl.savefig("DOS"+str(atomi)+".png")
      plt.close()

  if  (np.sum(par_element)==-4):
    return par_orbital,d_t2g,d_eg,All

  par_symbol=['s','p','d','f']
  if (symbol!=None):
    for i in range(4):
      par_symbol[i]=symbol[par_element[i]]+"_"+par_symbol[i]
  else:
    for i in range(4):
      par_symbol[i]="atomtype="+str(par_element[i])+" "+par_symbol[i]+"_orbital"
  par_color=['y','r','b','k']
  par_lim=np.zeros(8)

  plt.figure(0)
  for i in range(4):
    if (par_element[i] >=0):
      plt.plot(ef,par_orbital[i][:,0],par_color[i]+'-',label=par_symbol[i])
      plt.plot(ef,par_orbital[i][:,1],par_color[i]+'-')
  plt.legend()
  if (centerEf==True):
    plt.xlim(-7,7)
  if (zoomIn ==True):
    plt.xlim(zoomEmin,zoomEmax)
  pl.show()
  pl.savefig("pdf_orbital.png")
  plt.close()

  if (par_element[2] >=0):
    plt.figure(0)
    plt.plot(ef,d_t2g[:,0],'r-',label='_t2g')
    plt.plot(ef,d_t2g[:,1],'r-')
    plt.plot(ef,d_eg[:,0],'b-',label='_eg')
    plt.plot(ef,d_eg[:,1],'b-')
    data = np.hstack([d_t2g[:,0],d_eg[:,0],d_t2g[:,1],d_eg[:,1]]).reshape([4,nedos]).T
    if (centerEf==True):
      #plt.ylim(np.amin(data[loc_down:loc_up,2:]), np.amax(data[loc_down:loc_up,0:2]))
      plt.xlim(-7,7)
    if (zoomIn ==True):
      plt.xlim(zoomEmin,zoomEmax)
    plt.legend()
    pl.show()
    pl.savefig("d_decompose.png")
    plt.close()

  plt.figure(0)
  for i in range(ntype):
    if (symbol!=None):
      plt.plot(ef,All[:,i],label=symbol[i])
    else:
      plt.plot(ef,All[:,i],label="species "+str(i))
  #plt.ylim(np.amin(All), np.amax(All))
  if (centerEf==True):
    plt.xlim(-7,7)
  if (zoomIn==True):
    plt.xlim(zoomEmin,zoomEmax)
  plt.legend()
  pl.show()
  pl.savefig("All.png")
  plt.close()
  return par_orbital,d_t2g,d_eg,All

def write_pdos(nedos,ef,par_orbital,name):
  par_symbol=['s','p','d','f']
  if (symbol!=None):
    for i in range(4):
      par_symbol[i]=symbol[par_element[i]]+"_"+par_symbol[i]
  else:
    for i in range(4):
      par_symbol[i]=str(par_element[i])+"_"+par_symbol[i]

  for i in range(4):
    if (par_element[i]>=0):
      fdos = open(name+par_symbol[i]+"_PDOS.dat", 'w')
      for n in xrange(nedos):
          fdos.write('%10.3e %10.3e %10.3e\n' %(ef[n],par_orbital[i][n,0],par_orbital[i][n,1]))
      fdos.close()

if __name__ == '__main__':
  main()
