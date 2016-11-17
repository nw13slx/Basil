#!/usr/bin/env python

#usage: split_dos_atom argument-name argument
#argument name, parameter
# partial orbital: "s, p, d, f", # of element; i.e. p 1 f 0
# input element species: atom.ntype, [# of types, type 0 #, type 1 #, ...]
#  this setting overrides the one in poscar/contcar
# write_DOS0 yes/no
# write_PDOS yes/no
# centerEf yes/no
# zoomin emin emax

class control_knob:

  par_symbol=['s','p','d','f']

  def __init__(self,argv,atom,dos):
    if not argv:
        self.end="no input argument"
        return self.end
    self.name=""
    self.path=""
    self.run_pdos=True
    self.write_DOS0=True
    self.write_PDOS=True
    self.peratom=False
    self.perspecies=True
    self.wholerange=False
    self.centerEf=True
    self.zoomIn=False
    self.zoomEmax=None
    self.zoomEmin=None
    self.energyshift=None
    self.end=None
    self.doscar=""
    self.poscar=""
    self.atom=atom
    self.dos=dos

    narg=len(argv)
    i=1
    while (i<narg and (not self.end)):
      #atom.ntype: #of_types #of_type_1, #of_type_2 ...
      if (argv[i].lower() =="name"):
        i+=1
        if (i<narg):
          self.name=argv[i]
          i+=1
        else:
          self.end="there are not enough arguments for task name"
      elif (argv[i].lower() =="path"):
        i+=1
        if (i<narg):
          self.path=argv[i]
          i+=1
        else:
          self.end="there are not enough arguments for task path"
      elif (argv[i].lower() =="doscar"):
        i+=1
        if (i<narg):
          self.doscar=argv[i]
          i+=1
        else:
          self.end="there are not enough arguments for doscar name"
      elif (argv[i].lower() =="poscar"):
        i+=1
        if (i<narg):
          self.poscar=argv[i]
          i+=1
        else:
          self.end="there are not enough arguments for task name"
      elif (argv[i].lower() =="atom.ntype"):
        i+=1
        if (i<narg):
          self.atom.ntype=int(argv[i+1])
          i+=1
          if narg>(atom.ntype+i):
            self.atom.atomSpe = []
            for j in range(atom.ntype):
              self.atom.atomSpe += ([j]*int(argv[j+i]))
            i+=self.atom.ntype
            self.atom.natom=len(self.atom.atomSpe)
          else:
            self.end="there are not enough arguments for atom.ntype"
        else:
            self.end="there are not enough arguments for atom.ntype"
      elif (argv[i].lower()=="energyshift"):
        i+=1
        if (i<narg):
          self.energyshift=float(argv[i])
          i+=1
        else:
          self.end="there are not enough arguments for energy shift"
      elif (argv[i].lower()=="peratom"):
        i+=1
        if (i<narg):
          if ((argv[i].lower()=="yes") or (argv[i].lower()=="y")):
            self.peratom=True
            i+=1
          elif ((argv[i].lower()=="no") or (argv[i].lower()=="n")):
            self.peratom=False
            i+=1
          else:
            self.end="the argument for peratom has to be yes/no or y/n"
        else:
          self.end="there are not enough arguments for peratom"
      elif (argv[i].lower()=="zoomin"):
        self.zoomIn=True
        i+=1
        if ((i+1)<narg):
          self.zoomEmin=float(argv[i])
          i+=1
          self.zoomEmax=float(argv[i])
          i+=1
        else:
          self.end="there are not enough arguments for zoomin"
      elif (argv[i].lower()=="write_dos0"):
        i+=1
        if (i<narg):
          if ((argv[i].lower()=="yes") or (argv[i].lower()=="y")):
            self.write_DOS0=True
            i+=1
          elif ((argv[i].lower()=="no") or (argv[i].lower()=="n")):
            self.write_DOS0=False
            i+=1
        else:
          self.end="there are not enough arguments for write_dos0"
      elif (argv[i].lower()=="centeref"):
        i+=1
        if (i<narg):
          if ((argv[i].lower()=="yes") or (argv[i].lower()=="y")):
            self.centerEf=True
            i+=1
          elif ((argv[i].lower()=="no") or (argv[i].lower()=="n")):
            self.centerEf=False
            i+=1
        else:
          self.end="there are not enough arguments for centerEf"
      elif (argv[i].lower()=="wholerange"):
        i+=1
        if (i<narg):
          if ((argv[i].lower()=="yes") or (argv[i].lower()=="y")):
            self.wholerange=True
            i+=1
          elif ((argv[i].lower()=="no") or (argv[i].lower()=="n")):
            self.wholerange=False
            i+=1
        else:
          self.end="there are not enough arguments for wholerange"
      elif (argv[i].lower()=="write_pdos"):
        i+=1
        if (i<narg):
          if ((argv[i].lower()=="yes") or (argv[i].lower()=="y")):
            self.write_PDOS=True
            i+=1
          elif ((argv[i].lower()=="no") or (argv[i].lower()=="n")):
            self.write_PDOS=False
            i+=1
        else:
          self.end="there are not enough arguments for write_PDOS"
      elif (argv[i].lower() in self.par_symbol):
        j=self.par_symbol.index(argv[i].lower())
        i+=1
        if (i<narg):
          self.dos.par_element[j]=int(argv[i])
          i+=1
        else:
          self.end="there are not enough arguments for partial orbital"
      else:
        self.end=argv[i]+" is not an argument name"
        i+=1
      if  (all(x==-1 for x in self.dos.par_element) and (self.peratom==False)):
        self.run_pdos=False
