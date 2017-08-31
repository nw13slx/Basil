#!/usr/bin/env python

#usage: split_dos_atom argument-name argument
#argument name, parameter
# partial orbital: "s, p, d, f", # of element; i.e. p 1 f 0
# input element species: atom.ntype, [# of types, type 0 #, type 1 #, ...]
#  this setting overrides the one in poscar/contcar
# write_dos0 yes/no
# write_pdos yes/no
# center_ef yes/no
# zoomin emin emax

class control_knob:

  par_symbol=['s','p','d','f']

  def __init__(self,argv,atom,dos):
    if not argv:
        self.end="no input argument"
        return self.end

    self.name=""
    self.path=""
    self.doscar=""
    self.poscar=""

    self.run_pdos=True
    self.write_dos0=True
    self.write_pdos=True
    self.peratom=False
    self.perspecies=True
    self.whole_range=False
    self.center_ef=True

    self.zoom_in=False
    self.zoom_emax=None
    self.zoom_emin=None
    self.energyshift=None
    self.end=None
    self.atom=atom
    self.dos=dos

    yes_and_no_list=[]
    string_value_list=['name','path','doscar','poscar']
    for k, v in self.__dict__.items():
      if type(v) == bool:  #I changed it
          yes_and_no_list +=[k]

    narg=len(argv)
    i=1
    while (i<narg and (not self.end)):
      #atom.ntype: #of_types #of_type_1, #of_type_2 ...
      if (argv[i].lower() in string_value_list):
        i+=1
        if (i<narg):
          self.__dict__[argv[i-1].lower()]=argv[i]
          i+=1
        else:
          self.end="there are not enough arguments for %s"%argv[i-1]
      elif (argv[i].lower() =="atom.ntype"):
        i+=1
        if (i<narg):
          self.atom.ntype=int(argv[i+1])
          i+=1
          if narg>(atom.ntype+i):
            self.atom.species = []
            for j in range(atom.ntype):
              self.atom.species += ([j]*int(argv[j+i]))
            i+=self.atom.ntype
            self.atom.natom=len(self.atom.species)
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
      elif (argv[i].lower() in yes_and_no_list):
        argv_name=argv[i].lower()
        i+=1
        if (i<narg):
          if ((argv[i].lower()=="yes") or (argv[i].lower()=="y")):
            self.__dict__[argv_name]=True
            i+=1
          elif ((argv[i].lower()=="no") or (argv[i].lower()=="n")):
            self.__dict__[argv_name]=False
            i+=1
          else:
            self.end="the argument for %s has to be yes/no or y/n"%argv_name
        else:
          self.end="there are not enough arguments for %s"%argv_name
      elif (argv[i].lower()=="zoomin"):
        self.zoom_in=True
        i+=1
        if ((i+1)<narg):
          self.zoom_emin=float(argv[i])
          i+=1
          self.zoom_emax=float(argv[i])
          i+=1
        else:
          self.end="there are not enough arguments for zoomin"
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
      if  ((self.perspecies==False) and all(x==-1 for x in self.dos.par_element) and (self.peratom==False)):
        self.run_pdos=False
