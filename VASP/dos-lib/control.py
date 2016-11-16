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

class control_knob:

  _par_symbol=['s','p','d','f']

  def __init__(self,argv):
    if not argv:
        self.end="no input argument"
        return self.end
    self.name=""
    self.path=""
    self.write_DOS0=True
    self.write_PDOS=True
    self.peratom=False
    self.centerEf=True
    self.zoomIn=False
    self.zoomEmax=None
    self.zoomEmin=None
    self.energyshift=None
    self.par_element=[None,None,None,None] 
    self.ntype=None
    self.atomSpe=None
    self.end=None
    self.doscar=""
    self.poscar=""

    narg=len(argv)
    i=1
    while (i<narg and (not self.end)):
      #ntype: #of_types #of_type_1, #of_type_2 ...
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
      elif (argv[i].lower() =="ntype"):
        i+=1
        if (i<narg):
          self.ntype=int(argv[i+1])
          i+=1
          if narg>(ntype+i):
            self.atomSpe = []
            for j in range(ntype):
              self.atomSpe += ([j]*int(argv[j+i]))
            i+=self.ntype
          else:
            self.end="there are not enough arguments for ntype"
        else:
            self.end="there are not enough arguments for ntype"
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
          self.par_element[j]=int(argv[i])
          i+=1
          if (self.par_element[j]>=ntype):
            self.end="END, the partial orbital "+str(par_element[j])+"is higher than ntype "+str(ntype)
        else:
          self.end="there are not enough arguments for partial orbital"
      else:
        self.end=argv[i]+" is not an argument name"
        i+=1

  def dumpclean(self,obj0=None):
    '''this source code comes from http://stackoverflow.com/questions/15785719/how-to-print-a-dictionary-line-by-line-in-python'''
    if not obj0:
      obj=self.__dict__
    else:
      obj=obj0
    if type(obj) == dict:
      for k, v in obj.items():
        if type(v) == list:  #I changed it
            print k," : ",v  #I changed it
        elif hasattr(v, '__iter__'):
           print k
           self.dumpclean(v)
        else:
           if v:
             print '%s : %s' % (k, v)
    #elif type(obj) == list:
      #for v in obj:
      #  if hasattr(v, '__iter__'):
      #    self.dumpclean(v)
      #  else:
      #    print v
    else:
      print obj