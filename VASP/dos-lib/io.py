#data storage

#!/usr/bin/env python
import numpy as np
import os

class iodos:
  def __init__(self,control):
    ''' this function is to find out atomic coordinations and DOSCAR; the priority is command line input "ntype" > "poscar" input>CONTCAR >POSCAR , and "doscar" input >DOSCAR) '''
    self._cont=control
    self.end=None

    self._posf=None
    self._dosf=None

    self.natom=None
    self.ntype=None
    self.symbol=None
    self.atomSpe=None
    self.positions=None
    self.constraint=None

    self.nedos = None
    self.efermi = None
    self.Xenergy = None
    self.loc_down = None 
    self.loc_up = None
    self.fermiN = None
  def _find_file(self):
    control=self._cont
    if control.ntype:
      self.atomSpe=control.atomSpe
      self.natom=len(self.atomSpe)
    elif (os.path.isfile(control.path+control.poscar) and control.poscar):
      self._posf=open(control.path+control.poscar)
    elif os.path.isfile(control.path+"CONTCAR"):
      self._posf =open(control.path+ "CONTCAR")
    elif os.path.isfile(control.path+"POSCAR") :
      self._posf =open(control.path+"POSCAR")
    else:
      self.end="cannot find poscar"

    if os.path.isfile(control.path+control.doscar) and control.doscar:
      self._dosf=open(control.path+control.doscar)
    elif os.path.isfile(control.path+"DOSCAR"):
      self._dosf =open(control.path+"DOSCAR")
    else:
      self.end="cannot find doscar"
    return self.end

  def read_poscar(self):
    if self._find_file():
        return self.end
    if not self._posf:
        return self.end
    lines=self._posf.readlines()
    if not self.atomSpe:
      self.symbol=lines[5].strip().split()
      number=map(int,lines[6].strip().split())
      self.ntype=len(number)
      self.atomSpe=[]
      for i in range(len(number)):
        self.atomSpe += ([i]*number[i])
      self.natom=len(self.atomSpe)
    if ('S' in lines[7]) or ('s' in lines[7]):
        starting=9
        a=np.array([line.strip().split() for line in lines[starting:(starting+self.natom)]])
        self.positions=a[:,0:3]
        self.constraints=a[:,3:6]
    else:
        starting=8
        self.positions = np.loadtxt(lines[starting:(starting+self.natom)])
    self._posf.close()
  def read_tot_dosfile(self):
    control=self._cont
    line = self._dosf.readline()
    if (int(line.strip().split()[0]) != self.natom):
        self.end="doscar atom number is different from poscar"
        return 1
    for i in range(5):
      line = self._dosf.readline()
    split=line.strip().split()
    self.nedos = int(split[2])
    self.efermi = float(split[3])

    self.Xenergy = np.zeros(self.nedos)
    Current=np.zeros((self.nedos,2))
    #read the DOS0 file
    chunck = []
    for n in xrange(self.nedos):
        chunck.append(self._dosf.readline())
    data = np.loadtxt(chunck)
    self.Xenergy = data[:,0]
    Current[:,:] = data[:,1:2]

    if (control.centerEf==True):
      self.Xenergy -= self.efermi
      self.loc_down = np.argmin(abs(self.Xenergy+7))
      self.loc_up = np.argmin(abs(self.Xenergy-7))
      self.fermiN = np.argmin(abs(self.Xenergy))
    else:
      self.fermiN = np.argmin(abs(self.Xenergy-efermi))
    if (control.zoomIn==True):
      self.loc_down = np.argmin(abs(self.Xenergy-zoomEmin))
      self.loc_up = np.argmin(abs(self.Xenergy-zoomEmax))
  def dumpclean(self,obj0="everything"):
    '''this source code comes from http://stackoverflow.com/questions/15785719/how-to-print-a-dictionary-line-by-line-in-python'''
    if type(obj0)==str:
        if (obj0=="everything"):
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

#def read_pdos():
#  if ( dosf == None):
#      return "End, DOSCARfile is not opened yet"
#  global ncols 
#  start = dosf.tell()
#  line = dosf.readline()
#  line = dosf.readline().strip().split()
#  ncols= int(len(line))
#  dosf.seek(start)
#  if !(ncols==7 or ncols==19 or ncols==9 or ncols==33):
#      return "End, DOSCARfile do not have spin component"
#  return None
#def read_dos0():
#  #attempt to read the first line and analysis the number of column
#  ef = np.zeros(nedos)
#  Current=np.zeros((nedos,2))
#
#  #write the DOS0 file
#  chunck = []
#  for n in xrange(nedos):
#      chunck.append(f.readline())
#  global dos0 = np.loadtxt(chunck)
#
#def read_parorbit():
#  if  ((np.sum(par_element)==-4) and (peratom==False)):
#    return "END no partial orbital is needed"
#  nsites = (ncols -1)/2
#  for i in range(4):
#    if ((par_element[i]>=0) and (nsites<(2*i+1)):
#        return "END, nsites is smaller than the required orbital"
#
#  global d_t2g,d_eg, par_orbital,perelement
#  if (par_element[2]>=0):
#    d_t2g=np.zeros((nedos,2))
#    d_eg=np.zeros((nedos,2))
#  else:
#    d_t2g=None
#    d_eg=None
#  par_orbital=[] #store the partial s, p, d, f orbital
#  for i in range(4):
#    if par_element[i]>=0:
#      par_orbital+=[np.zeros((nedos,2))]
#    else:
#      par_orbital+=[None]
#  perelement=np.zeros((nedos,ntype))
#
#  for atomi in xrange(natoms):
#    #skip the first line, and loop over nedos
#    xEnergy,tot,partial=read_atomDOS(dosf)
#
#    if (atomSpe[atomi] == par_element[0]):
#      par_orbital[0][:,0]+=partial[:,0]
#      par_orbital[0][:,1]-=partial[:,1]
#    if (atomSpe[atomi] == par_element[1]):
#      par_orbital[1][:,0]+=partial[:,2]
#      par_orbital[1][:,1]-=partial[:,3]
#    if (atomSpe[atomi] == par_element[2]):
#      par_orbital[2][:,0]+=partial[:,4]
#      par_orbital[2][:,1]-=partial[:,5]
#      d_t2g[:,0]+=partial[:,8]
#      d_t2g[:,1]-=partial[:,9]
#      d_eg[:,0]+=partial[:,10]
#      d_eg[:,1]-=partial[:,11]
#    if (atomSpe[atomi] == par_element[3]):
#      par_orbital[3][:,0]+=partial[:,6]
#      par_orbital[3][:,1]-=partial[:,7]
#    All[:,atomSpe[atomi]] += tot[:,0] + tot[:,1]
#
#    #print data and png per atom
#    if (peratom==True):
#      matrix = np.hstack([xEnergy.reshape([len(xEnergy),1]),partial])
#      header_line="#"+name+"atom"+str(atomi)+"\n#E su sd pu pd du dd fu fd dt2gu dt2gd degu degd"
#      np.savetxt(name+'DOS'+str(atomi)+".dat",matrix,fmt='%15.8f',header=header_line)
#
#def read_atomDOS(pdosf)
#    line = pdosf.readline()
#    chunck=[]
#    for n in xrange(nedos):
#      chunck.append(pdosf.readline())
#    element = np.loadtxt(chunck)
#
#    xEnergy=np.zeros(nedos)
#    tot=np.zeros((nedos,2))
#    partial=np.zeros((nedos,12)) #s,p,d,f,dt2g,deg
#    for site in range(nsites):
#      tot[:,0] += element[:,site*2+1]
#      tot[:,1] += element[:,site*2+2]
#    xEnergy = element[:,0]
#    #collection, s, p, d, f, dt2g, deg
#    partial[:,0]=element[:,1] #s
#    partial[:,1]=element[:,2] #s
#    partial[:,2]=element[:,3] + element[:,5] + element[:,7] #p
#    partial[:,3]=element[:,4] + element[:,6] + element[:,8] #p
#    partial[:,4]=element[:,9] + element[:,11] + element[:,13] + element[:,15] + element[:,17]
#    partial[:,5]=element[:,10] + element[:,12] + element[:,14] + element[:,16] + element[:,18]
#    partial[:,6]=element[:,19] + element[:,21] + element[:,23] + element[:,25] + element[:,27]+element[:,29]+element[:,31]
#    partial[:,7]=(element[:,20] + element[:,22] + element[:,24] + element[:,26] + element[:,28]+element[:,30]+element[:,32])
#    partial[:,8]=element[:,9] + element[:,11] + element[:,15]
#    partial[:,9]=element[:,10] + element[:,12] + element[:,16]
#    partial[:,10]=element[:,13] + element[:,17]
#    partial[:,11]=element[:,14] + element[:,18]
#    return xEnergy,tot,partial
