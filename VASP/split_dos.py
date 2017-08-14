#
#The original script comes from

#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl
import sys

energyshift=0
if (len(sys.argv)>(int(sys.argv[1])+2)):
    energyshift=float(sys.argv[int(sys.argv[1])+2])
LaN=int(sys.argv[2])
CoN=int(sys.argv[3])
ON=int(sys.argv[4])
LaMIN=0
LaMAX=LaN
CoMIN=LaMAX
CoMAX=LaMAX+CoN
OMIN=CoMAX
OMAX=CoMAX+ON
print "input",LaN, CoN, ON,energyshift

### READ DOSCAR ###
def read_dosfile():
    f = open("DOSCAR", 'r')
    lines = f.readlines()
    f.close()
    #first line, the first number
    index = 0
    natoms = int(lines[index].strip().split()[0])
    index = 5
    #sixth line, the third and fourth number
    nedos = int(lines[index].strip().split()[2])
    efermi = float(lines[index].strip().split()[3])
    if (energyshift!=0):
        print "shift energy by", energyshift, "instead of efermi", efermi
        print "valence band difference ", efermi-energyshift
    #print natoms, nedos, efermi
    return lines, index, natoms, nedos, efermi

### WRITE DOS0 CONTAINING TOTAL DOS ###
def write_dos0(lines, index, nedos, efermi):

    fdos = open("DOS0", 'w')
    index +=1
    line = lines[index+2].strip().split()
    ncols = int(len(line))
    fdos.write('# %d \n' % (ncols))
    ef = np.zeros(nedos)
    Current=np.zeros((nedos,2))
    fermiN=0

    for n in xrange(nedos):
        e = float(lines[index].strip().split()[0])
        e_f = e-efermi
        ef[n] = e_f
        fdos.write('%15.8f ' % (e_f))
        for col in xrange(1, ncols):
            dos = float(lines[index].strip().split()[col])
            fdos.write('%15.8e ' % (dos))
            if ((e_f >-7) and (e_f<7)):
               if (col == 1) :
                  Current[n,0] = dos
               if (col == 2) :
                  Current[n,1] = dos
            if (e_f ==0):
                 fermiN=n
            if (n >1):
               if ((e_f > 0) and (ef[n-1] < 0)):
                 fermiN=n-1
            col +=1
        fdos.write('\n ')
        index +=1

    EGMIN_u=30
    EGMIN_d=30
    EGMAX_u=-30
    EGMAX_d=-30
    foundEg=0
    finfo = open("band-info", 'w')
    finfo.write('between 0 %15.8f %15.8f \n' % (ef[fermiN],ef[fermiN+1]))
    print ef[fermiN],ef[fermiN+1]
    for n in range(fermiN, 0, -1):
       if ((Current[n,0]>1e-8) and (EGMIN_u >0)):
          EGMIN_u=ef[n]
          Emin_u=Current[n,0]
       else:
          if ((Current[n,0]<=1e-8) and (Current[n,0]!=0) and (EGMIN_u>0)):
             print "exception egmin_u", ef[n],Current[n,0]
       if ((Current[n,1]>1e-8) and (EGMIN_d >0)):
          EGMIN_d=ef[n]
          Emin_d=Current[n,0]
       else:
          if ((Current[n,1]<=1e-8) and (Current[n,1]!=0) and (EGMIN_d>0)):
             print "exception egmin_d", ef[n],Current[n,1]
    for n in range(fermiN, nedos, 1):
       if ((Current[n,0]>1e-8) and (EGMAX_u <0)):
          EGMAX_u=ef[n]
          Emax_u=Current[n,0]
       else:
          if ((Current[n,0]<1e-8) and (Current[n,0]!=0) and (EGMAX_u<0)):
             print "exception egmax_u", ef[n],Current[n,0]
       if ((Current[n,1]>1e-8) and (EGMAX_d <0)):
          EGMAX_d=ef[n]
          Emax_d=Current[n,1]
       else:
          if ((Current[n,1]<1e-8) and (Current[n,1]!=0) and (EGMAX_d<0)):
             print "exception egmax_d", ef[n],Current[n,1]
    if ((Current[fermiN,0]>1e-8) and (Current[fermiN+1,0]>1e-8)):
        Eg_u=0
    else:
        Eg_u = EGMAX_u-EGMIN_u
    if ((Current[fermiN,1]>1e-8) and (Current[fermiN+1,1]>1e-8)):
        Eg_d=0
    else:
        Eg_d = EGMAX_d-EGMIN_d
    if (Eg_d < Eg_u):
        Egg=Eg_d
    else:
        Egg=Eg_u

    finfo.write('Eg_u %15.8f %15.8e %15.8f %15.8e %15.8f\n' % (EGMIN_u,Emin_u,EGMAX_u,Emax_u,Eg_u))
    finfo.write('Eg_d %15.8f %15.8e %15.8f %15.8e %15.8f\n' % (EGMIN_d,Emin_d,EGMAX_d,Emax_d,Eg_d))
    finfo.write('Eg %15.8f \n' % (Egg))
    print EGMIN_u,Emin_u,EGMAX_u,Emax_u,Eg_u
    print EGMIN_d,Emin_d,EGMAX_d,Emax_d,Eg_d
    print Egg
    plt.figure(0)
    if (energyshift !=0):
        delta = efermi-energyshift
    else:
        delta = 0
    plt.plot(ef+delta,Current[:,0],'r-',ef+delta,-Current[:,1],'b-')
    plt.fill_between((ef+delta)[:fermiN], 0, Current[:fermiN,0],facecolor='red')
    plt.fill_between((ef+delta)[:fermiN], 0, -Current[:fermiN,1],facecolor='blue')
    plt.xlim(-7,7)
    pl.show()
    pl.savefig("DOS0.png")
    plt.close()
    fdos.close()
    finfo.close()
    return index,ef

def write_spin(lines, index, nedos, natoms, ncols, efermi,ef):
    #pos=[]

    nsites = (ncols -1)/2
    Op=np.zeros((nedos,2))
    Cot2g=np.zeros((nedos,2))
    Coeg=np.zeros((nedos,2))
    Cod=np.zeros((nedos,2))
    Current=np.zeros((nedos,2))

    for i in xrange(natoms):
        chargei=0
        magi=0
        ## OPEN DOSi FOR WRITING ##
        #fdos = open("DOS"+str(i), 'w')
        index +=1
        if ((i>=OMIN) and (i<OMAX)):
            Oxygen = True
            #print "O",i
        else:
            Oxygen = False
        if ((i>=CoMIN) and (i<CoMAX)):
            Co = True
            #print "Co",i
        else:
            Co = False

        #fdos.write('# %d \n' % (ncols))
        #fdos.write('# %10.3f %10.3f %10.3f \n' % (pos[i,0], pos[i,1], pos[i,2]))

        ### LOOP OVER NEDOS ###
        for n in xrange(nedos):
            string_e = lines[index].strip().split()
            element = [float(x) for x in string_e]
            e = element[0]
            #e_f = e-efermi
            #fdos.write('%f ' % (e_f))

            #col=1   #edit by leleslx
            #for i in xrange(nsites):   #edit by leleslx
            #    dos_up = element[col]
            #    dos_down = element[col+1]*-1
            #    fdos.write('%10.3e %10.3e' % (dos_up, dos_down))
            #    col +=2  #edit by leleslx
            #fdos.write('\n ')

            #collection
            if (ef[n]<0):
               chargei +=element[1]
               chargei +=element[2]
               magi += element[1]
               magi -= element[2]
            if ((ef[n] >-7) and (ef[n]<7)):
               if (Oxygen == True):
                  Op[n,0]+=element[3]+element[5]+element[7]
                  Op[n,1]-=(element[4]+element[6]+element[8])
               if (Co == True):
                  Cod[n,0]+=element[9]+element[11]+element[13]+element[15]+element[17]
                  Cod[n,1]-=(element[10]+element[12]+element[14]+element[16]+element[18])
                  Cot2g[n,0]+=(element[9]+element[11]+element[15])
                  Cot2g[n,1]-=(element[10]+element[12]+element[16])
                  Coeg[n,0]+=(element[13]+element[17])
                  Coeg[n,1]-=(element[14]+element[18])
               Current[n,0] = element[1]
               Current[n,1] = -element[2]
            index +=1
        print "atom",i, chargei, magi
        #plt.figure(i)
        #plt.plot(ef,Current[:,0],'r-',ef,Current[:,1],'b-')
        #plt.xlim(-7,7)
        #pl.show()
        #pl.savefig(str(i)+"DOS.png")
        #plt.close()
    if (energyshift !=0):
        delta = efermi-energyshift
    else:
        delta = 0
    plt.figure(0)
    plt.plot(ef+delta,Op[:,0],'r-',ef+delta,Op[:,1],'r-',ef+delta,Cod[:,0],'b-',ef+delta,Cod[:,1],'b-')
    plt.xlim(-7,7)
    pl.show()
    pl.savefig("Op.png")
    plt.close()
    plt.figure(0)
    plt.plot(ef+delta,Cot2g[:,0],'r-',ef+delta,Cot2g[:,1],'r-',ef+delta,Coeg[:,0],'b-',ef+delta,Coeg[:,1],'b-')
    plt.xlim(-7,7)
    pl.show()
    pl.savefig("Co.png")
    plt.close()
        #fdos.close()
        #fdos.close()
        #fdos.close()
    return Op, Cot2g, Coeg, Cod

def write_PDOS(nedos,ef,Op,name):
    fdos = open(name+"PDOS.dat", 'w')
    for n in xrange(nedos):
        fdos.write('%10.3e %10.3e %10.3e\n' %(ef[n],Op[n,0],Op[n,1]))
    fdos.close()

#
if __name__ == '__main__':
 import sys
 import os
 import datetime
 import time
 import optparse

    #read in file, all lines store in 'lines' variable
    #index is the last line processed
    #nedos, number of lines for each DOS
 lines, index, natoms, nedos, efermi = read_dosfile()

    #write out the total DOS to DOS0
 index,ef = write_dos0(lines, index, nedos, efermi)

 ##Test if there a spin calculation was performed ##
 line = lines[index+2].strip().split()
 ncols = int(len(line))
 if ncols==7 or ncols==19 or ncols==9 or ncols==33:
    Op, Cot2g, Coeg, Cod = write_spin(lines, index, nedos, natoms, ncols, efermi,ef)
    is_spin=True

    write_PDOS(nedos,ef,Op,"OP")
    write_PDOS(nedos,ef,Cot2g,"Cot2g")
    write_PDOS(nedos,ef,Coeg,"Coeg")
    write_PDOS(nedos,ef,Cod,"Cod")
 else:
  #write_nospin(lines, index, nedos, natoms, ncols, efermi)
  is_spin=False
 print "Spin unrestricted calculation: ", is_spin
    #if is_spin:
    # write_spin(lines, index, ncols, natoms, nedos, efermi)
    #else:
   #  write_nospin(lines, index, ncols, natoms, nedos, efermi)
