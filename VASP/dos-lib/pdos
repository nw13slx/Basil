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

#default setting
PATH="../"
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
loc_down=None
loc_up=None
ferminN=None

#main function
def main():

  f = open(PATH+"DOSCAR", 'r')
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

def per_atom_plot():
  cm = pl.get_cmap('winter')
  rad=12
  scale=1/23.*2.
  zmin=np.amin(position[:,2])
  zmax=np.amax(position[:,2])
  lz=zmax-zmin
  zcenter=(zmin+zmax)*0.5
  Oplot=plt.figure(200)
  Ceplot=plt.figure(300)
  y=Current[loc_down:loc_up,0]+Current[loc_down:loc_up,1]
  color =cm(2.*abs(position[atomi,2]-zcenter)/lz)  # color will now be an RGBA tuple
  y_off=(position[atomi,2]-zmin)/lz
  if (atomSpe[atomi]==par_element[1]):
    plt.figure(200)
    ymax=np.amax(y)
    pl.plot(ef[loc_down:loc_up],y_off+y/ymax*scale,c=color)
  if (atomSpe[atomi]==par_element[3]):
    plt.figure(300)
    ymax=np.amax(y)
    pl.plot(ef[loc_down:loc_up],y_off+y/ymax*scale,c=color)

def simple_line_plot(name,x,y,l_d,l_u,line_color,line_label):
  a=plt.figure(0)
  for i in range(len(y)):
    if (y[i]!=None):
      if ((line_color[i]!=None) and line_label[i]!=None):
        plt.plot(x[l_d:l_u],y[i][l_d:l_u],c=line_color[i],label=line_label[i])
      elif (line_color[i]==None):
        plt.plot(x[l_d:l_u],y[i][l_d:l_u],label=line_label[i])
      else:
        plt.plot(x[l_d:l_u],y[i][l_d:l_u])
  plt.legend()
  pl.savefig(name+".png")
  plt.close()

# READ DOSCAR

  matrix = np.hstack([ef.reshape([len(e),1]),data[:,1:]])
  if write_DOS0:
    np.savetxt('DOS-all-dat.gz',matrix,fmt='%15.8f')

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

def write_spin(f, positions,atomSpe, nedos, natoms, ncols, efermi,ef):
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
    print "atom",atomi
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
      per_atom_plot()
#    rad=12
#    scale=1/23.*2.
#    Oplot=plt.figure(200)
#    Ceplot=plt.figure(300)
#    if ((peratom==True)):
#      if (atomSpe[atomi]==par_element[1]):
#        plt.figure(200)
#        color =cm(2.*abs(position[atomi,2]-zcenter)/lz)  # color will now be an RGBA tuple
#        y=Current[:,0]+Current[:,1]
#        ymax=np.amax(y)
#        y_off=(position[atomi,2]-zmin)/lz
#        pl.plot(ef,y_off+y/ymax*scale,c=color)
#      if (atomSpe[atomi]==par_element[3]):
#        if (zoomIn==True):
#        plt.figure(300)
#        color =cm(2.*abs(position[atomi,2]-zcenter)/lz)  # color will now be an RGBA tuple
#        y=Current[:,0]+Current[:,1]
#        ymax=np.amax(y)
#        print ymax
#        y_off=(position[atomi,2]-zmin)/lz
#        pl.plot(ef,y_off+y/ymax*scale,c=color)
#        #ax_Ce.plot(ef,y,Current[:,0],c=color)
#        #ax_Ce.plot(ef,y,-Current[:,1],c=color)
#      #matrix = np.hstack([ef.reshape([len(ef),1]),element[:,1:]])
#      #np.savetxt('DOS'+str(atomi)+".dat",matrix,fmt='%15.8f')
#      #plt.figure(i)
#      #plt.plot(ef,Current[:,0],'r-',ef,Current[:,1],'b-')
#      #if (centerEf==True):
#      #  plt.xlim(-7,7)
#      #if (zoomIn==True):
#      #  plt.xlim(zoomEmin,zoomEmax)
#      #  #plt.ylim(-np.amax(data[loc_down:loc_up,2]), np.amax(data[loc_down:loc_up,1]))
#      #pl.show()
#      #pl.savefig("DOS"+str(atomi)+".png")
#      #plt.close()
#
#  plt.figure(200)
#  print zoomIn, zoomEmin,zoomEmax
#  if (zoomIn==True):
#    plt.xlim(zoomEmin,zoomEmax)
#  plt.ylim(0,1+scale)
#  plt.show()
#  plt.savefig("DOS_O_INSIDE.png")
#  plt.close()
#  plt.figure(300)
#  if (zoomIn==True):
#    plt.xlim(zoomEmin,zoomEmax)
#    plt.ylim(-Ce_max[1],Ce_max[0])
#  plt.ylim(0,1+scale)
#  plt.savefig("DOS_Ce_INSIDE.png")
#  plt.close()

  if  (np.sum(par_element)==-4):
    return par_orbital,d_t2g,d_eg,All

  par_symbol=['s','p','d','f']
  if (symbol!=None):
    for i in range(4):
      par_symbol[i]=symbol[par_element[i]]+"_"+par_symbol[i]
  else:
    for i in range(4):
      par_symbol[i]="atomtype="+str(par_element[i])+" "+par_symbol[i]+"_orbital"

  y=[]
  line_label=[]
  line_color=['y','y','r','r','b','b','k','k']
  for i in range(4):
    y+=[par_orbital[i][:,0],par_orbital[i][:,1]]
    line_label+=[par_symbol[i],None]
  simple_line_plot("pdf_orbital.png",x,y,loc_down,loc_up,line_color,line_label):
  del y, line_label, line_color

  if (par_element[2] >=0):
    y=[]
    line_label=['_t2g',None,'_eg',None]
    line_color=['r','r','b','b']
    y=[d_t2g[:,0],d_t2g[:,1],d_eg[:,0],d_eg[:,1]]
    simple_line_plot("d_decompose.png",x,y,loc_down,loc_up,line_color,line_label):
    del y, line_label, line_color

  y=[]
  line_label=[]
  for i in range(ntype):
    y+=[All[:,i]]
    line_label+=["species "+str(i)]
  if (symbol!=None):
    label+=[symbol[i]]
  simple_line_plot("All.png",x,y,loc_down,loc_up,None,line_label):
  del y, line_label, line_color

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
