    if (write_PDOS==True):
      write_pdos(nedos,ef,par_orbital,name)
    f.close()
    del par_orbital,d_t2g,d_eg,All,ef,nedos

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
