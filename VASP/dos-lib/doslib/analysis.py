import numpy as np
class analysis:
  def __init__(self,atom,dos):
      self.atom=atom
      self.dos=dos

  def bandgap(self):
    #try to calculate the band gap
    ef=self.dos.Xenergy
    dos=self.dos
    dos0=dos.dos0
    fermiN=self.dos.fermiN
    nedos=self.dos.nedos

    EGMIN_u=None
    EGMIN_d=None
    EGMAX_u=None
    EGMAX_d=None
    foundEg=0
    finfo = open("band-info", 'w')
    finfo.write('#search between 0 %15.8f %15.8f \n' % (ef[fermiN],ef[fermiN+1]))
    for n in range(fermiN, 0, -1):
       if ((dos0[n,0]!=0) and (EGMIN_u is None)):
          EGMIN_u=ef[n]
       if ((dos0[n,1]!=0) and (EGMIN_d is None)):
          EGMIN_d=ef[n]
    for n in range(fermiN+2, nedos, 1):
       if ((dos0[n,0]!=0) and (EGMAX_u is None)):
          EGMAX_u=ef[n]
       if ((dos0[n,1]!=0) and (EGMAX_d is None)):
          EGMAX_d=ef[n]
    if (EGMIN_u >EGMIN_d):
        EGMIN = EGMIN_u
    else:
        EGMIN = EGMIN_d
    if (EGMAX_u < EGMAX_d):
        EGMAX = EGMAX_u
    else:
        EGMAX = EGMAX_d
    if ((dos0[fermiN,0]!=0) and (dos0[fermiN+1,0]!=0) and (dos0[fermiN+2,0]!=0 )):
        Eg_u=0
    else:
        Eg_u = EGMAX_u-EGMIN_u
    if ((dos0[fermiN,1]!=0) and (dos0[fermiN+1,1]!=0) and (dos0[fermiN+2,1]!=0 )):
        Eg_d=0
    else:
        Eg_d = EGMAX_d-EGMIN_d
    finfo.write('#spin     %15s %15s %15s\n' % ("VBM","CBM","gap"))
    finfo.write('spin_up   %15.8f %15.8f %15.8f\n' % (EGMIN_u,EGMAX_u,Eg_u))
    finfo.write('spin_down %15.8f %15.8f %15.8f\n' % (EGMIN_d,EGMAX_d,Eg_d))
    finfo.write('total     %15.8f %15.8f %15.8f\n' % (EGMIN,EGMAX,EGMAX-EGMIN))
    finfo.close()


  def peak_finder(self,center,span):
    ef=self.dos.Xenergy
    dos=self.dos
    dos0=dos.dos0
    fermiN=self.dos.fermiN
    nedos=self.dos.nedos

    window_d=center-span/2.
    window_u=center+span/2.
    win_down= np.argmin(abs(dos.Xenergy-window_d))
    win_up= np.argmin(abs(dos.Xenergy-window_u))
    win_span=(win_up-win_down)
    ana_dos=dos0[win_down:win_up,0]-dos0[win_down:win_up,1]

    from scipy.signal import find_peaks_cwt
    peak=find_peaks_cwt(ana_dos,np.arange(5,20),min_snr=0.8)
    sigma=(np.square(np.arange(win_span)-peak)).dot(ana_dos)/np.sum(ana_dos)
    sigma=np.sqrt(sigma)
    dx=dos.Xenergy[1]-x0
    #print avg*dx+x0,sigma*dx
    return dos.Xenergy[win_down+peak],sigma*dx

  def peak_weight_center(self,center,span):
    #try to calculate the band gap
    ef=self.dos.Xenergy
    dos=self.dos
    dos0=dos.dos0
    fermiN=self.dos.fermiN
    nedos=self.dos.nedos

    window_d=center-span/2.
    window_u=center+span/2.
    win_down= np.argmin(abs(dos.Xenergy-window_d))
    win_up= np.argmin(abs(dos.Xenergy-window_u))
    win_span=(win_up-win_down)
    ana_dos=dos0[win_down:win_up,0]-dos0[win_down:win_up,1]

    #assuming there is only one gaussian peak
    avg=(np.arange(win_span).dot(ana_dos))/np.sum(ana_dos)
    sigma=(np.square(np.arange(win_span)-avg)).dot(ana_dos)/np.sum(ana_dos)
    sigma=np.sqrt(sigma)

    x0=dos.Xenergy[0]
    dx=dos.Xenergy[1]-x0
    return dos.Xenergy[int(win_down+avg)],sigma*dx

  def find_atomic_plane(self,positions,axis):
    #bin all positions towrad one axis
    if (axis=='x'):
      direction=0
    elif (axis=='y'):
      direction=1
    elif (axis=='z'):
      direction=2
    else:
      print 'please input x, y, or z'
      return None,None,None

    natom=len(positions)
    #get the seperation part for each plane
    dz=0.1
    zmin=np.min(positions[:,direction])-3
    zmax=np.max(positions[:,direction])+3
    binz=np.zeros(int(np.ceil((zmax-zmin)/dz)))
    for i in range(natom):
      z=positions[i,direction]
      idz=int(np.floor((z-zmin)/dz))
      binz[idz]+=1
    from scipy.signal import find_peaks_cwt
    peak_atoms=find_peaks_cwt(binz,np.arange(2,10),min_snr=0.5)
    print "peak",peak_atoms*dz+zmin
    division=np.zeros(len(peak_atoms))
    nplane=len(division)
    for i in range(len(peak_atoms)-1):
      division[i]=zmin+((peak_atoms[i]+peak_atoms[i+1])/2.)*dz
      print -(peak_atoms[i]-peak_atoms[i+1])*dz
    division[nplane-1]=zmax
    #print nplane
    pid=np.array(range(natom))
    #for each atom, find the plane id
    for i in range(natom):
      z=positions[i,direction]
      planeid=0
      while ((planeid<(nplane-1)) and (z>division[planeid])):
        planeid+=1
      pid[i]=int(planeid)
    symbol=np.array(range(nplane))

    return pid, nplane, symbol

#    def locater(array):
#      length=len(array)
#      middle=length/2
#      if array[middle] < a[middle - 1]:
#          #only look at the left 1 ... n/2 - 1
#      elif array[middle] < a[middle + 1]:
#          #then only look at the right n/2 +1 ... n
#      else:
#          #n/2 is a peak
#
