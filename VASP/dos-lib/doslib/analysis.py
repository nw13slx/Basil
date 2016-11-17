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

    EGMIN_u=30
    EGMIN_d=30
    EGMAX_u=-30
    EGMAX_d=-30
    foundEg=0
    finfo = open("band-info", 'w')
    finfo.write('between 0 %15.8f %15.8f \n' % (ef[fermiN],ef[fermiN+1]))
    for n in range(fermiN, 0, -1):
       if ((dos0[n,0]!=0) and (EGMIN_u >0)):
          EGMIN_u=ef[n]
       if ((dos0[n,1]!=0) and (EGMIN_d >0)):
          EGMIN_d=ef[n]
    for n in range(fermiN, nedos, 1):
       if ((dos0[n,0]!=0) and (EGMAX_u <0)):
          EGMAX_u=ef[n]
       if ((dos0[n,1]!=0) and (EGMAX_d <0)):
          EGMAX_d=ef[n]
    if (EGMIN_u >EGMIN_d):
        EGMIN = EGMIN_u
    else:
        EGMIN = EGMIN_d
    if (EGMAX_u < EGMAX_d):
        EGMAX = EGMAX_u
    else:
        EGMAX = EGMAX_d
    if ((dos0[fermiN,0]!=0) and (dos0[fermiN+1,0]!=0)):
        Eg_u=0
    else:
        Eg_u = EGMAX_u-EGMIN_u
    if ((dos0[fermiN,1]!=0) and (dos0[fermiN+1,1]!=0)):
        Eg_d=0
    else:
        Eg_d = EGMAX_d-EGMIN_d
    finfo.write('Eg_u %15.8f %15.8f %15.8f\n' % (EGMIN_u,EGMAX_u,Eg_u))
    finfo.write('Eg_d %15.8f %15.8f %15.8f\n' % (EGMIN_d,EGMAX_d,Eg_d))
    finfo.write('Eg %15.8f %15.8f %15.8f\n' % (EGMIN,EGMAX,EGMAX-EGMIN))
    finfo.close()


  def peak_finder(self,center,span):
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
    ana_dos=dos0[win_down:win_up,0]+dos0[win_down:win_up,1]

    #assuming there is only one gaussian peak
    avg=(np.arange(win_span).dot(ana_dos))/np.sum(ana_dos)
    sigma=(np.arange(win_span)-avg)
    sigma=np.square(sigma)
    sigma=sigma.dot(dos0[win_down:win_up,0])/np.sum(ana_dos)
    sigma=np.sqrt(sigma)

    x0=dos.Xenergy[0]
    dx=dos.Xenergy[1]-x0
    print avg*dx+x0,sigma*dx
    return avg*dx+x0,sigma*dx

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
