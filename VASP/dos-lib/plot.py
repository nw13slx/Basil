import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl

class plot:
  cm = pl.get_cmap('winter')
  plot_atom_dict={}
  def __init__(self,control):
      self.control=control
      self.atom=control.atom
      self.dos=control.dos
      self.par_symbol=['s','p','d','f']
  def pdos(self):
      dos=self.dos
      atom=self.atom
      loc_up=dos.loc_up
      loc_down=dos.loc_down
      control=self.control
      simple_line_plot=self.simple_line_plot

      par_symbol=self.par_symbol
      if atom.symbol:
        for i in range(4):
          par_symbol[i]=atom.symbol[dos.par_element[i]]+"_"+par_symbol[i]
      else:
        for i in range(4):
          par_symbol[i]="atomtype="+str(dos.par_element[i])+" "+par_symbol[i]+"_orbital"
    
    #plot each orbital
      y=[]
      line_label=[]
      line_color=['y','y','r','r','b','b','k','k']
      for i in range(4):
        if (dos.par_orbital[i]!=None):
          y+=[dos.par_orbital[i][:,0],dos.par_orbital[i][:,1]]
        else:
          y+=[None,None]
        line_label+=[par_symbol[i],None]
      simple_line_plot("pdf_orbital",dos.Xenergy,y,loc_down,loc_up,line_color,line_label)
      del y, line_label, line_color
    
      if (dos.par_element[2] >=0):
        y=[]
        line_label=['_t2g',None,'_eg',None]
        line_color=['r','r','b','b']
        y=[d_t2g[:,0],d_t2g[:,1],d_eg[:,0],d_eg[:,1]]
        simple_line_plot("d_decompose",dos.Xenergy,y,loc_down,loc_up,line_color,line_label)
        del y, line_label, line_color
    
      y=[]
      line_label=[]
      line_color=[]
      for i in range(atom.ntype):
        y+=[dos.perspecies[:,i]]
        line_label+=["species "+str(i)]
        line_color+=[None]
      if atom.symbol:
        line_label+=[atom.symbol[i]]
      simple_line_plot("perspecies",dos.Xenergy,y,loc_down,loc_up,line_color,line_label)
      del y, line_label, line_color
    
  def simple_line_plot(self,name,x,y,l_d,l_u,line_color,line_label):
    a=plt.figure(0)
    for i in range(len(y)):
      if (y[i]!=None):
        if ((line_color[i]!=None) and line_label[i]!=None):
          plt.plot(x[l_d:l_u],y[i][l_d:l_u],c=line_color[i],label=line_label[i])
        elif ((line_color[i]==None) and line_label[i]!=None):
          plt.plot(x[l_d:l_u],y[i][l_d:l_u],label=line_label[i])
        elif ((line_color[i]!=None) and line_label[i]==None):
          plt.plot(x[l_d:l_u],y[i][l_d:l_u],c=line_color[i])
        else:
          plt.plot(x[l_d:l_u],y[i][l_d:l_u])
    plt.legend()
    pl.savefig(self.control.name+name+".png")
    plt.close()
  
  def plot_atom_start(self):
    rad=12
    scale=1/23.*2.
    atom=self.atom
    print atom.positions[:,2]
    zmin=np.amin(atom.positions[:,2])
    zmax=np.amax(atom.positions[:,2])
    lz=zmax-zmin
    zcenter=(zmin+zmax)*0.5
    Oplot=200
    Ceplot=300
    for i in ('zmin', 'zmax', 'lz','zcenter','Oplot','Ceplot','scale'):
       self.plot_atom_dict[i] = locals()[i]

  def plot_atom(self,atomi):
    Oplot=self.plot_atom_dict['Oplot']
    Ceplot=self.plot_atom_dict['Ceplot']
    zcenter=self.plot_atom_dict['zcenter']
    zmin=self.plot_atom_dict['zmin']
    scale=self.plot_atom_dict['scale']
    lz=self.plot_atom_dict['lz']

    atom=self.atom
    dos=self.dos

    y=dos.tot[dos.loc_down:dos.loc_up,0]+dos.tot[dos.loc_down:dos.loc_up,1]
    color =self.cm(2.*abs(atom.positions[atomi,2]-zcenter)/lz)  # color will now be an RGBA tuple
    y_off=(atom.positions[atomi,2]-zmin)/lz
    if (atom.species[atomi]==dos.par_element[1]):
      plt.figure(Oplot)
      ymax=np.amax(y)
      pl.plot(dos.Xenergy[dos.loc_down:dos.loc_up],y_off+y/ymax*scale,c=color)
    if (atom.species[atomi]==dos.par_element[3]):
      plt.figure(Ceplot)
      ymax=np.amax(y)
      pl.plot(dos.Xenergy[dos.loc_down:dos.loc_up],y_off+y/ymax*scale,c=color)

  def plot_atom_end(self):
    Oplot=self.plot_atom_dict['Oplot']
    Ceplot=self.plot_atom_dict['Ceplot']
    zcenter=self.plot_atom_dict['zcenter']
    zmin=self.plot_atom_dict['zmin']
    scale=self.plot_atom_dict['scale']

    plt.figure(Oplot)
    plt.ylim(0,1+scale)
    plt.show()
    plt.savefig("DOS_O_INSIDE.png")
    plt.close()

    plt.figure(Ceplot)
    plt.ylim(0,1+scale)
    plt.savefig("DOS_Ce_INSIDE.png")
    plt.close()
