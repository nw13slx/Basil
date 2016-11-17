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
      if all(x!=-1 for x in dos.par_element):
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
    
  def plot_tot_dos(self):
      y=[self.dos.dos0[:,0],self.dos.dos0[:,1]]
      line_label=[None,None]
      line_color=['red','blue']
      simple_line_plot("DOS-tot",dos.Xenergy,y,loc_down,loc_up,line_color,line_label)
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
    if (name=="DOS-tot"):
        dos0_extra()
    pl.savefig(self.control.name+name+".png")
    plt.close()
  
  def plot_atom_start(self):
      pass
  def plot_atom(self,atomi):
      pass

  def plot_atom_end(self):
      pass
  def dos0_extra(self):
    dos0=self.dos.dos0
    ef=self.dos.Xenergy
    fermiN=self.dos.fermiN
    plt.fill_between(ef[:fermiN], 0, dos0[:fermiN,0],facecolor='red')
    plt.fill_between(ef[:fermiN], 0, dos0[:fermiN,1],facecolor='blue')

