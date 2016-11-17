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
      if not self.control.run_pdos:
         return 0
      dos=self.dos
      atom=self.atom
      loc_up=dos.loc_up
      loc_down=dos.loc_down
      control=self.control
      simple=self.simple

    #plot each orbital
      if any(x!=-1 for x in dos.par_element):
        par_symbol=self.par_symbol
        y=[]
        line_label=[]
        line_color=['y','y','r','r','b','b','k','k']
        for i in range(4):
          if (dos.par_orbital[i]!=None):
            y+=[dos.par_orbital[i][:,0],dos.par_orbital[i][:,1]]
          else:
            y+=[None,None]
          if atom.symbol:
            line_label+=[atom.symbol[dos.par_element[i]]+"_"+par_symbol[i],None]
          else:
            line_label+=["atomtype="+str(dos.par_element[i])+" "+par_symbol[i]+"_orbital",None]
        simple("pdf_orbital",dos.Xenergy,y,loc_down,loc_up,line_color,line_label)
        del y, line_label, line_color
    
      if (dos.par_element[2] >=0):
        y=[]
        line_label=['_t2g',None,'_eg',None]
        line_color=['r','r','b','b']
        y=[d_t2g[:,0],d_t2g[:,1],d_eg[:,0],d_eg[:,1]]
        simple("d_decompose",dos.Xenergy,y,loc_down,loc_up,line_color,line_label)
        del y, line_label, line_color
    
      if self.control.perspecies:
        y=[]
        line_label=[]
        line_color=[]
        for i in range(atom.ntype):
          y+=[dos.perspecies[:,i]]
          line_label+=["species "+str(i)]
          line_color+=[None]
        if atom.symbol:
          line_label+=[atom.symbol[i]]
        simple("perspecies",dos.Xenergy,y,loc_down,loc_up,line_color,line_label)
        del y, line_label, line_color
    
  def tot_dos(self):
      y=[self.dos.dos0[:,0],self.dos.dos0[:,1]]
      line_label=[None,None]
      line_color=['r','b']
      self.simple("DOS-tot",self.dos.Xenergy,y,self.dos.loc_down,self.dos.loc_up,line_color,line_label)
      del y, line_label, line_color

  def simple(self,name,x,y,l_d,l_u,line_color,line_label):
    a=plt.figure(0)
    for i in range(len(y)):
      if (y[i]!=None):
        line=plt.plot(x[l_d:l_u],y[i][l_d:l_u])
        if (line_color[i]!=None):
          plt.setp(line,color=line_color[i])
        if (line_label[i]!=None):
          plt.setp(line,label=line_label[i])
    if any(x for x in line_label):
      plt.legend()
    if (name=="DOS-tot"):
        self.dos0_extra()
    pl.savefig(self.control.name+name+".png")
    plt.close()
  
  def atom_start(self):
      pass
  def atom(self,atomi):
      pass

  def atom_end(self):
      pass
  def dos0_extra(self):
    dos0=self.dos.dos0
    ef=self.dos.Xenergy
    fermiN=self.dos.fermiN
    plt.fill_between(ef[:fermiN], 0, dos0[:fermiN,0],facecolor='red')
    plt.fill_between(ef[:fermiN], 0, dos0[:fermiN,1],facecolor='blue')

