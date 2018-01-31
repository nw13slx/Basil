//this code is used to analyse the Print[ P_MOs ] 1
//output of ORCA. The density of state will be decomposed
//to different element and different angular momentum , s p d f
//this is called PDOS in the solid state dft community
//usage: MO <inputfile name> <outputfile name> 
//author: Lixin Sun nw13mifaso@gmail.com

#include "functions.h"

int main(int argc, char **argv){

    char temp[MAX_CHARACTER], * pch;
    string *content=new string[MAX_COLUMN];
    int column, line=0;
    bool Print_MO=false;
    
    if (argc<3){
      cout<<"need two arguments: orca_dos intput output"<<endl;
      return 1;
    }
    //locate the beginning of the print_MO
    ifstream In1(argv[1]);
    double *energy=new double[2*MAX_ENERGYLINE];
    double *occupancy=new double[2*MAX_ENERGYLINE];
    int *nstate=new int[2];
    double *homo=new double[2];
    double *lumo=new double[2];
    bool read=read_orbital(In1,nstate,energy,occupancy,homo,lumo);
    if (read==false){
      return 1;
    }

    ofstream out(argv[2]);
    // output the raw pdos without smearing
    out<<"# homo "<<homo[0]<<" "<<homo[1]<<" lumo "<<lumo[0]<<" "<<lumo[1]<<endl;;
    out<<"# energy up down sum " <<endl;

    //smearing
    double dE=0.1;
    double sigma=0.05;
    double sigma2=sigma*sigma;

    double base=homo[0]+EMIN;
    double up=homo[0]+EMAX;
    int grid=int((up-base)/dE);

    double *smearing=new double [2*grid];
    fill(smearing, smearing + 2*grid , 0);

    double * spreadE=new double[grid];
    double x[grid];
    for (int i=0;i<grid;i++) x[i]=base+i*dE;

    for (int spin=0; spin<2; spin++){
      double *p_energy=energy+spin*MAX_ENERGYLINE;
      double *p_smearing=smearing+spin*grid;
      for (int i=0;i<nstate[spin];i++){
        if ((p_energy[i]-homo[0])>=EMIN && (p_energy[i]-homo[0]) <=EMAX){
          gaussian(grid, x, spreadE,p_energy[i],sigma2);
          for (int igrid=0; igrid < grid; igrid++)
            p_smearing[igrid]+=spreadE[igrid];
        }
      }
    }
    cout<<"done smearing"<<endl;

    for (int i=0;i<grid;i++)
        out<<x[i]-homo[0]<<" "<<smearing[i]<<" "<<-smearing[i+grid]<<endl;
}
