//this code is used to turn the screen output of quantum espresso to json file format
//usage: qe_out2json <inputfile name> <outputfile name> <optional: header for wannier xyz file>
//author: Lixin Sun nw13mifaso@gmail.com

#include <algorithm>
#include <iostream>
#include <locale>
#include <fstream>          // file I/O suppport
#include <cstdlib>          // support for exit()
#include <stdio.h>
#include <sys/timeb.h>
#include <sys/types.h>
#include <time.h>
#include <cmath>
#include <malloc.h>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <sstream>
using namespace std;

#include "math.h"
#include "stdlib.h"       // for random
#include <vector>

//set up the size of array used for storage
#define MAX_ELEMENT 10
#define MAX_M   6
#define MAX_ENERGYLINE 4000
#define TEMP_ENERGYLINE 10

//set up the final plotting region
//the fermi level is shfited to zero
#define EMIN -30
#define EMAX 15

//set up the buffer size for reading
#define MAX_CHARACTER 1000
#define MAX_COLUMN 20

const double kbar2GPa=0.1;
const double Eh2eV=27.2113834;
const double Ryd2eV=27.2113834/2.;
const double bohr2a=0.529177249;
const double Ehau2eVa=Eh2eV/bohr2a;
const double Rydau2eVa=Ryd2eV/bohr2a;
const double kcal2eV=4.3363e-2;

double distance(double* cell,double *x1,double *x2,double *dx){
  dx[0]=x1[0]-x2[0];
  dx[1]=x1[1]-x2[1];
  dx[2]=x1[2]-x2[2];
  double hcell[3]={cell[0]/2.,cell[4]/2.,cell[8]/2.};
  double wcell[3]={cell[0],cell[4],cell[8]};
  for (int d=0;d<3;d++){
    while (dx[d]<-hcell[d]){
      dx[d]+=wcell[d];
    }
    while (dx[d]>hcell[d]){
      dx[d]-=wcell[d];
    }
  }
  return sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

}

bool isxyz(char c){
    if ((c=='x')||(c=='y')||(c=='z')) return true;
    else return false;
}

bool isdigit1(char c){
    if ((c>='0')&&(c<='9')) return true;
    else return false;
}

void gaussian(int grid, double *x, double *y, double x0, double sigma2){
    for (int i=0;i<grid;i++) y[i]=exp(-(x[i]-x0)*(x[i]-x0)/2./sigma2);
}

bool contain_alphabet(char * c){
    for (int i=0;i<strlen(c);i++){
        if (((c[i]>='a')&&(c[i]<='z'))||((c[i]>='A')&&(c[i]<='Z')))
            return true;
    }
    return false;
}

int find_pattern(ifstream &In1,const char *pattern,const char * mode=(const char *)"contains",bool back=false){
  char temp[MAX_CHARACTER], * pch;
  int beg_pos=In1.tellg();
  bool find=false;
  int beg_line=0;
  if (strcmp(mode,"contains")==0){
    while (In1.good() && !find){
      beg_line=In1.tellg();
      In1.getline(temp,MAX_CHARACTER); 
      pch=strstr(temp,pattern);
      if (pch!=NULL){
        find=true;
      }
    }
  }else{
    while (In1.good() && !find){
      beg_line=In1.tellg();
      In1.getline(temp,MAX_CHARACTER); 
      if (strcmp(temp,pattern)==0) find=true;
    }
  }
  if (!find){
      In1.clear();
      In1.seekg(beg_pos);
      beg_line=NULL;
  }
  if (back){
    In1.seekg(beg_line);
  }
  return beg_line;
}

int break_line(char *temp,string *content){
  char *pch = strtok (temp," ");
  int column=0;
  while (pch != NULL) {
      content[column]=string(pch);
      column++;
      pch = strtok (NULL, " ");
  }
  pch= NULL;
  return column;
}

int main(int argc, char **argv){

    char temp[MAX_CHARACTER], * pch;
    string* content=new string[MAX_COLUMN];
    int column, line=0;
    bool Print_MO=false;
    
    //locate the beginning of the print_MO
    ifstream In1(argv[1]);
    if ( !In1.good()){
        cout<< " the input file does not exist or is corrupted..."<<endl;
        return 1;
    }
    ofstream json_o(argv[2]);
    if ( !json_o.good()){
        cout<< " the output file does not exist or is corrupted..."<<endl;
        return 1;
    }

    char *pattern="Program PWSCF";
    int v_pos=find_pattern(In1,pattern,"contains",true);
    if (v_pos!=NULL){
    json_o.precision(5);
    json_o<<"{"<<endl;
    json_o<<"\"outputfile\":\""<<argv[1]<<"\","<<endl;
    json_o<<"\"software\":\"QE";
    In1.getline(temp,MAX_CHARACTER); 
    break_line(temp,content);
    json_o<<content[2]<<"_"<<content[3]<<"_"<<content[4]<<"_"<<content[5]<<"\","<<endl;
    }else{
      cout<<"This is not a QE output file"<<endl;
      return 0;
    }


    pattern="bravais-lattice index";
    int u_pos=find_pattern(In1,pattern,"contains",true);
    int ibrav=-1;
    if (u_pos!=NULL){
      In1.getline(temp,MAX_CHARACTER);
      column=break_line(temp,content);
      ibrav=atoi(content[3].c_str());
    }

    pattern="lattice parameter";
    u_pos=find_pattern(In1,pattern,"contains",true);
    double alat=-1;
    if (u_pos!=NULL){
      In1.getline(temp,MAX_CHARACTER);
      column=break_line(temp,content);
      alat=bohr2a*atof(content[4].c_str());
    }

    pattern="number of atoms/cell";
    u_pos=find_pattern(In1,pattern,"contains",true);
    In1.getline(temp,MAX_CHARACTER);
    column=break_line(temp,content);
    int atomn=atoi(content[4].c_str());

    pattern="celldm(1)";
    u_pos=find_pattern(In1,pattern,"contains",true);
    double cell[9]={0,0,0,0,0,0,0,0,0};
    if (u_pos!=NULL){
      In1.getline(temp,MAX_CHARACTER);
      column=break_line(temp,content);
      double a=atof(content[1].c_str());
      double a2=atof(content[3].c_str());
      double a3=atof(content[5].c_str());
      In1.getline(temp,MAX_CHARACTER);
      column=break_line(temp,content);
      double a4=atof(content[1].c_str());
      double a5=atof(content[3].c_str());
      double a6=atof(content[5].c_str());
      if (ibrav==1){
        cell[0]=a*bohr2a;
        cell[4]=a*bohr2a;
        cell[8]=a*bohr2a;
      }
      json_o<<"\"cell\":[";
      for (int i=0;i<9;i++){
        if (i>0){
          json_o<<", ";
        }
        json_o<<cell[i];
      }
      json_o<<"],"<<endl;
    }

    json_o<<"\"comp_info\":[";
    pattern="PseudoPot";
    int len_pat=strlen(pattern);
    u_pos=find_pattern(In1,pattern,"contains",true);
    string atom_species[100];
    int Zval[100];
    int ntype=0;

    while (u_pos!=NULL && In1.good()){
      In1.getline(temp,MAX_CHARACTER); 
      column=break_line(temp,content);
      atom_species[ntype]=content[4];
      In1.getline(temp,MAX_CHARACTER); 
      column=break_line(temp,content);
      if (ntype>0){
        json_o<<", ";
      }
      json_o<<"\""<<content[0]<<"\"";
      find_pattern(In1,"Zval","contains",true);
      In1.getline(temp,MAX_CHARACTER); 
      column=break_line(temp,content);
      Zval[ntype]=atoi(content[column-1].c_str());
      ntype++;
      u_pos=find_pattern(In1,pattern,"contains",true);
    }
    json_o<<"],"<<endl;

    pattern="LDA+U calculation ";
    u_pos=find_pattern(In1,pattern);
    if (u_pos!=NULL){
      In1.getline(temp,MAX_CHARACTER);
      column=break_line(temp,content);
      json_o<<"\"dft+u\":[";
      for (int i=0;i<column;i++){
        if (i>0){
          json_o<<", ";
        }
        json_o<<"\""<<content[i]<<"\"";
      }
      json_o<<"],"<<endl;
      In1.getline(temp,MAX_CHARACTER);
      column=break_line(temp,content);
      json_o<<"\"dft+u_para\":[";
      for (int i=0;i<column;i++){
        if (i>0){
          json_o<<", ";
        }
        json_o<<"\""<<content[i]<<"\"";
      }
      json_o<<"],"<<endl;
    }

    pattern="Cartesian axes";
    int c_pos=find_pattern(In1,pattern);
    double *atomx=new double [3*atomn];
    int *atomt=new int [atomn];
    if (c_pos!=NULL){
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      json_o<<"\"coordinates\":[";
      char species[100000]="\"species\":[";
      for (int i=0;i<atomn;i++){
        In1.getline(temp,MAX_CHARACTER);
        column=break_line(temp,content);
        if (i>0){
          json_o<<", ";
          strcat(species,",");
        }
        int i3=i*3;
        atomx[i3]=atof(content[6].c_str())*alat;
        atomx[i3+1]=atof(content[7].c_str())*alat;
        atomx[i3+2]=atof(content[8].c_str())*alat;
        json_o<<"["<<atomx[i3]<<","<<atomx[i3+1]<<","<<atomx[i3+2]<<"]";
        strcat(species,"\"");
        strcat(species,content[1].c_str());
        for (int t=0;t<ntype;t++){
          if (content[1]==atom_species[t]){
            atomt[i]=t;
          }
        }
        strcat(species,"\"");
      }
      json_o<<"],"<<endl;
      strcat(species,"],");
      json_o<<species<<endl;
      json_o<<"\"atomn\":"<<atomn<<","<<endl;
    }

    pattern="!";
    json_o<<std::fixed;
    int fe_pos=find_pattern(In1,pattern,"contains",true);
    if (fe_pos!=NULL){
      In1.seekg(fe_pos);
      In1.getline(temp,MAX_CHARACTER);
      break_line(temp,content);
      json_o<<"\"energy\":"<<atof(content[4].c_str())*Ryd2eV<<","<<endl;
    }

    pattern="total magnetization";
    int spin_pos=find_pattern(In1,pattern,"contains",true);
    if (spin_pos!=NULL){
      int s2_pos=find_pattern(In1,"Expectation value of");
      In1.getline(temp,MAX_CHARACTER);
      break_line(temp,content);
      json_o<<"\"S\":"<<content[3]<<","<<endl;
    }

    pattern="Forces acting on atoms";
    c_pos=find_pattern(In1,pattern);
    if (c_pos!=NULL){
      In1.getline(temp,MAX_CHARACTER);
      int catomn=0;
      json_o<<"\"forces\":[";
      for (int i=0;i<atomn;i++){
        In1.getline(temp,MAX_CHARACTER);
        column=break_line(temp,content);
        if (i>0){
          json_o<<", ";
        }
        json_o<<"["<<atof(content[6].c_str())*Rydau2eVa<<","<<atof(content[7].c_str())*Rydau2eVa<<","<<atof(content[8].c_str())*Rydau2eVa<<"]";
      }
      json_o<<"],"<<endl;
    }

    pattern="kbar";
    c_pos=find_pattern(In1,pattern);
    if (c_pos!=NULL){
      int catomn=0;
      json_o<<"\"stresses\":[";
      for (int i=0;i<3;i++){
        if (i>0) json_o<<", ";
        In1.getline(temp,MAX_CHARACTER);
        column=break_line(temp,content);
        json_o<<"["<<atof(content[3].c_str())*kbar2GPa<<","<<atof(content[4].c_str())*kbar2GPa<<","<<atof(content[5].c_str())*kbar2GPa<<"]";
      }
      json_o<<"],"<<endl;
    }

    //wannier part
    char wannier_xyz[100];
    if (argc>3){
      strcpy(wannier_xyz,argv[3]);
    }else{
      strcpy(wannier_xyz,&argv[1][6]);
    }
    char name[2][100];
    strcpy(name[0],wannier_xyz);
    strcat(name[0],"-u_centres.xyz");
    strcpy(name[1],wannier_xyz);
    strcat(name[1],"-d_centres.xyz");
    cout<<"wannier file"<<name[0]<<endl;

    bool winfail=false;
    double *dipole=new double [atomn*3];
    int *e=new int [atomn];
    int *spin=new int [atomn];
    for (int i=0;i<2;i++){
      ifstream win(name[i]);
      int electron=0;
      if (win.good()){
        win>>electron;
        electron-=atomn;
        win.getline(temp,MAX_CHARACTER);
        win.getline(temp,MAX_CHARACTER);
        cout<<temp<<endl;
        for (int el=0;el<electron;el++){
          win>>temp;
          double xx[3];
          win>>xx[0]>>xx[1]>>xx[2];
          int centeri=-1;
          double rmin=100;
          double mindx[3];
          for (int atom=0;atom<atomn;atom++){
            double dx[3];
            double dr=distance(cell,&atomx[atom*3],xx,dx);
            if (dr<rmin){
              rmin=dr;
              centeri=atom;
              mindx[0]=dx[0];
              mindx[1]=dx[1];
              mindx[2]=dx[2];
            }
          }
          if (centeri!=-1){
            int centeri3=centeri*3;
            dipole[centeri3]-=mindx[0];
            dipole[centeri3+1]-=mindx[1];
            dipole[centeri3+2]-=mindx[2];
            e[centeri]-=1;
            spin[centeri]+=(i-0.5)*2;
          }
        }
        win.close();
      }else{
        winfail=true;
      }
    }
    for (int atom=0;atom<atomn;atom++){
      e[atom]+=Zval[atomt[atom]];
    }
    if (!winfail){
      json_o<<"\"dipole_atom\":[";
      for (int i=0;i<atomn;i++){
        int i3=i*3;
        if (i>0){
          json_o<<", ";
        }
        json_o<<"["<<dipole[i3]<<","<<dipole[i3+1]<<","<<dipole[i3+2]<<"]";
        cout<<dipole[i3]<<","<<dipole[i3+1]<<","<<dipole[i3+2]<<" "<<e[i]<<endl;
      }
      json_o<<"],"<<endl;
      json_o<<"\"wannier_charge\":[";
      for (int i=0;i<atomn;i++){
        if (i>0){
          json_o<<", ";
        }
        json_o<<e[i];
      }
      json_o<<"],"<<endl;
      json_o<<"\"wannier_spin\":[";
      for (int i=0;i<atomn;i++){
        if (i>0){
          json_o<<", ";
        }
        json_o<<spin[i];
      }
      json_o<<"],"<<endl;
    }

    json_o<<"}"<<endl;
}

