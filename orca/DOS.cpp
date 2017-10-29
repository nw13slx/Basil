//this code is used to analyse the Print[ P_MOs ] 1
//output of ORCA. The density of state will be decomposed
//to different element and different angular momentum , s p d f
//this is called PDOS in the solid state dft community
//usage: MO <inputfile name> <outputfile name> 
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
#define MAX_ENERGYLINE 1000

//set up the final plotting region
//the fermi level is shfited to zero
#define EMIN -30
#define EMAX 15

//set up the buffer size for reading
#define MAX_CHARACTER 1000
#define MAX_COLUMN 20

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

int main(int argc, char **argv){

    char temp[MAX_CHARACTER], * pch,content[MAX_COLUMN][MAX_CHARACTER];
    int column, line=0;
    bool Print_MO=false;
    
    //locate the beginning of the print_MO
    ifstream In1(argv[1]);
    while (In1.good() && Print_MO==false ){
        line++;
        In1.getline(temp,MAX_CHARACTER); 

        column=0;
        pch = strtok (temp," ");
        while ((pch != NULL)&&(column<MAX_COLUMN)) {
            strcpy(content[column],pch);
            column++;
            pch = strtok (NULL, " ");
        }
        pch= NULL;

        if (column == 2 ){
            if (( strcmp(content[0],"ORBITAL") == 0 ) && (strcmp(content[1],"ENERGIES") == 0)){
                Print_MO=true;
            }
        }
    }
    if ((! In1.good()) || (Print_MO==false)){
        cout<<"ERROR: ORBITAL ENERGIES block not found"<<endl;
        return 1;
    }

    double *energy=new double[2*MAX_ENERGYLINE];
    double *occupancy=new double[2*MAX_ENERGYLINE];
    int nstate[2]={0,0};
    double efermi=-100000;
    
    //spin up and down
    In1.getline(temp,MAX_CHARACTER); //another buffering line
    for (int spin=0;spin<2;spin++){
      line++;
      In1.getline(temp,MAX_CHARACTER); //another buffering line
      line++;
      In1.getline(temp,MAX_CHARACTER); //another buffering line
      line++;
      double *p_energy=energy+spin*MAX_ENERGYLINE;
      double *p_occupancy=occupancy+spin*MAX_ENERGYLINE;
      bool ORBITAL=true;
      do{
          In1.getline(temp,MAX_CHARACTER); 
          line++;
          column=0;
          pch = strtok (temp," ");
          while (pch != NULL) {
              strcpy(content[column],pch);
              column++;
              pch = strtok (NULL, " ");
          }
          pch= NULL;

          if (column == 4){
            p_energy[nstate[spin]]=atof(content[3]);
            p_occupancy[nstate[spin]]=atof(content[1]);
            if ((p_occupancy[nstate[spin]]>0) && (p_energy[nstate[spin]]>efermi)){
                efermi=p_energy[nstate[spin]];
            }
            nstate[spin]++; 
          }
      }while (column==4 );
    }
    cout<<"done reading"<<endl;

    ofstream out(argv[2]);
    // output the raw pdos without smearing
    out<<"# efermi:"<<efermi<<endl;;
    out<<"# energy up down sum " <<endl;

    //smearing
    double dE=0.1;
    double sigma=0.05;
    double sigma2=sigma*sigma;

    double base=efermi+EMIN;
    double up=efermi+EMAX;
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
        if ((p_energy[i]-efermi)>=EMIN && (p_energy[i]-efermi) <=EMAX){
          gaussian(grid, x, spreadE,p_energy[i],sigma2);
          for (int igrid=0; igrid < grid; igrid++)
            p_smearing[igrid]+=spreadE[igrid];
        }
      }
    }
    cout<<"done smearing"<<endl;

    for (int i=0;i<grid;i++)
        out<<x[i]-efermi<<" "<<smearing[i]<<" "<<-smearing[i+grid]<<endl;
}
