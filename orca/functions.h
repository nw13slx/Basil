//this code is used to analyse 
//output of ORCA. 
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

#define MAX_COORD 1000

//set up the buffer size for reading
#define MAX_CHARACTER 1000
#define MAX_COLUMN 50
#define MAX_ELEMENT 10
#define DEFAULT_CUTOFF 3.2
#define MAX_ATOMS 1000

#define ANG2BOHR 0.52917721092
const double Eh2eV=27.2113834;
const double bohr2a=0.529177249;
const double Ehau2eVa=Eh2eV/bohr2a;
const double kcal2eV=4.3363e-2;

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

int find_pattern(ifstream &In1,string pattern,const char * mode=(const char *)"contains",bool back=false){
  char temp[MAX_CHARACTER], * pch;
  int beg_pos=In1.tellg();
  bool find=false;
  int beg_line=0;
  if (strcmp(mode,"contains")==0){
    while (In1.good() && !find){
      beg_line=In1.tellg();
      In1.getline(temp,MAX_CHARACTER); 
      pch=strstr(temp,pattern.c_str());
      if (pch!=NULL){
        find=true;
      }
    }
  }else{
    while (In1.good() && !find){
      beg_line=In1.tellg();
      In1.getline(temp,MAX_CHARACTER); 
      //if (strcmp(temp,pattern)==0) find=true;
      if (pattern.compare(temp)==0) find=true;
    }
  }
  if (back){
      In1.seekg(beg_line);
  }
  if (!find){
      In1.clear();
      In1.seekg(beg_pos);
      beg_line=-1;
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

void parse(char * temp, int & column, string *content){
  char temp0[MAX_CHARACTER];
  strcpy(temp0,temp);
  column=0;
  char * pch;
  pch = strtok (temp0," ");
  while ((pch != NULL)&&(column<MAX_COLUMN)) {
      content[column]=pch;
      column++;
      pch = strtok (NULL, " ");
  }
  pch= NULL;
}


bool read_orbital(ifstream &In1,int *nstate,double *energy,double *occupancy,double *homo,double *lumo,int *homo_i, int *lumo_i){
  char temp[MAX_CHARACTER], * pch;
  string *content=new string[MAX_COLUMN];
  int column, line=0;
  bool Print_MO=false;
  
  string pattern="ORBITAL ENERGIES";
  int beg_pos=In1.tellg();
  int v_pos=find_pattern(In1,pattern,"contains",true);
  if (v_pos!=-1){
    Print_MO=true;
  }else{
      cout<<"ERROR: ORBITAL ENERGIES block not found"<<endl;
      return false;
  }

  nstate[0]=0;
  nstate[1]=0;
  
  //spin up and down
  lumo[0]=INFINITY;
  lumo[1]=INFINITY;
  homo[0]=-INFINITY;
  homo[1]=-INFINITY;
  homo_i[0]=0;
  homo_i[1]=0;
  lumo_i[0]=0;
  lumo_i[1]=0;
  In1.getline(temp,MAX_CHARACTER); //another buffering line
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
        column=break_line(temp,content);
        if (column == 4){
          int number=atoi(content[0].c_str());
          double occ=atof(content[1].c_str());
          double e=atof(content[3].c_str());
          p_energy[nstate[spin]]=e;
          p_occupancy[nstate[spin]]=occ;
          if (occ>0 && e>homo[spin]){
            homo[spin]=e;
            homo_i[spin]=number;
          }else{
            if (occ==0 && homo[spin]<e && lumo[spin]==INFINITY){ //&& lumo[spin]>e){
              lumo[spin]=e;
              lumo_i[spin]=number;
            }
          }
          nstate[spin]++; 
        }
    }while (column==4 );
  }
  return true;
}
