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
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <map>
#include "atomic_mass.h"
using namespace std;

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

const double Eh2eV = 27.2113834;
const double bohr2a = 0.529177249;
const double Ehau2eVa = Eh2eV/bohr2a;
const double kcal2eV = 4.3363e-2;
const double autime2s = 6.582119514e-16 / Eh2eV ;
const double borh_autime2a_s = bohr2a / autime2s ; 
const double eV = 1.6021765e-19;
const double kin2eV = 1.660539040e-27*(1e-10*borh_autime2a_s)*(1e-10*borh_autime2a_s) / eV;
const double kin2T = kin2eV / 8.6173303e-5;
const map<string,double> AUM2kg::atomic_mass =  AUM2kg::create_map();


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


double * read_vel(ifstream &fin, int &natom){

  char temp[100];
  double *T;
  double x[3];
  double m;
  string symbol;

  if (!fin.good()){
    return NULL;
  }
  fin >> natom;
  if (!fin.good()){
    return NULL;
  }
  T = new double [natom];
  fin.getline(temp,100);
  fin.getline(temp,100);
  for (int i = 0; i < natom ; i++){
    int i3=i*3;
    fin >> temp >> x[0] >> x[1] >> x[2];
    symbol = temp;
    m  = AUM2kg::atomic_mass.find(symbol)->second;
    fin.getline(temp,100);
    T[i] = 1/2. * m * (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) * kin2T;
    // cout<<i << " " << symbol<<" "<<m << " "<< T[i]<<endl;
  } 
  return T;
}
