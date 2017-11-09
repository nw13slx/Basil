//average the local potential file along certain direction
//usage: average_loc input output direction<1-3>

//revised from a source code called vtotav.f. don't remember where did I find the code...
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

#define MAX_CHARACTER 1000
#define MAX_COLUMN 200

void parse(char * temp, int & column, string *content){
  column=0;
  char * pch;
  pch = strtok (temp," ");
  while ((pch != NULL)&&(column<MAX_COLUMN)) {
      content[column]=pch;
      column++;
      pch = strtok (NULL, " ");
  }
  pch= NULL;
}

int main(int argc, char **argv){
  if (argc < 2){
    cout<<" FAILED: need more input arguments"<<endl;
    return 1;
  }
  ifstream fin(argv[1]); 
  ofstream fout(argv[2]); 
  int ider=atoi(argv[3]);

  if (ider<1 || ider >3){
    cout<<" FAILD: direction has to be 1-3 for x-z"<<endl;
    return 1;
  } 
  ider -=1;

  char temp[MAX_CHARACTER];
  string content[MAX_COLUMN];
  int column, line=0;

  bool find_grid=false;
  do{
    fin.getline(temp,MAX_CHARACTER); 
    parse(temp,column,content);
    if (column==6){
        if ((content[0]=="Electrostatic") && (content[1]=="potential") && (content[4]=="grid")){
            find_grid=true;
        }
    }
  }while (find_grid==false &&fin.good());
  if (find_grid==false){
    cout<<"intput script need:\n  potgrid nx ny nz"<<endl;
    return 1;
  }

  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  
  int ng[3];
  fin>>temp>>temp>>temp>>ng[0]>>temp>>ng[1]>>temp>>ng[2];
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 


  int nlocal=ng[0]*ng[1]*ng[2];
  double * vlocal=new double[nlocal];
  for (int i=0;i<nlocal;i++) {
      fin >>temp>>temp>>temp>>vlocal[i];
      fin.getline(temp,MAX_CHARACTER); 
  }
  cout<<"read density"<<endl;

  double *vav = new double [ng[ider]];
  fill(vav, vav+ng[ider], 0);

  double scale=1./float(nlocal/ng[ider]);
  if (ider == 0) {
    for (int ix=0; ix<ng[0]; ix++){
      for (int iz=0; iz<ng[2]; iz++){
        for (int iy=0; iy<ng[1]; iy++){
          int ipl=iz+(iy+ix*ng[1])*ng[2];
          vav[ix] += vlocal[ipl]*scale;
        }
      }
    }
  }else if ( ider == 1){
    for (int iy=0; iy<ng[1]; iy++){
      for (int iz=0; iz<ng[2]; iz++){
        for (int ix=0; ix<ng[0]; ix++){
          int ipl=iz+(iy+ix*ng[1])*ng[2];
          vav[iy] += vlocal[ipl]*scale;
        }
      }
    }
  }else if (ider==2){
    for (int iz=0; iz<ng[2]; iz++){
      for (int iy=0; iy<ng[1]; iy++){
        for (int ix=0; ix<ng[0]; ix++){
          int ipl=iz+(iy+ix*ng[1])*ng[2];
          vav[iz] += vlocal[ipl]*scale;
        }
      }
    }
  }else { 
   cout<<"FAILED: ider should be 1-3"<<endl;
   return 1;
  } 

  fout<<ng[ider]<<" "<<ider+1<<endl;
  for (int i=0; i<ng[ider];i++){
    fout<<i<<" "<<vav[i]<<endl;
  }

}

