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

int main(int argc, char **argv){
  if (argc < 3){
    cout<<" FAILED: need more input arguments"<<endl;
    cout<<" locave input output direction" <<endl;
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

  char temp[MAX_CHARACTER], * pch,content[MAX_COLUMN][MAX_CHARACTER];
  int column, line=0;
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 

  fin.getline(temp,MAX_CHARACTER); 
  cout<<temp<<endl;
 
  column=0;
  pch = strtok (temp," ");
  while ((pch != NULL)&&(column<MAX_COLUMN)) {
      strcpy(content[column],pch);
      column++;
      pch = strtok (NULL, " ");
  }
  pch= NULL;

  double natom=0;
  for (int i=0;i<column;i++){
    natom+=atoi(content[i]);
  }
  cout << natom<<endl;;

  fin.getline(temp,MAX_CHARACTER); 
  column=0;
  pch = strtok (temp," ");
  while ((pch != NULL)&&(column<MAX_COLUMN)) {
      strcpy(content[column],pch);
      column++;
      pch = strtok (NULL, " ");
  }
  pch= NULL;
  if ((content[0][0]=='s') || (content[0][0]=='S') )
    fin.getline(temp,MAX_CHARACTER); 
  //read ions
  for (int i=0;i<natom;i++){
    fin.getline(temp,MAX_CHARACTER); 
  }
  fin.getline(temp,MAX_CHARACTER); //empty line
  cout<<"read ions"<<endl;

  //read in dimension 
  fin.getline(temp,MAX_CHARACTER); 
  column=0;
  pch = strtok (temp," ");
  while ((pch != NULL)&&(column<MAX_COLUMN)) {
      strcpy(content[column],pch);
      column++;
      pch = strtok (NULL, " ");
  }
  pch= NULL;
  if (column!=3){
    cout<<"FAILED: wrong format"<<endl;
    return 1;
  }
  int nlocal=1;
  int ng[3];
  for (int i=0;i<column;i++){
     ng[i]=atoi(content[i]);
     nlocal*=ng[i];
  }
  
  //read density
  double * vlocal=new double[nlocal];
  for (int i=0;i<nlocal;i++) fin >>vlocal[i];
  cout<<"read density"<<endl;

  double *vav = new double [ng[ider]];
  fill(vav, vav+ng[ider], 0);

  double scale=1./float(nlocal/ng[ider]);
  if (ider == 0) {
    for (int ix=0; ix<ng[0]; ix++){
      for (int iz=0; iz<ng[2]; iz++){
        for (int iy=0; iy<ng[1]; iy++){
          int ipl=ix+(iy+iz*ng[1])*ng[0];
          vav[ix] += vlocal[ipl]*scale;
        }
      }
    }
  }else if ( ider == 1){
    for (int iy=0; iy<ng[1]; iy++){
      for (int iz=0; iz<ng[2]; iz++){
        for (int ix=0; ix<ng[0]; ix++){
          int ipl=ix+(iy+iz*ng[1])*ng[0];
          vav[iy] += vlocal[ipl]*scale;
        }
      }
    }
  }else if (ider==2){
    for (int iz=0; iz<ng[2]; iz++){
      for (int iy=0; iy<ng[1]; iy++){
        for (int ix=0; ix<ng[0]; ix++){
          int ipl=ix+(iy+iz*ng[1])*ng[0];
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
