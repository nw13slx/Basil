//average the local potential file along certain direction
//usage: average_cube input output direction<1-3>

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
    cout<<" exec input output direction" <<endl;
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

  int natom;
  fin >> natom;
  fin.getline(temp,MAX_CHARACTER); 
  cout << "number of atoms " <<natom<<endl;;

  int nlocal, ng[3];
  fin >> ng[0];
  fin.getline(temp,MAX_CHARACTER); 
  fin >> ng[1];
  fin.getline(temp,MAX_CHARACTER); 
  fin >> ng[2];
  fin.getline(temp,MAX_CHARACTER); 
  cout << "number of grids " << ng[0]<< " "<<ng[1]<<" "<<ng[2]<<endl;
  nlocal = ng[0]*ng[1]*ng[2];
  cout << "number of values " << nlocal<<endl;

  cout<< "read ions"<<endl;
  for (int i=0;i<natom;i++){
    fin.getline(temp,MAX_CHARACTER); 
  }
  
  cout<< "read density"<<endl;
  double * vlocal=new double[nlocal];
  for (int i=0;i<nlocal;i++) {
      fin >>vlocal[i];
  }
  cout<<"finished reading density"<<endl;

  double *vav = new double [ng[ider]];
  fill(vav, vav+ng[ider], 0);

  cout<< "sum across direction "<<ider<<endl;
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
