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
  if (argc < 4){
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
  cout<<"average over ";
  switch (ider) {
      case 0: cout<<"x"<<endl; break;
      case 1: cout<<"y"<<endl; break;
      case 2: cout<<"z"<<endl; break;
  }

  char temp[MAX_CHARACTER], * pch,content[MAX_COLUMN][MAX_CHARACTER];
  int column, line=0;
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 

  int natom;
  fin >> natom;
  fin.getline(temp,MAX_CHARACTER); 
  cout << "number of atoms " <<natom<<endl;;

  int nlocal, ng[3];
  double ls[3][3];
  fin >> ng[0] >> ls[0][0] >> ls[0][1] >> ls[0][2];
  fin.getline(temp,MAX_CHARACTER); 
  fin >> ng[1] >> ls[1][0] >> ls[1][1] >> ls[1][2];
  fin.getline(temp,MAX_CHARACTER); 
  fin >> ng[2] >> ls[2][0] >> ls[2][1] >> ls[2][2];
  fin.getline(temp,MAX_CHARACTER); 
  cout << "number of grids " << ng[0]<< " "<<ng[1]<<" "<<ng[2]<<endl;
  nlocal = ng[0]*ng[1]*ng[2];
  cout << "number of values " << nlocal<<endl;

  cout<< "read ions"<<endl;
  double xx[3];
  double xmin=10000, xmax=-10000;
  for (int i=0;i<natom;i++){
    fin >> temp >> temp >> xx[0] >> xx[1] >> xx[2];
    fin.getline(temp,MAX_CHARACTER); 
    if (xx[ider] > xmax) xmax=xx[ider];
    if (xx[ider] < xmin) xmin=xx[ider];
  }
  
  cout<< "read density"<<endl;
  double vlocal[ng[0]][ng[1]][ng[2]];
  cout<<"finished reading density"<<endl;
  for (int ix=0; ix<ng[0]; ix++)
    for (int iy=0; iy<ng[1]; iy++)
      for (int iz=0; iz<ng[2]; iz++)
          fin >> vlocal[ix][iy][iz];

  double *vav = new double [ng[ider]];
  fill(vav, vav+ng[ider], 0);

  cout<< "sum across direction "<<ider<<endl;
  double scale=1./float(nlocal/ng[ider]);
  int index[3];
  for (index[0]=0; index[0]<ng[0]; index[0]++){
    for (index[1]=0; index[1]<ng[1]; index[1]++){
      for (index[2]=0; index[2]<ng[2]; index[2]++){
        vav[index[ider]] += vlocal[index[0]][index[1]][index[2]]*scale;
      }
    }
  }

  fout<<"# "<<xmin<<" "<<xmax<<endl;
  fout<<"# "<<ls[ider][0]<<" "<<ls[ider][1]<<" "<<ls[ider][2]<<endl;
  fout<<"# "<<ng[ider]<<" "<<ider+1<<endl;
  for (int i=0; i<ng[ider];i++){
    fout<<i<<" "<<vav[i]<<endl;
  }
}
