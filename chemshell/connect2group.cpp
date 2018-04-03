
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

#include "stdlib.h"       // for random
#include <vector>

#define MAX_CHARACTER 1000
#define MAX_COLUMN 200

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

int main(int argc, char **argv){
  double cutoff=3;
  if (argc < 2){
    cout<<" FAILED: need more input arguments"<<endl;
    return 1;
  }
  ifstream fin(argv[1]); 
  ofstream fout(argv[2]); 
  string sym1(argv[3]);
  string sym2(argv[4]);

  char temp[MAX_CHARACTER];
  string content[MAX_COLUMN];
  int column, line=0;

  /*block = fragment records = 0
   * block = title records = 1
   * molecule 1
   * block = coordinates records = 3627*/

  fin.getline(temp,MAX_CHARACTER); 
  fout<<temp<<endl;
  fin.getline(temp,MAX_CHARACTER); 
  fout<<temp<<endl;
  fin.getline(temp,MAX_CHARACTER); 
  fout<<temp<<endl;
  fin.getline(temp,MAX_CHARACTER); 
  fout<<temp<<endl;
  parse(temp,column,content);
  int natom=stof(content[5]);
  double * x= new double [3*natom];
  int *QM = new int [natom];
  int nQM=0;
  int *MM = new int [natom];
  int nMM=0;
  for (int i=0;i<natom;i++){
    fin.getline(temp,MAX_CHARACTER); 
    fout<<temp<<endl;
    parse(temp,column,content);
    x[i*3]=stof(content[1]);
    x[i*3+1]=stof(content[2]);
    x[i*3+2]=stof(content[3]);
    if (content[0].compare(sym1)==0){
      QM[nQM]=i;
      nQM++;
    }else if (content[0].compare(sym2)==0){
      MM[nMM]=i;
      nMM++;
    }
  }
  cout<<nQM<<" "<<nMM<<endl;

  fin.getline(temp,MAX_CHARACTER); 
  parse(temp,column,content);
  while (fin.good() and content[2]!="connectivity"){
    fout<<temp<<endl;
    fin.getline(temp,MAX_CHARACTER); 
    parse(temp,column,content);
  }

  int nconn=0;
  string *conn=new string[nQM*10];
  for (int i=0; i<nQM; i++){
    double *xQM=&x[QM[i]*3];
    int ncoord=0;
    for (int j=0; j<nMM; j++){
      double *xMM=&x[MM[j]*3];
      double dx[3];
      dx[0]=xQM[0]-xMM[0];
      dx[1]=xQM[1]-xMM[1];
      dx[2]=xQM[2]-xMM[2];
      double r=sqrt((dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]));
      if (r<cutoff/0.52917721092){
        ncoord++;
        conn[nconn]=to_string(QM[i]+1)+" "+to_string(MM[j]+1);
        nconn++;
      }
    }
    cout<<QM[i]<<" "<<ncoord<<endl;
  }

  fout<<"block = connectivity records = "<<nconn<<endl;
  for (int i=0;i<nconn;i++){
    fout<<conn[i]<<endl;
  }
  fin.close();
  fout.close();

}
