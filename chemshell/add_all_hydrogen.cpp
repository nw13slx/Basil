
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
#include <malloc.h>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <sstream>
using namespace std;

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

// only take xyz format

int main(int argc, char **argv){
  if (argc < 3){
    cout<<" FAILED: need more input arguments"<<endl;
    return 1;
  }
  
  ifstream fin(argv[1]); 
  ofstream fout(argv[2]); 
  //qm region number
  string qm(argv[3]);
  double cutoff=2;

  char temp[MAX_CHARACTER];
  string content[MAX_COLUMN];
  int column, line=0;

  fin.getline(temp,MAX_CHARACTER); 
  parse(temp,column,content);
  int natom=stoi(content[0]);

  fin.getline(temp,MAX_CHARACTER); 
  string *sym=new string[natom];
  double * x= new double [3*natom];
  double * q= new double [natom];
  bool * isQM = new bool [natom];
  int *QM = new int [natom];
  int nQM=0;
  int *MM = new int [natom];
  int nMM=0;
  for (int i=0;i<natom;i++){
    fin.getline(temp,MAX_CHARACTER); 
    parse(temp,column,content);
    sym[i]=content[0];
    x[i*3]=stof(content[1]);
    x[i*3+1]=stof(content[2]);
    x[i*3+2]=stof(content[3]);
    if (content[0].find(qm)!=std::string::npos){
      QM[nQM]=i;
      nQM++;
      isQM[i]=true;
    }else{ //if (content[0].compare(sym2)==0){
      MM[nMM]=i;
      nMM++;
      isQM[i]=false;
    }
    if (column>4){
      q[i]=stof(content[4]);
    }
  }
  fin.close();
  cout<<nQM<<" "<<nMM<<endl;

  int nh=0;
  double *h=new double[nQM*9];
  for (int i=0; i<nQM; i++){
    double *xQM=&x[QM[i]*3];
    int ncoord=0;
    for (int j=0; j<nMM; j++){
      double *xMM=&x[MM[j]*3];
      double dx[3];
      dx[0]=xMM[0]-xQM[0];
      dx[1]=xMM[1]-xQM[1];
      dx[2]=xMM[2]-xQM[2];
      double r=sqrt((dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]));
      if (r<cutoff){
        ncoord++;
        double dr=1.*r/1.91566;
        h[nh*3]=dx[0]/r*dr+xQM[0];
        h[nh*3+1]=dx[1]/r*dr+xQM[1];
        h[nh*3+2]=dx[2]/r*dr+xQM[2];
        nh++;
        if (q[MM[j]]>0){
          q[MM[j]]-=2.23/6.;
        }else{
          q[MM[j]]+=2.23/6.;
        }
      }
    }
  }
  
  fout<<nh+nQM<<endl<<endl;
  for (int i=0;i<nQM;i++) {
    fout<<sym[QM[i]]<<" "<<x[QM[i]*3]<<" "<<x[QM[i]*3+1]<<" "<<x[QM[i]*3+2]<<endl;
  }
  for (int i=0;i<nh;i++){
    fout<<"H "<<h[i*3]<<" "<<h[i*3+1]<<" "<<h[i*3+2]<<endl;
  }

  fout<<natom-nQM<<endl<<endl;
  for (int i=0;i<nMM;i++) {
    if (isQM[MM[i]]==false)
      fout<<sym[MM[i]]<<" "<<q[MM[i]]<<" "<<x[MM[i]*3]<<" "<<x[MM[i]*3+1]<<" "<<x[MM[i]*3+2]<<endl;
  }

  fout.close();
}
