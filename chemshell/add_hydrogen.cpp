//Usage: input outputname qm_region bd_sym h_sym true/fasle
//  input format : xyz
//author: Lixin Sun nw13mifaso@gmail.com
//
//it find all the possible atom labeled as "bd_sym" that is connected to the qm region
//  and assign them as part of the qm region
//
//add h to all the atoms with region number "qm_region" and label "bd_sym"
//

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
  if (argc < 5){
    cout<<" FAILED: need more input arguments"<<endl;
    cout<<" fin ffout qm_region bd_sym h_sym" <<endl;
    return 1;
  }
  
  ifstream fin(argv[1]); 
  char output[100];
  strcpy(output,argv[2]);
  //qm region number
  string qm(argv[3]);
  //boundary atoms symbol
  string bd_sym(argv[4]);
  //bond type to count
  string h_sym(argv[5]);
  bool reduce=atoi(argv[6]);
  double cutoff=3;


  char temp[MAX_CHARACTER];
  string content[MAX_COLUMN];
  int column, line=0;

  fin.getline(temp,MAX_CHARACTER); 
  parse(temp,column,content);
  int natom=atoi(content[0].c_str());

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
    x[i*3]=  atof(content[1].c_str());
    x[i*3+1]=atof(content[2].c_str());
    x[i*3+2]=atof(content[3].c_str());
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
      q[i]=atof(content[4].c_str());
    }
  }
  fin.close();
  cout<<"initial QM atoms: "<<nQM<<endl;
  cout<<"initial MM atoms: "<<nMM<<endl;

  if (reduce){
    int *new_QM = new int [natom];
    int new_nQM=0;
    cout<<"add all the surrounding oxygen"<<endl;
    for (int i=0; i<nQM; i++){
      if (sym[QM[i]].find(bd_sym)==std::string::npos){
        double *xQM=&x[QM[i]*3];
        for (int j=0; j<nMM; j++){
          if (isQM[MM[j]]==false && sym[MM[j]].find(bd_sym)!=std::string::npos){
            double *xMM=&x[MM[j]*3];
            double dx[3];
            dx[0]=xQM[0]-xMM[0];
            dx[1]=xQM[1]-xMM[1];
            dx[2]=xQM[2]-xMM[2];
            double r=sqrt((dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]));
            if (r<cutoff){
              cout<<new_nQM<<" find one more "<<bd_sym<<" "<<MM[j]<<" bond length: "<<r<<endl;
              new_QM[new_nQM]=MM[j];
              new_nQM++;
              isQM[MM[j]]=true;
            }
          }
        }
      }
    }
    for (int i=0;i<new_nQM;i++){
      QM[nQM]=new_QM[i];
      nQM++;
    }
    delete [] new_QM;
    cout<<"add "<<new_nQM<<" more "<<bd_sym<<" atoms to QM"<<endl;
  }

  int nh=0;
  double *h=new double[nQM*9];
  for (int i=0; i<nQM; i++){
    if (sym[QM[i]].find(bd_sym)!=std::string::npos){
      double *xQM=&x[QM[i]*3];
      int ncoord=0;
      for (int j=0; j<nMM; j++){
        if (isQM[MM[j]]==false && sym[MM[j]].find(h_sym)!=std::string::npos){
          double *xMM=&x[MM[j]*3];
          double dx[3];
          dx[0]=xMM[0]-xQM[0];
          dx[1]=xMM[1]-xQM[1];
          dx[2]=xMM[2]-xQM[2];
          double r=sqrt((dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]));
          if (r<cutoff){
            cout<<"add a hydrogen between "<<i<<" and "<<j<<" bond length: "<<r<<endl;
            ncoord++;
            double dr=1.*r/1.91566;
            h[nh*3]=dx[0]/r*dr+xQM[0];
            h[nh*3+1]=dx[1]/r*dr+xQM[1];
            h[nh*3+2]=dx[2]/r*dr+xQM[2];
            nh++;
            if (reduce)
              q[MM[j]]-=2.23/6.;
          }
        }
      }
    }
  }

  char name[100];
  sprintf(name,"%s.qmxyz",output);
  ofstream fout(name);
  fout<<nh+nQM<<endl<<endl;
  for (int i=0;i<nQM;i++) {
    fout<<sym[QM[i]]<<" "<<x[QM[i]*3]<<" "<<x[QM[i]*3+1]<<" "<<x[QM[i]*3+2]<<endl;
  }
  cout<<"add "<<nh<<" hydrogen atoms in total"<<endl;
  for (int i=0;i<nh;i++){
    fout<<"H "<<h[i*3]<<" "<<h[i*3+1]<<" "<<h[i*3+2]<<endl;
  }
  fout.close();

  sprintf(name,"%s.mmxyz",output);
  fout.open(name);
  fout<<natom-nQM<<endl<<endl;
  for (int i=0;i<nMM;i++) {
    if (isQM[MM[i]]==false)
      fout<<sym[MM[i]]<<" "<<q[MM[i]]<<" "<<x[MM[i]*3]<<" "<<x[MM[i]*3+1]<<" "<<x[MM[i]*3+2]<<endl;
  }
  fout.close();

  sprintf(name,"%s-overall.xyz",output);
  fout.open(name);
  fout<<nh+natom<<endl<<endl;
  for (int i=0;i<nQM;i++) {
    fout<<sym[QM[i]]<<" "<<x[QM[i]*3]<<" "<<x[QM[i]*3+1]<<" "<<x[QM[i]*3+2]<<endl;
  }
  cout<<"add "<<nh<<" hydrogen atoms in total"<<endl;
  for (int i=0;i<nh;i++){
    fout<<"H "<<h[i*3]<<" "<<h[i*3+1]<<" "<<h[i*3+2]<<endl;
  }
  for (int i=0;i<nMM;i++) {
    if (isQM[MM[i]]==false)
      fout<<sym[MM[i]]<<" "<<x[MM[i]*3]<<" "<<x[MM[i]*3+1]<<" "<<x[MM[i]*3+2]<<endl;
  }
  fout.close();
}
