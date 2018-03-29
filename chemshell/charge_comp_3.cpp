
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
  if (argc < 2){
    cout<<" FAILED: need more input arguments"<<endl;
    return 1;
  }
  ifstream fin(argv[1]); 
  ofstream fout(argv[2]); 
  double excess_q=atof(argv[3]);

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
  int nMM=0;
  int *MM=new int [natom];
  for (int i=0;i<natom;i++){
    fin.getline(temp,MAX_CHARACTER); 
    fout<<temp<<endl;
    parse(temp,column,content);
    char chk=content[0].back();
    if (chk=='3'){
      MM[i]=1;
      nMM++;
    }else{
      MM[i]=0;
    }
  }
  double delta=-excess_q/double(nMM);
  cout<<nMM<<" delta: "<<delta<<endl;

  fin.getline(temp,MAX_CHARACTER); 
  fout<<temp<<endl;
  fout<< std::setprecision(10);
  for (int i=0;i<natom;i++){
    double q;
    fin >>q;
    if (MM[i]==1){
      fout<<q+delta<<endl;
    }else{
      fout<<q<<endl;
    }
  }


  fin.getline(temp,MAX_CHARACTER); 
  while (fin.good() ){
    fout<<temp<<endl;
    fin.getline(temp,MAX_CHARACTER); 
  }

  fin.close();
  fout.close();

}
