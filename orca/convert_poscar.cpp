// this code is used to read the orca output
// and convert the atomic configuration to pdb format
// Q+ will be label as Qp and Q- will be label as Qn
// embedding potential is labeled as XXpo
// 
// usage: g++ convert_pdb.cpp -o orca_pdb
//        orca_pdb <inputfile name> <outputfile name> 
// author: Lixin Sun nw13mifaso@gmail.com

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

//set up the size of array used for storage
#define MAX_ELEMENT 10
#define MAX_COORD 1000
#define MAX_ENERGYLINE 1000

//set up the buffer size for reading
#define MAX_CHARACTER 1000
#define MAX_COLUMN 20

bool contain(const char *string, char c){
    int n=strlen(string);
    for (int i=0;i<n;i++){
        if (string[i]==c){
            return true;
        }
    }
    return false;
}

void gaussian(int grid, double *x, double *y, double x0, double sigma2){
    for (int i=0;i<grid;i++) y[i]=exp(-(x[i]-x0)*(x[i]-x0)/2./sigma2);
}

int main(int argc, char **argv){

    char temp[MAX_CHARACTER], * pch,content[MAX_COLUMN][MAX_CHARACTER];
    int column, line=0;
    bool Print_Conf=false;
    
    //locate the beginning of the print_MO
    ifstream In1(argv[1]);
    while (In1.good()){
        line++;
        In1.getline(temp,MAX_CHARACTER); 

        column=0;
        pch = strtok (temp," ");
        while ((pch != NULL)&&(column<MAX_COLUMN)) {
            strcpy(content[column],pch);
            column++;
            pch = strtok (NULL, " ");
        }
        pch= NULL;
        if (column == 3 ){
            if (( strcmp(content[0],"CARTESIAN") == 0 ) && (strcmp(content[1],"COORDINATES") == 0) && (strcmp(content[2],"(A.U.)")==0)){
                Print_Conf=true;
                break;
            }
        }
    }
    if ((! In1.good()) || (Print_Conf==false)){
        cout<<"ERROR: CARTESIAN COORDINATES (ANGSTROEM) block not found"<<endl;
        return 1;
    }
    In1.getline(temp,MAX_CHARACTER); //another buffering line
    In1.getline(temp,MAX_CHARACTER); //another buffering line
    line++;

    char element[MAX_ELEMENT][10];
    double energy[MAX_ENERGYLINE];
    int n_element=0;

    double *x=new double [MAX_ELEMENT*MAX_COORD*3];
    int *id=new int [MAX_ELEMENT*MAX_COORD];
    int *natom=new int [MAX_ELEMENT];
    for (int i=0;i<MAX_ELEMENT;i++) natom[i]=0;
    
    double boundary[6]={1000,0,1000,0,1000,0};
    int count=0;
    do{
      line++;
      In1.getline(temp,MAX_CHARACTER); 
      column=0;
      pch = strtok (temp," ");
      while (pch != NULL) {
          strcpy(content[column],pch);
          column++;
          pch = strtok (NULL, " ");
      }
      pch= NULL;

      if (column>7 && content[0][0]!='*' && content[0][0]!='>'){
        char ele[10];
        //second column is the label
        if (content[1][0]=='Q'){
            if (atof(content[2])>0) sprintf(ele," Qp ");
            else sprintf(ele," Qn ");
        } else if (contain(content[1],'>')==true){
            int idxToDel=strlen(content[1])-1;
            memmove(&content[1][idxToDel], &content[1][idxToDel + 1], strlen(content[1]) - idxToDel);
            sprintf(ele,"%2spo",content[1]);
        }else sprintf(ele,"%2s  ",content[1]);

        //recognize the element name
        int elementid=-1;
        for (int eid=0;eid<n_element;eid++){
            if (strcmp(ele,element[eid])==0) elementid=eid;
        }
        if (elementid==-1){
            elementid=n_element;
            strcpy(element[n_element],ele);
            n_element++;
        }

        double *xx=&x[elementid*MAX_COORD*3+natom[elementid]*3];
        count++;
        id[elementid*MAX_COORD+natom[elementid]]=count;

        xx[0]=atof(content[5])*0.52917721092;
        xx[1]=atof(content[6])*0.52917721092;
        xx[2]=atof(content[7])*0.52917721092;
        natom[elementid]++;

        if (xx[0]<boundary[0]) boundary[0]=xx[0];
        if (xx[1]<boundary[2]) boundary[2]=xx[1];
        if (xx[2]<boundary[4]) boundary[4]=xx[2];
        if (xx[0]>boundary[1]) boundary[1]=xx[0];
        if (xx[1]>boundary[3]) boundary[3]=xx[1];
        if (xx[2]>boundary[5]) boundary[5]=xx[2];
        
      }
    }while (column!=0);
    cout<<"done reading"<<endl;

    ofstream out(argv[2]);
    for (int j=0;j<n_element;j++){
      double *xx=&x[j*MAX_COORD*3];
      int *iid=&id[j*MAX_COORD];
      for (int k=0;k<natom[j];k++){
        double *xxx=&xx[k*3];
        int iiid=iid[k];
        xxx[0]-=(boundary[0]-5);
        xxx[1]-=(boundary[2]-5);
        xxx[2]-=(boundary[4]-5);
        sprintf(temp,"ATOM  %5d %4s              %8.3f%8.3f%8.3f%6.2f%6.2f          %2s ",iiid,element[j],xxx[0],xxx[1],xxx[2],1.0,1.0,"O");
        out<<temp<<endl;
      }
    }
    out<<"END"<<endl;
}
