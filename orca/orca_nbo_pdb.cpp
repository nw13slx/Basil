// this code is used to read the orca output
// and convert the atomic configuration to pdb format
// Q+ will be label as Qp and Q- will be label as Qn
// embedding potential is labeled as XXpo
// the partial charge is put at the b-factor column
// 
// usage: g++ convert_pdb.cpp -o orca_nbo_pdb
//        orca_nbo_pdb <orcaoutput> <nbooutput> <output filename>
//        orca_nbo_pdb orcaout.asdfa FILE.nbo output.pdb
// author: Lixin Sun nw13mifaso@gmail.com
//
// note: the results can be visualize by nbo.tcl

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
//#include <stdlib.h>
#include <string>
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
#define MAX_COLUMN 50


void gaussian(int grid, double *x, double *y, double x0, double sigma2){
    for (int i=0;i<grid;i++) y[i]=exp(-(x[i]-x0)*(x[i]-x0)/2./sigma2);
}


int split(char *temp, string * content){
   int column=0;
   stringstream ss(temp);
   while ((ss>>content[column])&&(column<MAX_COLUMN)) {
       column++;
   }
   return column;
}

int main(int argc, char **argv){

    char temp[MAX_CHARACTER], * pch;
    string content[MAX_COLUMN];
    int column, line=0;
    bool Print_Conf=false;
    
    //locate the beginning of the print_MO
    ifstream In1(argv[1]);
    while (In1.good()){
        line++;
        In1.getline(temp,MAX_CHARACTER); 

        column=split(temp,content);;
        if (column == 3 ){
            if (( content[0]=="CARTESIAN" ) && (content[1]=="COORDINATES") && (content[2]=="(A.U.)")){
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

    string element[10];
    double energy[MAX_ENERGYLINE];
    int n_element=0;

    double *x=new double [MAX_ELEMENT*MAX_COORD*3];
    int *id=new int [MAX_ELEMENT*MAX_COORD];
    int *natom=new int [MAX_ELEMENT];
    for (int i=0;i<MAX_ELEMENT;i++) natom[i]=0;
    
    double boundary[6]={1000,0,1000,0,1000,0};
    int natom_tot=0;
    do{
      line++;
      In1.getline(temp,MAX_CHARACTER); 
      column=split(temp,content);

      if (column>7 && content[0][0]!='*' && content[0][0]!='>'){
        string ele;
        //second column is the label
        if (content[1][0]=='Q'){
            if (stod(content[2])>0) ele=" Qp ";
            else ele=" Qn ";
        } else if ( content[1].find_first_of(">")!=string::npos){
            content[1].erase(content[1].length()-1,1);
            ele=content[1]+"po";
        }else ele=content[1]+"  ";
        while (ele.length()<4)
           ele=" "+ele;
        if (ele.length()>4)
            cout<<ele<<endl;

        //recognize the element name
        int elementid=-1;
        for (int eid=0;eid<n_element;eid++){
            if (ele==element[eid]) elementid=eid;
        }
        if (elementid==-1){
            elementid=n_element;
            element[n_element]=ele;
            n_element++;
        }

        double *xx=&x[elementid*MAX_COORD*3+natom[elementid]*3];
        natom_tot++;
        id[elementid*MAX_COORD+natom[elementid]]=natom_tot;

        xx[0]=stod(content[5])*0.52917721092;
        xx[1]=stod(content[6])*0.52917721092;
        xx[2]=stod(content[7])*0.52917721092;
        natom[elementid]++;

        if (xx[0]<boundary[0]) boundary[0]=xx[0];
        if (xx[1]<boundary[2]) boundary[2]=xx[1];
        if (xx[2]<boundary[4]) boundary[4]=xx[2];
        if (xx[0]>boundary[1]) boundary[1]=xx[0];
        if (xx[1]>boundary[3]) boundary[3]=xx[1];
        if (xx[2]>boundary[5]) boundary[5]=xx[2];
        
      }
    }while (column!=0);
    cout<<"done reading poscar"<<endl;
    In1.close();

    ifstream In2(argv[2]);
    while (In2.good()){
        line++;
        In2.getline(temp,MAX_CHARACTER); 
        column=split(temp,content);;
        if (column == 5 ){
            if (( content[0]=="Summary" ) && (content[1]=="of") && (content[2]=="Natural")){
                Print_Conf=true;
                break;
            }
        }
    }
    if ((! In2.good()) || (Print_Conf==false)){
        cout<<"ERROR: Summary of Natural Population Analysis block not found"<<endl;
        return 1;
    }
    In2.getline(temp,MAX_CHARACTER);
    In2.getline(temp,MAX_CHARACTER);
    In2.getline(temp,MAX_CHARACTER);
    In2.getline(temp,MAX_CHARACTER);
    In2.getline(temp,MAX_CHARACTER);
    double *q=new double[natom_tot];
    for (int i=0;i<natom_tot;i++){
        In2 >> temp>>temp>>q[i];
        In2.getline(temp,MAX_CHARACTER);
    }
/*
 Summary of Natural Population Analysis:

                                     Natural Population                 Natural
             Natural    ---------------------------------------------    Spin
  Atom No    Charge        Core      Valence    Rydberg      Total      Density
 -------------------------------------------------------------------------------
    O  1   -1.54331      1.99999     7.53196    0.01135     9.54331     0.00000
    */

    ofstream out(argv[3]);
    for (int j=0;j<n_element;j++){
      double *xx=&x[j*MAX_COORD*3];
      int *iid=&id[j*MAX_COORD];
      for (int k=0;k<natom[j];k++){
        double *xxx=&xx[k*3];
        int iiid=iid[k];
        double qq=q[iiid-1];
        xxx[0]-=(boundary[0]-5);
        xxx[1]-=(boundary[2]-5);
        xxx[2]-=(boundary[4]-5);
        sprintf(temp,"ATOM  %5d %4s              %8.3f%8.3f%8.3f%6.2f%6.2f          %2s ",iiid,element[j].c_str(),xxx[0],xxx[1],xxx[2],1.0,qq,"O");
        out<<temp<<endl;
      }
    }
    out<<"END"<<endl;
}
