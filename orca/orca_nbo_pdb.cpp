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

#include "functions.h"
#include "math.h"
#include "stdlib.h"       // for random
#include <vector>


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
    int n_element=0;

    double *x=new double [MAX_ELEMENT*MAX_COORD*3];
    int *id=new int [MAX_ELEMENT*MAX_COORD];
    int *natom=new int [MAX_ELEMENT];
    double *qtype=new double [MAX_ELEMENT*MAX_COORD];
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
        double qq=-50;
        if (content[1][0]=='Q'){
            qq=atof(content[2].c_str());
            if (qq>0) ele=" q+ ";
            else ele=" q- ";
        } else if ( content[1].find_first_of(">")!=string::npos){
            qq=atof(content[2].c_str());
            content[1].erase(content[1].length()-1,1);
            ele=content[1]+"qp";
        }else {
          ele=content[1]+"  ";
        }
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
        qtype[natom_tot]=qq;
        natom_tot++;
        id[elementid*MAX_COORD+natom[elementid]]=natom_tot;

        xx[0]=atof(content[5].c_str())*0.52917721092;
        xx[1]=atof(content[6].c_str())*0.52917721092;
        xx[2]=atof(content[7].c_str())*0.52917721092;
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
    string pattern="Summary of Natural";
    int beg_pos=In2.tellg();
    int v_pos=find_pattern(In2,pattern,"contains",true);
    if (v_pos<0){
      cout<<"ERROR: Summary of Natural Population Analysis block not found"<<endl;
      return 1;
    }else{
      Print_Conf=true;
    }
    In2.getline(temp,MAX_CHARACTER);
    In2.getline(temp,MAX_CHARACTER);
    In2.getline(temp,MAX_CHARACTER);
    In2.getline(temp,MAX_CHARACTER);
    In2.getline(temp,MAX_CHARACTER);
    In2.getline(temp,MAX_CHARACTER);
    double *q=new double[natom_tot];
    for (int i=0;i<natom_tot;i++){
      if (qtype[i]>-50){
        q[i]=qtype[i];
      }else{
        In2.getline(temp,MAX_CHARACTER);
        column=break_line(temp,content);
        if (column==8){
          q[i]=atof(content[2].c_str());
        }else if (column==7){
          q[i]=atof(content[1].c_str());
        }else{
          cout<<"???"<<column<<endl;
          q[i]=atof(content[1].c_str());
        }
      }
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
        //xxx[0]-=(boundary[0]-5);
        //xxx[1]-=(boundary[2]-5);
        //xxx[2]-=(boundary[4]-5);
        sprintf(temp,"ATOM  %5d %4s              %8.3f%8.3f%8.3f%6.2f%6.2f          %2s ",iiid,element[j].c_str(),xxx[0],xxx[1],xxx[2],1.0,qq,"O");
        out<<temp<<endl;
      }
    }
    out<<"END"<<endl;
    out.close();

   if (argc>5){
     int *coord=new int [MAX_ELEMENT*MAX_COORD];
     for (int i=0;i<(MAX_ELEMENT*MAX_COORD); i++){
       coord[i]=0;
     }

     ofstream out2(argv[4]);
     double cutoff=atof(argv[5]);
     for (int j=0;j<n_element;j++){
       double *xx=&x[j*MAX_COORD*3];
       int *iid=&id[j*MAX_COORD];
       int *cc=&coord[j*MAX_COORD];
       for (int k=0;k<natom[j];k++){
         double *xxx=&xx[k*3];
         int iiid=id[k];
         int *icc=&cc[k];
         for (int j1=j+1;j1<n_element;j1++){
           double *xx1=&x[j1*MAX_COORD*3];
           int *iid1=&id[j1*MAX_COORD];
           int *cc1=&coord[j1*MAX_COORD];
           for (int k1=0;k1<natom[j1];k1++){
             double *xxx1=&xx1[k1*3];
             int *icc1=&cc1[k1];
             int iiid1=iid1[k1];
             double dr=(xxx1[0]-xxx[0])*(xxx1[0]-xxx[0]);
             dr += (xxx1[1]-xxx[1])*(xxx1[1]-xxx[1]);
             dr += (xxx1[2]-xxx[2])*(xxx1[2]-xxx[2]);
             dr = sqrt(dr);
             if ( dr < cutoff ) {
               icc1[0]+=1;
               icc[0]+=1;
             }
           }
         }
       }
     }
      for (int j=0;j<n_element;j++){
        double *xx=&x[j*MAX_COORD*3];
        int *iid=&id[j*MAX_COORD];
        int *cc=&coord[j*MAX_COORD];
        for (int k=0;k<natom[j];k++){
          double *xxx=&xx[k*3];
          int iiid=iid[k];
          double qq=q[iiid-1];
          int icc=cc[k];
          out2 << iiid << " " <<element[j]<<" " <<qq<<" "<<icc<<endl;
       }
     }
     out2.close();

      
   }


}
