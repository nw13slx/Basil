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
#define MAX_COORD 1000
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
    string pattern="CARTESIAN COORDINATES (A.U.)";
    double *x=new double [MAX_ELEMENT*MAX_COORD*3];
    double *qecp=new double [MAX_ELEMENT*MAX_COORD];
    bool   *isecp = new bool [MAX_ELEMENT*MAX_COORD];
    int natom=0;
    int v_pos=find_pattern(In1,pattern);
    if (v_pos!=-1){
      In1.getline(temp,MAX_CHARACTER); //another buffering line
      In1.getline(temp,MAX_CHARACTER); //another buffering line
      line++;

      string element[10];
      double energy[MAX_ENERGYLINE];
      int n_element=0;
      do{
        line++;
        In1.getline(temp,MAX_CHARACTER); 
        column=split(temp,content);

        if (column>7 && content[0][0]!='*' && content[0][0]!='>'){
          string ele;
          //second column is the label
          if ((content[1][0]!='Q') && ( content[1].find_first_of(">")==string::npos)){
            double *xx=&x[natom*3];
            xx[0]=atof(content[5].c_str());
            xx[1]=atof(content[6].c_str());
            xx[2]=atof(content[7].c_str());
            isecp[natom]=false;
            natom++;
          }else if ( content[1].find_first_of(">")!=string::npos) {
            double *xx=&x[natom*3];
            qecp[natom]=atof(content[2].c_str());
            isecp[natom]=true;
            xx[0]=atof(content[5].c_str());
            xx[1]=atof(content[6].c_str());
            xx[2]=atof(content[7].c_str());
            natom++;
          }
        }
      }while (column!=0);
    } else{
      cout<<"ERROR: CARTESIAN COORDINATES (ANGSTROEM) block not found"<<endl;
      return 1;
    }
    In1.close();
    double *q=new double[natom];

    ifstream In2(argv[2]);
    pattern="Summary of Natural Population Analysis";
    v_pos=find_pattern(In2,pattern);
    if (v_pos!=-1){
      In2.getline(temp,MAX_CHARACTER);
      In2.getline(temp,MAX_CHARACTER);
      In2.getline(temp,MAX_CHARACTER);
      In2.getline(temp,MAX_CHARACTER);
      In2.getline(temp,MAX_CHARACTER);
      for (int i=0;i<natom;i++){
        if (isecp[i]==true){
          q[i]=qecp[i];
        } else{
          In2.getline(temp,MAX_CHARACTER); 
          column=split(temp,content);
          if (column >=7){
            q[i]=atof(content[2].c_str());
          }else{
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
    } else{
        cout<<"ERROR: Summary of Natural Population Analysis block not found"<<endl;
        return 1;
    }

    ofstream out(argv[3]);
    out<<"block = fragment records = 0"<<endl;
    out<<"block = title records = 1"<<endl;
    out<<"molecule 3"<<endl;
    out<<"block = coordinates records = "<<natom<<endl;
    for (int k=0;k<natom;k++){
      double *xx=&x[k*3];
      out<<"c "<<xx[0]<<" "<<xx[1]<<" "<<xx[2]<<endl;
    }
    out<<"block = atom_charges records = "<<natom<<endl;
    for (int k=0;k<natom;k++){
      double *xx=&x[k*3];
      out<<q[k]<<endl;
    }
//sprintf(temp,"ATOM  %5d %4s              %8.3f%8.3f%8.3f%6.2f%6.2f          %2s ",iiid,element[j].c_str(),xxx[0],xxx[1],xxx[2],1.0,qq,"O");
    out<<"block = connectivity records = 0"<<endl;
    delete [] x;
    delete [] q;
    delete [] qecp;
    delete [] isecp;
    return 0;
}
