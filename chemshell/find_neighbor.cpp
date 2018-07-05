// analyse average bond length between sym1-sym2, for inside, outside and across the QM/MM boundary
// input format, pun file  for chemshell
//usage: bondlength  input symbol1 symbol2 cutoff
//author: Lixin Sun nw13mifaso@gmail.com

#include "definition.h"

int main(int argc, char **argv){
  if (argc < 3){
    cout<<" FAILED: need more input arguments"<<endl;
    cout<<"usage: bondlength  input symbol1 symbol2 cutoff"<<endl;
    return 1;
  }
  ifstream fin(argv[1]); 
  int id0=atoi(argv[2]);
  double cutoff;
  if (argc < 4 ) cutoff=3;
  else cutoff=atof(argv[3]);
  double bondlength=1.5;

  char temp[MAX_CHARACTER];
  string content[MAX_COLUMN];
  int column, line=0;

  /*block = fragment records = 0
   * block = title records = 1
   * molecule 1
   * block = coordinates records = 3627*/


  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  parse(temp,column,content);
  int natom=atof(content[5].c_str());
  double * x= new double [3*natom];
  int *QM = new int [natom];
  int nQM=0;
  int *MM = new int [natom];
  int nMM=0;
  string * type=new string[natom];
  for (int i=0;i<natom;i++){
    fin.getline(temp,MAX_CHARACTER); 
    parse(temp,column,content);
    x[i*3]=atof(content[1].c_str());
    x[i*3+1]=atof(content[2].c_str());
    x[i*3+2]=atof(content[3].c_str());
    type[i]=content[0];
    if (content[0].find("1")!=std::string::npos){
      //cout<<"QM "<<content[0]<<endl;
      QM[nQM]=i;
      nQM++;
    } else if (content[0].find(string("2"))!=std::string::npos){
      //cout<<"MM "<<content[0]<<endl;
      MM[nMM]=i;
      nMM++;
    } else if (content[0].find(string("3"))!=std::string::npos){
      //cout<<"MM "<<content[0]<<endl;
      MM[nMM]=i;
      nMM++;
    }
  }
  cout<<nQM<<" "<<nMM<<endl;
  fin.close();

  double sx[3]={0,0,0};
  int ns=0;

  // QM-QM bond
  double *xQM=&x[QM[id0]*3];
  string sum="";
  for (int i=0; i<nQM; i++){
    string type1=type[QM[i]];
    double *xQM2=&x[QM[i]*3];
    double dx[3];
    dx[0]=xQM2[0]-xQM[0];
    dx[1]=xQM2[1]-xQM[1];
    dx[2]=xQM2[2]-xQM[2];
    double r=sqrt((dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]));
    if (r<cutoff/ANG2BOHR){
      cout<<i<<" "<<type1<<" "<<r*ANG2BOHR<<endl;
      sum=sum+" ";
      sum=sum+to_string(i);
      if (r!=0){
        sx[0]+=dx[0];
        sx[1]+=dx[1];
        sx[2]+=dx[2];
        ns++;
      }
    }
  }
  for (int i=0; i<nMM; i++){
    string type1=type[MM[i]];
    double *xMM=&x[MM[i]*3];
    double dx[3];
    dx[0]=xMM[0]-xQM[0];
    dx[1]=xMM[1]-xQM[1];
    dx[2]=xMM[2]-xQM[2];
    double r=sqrt((dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]));
    if (r<cutoff/ANG2BOHR){
      cout<<i<<" "<<type1<<" "<<r*ANG2BOHR<<endl;
      sum=sum+" ";
      sum=sum+to_string(i);
      if (r!=0){
        sx[0]+=dx[0];
        sx[1]+=dx[1];
        sx[2]+=dx[2];
        ns++;
      }
    }
  }
  sx[0]=sx[0]/double(ns);
  sx[1]=sx[1]/double(ns);
  sx[2]=sx[2]/double(ns);
  double ls=sqrt(sx[0]*sx[0]+sx[1]*sx[1]+sx[2]*sx[2]);
  cout<<"all "<<ns<<" neighbor "<<sum<<endl;
  cout<<"weight center "<<sx[0]<<" "<<sx[1]<<" "<<sx[2]<<endl;
  cout<<"distance "<<ls<<endl;
  cout<<"add atom "<<xQM[0]+sx[0]/ls*bondlength<<" "<<xQM[1]+sx[1]/ls*bondlength<<" "<<xQM[2]+sx[2]/ls*bondlength<<endl;

}
