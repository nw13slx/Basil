
//author: Lixin Sun nw13mifaso@gmail.com

#include "definition.h"

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
  int natom=atof(content[5].c_str());
  double * x= new double [3*natom];
  int *QM = new int [natom];
  int nQM=0;
  int *MM = new int [natom];
  int nMM=0;
  for (int i=0;i<natom;i++){
    fin.getline(temp,MAX_CHARACTER); 
    fout<<temp<<endl;
    parse(temp,column,content);
    x[i*3]=atof(content[1].c_str());
    x[i*3+1]=atof(content[2].c_str());
    x[i*3+2]=atof(content[3].c_str());
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
        sprintf(temp,"%d %d",QM[i]+1,MM[j]+1);
        conn[nconn]=temp;
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
