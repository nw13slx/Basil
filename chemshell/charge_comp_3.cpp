
//author: Lixin Sun nw13mifaso@gmail.com

#include "definition.h"

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
