//input format, pun file  for chemshell
//usage: add_shell input output coresymbol shellq
//default: cutoff=3 angstrom
//author: Lixin Sun nw13mifaso@gmail.com

#include "definition.h"

int main(int argc, char **argv){
  if (argc < 5){
    cout<<" FAILED: need more input arguments"<<endl;
    cout<<"usage: add_shell input output core_symbol shellq"<<endl;
    return 1;
  }
  ifstream fin(argv[1]); 
  ofstream fout(argv[2]);
  string sym=argv[3];
  double shellq=atof(argv[4]);

  char temp[MAX_CHARACTER];
  string content[MAX_COLUMN];
  int column, line=0;

  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  parse(temp,column,content);
  int natom = atof(content[5].c_str());
  double * x = new double [3*natom];
  double * q = new double [natom];
  bool *shell = new bool [natom];
  string * type = new string[natom];
  int nshell = 0;
  for (int i=0;i<natom;i++){
    fin.getline(temp,MAX_CHARACTER); 
    parse(temp,column,content);
    x[i*3]=atof(content[1].c_str());
    x[i*3+1]=atof(content[2].c_str());
    x[i*3+2]=atof(content[3].c_str());
    type[i]=content[0];
    if ((content[0].find(sym)!=string::npos) && (content[0].find("3")!=string::npos)) {
      shell[i]=true;
      nshell++;
    } else {
      shell[i]=false;
    } 
  }
  fin.getline(temp,MAX_CHARACTER); 
  for (int i=0;i<natom;i++){
    fin.getline(temp,MAX_CHARACTER); 
    parse(temp,column,content);
    q[i]=atof(content[0].c_str());
  }

  fout<<"block = fragment records = 0"<<endl;
  fout<<"block = title records = 1"<<endl;
  fout<<"molecule 1"<<endl;
  fout<<"block = coordinates records = "<<natom<<endl;
  for (int i=0;i<natom;i++){
    fout<<type[i]<<" "<<std::setprecision(14)<<x[i*3]<<" "<<x[i*3+1]<<" "<<x[i*3+2]<<endl;
  }
  fout<<"block = atom_charges records = "<<natom<<endl;
  for (int i=0;i<natom;i++){
    if (shell[i]){
      fout<<q[i]-shellq<<endl;
    }else{
      fout<<q[i]<<endl;
    }
  }
  fout<<"block = shells records = "<<nshell<<endl;
  for (int i=0;i<natom;i++){
    if (shell[i]){
      fout<<type[i]<<" "<<std::setprecision(14)<<x[i*3]<<" "<<x[i*3+1]<<" "<<x[i*3+2];
      fout<<" "<<std::setprecision(15)<<shellq<<" "<<i<<endl;
    }
  }
  fout.close();
  delete [] x;
  delete [] q;
  delete [] shell;
  delete [] type;

}
