//input format, pun file  for chemshell
//usage: add_ecp input output ecp_symbol cutoff
//default: cutoff=3 angstrom
//author: Lixin Sun nw13mifaso@gmail.com

#include "definition.h"

int main(int argc, char **argv){
  if (argc < 4){
    cout<<" FAILED: need more input arguments"<<endl;
    cout<<"usage: add_ecp input output ecp_symbol cutoff(default: 3)"<<endl;
    return 1;
  }
  ifstream fin(argv[1]); 
  ofstream fout(argv[2]);
  string sym=argv[3];
  double cutoff;
  if (argc < 5 ) cutoff=3/ANG2BOHR;
  else cutoff=atof(argv[4])/ANG2BOHR;

  char temp[MAX_CHARACTER];
  string content[MAX_COLUMN];
  int column, line=0;

  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  fin.getline(temp,MAX_CHARACTER); 
  parse(temp,column,content);
  int natom=atof(content[5].c_str());
  double * x= new double [3*natom];
  double * q= new double [natom];
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
    if (content[0].find("1")!=string::npos){
      QM[nQM]=i;
      nQM++;
    } else {
      MM[nMM]=i;
      nMM++;
    } 
  }
  fin.getline(temp,MAX_CHARACTER); 
  for (int i=0;i<natom;i++){
    fin.getline(temp,MAX_CHARACTER); 
    parse(temp,column,content);
    q[i]=atof(content[0].c_str());
  }
  fin.getline(temp,MAX_CHARACTER);
  int nshell=0;
  double * shellx, *shellq;
  int * core;
  string * shelltype;
  if (fin.good()){
    parse(temp,column,content);
    nshell=atoi(content[5].c_str());
    shellx=new double [3*nshell];
    shellq=new double [nshell];
    core=new int [nshell];
    shelltype=new string[natom];
    for (int i=0;i<nshell;i++){
      fin.getline(temp,MAX_CHARACTER); 
      parse(temp,column,content);
      shelltype[i]=content[0];
      shellx[i*3]=atof(content[1].c_str());
      shellx[i*3+1]=atof(content[2].c_str());
      shellx[i*3+2]=atof(content[3].c_str());
      shellq[i]=atof(content[4].c_str());
      cout<<std::fixed<<std::setprecision(9) << shellq[i]<<" "<<content[4]<<endl;
      core[i]=atoi(content[5].c_str());
    }
  }
  cout<<"nQM nMM ";
  cout<<nQM<<" "<<nMM<<endl;
  fin.close();

  double sx[3]={0,0,0};
  int ns=0;

  cout<<(sym+"3")<<endl;
  for (int i=0; i<nQM; i++){
    double *xQM=&x[QM[i]*3];
    for (int j=0; j<nMM; j++){
      string type1=type[MM[j]];
      if (type1.find(sym+"3")!=string::npos){
        double *xMM=&x[MM[j]*3];
        double dx[3];
        dx[0]=xMM[0]-xQM[0];
        dx[1]=xMM[1]-xQM[1];
        dx[2]=xMM[2]-xQM[2];
        double r=sqrt((dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]));
        if (r<cutoff){
          int pos=type1.find(string("3"));
          type[MM[j]][pos]='2';
          q[MM[j]]=0;
          cout<<MM[j]<<" "<<type[MM[j]]<<endl;
        }
      }
    }
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
    fout<<q[i]<<endl;
  }
  if (nshell>0){
    fout<<"block = shells records = "<<nshell<<endl;
    for (int i=0;i<nshell;i++){
      fout<<shelltype[i]<<" "<<std::setprecision(14)<<shellx[i*3]<<" "<<shellx[i*3+1]<<" "<<shellx[i*3+2];
      fout<<" "<<std::setprecision(15)<<shellq[i]<<" "<<core[i]<<endl;
    }
  }
  fout.close();
  delete [] x;
  delete [] q;

}
