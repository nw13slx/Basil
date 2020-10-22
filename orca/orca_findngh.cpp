//author: Lixin Sun nw13mifaso@gmail.com

#include "functions.h"

int main(int argc, char **argv){
  if (argc < 3){
    cout<<" FAILED: need more input arguments"<<endl;
    cout<<"usage: orca_findngh input id0 cutoff addatom_bondlength"<<endl;
    return 1;
  }
  ifstream fin(argv[1]);
  int id0=atoi(argv[2]);
  double cutoff;
  if (argc < 4 ) cutoff=3;
  else cutoff=atof(argv[3]);
  double bondlength=1.5;
  if (argc >4 ) bondlength=atof(argv[4]);

  char temp[MAX_CHARACTER];
  string content[MAX_COLUMN];
  int column, line=0;

  string mode="xyz";
  string pattern="%coords";
  int beg_pos=fin.tellg();
  int v_pos=find_pattern(fin,pattern,"contains",true);
  if (v_pos==-1){
    mode="none";
    fin.seekg(beg_pos);
  }else{
    mode="coords";
  }

  int natom=0;
  double *x=new double [3*MAX_ATOMS];
  string *type=new string[MAX_ATOMS];
  bool bohr=false;

  if (mode=="coords"){
    fin.getline(temp,MAX_CHARACTER); 
    beg_pos=fin.tellg();
    v_pos=find_pattern(fin,string("coords"),"contains",false);
    fin.getline(temp,MAX_CHARACTER); 
    parse(temp,column,content);
    while (content[0]!="end" && natom<MAX_ATOMS){
      type[natom]=content[0];
      x[natom*3]=atof(content[1].c_str());
      x[natom*3+1]=atof(content[2].c_str());
      x[natom*3+2]=atof(content[3].c_str());
      natom++;
      fin.getline(temp,MAX_CHARACTER); 
      parse(temp,column,content);
    }
    if (natom==MAX_ATOMS){
      cout<<"WARNING: the number of atoms is exceeding " <<MAX_ATOMS<<endl;
      cout<<"         anything larger than "<<MAX_ATOMS<<" is ignored"<<endl;
    }
    fin.seekg(beg_pos);
    v_pos=find_pattern(fin,string("Units bohrs"),"contains",false);
    if (v_pos>-1){
      bohr=true;
      cutoff/=bohr2a;
    }
  }

  double sx[3]={0,0,0};
  int ns=0;

  double *x0=&x[id0*3];
  string sum="";
  for (int i=0; i<natom; i++){
    string type1=type[i];
    double *x1=&x[i*3];
    double dx[3];
    dx[0]=x1[0]-x0[0];
    dx[1]=x1[1]-x0[1];
    dx[2]=x1[2]-x0[2];
    double r=sqrt((dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]));
    if (r<cutoff){
      cout<<i<<" "<<type1<<" "<<r*bohr2a<<endl;
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
  cout<<"add atom "<<x0[0]+sx[0]/ls*bondlength<<" "<<x0[1]+sx[1]/ls*bondlength<<" "<<x0[2]+sx[2]/ls*bondlength<<endl;
}
