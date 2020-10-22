
//author: Lixin Sun nw13mifaso@gmail.com

#include "definition.h"

int main(int argc, char **argv){
  if (argc < 2){
    cout<<" FAILED: need more input arguments"<<endl;
    return 1;
  }
  ifstream fin(argv[1]); 
  ofstream fout(argv[2]); 
  bool frac=true;
  if (argc >3){
    if ( strcmp (argv[3],"frac")==0){
      frac=true;
    }else{
      frac=false;
    }
  }

  char temp[MAX_CHARACTER];
  string content[MAX_COLUMN];
  int column, line=0;
  fin.getline(temp,MAX_CHARACTER); 

  double scale,b[9];
  fin>>scale;
  for (int i=0;i<9;i++) fin >>b[i];
  fin.getline(temp,MAX_CHARACTER); 

  double cell[6];
  for (int i=0; i<3;i++){
      int i3=i*3;
      cell[i]=double(sqrt(b[i3]*b[i3]+b[i3+1]*b[i3+1]+b[i3+2]*b[i3+2]))*scale;
  }
  cell[3]=acos((b[6]*b[3]+b[7]*b[4]+b[8]*b[5])/cell[2]/cell[1])/M_PI*180.;
  cell[4]=acos((b[0]*b[6]+b[1]*b[7]+b[2]*b[8])/cell[0]/cell[2])/M_PI*180.;
  cell[5]=acos((b[0]*b[3]+b[1]*b[4]+b[2]*b[5])/cell[0]/cell[1])/M_PI*180.;
  fout<<"cell"<<endl;
  for (int i=0; i<6;i++){
      fout<<cell[i]<<" ";
  }
  fout<<endl;

  //get name
  string element[10];
  fin.getline(temp,MAX_CHARACTER); 
  column=0;
  parse(temp,column,content);
  for (int i=0;i<column;i++){ element[i]=content[i];}

  //get natom
  fin.getline(temp,MAX_CHARACTER); 
  parse(temp,column,content);
  int nelement=column;
  double *natom=new double[column];
  for (int i=0;i<column;i++){
    natom[i]=atoi(content[i].c_str());
  }

  fin.getline(temp,MAX_CHARACTER); 
  parse(temp,column,content);
  if ((content[0][0]=='s') || (content[0][0]=='S') )
    fin.getline(temp,MAX_CHARACTER); 

  //read ions
  if (frac){
    fout<<"frac"<<endl;
    double x[3];
    for (int ele=0; ele<nelement;ele++){
      for (int i=0;i<natom[ele];i++){
        fin>>x[0]>>x[1]>>x[2];
        fin.getline(temp,MAX_CHARACTER); 
        fout<<element[ele]<<" "<<x[0]*cell[0]<<" "<<x[1]*cell[1]<<" "<<x[2]*cell[2]<<endl;
      }
    }
    fin.getline(temp,MAX_CHARACTER); //empty line
  }{
    fout<<"frac"<<endl;
    double x[3];
    for (int ele=0; ele<nelement;ele++){
      for (int i=0;i<natom[ele];i++){
        fin>>x[0]>>x[1]>>x[2];
        fin.getline(temp,MAX_CHARACTER); 
        fout<<element[ele]<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
      }
    }
    fin.getline(temp,MAX_CHARACTER); //empty line
  }

}
