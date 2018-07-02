// analyse average bond length between sym1-sym2, for inside, outside and across the QM/MM boundary
// input format, pun file  for chemshell
//usage: bondlength  input symbol1 symbol2 cutoff
//author: Lixin Sun nw13mifaso@gmail.com

#include "definition.h"

int main(int argc, char **argv){
  if (argc < 4){
    cout<<" FAILED: need more input arguments"<<endl;
    cout<<"usage: bondlength  input symbol1 symbol2 cutoff"<<endl;
    return 1;
  }
  ifstream fin(argv[1]); 
  string sym1(argv[2]);
  string sym2(argv[3]);
  double cutoff;
  if (argc < 5 ) cutoff=3;
  else cutoff=atof(argv[4]);

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

  // QM-MM bond
  double tally_12=0;
  double tally2_12=0;
  int n_12=0;
  for (int i=0; i<nQM; i++){
    double *xQM=&x[QM[i]*3];
    string type1=type[QM[i]];
    for (int j=0; j<nMM; j++){
      string type2=type[MM[j]];
      if ((type1.find(sym1)!=std::string::npos) && ( type2.find(sym2)!=std::string::npos) ||(type1.find(sym2)!=std::string::npos) && ( type2.find(sym1)!=std::string::npos)) {
        double *xMM=&x[MM[j]*3];
        double dx[3];
        dx[0]=xQM[0]-xMM[0];
        dx[1]=xQM[1]-xMM[1];
        dx[2]=xQM[2]-xMM[2];
        double r=sqrt((dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]));
        if (r<cutoff/ANG2BOHR){
          n_12++;
          tally_12+=r;
          tally2_12+=r*r;
        }
      }
    }
  }
  tally_12=tally_12*ANG2BOHR/double(n_12); 
  tally2_12=tally2_12*ANG2BOHR*ANG2BOHR/double(n_12); 
  cout<<"QM-MM "<< tally_12<<" +/- "<<sqrt(tally2_12-tally_12*tally_12)<<endl;

  // QM-QM bond
  double tally_11=0;
  double tally2_11=0;
  int n_11=0;
  for (int i=0; i<nQM; i++){
    double *xQM=&x[QM[i]*3];
    string type1=type[QM[i]];
    for (int j=i+1; j<nQM; j++){
      string type2=type[QM[j]];
      if ((type1.find(sym1)!=std::string::npos) && ( type2.find(sym2)!=std::string::npos) ||(type1.find(sym2)!=std::string::npos) && ( type2.find(sym1)!=std::string::npos)) {
        double *xQM2=&x[QM[j]*3];
        double dx[3];
        dx[0]=xQM[0]-xQM2[0];
        dx[1]=xQM[1]-xQM2[1];
        dx[2]=xQM[2]-xQM2[2];
        double r=sqrt((dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]));
        if (r<cutoff/ANG2BOHR){
          n_11++;
          tally_11+=r;
          tally2_11+=r*r;
        }
      }
    }
  }
  tally_11=tally_11*ANG2BOHR/double(n_11); 
  tally2_11=tally2_11*ANG2BOHR*ANG2BOHR/double(n_11); 
  cout<<"QM-QM "<< tally_11<<" +/- "<<sqrt(tally2_11-tally_11*tally_11)<<endl;

  // MM-MM bond
  double tally_22=0;
  double tally2_22=0;
  int n_22=0;
  for (int i=0; i<nMM; i++){
    double *xMM1=&x[MM[i]*3];
    string type1=type[MM[i]];
    for (int j=i+1; j<nMM; j++){
      string type2=type[MM[j]];
      if ((type1.find(sym1)!=std::string::npos) && ( type2.find(sym2)!=std::string::npos) ||(type1.find(sym2)!=std::string::npos) && ( type2.find(sym1)!=std::string::npos)) {
        double *xMM2=&x[MM[j]*3];
        double dx[3];
        dx[0]=xMM1[0]-xMM2[0];
        dx[1]=xMM1[1]-xMM2[1];
        dx[2]=xMM1[2]-xMM2[2];
        double r=sqrt((dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]));
        if (r<cutoff/ANG2BOHR){
          n_22++;
          tally_22+=r;
          tally2_22+=r*r;
        }
      }
    }
  }
  tally_22=tally_22*ANG2BOHR/double(n_22); 
  tally2_22=tally2_22*ANG2BOHR*ANG2BOHR/double(n_22); 
  cout<<"MM-MM "<< tally_22<<" +/- "<<sqrt(tally2_22-tally_22*tally_22)<<endl;

}
