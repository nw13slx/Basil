
//author: Lixin Sun nw13mifaso@gmail.com

#include "definition.h"

// only take xyz format

int main(int argc, char **argv){
  if (argc < 3){
    cout<<" FAILED: need more input arguments"<<endl;
    return 1;
  }
  
  ifstream fin(argv[1]); 
  ofstream fout(argv[2]); 
  //qm region number
  string qm(argv[3]);
  double cutoff=2;

  char temp[MAX_CHARACTER];
  string content[MAX_COLUMN];
  int column, line=0;

  fin.getline(temp,MAX_CHARACTER); 
  parse(temp,column,content);
  int natom=atoi(content[0].c_str());

  fin.getline(temp,MAX_CHARACTER); 
  string *sym=new string[natom];
  double * x= new double [3*natom];
  double * q= new double [natom];
  bool * isQM = new bool [natom];
  int *QM = new int [natom];
  int nQM=0;
  int *MM = new int [natom];
  int nMM=0;
  for (int i=0;i<natom;i++){
    fin.getline(temp,MAX_CHARACTER); 
    parse(temp,column,content);
    sym[i]=content[0];
    x[i*3]=atof(content[1].c_str());
    x[i*3+1]=atof(content[2].c_str());
    x[i*3+2]=atof(content[3].c_str());
    if (content[0].find(qm)!=std::string::npos){
      QM[nQM]=i;
      nQM++;
      isQM[i]=true;
    }else{ //if (content[0].compare(sym2)==0){
      MM[nMM]=i;
      nMM++;
      isQM[i]=false;
    }
    if (column>4){
      q[i]=atof(content[4].c_str());
    }
  }
  fin.close();
  cout<<nQM<<" "<<nMM<<endl;

  int nh=0;
  double *h=new double[nQM*9];
  for (int i=0; i<nQM; i++){
    double *xQM=&x[QM[i]*3];
    int ncoord=0;
    for (int j=0; j<nMM; j++){
      double *xMM=&x[MM[j]*3];
      double dx[3];
      dx[0]=xMM[0]-xQM[0];
      dx[1]=xMM[1]-xQM[1];
      dx[2]=xMM[2]-xQM[2];
      double r=sqrt((dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]));
      if (r<cutoff){
        ncoord++;
        double dr=1.*r/1.91566;
        h[nh*3]=dx[0]/r*dr+xQM[0];
        h[nh*3+1]=dx[1]/r*dr+xQM[1];
        h[nh*3+2]=dx[2]/r*dr+xQM[2];
        nh++;
        if (q[MM[j]]>0){
          q[MM[j]]-=2.23/6.;
        }else{
          q[MM[j]]+=2.23/6.;
        }
      }
    }
  }
  
  fout<<nh+nQM<<endl<<endl;
  for (int i=0;i<nQM;i++) {
    fout<<sym[QM[i]]<<" "<<x[QM[i]*3]<<" "<<x[QM[i]*3+1]<<" "<<x[QM[i]*3+2]<<endl;
  }
  for (int i=0;i<nh;i++){
    fout<<"H "<<h[i*3]<<" "<<h[i*3+1]<<" "<<h[i*3+2]<<endl;
  }

  fout<<natom-nQM<<endl<<endl;
  for (int i=0;i<nMM;i++) {
    if (isQM[MM[i]]==false)
      fout<<sym[MM[i]]<<" "<<q[MM[i]]<<" "<<x[MM[i]*3]<<" "<<x[MM[i]*3+1]<<" "<<x[MM[i]*3+2]<<endl;
  }

  fout.close();
}
