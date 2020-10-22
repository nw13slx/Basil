// analyse average bondangle between sym1-sym2-sym3, for inside, outside and across the QM/MM boundary
// input format, pun file  for chemshell
//usage: bondlength  input symbol1 symbol2 sym3 cutoff
//author: Lixin Sun nw13mifaso@gmail.com

#include "definition.h"

void angle(string *sym3, string *type3, double ** x3, double&tally, int &n,double cutoff,double cutoff2, double cutoff3,bool debug=false){
  if (debug){
    cout<<type3[0]<<" "<<type3[1]<<" "<<type3[2]<<" ";
  }
  if (type3[1].find(sym3[1])!=std::string::npos  && type3[2].find(sym3[2])!=std::string::npos && type3[0].find(sym3[0])!=std::string::npos){
    double dx21[3],dx23[3],r21,r23,r13,dx13[3];
    dx21[0]=x3[1][0]-x3[0][0];
    dx21[1]=x3[1][1]-x3[0][1];
    dx21[2]=x3[1][2]-x3[0][2];
    r21=sqrt((dx21[0]*dx21[0])+(dx21[1]*dx21[1])+(dx21[2]*dx21[2]));
    if (debug){
      cout<<r21<<" ";
    }
    if (r21<cutoff/ANG2BOHR){
      dx23[0]=x3[1][0]-x3[2][0];
      dx23[1]=x3[1][1]-x3[2][1];
      dx23[2]=x3[1][2]-x3[2][2];
      r23=sqrt((dx23[0]*dx23[0])+(dx23[1]*dx23[1])+(dx23[2]*dx23[2]));
      dx13[0]=x3[0][0]-x3[2][0];
      dx13[1]=x3[0][1]-x3[2][1];
      dx13[2]=x3[0][2]-x3[2][2];
      r13=sqrt((dx13[0]*dx13[0])+(dx13[1]*dx13[1])+(dx13[2]*dx13[2]));
      if (debug){
        cout<<r23<<" "<<r13<<" ";
      }
      if ((r23<cutoff/ANG2BOHR) && (r13>cutoff2/ANG2BOHR) && (r13<cutoff3/ANG2BOHR)) {
        double costheta=(dx21[0]*dx23[0]+dx21[1]*dx23[1]+dx21[2]*dx23[2])/r21/r23;
        double theta=acos(costheta)/M_PI*180;
        n++;
        tally+=theta;
      }
    }
  }
  if (debug){
    cout<<endl;
  }
}

int ** neighbor_list(double *x1, double *x2, int n1, int n2, double cutoff){
  int ** ngh_l = new int * [n1]; 
  for (int i=0;i<n1;i++){
    int ngh_temp[1000];
    int count=0;
    for (int j=0;j<n2;j++){
      double r,dx[3];
      dx[0]=x1[i*3]-x2[j*3];
      dx[1]=x1[i*3+1]-x2[j*3+1];
      dx[2]=x1[i*3+2]-x2[j*3+2];
      r=sqrt((dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]));
      if (r<cutoff/ANG2BOHR){
        ngh_temp[count]=j;
        count++;
      }
    }
    ngh_l[i]=new int [count+1];
    ngh_l[i][0]=count;
    for (int j=0;j<count;j++){
      ngh_l[i][j+1]=ngh_temp[j];
    }
  }
  return ngh_l;
}

int main(int argc, char **argv){
  if (argc < 5){
    cout<<" FAILED: need more input arguments"<<endl;
    return 1;
  }
  ifstream fin(argv[1]); 
  string sym3[3]={argv[2],argv[3], argv[4]};
  double cutoff,cutoff2,cutoff3;
  if (argc < 6 ) {
    cutoff=3;
    cutoff2=0;
    cutoff3=3;
  }else if (argc<7){
   cutoff=atof(argv[5]);
   cutoff2=0;
   cutoff3=cutoff;
  } else if (argc<8) {
   cutoff=atof(argv[5]);
   cutoff2=0;
   cutoff3=atof(argv[6]);
  } else{
   cutoff=atof(argv[5]);
   cutoff2=atof(argv[6]);
   cutoff3=atof(argv[7]);
  }

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
  double **x3=new double *[3];  
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

  double * MMx=new double[nMM*3];
  for (int i=0;i<nMM;i++){
    MMx[i*3]=x[MM[i]*3];
    MMx[i*3+1]=x[MM[i]*3+1];
    MMx[i*3+2]=x[MM[i]*3+2];
  }
  double * QMx=new double[nQM*3];
  for (int i=0;i<nQM;i++){
    QMx[i*3]=x[QM[i]*3];
    QMx[i*3+1]=x[QM[i]*3+1];
    QMx[i*3+2]=x[QM[i]*3+2];
  }
  int **ngh_list22=neighbor_list(MMx, MMx, nMM, nMM,cutoff);
  int **ngh_list12=neighbor_list(QMx, MMx, nQM, nMM,cutoff);
  int **ngh_list21=neighbor_list(MMx, QMx, nMM, nQM,cutoff);
  int **ngh_list11=neighbor_list(QMx, QMx, nQM, nQM,cutoff);
  delete [] MMx;
  delete [] QMx;

  double tally_111=0;
  int n_111=0;
  for (int i=0; i<nQM; i++){
    for (int j0=1; j0<=ngh_list11[i][0]; j0++){
      int j=ngh_list11[i][j0];
      for (int k0=j+1; k0<=ngh_list11[i][0]; k0++){
        int k=ngh_list11[i][k0];
        if (i!=j && i!=k){
          string type3[3];
          type3[0]=type[QM[j]];
          type3[1]=type[QM[i]];
          type3[2]=type[QM[k]];
          x3[0]=&x[QM[j]*3];
          x3[1]=&x[QM[i]*3];
          x3[2]=&x[QM[k]*3];
          angle(sym3, type3, x3, tally_111, n_111,cutoff,cutoff2,cutoff3);
          if (sym3[0]!=sym3[2]){
            type3[0]=type[QM[k]];
            type3[1]=type[QM[i]];
            type3[2]=type[QM[j]];
            x3[0]=&x[QM[k]*3];
            x3[1]=&x[QM[i]*3];
            x3[2]=&x[QM[j]*3];
            angle(sym3, type3, x3, tally_111, n_111,cutoff,cutoff2,cutoff3);
          } 
        }
      }
    }
  }
  tally_111=tally_111/double(n_111); 
  cout<<"QM-QM-QM "<< tally_111<<endl;

  double tally_211=0;
  int n_211=0;
  double tally_112=0;
  int n_112=0;
  for (int i=0; i<nQM; i++){
    for (int j0=1; j0<=ngh_list11[i][0]; j0++){
      int j=ngh_list11[i][j0];
      for (int k0=0; k0<=ngh_list12[i][0]; k0++){
        int k=ngh_list12[i][k0];
        string type3[3];
        type3[0]=type[QM[j]];
        type3[1]=type[QM[i]];
        type3[2]=type[MM[k]];
        x3[0]=&x[QM[j]*3];
        x3[1]=&x[QM[i]*3];
        x3[2]=&x[MM[k]*3];
        angle(sym3, type3, x3, tally_112, n_112,cutoff,cutoff2,cutoff3);
        if (sym3[0]!=sym3[2]){
          type3[0]=type[MM[k]];
          type3[2]=type[QM[j]];
          x3[0]=&x[MM[k]*3];
          x3[1]=&x[QM[i]*3];
          x3[2]=&x[QM[j]*3];
          angle(sym3, type3, x3, tally_211, n_211,cutoff,cutoff2,cutoff3);
        }
      }
    }
  }
  
  tally_112=tally_112/double(n_112); 
  cout<<"QM-QM-MM "<< tally_112<<endl;
  if (sym3[0]!=sym3[2]){
    tally_211=tally_211/double(n_211); 
    cout<<"MM-QM-QM "<< tally_211<<endl;
  }

  double tally_212=0;
  int n_212=0;
  for (int i=0; i<nQM; i++){
    for (int j0=1; j0<=ngh_list12[i][0]; j0++){
      int j=ngh_list12[i][j0];
      for (int k0=j+1; k0<=ngh_list12[i][0]; k0++){
        int k=ngh_list12[i][k0];
        string type3[3];
        type3[0]=type[MM[j]];
        type3[1]=type[QM[i]];
        type3[2]=type[MM[k]];
        x3[0]=&x[MM[j]*3];
        x3[1]=&x[QM[i]*3];
        x3[2]=&x[MM[k]*3];
        angle(sym3, type3, x3, tally_212, n_212,cutoff,cutoff2,cutoff3);
        if (sym3[0]!=sym3[2]){
          type3[0]=type[MM[k]];
          type3[2]=type[MM[j]];
          x3[0]=&x[MM[k]*3];
          x3[1]=&x[QM[i]*3];
          x3[2]=&x[MM[j]*3];
          angle(sym3, type3, x3, tally_212, n_212,cutoff,cutoff2,cutoff3);
        }
      }
    }
  }
  tally_212=tally_212/double(n_212); 
  cout<<"MM-QM-MM "<< tally_212<<endl;

  double tally_121=0;
  int n_121=0;
  for (int i=0; i<nMM; i++){
    for (int j0=1; j0<=ngh_list21[i][0]; j0++){
      int j=ngh_list21[i][j0];
      for (int k0=j+1; k0<=ngh_list21[i][0]; k0++){
        int k=ngh_list21[i][k0];
        string type3[3];
        type3[0]=type[QM[j]];
        type3[1]=type[MM[i]];
        type3[2]=type[QM[k]];
        x3[0]=&x[QM[j]*3];
        x3[1]=&x[MM[i]*3];
        x3[2]=&x[QM[k]*3];
        angle(sym3, type3, x3, tally_121, n_121,cutoff,cutoff2,cutoff3);
        if (sym3[0]!=sym3[2]){
          type3[0]=type[QM[k]];
          type3[2]=type[QM[j]];
          x3[0]=&x[QM[k]*3];
          x3[1]=&x[MM[i]*3];
          x3[2]=&x[QM[j]*3];
          angle(sym3, type3, x3, tally_121, n_121,cutoff,cutoff2,cutoff3);
        }
      }
    }
  }
  if (n_121>0){
   tally_121=tally_121/double(n_121); 
   cout<<"QM-MM-QM "<< tally_121<<endl;
  }

  double tally_221=0;
  int n_221=0;
  double tally_122=0;
  int n_122=0;
  for (int i=0; i<nMM; i++){
    for (int j0=1; j0<=ngh_list22[i][0]; j0++){
      int j=ngh_list22[i][j0];
      for (int k0=1; k0<=ngh_list21[i][0]; k0++){
        int k=ngh_list21[i][k0];
        if (i!=j){
          string type3[3];
          type3[0]=type[MM[j]];
          type3[1]=type[MM[i]];
          type3[2]=type[QM[k]];
          x3[0]=&x[MM[j]*3];
          x3[1]=&x[MM[i]*3];
          x3[2]=&x[QM[k]*3];
          angle(sym3, type3, x3, tally_221, n_221,cutoff,cutoff2,cutoff3);
          if (sym3[0]!=sym3[2]){
            type3[0]=type[QM[k]];
            type3[2]=type[MM[j]];
            x3[0]=&x[QM[k]*3];
            x3[1]=&x[MM[i]*3];
            x3[2]=&x[MM[j]*3];
            angle(sym3, type3, x3, tally_122, n_122,cutoff,cutoff2,cutoff3);
          }
        }
      }
    }
  }
  tally_221=tally_221/double(n_221); 
  cout<<"MM-MM-QM "<< tally_221<<endl;
  if (sym3[0]!=sym3[2]){
    tally_122=tally_122/double(n_122); 
    cout<<"QM-MM-MM "<< tally_122<<endl;
  }

  double tally_222=0;
  int n_222=0;
  for (int i=0; i<nMM; i++){
    for (int j0=1; j0<=ngh_list22[i][0]; j0++){
      int j=ngh_list22[i][j0];
      for (int k0=j+1; k0<=ngh_list22[i][0]; k0++){
        int k=ngh_list22[i][k0];
        if (i!=j){
          string type3[3];
          type3[0]=type[MM[j]];
          type3[1]=type[MM[i]];
          type3[2]=type[MM[k]];
          x3[0]=&x[MM[j]*3];
          x3[1]=&x[MM[i]*3];
          x3[2]=&x[MM[k]*3];
          angle(sym3, type3, x3, tally_222, n_222,cutoff,cutoff2,cutoff3);
        }
      }
    }
  }
  tally_222=tally_222/double(n_222); 
  cout<<"MM-MM-MM "<< tally_222<<endl;
  delete [] x;
  delete [] QM;
  delete [] MM;
  delete [] ngh_list11;
  delete [] ngh_list12;
  delete [] ngh_list21;
  delete [] ngh_list22;
}
