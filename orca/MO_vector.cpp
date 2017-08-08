//this code is used to analyse the Print[ P_GuessOrb ] 1
//the top three maximum component of each orbital will be output
//usage: orca_guessorb <inputfile name> <outputfile name> 
//author: Lixin Sun nw13mifaso@gmail.com

#include <algorithm>
#include <iostream>
#include <locale>
#include <fstream>          // file I/O suppport
#include <cstdlib>          // support for exit()
#include <stdio.h>
#include <sys/timeb.h>
#include <sys/types.h>
#include <time.h>
#include <cmath>
#include <malloc.h>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <sstream>
using namespace std;

#include "math.h"
#include "stdlib.h"       // for random
#include <vector>

//set up the size of array used for storage
#define MAX_ELEMENT 10
#define MAX_M   6
#define MAX_ENERGYLINE 4000
#define TEMP_ENERGYLINE 10

//set up the final plotting region
//the fermi level is shfited to zero
#define EMIN -30
#define EMAX 15

//set up the buffer size for reading
#define MAX_CHARACTER 1000
#define MAX_COLUMN 20

bool isxyz(char c){
    if ((c=='x')||(c=='y')||(c=='z')) return true;
    else return false;
}

bool isdigit1(char c){
    if ((c>='0')&&(c<='9')) return true;
    else return false;
}

void gaussian(int grid, double *x, double *y, double x0, double sigma2){
    for (int i=0;i<grid;i++) y[i]=exp(-(x[i]-x0)*(x[i]-x0)/2./sigma2);
}

bool contain_alphabet(char * c){
    for (int i=0;i<strlen(c);i++){
        if (((c[i]>='a')&&(c[i]<='z'))||((c[i]>='A')&&(c[i]<='Z')))
            return true;
    }
    return false;
}

int main(int argc, char **argv){

    char temp[MAX_CHARACTER], * pch,content[MAX_COLUMN][MAX_CHARACTER];
    int column, line=0;
    bool Print_MO=false;
    
    //locate the beginning of the print_MO
    ifstream In1(argv[1]);

    if ( !In1.good()){
        cout<< " the input file does not exist or is corrupted..."<<endl;
        return 1;
    }
    while (In1.good()){
        line++;
        In1.getline(temp,MAX_CHARACTER); 

        column=0;
        pch = strtok (temp," ");
        while ((pch != NULL)&&(column<MAX_COLUMN)) {
            strcpy(content[column],pch);
            column++;
            pch = strtok (NULL, " ");
        }
        pch= NULL;
        if (column == 2 ){
            if (( strcmp(content[0],"MOLECULAR") == 0 ) && (strcmp(content[1],"ORBITALS") == 0)){
                Print_MO=true;
                break;
            }
        }
    }
    if ((! In1.good()) || (Print_MO==false)){
        cout<<"ERROR: MOLECULAR ORBITALS block not found"<<endl;
        cout<<"please write below lines to your orca input script "<<endl;
        cout<<"\%output"<<endl;
        cout<<"Print[ P_Basis ] 2"<<endl;
        cout<<"Print[ P_MOs ] 1"<<endl;
        cout<<"end"<<endl;
        return 1;
    }
    In1.getline(temp,MAX_CHARACTER); //another buffering line
    line++;

    double energy[2*MAX_ENERGYLINE];
    double occupancy[2*MAX_ENERGYLINE];
    double max_coeff2[3*MAX_ENERGYLINE*2];
    string max_name[3*MAX_ENERGYLINE*2];
    int n_energy[2]={0,0};
    double efermi=-100000;

    fill(max_coeff2,max_coeff2+3*MAX_ENERGYLINE,0);

    double norm[TEMP_ENERGYLINE];
    
    //spin up and down
    for (int spin=0;spin<2;spin++){
      int energy_id=0;
      int energy_line=0;
      double *p_energy;
      double *p_occupancy;
      double *p_tally;
      double *p_maxcoeff2;
      string *p_maxname;
      bool start_sum=false;
      int current_column=0;
      do{
        line++;
        In1.getline(temp,MAX_CHARACTER); 
        if ( temp[0]=='-' && temp[1]=='-') break;

        column=0;
        pch = strtok (temp," ");
        while (pch != NULL) {
            strcpy(content[column],pch);
            column++;
            pch = strtok (NULL, " ");
        }
        pch= NULL;
        
        if (column!=0){
          if (contain_alphabet(content[0])==false) energy_line+=1;

          //the second line is for orbital energy
          if (energy_line == 2){
            p_energy=energy+spin*MAX_ENERGYLINE+energy_id;
            p_occupancy=occupancy+spin*MAX_ENERGYLINE+energy_id;
            p_maxcoeff2=max_coeff2+spin*MAX_ENERGYLINE*3+energy_id*3;
            p_maxname=max_name+spin*MAX_ENERGYLINE*3+energy_id*3;
            for (int columni=0;columni<column;columni++){
               p_energy[columni]=atof(content[columni])*27.2113834; 
            }

            //clear the temporary sum container. because this sum is not normalized.
            start_sum=true;
            current_column=column;
            fill(norm,norm+TEMP_ENERGYLINE,0);

            n_energy[spin]+=column;
            energy_id+=column;

          //it starts with non-space character
          } else if ((column>0) && (energy_line ==0)  && (contain_alphabet(content[0])==true)){

            //sum over 
            for (int columni = 2; columni < column; columni++){
              double coeff=atof(content[columni]);
              double coeff2=coeff*coeff;

              norm[columni-2]+=coeff2;
              int icolumn3=(columni-2)*3;
              double *max = p_maxcoeff2+icolumn3;
              string *name = p_maxname+icolumn3;
              bool find_max=false;
              for (int top=0;((top<3) &&(find_max==false)); top++){
                if (coeff2 > max[top] ){
                    for (int next=2;next<=(top+1);next++){
                       max[next] = max[next-1]; 
                       //sprintf(temp,"%s_%s",content[0],content[1]);
                       name[next]=content[0];
                       name[next]=name[next]+'_';
                       name[next]=name[next]+content[1];
                    }
                    max[top] = coeff2;
                    name[top] = content[0];
                    name[top] = name[top]+'_';
                    name[top] = name[top]+content[1];
                    find_max=true;
                }
              }
            }
          }else if (energy_line == 3){
            for (int columni=0;columni<column;columni++){
                p_occupancy[columni]=atof(content[columni]);
                if ((p_occupancy[columni]>0) && (p_energy[columni] > efermi)){
                        efermi=p_energy[columni];
                }
            }
          }else if (energy_line == 4 ) energy_line=0;
          else if ((energy_line ==1)||((energy_line==0) && (contain_alphabet(content[0])==false))){
            if (start_sum==true){
              for (int columni=0;columni<current_column;columni++){
                for (int top=0;top<3; top++){
                  p_maxcoeff2[columni*3+top]/=norm[columni];
                }
              }
              start_sum=false;
            }
          }
        }
      }while (column!=0);

      if (start_sum==true){
        for (int columni=0;columni<current_column;columni++){
          for (int top=0;top<3; top++){
            p_maxcoeff2[columni*3+top]/=norm[columni];
          }
        }
        start_sum=false;
      }
    }
    cout<<"done reading"<<endl;


    ofstream out(argv[2]);
    // output the raw pdos without smearing
    for (int spin=0; spin<2; spin++){
      double *p_energy=energy+spin*MAX_ENERGYLINE;
      double *p_occupancy=occupancy+spin*MAX_ENERGYLINE;
      double *p_maxcoeff2=max_coeff2+(spin*MAX_ENERGYLINE)*3;
      string *p_maxname=max_name+(spin*MAX_ENERGYLINE)*3;
      for (int i=0;i<n_energy[spin];i++){
        out<<i<<" "<<p_energy[i]<<" "<<p_occupancy[i];
        for (int top=0;top<3;top++){
            out<<" "<<p_maxname[i*3+top]<<" "<<p_maxcoeff2[i*3+top];
        }
        out<<endl;
      }
    }
}
