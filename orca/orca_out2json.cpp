//this code is used to analyse the Print[ P_MOs ] 1
//output of ORCA. The density of state will be decomposed
//to different element and different angular momentum , s p d f
//this is called PDOS in the solid state dft community
//usage: MO <inputfile name> <outputfile name> 
//author: Lixin Sun nw13mifaso@gmail.com

#include "functions.h"

int main(int argc, char **argv){

    char temp[MAX_CHARACTER], * pch;
    string* content=new string[MAX_COLUMN];
    int column, line=0;
    bool Print_MO=false;
    
    //locate the beginning of the print_MO
    ifstream In1(argv[1]);
    if ( !In1.good()){
        cout<< " the input file does not exist or is corrupted..."<<endl;
        return 1;
    }

    //first check whether the result is converge
    string pattern="ERROR";
    int beg_pos=In1.tellg();
    int v_pos=find_pattern(In1,pattern,"contains",true);
    if (v_pos!=-1){
      cout<<"the calculation is not converged"<<endl;
      return 1;
    }else{
      In1.seekg(beg_pos);
    }

    pattern="FINAL SINGLE POINT ENERGY";
    int fe_pos=find_pattern(In1,pattern,"contains",true);
    if (fe_pos!=-1){
      In1.getline(temp,MAX_CHARACTER);
      if ( strstr(temp,"fully converged!")){
        cout<<"the calculation is not converged"<<endl;
        return 1;
      }
      break_line(temp,content);
      In1.seekg(beg_pos);
    }else{
      cout<<"the calculation is not converged"<<endl;
      return 1;
    }

    double *energy=new double[2*MAX_ENERGYLINE];
    double *occupancy=new double[2*MAX_ENERGYLINE];
    int *nstate=new int[2];
    double *homo=new double[2];
    double *lumo=new double[2];
    int *homo_i=new int[2];
    int *lumo_i=new int[2];
    bool read=read_orbital(In1,nstate,energy,occupancy,homo,lumo,homo_i,lumo_i);
    if (read==false){
      return 1;
    }else{
      In1.seekg(beg_pos);
    }

    ofstream json_o(argv[2]);
    if ( !json_o.good()){
        cout<< " the output file does not exist or is corrupted..."<<endl;
        return 1;
    }
    json_o<<std::fixed;


    json_o.precision(5);
    json_o<<"{"<<endl;
    json_o<<"\"energy\":"<<atof(content[4].c_str())*Eh2eV<<","<<endl;
    json_o<<"\"homo\":["<<homo[0]<<","<<homo[1]<<"],"<<endl;
    json_o<<"\"lumo\":["<<lumo[0]<<","<<lumo[1]<<"],"<<endl;
    json_o<<"\"homo_id\":["<<homo_i[0]<<","<<homo_i[1]<<"],"<<endl;
    json_o<<"\"lumo_id\":["<<lumo_i[0]<<","<<lumo_i[1]<<"],"<<endl;
    json_o<<"\"gap\":["<<lumo[0]-homo[0]<<","<<lumo[1]-homo[1]<<"],"<<endl;

    json_o<<"\"outputfile\":\""<<argv[1]<<"\","<<endl;
    json_o<<"\"software\":\"ORCA";
    pattern="Program Version";
    v_pos=find_pattern(In1,pattern,"contains",true);
    In1.getline(temp,MAX_CHARACTER); 
    break_line(temp,content);
    json_o<<content[2]<<"\","<<endl;

    json_o<<"\"comp_info\":[";
    pattern="Your calculation utilizes ";
    int len_pat=pattern.length();
    int u_count=0;
    int u_pos=find_pattern(In1,pattern,"contains",true);
    while (u_pos!=-1 && In1.good()){
      In1.getline(temp,MAX_CHARACTER); 
      if (u_count>0){
        json_o<<", ";
      }
      json_o<<"\""<<&temp[len_pat]<<"\"";
      u_count++;
      u_pos=find_pattern(In1,pattern,"contains",true);
    }
    json_o<<"],"<<endl;

    pattern="CARTESIAN COORDINATES (ANGSTROEM)";
    int c_pos=find_pattern(In1,pattern);
    int atomn=0;
    double *x=new double[3000];
    int *type=new int[1000];
    char element[MAX_ELEMENT][10];
    int n_element=0;
    int eid=0;
    if (c_pos!=-1){
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      column=break_line(temp,content);
      json_o<<"\"coordinates\":[";
      char species[100000]="\"species\":[";
      while (column==4){
        if (atomn>0){
          json_o<<", ";
          strcat(species,",");
        }
        json_o<<"["<<content[1]<<","<<content[2]<<","<<content[3]<<"]";
        strcat(species,"\"");
        strcat(species,content[0].c_str());
        strcat(species,"\"");
        x[atomn*3]=atof(content[1].c_str());
        x[atomn*3+1]=atof(content[2].c_str());
        x[atomn*3+2]=atof(content[3].c_str());
        //recognize the element name
        int elementid=-1;
        for (int eid=0;eid<n_element;eid++){
            if (strcmp(content[0].c_str(),element[eid])==0) elementid=eid;
        }
        if (elementid==-1){
            strcpy(element[n_element],content[0].c_str());
            elementid=n_element;
            n_element++;
        }
        type[atomn]=elementid;

        atomn++;
        In1.getline(temp,MAX_CHARACTER);
        column=break_line(temp,content);
      }
      json_o<<"],"<<endl;
      strcat(species,"],");
      json_o<<species<<endl;
      json_o<<"\"atomn\":"<<atomn<<","<<endl;

      //coordination
      double cutoff=DEFAULT_CUTOFF;
      int *icc=new int [atomn];
      for (int i=0;i<atomn;i++){
        icc[i]=0;
      }
      for (int i=0;i<atomn;i++){
        double *xxx=&x[i*3];
        int t=type[i];
        for (int j=i;j<atomn;j++){
          double *xxx1=&x[j*3];
          int t1=type[j];
          if (t!=t1){
            double dr=(xxx1[0]-xxx[0])*(xxx1[0]-xxx[0]);
            dr += (xxx1[1]-xxx[1])*(xxx1[1]-xxx[1]);
            dr += (xxx1[2]-xxx[2])*(xxx1[2]-xxx[2]);
            dr = sqrt(dr);
            if ( dr < cutoff ) {
              icc[j]+=1;
              icc[i]+=1;
            }
          }
        }
      }
      json_o << "\"coord\":[";
      for (int i=0;i<atomn;i++){
        if (i>0){
          json_o<<",";
        }
        json_o<<"\""<<icc[i]<<"\"";
      }
      json_o<<"],"<<endl;
      delete [] icc;
      delete [] x;
      delete [] type;
    }

    pattern="Hamiltonian:";
    int ham_pos=find_pattern(In1,pattern);
    if (ham_pos!=-1){
      In1.getline(temp,MAX_CHARACTER);
      column=break_line(temp,content);
      json_o<<"\"hamiltonian\":[";
      int nha=0;
      while(column>0){
        if (temp[1]!=' '){
          if (nha>0) json_o<<",";
          json_o<<"\""<<content[0];
          if (content[column-1].compare("on")!=0){
            json_o<<" "<<content[column-1];
          }
          json_o<<"\"";
          nha++; 
        }
        In1.getline(temp,MAX_CHARACTER);
        column=break_line(temp,content);
      }
      json_o<<"],"<<endl;
    }



    pattern="SCF CONVERGENCE";
    int scf_pos=find_pattern(In1,pattern);
    In1.getline(temp,MAX_CHARACTER);
    In1.getline(temp,MAX_CHARACTER);

    In1.getline(temp,MAX_CHARACTER);
    break_line(temp,content);
    json_o<<"\""<<content[0]<<"_"<<content[1]<<"_"<<content[2]<<"\":"<<content[4]<<","<<endl;
    In1.getline(temp,MAX_CHARACTER);
    break_line(temp,content);
    json_o<<"\""<<content[0]<<"_"<<content[1]<<"_"<<content[2]<<"\":"<<content[4]<<","<<endl;
    In1.getline(temp,MAX_CHARACTER);
    break_line(temp,content);
    json_o<<"\""<<content[0]<<"_"<<content[1]<<"_"<<content[2]<<"\":"<<content[4]<<","<<endl;
    In1.getline(temp,MAX_CHARACTER);
    break_line(temp,content);
    json_o<<"\""<<content[0]<<"_"<<content[1]<<"_"<<content[2]<<"\":"<<content[4]<<","<<endl;

    pattern="UHF SPIN CONTAMINATION"; 
    int spin_pos=find_pattern(In1,pattern);
    if (spin_pos!=-1){
      int s2_pos=find_pattern(In1,"Expectation value of","contains",true);
      In1.getline(temp,MAX_CHARACTER);
      break_line(temp,content);
      json_o<<"\"S2\":"<<content[5]<<","<<endl;
    }

    pattern="MULLIKEN POPULATION ANALYSIS";
    int mq_pos=find_pattern(In1,pattern);
    if (mq_pos!=-1){
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      char mulliken_q[100000],mulliken_s[100000];
      strcpy(mulliken_q,"\"mulliken_charge\":[");
      strcpy(mulliken_s,"\"mulliken_spin\":[");
      for (int i=0;i<atomn;i++){
        In1.getline(temp,MAX_CHARACTER);
        break_line(temp,content);
        if (i>0) {
          strcat(mulliken_q,", ");
          strcat(mulliken_s,", ");
        }
        if (content[2]!=":"){
          strcat(mulliken_q,content[2].c_str());
          strcat(mulliken_s,content[3].c_str());
        }else{
          strcat(mulliken_q,content[3].c_str());
          strcat(mulliken_s,content[4].c_str());
        }
      }
      json_o<<mulliken_q<<"],"<<endl;
      json_o<<mulliken_s<<"],"<<endl;
    }

    pattern="LOEWDIN POPULATION ANALYSIS";
    mq_pos=find_pattern(In1,pattern);
    if (mq_pos!=-1){
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      char loewdin_q[100000],loewdin_s[100000];
      strcpy(loewdin_q,"\"loewdin_charge\":[");
      strcpy(loewdin_s,"\"loewdin_spin\":[");
      for (int i=0;i<atomn;i++){
        In1.getline(temp,MAX_CHARACTER);
        break_line(temp,content);
        if (i>0) {
          strcat(loewdin_q,", ");
          strcat(loewdin_s,", ");
        }
        if (content[2]!=":"){
          strcat(loewdin_q,content[2].c_str());
          strcat(loewdin_s,content[3].c_str());
        }else{
          strcat(loewdin_q,content[3].c_str());
          strcat(loewdin_s,content[4].c_str());
        }
      }
      json_o<<loewdin_q<<"],"<<endl;
      json_o<<loewdin_s<<"],"<<endl;
    }

    pattern="MAYER POPULATION ANALYSIS";
    mq_pos=find_pattern(In1,pattern);
    if (mq_pos!=-1){
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      char mayer_q[100000];
      strcpy(mayer_q,"\"mayer_charge\":[");
      for (int i=0;i<atomn;i++){
        In1.getline(temp,MAX_CHARACTER);
        break_line(temp,content);
        if (i>0) {
          strcat(mayer_q,", ");
        }
        strcat(mayer_q,content[4].c_str());
      }
      json_o<<mayer_q<<"],"<<endl;
    }

    pattern="Summary of Natural Population Analysis";
    int nbo_pos=find_pattern(In1,pattern);
    if (nbo_pos!=-1){
      char nbo_q[100000],nbo_s[100000];
      strcpy(nbo_q,"\"nbo_charge\":[");
      strcpy(nbo_s,"\"nbo_spin\":[");
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      for (int i=0;i<atomn;i++){
        In1.getline(temp,MAX_CHARACTER);
        break_line(temp,content);
        if (i>0) {
          strcat(nbo_q,", ");
          strcat(nbo_s,", ");
        }
        strcat(nbo_q,content[2].c_str());
        strcat(nbo_s,content[7].c_str());
      }
      json_o<<nbo_q<<"],"<<endl;
      json_o<<nbo_s<<"],"<<endl;
    }

    pattern="DFT-D V3";
    json_o<<std::fixed;
    int d3_pos=find_pattern(In1,pattern);
    if (d3_pos!=-1){
      d3_pos=find_pattern(In1,"Edisp/kcal","contains",true);
      In1.getline(temp,MAX_CHARACTER);
      break_line(temp,content);
      json_o<<"\"edisp\":"<<atof(content[2].c_str())*Eh2eV<<","<<endl;
      In1.getline(temp,MAX_CHARACTER);
      break_line(temp,content);
      json_o<<"\"e6\":"<<atof(content[3].c_str())*kcal2eV<<","<<endl;
      In1.getline(temp,MAX_CHARACTER);
      break_line(temp,content);
      json_o<<"\"e8\":"<<atof(content[3].c_str())*kcal2eV<<","<<endl;
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      break_line(temp,content);
      json_o<<"\"dispersion_corr\":"<<atof(content[2].c_str())*Eh2eV<<","<<endl;
    }


    pattern="CARTESIAN GRADIENT";
    c_pos=find_pattern(In1,pattern);
    if (c_pos!=-1){
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      column=break_line(temp,content);
      int catomn=0;
      json_o<<"\"forces\":[";
      while (column==6){
        if (catomn>0){
          json_o<<", ";
        }
        catomn++;
        json_o<<"["<<atof(content[3].c_str())*Ehau2eVa<<","<<atof(content[4].c_str())*Ehau2eVa<<","<<atof(content[5].c_str())*Ehau2eVa<<"]";
        In1.getline(temp,MAX_CHARACTER);
        column=break_line(temp,content);
      }
      json_o<<"],"<<endl;
    }

    pattern="DIPOLE MOMENT";
    int dip_pos=find_pattern(In1,pattern);
    if (dip_pos!=-1){
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      break_line(temp,content);
      json_o<<"\"dipole_elec\":["<<atof(content[2].c_str())*bohr2a<<","<<atof(content[3].c_str())*bohr2a<<","<<atof(content[4].c_str())*bohr2a<<"],"<<endl;
      In1.getline(temp,MAX_CHARACTER);
      break_line(temp,content);
      json_o<<"\"dipole_nuc\":["<<atof(content[3].c_str())*bohr2a<<","<<atof(content[4].c_str())*bohr2a<<","<<atof(content[5].c_str())*bohr2a<<"],"<<endl;
      In1.getline(temp,MAX_CHARACTER);
      In1.getline(temp,MAX_CHARACTER);
      break_line(temp,content);
      json_o<<"\"dipole\":["<<atof(content[4].c_str())*bohr2a<<","<<atof(content[5].c_str())*bohr2a<<","<<atof(content[6].c_str())*bohr2a<<"]"<<endl;
    }

    json_o<<"}"<<endl;
}

