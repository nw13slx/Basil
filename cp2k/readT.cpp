//this code is used to read the velocity output file from CP2K (with xyz format)
//it will output the average temperature of atom group for each frame, and the mean and deviation over time
//
//usage: readT filename index1(default: 0) index2(default: natom-1)
//author: Lixin Sun nw13mifaso@gmail.com

#include "functions.h"
#include "atomic_mass.h"

int main(int argc, char **argv){
  
  if (argc < 1){
    cout<<"need the input file name"<<endl;
    return 1;
  }
  int natom,index1 = 0, index2 = -1;
  ifstream fin(argv[1]);
  if (argc >= 4){
    index1 = atoi(argv[2]);
    index2 = atoi(argv[3]);
  }
  int frame = 0;
  double *T = read_vel(fin, natom);
  double Ttally = 0;
  double T2tally = 0;

  cout << "# frame T_avg" <<endl;

  while (T != NULL){
    double T_avg = 0;
    if (index2 == -1) index2 = natom-1;
    int count = 0;
    for (int i = index1; i <= index2 ; i++){
      T_avg += T[i];
      count++;
    }
    T_avg /= double(count);
    Ttally += T_avg ;
    T2tally += T_avg*T_avg ;
    cout << frame++ << " " << T_avg <<endl;
    delete [] T;
    T = read_vel(fin, natom);
  }
  Ttally = Ttally / double (frame);
  T2tally = T2tally / double (frame);
  cout << "# Ttally sqrt(-Ttally*Ttally + T2tally) " <<endl;
  cout << "# " << Ttally << " " << sqrt(-Ttally*Ttally + T2tally) <<endl;
  return 0;
}
