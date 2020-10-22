
//author: Lixin Sun nw13mifaso@gmail.com

#include <algorithm>
#include <iostream>
#include <locale>
#include <iomanip>
#include <fstream>          // file I/O suppport
#include <string>
#include <sys/timeb.h>
#include <sys/types.h>

#include <cstdlib>          // support for exit()
#include <cstdio>
#include <ctime>
#include <cmath>
#include <malloc.h>
#include <cstring>
#include <sstream>
using namespace std;

#include <vector>

#define ANG2BOHR 0.52917721092
#define MAX_CHARACTER 1000
#define MAX_COLUMN 200

void parse(char * temp, int & column, string *content){
  char temp0[MAX_CHARACTER];
  strcpy(temp0,temp);
  column=0;
  char * pch;
  pch = strtok (temp0," ");
  while ((pch != NULL)&&(column<MAX_COLUMN)) {
      content[column]=pch;
      column++;
      pch = strtok (NULL, " ");
  }
  pch= NULL;
}

