//Subroutine reads the input file.
//This subroutine reads the parameter from the input files which will 
//be used by the program. To include more parameters in the input file
//change the program apropriately

#include "functions.h"

void read_input()
{
ifstream fin;
string line;

fin.open("input_file.dat");

getline(fin,line);
getline(fin,line);
fin>>grid_name;    
getline(fin,line);
getline(fin,line);
fin>>Niter;
getline(fin,line);
fin>>tol;
}

