//Subroutine for reading the mesh file for this particular format
//New subrouinte need to be written for different format

#include "functions.h"

void reading_grid()
{
int temp=1;	
std::string line;
ifstream fin;
fin.open(grid_name.c_str());


while(temp<=3)
{
getline(fin,line);
temp++;
}

fin>>ndimn>>ntype;
getline(fin,line);
getline(fin,line);
fin>>nelem>>npoint>>nface;
getline(fin,line);
getline(fin,line);

cout<<"Number of Elements:"<<nelem<<"\nNumber of Points: "<<npoint<<"\nNumber of Boundary faces "<<nface<<endl;

//-----------------------------------------------------------------
//Initializing element connectivity array and reading from the file
inpoel = new double* [nelem];
for (int j=0; j < nelem; j++) inpoel[j] = new double[3];
 
for(int j=0;j<nelem;j++)
{
fin>>temp>>inpoel[j][0]>>inpoel[j][1]>>inpoel[j][2];
getline(fin,line);
}
getline(fin,line);

//Initializing coordinates of the points array and reading from the file
coord=new double* [npoint];
for (int j=0; j < npoint; j++) coord[j] = new double[2];
for(int j=0;j<npoint;j++)
{
fin>>temp>>coord[j][0]>>coord[j][1];
}
getline(fin,line);
getline(fin,line);

//End of the coordinate section
//--------------------------------------------------------------------
//Initializing the initial condition array and reading the values from file
unknown=new double* [npoint]; 
for (int j=0; j < npoint; j++) unknown[j] = new double[4];

for(int j=0;j<npoint;j++)
{
fin>>temp;
for(int i=0;i<4;i++) fin>>unknown[j][i]; 
getline(fin,line);
}
getline(fin,line);
//End of initial condition

//---------------------------------------------------------------------
// Initalizing the boundary faces vector and reading from the file
bface=new double* [nface]; 
for (int j=0; j < nface; j++) bface[j] = new double[5];
for(int j=0;j<nface;j++)
{
fin>>temp;
for(int i=0;i<5;i++) 
{
	fin>>bface[j][i];
//	getline(fin,line);
}
getline(fin,line);
}

fin.close();
}
