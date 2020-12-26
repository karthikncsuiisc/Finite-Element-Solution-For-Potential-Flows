//This function writes output file containing phi, pressure and velocity vector
//in a VTK format and wrties the x,y,velocity values into a text format for
//comparing the results with exact solution.  

#include "functions.h"


void results()
{
string filename;

ofstream fout;
filename=grid_name+"_results.vtk";
fout.open(filename.c_str(), ios::out | ios::trunc);
fout<<"# vtk DataFile Version 2.0"<<endl;
fout<<"2D Unstructured Grid of Linear Triangles"<<endl;
fout<<"ASCII\n"<<endl;
fout<<"DATASET UNSTRUCTURED_GRID"<<endl;
fout<<"POINTS "<<npoint<<" float"<<endl;
for(int j=0;j<npoint;j++)
{
fout<<coord[j][0]<<" "<<coord[j][1]<<" 0.00"<<endl;	
}
fout<<"\nCells "<<nelem<<" "<<4*nelem<<endl;
for(int j=0;j<nelem;j++)
{
fout<<"3 "<<inpoel[j][0]-1<<" "<<inpoel[j][1]-1<<" "<<inpoel[j][2]-1<<endl;	
}
fout<<"\nCELL_TYPES "<<nelem<<endl;
for(int j=0;j<nelem;j++)
{
fout<<"5"<<endl;	
}

fout<<"\nPOINT_DATA "<<npoint<<endl;
fout<<"SCALARS phi float"<<endl;
fout<<"LOOKUP_TABLE default"<<endl;

for(int i=0;i< npoint;i++)
fout<<unknown[i][0]<<endl;

//--------------------------------------------------------------------
//Calculating the elemental velocities
//-------------------------------------------------------------------

double **vele;   //elemental velocities
vele=new double* [nelem];
for(int i=0;i<nelem;i++) vele[i] =new double[2];
for (int i=0;i<nelem;i++)
for(int j=0; j<2;j++)
vele[i][j]=0;

for(int i=0; i<nelem; i++)
{			
    		int ind0,ind1,ind2;
			double a[3],b[3],c[3],Dtot,tempx,tempy;
            
            ind0=inpoel[i][0]-1;
            ind1=inpoel[i][1]-1;
			ind2=inpoel[i][2]-1;
			a[0]=coord[ind1][1]-coord[ind2][1];
			a[1]=coord[ind2][1]-coord[ind0][1];
			a[2]=coord[ind0][1]-coord[ind1][1];
    		b[0]=-coord[ind1][0]+coord[ind2][0];
			b[1]=-coord[ind2][0]+coord[ind0][0];
			b[2]=-coord[ind0][0]+coord[ind1][0];
			c[0]=coord[ind1][0]*coord[ind2][1]-coord[ind2][0]*coord[ind1][1];
			c[1]=coord[ind2][0]*coord[ind0][1]-coord[ind0][0]*coord[ind2][1];
			c[2]=coord[ind0][0]*coord[ind1][1]-coord[ind1][0]*coord[ind0][1];
			Dtot=(c[0]+c[1]+c[2]);
			
			vele[i][0]=(unknown[ind0][0]*a[0]+unknown[ind1][0]*a[1]+unknown[ind2][0]*a[2])/Dtot;
			vele[i][1]=(unknown[ind0][0]*b[0]+unknown[ind1][0]*b[1]+unknown[ind2][0]*b[2])/Dtot;
			

}

//--------------------------------------------------------------------
//End of elemental velocity calculation
//-------------------------------------------------------------------
//End od elemental velocity calculation


double **vel;
vel=new double* [npoint];
for(int i=0;i<npoint;i++) vel[i] =new double[2];

for (int i=0;i<npoint;i++)
for(int j=0; j<2;j++)
vel[i][j]=0;

//--------------------------------------------------------------------
//Calculating the nodal velocities from elemental values
//-------------------------------------------------------------------
double we;                                        //Weighting function =1/de, Distance of the node from element center
double *wetot;
wetot = new double[npoint];
for(int i=0; i<npoint; i++)
wetot[i]=0;

for(int i=0; i<nelem; i++)
{			
	        double cx,cy,dist;
	        int ind0, ind1,ind2;
            ind0=inpoel[i][0]-1;
            ind1=inpoel[i][1]-1;
			ind2=inpoel[i][2]-1;
			
			cx=coord[ind0][0]+coord[ind1][0]+coord[ind2][0];
			cx=cx/3;
			cy=coord[ind0][1]+coord[ind1][1]+coord[ind2][1];
			cy=cy/3;
			
			dist=sqrt(pow(coord[ind0][0]-cx,2)+pow(coord[ind0][1]-cy,2));
			we=1/(dist);
    		wetot[ind0]=wetot[ind0]+we;
			vel[ind0][0]=vel[ind0][0]+vele[i][0]*we;
			vel[ind0][1]=vel[ind0][1]+vele[i][1]*we;
			
            dist=sqrt(pow(coord[ind1][0]-cx,2)+pow(coord[ind1][1]-cy,2));
			we=1/(dist);
    		wetot[ind1]=wetot[ind1]+we;
			vel[ind1][0]=vel[ind1][0]+vele[i][0]*we;
			vel[ind1][1]=vel[ind1][1]+vele[i][1]*we;

            dist=sqrt(pow(coord[ind2][0]-cx,2)+pow(coord[ind2][1]-cy,2));
			we=1/(dist);
    		wetot[ind2]=wetot[ind2]+we;
			vel[ind2][0]=vel[ind2][0]+vele[i][0]*we;
			vel[ind2][1]=vel[ind2][1]+vele[i][1]*we;
}

for(int i=0;i< npoint;i++)
{
vel[i][0]=vel[i][0]/wetot[i];
vel[i][1]=vel[i][1]/wetot[i];
}
//--------------------------------------------------------------------
//End of nodal velocity calculations
//-------------------------------------------------------------------

//--------------------------------------------------------------------
//Calculating the Pressure inside the domain using Bernoli equation
// 0.5(u.u)+p/rho=const 
//---------------------------------------------------------------------
double Pref,Vref;
double rho,constant;                                  //Reference Gauge pressure and Reference Velocity
Pref=0;                                               //Reference Gauge Pressure is zero
Vref=1;                                               //Reference velocity is free stream velocity 
rho=1;                                                //Assuming density is one.
constant=Pref/rho+0.5*Vref*Vref;
double *pres;
pres=new double[npoint];

fout<<"SCALARS pressure float"<<endl;
fout<<"LOOKUP_TABLE default"<<endl;
for(int i=0;i<npoint;i++)
{
pres[i]=rho*constant-0.5*rho*(vel[i][0]*vel[i][0]+vel[i][1]*vel[i][1]);
fout<<pres[i]<<endl;
}


//fout<<"\nPOINT_DATA "<<npoint<<endl;
fout<<"VECTORS velocity float"<<endl;
//fout<<"LOOKUP_TABLE default"<<endl;

for(int i=0;i< npoint;i++)
{
fout<<vel[i][0]<<" "<<vel[i][1]<<" 0"<<endl;
}
fout.close();
//-------------------------------------------------------------
//End of Writing VKT file files for plotting
//-------------------------------------------------------------

//----------------------------------------------------------------------
//Outputting the x, y, vx, vy, cell size fvalues to a data file
//The mesh size is determined by first calculating the avearage of size 
//of each mesh element and than after taking the average of all elements.
//------------------------------------------------------------------------
   
double havg=0;
for(int i=0; i<nelem; i++)
{			
    		int ind0,ind1,ind2;
			double dy[3],dx[3],hele;

			ind0=inpoel[i][0]-1;
			ind1=inpoel[i][1]-1;
			ind2=inpoel[i][2]-1;
			dy[0]=coord[ind1][1]-coord[ind2][1];
			dy[1]=coord[ind2][1]-coord[ind0][1];
			dy[2]=coord[ind0][1]-coord[ind1][1];
			dx[0]=-coord[ind1][0]+coord[ind2][0];
			dx[1]=-coord[ind2][0]+coord[ind0][0];
			dx[2]=-coord[ind0][0]+coord[ind1][0];

			hele=sqrt(dx[0]*dx[0]+dy[0]*dy[0])
			    +sqrt(dx[1]*dx[1]+dy[1]*dy[1])
			    +sqrt(dx[2]*dx[2]+dy[2]*dy[2]);
            havg=havg+hele/3;             
}
havg=havg/nelem;


filename=grid_name+"_x_y_u_v_h.dat";
fout.open(filename.c_str(), ios::out | ios::trunc);

for(int i=0;i<npoint;i++)
fout<<i<<" "<<unknown[i][0]<<" "<<coord[i][0]<<" "<<coord[i][1]<<" "<<vel[i][0]<<" "<<vel[i][1]<<" "<<havg<<endl;
fout.close();

//------------------------------------------------------------------
//End of writing data file
//-----------------------------------------------------------------
}
