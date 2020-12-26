//Subroutine for assemblying stiffness matrix for potential equation using
//Lagrangian shape function and linear elements

#include "functions.h"

void system_matrices()
{

//---------------------------------------------------------------------
//Initializing the stiffness matrix(lhspo) and Assembling the Stiffness Matrix
//---------------------------------------------------------------------

lhspo = new double* [npoint];
for (int j=0; j < npoint; j++) lhspo[j] = new double[npoint];
for(int i=0; i<npoint; i++)	for(int j=0; j<npoint; j++) lhspo[i][j]=0;	


for(int i=0; i<nelem; i++)
{			
    		int ind0,ind1,ind2;
			double a[3],b[3],c[3],Dtot;

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
lhspo[ind0][ind0]=lhspo[ind0][ind0]+(a[0]*a[0]+b[0]*b[0])/2/Dtot;
lhspo[ind0][ind1]=lhspo[ind0][ind1]+(a[0]*a[1]+b[0]*b[1])/2/Dtot;
lhspo[ind0][ind2]=lhspo[ind0][ind2]+(a[0]*a[2]+b[0]*b[2])/2/Dtot;
lhspo[ind1][ind0]=lhspo[ind1][ind0]+(a[1]*a[0]+b[1]*b[0])/2/Dtot;
lhspo[ind1][ind1]=lhspo[ind1][ind1]+(a[1]*a[1]+b[1]*b[1])/2/Dtot;
lhspo[ind1][ind2]=lhspo[ind1][ind2]+(a[1]*a[2]+b[1]*b[2])/2/Dtot;
lhspo[ind2][ind0]=lhspo[ind2][ind0]+(a[2]*a[0]+b[2]*b[0])/2/Dtot;
lhspo[ind2][ind1]=lhspo[ind2][ind1]+(a[2]*a[1]+b[2]*b[1])/2/Dtot;
lhspo[ind2][ind2]=lhspo[ind2][ind2]+(a[2]*a[2]+b[2]*b[2])/2/Dtot;
}

//-------------------------------------------------------------------
//End of assembling stiffness matrix
//-------------------------------------------------------------------

//-------------------------------------------------------------------
//Initializing and assemblying  load(rhspo) matrix
//-------------------------------------------------------------------

rhspo = new double[npoint];
for(int i=0; i<npoint; i++)	rhspo[i]=0;	
for( int j=0; j<nface; j++)
{
	int ind0,ind1;
	double dx,dy,dl,cond,nx,ny,g;
	int cond_type;
	ind0=bface[j][0]-1;
	ind1=bface[j][1]-1;
	
	dx=coord[ind1][0]-coord[ind0][0];
	dy=coord[ind1][1]-coord[ind0][1];
	dl=dx*dx+dy*dy;
	dl=sqrt(dl);
	cond_type=bface[j][2];
	cond=bface[j][3];
    nx=(coord[ind1][1]-coord[ind0][1])/dl;
    ny=-(coord[ind1][0]-coord[ind0][0]/dl);
	
	g=cond*nx;
	
    rhspo[ind0]=rhspo[ind0]+dl*g/2;	
	rhspo[ind1]=rhspo[ind1]+dl*g/2;   

	
}
//-----------------------------------------
//Constraining \phi at a point
//-----------------------------------------
int dcp;
double Cbig=1e9;
for(int i=0;i<nface;i++)
{
double xbig=-1e10;

	if(bface[i][2]==4)
	{
		if(coord[i][0]>xbig) {xbig=coord[i][0]; dcp=i;} 
    }
}
        lhspo[dcp][dcp]=lhspo[dcp][dcp]*Cbig;
        rhspo[dcp]=lhspo[dcp][dcp]*Cbig*0;


//-------------------------------------------------------------------
// End of calculating the load matrix
//--------------------------------------------------------------------


//---------------------------------------------------------------------
//Modifying the Stifness-Matrix and load Matrix to implement Dirchilet Boundary condition
//--------------------------------------------------------------------- 
}
