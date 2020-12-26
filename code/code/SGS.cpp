//Subroutine to solve the system of equations Kx=b using Symmetric Gauss Sidel Method

#include "functions.h"

void SGS()
{
double res[npoint],dx2[npoint],dx1[npoint],dx[npoint];
double resnorm0,ratio;
double phi[npoint];

for(int i=0;i<npoint;i++) phi[i]=0.0; 


for(int i=0;i<npoint;i++){
res[i]=rhspo[i];
for(int j=0;j<npoint;j++)
res[i]=res[i]-lhspo[i][j]*phi[j];
}
//----------------------------------------------------------------
//SGS Algorithm Start
//----------------------------------------------------------------

//----------Start of forward fauss sidel----------
for(int iter=0;iter<Niter;iter++) 
{
cout<<"Iteration Number: "<<iter<<endl;
for(int i=0;i<npoint;i++){ 
dx2[i]=res[i];
for(int j=0;j<npoint;j++) if(j<i) dx2[i]=dx2[i]-lhspo[i][j]*dx2[j];
dx2[i]=dx2[i]/lhspo[i][i];
dx1[i]=dx2[i]*lhspo[i][i];
}
//----------end of forward gauss sidel----------

//--------Start of backward gauss sidel--------
for(int i=npoint-1;i>=0;i--)
{ 
dx[i]=dx1[i];
for(int j=0;j<npoint;j++) if(j>i) dx[i]=dx[i]-lhspo[i][j]*dx[j];
dx[i]=dx[i]/lhspo[i][i];
phi[i]=phi[i]+dx[i];
}
//--------End of backward gauss sidel--------

// Residue Calculation
double resnorm=0;
for(int i=0;i<npoint;i++)
{
res[i]=rhspo[i];
for(int j=0;j<npoint;j++)
res[i]=res[i]-lhspo[i][j]*phi[j];
resnorm=resnorm+res[i]*res[i];
}
resnorm=sqrt(resnorm);
if(iter==0) resnorm0=resnorm+1e-10;
ratio=resnorm/resnorm0;
cout<<"Absolute Error:"<<resnorm<<endl;
cout<<"Relative Error:"<<ratio<<endl;
cout<<"-----------------------------------------------------"<<endl;
if(ratio<tol) break;
}
if(ratio<tol)
cout<<"The solver converged"<<endl;
else
cout<<"Maximum number of iterations reached, N="<<Niter<<endl;

for(int i=0;i<npoint;i++) unknown[i][0]=phi[i];

}
