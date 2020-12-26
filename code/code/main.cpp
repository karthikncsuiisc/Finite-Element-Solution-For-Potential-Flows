// A C++ program for solving INCOMPRESSIBLE POTENTIAL FLOW equation. 
// The Program is arranged into different subroutines. All the functions
// are declared in the header file "functions.h". A make file is aviable
// to compile all the subroutines.
//

#include "functions.h"

string grid_name;       	 // Mesh file name
int nelem;              	 // Number of elements in the mesh
int npoint;                  // Number of nodes in the mesh
int nface;              	 // Number of boundary faces
int ndimn;              	 // Dimension of the problem
int ntype;              	 // Dimension of the problem
double **inpoel;		     // Node numbers of each element
double **coord; 	   		 // Coordinates of node points
double **unknown; 			 // Inital Conditions
double **bface; 			 // Boundary elements
double **lhspo;              // System matrix
double  *rhspo; 			 // Load matrix
double tol; 				 // Tolerance for SGS method
int Niter; 					 // Maximum number of iterations for SGS method.

int main()
{
read_input();
reading_grid();
system_matrices();
SGS();
results();
return 0;
}





