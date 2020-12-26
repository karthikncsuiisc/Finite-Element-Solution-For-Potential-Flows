#ifndef HEADER_H
#define HEADER_H
#include <math.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <cmath>

using namespace std;

extern string grid_name; //Mesh file name
extern int nelem; // Number of elements in the mesh
extern int npoint; // Number of nodes in the mesh
extern int nface; // Number of boundary faces
extern int ndimn; //Dimension of the problem
extern int ntype; //Dimension of the problem
extern double **inpoel; // Node numbers of each element
extern double **coord; // Coordinates of node points
extern double **unknown; // Inital Conditions
extern double **bface; // Boundary elements, corresponding points and 
extern double **lhspo; // System matrix
extern double *rhspo; // Load matrix
extern double tol; //Tolerance for SGS method
extern int Niter; //Maximum number of iterations for SGS method.
extern double havg; //Average cell size

void reading_grid();
void results();
void system_matrices();
void read_input();
void SGS();
#endif







