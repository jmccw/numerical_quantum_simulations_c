/*
  PY 4109: Advanced Computational Physics

  Program 2.2: Energy Eigenvalues of kinetic term and a potential

*/


// Header

#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
using namespace std;
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <complex.h>
#define complex _Complex
#include "gnuplot.cxx"

// Function from Press et. al.: Numerical Recipes in C
#include "num_recipes_tridiagonal.cxx"

const double pi = M_PI;

// Constants (max. size of matrix)
const int max_matrix = 10000;

// Input (will de distroyed during calculation of eigenvalues and eigenvectors)
double trimatrix_diag [max_matrix];
double trimatrix_subdiag [max_matrix];

// Output
double trimatrix_eigenvalue [max_matrix];
double trimatrix_eigenvector [max_matrix][max_matrix];

// only used during calculation of eigenvalues and eigenvectors
double trimatrix_result[max_matrix][max_matrix];
double* pointer_matrix [max_matrix];


// calculates eigenvalues and eigenvectors based on tqli funtion (do not change)

void eigenvalues_and_eigenvectors_tridiagonal_matrix (int nn)
{
	int n, m;
	for (n=0; n<nn; n++)
	{
		for (m=0; m<nn; m++)
		{
			if (n==m)
				trimatrix_result [n][m] = 1.0;
			else
				trimatrix_result [n][m] = 0.0;
		}
	}
	for (n=0; n<nn; n++)
		pointer_matrix [n] = &(trimatrix_result[n][0]);
	tqli (trimatrix_diag, trimatrix_subdiag, nn, pointer_matrix);
	for (n=0; n<nn; n++)
	{
		trimatrix_eigenvalue [n] = trimatrix_diag[n];
		for (m=0; m<nn; m++)
			trimatrix_eigenvector [n][m] = trimatrix_result[m][n];
	}
}

// constants
const double mass = 1.0;
const double hbar = 1.0;
const double length = 1.0;
const double potential_v0 = 400.0;

// potential (to be changed)
double potential (double x)
{
	// particle in a box with finite well
	if ((x >= 3.0*length/8.0) && (x <= 5.0 * length/8.0))
		return -potential_v0*sin((4.0*pi/length)*(x - 3.0*length/8));
	else
		return 0.0;
}


// main function (to be changed)
int main ()
{
	int nr;
	double delta_x, beta, x[max_matrix];
	int number_steps = 400;
	
	delta_x = length/(number_steps + 1.0);
	beta = 1.0/(delta_x * delta_x);
	
	for (nr=0; nr < number_steps; nr++){
		trimatrix_diag [nr] = 2.0/(delta_x*delta_x) + 2.0*mass/(hbar*hbar) * potential((nr+1)*delta_x);
		x[nr] = nr*delta_x;
	}
	for (nr=1; nr < number_steps; nr++)
		trimatrix_subdiag [nr] = -beta;
		
	eigenvalues_and_eigenvectors_tridiagonal_matrix(number_steps);
	
	cout << "Ground State Energy E_0: " << trimatrix_eigenvalue[number_steps-1]*hbar*hbar/(2.0*mass);
	cout << "\nExcited State Energy E_1: " << trimatrix_eigenvalue[number_steps-2]*hbar*hbar/(2.0*mass) << "\n";
	
	gnuplot_two_functions ("Numerical Eigenfunctions", "linespoints", "x", "Eigenfunction |E_n> (x)",
			x, trimatrix_eigenvector[number_steps-1], nr, "Ground State (n=0)", x, trimatrix_eigenvector[number_steps-2], nr, "First Excited State (n=1)");
	
	// Note: Energy eigenvalues = matrix eigenvalues *hbar^2/(2m)
	//for (nr=number_steps-1; nr>=0; nr--)
	//{
		//cout << "Energy-Eigenvalue " << number_steps-nr << ": " << trimatrix_eigenvalue[nr]*hbar*hbar/(2.0*mass) << "\n";
	//}
	
	return 0;
}
