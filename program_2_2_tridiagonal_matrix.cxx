/*
  PY 4109: Advanced Computational Physics

  Program 2.2: Eigenvalues and eigenvectors of tridiagonal matrix

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


// main function (to be changed)

int main ()
{
	int nr, nr2;
	
	trimatrix_diag [0] = 1.0;
	trimatrix_diag [1] = 2.0;
	trimatrix_diag [2] = 3.0;
	trimatrix_diag [3] = 4.0;
	
	trimatrix_subdiag [1] = 1.0;
	trimatrix_subdiag [2] = 1.0;
	trimatrix_subdiag [3] = 1.0;
	
	eigenvalues_and_eigenvectors_tridiagonal_matrix (4);
	
	for (nr=0; nr<4; nr++)
	{
		cout << "Eigenvalue " << nr << ": " << trimatrix_eigenvalue[nr] << "\n";
		cout << "Corresponding eigenvector:\n";
		for (nr2=0; nr2<4; nr2++)
			cout << trimatrix_eigenvector [nr][nr2] << "\n";
		cout << "\n";
	}
	
	return 0;
}
