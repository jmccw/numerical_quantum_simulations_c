// Q 2.1
// Jordan Walsh | 120387836

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

using namespace std;

int n; //Matrix size - user defined
const double pi = M_PI;
const int max_matrix = 10000;

// Constants (max. size of matrix)

// Input (will de distroyed during calculation of eigenvalues and eigenvectors)
double trimatrix_diag [max_matrix];
double trimatrix_subdiag [max_matrix];

// Output
double trimatrix_eigenvalue [max_matrix];
double trimatrix_eigenvector [max_matrix][max_matrix];

// only used during calculation of eigenvalues and eigenvectors
double trimatrix_result[max_matrix][max_matrix];
double* pointer_matrix [max_matrix];

// as defined in example program
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
 
void print_matrix(double M[][max_matrix]){ //Problem Set 1
	for(int i = 0; i < n; i++){	
		for(int j = 0; j < n; j++){ 
			cout << M[i][j];			//print element
			if(j != n-1) cout << ", ";	// (decoration*)
			else if (j == n-1 && i == n-1) cout << "]";
		}
		cout << "\n ";
	}
}

int main(int argc, char **argv)
{
	cout << "This program calculates and prints the eigenvalues and eigenvectors \
of a tridiagonal N x N matrix with subdiagonal elements = 1. \nEnter a positive \
integer N (max size n=100):\n\nN > ";
	cin >> n;					//user assignment of n.
	
	double M[n][max_matrix]; 	//declaration of print matrix "M"
	double sub_diagonals = 1;	//as specifed by question
	
	cout << "\nDeclare diagonal matrix elements of M, where d(i) represents \
the ith diagonal component M(i, i) [\"ith component of the ith row\" i = 0 -> N-1]:\n\n";
	
	//USER INPUT
	for(int i = 0; i < n; i++){
		cout << "d(" << i << ") > "; 
		cin >> trimatrix_diag[i];				//DIAGONAL ELEMENTS
		if(i != 0) trimatrix_subdiag[i] = sub_diagonals; // = 1, as specified by question.
			//in the example program trimatrix_subdiag[0] is undefined.
	}
	
	//DEFINING "PRINT MATRIX" - The purpose of this is to show the user their input (not a general purpose solution)
	for(int i = 0; i < n; i++){				
		M[i][i] = trimatrix_diag[i];	// Diagonals
		for(int j = 0; j < n; j++){ 								//
			if(j == i-1 || j == i+1) M[i][j] = sub_diagonals;		// Upper and Lower elements
			else if (j < i-1 || j > i+1) M[i][j] = 0;				//
		}
	}
	
	//SHOW USER INPUT
	cout << "\nYou have defined:\nM:\n["; 	
	print_matrix(M);							//see print_matrix() - formatting 
	cout << "\n";
	
	//CALCULATE + PRINT EIGENVALUES/EIGENVECTORS - verbatim example program.
	eigenvalues_and_eigenvectors_tridiagonal_matrix(n);
	
	for (int nr=0; nr<n; nr++)
	{
		cout << "Eigenvalue " << nr << ": " << trimatrix_eigenvalue[nr] << "\n";
		cout << "Corresponding eigenvector:\n";
		for (int nr2=0; nr2<n; nr2++)
			cout << trimatrix_eigenvector [nr][nr2] << "\n";
		cout << "\n";
	}


	return 0;
}

