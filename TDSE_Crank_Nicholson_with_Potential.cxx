/*
  PY 4109: Advanced Computational Physics

  Program 3.1: Solution of time-dependent free Schroedinger equation (Crank-Nicholson method)

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

const double pi = M_PI;

// Constants (max. size of matrix)
const int max_matrix = 10000;

// Input (will not be distroyed)
double complex trimatrix_diag [max_matrix];
// we assume the value 1 on the subdiagonal here
double complex vector_b [max_matrix];

// Output
double complex solution_psi [max_matrix];

// only used during solution_tridiagonal_matrixequation
double complex value_d [max_matrix], value_bx [max_matrix];


// Solution of a Matrix Equation with a Tridiagonal Matrix (do not change)

void solution_tridiagonal_matrixequation (int nn)
{
	int j;
	
	// First forward loop
	value_d [0] = trimatrix_diag [0];
	value_bx [0] = vector_b [0];
	for (j=1; j<nn; j=j+1)
	{
		value_d [j] = trimatrix_diag[j] - 1.0/value_d[j-1];
		value_bx [j] = vector_b [j] - value_bx [j-1]/value_d[j-1];
	}
	
	// Second backward loop
	solution_psi [nn-1] = value_bx[nn-1]/value_d[nn-1];
	for (j=nn-2; j>=0; j=j-1)
	{
		solution_psi [j] = (value_bx[j]-solution_psi [j+1])/value_d [j];
	}
}


double cabs2 (double complex z)
{
	return creal(z) * creal(z) + cimag(z) * cimag(z);
}


// constants
const double mass = 1.0;
const double hbar = 1.0;
const double length = 10.0;
double final_time = 0.4;
double const potential_v0 = 0.0;

double potential (double x)
{
	// particle in a box with finite well
	//if ((x >= 3.0*length/8.0) && (x <= 5.0 * length/8.0))
		//return -potential_v0*sin((4.0*pi/length)*(x - 3.0*length/8));
	//else
		//return potential_v0;
	if (x >= 5.0 * length/8.0)
		return potential_v0;
	else
		return 0.0;
}

double x0 = 5.0;  // center of initial wave packet
double p0 = 20.0;  // momentum of initial wave packet
double sigma = 0.2; // width parameter of initial wave packet 

const int number_x_steps = 1000;
const int number_t_steps = 1000;
double complex wavefunction [number_t_steps][number_x_steps];


// main function (to be changed)

int main ()
{
	int nt, nx;
	double delta_x, delta_t, x;
	double complex beta, betax[number_x_steps], alpha;
	double data_x [number_x_steps], data_y1 [number_x_steps], data_y2 [number_x_steps];
	
	delta_x = length/(number_x_steps + 1.0);
	delta_t = final_time/(1.0*number_t_steps);

	alpha = I * hbar * delta_t/(4.0*mass*delta_x*delta_x);
	//beta = -(1+potential(nx))/alpha - 2.0;
	//betax = -(1+potential(nx))/alpha + 2.0;

	// set trimatrix for Schroedinger evolution
	for (nx=0; nx<number_x_steps; nx=nx+1){
		//BETA IS NOW X DEPENDANT
		betax[nx] = -(1+potential(nx))/alpha + 2.0;
		beta = -(1+potential(nx))/alpha - 2.0;
		trimatrix_diag [nx] = beta;
	}

	// set initial wave function
	for (nx=0; nx < number_x_steps; nx=nx+1)
	{
		x = (nx+1) * delta_x;
		wavefunction [0][nx] = 1.0/sqrt(sqrt(pi)*sigma)*cexp(-(x-x0)*(x-x0)/(2.0*sigma*sigma)+ I * p0/hbar * x);
	}
	
	// time steps
	for (nt=0; nt < number_t_steps-1; nt=nt+1)
	{
		// calculate vector b
		vector_b [0] = betax[0] * wavefunction [nt][0] - wavefunction [nt][1];
		for (nx=1; nx < number_x_steps-1; nx=nx+1)
		{
			vector_b [nx] = -wavefunction [nt][nx-1] + betax[nx] * wavefunction [nt][nx] - wavefunction [nt][nx+1];
		}
		vector_b [number_x_steps-1] = -wavefunction [nt][number_x_steps-2] + betax[number_x_steps-1] * wavefunction [nt][number_x_steps-1];

		// solve matrixequation
		solution_tridiagonal_matrixequation (number_x_steps);
		
		// copy solution in wavefunction
		for (nx=0; nx<number_x_steps; nx=nx+1)
			wavefunction [nt+1][nx] = solution_psi [nx];
	}
	
	// set plotting data
	for (nx=0; nx < number_x_steps; nx=nx+1)
	{
		data_x [nx] = (nx+1) * delta_x;
		data_y1 [nx] = cabs2 (wavefunction [0][nx]);
		data_y2 [nx] = cabs2 (wavefunction [number_t_steps-1][nx]);
	}
	
	gnuplot_two_functions ("Free evolution", "linespoints", "x", "|psi|^2",
		data_x, data_y1, number_x_steps, "initial time",
		data_x, data_y2, number_x_steps, "final time");
		
	return 0;
}
