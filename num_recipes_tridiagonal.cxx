// Press et. al.: Numerical Recipes in C:
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void nrerror(string error_text)
/* Numerical Recipes standard error handler */
{
	cout << "Numerical Recipes run-time error...\n";
	cout << error_text << "\n";
	cout << "...now exiting to system...\n";
	exit(1);
}


double pythag(double a, double b)
// calculate (a^2 + b^2)^(1/2) without overflow
{
	double absa,absb;
	double r;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb)
	{
		r = absb/absa;
		return absa*sqrt(1.0+r*r);
	}
	else
	{
		if (absb == 0.0)
			return 0.0;
		else
		{
			r = absa/absb;
			return absb*sqrt(1.0+r*r);
		}
	}
}


// function tqli.c

void tqli(double dx[], double ex[], int n, double **zx)
// dx[]: diagonal elements of the tridiagonal matrix (contains eigenvalues at the end), 
// ex[]: subdiagonal elements (will be destroyed), note e[0] ignored
// zx: input: identity matrix
// rewritten here such that index 0..n-1 is used now
{
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;

	for (i=1;i<n;i++)
		ex[i-1]=ex[i];
	ex[n-1]=0.0;
	for (l=1;l<=n;l++)
	{
		iter=0;
		do
		{
			for (m=l;m<=n-1;m++)
			{
				dd=fabs(dx[m-1])+fabs(dx[m]);
				if ((fabs(ex[m-1])+dd) == dd)
					break;
			}
			if (m != l)
			{
				if (iter++ == 30)
					nrerror("Too many iterations in tqli");
				g=(dx[l]-dx[l-1])/(2.0*ex[l-1]);
				r=pythag(g,1.0);
				g=dx[m-1]-dx[l-1]+ex[l-1]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--)
				{
					f=s*ex[i-1];
					b=c*ex[i-1];
					ex[i]=(r=pythag(f,g));
					if (r == 0.0)
					{
						dx[i] -= p;
						ex[m-1]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=dx[i]-p;
					r=(dx[i-1]-g)*s+2.0*c*b;
					dx[i]=g+(p=s*r);
					g=c*r-b;
					for (k=1;k<=n;k++)
					{
						f=zx[k-1][i];
						zx[k-1][i]=s*zx[k-1][i-1]+c*f;
						zx[k-1][i-1]=c*zx[k-1][i-1]-s*f;
					}
				}
				if (r == 0.0 && i >= l)
					continue;
				dx[l-1] -= p;
				ex[l-1]=g;
				ex[m-1]=0.0;
			}
		}
		while (m != l);
	}
}

