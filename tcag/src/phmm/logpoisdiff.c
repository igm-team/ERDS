/*
**		Authors: Mingfu Zhu & David B. Goldstein
**		Organization: Center for Human Genome Variation, Duke School of Medicine
*/

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm.h"
#include <sys/types.h>
#include <unistd.h>

void Usage(char *name);
int main (int argc, char **argv)
{
	float	lam;
	float	n;
	int a, b;
	
	lam=atof(argv[1]);
	n=atof(argv[2]);
	a=atoi(argv[3]);
	b=atoi(argv[4]);
	double lpd=logpoisdiff(lam,n,a,b);
	fprintf(stdout, "%.2f", lpd);

}

double logpoisdiff(double lam, double n, int m, int k)
{
	double lpd, lama, lamb;
	double tiny=10e-30;
	lama=lam*m+tiny;
	lamb=lam*k+tiny;
	if (n==-1)
        return 1;
	else if (n==0)
        return lamb-lama;
	else{
		lpd=n*(log(lama)-log(lamb))-lama+lamb;
		return lpd;
	}
}
