/*
**		Authors: Mingfu Zhu & David B. Goldstein
**		Organization: Center for Human Genome Variation, Duke School of Medicine
**		Modified from Tapas Kanungo, "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite
**		State Models of Language," A. Kornai (editor), Cambridge University Press, 1999.																	    
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"
#include "hmm.h"


void ReadHMM(FILE *fp, HMM *phmm)
{
	int i, j, k;
	fscanf(fp, "N= %d\n", &(phmm->N)); 
	fscanf(fp, "A:\n");
	phmm->A = (double **) dmatrix(1, phmm->N, 1, phmm->N);
	for (i = 1; i <= phmm->N; i++) { 
		for (j = 1; j <= phmm->N; j++) {
			fscanf(fp, "%lf", &(phmm->A[i][j])); 
		}
		fscanf(fp,"\n");
	}

	fscanf(fp, "pi:\n");
	phmm->pi = (double *) dvector(1, phmm->N);
	for (i = 1; i <= phmm->N; i++) 
		fscanf(fp, "%lf", &(phmm->pi[i])); 

}

void FreeHMM(HMM *phmm)
{
	free_dmatrix(phmm->A, 1, phmm->N, 1, phmm->N);
	free_dvector(phmm->pi, 1, phmm->N);
}

int64_t countlines(FILE *fp){
	int64_t lines = 0;
	char ch, prev = '\n';

	while((ch = fgetc(fp)) != EOF){
		if (ch == '\n'){
			lines++;
		}
		prev = ch;
	}
	if (prev != '\n')
		lines++;
	return(lines);
}


void ReadSequence(FILE *fp, int T, float **pO, float **pL, int **pH)
{
	int *H;
	float *O, *L;
	int i;

	O = vector(1, T);
	L = vector(1, T);
	H = ivector(1, T);
	for (i=1; i <= T; i++)
		fscanf(fp,"%*d\t%*d\t%f\t%f\t%d\n", &O[i], &L[i], &H[i]);
	*pO = O;
	*pL = L;
	*pH = H;
}

void ReadSequenceshort(FILE *fp, int T, float **pO)
{
	float *O;
	int i;

	O = vector(1, T);
	for (i=1; i <= T; i++)
		fscanf(fp,"%f\n", &O[i]);
	*pO = O;
}

 
void PrintSequence(FILE *fp, int T, int *O, float *logp)
{
	int i;
	for (i=1; i <= T; i++) 
		fprintf(fp,"%d\t%f\n", O[i], logp[i]); 
}

double poisson(double lam, double n, int m)
{
	double p;
	lam*=m;
	if (n==-1)
        return 1;
	else if (n==0)
        return exp(-lam);
	else{
		p=exp(n*log(lam)-lam-gammln(n));
		if(p<1e-300)
			return 1e-300;
		else
			/*return (1-pow(10,-4))*p+pow(10,-6);*/
			return p;
	}
}


double gammln(double z)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
	24.01409824083091,-1.231739572450155,
	0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=z;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser);
}

float mean(float a[], int n) {
    if(n == 0)
        return 0.0;
    double sum = 0;
	int i;
    for( i= 1; i <= n; i++) {
       sum += a[i];
    }
    return sum / n;
}

float var(float a[], int n) {
    if(n <=2)
        return 0.0;
	int i;
	float mn;
    float sq_sum = 0;
    for( i= 1; i <= n; i++) {
       sq_sum += a[i] * a[i];
    }
    mn = mean(a, n);
    return sq_sum /(n-1) - mn * mn*n/(n-1);
}

double hetp(int n, int h, int N)
{
	if (h==0)
        return 1.0/N;
    if (h<10){
		if(n==1)
			return  pow(10,-4*h);
		if(n==2)
			return  pow(10,-2*h);
		else
			return 1.0/(N-2);
	}
    else{
		if(n==1|n==2)
			return  10e-30;
		else
			return 1.0/(N-2);
	}
}


