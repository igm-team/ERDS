/*
**		Authors: Mingfu Zhu & David B. Goldstein
**		Organization: Center for Human Genome Variation, Duke School of Medicine
**		Modified from Tapas Kanungo, "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite
**		State Models of Language," A. Kornai (editor), Cambridge University Press, 1999.																	    
*/

#include <stdio.h> 
#include "nrutil.h"
#include "hmm.h"
#include <math.h>
#define DELTA 0.01
#define NUMITR 10

void BaumWelch(HMM *phmm, int T, float *O, float *L, int *H, double **alpha, double **beta,
	double **gamma, int *pniter, double *plogprobinit, double *plogprobfinal)
{
	int	i, j, k;
	int	t, iteration=0;

	double	logprobf, logprobb,  threshold, temp, temp2;
	double	numeratorA, denominatorA;
	double	*c,*lam;

	double ***xi, *scale;
	double delta, deltaprev, logprobprev;
	

	deltaprev = 10e-70;

	xi = AllocXi(T, phmm->N);
	scale = dvector(1, T);

	lam = dvector(1, phmm->N);
	ForwardWithScale(phmm, T, O, L, H, alpha, scale, &logprobf);
	*plogprobinit = logprobf; /* log P(O |intial model) */
	BackwardWithScale(phmm, T, O, L, H, beta, scale, &logprobb);

	ComputeGamma(phmm, T, alpha, beta, gamma);
	ComputeXi(phmm, T, O, L, H, alpha, beta, xi);
	logprobprev = logprobf;

	do  {	
		for (i = 1; i <= phmm->N; i++) {
			phmm->pi[i] = gamma[1][i]*0.999995+0.000001;
		}
		
		/* reestimate transition matrix  and symbol prob in
		   each state */
		for (i = 1; i <= phmm->N-1; i++) { 
			denominatorA = 0.0;
			for (t = 1; t <= T - 1; t++)
				if(O[t]>=0){
					denominatorA += gamma[t][i];
				}
			for (j = 1; j <= phmm->N-1; j++) {
				numeratorA = 0.0;
				for (t = 1; t <= T - 1; t++) {
					if(O[t]>=0){
						numeratorA += xi[t][i][j];
					}
				}
				phmm->A[i][j] = .998*numeratorA/denominatorA;
			}
		}
		iteration++;		

		ForwardWithScale(phmm, T, O, L, H, alpha, scale, &logprobf);
		BackwardWithScale(phmm, T, O, L, H, beta, scale, &logprobb);
		ComputeGamma(phmm, T, alpha, beta, gamma);
		ComputeXi(phmm, T, O, L, H, alpha, beta, xi);

		delta = logprobf - logprobprev; 
		logprobprev = logprobf;
	}
	while (delta>DELTA & iteration<NUMITR); 
	
	*pniter = iteration;
	*plogprobfinal = logprobf; 
	FreeXi(xi, T, phmm->N);
	free_dvector(scale, 1, T);
}

void ComputeGamma(HMM *phmm, int T, double **alpha, double **beta, 
	double **gamma)
{
	int	i, j;
	int	t;
	double	denominator;

	for (t = 1; t <= T; t++) {
		denominator = 0.0;
		for (j = 1; j <= phmm->N; j++) {
			gamma[t][j] = alpha[t][j]*beta[t][j];
			denominator += gamma[t][j];
		}

		for (i = 1; i <= phmm->N; i++) 
			gamma[t][i] = gamma[t][i]/denominator;
	}
}

void ComputeXi(HMM* phmm, int T,float *O, float *L, int *H, double **alpha, double **beta, 
	double ***xi)
{
	int i, j;
	int t;
	double sum;

	for (t = 1; t <= T - 1; t++) {
		sum = 0.0;	
		for (i = 1; i <= phmm->N; i++) 
			for (j = 1; j <= phmm->N; j++) {
				xi[t][i][j] = alpha[t][i]*beta[t+1][j]
					*(phmm->A[i][j])
					* poisson(L[t+1],O[t+1],(j-1))
					* hetp(j,H[t+1],phmm->N);
				sum += xi[t][i][j];
			}

		for (i = 1; i <= phmm->N; i++) 
			for (j = 1; j <= phmm->N; j++) 
				xi[t][i][j]  /= sum;
	}
}

double *** AllocXi(int T, int N)
{
	int t;
	double ***xi;

	xi = (double ***) malloc(T*sizeof(double **));

	xi --;

	for (t = 1; t <= T; t++)
		xi[t] = dmatrix(1, N, 1, N);
	return xi;
}

void FreeXi(double *** xi, int T, int N)
{
	int t;
	for (t = 1; t <= T; t++)
		free_dmatrix(xi[t], 1, N, 1, N);

	xi ++;
	free(xi);

}
