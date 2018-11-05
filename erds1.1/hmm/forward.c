/*
**		Authors: Mingfu Zhu & David B. Goldstein
**		Organization: Center for Human Genome Variation, Duke School of Medicine
**		Modified from Tapas Kanungo, "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite
**		State Models of Language," A. Kornai (editor), Cambridge University Press, 1999.																	    
*/


#include <stdio.h>
#include "hmm.h"

void ForwardWithScale(HMM *phmm, int T, float *O, float *L, double **alpha, 
	double *scale, double *pprob)
/*  pprob is the LOG probability */
{
	int	i, j; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */

	/* 1. Initialization */

	scale[1] = 0.0;	
	for (i = 1; i <= phmm->N; i++) {
		alpha[1][i] = phmm->pi[i] * poisson(L[1],O[1],(i-1));
		scale[1] += alpha[1][i];
	}
	for (i = 1; i <= phmm->N; i++) 
		alpha[1][i] /= scale[1]; 

	/* 2. Induction */

	for (t = 1; t <= T - 1; t++) {
		scale[t+1] = 0.0;
		for (j = 1; j <= phmm->N; j++) {
			sum = 0.0;
			for (i = 1; i <= phmm->N; i++) 
				sum += alpha[t][i] * (phmm->A[i][j]); 

			alpha[t+1][j] = sum * poisson(L[t+1],O[t+1],(j-1));
			scale[t+1] += alpha[t+1][j];
		}
		for (j = 1; j <= phmm->N; j++) 
			alpha[t+1][j] /= scale[t+1]; 
	}

	/* 3. Termination */
	*pprob = 0.0;

	for (t = 1; t <= T; t++)
		*pprob += log(scale[t]);
	
}
