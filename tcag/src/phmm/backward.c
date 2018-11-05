/*
**		Authors: Mingfu Zhu & David B. Goldstein
**		Organization: Center for Human Genome Variation, Duke School of Medicine
**		Modified from Tapas Kanungo, "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite
**		State Models of Language," A. Kornai (editor), Cambridge University Press, 1999.																	    
*/

#include <stdio.h>
#include "hmm.h"

void BackwardWithScale(HMM *phmm, int T, float *O, float *L, int *H, double **beta, 
	double *scale, double *pprob)
{
        int     i, j;   /* state indices */
        int     t;      /* time index */
		double sum;
 
        /* 1. Initialization */
 
        for (i = 1; i <= phmm->N; i++)
			beta[T][i] = 1.0/scale[T]; 
 
        /* 2. Induction */
 
        for (t = T - 1; t >= 1; t--) {
			for (i = 1; i <= phmm->N; i++) {
				sum = 0.0;
				for (j = 1; j <= phmm->N; j++)
					sum += phmm->A[i][j] * 
					poisson(L[t+1],O[t+1],(j-1))*
					hetp(j,H[t+1],phmm->N)* 				beta[t+1][j];
					beta[t][i] = sum/scale[t];
			}
        }
}
