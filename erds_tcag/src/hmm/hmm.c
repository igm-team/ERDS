/*
**		Authors: Mingfu Zhu & David B. Goldstein
**		Organization: Center for Human Genome Variation, Duke School of Medicine
**		Modified from Tapas Kanungo, "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite
**		State Models of Language," A. Kornai (editor), Cambridge University Press, 1999.																	    
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
	int	t,	T;
	HMM  	hmm;
	int	N;
	int	M;
	double 	**alpha; 
	double	**beta;
	double	**gamma;
	float *O;	/* observation sequence for RD O[1..T] */
	float *L;	/* Lam sequence for L[1..T] */
	int	*q;	/* state sequence q[1..T] */
	int	niter;
	double	logprobinit, logprobfinal;
	FILE	*fp;
	float	*logp;
	double **delta;
	int	**psi;
	double 	proba, logproba; 
 
	if (argc != 3) {
		printf("Usage error \n");
		printf("Usage: esthmm <mod.hmm> <file.seq> \n");
		exit (1);
	}
	
 	/* read the hmm model */
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: HMM file not found \n");
				exit (1);
	}
	ReadHMM(fp, &hmm);
	fclose(fp);

	/* read the observed sequence */
	fp = fopen(argv[2], "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: File %s not found \n", argv[2]);
		exit (1);
	}
	T = countlines(fp);
	rewind(fp);
	ReadSequence(fp, T, &O, &L); 
	fclose(fp);
	
	/* allocate memory */
	alpha = dmatrix(1, T, 1, hmm.N);
	beta = dmatrix(1, T, 1, hmm.N);
	gamma = dmatrix(1, T, 1, hmm.N);

	/* call Baum Welch */
	BaumWelch(&hmm, T, O, L, alpha, beta, gamma, &niter, 
		&logprobinit, &logprobfinal);
		
	/* call Viterbi */
	q = ivector(1,T);
	delta = dmatrix(1, T, 1, hmm.N);
	psi = imatrix(1, T, 1, hmm.N);

	ViterbiLog(&hmm, T, O, L, delta, psi, q, &logproba); 

	/* find zero copy-number difference state */
	int normal= 2;
	/* write Viterbi and log odds sequences */
	logp = vector(1,T);	
	logp[1]=hmm.pi[q[1]] - hmm.pi[normal] +1*(hmm.A[q[1]][q[2]] - hmm.A[normal][q[2]])+ log(poisson(L[1],O[1],q[1]-1)) - log(poisson(L[1],O[1],normal-1));

	for (t = 2; t <= T; t++)
		logp[t]=1*(hmm.A[q[t - 1]][q[t]] - hmm.A[q[t - 1]][normal]) + log(poisson(L[t],O[t],q[t]-1)) - log(poisson(L[t],O[t],normal-1));

	PrintSequence(stdout, T, q, logp);

	/* free memory */
	free_dmatrix(alpha, 1, T, 1, hmm.N);
	free_dmatrix(beta, 1, T, 1, hmm.N);
	free_dmatrix(gamma, 1, T, 1, hmm.N);
	free_ivector(q, 1, T);
	free_vector(O, 1, T);
	free_vector(L, 1, T);
	free_vector(logp, 1, T);
	free_imatrix(psi, 1, T, 1, hmm.N);
	free_dmatrix(delta, 1, T, 1, hmm.N);
	FreeHMM(&hmm);
}
