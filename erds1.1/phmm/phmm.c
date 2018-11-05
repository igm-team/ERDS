/*
**		Authors: Mingfu Zhu & David B. Goldstein
**		Organization: Center for Human Genome Variation, Duke School of Medicine
**		Modified from Tapas Kanungo, "UMDphmm: Hidden Markov Model Toolkit," in "Extended Finite
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
	HMM  	phmm;
	int	N;
	int	M;
	double 	**alpha; 
	double	**beta;
	double	**gamma;
	float *O;	/* observation sequence for RD O[1..T] */
	int	*H;	/* observation sequence for Hets H[1..T] */
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
		printf("Usage: phmm <mod.hmm> <file.seq> \n");
		exit (1);
	}
	
 	/* read the phmm model */
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: phmm file not found \n");
				exit (1);
	}
	ReadHMM(fp, &phmm);
	fclose(fp);

	/* read the observed sequence */
	fp = fopen(argv[2], "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: File %s not found \n", argv[2]);
		exit (1);
	}
	T = countlines(fp);
	rewind(fp);
	ReadSequence(fp, T, &O, &L, &H);
	fclose(fp);
	/* allocate memory */
	alpha = dmatrix(1, T, 1, phmm.N);
	beta = dmatrix(1, T, 1, phmm.N);
	gamma = dmatrix(1, T, 1, phmm.N);

	/* call Baum Welch */
	BaumWelch(&phmm, T, O, L, H, alpha, beta, gamma, &niter, 
		&logprobinit, &logprobfinal);		
	/* call Viterbi */
	q = ivector(1,T);
	delta = dmatrix(1, T, 1, phmm.N);
	psi = imatrix(1, T, 1, phmm.N);
	logp=vector(1,T);

	ViterbiLog(&phmm, T, O, L, H, delta, psi, q, &logproba, logp); 	
	PrintSequence(stdout, T, q, logp);
	
	/* free memory */
	free_dmatrix(alpha, 1, T, 1, phmm.N);
	free_dmatrix(beta, 1, T, 1, phmm.N);
	free_dmatrix(gamma, 1, T, 1, phmm.N);
	free_ivector(q, 1, T);
	free_vector(O, 1, T);
	free_vector(L, 1, T);
	free_ivector(H, 1, T);
	free_vector(logp, 1, T);
	free_imatrix(psi, 1, T, 1, phmm.N);
	free_dmatrix(delta, 1, T, 1, phmm.N);
	FreeHMM(&phmm);
}
