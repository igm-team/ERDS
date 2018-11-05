/*
**		Authors: Mingfu Zhu & David B. Goldstein
**		Organization: Center for Human Genome Variation, Duke School of Medicine
**		Modified from Tapas Kanungo, "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite
**		State Models of Language," A. Kornai (editor), Cambridge University Press, 1999.																	    
*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct {
	int N;		/* number of states;  Q={1,2,...,N} */
	double	**A;	/* A[1..N][1..N]. a[i][j] is the transition prob
			   of going from state i at time t to state j
			   at time t+1 */
	double	*pi;	/* pi[1..N] pi[i] is the initial state distribution. */
} HMM;
void ReadHMM(FILE *fp, HMM *phmm);
void FreeHMM(HMM *phmm);

int64_t countlines(FILE *fp);
void ReadSequence(FILE *fp, int T, float **pO, float **pL);
void PrintSequence(FILE *fp, int T, int *O, float *logp);
 
void ForwardWithScale(HMM *phmm, int T, float *O, float *L, double **alpha,
        double *scale, double *pprob);
void BackwardWithScale(HMM *phmm, int T, float *O, float *L, double **beta,
        double *scale, double *pprob);
void BaumWelch(HMM *phmm, int T, float *O, float *L, double **alpha, double **beta,
        double **gamma, int *niter, double *plogprobinit, double *plogprobfinal);

double poisson(double lam, double n, int m);
double gammln(double xx);
float mean(float a[], int n);
float std_dev(float a[], int n);
double *** AllocXi(int T, int N);
void FreeXi(double *** xi, int T, int N);
void ComputeGamma(HMM *phmm, int T, double **alpha, double **beta,
        double **gamma);
void ComputeXi(HMM* phmm, int T, float *O, float *L, double **alpha, double **beta,
        double ***xi);
void ViterbiLog(HMM *phmm, int T, float *O, float *L, double **delta, int **psi,
        int *q, double *pprob);

double gln(double z);
 
#define MAX(x,y)        ((x) > (y) ? (x) : (y))
#define MIN(x,y)        ((x) < (y) ? (x) : (y))
