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
void ReadSequence(FILE *fp, int T, float **pO, float **pL, int **pH);
void ReadSequenceshort(FILE *fp, int T, float **pO);
void PrintSequence(FILE *fp, int T, int *O, float *logp);
 
void ForwardWithScale(HMM *phmm, int T, float *O, float *L, int *H, double **alpha,
        double *scale, double *pprob);
void BackwardWithScale(HMM *phmm, int T, float *O, float *L, int *H, double **beta,
        double *scale, double *pprob);
void BaumWelch(HMM *phmm, int T, float *O, float *L, int *H, double **alpha, double **beta,
        double **gamma, int *niter, double *plogprobinit, double *plogprobfinal);

double poisson(double lam, double n, int m);
double logpoisdiff(double lam, double n, int m, int k);
double hetp(int n, int h, int N);
double gammln(double xx);
float mean(float a[], int n);
float var(float a[], int n);
int floatcomp(const void* elem1, const void* elem2);
double *** AllocXi(int T, int N);
void FreeXi(double *** xi, int T, int N);
void ComputeGamma(HMM *phmm, int T, double **alpha, double **beta,
        double **gamma);
void ComputeXi(HMM* phmm, int T, float *O, float *L, int *H, double **alpha, double **beta,
        double ***xi);
void ViterbiLog(HMM *phmm, int T, float *O, float *L, int *H, double **delta, int **psi,
        int *q, double *pprob, float *logp);

double gln(double z);
 
#define MAX(x,y)        ((x) > (y) ? (x) : (y))
#define MIN(x,y)        ((x) < (y) ? (x) : (y))
 

