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
#define DELTA 0.01
#define NUMITR 10


void Usage(char *name);
int main (int argc, char **argv)
{
	int	j=0, K, t, T, ITR,h_idx,l_idx;
	float	*c, *lam, *mu, *sum1, *sum2, *sum3, *sum4;
	float	mn, vr, h_bd, l_bd, ratio=5, denominator, prev_lam;
	double tiny=10e-30;
	float	*O, *OO;	/* observation sequence for RD O[1..T] */
	FILE	*fp;
	
	int N=3;
	lam = vector(1, N);
	mu = vector(1, N);
	c = vector(1, N);
	sum1 = vector(1, N);
	sum2 = vector(1, N);
	sum3 = vector(1, N);
	sum4 = vector(1, N);

	c[1]=0.001;
	c[2]=0.998;
	c[3]=0.001;

	if (argc != 2) {
		printf("Usage error \n");
		printf("Usage: em <file.seq> \n");
		exit (1);
	}
	
	/* read the observed sequence */
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: File %s not found \n", argv[1]);
		exit (1);
	}
	T = countlines(fp);
	rewind(fp);
	
	ReadSequenceshort(fp, T, &O); 
	fclose(fp);

	qsort(O,T, sizeof(float), floatcomp);
	h_idx=MIN(floor(T*(1-c[3])),T-2);
	l_idx=ceil(T*c[1]);
	K=h_idx-l_idx+1;
	OO = vector(1, K);
	j=1;
	for (t =l_idx; t <= h_idx; t++) {
		OO[j]=O[t];
		j++;
	}
	free_vector(O, 1, T);
	mn=mean(OO,K);
	vr=var(OO,K);
	if (T <500){
		lam[2]=0;
		if(T>10){
			lam[2]=mn;
			if(mn>1 & vr >1){
				ratio=vr/mn;		
			}
		}
	}
	else{
		lam[2]=mn;
		lam[1]=lam[2]/2;
		lam[3]=lam[2]*1.5;
		do{
			ITR++;
			sum1 = vector(1, N);
			sum2 = vector(1, N);
			sum3 = vector(1, N);
			sum4 = vector(1, N);
			for(t=1; t<=K; t++){
				denominator=tiny;
				for(j=1; j<=N; j++){
					denominator+=poisson(lam[1]/ratio,OO[t]/ratio,j)*c[j];
				}
				for(j=1; j<=N; j++){
					sum1[j]+=(poisson(lam[1]/ratio,OO[t]/ratio,j)*c[j]/denominator)*OO[t];
					sum2[j]+=(poisson(lam[1]/ratio,OO[t]/ratio,j)*c[j]/denominator)*pow((OO[t]-lam[j]),2);
					sum3[j]+=(poisson(lam[1]/ratio,OO[t]/ratio,j)*c[j]/denominator);
				}
			}
			prev_lam=lam[2];
			for(j=1; j<=N; j++){
				lam[j]=sum1[j]/(sum3[j]+tiny);
				mu[j]=sum2[j]/(sum3[j]+tiny);
			}
			lam[1]=lam[2]/2;
			lam[3]=lam[2]*1.5;
			if(lam[2]>0){
				ratio=mu[2]/lam[2];
			}
			
			for(t=1; t<=K; t++){
				denominator=tiny;
				for(j=1; j<=N; j++){
					denominator+=poisson(lam[1]/ratio,OO[t]/ratio,j)*c[j];
				}
				for(j=1; j<=N; j++){
					sum4[j]+=(poisson(lam[1]/ratio,OO[t]/ratio,j)*c[j]/denominator);
				}
			}
			for(j=1; j<=N; j++){
				c[j]=sum4[j]/K;
			}
		}
		while(ITR<=NUMITR & (prev_lam-lam[2]>DELTA | prev_lam-lam[2]<-DELTA));
	}
	fprintf(stdout, "%.2f\t%.2f\n", lam[2],ratio);
	free_vector(OO, 1, K);
}


int floatcomp(const void* elem1, const void* elem2)
{
    if(*(const float*)elem1 < *(const float*)elem2)
        return -1;
    return *(const float*)elem1 > *(const float*)elem2;
}
