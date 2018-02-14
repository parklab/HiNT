/*
Written by Ruibin Xi
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>
#include "../lib/read.h"
#include <getopt.h>
#include <math.h>
#include <unistd.h>

int main(int argc, char **argv)
{	FILE **infiles=NULL;
	int i, nfiles;
	double ** data;
	int *ncols,*nrows, *nrows_total,k,binTotal;
	double varMeanRatio;
	double sum_obs, sum_exp, sum_expsq, theta, obs_exp_diff;

	if(argc<2){
		fprintf(stderr,"%s <file1> <file2> ...\n",argv[0]);
		exit(1);
		}

	infiles = (FILE **) malloc(sizeof(FILE *)*argc);
	nfiles = argc - 1;
	for(i=0;i<argc-1;i++){
		infiles[i] = fopen(argv[i+1],"r");
		if(infiles[i]==NULL){
			fprintf(stderr,"fopen: no such file or diectory: %s\n", argv[i+1]);
			exit(1);
			}
		}
	ncols = (int *) malloc(sizeof(int)*(nfiles+1));
	nrows = (int *) malloc(sizeof(int)*(nfiles+1));
	nrows_total = (int *) malloc(sizeof(int)*(nfiles+1));
	data = (double **) malloc(sizeof(double *)*(nfiles+1));

	for(i=0;i<nfiles;i++){
		nrows_total[i] = 0;
		}

	varMeanRatio = 0.0;
	binTotal = 0;
	sum_obs = 0.0;
	sum_exp = 0.0;
	sum_expsq = 0.0;
	theta = 0.0;
	for(i=0;i<nfiles;i++){
		data[i] = read_table(infiles[i], &nrows[i],&ncols[i],-1,0);
		fprintf(stderr,"read %d lines\n",nrows[i]);
		if(ncols[i]<4){
			fprintf(stderr,"Error: %s is expected to have at least 4 columns, but has %d columns\n",argv[i+1],ncols[i]);
			exit(1);
			}
		if(data[i]!=NULL && nrows[i]>0){
			for(k=0;k<nrows[i];k++){
				if(data[i][k*ncols[i]+3]>0.0){
					varMeanRatio += data[i][k*ncols[i]+4]/data[i][k*ncols[i]+3];
					binTotal++;
					//sum_obs += data[i][k*ncols[i]+2];
					//sum_exp += data[i][k*ncols[i]+3];
					//sum_expsq += data[i][k*ncols[i]+3]*data[i][k*ncols[i]+3];
					obs_exp_diff = (data[i][k*ncols[i]+2]-data[i][k*ncols[i]+3]);
					theta += (obs_exp_diff*obs_exp_diff - data[i][k*ncols[i]+3])/(data[i][k*ncols[i]+3]*data[i][k*ncols[i]+3]);
					}
				
				}
			free(data[i]); data[i] = NULL;
			}

		}

	theta = theta/binTotal;
	varMeanRatio = varMeanRatio/binTotal;
	fprintf(stdout,"theta\tvarMeanRatio\n%g\t%g\n",theta,varMeanRatio);

	return 0;
}
