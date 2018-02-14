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
	int i,j, nfiles;
	double ** data;
	int *ncols,*nrows, *nrows_total,k,s;
	int mline_rd = 10000;
	int flag = 0;

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

	flag = 0; k = 0;
	while(flag!=1){
		for(i=0;i<nfiles;i++){
			data[i] = read_table(infiles[i], &nrows[i],&ncols[i],mline_rd,0);
			nrows_total[i] += nrows[i];
			if(nrows_total[i]!=nrows_total[0]){
				fprintf(stderr,"Erorr, %s has %d rows of numeric data, but %s has %d rows of numeric data\n",argv[1],nrows_total[0],argv[i+1],nrows_total[i]);
				exit(1);
				}
			if(ncols[i]<4) {
				fprintf(stderr,"Error, the files must have at least 4 columns. %s has only %d columns\n",argv[i+1],ncols[i]);
				exit(1);
				}
			if(nrows[i] < mline_rd) {flag = 1;}
			}

		for(j=0;j<nrows[0];j++){
			s = j*ncols[0];
			fprintf(stdout,"%.0f\t%.0f\t%g\t%g",data[0][s],data[0][s+1],data[0][s+2],data[0][s+3]);
			for(i=1;i<nfiles;i++){
				if(data[i][j*ncols[i]] - data[0][j*ncols[0]]!=0.0 || data[i][j*ncols[i]+1] - data[0][j*ncols[0]+1]!=0.0){
					fprintf(stderr,"Error: the start/end positions of the %dth row of %s and %s are different\n", k,argv[i+1],argv[0]);
					exit(1);
					}

				fprintf(stdout,"\t%g\t%g",data[i][j*ncols[i]+2],data[i][j*ncols[i]+3]);
				}
			fprintf(stdout,"\n");
			}
		
		/*free the space allocated to data[i]*/
		for(i=0;i<nfiles;i++){
			if(data[i]!=NULL) {free(data[i]);data[i]=NULL;}
			ncols[i] = 0;
			nrows[i] = 0;
			}
		}


	return 0;
}
