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
#include "ITree.h"


typedef struct args_countRead_t{
        FILE *readFile1, *readFile2, *regionFile, *output; /*In one sample cae, readFile1 should be provided; in two sample cae, both readFile1 and readFile2 should be provided*/
	char *infile1, *infile2, *region;
	char *chrom;
	int numSample; /*numSample should be 1 or 2*/
        } args_countRead;


static void explain_command(char **argv);
static args_countRead option_assign(int argc, char **argv);


int main(int argc, char **argv)
{	args_countRead args;
	double *region=NULL, *read1=NULL, *read2=NULL;
	int nrowRegion, ncolRegion, nrowRead1=0, ncolRead1=0, nrowRead2=0,ncolRead2=0;
	int i,j;
	ITree itree=NULL;
	ITreeIntervalSet Iset;
	double ** readcount = NULL; /*This will store the number of reads (observed and expected; case and control) in the given region*/
	double *tmp;

	args = option_assign(argc,argv);

	read1 = read_table(args.readFile1, &nrowRead1, &ncolRead1, -1,0);
	if(ncolRead1<4){
		fprintf(stderr,"%s in incorrect format: expecting at least 4 columns, get %d columns\n",args.infile1,ncolRead1);
		exit(1);
		}

	if(args.readFile2!=NULL){
		read2 = read_table(args.readFile2, &nrowRead2, &ncolRead2, -1,0);
		if(ncolRead2<4){
			fprintf(stderr,"%s in incorrect format: expecting at least 4 columns, get %d columns\n",args.infile2,ncolRead2);
			exit(1);
			}
		if(nrowRead1 != nrowRead2){
			fprintf(stderr,"Case and control data have different number of rows: %d rows in %s, %d rows in %s\n", nrowRead1, args.infile1, nrowRead2, args.infile2);
			exit(1);
			}
		}
	
        region = read_table(args.regionFile, &nrowRegion, &ncolRegion, -1,0);
        if(ncolRegion!=2){
                fprintf(stderr,"%s in incorrect Format: expecting 2 columns, got %d columns\n", args.region,ncolRegion);
                exit(1);
                }

	readcount = (double **) malloc(sizeof(double *)*(nrowRegion+10)); /*Each row correspond to one region*/
	if(readcount==NULL){
		fprintf(stderr,"Error, memory allocation failed\n"); exit(1);
		}
	for(i=0;i<nrowRegion;i++){
		readcount[i] = (double *) malloc(sizeof(double)*10);
		if(readcount[i]==NULL) {fprintf(stderr,"Error: memory allocation failed\n"); exit(1);}
		readcount[i][0] = readcount[i][1] = readcount[i][2] = readcount[i][3] = readcount[i][4] = readcount[i][5] = 0.0; 
		/*readcount[i][0] is the number of bins overlapping with the ith region
		  readcount[i][1] is the total size of the bins overlapping with the ith region
		  readcount[i][2] the observed case count
		  readcount[i][3] the expected case count
		  readcount[i][4] the observed control count
		  readcount[i][5] the expected control count
   		*/
		}

	Iset.size = 0;
	Iset.ll = 0;
	Iset.Is = NULL;
	if(nrowRegion>0){
		/*First create an interval tree*/
		itree = ITree_create(nrowRegion);
		for(i=0;i<nrowRegion;i++){
			ITree_addInterval(itree,region[ncolRegion*i],region[ncolRegion*i+1],(void *) readcount[i]);
			}
		ITreeConstruct(itree);	

		/*now handle the queries*/
		for(i=0;i<nrowRead1;i++){
			ITreeQuery(itree,read1[ncolRead1*i],read1[ncolRead1*i+1],&Iset);
			for(j=0;j<Iset.ll;j++){
				tmp = (double *) Iset.Is[j]->pt;
				tmp[0] += 1.0;
				tmp[1] += read1[ncolRead1*i+1] - read1[ncolRead1*i] + 1;
				tmp[2] += read1[ncolRead1*i+2];
				tmp[3] += read1[ncolRead1*i+3];
				}
			Iset.ll = 0;
			}

		Iset.ll = 0;
		for(i=0;i<nrowRead2;i++){ /*if nrowRead2 = 0 , i.e. read2=NULL, this part will not be performed*/
			ITreeQuery(itree,read2[ncolRead2*i],read2[ncolRead2*i+1],&Iset);
			for(j=0;j<Iset.ll;j++){
				tmp = (double *) Iset.Is[j]->pt;
				tmp[4] += read2[ncolRead2*i+2];
				tmp[5] += read2[ncolRead2*i+3];
				}
			Iset.ll = 0;
			}
		}


	/*Now print the result*/
	for(i=0;i<nrowRegion;i++){
		if(args.chrom!=NULL) fprintf(args.output,"%s\t",args.chrom);
		fprintf(args.output,"%.0f\t%.0f\t",region[ncolRegion*i],region[ncolRegion*i+1]);
		fprintf(args.output,"%.0f\t",readcount[i][0]);
		if(readcount[i][1] > 5.0*(region[ncolRegion*i+1]-region[ncolRegion*i])){ 
			/*overall size of the bins overlapping with this region is too LARGE, 
			  This means that the bin size in the normalized data is too large or the region has low mappability
			  In any case, the genotyping in this region would be not reliable, so I will set the read counts in these regions as 0
			*/
			readcount[i][2] = readcount[i][3] = readcount[i][4] = readcount[i][5] = 0.0;
			}

		fprintf(args.output,"%g\t%g",readcount[i][2],readcount[i][3]);
		if(args.numSample==2){
			fprintf(args.output,"\t%g\t%g",readcount[i][4],readcount[i][5]);
			}
		fprintf(args.output,"\n");
		}

	return 0;
}


static void explain_command(char **argv)
{       fprintf(stderr,"Usage: %s [options] <regionFile>\n", argv[0]);
        fprintf(stderr,"Options:\n");
        fprintf(stderr,"       -t <string>: the normalized data for the case genome; This is is an required option\n");
        fprintf(stderr,"       -c <string>: the normalized data for the control genome; This is optional.\n");
	fprintf(stderr,"       -o <string>: the output file; Default is <stdout>\n");
	fprintf(stderr,"       -h: print this message\n");
        fprintf(stderr,"       --chrom <string>: the chromosome name\n");
        return;
}


static args_countRead option_assign(int argc, char **argv)
{       args_countRead args;
        int c,option_index;
	char *out = NULL;
        static struct option long_options[] =
                { {"chrom", required_argument,0,0},
                  {0,0,0,0}
                };

	args.infile1 = NULL;
	args.infile2 = NULL;
	args.region = NULL;
        args.readFile1 = NULL;
        args.readFile2 = NULL;
	args.regionFile = NULL;
        args.output = stdout;
	args.chrom = NULL;
	args.numSample = 1;

        while((c=getopt_long(argc,argv,"t:c:o:h",long_options,&option_index))!=-1){
                switch(c){
                        case 0:
                                args.chrom = strdup(optarg);
                                break;
                        case 't':
                                args.infile1 = strdup(optarg);
                                break;
                        case 'c':
                                args.infile2 = strdup(optarg);
                                break;
			case 'o':
				out = strdup(optarg);
				break;
                        case 'h':
                                explain_command(argv);
                                exit(0);
                                break;
                        case '?':
                                exit(1);
                                break;
                        default:
                                abort();
                        }

                }
        if (argc - optind!=1){
                explain_command(argv);
                exit(1);
                }

	args.region = strdup(argv[optind]);
	args.regionFile = fopen(args.region,"r");
	if(args.regionFile==NULL){
		fprintf(stderr,"fopen %s: %s\n", args.region, strerror(errno));
		}

	if(out!=NULL){
		args.output = fopen(out,"w");
		if(args.output==NULL){
			fprintf(stderr,"fopen %s: %s\n",out,strerror(errno));
			}
		}


	if(args.infile1==NULL){
		fprintf(stderr,"Error, the option -t is required\n");
		exit(1);
		}else{
		args.readFile1 = fopen(args.infile1,"r");
		if(args.readFile1 == NULL){
			fprintf(stderr,"fopen %s: %s\n", args.infile1, strerror(errno));
			}
		}

	if(args.infile2!=NULL){
		args.numSample = 2;
		args.readFile2 = fopen(args.infile2, "r");
		if(args.readFile2 == NULL){
			 fprintf(stderr,"fopen %s: %s\n", args.infile2, strerror(errno));
			}
		}

        return args;
}
    
