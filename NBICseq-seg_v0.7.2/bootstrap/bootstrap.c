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
#include "../lib/gamma.h"
#include <getopt.h>
#include <math.h>
#include <unistd.h>

static int HghLevelBinCount = 100; /*For computational efficiency, if a segment has too many bins (i.e. the segment is too large), I will bootstrap on bins of HghLevelBinCount*BinSize*/

static double * BinMerge2HhgLevel(double *read, int nrowRead, int ncolRead, int *nrowReadHghLevel);
static double calcuateZscore(double *read, int nrowRead, int ncolRead, double obs_freq, int binNum, int nsample, int control);
/*control is a boolean variable, 1 means calculating the z-score with a control, 0 means calculating the z-score without a control (even if there is a control)
  Note that if there is no control in read (i.e.e ncolRead==2), the argument control has no effect
*/

int main(int argc, char **argv)
{	FILE *segfile=NULL, *readfile = NULL;
	int nrowSeg,ncolSeg, nrowRead, ncolRead;
	double *seg=NULL, *read=NULL;
	double *readHghLevel = NULL;
	int nrowReadHghLevel;
	int segtotalbin;
	int i, j;
	double zscore, pvalue;
	int nSample = 100;
	double obs_freq;

	if(argc!=3){
		fprintf(stderr,"%s <segFile> <ReadFile>\n",argv[0]);
		exit(1);
		}


	segfile = fopen(argv[1],"r");
	readfile = fopen(argv[2],"r");
	if(segfile==NULL){
		fprintf(stderr,"No such file or directory: %s\n",argv[1]);
		exit(1);
		}
	if(readfile==NULL){
		fprintf(stderr,"No such file or directory: %s\n",argv[2]);
		exit(1);
		}


	seg = read_table(segfile, &nrowSeg, &ncolSeg, -1,0);
	if(ncolSeg!=6 && ncolSeg!=9){
		fprintf(stderr,"%s in incorrect Format: expecting 6 or 9 columns, got %d columns\n", argv[1],ncolSeg);
		exit(1);
		}
	fclose(segfile);
	
	read = read_table(readfile, &nrowRead, &ncolRead, -1,0);
	if(ncolRead!=2 && ncolRead!=4){
		fprintf(stderr,"%s in incorrect format: expecting 2 or 4 columns, get %d columns\n",argv[2],ncolRead);
		exit(1);
		}

	segtotalbin = 0;
	for(i = 0; i < nrowSeg; i++){
		segtotalbin += (int) seg[i*ncolSeg+2];
		}
	if(segtotalbin > nrowRead){
		fprintf(stderr,"Error: more bins in %s (%d) than in %s (%d)\n",argv[1],segtotalbin,argv[2],nrowRead);
		}

	readHghLevel = BinMerge2HhgLevel(read, nrowRead, ncolRead, &nrowReadHghLevel);

	segfile = fopen(argv[1],"w");
        if(segfile==NULL){
                fprintf(stderr,"No such file or directory: %s\n",argv[1]);
                exit(1);
                }
	if(ncolSeg==6){
		fprintf(segfile,"start\tend\tbinNum\tobserved\texpected\tlog2.copyRatio\tpvalue\n");
		}
	if(ncolSeg==9){
		fprintf(segfile,"start\tend\tbinNum\ttumor\ttumor_expect\tnormal\tnormal_expect\tlog2.copyRatio\tlog2.TumorExpectRatio\tpvalue\tpvalue.TumorVsExpected\n");
		}
	for(i=0;i<nrowSeg;i++){
		for(j=0;j<ncolSeg;j++){
			if(j>3){
				fprintf(segfile,"%g\t",seg[i*ncolSeg+j]);
				}else{
				fprintf(segfile,"%.0f\t",seg[i*ncolSeg+j]);
				}
			}
		if(ncolRead==2) {obs_freq = seg[i*ncolSeg+3]/(seg[i*ncolSeg+3]+seg[i*ncolSeg+4]+1e-10);}
		else {
			double tmp1, tmp2;
			tmp1 = seg[i*ncolSeg+3]/(seg[i*ncolSeg+4]+1e-10);
			tmp2 = seg[i*ncolSeg+5]/(seg[i*ncolSeg+6]+1e-10);
			obs_freq = tmp1/(tmp1+tmp2+1e-10);
			}
		if(seg[i*ncolSeg+2]>50*HghLevelBinCount){
			zscore = calcuateZscore(readHghLevel, nrowReadHghLevel, ncolRead,obs_freq,(int) seg[i*ncolSeg+2]/HghLevelBinCount, nSample, 1);
			}else{
			zscore = calcuateZscore(read, nrowRead, ncolRead,obs_freq,(int) seg[i*ncolSeg+2], nSample, 1);
			}
			
		//fprintf(segfile,"%g\t",zscore);
		pvalue = pchisq(1,zscore*zscore,0);
		fprintf(segfile,"%g",pvalue);

		if(ncolRead==4){ /*calculate p-value for tumor without a control*/
			obs_freq = seg[i*ncolSeg+3]/(seg[i*ncolSeg+3]+seg[i*ncolSeg+4]+1e-10);
			if(seg[i*ncolSeg+2]>50*HghLevelBinCount){
				zscore = calcuateZscore(readHghLevel, nrowReadHghLevel, ncolRead,obs_freq,(int) seg[i*ncolSeg+2]/HghLevelBinCount, nSample, 0);
				}else{
				zscore = calcuateZscore(read, nrowRead, ncolRead,obs_freq,(int) seg[i*ncolSeg+2], nSample, 0);
				}
			pvalue = pchisq(1,zscore*zscore,0);
			fprintf(segfile,"\t%g",pvalue);
			}
		fprintf(segfile,"\n");
		}

	return 0;
}



static double * BinMerge2HhgLevel(double *read, int nrowRead, int ncolRead, int *nrowReadHghLevel)
{	int k,i,j, nrowNewBin;
	double *newBin;

	nrowNewBin = (nrowRead+1)/HghLevelBinCount;

	newBin = (double *) malloc(sizeof(double)*(ncolRead*nrowNewBin+10));
	if(newBin==NULL){
		fprintf(stderr,"Error in BinMerge2HhgLevel: memory allocation failed\n");
		exit(1);
		}

	for(k = 0;k < nrowNewBin; k++){
		for(j=0;j<ncolRead;j++){
			newBin[ncolRead*k+j] = 0.0;
			}

		for(i=HghLevelBinCount*k; i < HghLevelBinCount*(k+1) && i<nrowRead; i++){
			for(j=0;j<ncolRead;j++){
				newBin[ncolRead*k+j] += read[ncolRead*i+j];
				}
			}
		}

	*nrowReadHghLevel = nrowNewBin;


	return newBin;

}


static double calcuateZscore(double *read, int nrowRead, int ncolRead, double obs_freq,int binNum, int nsample_in, int control)
{	int index,i,k, index_start;
	static double *obs_freq_sample=NULL;
	static int size=0;
	double obs, expected, normal_obs, normal_expected;
	double meanFreq, varFreq, sdFreq;
	int nsample;

	nsample = nsample_in < 20 ? 20 : nsample_in;
	if(size < nsample){
		size = nsample+100;
		obs_freq_sample = (double *) realloc(obs_freq_sample,sizeof(double)*size);
		if(obs_freq_sample==NULL){
			fprintf(stderr,"Errorin calcuateZscore: memory allocation failed\n");
			exit(1);
			}
		}
	if(binNum<=0){
		fprintf(stderr,"Error: binNum=%d is nonpositive\n",binNum);
		exit(1);
		}

	meanFreq = 0.0;
	for(i=0;i<nsample;i++){
		obs = 0.0;
		expected = 0.0;
		normal_obs = 0.0;
		normal_expected = 0.0;
		index_start = (int) floor(rand_lp()*nrowRead);
		for(k=0;k<binNum;k++){
			index = (int) floor(rand_lp()*nrowRead);
			//index = (index_start+k)%nrowRead;
			obs += read[ncolRead*index];
			expected += read[ncolRead*index+1];
			if(ncolRead==4 && control==1){
				normal_obs += read[ncolRead*index+2];
				normal_expected += read[ncolRead*index+3];
				}
			}
		if(ncolRead==2 || (ncolRead==4 && control == 0)){
			obs_freq_sample[i] = obs/(obs+expected+1e-10); /*make sure the denominator is nonzero*/
			meanFreq += obs_freq_sample[i];
			}else if(ncolRead==4 && control == 1){
			double tmp1, tmp2;
			tmp1 = obs/(expected+1e-10); /*normalized tumor*/
			tmp2 = normal_obs/(normal_expected+1e-10); /*normalized normal*/
			obs_freq_sample[i] = tmp1/(tmp1+tmp2+1e-10);
			meanFreq += obs_freq_sample[i];
			}
		}
	meanFreq = meanFreq/nsample;

	varFreq = 0.0;
	for(i=0;i<nsample;i++){
		varFreq += (obs_freq_sample[i]-meanFreq)*(obs_freq_sample[i]-meanFreq);
		}
	varFreq = varFreq/(nsample-1);
	sdFreq = sqrt(varFreq);
	//fprintf(stdout,"mean = %g, sd = %g\n", meanFreq, sdFreq);
	if(sdFreq<0.01) sdFreq = 0.01; /*Make sure sdFreq not too small*/

	return (obs_freq-meanFreq)/sdFreq;
}

