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
static double calcuatePvalue(double *read, int nrowRead, int ncolRead, double obs_freq, int binNum, int nsample, int control, int normalApproximate);
/*control is a boolean variable, 1 means calculating the pvalue with a control, 0 means calculating the pvalue without a control (even if there is a control)
  Note that if there is no control in read (i.e.e ncolRead==2), the argument control has no effect
  if binNum = 0, reutrn 0;
  
  normalApproximate is a boolean variable, 1 means use normal approximation for calucating the pvalue with normal approximationg, 0 means no.
*/
static double calucateObsFreq(double tumor, double tumor_exp, double normal, double normal_exp, int control);
/*if control = 1, use the control for caluclating the observed Frequency
  Otherwise, only use the tumor and tumor_expected. 
 */

static void getTotalRead(double *read, int nrowRead, int ncolRead, double *totalRead);

int main(int argc, char **argv)
{	FILE *segfile=NULL, *readfile = NULL;
	int nrowSeg,ncolSeg, nrowRead, ncolRead;
	double *seg=NULL, *read=NULL;
	double *readHghLevel = NULL;
	int nrowReadHghLevel;
	int segtotalbin;
	int i, j;
	double pvalue;
	int nSample = 200;
	double obs_freq, tumorExp_overall_freq=0.5, tumorNormal_overall_Freq=0.5;
	double total_read[4];
	int normalApproximate;

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
	if(ncolSeg!=5 && ncolSeg!=7){
		fprintf(stderr,"%s in incorrect Format: expecting 5 or 7 columns, got %d columns\n", argv[1],ncolSeg);
		exit(1);
		}
	fclose(segfile);
	
	read = read_table(readfile, &nrowRead, &ncolRead, -1,0);
	if(ncolRead!=2 && ncolRead!=4){
		fprintf(stderr,"%s in incorrect format: expecting 2 or 4 columns, get %d columns\n",argv[2],ncolRead);
		exit(1);
		}
	if((ncolRead==2&&ncolSeg!=5)||(ncolRead==4&&ncolSeg!=7)){
		fprintf(stderr,"Expecting 2 and 5 (or 4 and 7) columns in %s and %s, got %d and %d columns\n", argv[2], argv[1], ncolRead, ncolSeg);
		exit(1);
		}

	getTotalRead(read, nrowRead, ncolRead, total_read);
	tumorExp_overall_freq = calucateObsFreq(total_read[0], total_read[1], 0.0, 0.0, 0);
	if(ncolRead==4){
		tumorNormal_overall_Freq = calucateObsFreq(total_read[0], total_read[1], total_read[2], total_read[3], 1);
		}

	segtotalbin = 0;
	for(i = 0; i < nrowSeg; i++){
		segtotalbin += (int) seg[i*ncolSeg+2];
		}
	if(segtotalbin > nrowRead){
		fprintf(stderr,"Error: more bins in %s (%d) than in %s (%d)\n",argv[1],segtotalbin,argv[2],nrowRead);
		}

	readHghLevel = BinMerge2HhgLevel(read, nrowRead, ncolRead, &nrowReadHghLevel);

	//segfile = fopen(argv[1],"w");
        if(segfile==NULL){
                fprintf(stderr,"No such file or directory: %s\n",argv[1]);
                exit(1);
                }
	if(ncolSeg==5){
		fprintf(stdout,"start\tend\tbinNum\tobserved\texpected\tlog2.copyRatio\tpvalue\n");
		}
	if(ncolSeg==7){
		fprintf(stdout,"start\tend\tbinNum\ttumor\ttumor_expect\tnormal\tnormal_expect\tlog2.copyRatio\tlog2.TumorExpectRatio\tpvalue\tpvalue.TumorVsExpected\n");
		}
	for(i=0;i<nrowSeg;i++){
		for(j=0;j<ncolSeg;j++){
			if(j>3){
				fprintf(stdout,"%g\t",seg[i*ncolSeg+j]);
				}else{
				fprintf(stdout,"%.0f\t",seg[i*ncolSeg+j]);
				}
			}
		if(ncolSeg==5){
			if(seg[i*ncolSeg+3]+seg[i*ncolSeg+4]==0.0){
				fprintf(stdout,"0\t");
				}else{
				fprintf(stdout,"%g\t",log2((seg[i*ncolSeg+3]+1e-10)/(seg[i*ncolSeg+4]+1e-10) + 1e-10) - log2(total_read[0]/(total_read[1]+1e-10) + 1e-10)); /*log2.copyRatio*/
				}
			}
		else if(ncolSeg==7){
			double tmp1, tmp2;
			if(seg[i*ncolSeg+3]+seg[i*ncolSeg+4]+seg[i*ncolSeg+5]+seg[i*ncolSeg+6]==0.0){
				fprintf(stdout,"0\t0\t");
				}else{
				tmp1 = ((seg[i*ncolSeg+3]+1e-10)/(seg[i*ncolSeg+4]+1e-10))/(total_read[0]/(total_read[1]+1e-10));
				tmp2 = ((seg[i*ncolSeg+5]+1e-10)/(seg[i*ncolSeg+6]+1e-10))/(total_read[2]/(total_read[3]+1e-10));
				fprintf(stdout,"%g\t",log2((tmp1+1e-10)/(tmp2+1e-10))+1e-10);/*log2.copyRatio*/
				fprintf(stdout,"%g\t",log2((seg[i*ncolSeg+3]+1e-10)/(seg[i*ncolSeg+4]+1e-10) + 1e-10) - log2(total_read[2]/(total_read[3]+1e-10) + 1e-10)); /*log2.TumorExpectRatio*/
				}
			}

		if(ncolRead==2) {
			obs_freq = calucateObsFreq(seg[i*ncolSeg+3], seg[i*ncolSeg+4], 0.0, 0.0, 0);
			obs_freq = seg[i*ncolSeg+3]+seg[i*ncolSeg+4] == 0 ? tumorExp_overall_freq: obs_freq;
			normalApproximate = (seg[i*ncolSeg+4] > 100) ? 1 : 0;
			}
		else {
			obs_freq = calucateObsFreq(seg[i*ncolSeg+3], seg[i*ncolSeg+4], seg[i*ncolSeg+5], seg[i*ncolSeg+6], 1);
			obs_freq = seg[i*ncolSeg+3]+seg[i*ncolSeg+4]+seg[i*ncolSeg+5]+seg[i*ncolSeg+6] ==0 ? tumorNormal_overall_Freq : obs_freq;
			normalApproximate = (seg[i*ncolSeg+4] > 100 && seg[i*ncolSeg+6] > 100) ? 1: 0;
			}

		if(normalApproximate==1){
			nSample = 1000;
			}else{
			nSample = 100000;
			}
		if(seg[i*ncolSeg+2]>50*HghLevelBinCount){
			pvalue = calcuatePvalue(readHghLevel, nrowReadHghLevel, ncolRead,obs_freq,(int) seg[i*ncolSeg+2]/HghLevelBinCount, nSample, 1, normalApproximate);
			}else{
			pvalue = calcuatePvalue(read, nrowRead, ncolRead,obs_freq,(int) seg[i*ncolSeg+2], nSample, 1, normalApproximate);
			}
		fprintf(stdout,"%g",pvalue);

		if(ncolRead==4){ /*calculate p-value for tumor without a control*/
			obs_freq = calucateObsFreq(seg[i*ncolSeg+3], seg[i*ncolSeg+4], 0.0, 0.0, 0);;
			obs_freq = seg[i*ncolSeg+3]+seg[i*ncolSeg+4] == 0 ? tumorExp_overall_freq : obs_freq;
			if(seg[i*ncolSeg+2]>50*HghLevelBinCount){
				pvalue = calcuatePvalue(readHghLevel, nrowReadHghLevel, ncolRead,obs_freq,(int) seg[i*ncolSeg+2]/HghLevelBinCount, nSample, 0, normalApproximate);
				}else{
				pvalue = calcuatePvalue(read, nrowRead, ncolRead,obs_freq,(int) seg[i*ncolSeg+2], nSample, 0, normalApproximate);
				}
			fprintf(stdout,"\t%g",pvalue);
			}
		fprintf(stdout,"\n");
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


static void getTotalRead(double *read, int nrowRead, int ncolRead, double *totalRead)
{	int i, j;
	
	for(j=0;j<ncolRead;j++){
		totalRead[0] = 0.0;
		}

	for(i=0;i<nrowRead;i++){
		for(j=0;j<ncolRead;j++){
			totalRead[j] += read[i*ncolRead + j];
			}
		}

	return;
}


static double calucateObsFreq(double tumor_obs, double tumor_exp, double normal_obs, double normal_exp, int control)
{	double obsFreq;

	if(control==0){
		obsFreq = (tumor_obs+1e-10)/(tumor_obs+tumor_exp+1e-10);
		}else{
		double tmp1, tmp2;
		tmp1 = (tumor_obs+1e-10)/(tumor_exp+1e-10);
		tmp2 = (normal_obs+1e-10)/(normal_exp + 1e-10);
		obsFreq = (tmp1+1e-10)/(tmp1+tmp2+1e-10);
		}
	return obsFreq;
}


static double calcuatePvalue(double *read, int nrowRead, int ncolRead, double obs_freq,int binNum, int nsample_in, int control, int normalApproximate)
{	int index,i,k, index_start;
	static double *obs_freq_sample=NULL;
	static int size=0;
	double obs, expected, normal_obs, normal_expected;
	double meanFreq, varFreq, sdFreq, pvalue, zscore;
	double cnt_lowtail, cnt_hightail;
	int nsample;

	nsample = nsample_in < 100 ? 100 : nsample_in;
	if(size < nsample){
		size = nsample+100;
		obs_freq_sample = (double *) realloc(obs_freq_sample,sizeof(double)*size);
		if(obs_freq_sample==NULL){
			fprintf(stderr,"Error in calcuateZscore: memory allocation failed\n");
			exit(1);
			}
		}
	if(binNum<0){
		fprintf(stderr,"Error: binNum=%d is nonpositive\n",binNum);
		exit(1);
		}
	if(binNum==0){
		return 1.0;
		}

	meanFreq = 0.0;
	cnt_lowtail = 0.0;
	cnt_hightail = 0.0;
	for(i=0;i<nsample;i++){
		obs = 0.0;
		expected = 0.0;
		normal_obs = 0.0;
		normal_expected = 0.0;
		index_start = (int) floor(rand_lp()*nrowRead);
		for(k=0;k<binNum;k++){
			//index = (int) floor(rand_lp()*nrowRead);
			index = (index_start+k)%nrowRead;
			obs += read[ncolRead*index];
			expected += read[ncolRead*index+1];
			if(ncolRead==4 && control==1){
				normal_obs += read[ncolRead*index+2];
				normal_expected += read[ncolRead*index+3];
				}
			}
		if(ncolRead==2 || (ncolRead==4 && control == 0)){
			obs_freq_sample[i] = calucateObsFreq(obs,expected,0.0,0.0, 0);
			meanFreq += obs_freq_sample[i];
			}else if(ncolRead==4 && control == 1){
			obs_freq_sample[i] = calucateObsFreq(obs,expected,normal_obs,normal_expected, 1);
			meanFreq += obs_freq_sample[i];
			}
		cnt_lowtail += (obs_freq_sample[i] <= obs_freq+0.01);
		cnt_hightail += (obs_freq_sample[i] >= obs_freq-0.01);
		}
	meanFreq = meanFreq/nsample;

	varFreq = 0.0;
	for(i=0;i<nsample;i++){
		varFreq += (obs_freq_sample[i]-meanFreq)*(obs_freq_sample[i]-meanFreq);
		}
	varFreq = varFreq/(nsample-1);
	sdFreq = sqrt(varFreq);
	if(sdFreq<0.01) sdFreq = 0.01; /*Make sure sdFreq not too small*/

	zscore =  (obs_freq-meanFreq)/sdFreq;
	if(normalApproximate){
		pvalue = pchisq(1,zscore*zscore,0);
		}else{
		pvalue = cnt_lowtail < cnt_hightail ? 2*(cnt_lowtail)/nsample : 2*(cnt_hightail)/nsample;
		pvalue = pvalue < 1.0 ? pvalue : 1.0;
		pvalue = pvalue > 0.0 ? pvalue : 1.0/nsample;
		}

	return pvalue;
}

