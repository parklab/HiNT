/* Written by ruibin Xi
 * The input file should be of the form
 *
 * bin start | bin end | tumor1 count |matched normal1 count| ... | tumorK count
| matched normalK count|

Modified on Feb 18 2009 so that the tumor or normal count could be non-integer
 * */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include "bin.h"
#include <getopt.h>


typedef struct args_mbic_t{
	FILE *input;
	int remove_commonCNV; /*whether or not remove the candidate CNVs common to all samples*/
	double lambda;
	} args_mbic;


bin_list *remove_false_CNV(bin_list *head, double cutoff); /*remove the bins with log2 copy ratios of all samples less than -cutoff or greater than cutoff, as these are highly likely to be false positives*/


static void explain_command(char **argv);
static args_mbic option_assign(int argc, char **argv);

int main(int argc, char **argv)
{
	bin_list *head, *bin_tmp;
	int n_bins,n_merged,level,flag=0;
	args_mbic args;
	
	args = option_assign(argc,argv);
	set_lambda(args.lambda);
	
	head = bin_read(args.input,&n_bins);
	fprintf(stderr,"Loaded %d bins\n",n_bins);

	/*merge the neighboring bins*/
	flag = 0;level = 2;
        while(flag==0){
                do {
                        n_merged = MBIC_seq(head,level);
                        fprintf(stderr,"Level %d: merged %d bins, remaining %d\n",level,n_merged, ll_bins(head));
                        ++level;
                } while (n_merged>0);

                level = level -1;
                if(level<=2&&n_merged==0) flag = 1; else level = level-1;
        }

	/*remove the likely false positive and further merge the bins*/
	if(args.remove_commonCNV==1){
		head = remove_false_CNV(head,0.5);
		if(head!=NULL){
		        flag = 0;level = 2;
			while(flag==0){
		                do {
		                        n_merged = MBIC_seq(head,level);
		                        fprintf(stderr,"Level %d: merged %d bins, remaining %d\n",level,n_merged, ll_bins(head));
	        	                ++level;
		                } while (n_merged>0);

	        	        level = level -1;
		                if(level<=2&&n_merged==0) flag = 1; else level = level-1;
	        		}
			}
		}	

	/*print the merged bins*/
	bin_list_print(head);
	/*free the space allocated to the list head*/
	bin_tmp=head;
	while(bin_tmp!=NULL)
		{ head = head->next;
		  free(bin_tmp);
		  bin_tmp = head==NULL ? NULL : head->next;
		}

	return 0;
}


bin_list *remove_false_CNV(bin_list *head_in, double cutoff)
{	bin_list *node, *node1, *head;
	int i, flag1,k=0;
	double *tumor_total, *normal_total, *log2ratio, tmp;

	if(head_in->n_sample<=1) return head_in;
	head = head_in;

	/*first calculate the overall log2 copy ratio for each sample*/
	tumor_total = (double *) malloc(sizeof(double)*head->n_sample); 
	normal_total = (double *) malloc(sizeof(double)*head->n_sample);
	log2ratio = (double *) malloc(sizeof(double)*head->n_sample);

	if(tumor_total==NULL||normal_total==NULL||log2ratio==NULL){
		fprintf(stderr,"Error in remove_false_CNV: memory allocation failed\n");
		exit(1);
		}

	for(i=0;i<head->n_sample;i++){
		tumor_total[i] = 0.0;
		normal_total[i] = 0.0;
		}

	node=head;
	while(node!=NULL){
		for(i=0;i<node->n_sample;i++){
			tumor_total[i] += node->tumor[i];
			normal_total[i] += node->normal[i];
			}
		node = node->next;
		}
	for(i=0;i<head->n_sample;i++){
		log2ratio[i] = log(tumor_total[i]/normal_total[i]+0.001)/log(2);
		}


	node = head;
	k = 0;
	while(node!=NULL){
		flag1 = 1; i = 0;
		while(flag1==1&&i<node->n_sample){
			tmp = log(node->tumor[i]/node->normal[i]+0.001)/log(2);
			if(tmp < fabs(cutoff) && tmp < log2ratio[i] + fabs(cutoff)){ flag1 = 0;}
			i++;
			}


		if(flag1==1){ /*all log2 copy ratio >= fabs(cutoff)*/
			node1 = node;
			node = node->next;
			if(node1==head) head = node1->next; /*head will be removed, update head*/
			node_del(node1);
			k++;
			}
		else{ node = node->next;}
		}
	fprintf(stderr,"Removed %d CNV regions that are likely to be germline CNVs\n",k);

	return head;
}


static void explain_command(char **argv)
{	fprintf(stderr,"Usage: %s [options]\n", argv[0]);
	fprintf(stderr,"Options:\n");
	fprintf(stderr,"       -i <string>: the input file name; default stdin\n");
	fprintf(stderr,"       -l <float>: the penalty lambda of MBIC-seq; default 1.2\n");
	fprintf(stderr,"       --rm: remove the candidate CNV regions that are common to all samples;\n");
	fprintf(stderr,"           use this only if the reference is the expected and both tumor and normal present in the binned data,\n");
	fprintf(stderr,"           in which case these regions are likely to be false positives or germline CNVs\n");

	return;	
}


static args_mbic option_assign(int argc, char **argv)
{	args_mbic args;
	char *infile=NULL;
	int c,option_index;
	static struct option long_options[] =
		{ {"rm", no_argument ,0,0},
		  {0,0,0,0}
		};

	args.lambda = 1.2;
	args.remove_commonCNV = 0;
	args.input = stdin;

	while((c=getopt_long(argc,argv,"i:hl:",long_options,&option_index))!=-1){
		switch(c){
			case 0: 
				args.remove_commonCNV = 1;
				break;
			case 'i':
				infile = strdup(optarg);
				break;
			case 'l':
				args.lambda = atof(optarg);
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
	if (argc - optind!=0){
		explain_command(argv);
		exit(1);
		}

	if(args.lambda<=0.0){
		fprintf(stderr,"lambda must be positive\n");
		exit(1);
		}
	if(infile!=NULL){
		args.input = fopen(infile,"r");
		if(args.input==NULL){
			fprintf(stderr,"fopen %s: %s\n",infile,strerror(errno));
			exit(1);
			}
		}

	return args;
}
