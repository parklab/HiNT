#include <math.h>
#include "../lib/rbtree.h"
#include "../lib/read.h"

struct BIN {
	int start, end; /*start and end positions of the bin*/
	int num_bins; /*how many bins merged from the initial bins*/
	int n_sample; /*number of samples*/
	double *tumor, *normal; /*numbers of tumor and normal reads in the bin*/
	double *prob; /*percentages of the tumor reads in the bin*/
	double bic_diff; /*bic_diff of the bin, bic_diff = bic_new - bic_old*/
	struct BIN *next, *pre;
};

typedef struct BIN bin_list;
/* log_lk is the log likelihood of this bin
 * bic_diff is the bic difference if merge this bin with the next level (level>=2) bins
 * level is defined in the file bin.c
 * */

void set_lambda(double x);/*set lambda as x*/

void bin_free(bin_list * bin); /*free the memory allocated to bin*/

bin_list * bin_read(FILE *bin_file, int *n_bins); 
/* read the bin_file
 * n_bins: number of bins loaded.
 * bin_file: the file that store the bin information
 * This function also calculate the sample size while loading the data.
 *
 * Note: this function does NOT compute the bic_diff of the bins.
 * */

void node_del(bin_list *node); /*delete the given node in the bin list*/
void bin_merge(bin_list *bin, int k); /*merge the k nodes after the current bin (inclusive)*/

void bin_list_print(bin_list *head);/*print bin start from head*/

int ll_bins(bin_list *head); /*give the number of bins left*/

double BIC_diff(bin_list *bin, int k); 
/* Calculate the BIC difference 
 * bin: the current bin
 * k: the number of the bins after bin (inclusive) that would be merged.
 *
 * the returned value is (BIC after merging) - (BIC before merging) 
 *
 * if there are less than k bins after bin (inclusive), the BIC difference will be 0;
 * */

int bin_compare(void *left, void *right);
/* compare two bin based on their bic_diff 
 * if left < right return -1
 * if left > right return 1
 * if left = right return 0
 *
 * if their bic_diff are the same, use their start position to determine which is larger.
 * */




rbtree_node min_tree_node(rbtree tree);/*Get the minimum node of the tree*/

void rbtree_free(rbtree tree, bin_list *head);/*free the memory used by the rbtree*/

rbtree rbtree_ini(bin_list *head, int k); 
/* create a rbtree corresponding to the bin list head 
 * head: the bin list
 * k: number of bins attempting to merge.
 * */


int MBIC_seq(bin_list *head, int level); 
/* The multi-sample version of BIC-seq, the returned value is the number of level-tuples of bins that were merged.
 * the major function of this programe.
 * */
