#include "bin.h"


static double lambda = 1;
static double N = 0.0,logN = 0.0; /*total number of reads and log total number of reads*/
static int num_sample=1; /*number of samples*/
static int flag_BIC_diff = 1; /*flag if it's the first time to call the function BIC_diff*/


static bin_list * bin_from_array(double *line, int ncol); 
/* creat a bin from a double array line with ncol elements
 * the array line is assumed of the following form
 * start| end | tumor 1 | matched normal 1 | .... | tumor n | matched normal n
 * */


void bin_print(bin_list *bin);/*print one bin*/

static bin_list * bin_from_array(double *line, int ncol)
	{ bin_list *bin;
	  int i, n_sample;

	  bin = (bin_list *) malloc(sizeof(bin_list));
	  bin->n_sample = (ncol-2)/2;
	  n_sample = bin->n_sample;

	  bin->normal = (double *) malloc(sizeof(double)*(n_sample));
	  bin->tumor = (double *) malloc(sizeof(double)*(n_sample));
	  bin->prob = (double *) malloc(sizeof(double)*n_sample);

	  bin->start = (int) line[0];
	  bin->end = (int) line[1];
	  bin->num_bins = 1;
	  bin->pre = NULL;
	  bin->next = NULL;

	  for(i=0;i<n_sample;i++)
		{ bin->tumor[i] =  line[i*2+2];
		  bin->normal[i] =  line[i*2+3];
		  if(bin->tumor[i]<0.0 || bin->normal[i]< 0.0) {
			fprintf(stderr,"Error, \"read count\" is negative (tumor_%d=%g, normal_%d = %g)\n",i,bin->tumor[i],i,bin->normal[i]);
			exit(1);
			}
		  N += line[i*2+2]+line[i*2+3];
		  if(line[i*2+2]+ line[i*2+3] > 0) bin->prob[i] = line[i*2+2]/(line[i*2+2]+line[i*2+3]);
		  else bin->prob[i] = 0.0;
		}

	  return bin;
	
	}


void set_lambda(double x)
	{ lambda = x;
	  if(lambda<=0)
          	{ fprintf(stderr,"Error: lambda musth be a positive number.\n");
                  exit(1);
                }
	}

void bin_free(bin_list * bin)
	{ free(bin->tumor);
	  free(bin->normal);
	  free(bin->prob);
	  free(bin);
	}

bin_list * bin_read(FILE *bin_file, int *n_bins)
	{ double * line;
	  int ncol,nrow,n_sample;
	  bin_list *bins, *bin_tmp, *bin_curr; /*bins: the head of the list*/
	  int j;

	  if(bin_file==NULL)
		{ fprintf(stderr,"Error in bin_read: empty FILE stream.\n");
		  exit(1);
		}

	  
	  bins=NULL;
	  bin_tmp=NULL;
	  bin_curr = NULL;

	  line = read_table(bin_file,&nrow,&ncol,1,0);
	  n_sample = ncol-2;
	  if(n_sample<=0||n_sample%2!=0)
		{ fprintf(stderr,"Error in bin_read: bin file is in incorrect format.\n");
		  exit(1);
		}
	  num_sample = n_sample/2;

	  j = 0;
	  while(nrow>0&&line!=NULL)
		{ if(ncol-2!=n_sample)
               		 { fprintf(stderr,"Error in bin_read: unequal number of columns.\n");
                  	   exit(1);
               		 }
		  bin_tmp = bin_from_array(line,ncol);
		  bin_tmp -> next = NULL;
		  if(bins==NULL)
			{ bins = bin_tmp;
			  bins->pre = NULL;
			  bin_curr = bins;
			}
		  else
			{ bin_curr -> next = bin_tmp;
			  bin_tmp -> pre = bin_curr;
			  bin_curr = bin_tmp;
			}
		  
		  j++; 
  		  /*printf("j= %d\n",j);
		  for(i=0;i<ncol;i++)
			{ printf("\t %10.0f", line[i]);
			}
		  printf("\n");*/


                  free(line);/*free the memory allocated to line*/
		  line = read_table(bin_file,&nrow,&ncol,1,0);
		  //printf("nrow=%d,ncol=%d\n",nrow,ncol);
		}

	  logN = N==0.0 ? 0.0 : log(N);
	  //printf("loop done!\n");
	  *n_bins = j;
	  return bins;
	}

void node_del(bin_list *node)
	{ bin_list *bin_tmp;
	  
	  if(node==NULL) return;
	  bin_tmp = node;
	  if(node->pre!=NULL&&node->next!=NULL) /*node is not the head and the tail of the list*/
		{ node->pre->next = node->next;
		  node->next->pre = node->pre;
		}
	  else if(node->pre==NULL&&node->next!=NULL) /*node is the head but not the tail*/
		{ node->next->pre = NULL;
		}
	  else if((node->pre!=NULL&&node->next==NULL))/*tail but not head*/
		{ node->pre->next=NULL;}

	  bin_free(bin_tmp); /*free the space allocated to the current node*/
	
	  //printf("node deleted.\n\n");
	  return;
	}

void bin_merge(bin_list *bin, int k)
	{ bin_list *bin_tmp,*bin_tmp2;
	  int i,j;
	
	  if(bin==NULL)
		{ fprintf(stderr,"Error in bin_merge: the bin is empty.\n");
		  exit(1);
		}
	  bin_tmp = bin->next;
	  j=0;i=0;
	  while(bin_tmp!=NULL&&j<k-1) /*update the normal and tumor read in the current bin*/
		{ for(i=0;i<num_sample;i++)
			{ bin->normal[i] += bin_tmp->normal[i];
			  bin->tumor[i] += bin_tmp->tumor[i];
			}
		  bin->end = bin_tmp->end;/*update bin position*/
		  bin->num_bins += bin_tmp->num_bins;
		  bin_tmp2 = bin_tmp;
		  bin_tmp = bin_tmp->next;
		  node_del(bin_tmp2); /*delete this bin*/
		  j++;
		}
	  if(j<k-1) 
		{ fprintf(stderr,"Error in bin_merge: asked to merged %d bins, but only %d bins left after the current bin.\n",k,j);
		  exit(1);
		}

	  for(i=0;i<num_sample;i++) /*update the prob in the current bin*/
		{ if(bin->tumor[i]+bin->normal[i]>0)
			{ bin->prob[i] = ((double) bin->tumor[i])/((double) bin->tumor[i]+bin->normal[i]);}
		  else bin->prob[i] = 0.0;
		}

	  bin_tmp = bin;
	  j=0;
	  while(bin_tmp!=NULL&&j<k) /*update the BIC difference*/
		{ bin_tmp->bic_diff = BIC_diff(bin_tmp,k);
		  bin_tmp = bin_tmp->pre;
		  j++;
		}

	  return;
	  
	}


void bin_print(bin_list *bin)
{	int i;

	if(bin==NULL) return;
	printf("start\tend\tbinNum\tbic_dif\t");
	for(i=0;i<num_sample;i++)
		{ printf("tumor_%d\tnormal_%d\t",i,i);}
	printf("\n");

	printf("%d\t%d\t%d\t\t%g",bin->start,bin->end,bin->num_bins,bin->bic_diff);
	for(i=0;i<bin->n_sample;i++){
		printf("%g\t%g\t",bin->tumor[i],bin->normal[i]);
		}
	printf("\n");
}

void bin_list_print(bin_list *head)
        { int i;
          bin_list *bins;

	  if(head==NULL) return;
          bins = head;
	  fprintf(stderr,"lambda is %g\n",lambda);
	  /*print the header of the file*/
	  printf("start\tend\tbinNum\t");
	  for(i=0;i<num_sample;i++)
		{ printf("tumor_%d\tnormal_%d\t",i,i);}
	  printf("\n");

          while(bins!=NULL)
                { printf("%d\t%d\t%d\t",bins->start,bins->end,bins->num_bins);
                  for(i=0;i<bins->n_sample;i++)
                        { printf("%g\t%g\t",bins->tumor[i],bins->normal[i]);
                        }
		  //printf("%g",bins->bic_diff);
		  printf("\n");

                  /*for(i=0;i<bins->n_sample;i++)
                        { printf("%g\t",bins->prob[i]);}
		  printf("\n");*/
		  bins = bins->next;
                }
	  return;
        }


int ll_bins(bin_list *head)
	{ bin_list *bin;
	  int i;

	  i=0;
	  bin = head;
	  while(bin!=NULL)
		{ bin = bin->next; i++;
		}
	  return i;
	}


double BIC_diff(bin_list *bin, int k)
	{ bin_list *bin_tmp;
	  static double *n_tumor, *n_normal; 
	  /*The number of normal and tumor reads in the merged bin. */
	  static double *prob;/*prob: the frequency of the tumor reads in the merged bin*/
	  int i,j;
	  double BIC_old=0.0, BIC_new=0.0, lgL;

	  if(flag_BIC_diff==1) /*it's the first time to call the function BIC_diff*/
		{ n_tumor = (double *) malloc(sizeof(double)*(num_sample+1));
		  n_normal = (double *) malloc(sizeof(double)*(num_sample+1));
		  prob = (double *) malloc(sizeof(double)*(num_sample+1));
		  flag_BIC_diff = 0;
		  /*printf("memory allocated\n");*/
		}

	  for(i=0;i<num_sample;i++)
		{ n_tumor[i] = 0.0;
		  n_normal[i] = 0.0;
		}

	  j = 0;
	  bin_tmp = bin;
	  while(j<k&&bin_tmp!=NULL)
		{ for(i=0;i<num_sample;i++)
                	{ n_tumor[i] += bin_tmp->tumor[i];
                 	  n_normal[i] += bin_tmp->normal[i];
			  if(bin_tmp->tumor[i]==0.0||bin_tmp->normal[i]==0.0)
				{ lgL = 0.0;}
			  else lgL = bin_tmp->tumor[i]*log(bin_tmp->prob[i]) + bin_tmp->normal[i]*log(1.0-bin_tmp->prob[i]);

			  //if(isnan(lgL)) {fprintf(stderr,"Erorr, nan produced (n_tumor=%g, n_normal = %g)\n", bin_tmp->tumor[i],bin_tmp->normal[i]);}

			  BIC_old +=  -2*lgL + lambda*logN;
               	 	} 
		  j++;
		  bin_tmp = bin_tmp->next;
		}
	  //printf("BIC old is %10.7f, lambda = %g\n",BIC_old,lambda);

          if(j<k) return 0.0;
	  
	  for(i=0;i<num_sample;i++)
		{ if(n_tumor[i]!=0.0&&n_normal[i]!=0.0)
			{ prob[i] = n_tumor[i]/(n_tumor[i]+n_normal[i]);
			  lgL = n_tumor[i]*log(prob[i]) + n_normal[i]*log(1-prob[i]);
			}
		  else lgL = 0.0;
		  //if(isnan(lgL)) fprintf(stderr,"Error, nan produced (n_tumor=%g, n_normal=%g)\n",n_tumor[i],n_normal[i]);

		  BIC_new += -2*lgL + lambda*logN;
		}
	  //printf("BIC new is %10.7f, lambda = %g\n\n",BIC_new,lambda);	  

	  //printf("there are %d samples and total number of reads is %g, logN is %10.7f\n",num_sample,N,logN);
	  /*if(isnan(BIC_new-BIC_old)) 
		{ printf("Error in BIC_diff: NAN produced.\n");
		  printf("BIC_old = %g\t BIC_new = %g\n",BIC_old,BIC_new);
		  exit(1);
		}*/
	  //if(isnan(BIC_new-BIC_old)) fprintf(stderr,"Error, nan produced. BIC_new = %g, BIC_old=%g\n",BIC_new, BIC_old);
	  return BIC_new-BIC_old;

	}




int bin_compare(void *left, void *right)
        { int flag=0;
          bin_list *left_bin,*right_bin;

          if(left==NULL||right==NULL)
                { fprintf(stderr,"Error in compare_cluster: empty pointer\n");
                  exit(1);
                }

          left_bin = (bin_list *) left;
          right_bin = (bin_list *) right;

          if(left_bin->bic_diff-right_bin->bic_diff<0.0) flag = -1;
          else if(left_bin->bic_diff-right_bin->bic_diff>0.0) flag = 1;
          else
                { if(left_bin->start<right_bin->start) flag = -1;
                  else if(left_bin->start>right_bin->start) flag = 1;
                  else flag=0;
                }
          return flag;
        }




rbtree_node min_tree_node(rbtree tree)
        { rbtree_node tree_node;
          if(tree==NULL)
                { fprintf(stderr,"Error in min_node: tree does not exit.\n");
                }
          if(tree->root==NULL) return NULL;
          tree_node = tree->root;
          while(tree_node!=NULL&&tree_node->left!=NULL)
                { tree_node = tree_node->left;
                }
          return tree_node;
        }


void rbtree_free(rbtree tree,bin_list *head)
                { bin_list *nd;
                  nd = head;
                  while(nd!=NULL)
                        { if(rbtree_delete(tree, (void*) nd, bin_compare)==1){
				bin_print(nd);
				}
                          nd = nd->next;
                        }
                if(tree->root!=NULL) {fprintf(stderr,"Warining in rbtree_free: tree is still not empty.\n");}
                free(tree);

                return;
                }

rbtree rbtree_ini(bin_list *head, int k)
	{ bin_list *bin_tmp;
	  rbtree tree;

	  tree = rbtree_create(); /*intitilize the tree*/

	  bin_tmp = head;
	  while(bin_tmp!=NULL)
		{ bin_tmp->bic_diff = BIC_diff(bin_tmp,k);
		  rbtree_insert(tree,(void *)bin_tmp,(void *)bin_tmp,bin_compare); 
		  bin_tmp = bin_tmp->next;
		}

	 return tree;
	}

int MBIC_seq(bin_list *head, int level)
	{ rbtree tree;
	  rbtree_node tree_node;
	  bin_list *bin_min,*bin_tmp;
	  int i,k;

	  tree = rbtree_ini(head,level); /*initialize the tree*/

	  tree_node = min_tree_node(tree);/*the minimum tree node*/
          if(tree_node==NULL) return 0; /*empty tree, no bin merged*/
          bin_min = (bin_list *) tree_node->key;

	  //printf("\nstart merging\n");
	  //printf("\n");

	  k =0;
          while(bin_min->bic_diff<0)
		{ 
		  /*before merging the bins, we need to delete the tree nodes corresponding to these bins
  		    (also include the bins before the bin bin_min)*/
		   bin_tmp = bin_min; i= 0;
		   while(i<level&&bin_tmp!=NULL) /*delete nodes corresponding to the bins after bin_min (inclusive)*/
			{ //if(rbtree_lookup(tree, (void*) bin_tmp, bin_compare)==NULL)  printf("Node not found1.\n");
                  	  if(rbtree_delete(tree,(void *)bin_tmp,bin_compare)==1){
				bin_print(bin_tmp);
				}
			  bin_tmp = bin_tmp->next;
			  i++;
			}
		   bin_tmp = bin_min->pre;i=0;
		   while(i<level-1&&bin_tmp!=NULL) /*delete the nodes corresponding to the bins befor bin_min(exclusive)*/
                        { //if(rbtree_lookup(tree, (void*) bin_tmp, bin_compare)==NULL) printf("Node not found2.\n");
                          if(rbtree_delete(tree,(void *)bin_tmp,bin_compare)==1){
				bin_print(bin_tmp);
				}
                          bin_tmp = bin_tmp->pre;
                          i++;
                        }

		  bin_merge(bin_min,level); /*merge the bins and update the bic_diff of the bins*/

		  /*insert the nodes corresponding to the bins before bin_min(inclusive)*/
		  bin_tmp = bin_min; i = 0;
		  while(i<level&&bin_tmp!=NULL)
			{ rbtree_insert(tree,(void *)bin_tmp,(void *)bin_tmp,bin_compare);
			  bin_tmp = bin_tmp->pre; i++;
			}

		
		  tree_node = min_tree_node(tree);/*get the the minimum tree node*/
		  bin_min = (bin_list *) tree_node->key;
		  k++;

		  //printf("\n");
		}

	  rbtree_free(tree, head); /*free the memory used by tree.*/

	  return k;
	}
