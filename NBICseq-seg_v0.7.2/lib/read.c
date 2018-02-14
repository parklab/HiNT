/* This is the read data program written by Ruibin Xi Oct 12, 2009.
 * Modify the function read_table. The function read_table will ignore any lines that doesn't start from numbers
 * Dec 27 2009
 * Modified on Feb 1 2012 (In the earliear version, when the input file contains no numeric value, this program will freeze)
 * Modified on Mar 21 2012, fixed a bug in read_table (when lines = 1, read_table of the version on Feb 1 2012 will read 2 lines)
*/

#include "read.h"


static char *ll_buffer = NULL; /*a string that is used to store a line read from a file*/
static int buffer_size = LR_BUFFER; /*buffer_size+10 is The actual size of ll_buffer*/

static int is_number(char c);
/*test whether a character c is '0','1',..,'9' or '-' (minus sign)
 *  *if c is '0','1',..,'9' return 1
 *   *if c is '-' return -1
 *    *otherwise return 0
 *     * */

static int is_numeric(char *ll);
/*test if a character array starts with numeric values
 *return 1 if true, 0 otherwise
 * */


static int is_number(char c)
        { if(c<='9'&&c>='0') return 1;
          else if(c=='-') return -1;
          else return 0;
        }

static int is_numeric(char *ll)
	{ int len,flag;
	  len = strlen(ll);
	  
	  if(len==0) return 0;

	  flag = is_number(ll[0]);
	  if(flag==0) return 0;
	  else if(flag==-1)
		{ if(len==1||is_number(ll[1])!=1) return 0;
		  else return 1;
		}	
	  else return 1;
	}



char * lr_read(FILE *input)
{	char chr;
	int i=0;
	
	if(input == NULL) return NULL; 
	if(ll_buffer==NULL) ll_buffer = (char *) malloc(sizeof(char)*(buffer_size+10));
	if(ll_buffer==NULL) return NULL;

	do
		{ chr = fgetc(input);
		  ll_buffer[i] = chr; ll_buffer[i+1] = '\0'; i++;
		  if(strlen(ll_buffer)>= buffer_size) 
			{ buffer_size = buffer_size + LR_BUFFER;
			  ll_buffer = (char *)  realloc(ll_buffer,sizeof(char)*(buffer_size+10));
			  if(ll_buffer==NULL) return NULL;
			}
		}while(chr!='\n'&&chr!=EOF);

	return ll_buffer;
}



double * read_table(FILE *input, int *nrow, int *ncol, int lines,int skip)
{	double * mtrx;
	int i=0,j=0,flag=1,low_mmry = 40000, is_num=0,k;
	char *ll,*sub_str;
	//char chr;
	double tmp;
	
        (*nrow) = 0;
        (*ncol) = 2;
	for(i=0;i<skip;i++)
     		{ ll = lr_read(input);
		  if(strlen(ll)==0) return NULL;
		  //printf("%s\n",ll);
		}
	
        i=0;j=0;flag=1;k = skip;
	mtrx = (double *) malloc(sizeof(double)*(low_mmry+10));
	if(mtrx==NULL){
		fprintf(stderr,"Error: memory allocation failed\n");
		exit(1);
		}

	if(lines==0) return mtrx;

	
	ll = lr_read(input); k++; /*k indicate which line is under consideration*/
	if(strlen(ll)==0) return NULL;
	sub_str = strtok(ll," \t");
	is_num = is_numeric(sub_str);
	while(is_num==0) /*The current line is not numeric, read the next line*/
		{ ll = lr_read(input); k++;
		  if(strlen(ll)==0) return NULL;
	          sub_str = strtok(ll," \t");
       		  is_num = is_numeric(sub_str);
		if(feof(input)) return NULL;
		}

	while(sub_str!=NULL)
		{ 
		  tmp = atof(sub_str);
		  mtrx[j] = tmp;
		  j++;
		  sub_str=strtok(NULL," \t");
		  if(j>=low_mmry)
			 { low_mmry += 4000; 
			   mtrx = (double *) realloc(mtrx,sizeof(double)*(low_mmry+10));
			 }
		  if(sub_str!=NULL && is_numeric(sub_str)==0)
			{ fprintf(stderr,"Error in read_table: expecting numeric value, got \"%s\" (row %d)\n",sub_str,k);
			  exit(1);
			}
		}
	i++;
	(*ncol) = j;

	if(i< lines || lines < 0){ /*lines < 0 means read all data */
		while(!feof(input))
			{ ll = lr_read(input); k++;		  
	 		  sub_str = strtok(ll," \t");
	        	  is_num = is_numeric(sub_str);
		  
			  if(is_num==1)
				{ j = 0;
			  	  while(sub_str!=NULL&&j<=(*ncol))
					{ tmp = atof(sub_str);
					  mtrx[i*(*ncol)+j] = tmp; 
					  j++;
					  sub_str=strtok(NULL," \t");
        	          		  if(sub_str!=NULL && is_numeric(sub_str)==0)
                	    			    { fprintf(stderr,"Error in read_table: expecting numeric value, got \"%s\" (row %d)\n",sub_str,k);
                       			  	      exit(1);
                      				    }
					}

				  i++;
				  if(j<(*ncol)||j>(*ncol))  i--; /*skip any line that has different number of columns from the first line*/
				  if((i+1)*(*ncol)>=low_mmry)
					{ low_mmry += 200*(*ncol); /*not enough space, create more*/
					  mtrx = (double *) realloc(mtrx,sizeof(double)*(low_mmry+10));
					}
		
				  if(lines>0&&i>=lines) break; /*alread read lines lines, terminate the while loop*/
				}
			}
		}
	(*nrow) = i;
	return mtrx;
}


