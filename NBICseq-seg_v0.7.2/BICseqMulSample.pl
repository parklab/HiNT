#!/usr/bin/env perl
use strict;

use FindBin qw($Bin);
my $path = $Bin;

my $bicseq = "$path/MBIC-seq_v0.1.5/MBICseq";
my $combinefile = "$path/combineFile/combineFile";
my $rreport = "$path/R/reportOneSample.R";
my $estimatLambdaFct = "$path/EstimateLambdaFactor/EstLambdaFct";
my $bootstrap = "$path/bootstrap/bootstrapTest";
my $rcombineBoot = "$path/R/combineSegBoostrap.R";
my $plotSeg = "$path/R/plotProfile.R";

if(!(-e $plotSeg)) {
	die("Cannot find $plotSeg\n");
	}

if(!(-e $rcombineBoot)){
        die("Cannot find $rcombineBoot\n");
        }


if(!(-e $bicseq)){
	die("Cannot find $bicseq\n");
	}
if(!(-e $combinefile)){
	die("Cannot find $combinefile\n");
	}
if(!(-e $rreport)){
	die("Cannot find $rreport\n");
	}


my $help;
my $lambda=2;
my $configfile;
my $output;
my $tmpdir="$path/tmp/";
my $figname="";
my $plot = 0;
my $figtitle="";
my $rm = "Y";
my $nrm = "";
my $noscale;
my $permutation;
my $strict;
my $detail_MS;

use Getopt::Long;

my $invalid = GetOptions("help"=>\$help,"lambda=f"=>\$lambda,"tmp=s"=>\$tmpdir,"fig=s"=>\$figname,"title=s"=>\$figtitle, "nrm"=>\$nrm, "noscale"=>\$noscale, "bootstrap"=>\$permutation, "strict"=>\$strict,"detail"=>\$detail_MS);

my $size = $#ARGV+1;

if($help|!$invalid||$size!=2||$lambda<=0) {
	print "BICseqMulSample.pl [options] <config file> <output>\n";
	print "Options:\n";
	print "        --lambda=<float>: the (positive) penalty used for BIC-seq\n";
	print "        --tmp=<string>: the tmp directory; If unspecified, use $path/tmp/\n";
	print "        --help: pring this message\n";
	print "        --fig=<string>: plot the CNV profile in a png file\n";
	print "        --title=<string>: the title of the figure\n";
	print "        --bootstrap: perform bootstrap test to assign confidence\n";
	print "        --noscale: do not automatically adjust the lambda parameter according to the noise level in the data\n";
	print "        --strict: if specified, use a more stringent method to ajust the lambda parameter\n";
	print "        --detail: print the detailed information of the segmentaion result\n";
	die("\n");
	}

$configfile = $ARGV[0];
$output = $ARGV[1];

if($tmpdir !~/\//){$tmpdir = $tmpdir."\/";}
if(!(-e $tmpdir)){mkdir $tmpdir or die "Cannot create the dirctory $tmpdir\n";}

if($figname){
	$plot = 1;
	open(TEST,">$figname") or die("Cannot open $figname\n");
	close(TEST);
	unlink $figname or die("Cannot write to $figname\n");
	}

if($nrm){ $rm = "N";
        } else{
          $rm = "Y";
        }


open(OUTPUT,">$output") or die("Cannot open the file: $output\n");
close(OUTPUT);


open(CONFIG, "<$configfile") or die("No such file or directory: $configfile\n");

use File::Temp qw/tempfile/;


#### first check if files exist
my $i=0;
my %chrom_names = {};
my @chroms_array = ();
#my $cmd_estLambdaFct = "$estimatLambdaFct ";
my @cmd_estLambdaFct = ();
my $num_Sample = 1;
while(<CONFIG>){
	chomp;
	my @row = split(/\t/);
	my $num_items = $#row+1;
	if($num_items < 2){die("Config file must have at least of 2 columns\n");}
	if($i==0){
		$num_Sample = $num_items-1;
		}else{
		my $rownum_tmp = $i+1;
		if($num_Sample!=$num_items-1) {die("Error in config file: row $rownum_tmp has different number of columns than row 1\n");}
		}
	my @Indx_BinFile = (1 .. $num_items-1);

	if($i>0){
		my $chrom = $row[0];
		foreach my $indx (@Indx_BinFile){
			my $binFile = $row[$indx];
			if(!(-e $binFile)){
				die("No such file or directory $binFile\n");
				}
			$cmd_estLambdaFct[$indx-1] = $cmd_estLambdaFct[$indx-1]." $binFile";
			}
		if(exists $chrom_names{$chrom}){
			die("Error in the configure file: two chromosome $chrom\n");
			}else{
			$chrom_names{$chrom} = "";
			push(@chroms_array, $chrom);
			}
		}else{
		foreach my $indx (@Indx_BinFile){
			push(@cmd_estLambdaFct, "$estimatLambdaFct ");
			}
		}
	$i = $i+1;
	}

close(CONFIG);



#### rescale the lambda paramter if required
use List::Util qw(max);
my @lambda_scale_theta = ();
my @lambda_scale_varMeanRatio = ();
if(!$noscale){
	foreach my $cmd_estLF (@cmd_estLambdaFct){
	        my $tmp_lambda_factor;
	        (undef,$tmp_lambda_factor) = tempfile("lambdaFactor"."_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);
	        print "$cmd_estLF > $tmp_lambda_factor\n";
	        system("$cmd_estLF > $tmp_lambda_factor");
	        open(IN, "<$tmp_lambda_factor") or die("No such file or directory: $tmp_lambda_factor\n");
        	my $i = 0;
	        while(<IN>){
        	        chomp;
                	my @row = split(/\s+/);
	                my $num_items = $#row+1;

        	        if($i==1){
                	        push(@lambda_scale_theta,$row[0]);
                        	push(@lambda_scale_varMeanRatio, $row[1]);
	                        }
        	        $i = $i+1;
                	}
	        close(IN);
		if(-e $tmp_lambda_factor) {unlink "$tmp_lambda_factor" or print "Failed to delete $tmp_lambda_factor\n";}
		}

	my $lambda_scale_theta_max = max @lambda_scale_theta;
        if($lambda_scale_theta_max>0){
		if($strict){
			$lambda = ($lambda_scale_theta_max+1)*($lambda_scale_theta_max+1)*$lambda;
			}else{
			$lambda = ($lambda_scale_theta_max+1)*$lambda;
			}
                }

	}

### then combine the data and run BICseq

open(CONFIG, "<$configfile") or die("No such file or directory: $configfile\n");
my $tmp_readdata;
(undef,$tmp_readdata) = tempfile("readdata_all"."_XXXXXXX",SUFFIX=>".bin",DIR => $tmpdir,OPEN=>0);
my $flag = 0;
my $i=0;
my %tmp_seg_files = {};
while(<CONFIG>){
        chomp;
        my @row = split(/\t/);
        my $num_items = $#row+1;
        if($num_items < 2){die("Config file must have at least 2 columns\n");}
	my @Indx_BinFile = (1 .. $num_items-1);
        if($i>0){
                my $chrom = $row[0];
		my $cmd_combinefile = $combinefile;
                foreach my $indx (@Indx_BinFile){
                        my $binFile = $row[$indx];
                        if(!(-e $binFile)){
                                die("No such file or directory $binFile\n");
                                }
                        $cmd_combinefile = $cmd_combinefile." $binFile";
                        }

		my $tmpfile;
		(undef, $tmpfile) = tempfile("seg_".$chrom."_XXXXXXX",SUFFIX=>".mbic",DIR => $tmpdir,OPEN=>0);
		if(!(exists $tmp_seg_files{$chrom})){
			$tmp_seg_files{$chrom} = $tmpfile;
			}
		my $cmd = "$cmd_combinefile | $bicseq -l $lambda > $tmpfile";
		print "$cmd\n";
		if(system($cmd) !=0 ) {die("\n");};
		if(!(-e $tmpfile)){
			die("Error: failed to create the temp file: $tmpfile\n");
			}

                if($permutation){
                        my $cmd1;
                        if($flag==0) {
                                $cmd1 = "$cmd_combinefile | cut -f1-2 --complement > $tmp_readdata";
                                $flag = 1;
                                }else{
				$cmd1 = "$cmd_combinefile | cut -f1-2 --complement >> $tmp_readdata";
                                }
			print "$cmd1\n";
                        if(system($cmd1)!=0){die("\n");}
                        if(!(-e $tmp_readdata)){
                                die("Error: failed to create the temp file $tmp_readdata\n");
                                }
                        }

                }
        $i = $i+1;
        }

close(CONFIG);



### calculate some parameters used for reformated data
my $total_bin = 0;
my $tumor_total = 0;
my $normal_total = 0;

my $combined_mbicfile;
(undef, $combined_mbicfile) = tempfile("combined.seg.1st._XXXXXXX",SUFFIX=>".mbic",DIR => $tmpdir,OPEN=>0);
open(OUTPUT1,">$combined_mbicfile") or die("Cannot create the file: $combined_mbicfile\n");
print OUTPUT1 "chrom\tstart\tend\tbinNum";
my $index_sample = 1;
while ($index_sample <= $num_Sample){
	print OUTPUT1 "\tobserved_$index_sample\texpected_$index_sample";
	$index_sample = $index_sample + 1;
	}
print OUTPUT1 "\n";
#foreach my $key (keys %tmp_seg_files){
foreach my $key (@chroms_array){ ## here key is actually the chromosome names
	if($tmp_seg_files{$key}){
		open(SEG_FILE,"<$tmp_seg_files{$key}") or die("Error: cannot open the temp file: $tmp_seg_files{$key}\n");

		my $i = 0;
		while(<SEG_FILE>){
			chomp;
			my @row = split(/\t/);
			my $num_item = $#row + 1;
			if($i>0 && $num_item== 2*$num_Sample + 3){
				print OUTPUT1 "$key\t$_\n";
				$total_bin = $total_bin +  $row[2];
				}
			$i = $i+1;
			}
		close(SEG_FILE);
		}
	}


close(OUTPUT1);

print "$combined_mbicfile\n";

my $output_detail;
(undef, $output_detail) = tempfile("detailed.seg._XXXXXXX",SUFFIX=>".seg",DIR => $tmpdir,OPEN=>0);
### report the final segmentation
my $cmd = "R --slave --args $combined_mbicfile $rm $output_detail < $rreport";
print "$cmd\n";
if(system($cmd)!=0) {die("\n");};

if($detail_MS){
	my $cmd = "cat $output_detail > $output";
	print "$cmd\n";
	if(system($cmd)!=0) {die("\n");};
	}else{
	my $cmd = "cut -f 1-7 $output_detail > $output";
	print "$cmd\n";
	if(system($cmd)!=0) {die("\n");};
	}

#### perform a permuation test to assign confidence
#if($permutation){
#        my $seg_nochrom_tmp;
#        (undef, $seg_nochrom_tmp) = tempfile("seg_nochrom_XXXXXXX",SUFFIX=>".seg",DIR => $tmpdir,OPEN=>0);
#        my $cmd = "cut -f2-10 $output > $seg_nochrom_tmp";
#        if(system("$cmd")!=0){die("\n");};
#        my $cmd = "$bootstrap $seg_nochrom_tmp $tmp_readdata";
#        print "$cmd\n";
#        if(system($cmd)!=0) {die("\n")};
#        my $cmd = "R --slave --args $output $seg_nochrom_tmp < $rcombineBoot";
#        print "$cmd\n";
#        if(system($cmd)!=0) {die("\n");};

#        if(-e $seg_nochrom_tmp) {
#                unlink($seg_nochrom_tmp) or print "Failed to delete the file $seg_nochrom_tmp\n";
#                }
#        }


##### now plot the profile if required
if($plot!=0){
	$cmd = "R --slave --args $output $figname $figtitle < $plotSeg";
	print "$cmd\n";
	if(system($cmd)!=0) {print "Error while plotting the segmentation profile\n";}
	}

### remove the temp files
foreach my $key (keys %tmp_seg_files){
        if($tmp_seg_files{$key}){
		unlink($tmp_seg_files{$key}) or print "Failed to delete the file $tmp_seg_files{$key}\n";
                }
        }

unlink($output_detail) or print "Failed to delete the file $output_detail\n";
unlink($combined_mbicfile) or print "Failed to delete the file $combined_mbicfile\n";
if(-e $tmp_readdata){
        unlink($tmp_readdata) or print "Failed to delete the file $tmp_readdata\n";
        }

