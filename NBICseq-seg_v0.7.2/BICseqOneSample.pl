#!/usr/bin/env perl
use strict;

use FindBin qw($Bin);
my $path = $Bin;

my $bicseq = "$path/MBIC-seq_v0.1.5/MBICseq";
my $combinefile = "$path/combineFile/combineFile";
my $report = "$path/R/reportOneSample.R";
my $bootstrap = "$path/bootstrap/bootstrapTest";
my $estimatLambdaFct = "$path/EstimateLambdaFactor/EstLambdaFct";

my $plotSeg = "$path/R/plotProfile.R";

if(!(-e $plotSeg)) {
        die("Cannot find $plotSeg\n");
        }


if(!(-e $estimatLambdaFct)){
	die("Cannot find $estimatLambdaFct\n");
	}

if(!(-e $bicseq)){
	die("Cannot find $bicseq\n");
	}
if(!(-e $combinefile)){
	die("Cannot find $combinefile\n");
	}
if(!(-e "$report")){
	die("Cannot find $report\n");
	}
if(!(-e $bootstrap)){
        die("Cannot find $bootstrap\n");
        }

my $rcombineBoot = "$path/R/combineSegBoostrap.R";
if(!(-e $rcombineBoot)){
        die("Cannot find $rcombineBoot\n");
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
my $permutation;
my $noscale;
my $strict;

use Getopt::Long;

my $invalid = GetOptions("help"=>\$help,"lambda=f"=>\$lambda,"tmp=s"=>\$tmpdir,"fig=s"=>\$figname,"title=s"=>\$figtitle, "nrm"=>\$nrm, "bootstrap"=>\$permutation, "noscale"=>\$noscale, "strict"=>\$strict);

my $size = $#ARGV+1;

if($help|!$invalid||$size!=2||$lambda<=0) {
	print "BICseqOneSample.pl [options] <config file> <output>\n";
	print "Options:\n";
	print "        --lambda=<float>: the (positive) penalty used for BIC-seq\n";
	print "        --tmp=<string>: the tmp directory; If unspecified, use $path/tmp/\n";
	print "        --help: pring this message\n";
	print "        --fig=<string>: plot the CNV profile in a png file\n";
	print "        --title=<string>: the title of the figure\n";
	print "        --nrm: Do not remove segments with bad mappability\n";
	print "        --bootstrap: perform bootstrap test to assign confidence\n";
	print "        --noscale: do not automatically adjust the lambda parameter according to the noise level in the data\n";
	print "        --strict: if specified, use a more stringent method to ajust the lambda parameter\n";
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
my $cmd_estLambdaFct = "$estimatLambdaFct ";
while(<CONFIG>){
	chomp;
	my @row = split(/\t/);
	my $num_items = $#row+1;
	if($num_items!=2){die("Config file must be 2 column file\n");}
	if($i>0){
		my $chrom = $row[0];
		my $tumor = $row[1];
		if(!(-e $tumor)){
			die("No such file or directory $tumor\n");
			}
		if(exists $chrom_names{$chrom}){
			die("Error in the configure file: two chromosome $chrom\n");
			}else{
			$chrom_names{$chrom} = "";
			push(@chroms_array, $chrom);
			}
		$cmd_estLambdaFct = $cmd_estLambdaFct." $tumor";
		}
	$i = $i+1;
	}

close(CONFIG);


### automatically adjust the lambda parameter if required
if(!$noscale){
	my $tmp_lambda_factor;
	(undef,$tmp_lambda_factor) = tempfile("lambdaFactor"."_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);
	print "$cmd_estLambdaFct\n";
	system("$cmd_estLambdaFct > $tmp_lambda_factor");
	open(IN, "<$tmp_lambda_factor") or die("No such file or directory: $tmp_lambda_factor\n");
	my $lambda_scale_theta;
	my $lambda_scale_varMeanRatio;
	my $i = 0;
	while(<IN>){
		chomp;
		my @row = split(/\s+/);
		my $num_items = $#row+1;

		if($i==1){
			$lambda_scale_theta = $row[0];
			$lambda_scale_varMeanRatio = $row[1];
			}
		$i = $i+1;
		}
	close(IN);
	if($lambda_scale_theta>0){
		if($strict){
        		$lambda = ($lambda_scale_theta+1)*($lambda_scale_theta+1)*$lambda;
			}else{
			$lambda = ($lambda_scale_theta+1)*$lambda;
			}
        	}
	if(-e $tmp_lambda_factor){
		unlink $tmp_lambda_factor or print "Failed to delete $tmp_lambda_factor\n";
		}
	}

### then combine the data and run BICseq

open(CONFIG, "<$configfile") or die("No such file or directory: $configfile\n");

my $i=0;
my %tmp_seg_files = {};
my $flag = 0;
my $tmp_readdata;
(undef,$tmp_readdata) = tempfile("readdata_all"."_XXXXXXX",SUFFIX=>".bin",DIR => $tmpdir,OPEN=>0);
while(<CONFIG>){
        chomp;
        my @row = split(/\t/);
        my $num_items = $#row+1;
        if($num_items!=2){die("Config file must be 3 column file\n");}
        if($i>0){
                my $chrom = $row[0];
                my $tumor = $row[1];
                if(!(-e $tumor)){
                        die("No such file or directory $tumor\n");
                        }

		my $tmpfile;
		(undef, $tmpfile) = tempfile("seg_".$chrom."_XXXXXXX",SUFFIX=>".mbic",DIR => $tmpdir,OPEN=>0);
		if(!(exists $tmp_seg_files{$chrom})){
			$tmp_seg_files{$chrom} = $tmpfile;
			}
		my $cmd = "$combinefile $tumor | $bicseq -l $lambda > $tmpfile";
		print "$cmd\n";
		if(system($cmd) !=0 ) {die("\n");};
		if(!(-e $tmpfile)){
			die("Error: failed to create the temp file: $tmpfile\n");
			}

		if($permutation){
			my $cmd1;
			if($flag==0) {
				$cmd1 = "cut -f3,4 $tumor > $tmp_readdata";
				$flag = 1;
				}else{
				$cmd1 = "cut -f3,4 $tumor | tail -n +2 >> $tmp_readdata";
				}
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
my $observed_total = 0;
my $expected_total = 0;

my $combined_mbicfile;
(undef, $combined_mbicfile) = tempfile("combined.seg_XXXXXXX",SUFFIX=>".mbic",DIR => $tmpdir,OPEN=>0);
open(OUTPUT1,">$combined_mbicfile") or die("Cannot create the file: $combined_mbicfile\n");
print OUTPUT1 "chrom\tstart\tend\tbinNum\tobserved\texpected\n";
#foreach my $key (keys %tmp_seg_files){
foreach my $key (@chroms_array){ ## here key is actually the chromosome names
	if($tmp_seg_files{$key}){
		open(SEG_FILE,"<$tmp_seg_files{$key}") or die("Error: cannot open the temp file: $tmp_seg_files{$key}\n");

		my $i = 0;
		while(<SEG_FILE>){
			chomp;
			my @row = split(/\t/);
			my $num_item = $#row + 1;
			if($i>0 && $num_item==5){
				print OUTPUT1 "$key\t$_\n";
				$total_bin = $total_bin +  $row[2];
				$observed_total = $observed_total + $row[3];
				$expected_total = $expected_total + $row[4];
				}
			$i = $i+1;
			}
		close(SEG_FILE);
		}
	}


close(OUTPUT1);


### report the final segmentation
my $cmd = "R --slave --args $combined_mbicfile $rm $output < $report";
print "$cmd\n";
if(system($cmd)!=0) {die("\n");};

#### perform a permuation test to assign confidence
if($permutation){
	my $seg_nochrom_tmp;
	(undef, $seg_nochrom_tmp) = tempfile("seg_nochrom_XXXXXXX",SUFFIX=>".seg",DIR => $tmpdir,OPEN=>0);
	my $cmd = "cut -f2-7 $output > $seg_nochrom_tmp";
	if(system("$cmd")!=0) {die("\n")};
	my $cmd = "$bootstrap $seg_nochrom_tmp $tmp_readdata";
	print "$cmd\n";
	if(system($cmd)!=0) {die("\n");};
	my $cmd = "R --slave --args $output $seg_nochrom_tmp < $rcombineBoot";
	print "$cmd\n";
	if(system($cmd)!=0) {die("\n");};

	if(-e $seg_nochrom_tmp) {
		unlink($seg_nochrom_tmp) or print "Failed to delete the file $seg_nochrom_tmp\n";
		}
	}


##### now plot the profile if required
if($plot!=0){
        $cmd = "R --slave --args $output $figname $figtitle < $plotSeg";
        if(system($cmd)!=0) {print "Error while plotting the segmentation profile\n";}
        }


### remove the temp files
foreach my $key (keys %tmp_seg_files){
        if($tmp_seg_files{$key}){
		unlink($tmp_seg_files{$key}) or print "Failed to delete the file $tmp_seg_files{$key}\n";
                }
        }

unlink($combined_mbicfile) or print "Failed to delete the file $combined_mbicfile\n";

if(-e $tmp_readdata){
	unlink($tmp_readdata) or print "Failed to delete the file $tmp_readdata\n";
	}
