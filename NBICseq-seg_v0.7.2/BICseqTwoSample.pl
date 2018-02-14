#!/usr/bin/env perl
use strict;

use FindBin qw($Bin);
my $path = $Bin;

my $bicseq = "$path/MBIC-seq_v0.1.5/MBICseq";
my $combinefile = "$path/combineFile/combineFile";
my $rreport = "$path/R/report.R";
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
my $rm = "--rm";
my $nrm = "";
my $noscale;
my $permutation;
my $strict;

use Getopt::Long;

my $invalid = GetOptions("help"=>\$help,"lambda=f"=>\$lambda,"tmp=s"=>\$tmpdir,"fig=s"=>\$figname,"title=s"=>\$figtitle, "nrm"=>\$nrm, "noscale"=>\$noscale, "bootstrap"=>\$permutation, "strict"=>\$strict);

my $size = $#ARGV+1;

if($help|!$invalid||$size!=2||$lambda<=0) {
	print "BICseqTwoSample.pl [options] <config file> <output>\n";
	print "Options:\n";
	print "        --lambda=<float>: the (positive) penalty used for BIC-seq\n";
	print "        --tmp=<string>: the tmp directory; If unspecified, use $path/tmp/\n";
	print "        --help: pring this message\n";
	print "        --fig=<string>: plot the CNV profile in a png file\n";
	print "        --title=<string>: the title of the figure\n";
	print "        --nrm: do not remove likely germline CNVs\n";
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

if($nrm){ $rm = "";
	}

open(OUTPUT,">$output") or die("Cannot open the file: $output\n");
close(OUTPUT);


open(CONFIG, "<$configfile") or die("No such file or directory: $configfile\n");

use File::Temp qw/tempfile/;


#### first check if files exist
my $i=0;
my %chrom_names = {};
my @chroms_array = ();
my $cmd_estLambdaFct_tumor = "$estimatLambdaFct ";
my $cmd_estLambdaFct_normal = "$estimatLambdaFct ";
while(<CONFIG>){
	chomp;
	my @row = split(/\t/);
	my $num_items = $#row+1;
	if($num_items!=3){die("Config file must be 3 column file\n");}
	if($i>0){
		my $chrom = $row[0];
		my $tumor = $row[1];
		my $normal = $row[2];
		if(!(-e $tumor)){
			die("No such file or directory $tumor\n");
			}
		if(!(-e $normal)){
			die("No such file or directory $normal\n");
			}
		if(exists $chrom_names{$chrom}){
			die("Error in the configure file: two chromosome $chrom\n");
			}else{
			$chrom_names{$chrom} = "";
			push(@chroms_array, $chrom);
			}
		$cmd_estLambdaFct_tumor = $cmd_estLambdaFct_tumor." $tumor";
		$cmd_estLambdaFct_normal = $cmd_estLambdaFct_normal." $normal";
		}
	$i = $i+1;
	}

close(CONFIG);


#### rescale the lambda paramter if required
if(!$noscale){
        my $tmp_lambda_factor;
        (undef,$tmp_lambda_factor) = tempfile("lambdaFactor"."_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);
        print "$cmd_estLambdaFct_tumor > $tmp_lambda_factor\n";
        system("$cmd_estLambdaFct_tumor > $tmp_lambda_factor");
        open(IN, "<$tmp_lambda_factor") or die("No such file or directory: $tmp_lambda_factor\n");
        my $lambda_scale_theta_tumor;
        my $lambda_scale_varMeanRatio_tumor;
        my $i = 0;
        while(<IN>){
                chomp;
                my @row = split(/\s+/);
                my $num_items = $#row+1;

                if($i==1){
                        $lambda_scale_theta_tumor = $row[0];
                        $lambda_scale_varMeanRatio_tumor = $row[1];
                        }
                $i = $i+1;
                }
        close(IN);
	if(-e $tmp_lambda_factor) {unlink "$tmp_lambda_factor" or print "Failed to delete $tmp_lambda_factor\n";}

	print "$cmd_estLambdaFct_normal > $tmp_lambda_factor\n";
        system("$cmd_estLambdaFct_normal > $tmp_lambda_factor");
        open(IN, "<$tmp_lambda_factor") or die("No such file or directory: $tmp_lambda_factor\n");
        my $lambda_scale_theta_normal;
        my $lambda_scale_varMeanRatio_normal;
        my $i = 0;
        while(<IN>){
                chomp;
                my @row = split(/\s+/);
                my $num_items = $#row+1;
        
                if($i==1){
                        $lambda_scale_theta_normal = $row[0];
                        $lambda_scale_varMeanRatio_normal = $row[1];
                        }
                $i = $i+1;
                }
        close(IN);
	if(-e $tmp_lambda_factor) {unlink "$tmp_lambda_factor" or print "Failed to delete $tmp_lambda_factor\n";}


	my $lambda_scale_theta = $lambda_scale_theta_tumor;
	if($lambda_scale_theta_tumor < $lambda_scale_theta_normal) {$lambda_scale_theta = $lambda_scale_theta_normal;}
        if($lambda_scale_theta>0){
		if($strict){
			$lambda = ($lambda_scale_theta+1)*($lambda_scale_theta+1)*$lambda;
			}else{
			$lambda = ($lambda_scale_theta+1)*$lambda;
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
        if($num_items!=3){die("Config file must be 3 column file\n");}
        if($i>0){
                my $chrom = $row[0];
                my $tumor = $row[1];
                my $normal = $row[2];
                if(!(-e $tumor)){
                        die("No such file or directory $tumor\n");
                        }
                if(!(-e $normal)){
                        die("No such file or directory $normal\n");
                        }

		my $tmpfile;
		(undef, $tmpfile) = tempfile("seg_".$chrom."_XXXXXXX",SUFFIX=>".mbic",DIR => $tmpdir,OPEN=>0);
		if(!(exists $tmp_seg_files{$chrom})){
			$tmp_seg_files{$chrom} = $tmpfile;
			}
		my $lambda_1st = $lambda/2.0;
		my $cmd = "$combinefile $tumor $normal | $bicseq -l $lambda_1st $rm > $tmpfile";
		print "$cmd\n";
		if(system($cmd) !=0 ) {die("\n");};
		if(!(-e $tmpfile)){
			die("Error: failed to create the temp file: $tmpfile\n");
			}

                if($permutation){
                        my $cmd1;
                        if($flag==0) {
                                $cmd1 = "$combinefile $tumor $normal | cut -f3-6 > $tmp_readdata";
                                $flag = 1;
                                }else{
                                $cmd1 = "$combinefile $tumor $normal | cut -f3-6 >> $tmp_readdata";
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
print OUTPUT1 "chrom\tstart\tend\tbinNum\ttumor\ttumor_expected\tnormal\tnormal_expected\n";
#foreach my $key (keys %tmp_seg_files){
foreach my $key (@chroms_array){ ## here key is actually the chromosome names
	if($tmp_seg_files{$key}){
		open(SEG_FILE,"<$tmp_seg_files{$key}") or die("Error: cannot open the temp file: $tmp_seg_files{$key}\n");

		my $i = 0;
		while(<SEG_FILE>){
			chomp;
			my @row = split(/\t/);
			my $num_item = $#row + 1;
			if($i>0 && $num_item==7){
				print OUTPUT1 "$key\t$_\n";
				$total_bin = $total_bin +  $row[2];
				$tumor_total = $tumor_total + $row[3];
				$normal_total = $normal_total + $row[5];
				}
			$i = $i+1;
			}
		close(SEG_FILE);
		}
	}


close(OUTPUT1);

### reformate the output from the 1st step of MBIC-seq
my $combined_mbicfile2;
(undef, $combined_mbicfile2) = tempfile("combined.seg.2nd._XXXXXXX",SUFFIX=>".mbic",DIR => $tmpdir,OPEN=>0);

open(OUTPUT2,">$combined_mbicfile2") or die("Cannot create the file: $combined_mbicfile2\n");
print OUTPUT2 "chrom\tstart\tend\tbinNum\ttumor\tnormal\n";

#foreach my $key (keys %tmp_seg_files){
foreach my $key (@chroms_array){
	if($tmp_seg_files{$key}){
	        open(SEG_FILE,"<$tmp_seg_files{$key}") or die("Error: cannot open the temp file: $tmp_seg_files{$key}\n");
		my $tmpfile;
		(undef, $tmpfile) = tempfile("tmp_XXXXXXX",SUFFIX=>".bin",DIR => $tmpdir,OPEN=>0);

		open(TEMP, ">$tmpfile") or die("Faile to open the file $tmpfile\n");
		my $i = 0;
	        while(<SEG_FILE>){
        	        chomp;
              		my @row = split(/\t/);
                	my $num_item = $#row + 1;
                	if($i>0 && $num_item==7){
				my $binNum = $row[2];
				my $weight = $binNum*($tumor_total/$total_bin);
				my $tumor_normalized = ($row[3]/$row[4])*$weight;
				my $weight = $binNum*($normal_total/$total_bin);
				my $normal_normalized = ($row[5]/$row[6])*$weight;
				print TEMP "$row[0]\t$row[1]\t$tumor_normalized\t$normal_normalized\n";
                        	}
			$i = $i +1;
                	}
        	close(SEG_FILE);
		close(TEMP);

		my $tmpfile2;
		(undef, $tmpfile2) = tempfile("seg_".$key."\.2nd"."_XXXXXXX",SUFFIX=>".bin",DIR => $tmpdir,OPEN=>0);
		my $cmd = "cat $tmpfile | $bicseq -l $lambda > $tmpfile2";
		if(system($cmd)!=0) {die("\n");};

		$i = 0;
		open(TEMP,"<$tmpfile2") or die("Error, Failed to open the temp file: $tmpfile2\n");
		while(<TEMP>){
			chomp;
			if($i>0){
				print OUTPUT2 "$key\t$_\n";
				}
			$i = $i + 1;
			}
	
		close(TEMP);

		unlink($tmpfile) or print("Failed to delete the file $tmpfile\n");
		unlink($tmpfile2) or print("Failed to delete the file $tmpfile2\n");
		}
        }

close(OUTPUT2);


### report the final segmentation
my $cmd = "R --slave --args $combined_mbicfile $combined_mbicfile2 $output < $rreport";
print "$cmd\n";
if(system($cmd)!=0) {die("\n");};


#### perform a permuation test to assign confidence
if($permutation){
        my $seg_nochrom_tmp;
        (undef, $seg_nochrom_tmp) = tempfile("seg_nochrom_XXXXXXX",SUFFIX=>".seg",DIR => $tmpdir,OPEN=>0);
        my $cmd = "cut -f2-10 $output > $seg_nochrom_tmp";
        if(system("$cmd")!=0){die("\n");};
        my $cmd = "$bootstrap $seg_nochrom_tmp $tmp_readdata";
        print "$cmd\n";
        if(system($cmd)!=0) {die("\n")};
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
	print "$cmd\n";
	if(system($cmd)!=0) {print "Error while plotting the segmentation profile\n";}
	}

### remove the temp files
foreach my $key (keys %tmp_seg_files){
        if($tmp_seg_files{$key}){
		unlink($tmp_seg_files{$key}) or print "Failed to delete the file $tmp_seg_files{$key}\n";
                }
        }

unlink($combined_mbicfile) or print "Failed to delete the file $combined_mbicfile\n";
unlink($combined_mbicfile2) or print "Failed to delete the file $combined_mbicfile2\n";
if(-e $tmp_readdata){
        unlink($tmp_readdata) or print "Failed to delete the file $tmp_readdata\n";
        }

