#!/usr/bin/env perl
use strict;

use FindBin qw($Bin);
my $path = $Bin;

my $bicseq = "$path/MBIC-seq_v0.1.5/MBICseq";
my $combinefile = "$path/combineFile/combineFile";
my $rreport = "$path/R/report.R";
my $bicseqOneSample = "$path/BICseqOneSample.pl";
my $bidsegTwoSample = "$path/BICseqTwoSample.pl";
my $bidsegMulSample = "$path/BICseqMulSample.pl";
my $bootstrap = "$path/bootstrap/bootstrapTest";
my $plotSeg = "$path/R/plotProfile.R";

if(!(-e $plotSeg)) {
        die("Cannot find $plotSeg\n");
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
if(!(-e $bootstrap)){
	die("Cannot find $bootstrap\n");
	}


my $reportOnesample = "$path/R/reportOneSample.R";
my $rcombineBoot = "$path/R/combineSegBoostrap.R";
if(!(-e "$reportOnesample")){
        die("Cannot find $reportOnesample\n");
        }
if(!(-e $rcombineBoot)){
        die("Cannot find $rcombineBoot\n");
        }

if(!(-e $bicseqOneSample)){
	die("Cannot find $bicseqOneSample\n");
	}

if(!(-e $bidsegTwoSample)){
	die("Cannot find $bidsegTwoSample\n");
	}


my $help;
my $lambda=2;
my $configfile;
my $output;
my $tmpdir="$path/tmp/";
my $figname="";
my $plot = 0;
my $figtitle="";
my $nrm = "";
my $permutation = "";
my $noscale;
my $strict;
my $control;
my $detail_MS;

use Getopt::Long;

my $invalid = GetOptions("help"=>\$help,"lambda=f"=>\$lambda,"tmp=s"=>\$tmpdir,"fig=s"=>\$figname,"title=s"=>\$figtitle, "nrm"=>\$nrm, "bootstrap"=>\$permutation, "noscale"=>\$noscale, "strict"=>\$strict,"control"=>\$control, "detail"=>\$detail_MS);

my $size = $#ARGV+1;

if($help|!$invalid||$size!=2||$lambda<=0) {
	print "Usage: NBICseq-seg.pl [options] <configFile> <output>\n";
	print "Options:\n";
	print "        --lambda=<float>: the (positive) penalty used for BIC-seq\n";
	print "        --tmp=<string>: the tmp directory; If unspecified, use $path/tmp/\n";
	print "        --help: pring this message\n";
	print "        --fig=<string>: plot the CNV profile in a png file\n";
	print "        --title=<string>: the title of the figure\n";
	print "        --nrm: do not remove likely germline CNVs (with a matched normal) or segments with bad mappability (without a matched normal)\n";
	print "        --bootstrap: perform bootstrap test to assign confidence (only for one sample case)\n";
	print "        --noscale: do not automatically adjust the lambda parameter according to the noise level in the data\n";
	print "        --strict: if specified, use a more stringent method to ajust the lambda parameter\n";
	print "        --control: the data has a control genome\n";
	print "        --detail: if specified, print the detailed segmentation result (for multiSample only)\n";
	die("\n");
	}

$configfile = $ARGV[0];
$output = $ARGV[1];

if($noscale){
	$noscale = "--noscale";
	}else{
	$noscale = "";
	}
if($strict){
	$strict = "--strict";
	}else{
	$strict = "";
	}

if($detail_MS){
	$detail_MS = "--detail";
	}else{
	$detail_MS = "";
	}

if($tmpdir !~/\//){$tmpdir = $tmpdir."\/";}
if(!(-e $tmpdir)){mkdir $tmpdir or die "Cannot create the dirctory $tmpdir\n";}

my $figoption = "";
if($figname){
	$plot = 1;
	open(TEST,">$figname") or die("Cannot open $figname\n");
	close(TEST);
	unlink $figname or die("Cannot write to $figname\n");
	$figoption = "--fig=\"$figname\"";
	}
my $titleopt = "";
if($figtitle){
	$titleopt = "--title=\"$figtitle\"";
	}

if($nrm){ $nrm = "--nrm";
	} else{
	  $nrm = "";
	}
my $permutationOpt = "";
if($permutation){
	$permutationOpt = "-bootstrap";
	}


open(OUTPUT,">$output") or die("Cannot open the file: $output\n");
close(OUTPUT);


open(CONFIG, "<$configfile") or die("No such file or directory: $configfile\n");

use File::Temp qw/tempfile/;


#### first check if files exist
my $i=0;
my %chrom_names = {};
my @chroms_array = ();
my $columNum_FirstLine = 0;
my $flag = 0;
while(<CONFIG>){
	chomp;
	my @row = split(/\t/);
	my $num_items = $#row+1;
	if($num_items<2){die("Config file must have at leaste 2 columns\n");}
	if($flag==0){
		$columNum_FirstLine = $num_items;
		$flag = 1;
		}else{
		if($columNum_FirstLine != $num_items){
			my $rownum_tmp = $i+1;
			die("Error in config file: there are $num_items columns in the 1st row but $columNum_FirstLine columns in the $rownum_tmp"."th row\n");
			}
		}
	$i = $i+1;
	}

close(CONFIG);

my $cmd;
if($control){
	$cmd = "$bidsegTwoSample $noscale $strict --tmp=$tmpdir --lambda=$lambda $nrm $figoption $titleopt $permutationOpt $configfile $output";
	}else{
	if($columNum_FirstLine==2){
		$cmd = "$bicseqOneSample $noscale $strict --tmp=$tmpdir --lambda=$lambda $nrm $figoption $titleopt $permutationOpt  $configfile $output";
		}else{
		$cmd = "$bidsegMulSample $noscale $strict --tmp=$tmpdir --lambda=$lambda $nrm $figoption $titleopt $permutationOpt $detail_MS $configfile $output";
		}
	}

print "$cmd\n";
if(system($cmd)!=0){
	die("\n");
	}

