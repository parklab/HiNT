#!/usr/bin/env perl
use strict;

use FindBin qw($Bin);
my $path = $Bin;

my $combinefile = "$path/combineFile/combineFile";
my $countRead = "$path/genotype/countRead";
my $genotype = "$path/genotype/genotype";

if(!(-e $countRead)) {
	die("Cannot find $countRead\n");
	}

if(!(-e $genotype)){
        die("Cannot find $genotype\n");
        }

if(!(-e $combinefile)){
	die("Cannot find $combinefile\n");
	}


my $help;
my $configfile;
my $regionfile;
my $output;
my $tmpdir="$path/tmp/";

use Getopt::Long;

my $invalid = GetOptions("help"=>\$help,"tmp=s"=>\$tmpdir);

my $size = $#ARGV+1;

if($help|!$invalid||$size!=3) {
	print "genotype.pl [options] <configFile> <RegionFile> <Output>\n";
	print "Options:\n";
	print "        --tmp=<string>: the tmp directory; If unspecified, use $path/tmp/\n";
	print "        --help: pring this message\n";
	die("\n");
	}

$configfile = $ARGV[0];
$regionfile = $ARGV[1];
$output = $ARGV[2];

if($tmpdir !~/\//){$tmpdir = $tmpdir."\/";}
if(!(-e $tmpdir)){mkdir $tmpdir or die "Cannot create the dirctory $tmpdir\n";}

open(REGIONFILE, "<$regionfile") or die("Cannot open the file: $regionfile\n");
close(REGIONFILE);

open(OUTPUT,">$output") or die("Cannot open the file: $output\n");
close(OUTPUT);


use File::Temp qw/tempfile/;


#### first check if files exist
open(CONFIG, "<$configfile") or die("No such file or directory: $configfile\n");
my $i=0;
my %normalized_data = {};
my @chroms_array = ();
while(<CONFIG>){
	chomp;
	my @row = split(/\t/);
	my $num_items = $#row+1;
	if($num_items!=3 && $num_items!=2){die("Config file must be a 2 or 3 column file\n");}
	if($i>0){
		my $chrom = $row[0];
		my $tumor = $row[1];
		my $normal = "";
		if(!(-e $tumor)){
			die("No such file or directory $tumor\n");
			}
		if($num_items==3){
			$normal = $row[2];
			if(!(-e $normal)){
				die("No such file or directory $normal\n");
				}
			}
		if(exists $normalized_data{$chrom}){
			die("Error in the configure file: two chromosome $chrom\n");
			}else{
			$normalized_data{$chrom} = "-t $tumor";
			if($num_items==3) {$normalized_data{$chrom} = $normalized_data{$chrom}." -c $normal";}
			push(@chroms_array, $chrom);
			}
		}
	$i = $i+1;
	}

close(CONFIG);


##########  check the region file ############

my %regionsbychrom;
my @region_Chroms = ();
open(REGIONFILE, "<$regionfile") or die("Cannot open the file: $regionfile\n");
my $i = 0;
while(<REGIONFILE>){
	chomp;
	my @row = split(/\t/);
	my $num_items = $#row+1;
	if($num_items!=3){die("Region file must be a 3 column file\n");}
	my $chrom = $row[0];
	my $start = $row[1];
	my $end = $row[2];
	if($i > 0){ ## the first line is the header
		if(!(exists $normalized_data{$chrom})){
			die("No such chromosome in the config file: $chrom\n");
			}
		if(!exists $regionsbychrom{$chrom}){
			@{$regionsbychrom{$chrom}} = ("$start\t$end");
			}else{
			push(@{$regionsbychrom{$chrom}}, "$start\t$end");
			}
		push(@region_Chroms, $chrom);
		}
	$i = $i+1;
	}



### then combine the normized data for bootstrap test

open(CONFIG, "<$configfile") or die("No such file or directory: $configfile\n");
my $tmp_readdata;
(undef,$tmp_readdata) = tempfile("readdata_all"."_XXXXXXX",SUFFIX=>".bin",DIR => $tmpdir,OPEN=>0);
my $flag = 0;
my $i=0;
while(<CONFIG>){
        chomp;
        my @row = split(/\t/);
        my $num_items = $#row+1;
	if($num_items!=3 && $num_items!=2){die("Config file must be a 2 or 3 column file\n");}
        if($i>0){
                my $chrom = $row[0];
                my $tumor = $row[1];
		my $normal = "";
		my $cut_range = "3-4";
		if($num_items==3){
	                $normal = $row[2];
			$cut_range = "3-6";
			}

                my $cmd1;
                if($flag==0) {
                	$cmd1 = "$combinefile $tumor $normal | cut -f$cut_range > $tmp_readdata";
                        $flag = 1;
                        }else{
                        $cmd1 = "$combinefile $tumor $normal | cut -f$cut_range >> $tmp_readdata";
                        }
		print "$cmd1\n";
                if(system($cmd1)!=0){die("\n");}
                if(!(-e $tmp_readdata)){
                	die("Error: failed to create the temp file $tmp_readdata\n");
                        }

                }
        $i = $i+1;
        }

close(CONFIG);

#### now perform CNV genotyping
#### firstly count the number of reads in the bins overlapping with the given regions
my $tmp_region_count_file;
(undef,$tmp_region_count_file) = tempfile("region_count_file"."_XXXXXXX",SUFFIX=>".bin",DIR => $tmpdir,OPEN=>0);
my $flag = 0;
foreach my $chrom (keys %regionsbychrom){
	if(exists $regionsbychrom{$chrom}){
		my $tmp_region_file;
		(undef,$tmp_region_file) = tempfile("region_file_$chrom"."_XXXXXXX",SUFFIX=>".bin",DIR => $tmpdir,OPEN=>0);
		open(REGION_CHROM, ">$tmp_region_file") or die("Failed to open $tmp_region_file\n");
		foreach my $interval (@{$regionsbychrom{$chrom}}){
			print REGION_CHROM "$interval\n";
			}
		close(REGION_CHROM);

		my $cmd = "";
		if($flag==0){
			$cmd = "$countRead $normalized_data{$chrom} --chrom $chrom -o $tmp_region_count_file $tmp_region_file\n";
			$flag = 1;
			}else{
			$cmd = "$countRead $normalized_data{$chrom} --chrom $chrom $tmp_region_file >> $tmp_region_count_file\n";
			}
		print "$cmd\n";
		if(system($cmd)!=0) {die("\n");}
		unlink($tmp_region_file) or print "Failed to delete the file $tmp_region_file\n";
		}
	}


my $tmp_region_count_file_nochrom;
(undef,$tmp_region_count_file_nochrom) = tempfile("region_count_file_nochrom"."_XXXXXXX",SUFFIX=>".bin",DIR => $tmpdir,OPEN=>0);
open(OUT_REGION_NOCHROM, ">$tmp_region_count_file_nochrom") or die("Failed to open $tmp_region_count_file_nochrom\n");
open(IN_REGION,  "<$tmp_region_count_file") or die("Failed to open  $tmp_region_count_file\n");
my @region_chroms_by_order = ();
while (<IN_REGION>){
	chomp;
	my @row = split(/\t/);
	push (@region_chroms_by_order, $row[0]);
	my $j=0;
	foreach my $item(@row){
		if($j>1) {print OUT_REGION_NOCHROM "\t";}
		if($j>0) {print OUT_REGION_NOCHROM "$item";}
		$j = $j+1;
		}
	print OUT_REGION_NOCHROM "\n";
	#print "$row[0]\n";
	}
close(OUT_REGION_NOCHROM);
close(IN_REGION);
#print "OK\n";

### now assign pvalues with bootstrap for the given regions.
my $tmp_output;
(undef,$tmp_output) = tempfile("tmp_output"."_XXXXXXX",SUFFIX=>".bin",DIR => $tmpdir,OPEN=>0);
open(OUT_TMP, ">$tmp_output") or die("Failed to open $tmp_output\n");
close(OUT_TMP);
my $cmd = "$genotype $tmp_region_count_file_nochrom $tmp_readdata > $tmp_output\n";
print "$cmd\n";
if(system($cmd)!=0){die("\n");}


### finally report the result
open(OUTPUT, ">$output") or die("Failed to open $output\n");
open(OUT_TMP, "<$tmp_output") or die("Failed to open $tmp_output\n");
my $i = 0;
while(<OUT_TMP>){
	chomp;
	if($i==0){
		print OUTPUT "chrom\t$_\n";
		}else{
		print OUTPUT $region_chroms_by_order[$i-1]."\t$_\n";
		}
	$i = $i+1;
	}
close(OUTPUT);
close(OUT_TMP);

unlink $tmp_output or print "Failed to delete the file $tmp_output\n";
unlink($tmp_region_count_file_nochrom) or print "Failed to delete the file $tmp_region_count_file_nochrom\n";
unlink($tmp_region_count_file) or print "Failed to delete the file $tmp_region_count_file\n";
if(-e $tmp_readdata){
        unlink($tmp_readdata) or print "Failed to delete the file $tmp_readdata\n";
        }

