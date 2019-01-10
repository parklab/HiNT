#!/usr/bin/perl

# Author:
# Dhawal Jain, Dept of Biomedical Informatics, Harvard
#
# This script identifies non-Hi-C chimeric reads and writes them down into 2 output files
# Input: 
#  Pairsam file
# Enzyme:
#  Mbo1/Dpn2/HindIII
# All output options: 
#  FASTQ file: Softclipped part of the non Hi-C chmeric reads 
#  {OUTPUT_PREFIX}_2.txt  File contains non Hi-C chimeric reads (softclipped).\n";
#  {OUTPUT_PREFIX}_3.txt  File contains non Hi-C chimeric reads (reads with MAPQ<30, whole read sequence is reported).\n";
#  {OUTPUT_PREFIX}_2dedup.txt  File contains DUPLICATE REMOVED non Hi-C chimeric reads (softclipped).\n";
#  {OUTPUT_PREFIX}_3dedup.txt  File contains DUPLICATE REMOVED non Hi-C chimeric reads (reads with MAPQ<30, whole read sequence is reported).\n";
#  Note: Not all supplementary alinments are produced by BWA 
#  Remap Fastq file for further analysing the repeat or SV breakpoints
#

use warnings FATAL => "all";
use strict;
use POSIX;
use Getopt::Long;
use Data::Dumper;
use constant { true => 1, false => 0 };
use experimental 'smartmatch';
use utf8;
use open qw(:std :utf8);


######----------------------------------------------------------------------------------------------------------------
###    inputs
######----------------------------------------------------------------------------------------------------------------
my $psam = "";
my $outprefix = "project";
my $re = "MboI";
my $sites = "";
my $cutoff = 2;
my $min_clip = 28;
my $help = 0;

Getopt::Long::GetOptions(
  'psam=s'          => \$psam,
  're=s'            => \$re,
  'sites=s'         => \$sites,
  'outprefix:s'     => \$outprefix,
  'margin:s'        => \$cutoff,
  'min_clip:s'      => \$min_clip,
  'help'            => \$help,
  'h'               => \$help,
) or die "Incorrect input! Use -h for usage.\n";

if ($help) {
  print "\nUsage: perl extract_nonHiC.pl -psam psam.gz -re [STRING] -sites sites.txt -outprefix [STRING] -margin [INT,2] -min_clip [INT,28]\n\n\n";
  print "required options:\n";
  print "  -psam            The Input psam file\n";
  print "  -sites           RE sites along the genome (produced using 'generate_site_positions.py' script from Juicer package)\n";
  print "  -re              Name of the RE [MboI, DpnII, HindIII]\n";
  print "optional options:\n";
  print "  -outprefix       Output file PREFIX [default: project]\n";
  print "  -margin          Upstream Margin for finding genomic RE site [default: 2]\n";
  print "  -min_clip        Minimum length of softclipped part on the read [default : 28]\n";
  print "  -help|-h         Display usage information.\n\n\n";
  print "Default outputs:\n";
  print "  {OUTPUT_PREFIX}_2.txt  File contains non Hi-C chimeric reads (softclipped).\n";
  print "  {OUTPUT_PREFIX}_3.txt  File contains non Hi-C chimeric reads (reads with MAPQ<30, whole read sequence is reported).\n";
  print "  {OUTPUT_PREFIX}_2dedup.txt  File contains DUPLICATE REMOVED non Hi-C chimeric reads (softclipped).\n";
  print "  {OUTPUT_PREFIX}_3dedup.txt  File contains DUPLICATE REMOVED non Hi-C chimeric reads (reads with MAPQ<30, whole read sequence is reported).\n\n\n";
  print "Other hardcoded parameters:\n";
  print "  Enzyme collection: Mbo1/Dpn2/HindIII\n";
  print "  Proximity filter for mapping two mates of the pair : 500bp \n\n\n";
  exit 0;
}

######### Global varialbes
my %counts;
my $motif;
my $len_RE = 8;

if($re eq "DpnII" or $re eq "MboI"){
  $motif = 'GATC';
}elsif($re eq "HindIII"){
  $motif = "AAGCTT";
}else{
  print " Restriction enzyme not specified. Exiting!\n";
  exit 1;
}
$len_RE = 8 if($motif eq "GATC");
$len_RE = 10 if($motif eq "AAGCTT");
#$sites = "/n/data1/hms/dbmi/park/Dhawal/Genomes/hg19_annotations/hg19_DpnII.txt";

######----------------------------------------------------------------------------------------------------------------
###    checks and primary readings
######----------------------------------------------------------------------------------------------------------------
if($psam eq ""){
  print "ERROR: Required input psam file not provided\n";
  exit 1;   
} 
if($sites eq ""){
  print "ERROR: Required genomic TE location file not provided\n";
  exit 1;   
} 
if($motif eq ""){
  print "ERROR: Required value for RE motif not provided\n";
  exit 1;   
} 


### Read Juicer formatted RE motif locations as a hash
open FILE, $sites or die $!;
my %chromosomes;
while (<FILE>) {
  my @locs = split(" ",$_);
  my $key = shift(@locs);
  my $ref = \@locs;
  $chromosomes{$key} = $ref;
}
close(FILE);


######----------------------------------------------------------------------------------------------------------------
###    main
######----------------------------------------------------------------------------------------------------------------
if($psam =~ /.gz/){
  open IN,"zcat $psam|" or next "Can't open file $psam";
}else{
  open IN,"$psam" or next "Can't open file $psam";
}
#open( OUT, ">", $outprefix."_1.fastq") or die "Error: Coudln't open file for writing: $!";
#open( OUT1, ">", $outprefix."_2.txt") or die "Error: Coudln't open file for writing: $!";
open( OUT2, ">", $outprefix."_3.txt") or die "Error: Coudln't open file for writing: $!";


while(<IN> ) {
  next if(/^(\#)/); 
  s/\n//; 
  s/\r//;  
  my @psam = split(/\t+/);
  $counts{flag}{$psam[7]}++;  #flag count
  
  $counts{pairs}{total}++;
  if($counts{pairs}{total} % 100000 == 0 ){
    print "Read-pairs: $counts{pairs}{total} \t";
    print "non-HiC chimeric mates: $counts{mates}{nonHic_chimera} \n";
  }
  ## omit UU,UR,RU reads 
  next if($psam[7] eq "UU" or $psam[7] eq "UR" or $psam[7] eq "RU");
  next if($psam[7] eq "UU" or $psam[7] eq "UR" or $psam[7] eq "RU");

  my($a,$b,$c,$d)=parse_psam_reads($psam[8], $psam[9]);
  my @read1 = @{$a};
  my @read2 = @{$b};
  my @read3 = @{$c};
  my @read4 = @{$d};

  $counts{pairs}{used_pairs}++; 
 
  next if(!@read1 && !@read2);
  ## Get rid of all reads that have atleast 2 N bases
  if(@read1){
    my $num_N1 = $read1[9] =~ tr/N//;
    undef(@read1) if($num_N1>1);
  }
  if(@read2){
    my $num_N2 = $read2[9] =~ tr/N//;
    undef(@read2) if($num_N2>1);
  }
 
  ############################# Read1 present
  if(@read1){
    $counts{mates}{used_mates}++; 
    my $cnt = nonHiC_by_LigationSignature(\*OUT, \*OUT1, \*OUT2, \@read1, \@read3, 1);
    if($cnt==2){
      print "MAJOR ERROR\n";
      exit 1; 
    }
  }
  ############################# Read2 present
  if(@read2){ 
    $counts{mates}{used_mates}++; 
    my $cnt = nonHiC_by_LigationSignature(\*OUT, \*OUT1, \*OUT2, \@read2, \@read4, 2); 
    if($cnt==2){
      print "MAJOR ERROR\n";
      exit 1; 
    }
  }
  undef(@read1);
  undef(@read2);
  undef(@read3);
  undef(@read4);
  undef(@psam);
}

close(IN);
close(OUT);
close(OUT1);
close(OUT2);


### Remove duplicates
my $z;
#$z = remove_duplicates($outprefix."_2.txt",$outprefix."_2dedup.txt");
$z = remove_duplicates($outprefix."_3.txt",$outprefix."_3dedup.txt");


my $hic = $counts{mates}{used_mates} - $counts{mates}{nonHic_chimera};

print "###---------------------------Summary------------------------------------------------------------\n";
print "Input psam file:                                          $psam\n";
print "RE enzyme:                                                $re\n";
print "RE motif:                                                 $motif\n";
print "Input file with genomic location of RE motifs:            $sites\n";
print "Total read pairs in the file:                             $counts{pairs}{total} \n";
print "Flag Summary\n";
foreach my $cl (sort keys %{$counts{flag}}){
  print "  $cl:  $counts{flag}{$cl}\n";
}
print "Total non-HiC read pairs used:                            $counts{pairs}{used_pairs}\n";
print "Total non-HiC mates used:                                 $counts{mates}{used_mates}\n";
print "Total Hi-C chimeric/unmapped etc. mates:                  $hic\n";
print "Total unmapped mates                                      $counts{mates}{count_unmapped}\n";
print "Total unmapped mates with RE motif                        $counts{mates}{motif_in_unmapped}\n";
print "Mates with RE motif at clip junction:                     $counts{mates}{motif_spans}\n";
print "Mates with RE motif within unmapped clip sequence:        $counts{mates}{motif_hardClip}\n";
print "Mates with reference RE site annotation:                  $counts{mates}{motif_ref}\n";
print "Mates with soft-clipping match at proximity <500bp:       $counts{mates}{motif_proximity}\n";
print "Mates with both 5' and 3' clips (>$min_clip bp):          $counts{mates}{double_clip}\n";
print "Filtered mates with both 5' and 3' clips (>$min_clip bp): $counts{mates}{motif_multiple}\n";
print "Total non-Hi-C chimeric mates:                            $counts{mates}{nonHic_chimera}\n";
print "Total non-Hi-C chimeric mates after duplicate removal:    $z \n";


######----------------------------------------------------------------------------------------------------------------
###    subroutines
######----------------------------------------------------------------------------------------------------------------
sub remove_duplicates {
  my ($file_in, $file_out) = @_;
  open( IN1, "<", $file_in) or die "Can't open file $file_in";
  open( OUT2, ">", $file_out) or die "Error: Coudln't open file for writing: $!";

  my %hash;
  while(<IN1>){
    next if(/^(\#)/); 
    my @psam = split(/\t+/);
    $psam[15] =~ s/\n//;

  if($hash{$psam[1]." ".$psam[2]." ".$psam[3]." ".$psam[5]." ".$psam[8]." ".$psam[9]." ".$psam[10]." ".$psam[12]}){
      my @temp = split(" ", $hash{$psam[1]." ".$psam[2]." ".$psam[3]." ".$psam[5]." ".$psam[8]." ".$psam[9]." ".$psam[10]." ".$psam[12]});
      if($temp[4] < $psam[4]){
          $hash{$psam[1]." ".$psam[2]." ".$psam[3]." ".$psam[5]." ".$psam[8]." ".$psam[9]." ".$psam[10]." ".$psam[12]} = $psam[0]." ".$psam[1]." ".$psam[2]." ".$psam[3]." ".$psam[4]." ".$psam[5]." ".$psam[6]." ".$psam[7]." ".$psam[8]." ".$psam[9]." ".$psam[10]." ".$psam[11]." ".$psam[12]." ".$psam[13]." ".$psam[14]." ".$psam[15];
      }else{
        next;
      }
    }else{
      $hash{$psam[1]." ".$psam[2]." ".$psam[3]." ".$psam[5]." ".$psam[8]." ".$psam[9]." ".$psam[10]." ".$psam[12]} = $psam[0]." ".$psam[1]." ".$psam[2]." ".$psam[3]." ".$psam[4]." ".$psam[5]." ".$psam[6]." ".$psam[7]." ".$psam[8]." ".$psam[9]." ".$psam[10]." ".$psam[11]." ".$psam[12]." ".$psam[13]." ".$psam[14]." ".$psam[15];
    }
  }

  while(my($i, $j) = each %hash){
    my @temp = split(" ",$j);
    print OUT2 join("\t",@temp),"\n";
  }
  close(IN1);
  close(OUT2);
  return(scalar(keys(%hash)));
}

sub nonHiC_by_LigationSignature {
  my ($fh1, $fh2, $fh3, $a, $b, $readnum ) = @_;
  my @read1 = @{$a};
  my @read3 = @{$b};
  
  ## uses following global variables
  # $len_RE, $min_clip, %chromosomes, $cutoff, $motif, %counts

  ## Three pronged approach:
  # 1. Check if RE ligation motif spans clip junction
  # 2. If hardclipped part of the read (i.e. unmapped) bears RE ligation motif (this represents likely indels in the analyzed genome)
  # 3. If known reference-genome RE annotation overlaps the clip location for either primary and supplement alignments (with/without RE motif)


  ################## 
  ## Nonmapping reads do not have CIGAR, check if they have RE motif. If not, print them and move on
  ################## 
  if($read1[5] eq "*"){
      $counts{mates}{count_unmapped}++;
      if(re_search($read1[9],$motif,-1)==0){
          #print $fh1 "@","$read1[0]_$readnum 1:N:0:0 $read1[-2]\n$read1[9]\n+\n$read1[10]\n";
          #print $fh2 "$read1[0]_$readnum\t$read1[2]\t$read1[3]\t$read1[-1]\t$read1[4]\t$read1[5]\t$read1[9]\t$read1[10]\t";
          print $fh3 "$read1[0]_$readnum\t$read1[2]\t$read1[3]\t$read1[-1]\t$read1[4]\t$read1[5]\t$read1[9]\t$read1[10]\t";
          if(@read3){
             #print $fh2 "$read3[2]\t$read3[3]\t$read3[-1]\t$read3[4]\t$read3[5]\t$read3[9]\t-\t-\n";
             print $fh3 "$read3[2]\t$read3[3]\t$read3[-1]\t$read3[4]\t$read3[5]\t$read3[9]\t-\t-\n";
          }else{
             #print $fh2 "-\t-\t-\t-\t-\t-\t-\t-\n";
             print $fh3 "-\t-\t-\t-\t-\t-\t-\t-\n";
          }
      }else{
          $counts{mates}{motif_in_unmapped}++; 
      }
      return(0); #success
  }

  ################## 
  ## CIGAR splitting 
  ################## 
  my @len1 = split (/\D+/,$read1[5]); # storing the length per operation
  my @ops1 = split (/\d+/,$read1[5]); # storing the operation
  shift @ops1; # remove the empty first element
  return(1) unless (scalar @len1 == scalar @ops1);
  return(1) if(!@ops1);

  ################## 
  ## Non-HiC chimera test
  ################## 
  my $motif_len = length($motif);
  
  my ($seq1, $qual1) = "";
  my $olap = "-";
  my $xalign = "-";
  
  my %localcount;
  $localcount{motif_spans}=0;
  $localcount{motif_hardClip}=0;
  $localcount{motif_ref}=0;
  $localcount{count_prox}=0;
  $localcount{nonHic_chimera}=0;

  ## 5' clip of read1
  my $is_match_a=0;
  if($ops1[0] eq "S" && $len1[0] >=  $min_clip){    
    if(re_search($read1[9],$motif, $len1[0])>0){
       $is_match_a++; 
       $localcount{motif_spans}++;
    }elsif(RE_within_HardclippedSeq(\@read1,\@len1,\@ops1,\@read3)>0){
       $is_match_a++; 
       $localcount{motif_hardClip}++; 
    }elsif(@read3){
      my $x1= $read1[3];
      my $x2= $read3[3];
      if($read3[-1] eq "+"){ $x2 = $read3[3]+length($read3[9]); }
      my $pos1 = &binarysearch( $x1, $chromosomes{$read1[2]});
      my $pos2 = &binarysearch( $x2, $chromosomes{$read3[2]});      
      if( (($pos1-$cutoff)<= $x1 and ($pos1+$motif_len+$cutoff)>=$x1) or (($pos2-$cutoff)<=$x2 and ($pos2+$motif_len+$cutoff)>=$x2 ) ){
        $is_match_a++; 
        $localcount{motif_ref}++; 
      }
    }

    if($is_match_a==0){
      if(linearProximity_bw_mapped_softclip(\@read1,\@read3)==0){
            $localcount{nonHic_chimera}++;
            $seq1 = substr($read1[9], 0, $len1[0]);
            $qual1 = substr($read1[10], 0, $len1[0]);
        
            #print $fh1 "@","$read1[0]_$readnum 1:N:0:0 $read1[-2]\n$seq1\n+\n$qual1\n";
            #print $fh2 "$read1[0]_$readnum\t$read1[2]\t$read1[3]\t$read1[-1]\t$read1[4]\t$read1[5]\t$seq1\t$qual1\t";
                if($read1[4]<30){
                    print $fh3 "$read1[0]_$readnum\t$read1[2]\t$read1[3]\t$read1[-1]\t$read1[4]\t$read1[5]\t$read1[9]\t$read1[10]\t";
                }else{
                    print $fh3 "$read1[0]_$readnum\t$read1[2]\t$read1[3]\t$read1[-1]\t$read1[4]\t$read1[5]\t$seq1\t$qual1\t";
                }
                if(@read3){
                  if($read3[-3] =~ /^XA:Z:/){
                      $xalign = $read3[-3];
                  }
                  $olap = find_olap_in_chimera_mapping(\@read1,\@len1,\@ops1,\@read3);
                  if(!$olap){$olap="-";}
                  #print $fh2 "$read3[2]\t$read3[3]\t$read3[-1]\t$read3[4]\t$read3[5]\t$read3[9]\t$olap\t$xalign\n";
                  print $fh3 "$read3[2]\t$read3[3]\t$read3[-1]\t$read3[4]\t$read3[5]\t$read3[9]\t$olap\t$xalign\n";
                }else{
                  #print $fh2 "-\t-\t-\t-\t-\t-\t-\t-\n";
                  print $fh3 "-\t-\t-\t-\t-\t-\t-\t-\n";
                }
      }else{
        $localcount{count_prox}++;
      }
    }
  }


  ## 3' clip of read1
  my $is_match_b=0;
  if($ops1[-1] eq "S" && $len1[-1] >=  $min_clip){
    if(re_search($read1[9],$motif, (length($read1[9])-$len1[-1]) )>0){
       $is_match_b++; 
       $localcount{motif_spans}++;
    }elsif(RE_within_HardclippedSeq(\@read1,\@len1,\@ops1,\@read3)>0){
       $is_match_b++; 
       $localcount{motif_hardClip}++; 
    }elsif(@read3){
      my $x1= $read1[3]+length($read1[9])-$len1[-1];
      my $x2= $read3[3];
      if($read3[-1] eq "-"){ $x2 = $read3[3]+length($read3[9]); }
      my $pos1 = &binarysearch( $x1, $chromosomes{$read1[2]});
      my $pos2 = &binarysearch( $x2, $chromosomes{$read3[2]});      
      if( (($pos1-$cutoff)<= $x1 and ($pos1+$motif_len+$cutoff)>=$x1) or (($pos2-$cutoff)<=$x2 and ($pos2+$motif_len+$cutoff)>=$x2 )){
        $is_match_b++; 
        $localcount{motif_ref}++; 
      }
    }

    if($is_match_b==0){
      if(linearProximity_bw_mapped_softclip(\@read1,\@read3)==0){
            $localcount{nonHic_chimera}++;

            $seq1 = substr($read1[9],-$len1[-1],$len1[-1]);
            $qual1 = substr($read1[10],-$len1[-1],$len1[-1]);           
            #print $fh1 "@","$read1[0]_$readnum 1:N:0:0 $read1[-2]\n$seq1\n+\n$qual1\n";
            #print $fh2 "$read1[0]_$readnum\t$read1[2]\t$read1[3]\t$read1[-1]\t$read1[4]\t$read1[5]\t$seq1\t$qual1\t"; 
                if($read1[4]<30){
                    print $fh3 "$read1[0]_$readnum\t$read1[2]\t$read1[3]\t$read1[-1]\t$read1[4]\t$read1[5]\t$read1[9]\t$read1[10]\t";
                }else{
                    print $fh3 "$read1[0]_$readnum\t$read1[2]\t$read1[3]\t$read1[-1]\t$read1[4]\t$read1[5]\t$seq1\t$qual1\t";
                }
                if(@read3){
                  if($read3[-3] =~ /^XA:Z:/){
                      $xalign = $read3[-3];
                  }
                    $olap = find_olap_in_chimera_mapping(\@read1,\@len1,\@ops1,\@read3);
                  if(!$olap){$olap="-";}
                    #print $fh2 "$read3[2]\t$read3[3]\t$read3[-1]\t$read3[4]\t$read3[5]\t$read3[9]\t$olap\t$xalign\n";
                    print $fh3 "$read3[2]\t$read3[3]\t$read3[-1]\t$read3[4]\t$read3[5]\t$read3[9]\t$olap\t$xalign\n";
               }else{
                    #print $fh2 "-\t-\t-\t-\t-\t-\t-\t-\n";
                    print $fh3 "-\t-\t-\t-\t-\t-\t-\t-\n";
                }
      }else{
         $localcount{count_prox}++;
      }      
    }
  }


  if($ops1[-1] eq "S" && $len1[-1] >=  $min_clip && $ops1[0] eq "S" && $len1[0] >=  $min_clip){
    $counts{mates}{double_clip}++;
  }
  
  ## counters
  if($localcount{motif_spans}>0 && $localcount{motif_hardClip}==0 && $localcount{motif_ref}==0 && $localcount{count_prox}==0 && $localcount{nonHic_chimera}==0){
    $counts{mates}{motif_spans}++;
  }elsif($localcount{motif_spans}>0 && $localcount{motif_hardClip}==0 && $localcount{motif_ref}==0 && $localcount{count_prox}==0 && $localcount{nonHic_chimera}==0){
    $counts{mates}{motif_spans}++;
  }elsif($localcount{motif_spans}==0 && $localcount{motif_hardClip}>0 && $localcount{motif_ref}==0 && $localcount{count_prox}==0 && $localcount{nonHic_chimera}==0){
    $counts{mates}{motif_hardClip}++;
  }elsif($localcount{motif_spans}==0 && $localcount{motif_hardClip}==0 && $localcount{motif_ref}>0 && $localcount{count_prox}==0 && $localcount{nonHic_chimera}==0){
    $counts{mates}{motif_ref}++;
  }elsif($localcount{motif_spans}==0 && $localcount{motif_hardClip}==0 && $localcount{motif_ref}==0 && $localcount{count_prox}>0 && $localcount{nonHic_chimera}==0){
    $counts{mates}{motif_proximity}++;
  }elsif($localcount{motif_spans}==0 && $localcount{motif_hardClip}==0 && $localcount{motif_ref}==0 && $localcount{count_prox}==0 && $localcount{nonHic_chimera}>0){
    $counts{mates}{nonHic_chimera}++;
  }elsif($localcount{nonHic_chimera}>0){
    $counts{mates}{nonHic_chimera}++;
  }elsif($localcount{nonHic_chimera}==0 && ($is_match_a>0 and $is_match_b>0)){
    $counts{mates}{motif_multiple}++;
  }

  
  if($localcount{motif_spans}==0 && $localcount{motif_hardClip}==0 && $localcount{motif_ref}==0 && $localcount{count_prox}==0 && $localcount{nonHic_chimera}==0){
    if($is_match_a>0 or $is_match_b>0){
     return(2);
    }else{
      return(0);
    }

  }else{
    return(0);
  }
}

sub linearProximity_bw_mapped_softclip {
  my ($a,$b) = @_;
  my @read1 = @{$a};
  my @read3 = @{$b};
  
  if(!@read3 or !@read1){   ## If suppl alignment is absent, return 0
    return(0);
  } 
  elsif($read3[4]<30){   ## If suppl alignment has alignment score < 30, return 0
    return(0);
  }
  elsif($read1[2] ne $read3[2] ){   ## If suppl alignment maps to another chromosome, return 0
    return(0);
  }
  elsif($read1[2] eq "*"){ ## If read is umapped, return 0
    return(0);
  }
  elsif( abs($read1[3]-$read3[3]) <= 500 && $read1[2] eq $read3[2] ){  
  ## Tetracutter cuts every 256bp on an average, hence, minimum resolution 256bp. Any two fragments within 2times the distance are likely to be non-informative
    return(1);
  }
  else{ ## for any other unforseen cases, return 0
    return(0);
  }
}

sub find_olap_in_chimera_mapping {  #find_olap_in_chimera_mapping(\@read1,\@len1,\@ops1,\@read3)
  my ($a,$b,$c,$d) =@_;
  my @read1 = @{$a};
    my @len1 = @{$b};
    my @ops1 = @{$c};
    my @read3 = @{$d};

  if(!@read3){
      return("");
  }

  my @len3 = split (/\D+/,$read3[5]); # storing the length per operation
  my @ops3 = split (/\d+/,$read3[5]); # storing the operation
  shift @ops3; # remove the empty first element
  unless (scalar @len3 == scalar @ops3){
    return("");
  }
  if(!@ops3){
   return("");
  }

  my $seq = $read1[9];
  if($read1[-1] eq "-") {
    @len1 = reverse(@len1);
    @ops1 = reverse(@ops1);
    $seq = reverse_complement($seq);
  }
  if($read3[-1] eq "-") {
    @len3 = reverse(@len3);
    @ops3 = reverse(@ops3);
  }

    my $j=0;
    my $i=0;
    for($i=0;$i<scalar(@len1);$i++){
      if($ops1[$i] eq "D"){
        next;
      } 
      $j = $j+$len1[$i];
      if($ops1[$i] eq "S" or $ops1[$i] eq "H"){
        substr($seq,$j-$len1[$i],$len1[$i]) = "-" x $len1[$i];
      }
    }
    $j=0;
    $i=0;
    for($i=0;$i<scalar(@len3);$i++){
      if($ops3[$i] eq "D"){
        next;
      } 
      $j = $j+$len3[$i];
      if($ops3[$i] eq "S" or $ops3[$i] eq "H"){
        substr($seq,$j-$len3[$i],$len3[$i]) = "-" x $len3[$i];
      }
    }

   $seq =~ s/^-{1,}//s;
   $seq =~ s/-{1,}$//s;
   return($seq);
}

sub re_search{  #re_search($seq,$motif,-1)
  my ($seq,$motif,$clip) = @_;
  ## if CIGAR is unavailable, use clip as -1
  
  if($motif eq "" or !defined $motif){
    print " Motif not defined in Dist2RE_fromSeq sub!! exiting! \n";
    exit 1;
  }elsif($motif eq "GATC"){  ########## DpnII/ Mbo1
     my $c = index($seq, "GATCGATC",0);
     return(0) if($c == -1);
     return(1) if($c != -1 and $clip == -1);
     $clip -=1 if($clip !=0); 
     if($c <= $clip and ($c+8) >= $clip ){
        return(1);     
     }else{
        return(0);
     }     
  }elsif($motif eq "AAGCTT"){  ######## HindIII
     my $c = index($seq, "AAGCTAGCTT",0);
     return(0) if($c == -1);
     return(1) if($c != -1 and $clip == -1);
     $clip -=1 if($clip !=0); 
     if($c <= $clip and ($c+10) >= $clip ){
        return(1);     
     }else{
        return(0);
     }     
  }
}

sub re_search_stringent {  #re_search_stringent($seq)
  my $seq= shift;
  my @patterns;

  if($motif eq "GATC"){
     @patterns = (qr/GATCGATC/);
  }
  if($motif eq "AAGCTT"){
     @patterns = (qr/AAGCTAGCTT/);
  }

  my $is_match=0;
  $is_match++ if( $seq ~~ @patterns );
  return($is_match);
}

sub RE_within_HardclippedSeq { #RE_within_HardclippedSeq(\@read1,\@len1,\@ops1,\@read3)
  my ($a,$b,$c,$d) =@_;
  my @read1 = @{$a};
  my @len1 = @{$b};
  my @ops1 = @{$c};
  my @read3 = @{$d};

  if(!@read3){
      return(0)
  }
  if($read3[5] eq "*"){
      return(0)
  }

  my @len3 = split (/\D+/,$read3[5]); # storing the length per operation
  my @ops3 = split (/\d+/,$read3[5]); # storing the operation
  shift @ops3; # remove the empty first element
  unless (scalar @len3 == scalar @ops3){
    return(0);
  }
  if(!@ops3){
   return(0);
  }

  my $seq = $read1[9];
  if($read1[-1] eq "-") {
      @len1 = reverse(@len1);
      @ops1 = reverse(@ops1);
      $seq = reverse_complement($seq);
  }
  if($read3[-1] eq "-") {
      @len3 = reverse(@len3);
      @ops3 = reverse(@ops3);
  }



  my $j=0;
  my $i=0;
  for($i=0;$i<scalar(@len1);$i++){
      if($ops1[$i] eq "D"){
        next;
      } 
      $j = $j+$len1[$i];
      if($ops1[$i] eq "M" or $ops1[$i] eq "I"){
        substr($seq,$j-$len1[$i],$len1[$i]) = "-" x $len1[$i];
      }
  }
  $j=0;
  $i=0;
  for($i=0;$i<scalar(@len3);$i++){
      if($ops3[$i] eq "D"){
        next;
      } 
      $j = $j+$len3[$i];
      if($ops3[$i] eq "M" or $ops3[$i] eq "I"){
        substr($seq,$j-$len3[$i],$len3[$i]) = "-" x $len3[$i];
      }
  }

  $seq =~ s/^[ATGC]+/-/s;
  $seq =~ s/[ATGC]+$/-/g;
  $seq =~ s/-//g;
  if(length($seq)>= $len_RE){
     if(re_search_stringent($seq)==0){
        return(0);
      }else{
        return(1);
      }
  }else{
      return(0);
  }
}

# Returns distance between x and it's nearest RE site
sub binarysearch { #binarysearch($pos, \@locations)
  my ($x, $a) = @_;            # search for x in array a
  my ($l, $u) = (0, @$a - 1);          # lower, upper end of search interval
  my $i;                               # index of probe
  my $q1=0;
  my $q2=0;


  while ($l <= $u) {
    $i = int(($l + $u)/2);
    if ($a->[$i] < $x) {
      $l = $i+1;
    }
    elsif ($a->[$i] > $x) {
      $u = $i-1;
    } 
    else {
      # return $i+1; #original 
      $q1 = abs($a->[$i] - $x);
      $q2 = $q1;  
      if($a->[$i+1]){
        $q2 = abs($a->[$i+1] - $x);
      }
      if($q1 > $q2){
        return($a->[$i+1]); # found
      }else{
        return($a->[$i]); # found
      }
      #if($upper==1){
      #  return($q2);
      #}else{
      #  return($q1);
      #}
      
    }
  }
  $q2 = abs($a->[$l-1] - $x);
  $q1 = $q2;
  if($a->[$l]){
    $q1 = abs($a->[$l] - $x);  
  }
  if($q1 < $q2){
    return($a->[$l]); # found
  }else{
    return($a->[$l-1]); # found
  }
  #if($upper==1){
  #  return($q1);
  #}else{
  #  return($q2);
  #}
  #return($l); #original         
}

sub flag2binary { #flag2binary($num)
  my $num = shift;
  my $leftover;
  my $res="";
  while ($num>= 1) {
      $leftover = $num % 2;  
      $res = $leftover." ".$res;
      $num = $num / 2;   
    }

    my @res = split(" ",reverse($res));
    my $i=scalar @res;

    for($i=scalar@res,$i <= 10, $i++){
      push(@res,0);
    }

  return(@res);
}

sub reverse_complement {  #reverse_complement($seq)
        my $seq = shift;

        my $rcseq = reverse($seq);
        $rcseq =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;

        return $rcseq;
}

sub parse_psam_reads { #parse_psam_reads($side1,$side2)
  my ($side1, $side2) =@_;
  $side1=~ s/\n//; 
  $side1=~ s/\r//;  
  $side2=~ s/\n//; 
  $side2=~ s/\r//;  

  my @a = split("NEXT_SAM",$side1);
  my @b = split("NEXT_SAM",$side2);
  push(@a, @b);
  #$cq{scalar(@a)}++;

  my (@r1, @r2, @r3, @r4);
  #0  1 0x1 template having multiple segments in sequencing
  #1  2 0x2 each segment properly aligned according to the aligner
  #2  4 0x4 segment unmapped
  #3  8 0x8 next segment in the template unmapped
  #4  16 0x10 SEQ being reverse complemented
  #5  32 0x20 SEQ of the next segment in the template being reverse complemented
  #6  64 0x40 the first segment in the template
  #7  128 0x80 the last segment in the template
  #8  256 0x100 secondary alignment
  #9  512 0x200 not passing filters, such as platform/vendor quality controls
  #10  1024 0x400 PCR or optical duplicate
  #11 2048 0x800 supplementary alignment

  foreach my $l (@a){
    my @temp;
    $l =~ s/^\x{019}*//;
    $l =~ s/\x{019}*$//;

    @temp = split("\x{019}",$l);  ## check the separator
  
    my @cigar = flag2binary($temp[1]);
    if($cigar[4]==1){
      push(@temp,"-");
     }else{
       push(@temp,"+")
    }
    @r1 = @temp if($cigar[6]==1 && $cigar[8]==0);
    @r2 = @temp if($cigar[7]==1 && $cigar[8]==0);
    
    ## replace supl. alignment for 1st read only if the mapping quality of the alternate alignment is higher 
    if($cigar[6]==1 && $cigar[8]==1){
      if(@r3){
        if($r3[4] < $temp[4]){
          @r3=@temp;
        }
      }else{
        @r3=@temp;
      }
    }

    ## replace supl. alignment for 2nd read only if the mapping quality of the alternate alignment is higher 
    if($cigar[7]==1 && $cigar[8]==1){
      if(@r4){
        if($r4[4] < $temp[4]){
          @r4=@temp;
        }
      }else{
        @r4=@temp;
      }
    }

  }

  return(\@r1,\@r2, \@r3, \@r4);
}

#####----------------------------------------- END --------------------------- 