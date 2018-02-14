#!/usr/bin/perl

# Author:
# Dhawal Jain, Dept of Biomedical Informatics, Harvard

# This script identifies non-Hi-C chimeric reads and writes them down into following output files
# The script uses pairsam file produced using pairsamtools as an input and requires pre-classification of the reads 
# At present, the script is customized for Mbo1/Dpn2 RE

# FASTQ file: Softclipped part of the non Hi-C chmeric reads 
# {OUTPUT_PREFIX}_2.txt  File contains non Hi-C chimeric reads (softclipped).\n";
# {OUTPUT_PREFIX}_3.txt  File contains non Hi-C chimeric reads (reads with MAPQ<30, whole read sequence is reported).\n";
# {OUTPUT_PREFIX}_2dedup.txt  File contains DUPLICATE REMOVED non Hi-C chimeric reads (softclipped).\n";
# {OUTPUT_PREFIX}_3dedup.txt  File contains DUPLICATE REMOVED non Hi-C chimeric reads (reads with MAPQ<30, whole read sequence is reported).\n";
# Note: Not all supplementary alinments are produced by BWA
# Remap Fastq file for further analysing the repeat or SV breakpoints


use warnings FATAL => "all";
use strict;
use POSIX;
use Getopt::Long;
use Data::Dumper;
use List::MoreUtils qw(uniq);
use constant { true => 1, false => 0 };
#use experimental 'smartmatch';
use utf8;
use open qw(:std :utf8);


#=head
#########################################################################################################################
######## inputs, modified for parsing pairsam output
#########################################################################################################################
my $psamfile = "";
my $output_file = "stdout";
my $help = 0;
my $cutoff = 12;
my $leg = 12;
my $min_olap = 28;
my $sites = "/n/data1/hms/dbmi/park/Dhawal/Genomes/hg19_annotations/hg19_DpnII.txt";

Getopt::Long::GetOptions(
  'psam=s'          => \$psamfile,
  'sites=s'         => \$sites,
  'output_file:s'            => \$output_file,
  'margin_up:s'        => \$cutoff,
  'margin_down:s'      => \$leg,
  'min_olap:s'         => \$min_olap,
  'help'             => \$help,
  'h'                => \$help,
  ) or die "Incorrect input! Use -h for usage.\n";

if ($help) {
  print "\nUsage: perl extract_nonHiC.pl -psam psamfile.gz -sites sites.txt -output_file [STRING] -margin_up [INT,12] -margin_down [INT,12] -min_olap [INT,28]\n";
  print "Options:\n";
  print "  -psam             The Input psam file\n";
  print "  -sites     RE sites along the genome\n";
  print "  -output_file     Output file PREFIX\n";
  print "  -margin_up     Upstream Margin for finding RE site\n";
  print "  -margin_down     Upstream Margin for finding RE site\n";
  print "  -min_olap     Minimum length of softclipped part on the read\n";
  print "  -help|-h         Display usage information.\n\n\n";
  print "   Default outputs:\n";
  print "  {OUTPUT_PREFIX}_2.txt  File contains non Hi-C chimeric reads (softclipped).\n";
  print "  {OUTPUT_PREFIX}_3.txt  File contains non Hi-C chimeric reads (reads with MAPQ<30, whole read sequence is reported).\n";
  print "  {OUTPUT_PREFIX}_2dedup.txt  File contains DUPLICATE REMOVED non Hi-C chimeric reads (softclipped).\n";
  print "  {OUTPUT_PREFIX}_3dedup.txt  File contains DUPLICATE REMOVED non Hi-C chimeric reads (reads with MAPQ<30, whole read sequence is reported).\n";
  print "  The script only works for Mbo1/Dpn2 REs at the momemnt. Dhawal needs to include additional RE sites in fuzzy search\n";
  print "  You need to input parsed pairsam file\n";
  print "  Presently, one parameter, proximity of mapping, is hardcoded. It represents proximity between primary and secondary alignments for a mate\n\n";
  exit 0;
}

## Since the code does not have >= comparison
$cutoff = $cutoff+1;
$leg = $leg+1;


#########################################################################################################################
######## RE site maps
#########################################################################################################################
# read in restriction site file and store as multidimensional array
# array are referenced for quik searching
# RE site locations are produced using generate_site_positions.py script from Juicer
#my $site_file = "A:/work/hg19_annotations/hg19_DpnII.txt";
open FILE, $sites or die $!;

my %chromosomes;
while (<FILE>) {
  my @locs = split(" ",$_);
  my $key = shift(@locs);
  my $ref = \@locs;
  $chromosomes{$key} = $ref;
}
close(FILE);



#########################################################################################################################
######## Extracting chimeric reads
#########################################################################################################################
#open( OUT, ">", $output_file."_1.fastq") or die "Error: Coudln't open file for writing: $!";
open( OUT1, ">", $output_file."_2.txt") or die "Error: Coudln't open file for writing: $!";
open( OUT2, ">", $output_file."_3.txt") or die "Error: Coudln't open file for writing: $!";
open IN,"zcat $psamfile|" or die "Can't open file $psamfile";


my $RE_site = "GATC";  ## RE site
my $len_RE = length($RE_site); ## for 4bp long RE site, search space is 24bp centered on the softclipping junction, of which 8bp need to match
my $total_read_pairs =0; ## Total number of reads
my $total = 0;  ## Total number of reported alignemnts (include supplimentary alignments)
my $nonHic_chimera = 0; 
my @counts=0;

while(<IN> ) {
  next if(/^(\#)/); 
  s/\n//; 
  s/\r//;  
  my @psam = split(/\t+/);
  
  $total++;
  if($total % 100000 == 0 ){
    print "Read-pairs: $total \t";
    print "non-HiC chimeric reads: $counts[4] \n";
  }

  ## Neglect LL and CX reads 
  next if($psam[7] eq "LL" or $psam[7] eq "CX");

  my($a,$b,$c,$d)=parse_psam_reads($psam[8], $psam[9]);
  my @read1 = @{$a};
  my @read2 = @{$b};
  my @read3 = @{$c};
  my @read4 = @{$d};


  $total_read_pairs++;  ## counts total reads. Each counter represents paired read with primary alignment 

  ## QUICK CHECKS!! 
  ## Get rid of all reads that have atleast 2 consecutive N bases
  if(!@read1 && !@read2){
    next;
  }
  if($read1[9]=~ m/NN/){
    undef(@read1);
  }
  if($read2[9]=~ m/NN/){
    undef(@read2);
  }
 
  ############################# Both primary alignments are reported
  if(@read1 && @read2){
    @counts = nonHiC_by_LigationSignature(\*OUT, \*OUT1, \*OUT2, \@read1, \@read3, $counts[0], $counts[1],$counts[2],$counts[3],$counts[4],$counts[5], 1);
    @counts = nonHiC_by_LigationSignature(\*OUT, \*OUT1, \*OUT2, \@read2, \@read4, $counts[0], $counts[1],$counts[2],$counts[3],$counts[4],$counts[5], 2);
  }

  ############################# Read2 present but not Read1
  if(!@read1 && @read2){ 
     @counts = nonHiC_by_LigationSignature(\*OUT, \*OUT1, \*OUT2, \@read2, \@read4, $counts[0], $counts[1],$counts[2],$counts[3],$counts[4],$counts[5], 2); 
  }

  ############################# Read1 present but not Read2
  if(!@read2 && @read1){   
     @counts = nonHiC_by_LigationSignature(\*OUT, \*OUT1, \*OUT2, \@read1, \@read3, $counts[0], $counts[1],$counts[2],$counts[3],$counts[4],$counts[5], 1);
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

### Remove duplicates
my $z = remove_duplicates($output_file."_2.txt",$output_file."_2dedup.txt");
$z = remove_duplicates($output_file."_3.txt",$output_file."_3dedup.txt");

print "*********************************************\n***Summary of the psam file ",$psamfile,"\n*********************************************\n";
print "Total Alignments processed: $total \n";
print "Total read pairs: $total_read_pairs\n";
my $hic = $total_read_pairs*2 - $counts[4];
print "Total non-Hi-C chimeric mates: $counts[4]\n";
print "Total Hi-C chimeric/unmapped mates: $hic\n";
print "Mates filtered after RE filtering: $counts[0]\n";
print "Mates filtered after fuzzy RE filtering in unmapped region: $counts[1]\n";
print "Mates filtered after fuzzy RE filtering in 12bp window: $counts[2]\n";
print "Mates filtered after Soft-clipping match proximity filtering: $counts[3]\n";
print "RE site: $RE_site\n";
print "Total non-Hi-C chimeric mates after duplicate removal: $z \n";


##################################################################################################################################
########### Subroutines
########################################################################################### #######################################
sub remove_duplicates {
  my ($file_in, $file_out) = @_;
  open( IN1, "<", $file_in) or die "Can't open file $psamfile";
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


sub parse_psam_reads {
  my ($side1, $side2) =@_;
  $side1=~ s/\n//; 
  $side1=~ s/\r//;  
  $side2=~ s/\n//; 
  $side2=~ s/\r//;  

  my @a = split("NEXT_SAM",$side1);
  my @b = split("NEXT_SAM",$side2);
  push(@a, @b);

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

sub nonHiC_by_LigationSignature {
  my ($fh1, $fh2, $fh3, $a, $b, $count_re,$count_restr,$count_restrSC,$count_prox, $nonHic_chimera,$count_unmapped,$readnum ) = @_;
  my @read1 = @{$a};
  my @read3 = @{$b};

  ################## 
  ## Nonmapping reads do not have CIGAR, check if they have RE. If not, print them and move on
  ################## 
  if($read1[5] eq "*"){
      if(re_fuzzy_search($read1[9])==0){
             #print $fh1 "@","$read1[0]_$readnum 1:N:0:0 $read1[-2]\n$read1[9]\n+\n$read1[10]\n";
             print $fh2 "$read1[0]_$readnum\t$read1[2]\t$read1[3]\t$read1[-1]\t$read1[4]\t$read1[5]\t$read1[9]\t$read1[10]\t";
             print $fh3 "$read1[0]_$readnum\t$read1[2]\t$read1[3]\t$read1[-1]\t$read1[4]\t$read1[5]\t$read1[9]\t$read1[10]\t";
             if(@read3){
               print $fh2 "$read3[2]\t$read3[3]\t$read3[-1]\t$read3[4]\t$read3[5]\t$read3[9]\t-\t-\n";
              print $fh3 "$read3[2]\t$read3[3]\t$read3[-1]\t$read3[4]\t$read3[5]\t$read3[9]\t-\t-\n";
             }else{
               print $fh2 "-\t-\t-\t-\t-\t-\t-\t-\n";
               print $fh3 "-\t-\t-\t-\t-\t-\t-\t-\n";
             }
      }else{
        $count_re++; 
      }
      $count_unmapped++;

      return($count_re,$count_restr,$count_restrSC,$count_prox, $nonHic_chimera,$count_unmapped);
  }

  ################## 
  ## CIGAR splitting 
  ################## 
  my @len1 = split (/\D+/,$read1[5]); # storing the length per operation
  my @ops1 = split (/\d+/,$read1[5]); # storing the operation
  shift @ops1; # remove the empty first element
  unless (scalar @len1 == scalar @ops1){
    return($count_re,$count_restr,$count_prox, $nonHic_chimera,$count_unmapped);
  }
  if(!@ops1){
    return($count_re,$count_restr,$count_restrSC,$count_prox, $nonHic_chimera,$count_unmapped);
  }

  ################## 
  ## Non-HiC chimera test
  ################## 
  my ($seq1, $qual1) = "";
  my $olap = "-";
  my $xalign = "-";

  ## check if the read has minimum 28bp softclipping at 5'/3' end 
  ## if yes, check the RE motif around the softclipping junction, usng 3 pronged approach
  my @is_re_match = re_comprehensive_search_V1(\@read1,\@len1,\@ops1,\@read3,\%chromosomes,$len_RE,$min_olap,$count_re,$count_restr,$count_restrSC);
  $count_re = $is_re_match[0];
  $count_restr = $is_re_match[1];
  $count_restrSC = $is_re_match[2];

  if($is_re_match[3]==0){
   ## No RE present  
   #$count_restr++;

  ################### 5' end of the first read
    if($ops1[0] eq "S" && $len1[0] >=  $min_olap){
      if(linearProximity_bw_mapped_softclip(\@read1,\@read3)==0){
            $nonHic_chimera++;

            $seq1 = substr($read1[9], 0, $len1[0]);
            $qual1 = substr($read1[10], 0, $len1[0]);
            #print $fh1 "@","$read1[0]_$readnum 1:N:0:0 $read1[-2]\n$seq1\n+\n$qual1\n";
        
            ## Output to text file with information about the main and softclipped read 
            print $fh2 "$read1[0]_$readnum\t$read1[2]\t$read1[3]\t$read1[-1]\t$read1[4]\t$read1[5]\t$seq1\t$qual1\t";
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
                  print $fh2 "$read3[2]\t$read3[3]\t$read3[-1]\t$read3[4]\t$read3[5]\t$read3[9]\t$olap\t$xalign\n";
                  print $fh3 "$read3[2]\t$read3[3]\t$read3[-1]\t$read3[4]\t$read3[5]\t$read3[9]\t$olap\t$xalign\n";
                }else{
                  print $fh2 "-\t-\t-\t-\t-\t-\t-\t-\n";
                  print $fh3 "-\t-\t-\t-\t-\t-\t-\t-\n";
                }
      }
      else{
        $count_prox++;
      }
    }

  ################### 3' end of the first read
    if($ops1[-1] eq "S" and $len1[-1] >=  $min_olap){
      if(linearProximity_bw_mapped_softclip(\@read1,\@read3)==0){
            $nonHic_chimera++;

            $seq1 = substr($read1[9],-$len1[-1],$len1[-1]);
            $qual1 = substr($read1[10],-$len1[-1],$len1[-1]);
            #print $fh1 "@","$read1[0]_$readnum 1:N:0:0 $read1[-2]\n$seq1\n+\n$qual1\n";

            ## Generate text file with information about the main and softclipped read 
              print $fh2 "$read1[0]_$readnum\t$read1[2]\t$read1[3]\t$read1[-1]\t$read1[4]\t$read1[5]\t$seq1\t$qual1\t"; 
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
                    print $fh2 "$read3[2]\t$read3[3]\t$read3[-1]\t$read3[4]\t$read3[5]\t$read3[9]\t$olap\t$xalign\n";
                    print $fh3 "$read3[2]\t$read3[3]\t$read3[-1]\t$read3[4]\t$read3[5]\t$read3[9]\t$olap\t$xalign\n";
               }else{
                    print $fh2 "-\t-\t-\t-\t-\t-\t-\t-\n";
                    print $fh3 "-\t-\t-\t-\t-\t-\t-\t-\n";
                }
      }else{
         $count_prox++;
      }
    }

  }

  return($count_re,$count_restr,$count_restrSC,$count_prox, $nonHic_chimera,$count_unmapped); ## counter
}

sub re_comprehensive_search_V1 {
  my($a,$b,$c,$d,$e,$len_RE, $min_olap,$count_re,$count_restr,$count_restrSC) =@_;
  my @read1 = @{$a};
  my @len1 = @{$b};
  my @ops1 = @{$c};
  my @read3 = @{$d};
  my %chromosomes = %{$e};

  my $is_match=0;
  my $seq="";
  
  if($ops1[0] eq "S" && $len1[0] >=  $min_olap){
    my $x=0;
    my $pos=0;

    if($read1[-1] eq "+"){
      $x= $read1[3];
      $pos = &binarysearch( $x, $chromosomes{$read1[2]});
      if($pos > ($x-$leg) && $pos < ($x+$cutoff)){
        $count_re++;
        return($count_re,$count_restr,$count_restrSC,1);
      }
  
      if(@read3 && $read3[-1] eq "+"){
          $x = $read3[3]+length($read3[9]);
          $pos = &binarysearch( $x , $chromosomes{$read3[2]});
          if($pos > ($x-$cutoff) && $pos < ($x+$leg) ){
            $count_re++;
            return($count_re,$count_restr,$count_restrSC,1);
          }
      }
      if(@read3 && $read3[-1] eq "-"){
          $x= $read3[3];
          $pos = &binarysearch( $x , $chromosomes{$read3[2]});
          if($pos > ($x-$leg) && $pos < ($x+$cutoff) ){
            $count_re++;
            return($count_re,$count_restr,$count_restrSC,1);
          }
      }
    }

    if($read1[-1] eq "-"){
      #$x = $read1[3]+length($read1[9])-$len1[0];
      $x = $read1[3];
      $pos = &binarysearch( $x , $chromosomes{$read1[2]});
      if( $pos > ($x-$leg) && $pos < ($x+$cutoff) ){
        $count_re++;
        return($count_re,$count_restr,$count_restrSC,1);
      }

      if(@read3 && $read3[-1] eq "+"){
          $x = $read3[3];
          $pos = &binarysearch( $x , $chromosomes{$read3[2]});
          if($pos > ($x-$leg) && $pos < ($x+$cutoff) ){
            $count_re++;
            return($count_re,$count_restr,$count_restrSC,1);
          }
      }
      if(@read3 && $read3[-1] eq "-"){
          $x= $read3[3]+length($read3[9]);
          $pos = &binarysearch( $x , $chromosomes{$read3[2]});
          if($pos > ($x-$cutoff) && $pos < ($x+$leg) ){
            $count_re++;
            return($count_re,$count_restr,$count_restrSC,1);
          }
      }
    }

    if(RE_within_HardclippedSeq(\@read1,\@len1,\@ops1,\@read3)==1){
      $count_restr++;
      return($count_re,$count_restr,$count_restrSC,1);
    }

    $seq = substr($read1[9], $len1[0]- ($len_RE+2), 2*($len_RE+2));
    if(re_fuzzy_search_stringent($seq)==1){
      $count_restrSC++;
      return($count_re,$count_restr,$count_restrSC,1);
    }else{
      return($count_re,$count_restr,$count_restrSC,0);
    }
  }



  if($ops1[-1] eq "S" and $len1[-1] >=  $min_olap){
    my $x=0;
    my $pos=0;

    if($read1[-1] eq "+"){
      $x= $read1[3]+length($read1[9])-$len1[-1];
      $pos = &binarysearch( $x, $chromosomes{$read1[2]});
      if($pos > ($x-$cutoff) && $pos < ($x+$leg)){
        $count_re++;
        return($count_re,$count_restr,$count_restrSC,1);
      }
      if(@read3 && $read3[-1] eq "+"){
        $x = $read3[3];
        $pos = &binarysearch( $x, $chromosomes{$read3[2]});
        if($pos > ($x-$leg) && $pos < ($x+$cutoff)){
          $count_re++;
          return($count_re,$count_restr,$count_restrSC,1);
        }
      }
      if(@read3 && $read3[-1] eq "-"){
        $x= $read3[3]+length($read3[9]);
        $pos = &binarysearch( $x, $chromosomes{$read3[2]});
        if($pos > ($x-$cutoff) && $pos < ($x+$leg)){
          $count_re++;
          return($count_re,$count_restr,$count_restrSC,1);
        }      
      }
    }

    if($read1[-1] eq "-"){
      $x= $read1[3]+length($read1[9])-$len1[-1];
      $pos = &binarysearch( $x, $chromosomes{$read1[2]});
      if($pos > ($x-$leg) && $pos < ($x+$cutoff)){
        $count_re++;
        return($count_re,$count_restr,$count_restrSC,1);
      }
      if(@read3 && $read3[-1] eq "+"){
        $x = $read3[3]+length($read3[9]);
        $pos = &binarysearch( $x, $chromosomes{$read3[2]});
        if($pos > ($x-$leg) && $pos < ($x+$cutoff)){
          $count_re++;
          return($count_re,$count_restr,$count_restrSC,1);
        }
      }
      if(@read3 && $read3[-1] eq "-"){
        $x= $read3[3];
        $pos = &binarysearch( $x, $chromosomes{$read3[2]});
        if($pos > ($x-$cutoff) && $pos < ($x+$leg)){
          $count_re++;
          return($count_re,$count_restr,$count_restrSC,1);
        }      
      }
    }

    if(RE_within_HardclippedSeq(\@read1,\@len1,\@ops1,\@read3)==1){
      $count_restr++;
      return($count_re,$count_restr,$count_restrSC,1);
    }

    $seq = substr($read1[9], -($len1[-1]+($len_RE+2)),2*($len_RE+2));
    if(re_fuzzy_search_stringent($seq)==1){
      $count_restrSC++;
       return($count_re,$count_restr,$count_restrSC,1);
    }else{
      return($count_re,$count_restr,$count_restrSC,0);
    }
  }

  return($count_re,$count_restr,$count_restrSC,1);
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
  elsif( abs($read1[3]-$read3[3]) <= 512 && $read1[2] eq $read3[2] ){  
  ## Tetracutter cuts every 256bp on an average, hence, minimum resolution 256bp. Any two fragments within 2times the distance are likely to be non-informative
    return(1);
  }
  else{ ## for any other unforseen cases, return 0
    return(0);
  }
}


# binary search of the nearest RE site along the direction of sequencing
# Returns distance between x and it's nearest RE site
sub binarysearch {
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


sub re_fuzzy_search {
  my $seq= shift;

    ###  Fuzzy pattern search!
    ###  Carry out fuzzy search by allowing maximum 3 insertions between the RE sites

   my @patterns = (qr/GATCGATC/,
          qr/GATC[GATC]{0,3}GATC/,
          qr/GATC[ATCG]{0,1}ATC/,
          qr/GAT[ATGC]{0,1}GATC/
   );

  my $is_match=0;

  if( $seq ~~ @patterns ) {
    $is_match++;
    }; 
    return($is_match);
}

sub re_fuzzy_search_stringent {
  my $seq= shift;

  ##  Fuzzy pattern search!
  ##  Dpn2/ Mbo1 cuts to create 4bp overhang. At ligation junction, there should be 8bp signature
  ##  For fuzzy analysis, allow deletion of maximum 2 bases
  ##  Carry out fuzzy search over a 12bp window centered on softclipped junction
  my @patterns = ( qr/GATCGATC/,
                   qr/GATC[GATC]{0,3}GATC/,
                   qr/GATC/,
                   qr/GATC[ATCG]{0,1}ATC/,
                   qr/GAT[ATGC]{0,1}GATC/
  );

  my $is_match=0;

  if( $seq ~~ @patterns ) {
    $is_match++;
    };
    return($is_match);
}


sub flag2binary {
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


sub reverse_complement {
        my $seq = shift;

        my $rcseq = reverse($seq);
        $rcseq =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;

        return $rcseq;
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
  if(length($seq)>3){
     if(re_fuzzy_search_stringent($seq)==0){
        return(0);
      }else{
        return(1);
      }
  }else{
      return(0);
  }
}



exit;

#######################################################################################################################################################################
############## additional subroutines, not a part of the pipeline written above!
#######################################################################################################################################################################
=head


# binary search of the nearest RE site along the direction of sequencing
# Returns distance between x and it's nearest RE site
sub binarysearch {
  my ($x, $a, $upper) = @_;            # search for x in array a
  my ($l, $u) = (0, @$a - 1);          # lower, upper end of search interval
  my $i;                               # index of probe
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
      my $q1 = abs($a->[$i] - $x);
      my $q2 = abs($a->[$i+1] - $x);
      #if($q1 > $q2){
      #  return($q2); # found
      #  }else{
      #    return($q1); # found
      #  }
      if($upper==1){
        #return($q2);
        return($a->[$i+1]);
      }else{
        #return($q1);
        return($a->[$i+1]);
      }
      
    }
  }

  my $q1 = abs($a->[$l] - $x);
  my $q2 = abs($a->[$l-1] - $x);
  #if($q1 > $q2){
  #  return($q2); # found
  #}else{
  #  return($q1); # found
  #}
  if($upper==1){
    #return($q1);
    return($a->[$l]);
  }else{
    #return($q2);
    return($a->[$l-1]);
  }


  #return($l); #original         
}

sub lc_substr {
  my ($str1, $str2) = @_; 
  my $l_length = 0; # length of longest common substring
  my $len1 = length $str1; 
  my $len2 = length $str2;

  # $str1 : longest, $str2 shorter
  my @char1 = (undef, split(//, $str1)); # indexed from 1
  my @char2 = (undef, split(//, $str2)); 
  my @lc_suffix; # check table
  for my $n1 ( 1 .. $len1 ) { 
    for my $n2 ( 1 .. $len2 ) { 
      if ($char1[$n1] eq $char2[$n2]) {
         $lc_suffix[$n1-1][$n2-1] ||= 0;
         $lc_suffix[$n1][$n2] = $lc_suffix[$n1-1][$n2-1] + 1;
         if ($lc_suffix[$n1][$n2] > $l_length) {
           $l_length = $lc_suffix[$n1][$n2];
         }
         if($l_length >= 14){
          return(1);
         }  ## to reduce computation time, break iif the common minimum length is >14. This parameter is not optimized

      }
    }
  }   
  
  if($l_length < 14){
      # $str1 : longest, $str2 shorter
      # reverse complement $str2
      $str2 = reverse_complement($str2);
      $l_length = 0; # length of longest common substring
      $len1 = length $str1; 
      $len2 = length $str2;
      @char1 = (undef, split(//, $str1)); # indexed from 1
      @char2 = (undef, split(//, $str2)); 
      undef(@lc_suffix); # check table
      for my $n1 ( 1 .. $len1 ) { 
        for my $n2 ( 1 .. $len2 ) { 
          if ($char1[$n1] eq $char2[$n2]) {
             $lc_suffix[$n1-1][$n2-1] ||= 0;
             $lc_suffix[$n1][$n2] = $lc_suffix[$n1-1][$n2-1] + 1;
             if ($lc_suffix[$n1][$n2] > $l_length) {
               $l_length = $lc_suffix[$n1][$n2];
             }
             if($l_length >= 14){
              return(1);
             }  ## to reduce computation time, break iif the common minimum length is >14. This parameter is not optimized
          }
        }
      }   
  }
  if($l_length<14){
    return(0);
  }
}

sub re_comprehensive_search {
  my($a,$b,$c,$d,$e,$len_RE, $min_olap) =@_;
  my @read1 = @{$a};
  my @len1 = @{$b};
  my @ops1 = @{$c};
  my @read3 = @{$d};
  my %chromosomes = %{$e};

  my $cutoff = 101;
  my $is_match=0;
  my $seq="";

  ## Primary and secondary pattern search is used here with hope to reduce computational time in mapping nearest upstream RE site

  if($ops1[0] eq "S" && $len1[0] >=  $min_olap){
    ## stringent fuzzy search
    $seq = substr($read1[9], ($len1[0]-($len_RE+2)),2*($len_RE+2));
    if(re_fuzzy_search_stringent($seq)==1){
      return(1);
    }

    ## comprehensive CIGAR based search
    my $upper=0;
    my $dist=200;
    
    if(@read3){
     ## 1) Is there any GATC site within last 2*$len_RE bases of softclipped read ?
          $seq = substr($read1[9], $len1[0]-2*$len_RE,2*$len_RE);
          if(re_fuzzy_search_stringent($seq)==1){
              return(1);
          }
     ## 2) Is there any GATC site within last 2*$len_RE bases of supplimentary alignment of softclipped read ?
          $seq = $read1[9];
          if($read3[-1] eq "-"){
            $seq = reverse_complement($seq);
          }
          $seq = substr($seq, length($seq)-(2*$len_RE),2*$len_RE);
          if(re_fuzzy_search_stringent($seq)==1){
              return(1);
          }


     ## 3) Is there any GATC site within 100bases of the softclipped read on reference genome, in the direction of softclipping junction?
          $upper=1;
          $dist = &binarysearch( ($read3[3]+length($read3[9])), $chromosomes{$read3[2]},$upper); ## for both orientations, this holds true
          if($dist <= $cutoff){
            return(1);
          }

          $upper=0;
          $dist = &binarysearch( $read1[3], $chromosomes{$read1[2]},$upper);
          if($dist <= $cutoff){
            return(1);
          }else{  #  print "$read1[2] $read1[3] $read1[4] $read1[5] $dist $read1[9] $read3[9] $read3[2] $read3[3]\n";
            return(0);
          }
    }else{
          $upper=0;
          $dist = &binarysearch( $read1[3], $chromosomes{$read1[2]},$upper);
          if($dist < $cutoff){
            return(1);
          }else{
            return(0);
          }
    }
  }


  if($ops1[-1] eq "S" and $len1[-1] >=  $min_olap){
   ## stringent fuzzy search
    $seq = substr($read1[9], -($len1[-1]+($len_RE+2)),2*($len_RE+2));
    if(re_fuzzy_search_stringent($seq)==1){
      return(1);
    }

    ## comprehensive CIGAR based search
    my $upper=1;
    my $dist=200;
    
    if(@read3){
     ## 1) Is there any GATC site within first 2*$len_RE bases of softclipped read ?
         $seq = substr($read1[9], -$len1[-1],2*$len_RE);
         if(re_fuzzy_search_stringent($seq)==1){
              return(1);
          }
     ## 2) Is there any GATC site within first 2*$len_RE bases of supplimentary alignment of softclipped read ?
          $seq = $read3[9];
          if($read3[-1] eq "-"){
            $seq = reverse_complement($seq);
          }
          $seq = substr($seq, 0,2*$len_RE);
          if(re_fuzzy_search_stringent($seq)==1){
              return(1);
          }

     ## 3) Is there any GATC site within 100bases of the softclipped read on reference genome?
          $upper=0;
          $dist = &binarysearch( $read3[3], $chromosomes{$read3[2]},$upper); ## for both orientations, this holds true
          if($dist <= $cutoff){
            return(1);
          }

          $upper=1;
          $dist = &binarysearch( ($read1[3]+length($read1[9])-$len1[-1]), $chromosomes{$read1[2]},$upper);   ## get the softclipping junction point 
          if($dist <= $cutoff){
            return(1);
          }else{
            return(0);
          }
    }else{
          $upper=1;
          $dist = &binarysearch( ($read1[3]+length($read1[9])-$len1[-1]), $chromosomes{$read1[2]},$upper);
          if($dist < $cutoff){
            return(1);
          }else{
            return(0);
          }
    } 
  }

  return(1);
}

sub re_comprehensive_search_V1 {
  my($a,$b,$c,$d,$e,$len_RE, $min_olap,$count_re,$count_restr,$count_restrSC) =@_;
  my @read1 = @{$a};
  my @len1 = @{$b};
  my @ops1 = @{$c};
  my @read3 = @{$d};
  my %chromosomes = %{$e};

  my $cutoff = 12;
  my $leg = 50;
  my $is_match=0;
  my $seq="";
  
  if($ops1[0] eq "S" && $len1[0] >=  $min_olap){
    my $x=0;
    my $pos=0;

    if($read1[-1] eq "+"){
      $x= $read1[3];
      $pos = &binarysearch( $x, $chromosomes{$read1[2]});
      if($pos > ($x-$leg) && $pos < ($x+$cutoff)){
        $count_re++;
        return($count_re,$count_restr,$count_restrSC,1);
      }
    }
    if($read1[-1] eq "-"){
      #$x = $read1[3]+length($read1[9])-$len1[0];
      $x = $read1[3];
      $pos = &binarysearch( $x , $chromosomes{$read1[2]});
      if( $pos > ($x-$leg) && $pos < ($x+$cutoff) ){
        $count_re++;
        return($count_re,$count_restr,$count_restrSC,1);
      }
    }
    if(@read3){
      if($read3[-1] eq "+"){
        $x = $read3[3]+length($read3[9]);
        $pos = &binarysearch( $x , $chromosomes{$read3[2]});
        if($pos > ($x-$cutoff) && $pos < ($x+$leg) ){
          $count_re++;
          return($count_re,$count_restr,$count_restrSC,1);
        }
      }
      if($read3[-1] eq "-"){
        $x= $read3[3]+length($read3[9]);
        $pos = &binarysearch( $x , $chromosomes{$read3[2]});
        if($pos > ($x-$cutoff) && $pos < ($x+$leg) ){
          $count_re++;
          return($count_re,$count_restr,$count_restrSC,1);
        }
      }
    }
    $seq = substr($read1[9], $len1[0]- ($len_RE+2), 2*($len_RE+2));
    if(re_fuzzy_search_stringent($seq)==1){
      $count_restr++;
      if($count_restr<5){
        print join("\t",@read1),"\n";
        print join("\t",@read3),"\n\n";
      }
      return($count_re,$count_restr,$count_restrSC,1);
    }
    $seq = substr($read1[9], $len1[0]-12, 12);
    if(re_fuzzy_search_stringent($seq)==1){
      $count_restrSC++;
      #if($count_restrSC==2){
      #  print join("\t",@read1),"\n";
      #  print join("\t",@read3),"\n\n";
      #}
      return($count_re,$count_restr,$count_restrSC,1);
    }else{
      return($count_re,$count_restr,$count_restrSC,0);
    }
  }

  if($ops1[-1] eq "S" and $len1[-1] >=  $min_olap){
    my $x=0;
    my $pos=0;

    if($read1[-1] eq "+"){
      $x= $read1[3]+length($read1[9])-$len1[-1];
      $pos = &binarysearch( $x, $chromosomes{$read1[2]});
      if($pos > ($x-$cutoff) && $pos < ($x+$leg)){
        $count_re++;
        return($count_re,$count_restr,$count_restrSC,1);
      }
    }
    if($read1[-1] eq "-"){
      $x= $read1[3]+length($read1[9])-$len1[-1];
      $pos = &binarysearch( $x, $chromosomes{$read1[2]});
      if($pos > ($x-$cutoff) && $pos < ($x+$leg)){
        $count_re++;
        return($count_re,$count_restr,$count_restrSC,1);
      }
    }
   if(@read3){
      if($read3[-1] eq "+"){
        $x = $read3[3];
        $pos = &binarysearch( $x, $chromosomes{$read3[2]});
        if($pos > ($x-$leg) && $pos < ($x+$cutoff)){
          $count_re++;
          return($count_re,$count_restr,$count_restrSC,1);
        }
      }
      if($read3[-1] eq "-"){
        $x= $read3[3]+length($read3[9]);
        #//$x= $read3[3];
        $pos = &binarysearch( $x, $chromosomes{$read3[2]});
        if($pos > ($x-$cutoff) && $pos < ($x+$leg)){
          $count_re++;
          return($count_re,$count_restr,$count_restrSC,1);
        }      
      }
    }
    $seq = substr($read1[9], -($len1[-1]+($len_RE+2)),2*($len_RE+2));
    if(re_fuzzy_search_stringent($seq)==1){
      $count_restr++;
      if($count_restr<5){
        print join("\t",@read1),"\n";
        print join("\t",@read3),"\n\n";
      }
       return($count_re,$count_restr,$count_restrSC,1);
    }
    $seq = substr($read1[9], -($len1[-1]),12);
    if(re_fuzzy_search_stringent($seq)==1){
      $count_restrSC++;
      #if($count_restrSC==1){
      #  print join("\t",@read1),"\n";
      #  print join("\t",@read3),"\n\n";
      #}

      return($count_re,$count_restr,$count_restrSC,1);
    }else{
      return($count_re,$count_restr,$count_restrSC,0);
    }
  }

  return($count_re,$count_restr,$count_restrSC,1)
}

=cut













