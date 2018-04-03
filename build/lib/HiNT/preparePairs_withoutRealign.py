#This script is trying to extract the information of nonHi-C related chrimeric reads that gnerated from extract_nonHiC.pl (perl script) to pairs format.
import os,sys
from multiprocessing import Pool

chimericReads = sys.argv[1]
workdir = sys.argv[2]
jobName = sys.argv[3]

WORKING_FOlDER='%s/'%workdir
print WORKING_FOlDER


def makePairFormat(read):
	outres = False

	line = read.strip().split('\t')
	if len(line) == 16:
		readID,chrom1,pos1,strand1,score1,cigar1,sequence1,seq_qual1,chrom2,pos2,strand2,score2,cigar2,sequence2,flag,alternatives = line
	else:
		readID,chrom1,pos1,strand1,score1,cigar1,sequence1,chrom2,pos2,strand2,score2,cigar2,sequence2 = line
	if cigar1 != '*' and cigar2 != '-':
		start1 = pos1
		end1 = int(pos1) + 1
		start2 = pos2
		end2 = int(pos2) + 1
		res = [readID,chrom1,start1,str(end1),chrom2,start2,str(end2),strand1,strand2,cigar1,cigar2]
		outres = '\t'.join(res) + '\n'
		#outf.write('\t'.join(res) + '\n')

		if alternatives != '-':
			alternatives = alternatives.lstrip('XA:Z:').rstrip(';').split(';')
			for i in range(len(alternatives)):
				[chrom2,strand_pos2,cigar2,mismatch2] = alternatives[i].split(',')
				stand2 = strand_pos2[0]
				start2 = int(strand_pos2[1:])
				end2 = start2 + 1
				res = [readID,chrom1,start1,str(end1),chrom2,str(start2),str(end2),strand1,strand2,cigar1,cigar2]
				outres += '\t'.join(res) + '\n'
	return outres

print "1, get regions Infomation..."
#regions = get_regions(chromlengthf,breakpointsf,windowsize=5)
#print regions

#regions_fas = gnrt_region_ref(genome_ref, regions)
#index_region_ref(regions_fas)

print "2, writing fastq files, parallele processing..."
inf = open(chimericReads).xreadlines()
Infos = []

output = WORKING_FOlDER + (jobName + '.pairs')
pool = Pool(16)
with open(output,'w') as outf:
	for res in pool.imap(makePairFormat, inf):
		if res != False:
			outf.write(res)
inf.close()
