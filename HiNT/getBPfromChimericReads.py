import os,sys
from pkg_resources import resource_filename
from multiprocessing import Pool
import numpy as np
import pypairix
from HiNT.corelib import *

def extrac_nonHiCchimeric(chimericPairsam,restrictionSites,restrictionEnzyme,outputPrefix):
	extractionScript = resource_filename('HiNT', 'externalScripts/extract_nonHiC.pl')
	command = "perl %s -psam %s -sites %s -re %s -outprefix %s"%(extractionScript, chimericPairsam, restrictionSites, restrictionEnzyme, outputPrefix)
	run_cmd(command)
	outputChimeras = outputPrefix+'_3dedup.txt'
	return outputChimeras

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

def getchromsize(chromlengthf):
	chromSizeInfo = {}
	inf = open(chromlengthf)
	for line in inf:
		line = line.strip().split('\t')
		chrom,size = line
		chromSizeInfo[chrom] = int(size)
	return chromSizeInfo

def readFilteredBPs(breakpointsf,chromSizeInfo):
	bpInfo = {}
	inf = open(breakpointsf)
	for line in inf:
		if line.startswith('chrom'):
			pass
		else:
			line = line.strip().split('\t')
			chrompair,chrom1,bp1,chrom2,bp2 = line
			try:
				bpInfo[chrompair].append([(chrom1,bp1),(chrom2,bp2)])
			except:
				bpInfo[chrompair] = [[(chrom1,bp1),(chrom2,bp2)]]
	#print bpInfo
	return bpInfo

def is_qualified_clipped(cigar, cutoff_len):
	l=len(cigar)
	signal=[]
	lenth=[]
	temp=""
	for i in range(l):
		if cigar[i]>="0" and cigar[i]<="9":
			temp=temp+cigar[i]
		else:
			signal.append(cigar[i])
			try:
				lenth.append(int(temp))
			except ValueError:
				print("Error: ", cigar, temp)
			temp=""

	b_hardclip=False
	cnt_m=0
	clip_flag=0
	left_clip_len=0
	right_clip_len=0

	for i in range(len(signal)):
		if signal[i]=="M":
			cnt_m=cnt_m+lenth[i]

	if (signal[0]=="S" or signal[0]=="H") and lenth[0]>=cutoff_len:#left-clip
		clip_flag=1
		left_clip_len=lenth[0]
		if signal[0]=="H":
			b_hardclip=True
	if (signal[len(signal)-1]=="S" or signal[len(signal)-1]=="H") and lenth[len(signal)-1]>=cutoff_len: #right-clip
		clip_flag=clip_flag+2 #if this is 3, then both side clipped
		right_clip_len=lenth[len(signal)-1]
		if signal[len(signal)-1]=="H":
			b_hardclip=True

	return b_hardclip, cnt_m, clip_flag, left_clip_len, right_clip_len

def getclipPos_pair(read,m_clip_pos,reverse=False):
	CLIP_CUTOFF=13
	if reverse == True:
		ID,chromb,posb,endb,chroma,posa,enda,strandb,stranda,cigarb,cigara = read
	else:
		ID,chroma,posa,enda,chromb,posb,endb,stranda,strandb,cigara,cigarb = read
	b_hardclip_a, cnt_m_a, clip_flag_a, left_clip_len_a, right_clip_len_a = is_qualified_clipped(cigara, CLIP_CUTOFF)
	b_hardclip_b, cnt_m_b, clip_flag_b, left_clip_len_b, right_clip_len_b = is_qualified_clipped(cigarb, CLIP_CUTOFF)
	clip_pos_a=int(posa)
	if right_clip_len_a>0:#right clip
		clip_pos_a=int(posa)+cnt_m_a-1
	clip_pos_b = int(posb)
	if right_clip_len_b>0:#right clip
		clip_pos_b=int(posb)+cnt_m_b-1
	pair_pos = (clip_pos_a,clip_pos_b)
	if pair_pos not in m_clip_pos:
		m_clip_pos[pair_pos]=0
	m_clip_pos[pair_pos]=m_clip_pos[pair_pos]+1

	return m_clip_pos

def singleSide_clip_pos_calculation(m_clip_pos_forward,m_clip_pos_reverse,window_size=10):
	#clip_pos_a and clip_pos_b will be validated only when reads that chimeric in chr1~chr8 has such clip pos, and reads that chimeric in chr8~chr1 have such clip position too.
	clip_Pos_statistics = {}

	forward_keys = list(m_clip_pos_forward.keys())
	reverse_keys = list(m_clip_pos_reverse.keys())
	for i in range(len(forward_keys)):
		posa,posb = forward_keys[i]
		t = 0
		for j in range(-window_size,window_size+1):
			temp_posa = posa + j
			for l in range(-window_size,window_size+1):
				temp_posb = posb + l
				temp_pos = (temp_posa,temp_posb)
				if temp_pos in reverse_keys:
					t += 1
					if temp_pos not in clip_Pos_statistics:
						clip_Pos_statistics[temp_pos] = [m_clip_pos_forward[(posa,posb)],m_clip_pos_reverse[temp_pos],m_clip_pos_forward[(posa,posb)]+m_clip_pos_reverse[temp_pos]]
					else:
						clip_Pos_statistics[temp_pos] = [x + y for x,y in zip(clip_Pos_statistics[temp_pos], [m_clip_pos_forward[(posa,posb)],m_clip_pos_reverse[temp_pos],m_clip_pos_forward[(posa,posb)]+m_clip_pos_reverse[temp_pos]])]
		if t == 0:
			if (posa, posb) not in clip_Pos_statistics:
				clip_Pos_statistics[(posa,posb)] = [m_clip_pos_forward[(posa,posb)],0,m_clip_pos_forward[(posa,posb)]+0]
			else:
				clip_Pos_statistics[(posa,posb)] = [x + y for x,y in zip(clip_Pos_statistics[(posa,posb)], [m_clip_pos_forward[(posa,posb)],0,m_clip_pos_forward[(posa,posb)]+0])]

	for i in range(len(reverse_keys)):
		posa,posb = reverse_keys[i]
		m = 0
		for j in range(-window_size,window_size+1):
			temp_posa = posa + j
			for l in range(-window_size,window_size+1):
				temp_posb = posb + l
				temp_pos = (temp_posa,temp_posb)
				if temp_pos in forward_keys:
					m += 1
					pass
		if m == 0:
			if (posa,posb) not in clip_Pos_statistics:
				clip_Pos_statistics[(posa,posb)] = [0, m_clip_pos_reverse[(posa,posb)], m_clip_pos_reverse[(posa,posb)]+0]
			else:
				clip_Pos_statistics[(posa,posb)] = [x + y for x,y in zip(clip_Pos_statistics[(posa,posb)], [0,m_clip_pos_reverse[(posa,posb)],m_clip_pos_reverse[(posa,posb)]+0])]

	return clip_Pos_statistics

def further_filtering(clip_Pos_statistics,bp1,bp2):
	#The input clip_pos_dic can be either double validated clip positions or single side validated clip positions
	resolution = 100000
	filtered_clip_pos_statistics = {}
	print(clip_Pos_statistics)
	for cp in clip_Pos_statistics:
		if ((clip_Pos_statistics[cp][-2] >= 1) and (clip_Pos_statistics[cp][-3] >= 1)) or (clip_Pos_statistics[cp][-1] >= 2):
			filtered_clip_pos_statistics[cp] = clip_Pos_statistics[cp]
		elif (abs(cp[0] - int(bp1)*resolution) <= 2*resolution) and (abs(cp[1]-int(bp2)*resolution) <= 2*resolution):
			filtered_clip_pos_statistics[cp] = clip_Pos_statistics[cp]
		else:
			pass
	return filtered_clip_pos_statistics

def getReads(bpInfo,chromSizeInfo,chimericReadPairs,outdir,outname,windowsize=5,window_size=10):
	#windowsize is for the break point region detected from Hi-C, window_size is for the sliding window searching size
	tb = pypairix.open(chimericReadPairs)
	output = os.path.join(outdir,outname + '.supportedReads.txt')
	outf = open(output,'w')
	output3 = os.path.join(outdir,outname + '.singleValidated_pairs_BPs.txt')
	outf3 = open(output3,'w')
	resolution = 100000
	for chrompair in bpInfo:
		print(chrompair)
		chrompair_res1 = []
		chrompair_res2 = []
		chrompair_res3 = []
		bps = bpInfo[chrompair]
		supportedReads =[]
		for bp_detail in bps:
			chrom1Info,chrom2Info = bp_detail
			chr1,bp1 = chrom1Info
			chr2,bp2 = chrom2Info
			start = max(0,(int(bp1)-windowsize)*resolution)
			end = min((int(bp1)+windowsize)*resolution,chromSizeInfo[chrom1Info[0]])
			start2 = max(0,(int(bp2)-windowsize)*resolution)
			end2 = min((int(bp2)+windowsize)*resolution,chromSizeInfo[chrom2Info[0]])
			it1 = tb.query2D(chrom1Info[0],start,end,chrom2Info[0],start2,end2)
			m_clip_pos_forward = {}
			m_clip_pos_reverse = {}
			tempreads = []
			for read in it1:
				readinfo = read[1:]
				if readinfo not in tempreads:
					m_clip_pos_forward = getclipPos_pair(read,m_clip_pos_forward)
					tempreads.append(readinfo)
				else:
					pass
				if readinfo not in supportedReads:
					supportedReads.append(readinfo)
			it2 = tb.query2D(chrom2Info[0],start2,end2,chrom1Info[0],start,end)
			for read in it2:
				readinfo = read[1:]

				if readinfo not in tempreads:
					m_clip_pos_reverse = getclipPos_pair(read,m_clip_pos_reverse,reverse=True)
					tempreads.append(readinfo)
				else:
					pass
				if readinfo not in supportedReads:
					supportedReads.append(readinfo)
			clip_pos_statistics = singleSide_clip_pos_calculation(m_clip_pos_forward,m_clip_pos_reverse,window_size=10)
			filtered_clip_pos_statistics = further_filtering(clip_pos_statistics,bp1,bp2)
			reg_a = '_'.join([str(start),str(end)])
			reg_b = '_'.join([str(start2),str(end2)])
			for cp in filtered_clip_pos_statistics:
				#print chrompair_res
				pos_a,pos_b = cp
				hit_forward,hit_reverse,hit_all = filtered_clip_pos_statistics[cp]
				if hit_all == 0:
					pass
				else:
					if cp not in chrompair_res1:
						chrompair_res1.append(cp)
						res = [chrom1Info[0],reg_a,str(pos_a),chrom2Info[0],reg_b,str(pos_b)] + [str(hit_forward),str(hit_reverse),str(hit_all)]
						outf3.write('\t'.join(res)+"\n")
					else:
						pass
		for sr in supportedReads:
			outf.write('\t'.join(sr) + '\n')

	outf3.close()
	outf.close()
	return output3

def getBPfromChimeras(chimericPairsam,restrictionSites,restrictionEnzyme,outdir,name,chromlengthf,validbreakpointsf,pairixpath):
	chimericDir = os.path.join(outdir,'BPfromChimeric')
	if not os.path.isdir(chimericDir):
		os.mkdir(chimericDir)
	chimerasPrefix = os.path.join(chimericDir,name+'_nonHiC_chimeric')
	outputChimeras = extrac_nonHiCchimeric(chimericPairsam,restrictionSites,restrictionEnzyme,chimerasPrefix)
	inf = open(outputChimeras)
	Infos = []
	outputChimeraPairs = os.path.join(chimericDir,name+'_nonHiC_chimeric_3dedup.pairs')
	pool = Pool(16)
	with open(outputChimeraPairs,'w') as outf:
		for res in pool.imap(makePairFormat, inf):
			if res != False:
				outf.write(res)
	inf.close()
	sortedChimeras = os.path.join(chimericDir,name+'_nonHiC_chimeric_3dedup.bsorted.pairs.gz')
	command1 = "sort -k2,2 -k5,5 -k3,3n -k 6,6n %s | bgzip -c > %s"%(outputChimeraPairs,sortedChimeras)
	print(command1)
	run_cmd(command1)
	command2 = "%s -f -s2 -d5 -b3 -e4 -u6 -v7 %s"%(pairixpath,sortedChimeras)
	print(command2)
	run_cmd(command2)
	chromInfo = getchromsize(chromlengthf)
	bpInfo = readFilteredBPs(validbreakpointsf,chromInfo)
	BP_fromChimerasfile = getReads(bpInfo,chromInfo,sortedChimeras,outdir,name,windowsize=5,window_size=10)
	return BP_fromChimerasfile
