#This script will identify breakpoints with single base pair resolution by using chimeric reads. singleValidated_pairs_BPs will be used to summarize and be as the candidate BPs
import os,sys
import pypairix
from random import randint, choice,seed


chromlengthf = sys.argv[1] #chromosome length information, 2 columns: chromosome, length
breakpointsf = sys.argv[2] #breakpoints detected from Hi-C inter-chromosomal interaction matrices with 100kb resolution
chimericReadPairs = sys.argv[3] #nonHi-C related chimeric reads with pairs format
outdir = sys.argv[4] #output directory
outname = sys.argv[5] #name out output file
resolution = int(sys.argv[6]) #resolution of breakpoints from Hi-C, default should be 100kb


def getchromsize(chromlengthf):
	chromSizeInfo = {}
	inf = open(chromlengthf)
	for line in inf:
		line = line.strip().split('\t')
		chrom,size = line
		chromSizeInfo[chrom] = int(size)
	return chromSizeInfo

def readFilteredBPs(breakpointsf,chromSizeInfo):
	#breakpointsf is a list of breakpoints in 100kb resolution, format like: chr10_chr18; chr10; 424; chr18; 251
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

	print bpInfo
	return bpInfo

##Check whether reads is qualified clipped reads
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
                print "Error: ", cigar, temp
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

def getclipPos(read,m_clip_pos_a,m_clip_pos_b):
	CLIP_CUTOFF=13
	ID,chroma,posa,enda,chromb,posb,endb,stranda,strandb,cigara,cigarb = read
	b_hardclip_a, cnt_m_a, clip_flag_a, left_clip_len_a, right_clip_len_a = is_qualified_clipped(cigara, CLIP_CUTOFF)
	b_hardclip_b, cnt_m_b, clip_flag_b, left_clip_len_b, right_clip_len_b = is_qualified_clipped(cigarb, CLIP_CUTOFF)

	clip_pos_a=int(posa)
	if right_clip_len_a>0:#right clip
		clip_pos_a=int(posa)+cnt_m_a-1

	if clip_pos_a not in m_clip_pos_a:
		m_clip_pos_a[clip_pos_a]=0
	m_clip_pos_a[clip_pos_a]=m_clip_pos_a[clip_pos_a]+1


	clip_pos_b = int(posb)
	if right_clip_len_b>0:#right clip
		clip_pos_b=int(posb)+cnt_m_b-1

	if clip_pos_b not in m_clip_pos_b:
		m_clip_pos_b[clip_pos_b]=0
	m_clip_pos_b[clip_pos_b]=m_clip_pos_b[clip_pos_b]+1

	return m_clip_pos_a,m_clip_pos_b

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

def doubleValidation_pair(m_clip_pos_forward,m_clip_pos_reverse,window_size=10):
	#clip_pos_a and clip_pos_b will be validated only when reads that chimeric in chr1~chr8 has such clip pos, and reads that chimeric in chr8~chr1 have such clip position too.
	common_clip_Pos = {}

	forward_keys = m_clip_pos_forward.keys()
	reverse_keys = m_clip_pos_reverse.keys()
	for i in range(len(forward_keys)):
		posa,posb = forward_keys[i]
		for j in range(-window_size,window_size+1):
			temp_posa = posa + j
			for l in range(-window_size,window_size+1):
				temp_posb = posb + l
				temp_pos = (temp_posa,temp_posb)
				if temp_pos in reverse_keys:
					if temp_pos not in common_clip_Pos:
						common_clip_Pos[temp_pos] = [m_clip_pos_forward[(posa,posb)],m_clip_pos_reverse[temp_pos],m_clip_pos_forward[(posa,posb)]+m_clip_pos_reverse[temp_pos]]
					else:
						common_clip_Pos[temp_pos] = [x + y for x,y in zip(common_clip_Pos[temp_pos], [m_clip_pos_forward[(posa,posb)],m_clip_pos_reverse[temp_pos],m_clip_pos_forward[(posa,posb)]+m_clip_pos_reverse[temp_pos]])]
	#print common_clip_Pos
	return common_clip_Pos

def singleSide_clip_pos_calculation(m_clip_pos_forward,m_clip_pos_reverse,window_size=10):
	#clip_pos_a and clip_pos_b will be validated only when reads that chimeric in chr1~chr8 has such clip pos, and reads that chimeric in chr8~chr1 have such clip position too.
	clip_Pos_statistics = {}

	forward_keys = m_clip_pos_forward.keys()
	reverse_keys = m_clip_pos_reverse.keys()
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
	filtered_clip_pos_statistics = {}
	print clip_Pos_statistics
	for cp in clip_Pos_statistics:
		if ((clip_Pos_statistics[cp][-2] >= 1) and (clip_Pos_statistics[cp][-3] >= 1)) or (clip_Pos_statistics[cp][-1] >= 2):
			filtered_clip_pos_statistics[cp] = clip_Pos_statistics[cp]
		elif (abs(cp[0] - int(bp1)*resolution) <= 2*resolution) and (abs(cp[1]-int(bp2)*resolution) <= 2*resolution):
			filtered_clip_pos_statistics[cp] = clip_Pos_statistics[cp]
		else:
			pass

	#### The following are for further filtering, it is a little bit strigent, if the reads coverage is not high enough, not recommend to use
	'''
	need_to_delete = []

	if len(filtered_clip_pos_statistics) > 1:
		for cp in filtered_clip_pos_statistics:
			for c in range(1,6):
				if (filtered_clip_pos_statistics[cp][-2]*filtered_clip_pos_statistics[cp][-3] == 0) and (filtered_clip_pos_statistics[cp][-1]<=c) and (len(need_to_delete) < (len(filtered_clip_pos_statistics)-1)):
					if cp not in need_to_delete:
						need_to_delete.append(cp)
						break
				elif len(need_to_delete) >= (len(filtered_clip_pos_statistics)-1):
					break
				else:
					continue

	for ncp in need_to_delete:
		del filtered_clip_pos_statistics[ncp]
	'''
	return filtered_clip_pos_statistics

#use pairix to get the valid chimeric reads wihtin the defined regions
def getReads(bpInfo,chromSizeInfo,chimericReadPairs,outdir,outnam,windowsize=5,window_size=10):
	#windowsize is for the break point region detected from Hi-C, window_size is for the sliding window searching size
	tb = pypairix.open(chimericReadPairs)
	output = os.path.join(outdir,outname + '.supportedReads.txt') #all supported reads 
	outf = open(output,'w')
	output2 = os.path.join(outdir,outname + '.doubleValidated_pairs_BPs.txt')
	outp2 = open(output2,'w')
	output3 = os.path.join(outdir,outname + '.singleValidated_pairs_BPs.txt')
	outp3 = open(output3,'w')
	output4 = os.path.join(outdir,outname + '.all_pairs_BPs.txt')
	outp4 = open(output4,'w')

	for chrompair in bpInfo:
		print chrompair
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
			#res = [chrom1Info[0],str(start),str(end),chrom2Info[0],str(start2),str(end2)]
			#outf.write('\t'.join(res) + '\n')
			it1 = tb.query2D(chrom1Info[0],start,end,chrom2Info[0],start2,end2)
			#print chrom1Info[0],start,end,chrom2Info[0],start2,end2
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
			common_clip_Pos = doubleValidation_pair(m_clip_pos_forward,m_clip_pos_reverse,window_size=10)
			#print bp1,bp2
			reg_a = '_'.join([str(start),str(end)])
			reg_b = '_'.join([str(start2),str(end2)])
			#print reg_a,reg_b
			#print common_clip_Pos
			for cp in common_clip_Pos:
				#print cp
				#print common_clip_Pos[cp]
				if cp not in chrompair_res2:
					chrompair_res2.append(cp)
					pos_a,pos_b = cp
					hit_forward,hit_reverse,hit_all = common_clip_Pos[cp]
					res = [chrom1Info[0],reg_a,str(pos_a),chrom2Info[0],reg_b,str(pos_b)] + [str(hit_forward),str(hit_reverse),str(hit_all)]
					outp2.write('\t'.join(res)+"\n")
				else:
					pass
			clip_pos_statistics = singleSide_clip_pos_calculation(m_clip_pos_forward,m_clip_pos_reverse,window_size=10)
			filtered_clip_pos_statistics = further_filtering(clip_pos_statistics,bp1,bp2)

			#filtered_clip_pos_statistics = clip_pos_statistics
			#print filtered_clip_pos_statistics

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
						outp3.write('\t'.join(res)+"\n")
					else:
						pass
			for cp in clip_pos_statistics:
				#print chrompair_res
				pos_a,pos_b = cp
				hit_forward,hit_reverse,hit_all = clip_pos_statistics[cp]

				if cp not in chrompair_res3:
					chrompair_res3.append(cp)
					res = [chrom1Info[0],reg_a,str(pos_a),chrom2Info[0],reg_b,str(pos_b)] + [str(hit_forward),str(hit_reverse),str(hit_all)]
					outp4.write('\t'.join(res)+"\n")
				else:
					pass

		for sr in supportedReads:
			outf.write('\t'.join(sr) + '\n')

	outp2.close()
	outp3.close()
	outp4.close()
	outf.close()

def main():
	chromInfo = getchromsize(chromlengthf)
	#print chromInfo
	bpInfo = readFilteredBPs(breakpointsf,chromInfo)
	#print bpInfo
	getReads(bpInfo,chromInfo,chimericReadPairs,outdir,outname,windowsize=5,window_size=10)
if __name__ == '__main__':
	main()
