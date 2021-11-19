#!/usr/bin/python
import pysam
from collections import defaultdict, Counter
from HTSeq import GenomicInterval
from time import time
import itertools
import numpy as np ; np.seterr(all = "ignore")



def filter_reads(pysam_read, P_minMapQual):
	tags = dict(pysam_read.tags)

	if pysam_read.mapping_quality < P_minMapQual:
		return True, 'mapq'
	elif len(pysam_read.query_alignment_sequence) < 15:
		return True, 'short'
	elif "N" in pysam_read.query_alignment_sequence:
		return True, 'N'
	elif pysam_read.is_unmapped:
		return True, 'Unmap'
	elif pysam_read.is_secondary:
		return True, 'Second'
	elif pysam_read.is_duplicate:
		return True, 'Dup'
	elif pysam_read.is_qcfail:
		return True, 'qcfail'
	elif "NH" in tags and tags["NH"] > 1:
		return True, 'nh>1'
	elif pysam_read.is_paired and not pysam_read.is_proper_pair:
		return True, 'unpair'
	elif pysam_read.is_supplementary:
		return True, 'supp'
	else:
		return False, 'good'



def get_read_strand(pysam_read, strandedness):
	
	if strandedness == 'fr-firststrand':
		strand_1, strand_2 = '+', '-'
	elif strandedness == 'fr-secondstrand':
		strand_1, strand_2 = '-', '+'
	elif strandedness == 'fr-unstrand':
		strand_1, strand_2 = '.', '.'

	if pysam_read.is_read2 and pysam_read.is_reverse:      
		r_strand = strand_2			
	elif pysam_read.is_read1 and pysam_read.is_reverse:	     
		r_strand = strand_1
	elif pysam_read.is_read1 and not pysam_read.is_reverse:  
		r_strand = strand_2
	elif pysam_read.is_read2 and not pysam_read.is_reverse:  
		r_strand = strand_1
	else:
		if pysam_read.is_reverse:
			r_strand = strand_1
		else:
			r_strand = strand_2
	return r_strand


def process_CIGAR(pos, cigar):
	block_lengths  = [] 
	soft_clippings = []
	blocks = 0
	for (ctype, length) in cigar:
		if ctype == 0: #M
			blocks += length
		elif ctype == 2: #D
			blocks += length
		elif ctype == 3: #N
			block_lengths.append((pos, pos + blocks))
			pos += (blocks + length)
			blocks = 0
		elif ctype == 4: #S
			soft_clippings.append(length)
	block_lengths.append((pos, pos + blocks))
	return block_lengths, soft_clippings


def get_chroms_from_bam(bam_file):
	bam_handle = pysam.Samfile(bam_file, 'rb')
	bam_chroms = bam_handle.references
	bam_handle.close()
	return bam_chroms



def get_sample_name(bam_file):
	bam_handle = pysam.Samfile(bam_file, 'rb')	
	if "RG" not in bam_handle.header:
		sample_name = bam_file
	elif len(bam_handle.header['RG']) == 1 and 'SM' in bam_handle.header['RG'][0]:
		sample_name = bam_handle.header['RG'][0]['SM']
	else:
		sample_name = bam_file
	return sample_name



def find_read_clusters(bam_file, strandedness, CUTOFFS, chrom, start, end):
	
	read_block_coords = defaultdict(dict)

	bam_handle = pysam.Samfile(bam_file, 'rb')

	for read in bam_handle.fetch(chrom, start = start, end = end, until_eof = True):
		chrom = read.reference_name

		if filter_reads(read, CUTOFFS["P_minMapQual"])[0]:
			continue

		r_blocks, soft_clippings = process_CIGAR(read.pos, read.cigar)

		r_strand = get_read_strand(read, strandedness)

		for i, (sB, eB) in enumerate(r_blocks):
			try:
				read_block_coords[(chrom, r_strand)][sB] += 1
			except:
				read_block_coords[(chrom, r_strand)][sB] = 1
			try:
				read_block_coords[(chrom, r_strand)][eB] -= 1
			except:
				read_block_coords[(chrom, r_strand)][eB] = -1

	bam_handle.close()

	P_minCountTotal = CUTOFFS["DP"]

	Read_Clusters = list()
	RC_cov_cutoff = max(P_minCountTotal/2, 1)
	i = 0

	for (chrom, strand), coords in read_block_coords.items():
		Cum_Cov_block = 0
		Max_Cov_block = 0

		RC_started = False
		prev_coord = 0
		Area = 0.0

		_read_clusters = list()

		for coord in sorted(coords): 
			Cum_Cov_block += read_block_coords[(chrom, strand)][coord]
		
			if Cum_Cov_block > Max_Cov_block:
				Max_Cov_block = Cum_Cov_block

			if Cum_Cov_block >= RC_cov_cutoff and not RC_started:
				start_coord = coord
				Area = 0.0
				RC_started = True
				

			elif Cum_Cov_block < RC_cov_cutoff and RC_started:
				dx = coord - prev_coord
				dy = Cum_Cov_block  
				Area += dx*dy 

				_read_clusters.append((start_coord, coord, Max_Cov_block, Area))
				RC_started = False
				Max_Cov_block = 0
				Area = 0

			else:
				dx = coord - prev_coord
				dy = Cum_Cov_block  
				Area += dx*dy 

			prev_coord = coord 

		for (start, end, Max_Cov, Area) in _read_clusters:
			RC_length = end - start		

			v = (end - start)*Max_Cov
			
			if Max_Cov < P_minCountTotal:
				continue

			if RC_length > 5000:
				sub_RC_length = 3000
				remainder = RC_length % sub_RC_length
				while remainder < 350 and remainder > 0:
					sub_RC_length += 5
					remainder = RC_length % sub_RC_length

				for s_start in range(start, end, sub_RC_length):
					s_end = min(s_start + sub_RC_length, end)
					RC = [i, chrom, strand, s_start, s_end, Max_Cov]
					Read_Clusters.append(RC)
					i += 1
			else:
				RC = [i, chrom, strand, start, end, Max_Cov]
				Read_Clusters.append(RC)
				i += 1

		del _read_clusters
		
	del read_block_coords

	return Read_Clusters


	

class genome_ref_pos():
	
	def __init__(self, ranges, chrom, RC_Start, RC_End, RC_index, Seq, faidx_obj, strand, RC_read_overlap):
		self.ranges = ranges
		self.chrom  = chrom
		self.pos0   = ranges[0][0]
		self.end1   = ranges[-1][1]
		self.RC     = GenomicInterval(chrom, RC_Start, RC_End)
		self.RC_i   = RC_index
		self.Seq    = Seq.upper()
		self.Strand = strand
		self.Read_OV = RC_read_overlap
		self.Faidx  = str(faidx_obj.fetch(chrom, ranges[0][0] + 1, ranges[-1][1])).upper()

	def rc_pos(self, genome_pos):
		RC_Pos = None
		CumLen = 2
		for (s, e) in self.ranges:
			if s <= genome_pos and genome_pos < e:
				RC_Pos = CumLen + genome_pos - s
			CumLen += (e - s)  + 1
		return RC_Pos

	def get_sequence(self, g_start, g_end):
		off = g_start - self.pos0
		assert off >= 0 
		return self.Faidx[off : off + (g_end - g_start)]

	def genome_pos(self, relative_pos):
		CumLen = 2
		RealPos = None
		for (s, e) in self.ranges:
			ExonLen = e - s 
			if relative_pos > CumLen + ExonLen:
				CumLen += ExonLen + 1
			else:
				RealPos = s + relative_pos - CumLen
				break
		if RealPos is None:
			RealPos = self.end1
		return RealPos


def remove_homopolymer_ends(read, blocks):
	READSEQ = read.query_alignment_sequence
	
	i = 0
	lim = len(READSEQ) - 12
	while i < lim and READSEQ[i].upper() == READSEQ[i + 1].upper():
		i += 1

	j = len(READSEQ)
	lim = 12
	while j > lim and READSEQ[j - 1].upper() == READSEQ[j - 2].upper():
		j -= 1

	prefix = i if i > 10 else 0
	sufix  = j if j < (len(READSEQ) - 10) else len(READSEQ)

	if sufix - prefix <= 15:
		prefix, sufix = 0, len(READSEQ)

	new_blocks = []
	rlen = 0
	for (b1, b2) in blocks:
		blen = (b2 - b1)
		b1_new, b2_new = b1, b2

		if rlen + blen <= prefix:
			rlen += blen
			continue
		else:
			b1_new = b1 + max(prefix - rlen, 0)

		if rlen + blen <= sufix:
			pass
		elif rlen > sufix:
			continue
		else:
			b2_new = b1 + (sufix - rlen)

		rlen += blen
		new_blocks.append((b1_new, b2_new))

	aligned_qual = read.query_alignment_qualities[prefix : sufix]
	aligned_seq  = read.query_alignment_sequence[prefix : sufix]	

	return aligned_seq, aligned_qual, new_blocks




def get_RC_reads(bam_file, RC_info, genome, strandedness, CUTOFFS):

	RC_index, chrom, strand, RC_Start, RC_End, MaxCov = RC_info

	SetOfReadNames = defaultdict(dict)	

	_read_blocks = defaultdict(int)
	_unique_seqs = defaultdict(dict)
	_unique_pos  = defaultdict(dict)


	if MaxCov > 2000:
		subsample_rate = int(MaxCov/1000)
	else:
		subsample_rate = 1 


	bam_handle = pysam.Samfile(bam_file, 'rb')

	for ri, read in enumerate(bam_handle.fetch(chrom, RC_Start, RC_End)):

		if ri % subsample_rate != 0:
			continue

		read_name = "{}:{}".format(read.query_name, read.is_read1)

		if get_read_strand(read, strandedness) != strand:
			continue

		if filter_reads(read, CUTOFFS["P_minMapQual"])[0]:
			continue
		
		blocks, soft_clippings = process_CIGAR(read.pos, read.cigar)

		if not any(b1 <= RC_End and RC_Start <= b2 for (b1, b2) in blocks):
			continue

		if any(SC > CUTOFFS["P_minSoftClipLen"] for SC in soft_clippings):
			continue

		aligned_seq, aligned_qual, new_blocks = remove_homopolymer_ends(read, blocks)

		for (b1, b2) in new_blocks:
			_read_blocks[b1] += 1
			_read_blocks[b2] -= 1

		SetOfReadNames[read_name]['Dir']    = read.is_reverse
		SetOfReadNames[read_name]['Quals']  = aligned_qual
		SetOfReadNames[read_name]['seq']    = aligned_seq
		SetOfReadNames[read_name]['blocks'] = new_blocks
		
		_unique_pos[(tuple(new_blocks), aligned_seq)][read_name] = sum(aligned_qual)


	bam_handle.close()

	RC_read_overlap = 0.0

	## read mono-clusters

	if all(len(v['blocks']) == 1 for (r, v) in SetOfReadNames.items()) and len(SetOfReadNames) < 500:
		read_pairwise_overlap = []
		for (r1, r2) in itertools.combinations(list(SetOfReadNames), 2):
			(r1_bs, r1_be) = SetOfReadNames[r1]['blocks'][0]
			(r2_bs, r2_be) = SetOfReadNames[r2]['blocks'][0]
			r1_len = r1_be - r1_bs
			r2_len = r2_be - r2_bs
			ov = 0.0	
			if r1_bs <= r2_be and r2_bs <= r1_be:
				ov = min([r2_len, r1_len, r1_be - r2_bs, r2_be - r1_bs])
			ov = ov/min(r1_len, r2_len)
			read_pairwise_overlap.append(ov)

		RC_read_overlap = np.median(read_pairwise_overlap) 


	gqv_reads_total = len(SetOfReadNames)
	gqv_reads_lost  = 0

	SetOfReadGroups = defaultdict(dict)
	list_index = 1

	for (unique_block, unique_seq), read_dict in _unique_pos.items():
		SetOfReadGroups[list_index]['seq']    = unique_seq
		SetOfReadGroups[list_index]['blocks'] = unique_block
		SetOfReadGroups[list_index]['reads']  = list(read_dict.keys())
		list_index += 1


	RC_cov_cutoff = 1
	Cum_Cov    = 0
	RC_started = False
	Covered_Ranges = []
	for coord in sorted(_read_blocks.keys()):
		Cum_Cov += _read_blocks[coord]

		if Cum_Cov >= RC_cov_cutoff and not RC_started:
			start_coord = coord
			RC_started = True
		elif Cum_Cov < RC_cov_cutoff and RC_started:
			if Covered_Ranges:
				if start_coord - Covered_Ranges[-1][1] > 20:
					Covered_Ranges.append([start_coord, coord])
				else:
					Covered_Ranges[-1][1] = coord
			else:
				Covered_Ranges.append([start_coord, coord])
			RC_started = False


	adapter1, adapter2 = "XX", "YY"

	Genome_reference_sequence = adapter1
	for (s, e) in Covered_Ranges:
		Genome_reference_sequence += str(genome.fetch(chrom, s + 1, e)).upper() + "I"
	Genome_reference_sequence = Genome_reference_sequence.rstrip("I")	
	Genome_reference_sequence += adapter2

	del _read_blocks, _unique_seqs

	
	if not Covered_Ranges or len(Genome_reference_sequence) < 50:
		return None
	else:
		RC_ext_Start = Covered_Ranges[0][0] 
		RC_ext_End   = Covered_Ranges[-1][1]
			
		REFPOS = genome_ref_pos(Covered_Ranges, chrom, RC_Start, RC_End, RC_index, 
					Genome_reference_sequence, genome, strand, RC_read_overlap)
					
		Ref_Blocks = [(RC_ext_Start - len(adapter1), RC_ext_End + len(adapter2))]

		SetOfReadNames[0]['Blocks'] = Ref_Blocks		
		SetOfReadGroups[0] = {'seq': Genome_reference_sequence, 'reads': [0], 'blocks': Ref_Blocks}

		return SetOfReadNames, REFPOS, SetOfReadGroups, gqv_reads_total, gqv_reads_lost



