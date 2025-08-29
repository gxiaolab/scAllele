#!/usr/bin/python
import pysam
from collections import defaultdict, Counter
from HTSeq import GenomicInterval
from time import time
import itertools
import random
import numpy as np ; np.seterr(all = "ignore")
import os 


def read_bam_input(option_bam):
	if "," in option_bam:
		bam_list = option_bam.split(',')
	elif option_bam.endswith('.bam'):
		bam_list = [option_bam]
	else:
		bam_list = []
		with open(option_bam) as f:
			for line in f:
				bam_file = line.strip() 
				bam_list.append(bam_file)

	for bam_file in bam_list:
		assert os.path.exists(bam_file) & os.path.exists(bam_file + '.bai'), f"invalid bam {bam_file}"

	return bam_list



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
	try:
		sample_name = bam_handle.header['RG'][0]['SM']
	except:
		sample_name = bam_file
	
	return sample_name



def find_super_read_clusters(bam_file_list, options, chrom = None, f_start = None, f_end = None):
	read_block_coords = defaultdict(dict)

	for bam_file in bam_file_list:

		bam_handle = pysam.Samfile(bam_file, 'rb')

		for read in bam_handle.fetch(chrom, start = f_start, end = f_end):
			chrom = read.reference_name

			if filter_reads(read, options.minMapQ)[0]:
				continue

			sB, eB = read.positions[0], read.positions[-1] 

			try:
				read_block_coords[chrom][sB] += 1
			except:
				read_block_coords[chrom][sB] = 1
			try:
				read_block_coords[chrom][eB] -= 1
			except:
				read_block_coords[chrom][eB] = -1

		bam_handle.close()


	Read_Clusters = []
	RC_cov_cutoff = options.minCoverage 

	for chrom, coords in read_block_coords.items():
		Cum_Cov_block = 0
		cluster_count = 0 

		RC_started = False

		coords = sorted(coords)

		for coord in coords: 
			Cum_Cov_block += read_block_coords[chrom][coord]
		
			if Cum_Cov_block >= RC_cov_cutoff and not RC_started:
				start_coord = coord
				RC_started  = True
			elif Cum_Cov_block < RC_cov_cutoff and RC_started:
				cluster_count += 1
				if cluster_count >= 1000 or coord == coords[-1]:
					Read_Clusters.append([chrom, start_coord, coord])
					RC_started = False
					cluster_count = 0 

	return Read_Clusters



def find_read_clusters(bam_file, options, chrom, start, end):
	
	read_block_coords = defaultdict(dict)

	bam_handle = pysam.Samfile(bam_file, 'rb')

	for read in bam_handle.fetch(chrom, start = start, end = end, until_eof = True):
		chrom = read.reference_name

		if filter_reads(read, options.minMapQ)[0]:
			continue

		r_blocks, soft_clippings = process_CIGAR(read.pos, read.cigar)

		if any(SC > options.maxSoftClip for SC in soft_clippings):
			continue

		r_strand = get_read_strand(read, options.strandedness)

		new_read = remove_homopolymer_ends(read, r_blocks)
		
		if new_read is None:
			continue

		aligned_seq, aligned_qual, r_blocks = new_read


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


	Read_Clusters = list()
	RC_cov_cutoff = options.minCoverage 
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
			
			if Max_Cov < options.minCoverage:
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


	
adapter1, adapter2 = "XX", "YY"
special_bases = set(['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'])

class genome_ref_pos():
	
	def __init__(self, CHROM, Covered_Ranges, RC_Range, RC_index, strand, genome):

		RC_ext_Start = Covered_Ranges[0][0] 
		RC_ext_End   = Covered_Ranges[-1][1]

		RC_ext_seq_compressed = adapter1
		for (s, e) in Covered_Ranges:
			
			block_seq = str(genome.fetch(CHROM, s + 1, e)).upper()
			for b in special_bases:
				block_seq = block_seq.replace(b, 'N')

			RC_ext_seq_compressed += block_seq  + "I"
		RC_ext_seq_compressed = RC_ext_seq_compressed.rstrip("I")	
		RC_ext_seq_compressed += adapter2

		RC_ext_seq_full = str(genome.fetch(CHROM, RC_ext_Start + 1, RC_ext_End)).upper()
		for b in special_bases:
				RC_ext_seq_full = RC_ext_seq_full.replace(b, 'N')
		RC_ext_seq_full = adapter1 + RC_ext_seq_full + adapter2


		self.Seq_ext      = RC_ext_seq_compressed
		self.Seq_ext_full = RC_ext_seq_full 

		RC_Start, RC_End  = RC_Range
 
		self.gi     = GenomicInterval(CHROM, RC_Start, RC_End)
		self.ext_gi = GenomicInterval(CHROM, RC_ext_Start, RC_ext_End)

		self.ranges = Covered_Ranges
		self.strand = strand 
		self.RC_i   = RC_index

	def rc_pos(self, genome_pos):
		RC_Pos = None
		CumLen = 2
		for (s, e) in self.ranges:
			if s <= genome_pos and genome_pos < e:
				RC_Pos = CumLen + genome_pos - s
			CumLen += (e - s)  + 1
		return RC_Pos

	def get_sequence(self, g_start, g_end):
		off = g_start - self.ext_gi.start + 2
		assert off >= 0 , f"VAR start {g_start} off {off} GI {self.gi} GI_ext {self.ext_gi}"
		return self.Seq_ext_full[off : off + (g_end - g_start)]

	def genome_pos(self, relative_pos):
		CumLen = len(adapter1)
		RealPos = None
		for (s, e) in self.ranges:
			ExonLen = e - s 
			if relative_pos > CumLen + ExonLen:
				CumLen += ExonLen + 1
			else:
				RealPos = s + relative_pos - CumLen
				break
		if RealPos is None:
			RealPos = self.ext_gi.end
		return RealPos


def remove_homopolymer_ends(read, blocks):
	READSEQ = read.query_alignment_sequence
	if "A"*50 in READSEQ or "C"*50 in READSEQ or "G"*50 in READSEQ or "T"*50 in READSEQ:
		return None
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



def check_read_pileup(SetOfReadNames):

	## read mono-clusters

	if any(len(v['blocks']) != 1 for v in SetOfReadNames.values()):
		return 0.0

	if len(SetOfReadNames) > 300:
		r_names = random.sample(list(SetOfReadNames), 300)
	else:
		r_names = list(SetOfReadNames)

	read_pairwise_overlap = []
	for (r1, r2) in itertools.combinations(r_names, 2):
		(r1_bs, r1_be) = SetOfReadNames[r1]['blocks'][0]
		(r2_bs, r2_be) = SetOfReadNames[r2]['blocks'][0]
		r1_len = r1_be - r1_bs
		r2_len = r2_be - r2_bs
		
		if r1_bs <= r2_be and r2_bs <= r1_be:
			ov = min([r2_len, r1_len, r1_be - r2_bs, r2_be - r1_bs])
		else:
			ov = 0.0
		ov = ov/min(r1_len, r2_len)
		read_pairwise_overlap.append(ov)

	RC_read_overlap = np.median(read_pairwise_overlap) 
	return RC_read_overlap 



def get_RC_reads(bam_file, RC_info, genome, SM, options):

	RC_index, chrom, strand, RC_Start, RC_End, MaxCov = RC_info

	SetOfReadNames = defaultdict(dict)	

	_read_blocks = defaultdict(int)
	_unique_seqs = defaultdict(dict)
	_unique_pos  = defaultdict(dict)


	subsample_rate = 2000/float(MaxCov) 

	bam_handle = pysam.Samfile(bam_file, 'rb')

	random.seed(0)
	
	for ri, read in enumerate(bam_handle.fetch(chrom, RC_Start, RC_End)):	

		if random.random() > subsample_rate:
			continue

		read_name = "{}:{}:{}".format(SM, read.query_name, read.is_read1)

		if get_read_strand(read, options.strandedness) != strand:
			continue

		if filter_reads(read, options.minMapQ)[0]:
			continue
		
		blocks, soft_clippings = process_CIGAR(read.pos, read.cigar)

		if not any(b1 <= RC_End and RC_Start <= b2 for (b1, b2) in blocks):
			continue

		if any(SC > options.maxSoftClip for SC in soft_clippings):
			continue

		new_read = remove_homopolymer_ends(read, blocks)
		if new_read is None:
			continue
		aligned_seq, aligned_qual, new_blocks = new_read
		for (b1, b2) in new_blocks:
			_read_blocks[b1] += 1
			_read_blocks[b2] -= 1

		SetOfReadNames[read_name]['Dir']    = read.is_reverse
		SetOfReadNames[read_name]['Quals']  = aligned_qual
		SetOfReadNames[read_name]['seq']    = aligned_seq
		SetOfReadNames[read_name]['blocks'] = new_blocks
		SetOfReadNames[read_name]['PcrDup'] = read.is_duplicate
		
		_unique_pos[(tuple(new_blocks), aligned_seq)][read_name] = sum(aligned_qual)


	bam_handle.close()

	RC_overlap_median = check_read_pileup(SetOfReadNames)

	if RC_overlap_median >= 0.95:
		return None


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

	del _read_blocks, _unique_seqs
	
	if not Covered_Ranges or sum((b - a) for a, b in Covered_Ranges) < options.kmer*2:
		return None
	else:
		REFPOS = genome_ref_pos(chrom, Covered_Ranges, (RC_Start, RC_End), RC_index, strand, genome)
		Ref_Blocks = [(REFPOS.ext_gi.start - len(adapter1), REFPOS.ext_gi.end + len(adapter2))]
		
		SetOfReadNames[0]['Blocks'] = Ref_Blocks		
		SetOfReadGroups[0] = {'seq': REFPOS.Seq_ext, 'reads': [0], 'blocks': Ref_Blocks}

		return SetOfReadNames, REFPOS, SetOfReadGroups



