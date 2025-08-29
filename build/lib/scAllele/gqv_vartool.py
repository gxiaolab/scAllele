#!/usr/bin/python
import re
import sys
import numpy as np ; np.seterr(all = "ignore")
import itertools 
import HTSeq
from time import time
from collections import defaultdict, abc 
# from Bio import Align.PairwiseAligner
from Bio import Align
from Bio.Align.substitution_matrices import Array

import functools

sys.setrecursionlimit(100000)

## General scoring

_Intron_Gap = 500

_Gap = 5

_min_complex_length = 4

## Short sequence scoring

_gap_open, _gap_ext = -0.6*_Gap, -0.4*_Gap  

bases = ["A", "C", "G", "T", "I", "N", "Y", "X"]

_base_pair_score = Array(''.join(bases), dims = 2)

for b1, b2 in itertools.product(bases, bases):
	if b1 == b2:
		_base_pair_score[b1, b2] = 0.0
	elif b1 == "I" or b2 == "I":
		_base_pair_score[b1, b2] = -1*_Intron_Gap
	elif b1 == "N" or b2 == "N":
		_base_pair_score[b1, b2] = -1
	else:
		_base_pair_score[b1, b2] = -1


aligner = Align.PairwiseAligner() 
aligner.mode = "global"
aligner.extend_gap_score = _gap_ext
aligner.open_gap_score   = _gap_open
aligner.substitution_matrix = _base_pair_score


class memoized(object):
	def __init__(self, function_name):
		self.function_name = function_name
		self.cache = {}

	def __call__(self, *args):
		if not isinstance(args, abc.Hashable):
			return self.function_name(*args)
		if args in self.cache:
			return self.cache[args]
		else:
			value = self.function_name(*args)
			self.cache[args] = value
			return value

	def __get__(self, obj, objtype):
		return functools.partial(self.__call__, obj)


class variant:
	def __init__(self, g_start, g_end, REF, ALT, r_pos):
		self.g_start = g_start
		self.g_end   = g_end
		self.REF = REF.upper()
		self.ALT = ALT.upper()
		self.RIPos = r_pos
		self.RGPos = {}
		self.match_intron = False
		self.rep_count  = 0
		self.rep_region = None
    
	def __hash__(self):
		return hash(str(self))

	def __eq__(self, other):
		return str(self) == str(other)

	def L_ref(self): 
		return len(self.REF)

	def L_alt(self): 
		return len(self.ALT)

	def __str__(self):
		if self.L_ref() > 20:
			REFSEQ = "{}__({})__{}".format(self.REF[:10], self.L_ref(), self.REF[-10:])
		else:
			REFSEQ = self.REF
		return "{}-{}:{}>{}".format(self.g_start, self.g_end, REFSEQ, self.ALT)

	def length(self):
		return self.L_alt() - self.L_ref()

	def shift_Rpos(self, offset):
		for ri, (ri_s, ri_e) in self.RIPos.items():
			self.RIPos[ri] = (ri_s + offset, ri_e + offset)

	def get_variant_type(self):
		if   self.L_ref() == 0 and self.L_alt() == 0:
			variant_type = "NO_VARIANT"
		elif self.L_ref() == 1 and self.L_alt() == 1:
			variant_type = "SNV"
		elif self.L_ref() == 0 and self.L_alt() > 0:
			variant_type = "INSERTION"
		elif self.L_ref() > 0 and self.L_alt() == 0:
			if self.match_intron:
				variant_type = "INTRON"
			else:
				variant_type = "DELETION"
		elif self.L_ref() * self.L_alt() > 1:
			if self.REF[0] == self.ALT[0]:
				if self.L_ref() == 1 and self.L_alt() > 1:
					variant_type = "INSERTION"
				elif self.L_ref() > 1 and self.L_alt() == 1:
					variant_type = "DELETION"
				else:
					variant_type = "COMPLEX"
			else:
				variant_type = "COMPLEX"
		return variant_type


	def trim(self, first = "sufix", min_match = 3):
		Pfx, Sfx, ALT, REF = trim_sequences(self.ALT, self.REF, first, min_match)
		self.REF = REF 
		self.ALT = ALT
		self.g_start = self.g_start + Pfx
		self.g_end   = self.g_end - Sfx
		for read_name, (r_start, r_end) in self.RIPos.items():
			self.RIPos[read_name] = (r_start + Pfx, r_end - Sfx)		


	def same_variant(self, var2, REFPOS):
		min_s = min(self.g_start, var2.g_start)
		max_e = max(self.g_end, var2.g_end)
		hap1 = REFPOS.get_sequence(min_s, self.g_start) + self.ALT + REFPOS.get_sequence(self.g_end, max_e)
		hap2 = REFPOS.get_sequence(min_s, var2.g_start) + var2.ALT + REFPOS.get_sequence(var2.g_end, max_e)
		return hap1 == hap2


	def leftalign(self, REFPOS, max_shift = 100):
		max_shift = max(REFPOS.ext_gi.start + 1, self.g_start - max_shift)
		seq = REFPOS.get_sequence(max_shift, self.g_start + len(self.REF))
		assert seq.endswith(self.REF)

		if self.L_ref():
			seq = seq[:-self.L_ref()]
		
		ref, alt = self.REF, self.ALT 
		offset = 0
		quit = False
		while offset < len(seq) and not quit:
			quit = True
			if (ref and alt) and ref[-1] == alt[-1]:
				ref, alt = ref[:-1], alt[:-1]
				quit = False
			if len(ref) == 0 or len(alt) == 0:
				offset += 1
				ref = seq[-offset] + ref
				alt = seq[-offset] + alt
				quit = False

		self.REF = ref.upper()
		self.ALT = alt.upper() 
		self.g_start = self.g_start - offset
		self.g_end   = self.g_end - offset + 1
		for RI, (r_start, r_end) in self.RIPos.items():
			self.RIPos[RI] = (r_start - offset, r_end - offset + 1)


@memoized
def trim_sequences(rid_seq, ref_seq, first = 'prefix', min_match = 3):
	'''	
	Trim common sequences and label the variant.
	'''
	def _find_prefix(rid_seq, ref_seq):
		minLength = min(len(rid_seq), len(ref_seq))
		prefix = 0
		while prefix + 1 <= minLength and rid_seq[:prefix + 1] == ref_seq[:prefix + 1]:
			prefix += 1
		return prefix
	
	def _find_sufix(rid_seq, ref_seq):
		minLength = min(len(rid_seq), len(ref_seq))
		sufix = 0
		while sufix + 1 <= minLength and rid_seq[-(sufix + 1):] == ref_seq[-(sufix + 1):]:
			sufix += 1   
		return sufix

	sufix_final  = 0
	prefix_final = 0

	if first == 'prefix':

		prefix = _find_prefix(rid_seq, ref_seq)

		if prefix >= min_match or prefix == min(len(rid_seq), len(ref_seq)):
			rid_seq = rid_seq[prefix:]
			ref_seq = ref_seq[prefix:]
			prefix_final = prefix
		
		sufix = _find_sufix(rid_seq, ref_seq)		

		if sufix >= min_match or sufix == min(len(rid_seq), len(ref_seq)):
			rid_seq = rid_seq[:(-sufix or None)]
			ref_seq = ref_seq[:(-sufix or None)]
			sufix_final = sufix
	else:

		sufix = _find_sufix(rid_seq, ref_seq)

		if sufix >= min_match or sufix == min(len(rid_seq), len(ref_seq)):
			rid_seq = rid_seq[:(-sufix or None)]
			ref_seq = ref_seq[:(-sufix or None)]
			sufix_final = sufix

		prefix = _find_prefix(rid_seq, ref_seq)
		
		if prefix >= min_match or prefix == min(len(rid_seq), len(ref_seq)):
			rid_seq = rid_seq[prefix:]
			ref_seq = ref_seq[prefix:]
			prefix_final = prefix

	return prefix_final, sufix_final, rid_seq, ref_seq



def convert_to_var(REF, ALT, Ref_S, Ref_E, Rid_S, Rid_E, RI, REFPOS):

	Ref_S = REFPOS.genome_pos(Ref_S)
	Ref_E = REFPOS.genome_pos(Ref_E)
	if "I" in REF:
		REF = REFPOS.get_sequence(Ref_S, Ref_E)
	var = variant(Ref_S, Ref_E, REF, ALT, {RI : (Rid_S, Rid_E)})
	return var



def calculate_edit_distance(var, Junctions, REFPOS):

	is_intron = False 
	for (jxn_s, jxn_e) in sorted(Junctions):
		if jxn_s <= var.g_end and var.g_start <= jxn_e:
			jxn_var = variant(jxn_s, jxn_e, '', '', {})
			is_intron = var.same_variant(jxn_var, REFPOS)

	var.match_intron = is_intron

	if var.get_variant_type() == "INTRON":
		editd = _Intron_Gap
		var_set = [var]

	elif var.get_variant_type() == "NO_VARIANT":
		editd = 0.0 
		var_set = []

	elif var.get_variant_type() == "DELETION":
		editd = var.L_ref()*_Gap
		var_set = [var]

	elif var.get_variant_type() == "INSERTION":
		editd = var.L_alt()*_Gap
		var_set = [var]

	elif var.get_variant_type() == "SNV":
		editd = abs(_base_pair_score[(var.ALT, var.REF)])
		var_set = [var]

	elif var.get_variant_type() == "COMPLEX":
		editd, var_set = split_complex_var(var, Junctions, REFPOS)

	editd_wo_intron = editd
	n_vars = 0

	for v in var_set:
		if v.match_intron:
			editd_wo_intron -= _Intron_Gap
		if v.get_variant_type() in ('INSERTION', 'DELETION'):
			n_vars += 1
		elif v.get_variant_type() == 'COMPLEX':
			n_vars += 2

	return editd, editd_wo_intron, var_set, n_vars




def split_complex_var(var, Junctions, REFPOS):

	Matrix_Product_Cutoff_1 = 2000
	LenSim = np.log2(var.L_ref()/var.L_alt())

	if var.L_alt() <= _min_complex_length and var.L_ref() <= _min_complex_length:
		complex_var_type = "BOTH_SHORT"
		if var.L_ref() == var.L_alt():
			var_set, editd = edit_distance_function_msnp(var)
		else:
			editd = edit_distance_function_WA(var)
			var_set  = [var]

	elif var.L_alt() * var.L_ref() < Matrix_Product_Cutoff_1 and abs(LenSim) < 0.6:
		complex_var_type = "BOTH_MEDIUM"
		var_set, editd = edit_distance_function_biopython(var, Junctions, REFPOS)

	else:
		complex_var_type = "LONG"
		var_set, editd = edit_distance_function_internal_aln(var, Junctions, REFPOS)

	return editd, var_set



def edit_distance_lazy_calc(var):
	delta  = abs(var.L_ref() - var.L_alt())*_Gap
	common = min(var.L_ref(), var.L_alt())
	return common + delta



def edit_distance_function_msnp(var):
	editd = 0.0
	var_set = []
	for i in range(var.L_ref()):
		snp_ReadPos = {}
		for r_name, (r_start, r_end) in var.RIPos.items():
			snp_ReadPos[r_name] = (r_start + i, r_start + i + 1)
		snp_var = variant(var.g_start + i, var.g_start + i + 1, var.REF[i], var.ALT[i], snp_ReadPos)
		editd  += abs(_base_pair_score[(var.REF[i], var.ALT[i])])
		var_set.append(snp_var)
	return var_set, editd



def edit_distance_function_biopython(var, Junctions, REFPOS):

	alignments = aligner.align(var.ALT, var.REF) 

	# ALT_aln, REF_aln, Edit_Distance, Start, End = alignments[0]
	ALT_aln, REF_aln = alignments[0] 
	assert len(ALT_aln) == len(REF_aln)
	
	# --TGT-TTGT 
	# GGTGTG--G-


	total_editd = 0
	coord_list = []

	ref_indexes, rid_indexes = [], []
	ref_i, rid_i = var.g_start, 0

	for (b_f, b_r) in zip(REF_aln + 'X', ALT_aln + 'X'):
		ref_indexes.append(ref_i)
		rid_indexes.append(rid_i)
		if b_f != "-" and b_r != "-":
			ref_i += 1
			rid_i += 1
		elif b_f == "-":
			rid_i += 1
		elif b_r == "-":
			ref_i += 1

	IN_started, DL_started = False, False

	for i, (b_f, b_r) in enumerate(zip(REF_aln, ALT_aln)):

		if b_f == "-" and b_r == "-":
			raise ValueError("Weird alignment {} {}\n".format(REF_aln, ALT_aln))

		if b_f == "-" and not IN_started:
			var_ins_start = i
			IN_started = True
		elif b_f != "-" and IN_started:
			coord_list.append((var_ins_start, i))
			IN_started = False

		if b_r == "-" and not DL_started:
			var_del_start = i
			DL_started = True
		elif b_r != "-" and DL_started:
			coord_list.append((var_del_start, i))
			DL_started = False

		if b_r != b_f and b_f != "-" and b_r != "-":
			coord_list.append((i, i + 1))
			 
	if IN_started:
		coord_list.append((var_ins_start, i + 1))
	if DL_started:
		coord_list.append((var_del_start, i + 1))

	final_var_set = []
	total_n_vars = 0
	for (s, e) in coord_list:
		tmp_Rid_S, tmp_Rid_E = rid_indexes[s], rid_indexes[e]
		tmp_Ref_S, tmp_Ref_E = ref_indexes[s], ref_indexes[e]
		ALT = ALT_aln[s : e].strip("-")
		REF = REF_aln[s : e].strip("-")
		_ReadPos = {}
		for r_name, (r_start, r_end) in var.RIPos.items():
			_ReadPos[r_name] = (r_start + tmp_Rid_S, r_start + tmp_Rid_E)
		tmpvar = variant(tmp_Ref_S, tmp_Ref_E, REF, ALT, _ReadPos)

		editd, editd_wo_intron, var_set, n_vars = calculate_edit_distance(tmpvar, Junctions, REFPOS)
		total_editd += editd
		total_n_vars += n_vars 
		final_var_set.extend(var_set)

	return final_var_set, total_editd



def edit_distance_function_WA(var):

	mat = np.zeros((var.L_ref() + 1, var.L_alt() + 1)).astype(int)
	for i, char in enumerate(var.REF):
		mat[i + 1, 0] = (i + 1)*_Gap
	for j, char in enumerate(var.ALT):
		mat[0, j + 1] = (j + 1)*_Gap

	for i, char_ref in enumerate(var.REF):
		for j, char_rid in enumerate(var.ALT):
			subs = abs(_base_pair_score[(char_rid, char_ref)])
			mat[i + 1, j + 1] = min(mat[i, j] + subs, 
									mat[i, j + 1] + _Gap, 
									mat[i + 1, j] + _Gap)

	return mat[-1][-1]




def edit_distance_function_internal_aln(var, Junctions, REFPOS):

	def decompose_seq(seq, K, _junctions):

		intronic_seqs = {}
		for jxn_s, jxn_e in sorted(list(_junctions)):
			if var.g_start <= jxn_e and var.g_end >= jxn_s:
				pos = jxn_s - var.g_start
				intronic_seqs[jxn_s - var.g_start] = var.REF[pos : pos + (jxn_e - jxn_s)]
		ov_nodes = []

		j = 0
		while j < len(seq) - (K - 2):
			if j in intronic_seqs:
				intronic_seq = intronic_seqs[j]
				ov_nodes.append(intronic_seq)
				j = j + len(intronic_seq) - (K - 2)
			else:
				node = seq[j : j + (K - 1)]
				ov_nodes.append(node)
				j += 1
			if len(ov_nodes) > int(1e3):
				j = len(seq)
				ov_nodes = []

		return ov_nodes

	ALT_path = decompose_seq(var.ALT, _min_complex_length, [])
	REF_path = decompose_seq(var.REF, _min_complex_length, Junctions) 

	#### Too complex
	if len(REF_path)*len(ALT_path) > int(1e4) or len(REF_path)*len(ALT_path) == 0:
		return [var], edit_distance_lazy_calc(var)


	mat = np.zeros((len(ALT_path), len(REF_path))).astype(int)
	max_scores = {}

	for i in range(mat.shape[0]):
		cum_j = 0
		for j in range(mat.shape[1]):
			if ALT_path[i] == REF_path[j]:
				if i > 0 and j > 0:
					mat[i][j] = 1 + mat[i - 1][j - 1] 
				else:
					mat[i][j] = 1 

				mat_score = mat[i][j]
				
				if mat_score >= 2:
					origin = (i - mat_score + 1, cum_j - mat_score + 1)
					if origin not in max_scores or mat_score > max_scores[origin][0]:
						max_scores[origin] = (mat_score, i, cum_j)

			cum_j += len(REF_path[j]) - _min_complex_length + 2 


	del mat

	all_var_sets = []

	top = 1
	while max_scores and top <= 3:
		top += 1
		
		origin_max_score = max(max_scores, key = max_scores.get)

		RidS, RefS = origin_max_score

		(max_score, RidE, RefE) = max_scores[origin_max_score]
		max_scores.pop(origin_max_score)

		var_set = []

		### upstream of the aligned region
		Pfx, Sfx, S2, S1 = trim_sequences(var.ALT[:RidS], var.REF[:RefS])
		assert len(S1) + Pfx + Sfx == len(var.REF[:RefS]), "PP1 {} {} {}".format(S1, Pfx, Sfx, var.REF[:RefS])
		assert len(S2) + Pfx + Sfx == len(var.ALT[:RidS]), "PP2 {} {} {}".format(S2, Pfx, Sfx, var.ALT[:RidS])

		P1s, P1e = var.g_start + Pfx, var.g_start + RefS - Sfx
		P2s, P2e = Pfx, RidS - Sfx

		assert (P1e - P1s) == len(S1)
		assert (P2e - P2s) == len(S2)

		_ReadPos = {}
		for RI, (s, e) in var.RIPos.items():
			_ReadPos[RI] = (s + P2s, s + P2s + len(S2))
		_tmp_var = variant(P1s, P1e, S1, S2, _ReadPos)
		var_set.append(_tmp_var)

		
		### aligned region
		S1, S2 = var.REF[RefS : RefE], var.ALT[RidS : RidE]
		assert S1 == S2, "{} != {}".format(S1, S2)
		
		### downstream of the aligned region
		Pfx, Sfx, S2, S1 = trim_sequences(var.ALT[RidE:], var.REF[RefE:])
		assert len(S1) + Pfx + Sfx == len(var.REF[RefE:]), "PP3 {} {} {}".format(S1, Pfx, Sfx, var.REF[RefE:])
		assert len(S2) + Pfx + Sfx == len(var.ALT[RidE:]), "PP4 {} {} {}".format(S2, Pfx, Sfx, var.ALT[RidE:])

		P1s, P1e = var.g_start + RefE + Pfx, var.g_end - Sfx
		P2s, P2e = RidE + Pfx, len(var.ALT) - Sfx

		assert (P1e - P1s) == len(S1)
		assert (P2e - P2s) == len(S2)
		
		_ReadPos = {}
		for RI, (s, e) in var.RIPos.items():
			_ReadPos[RI] = (s + P2s, s + P2s + len(S2))
		_tmp_var = variant(P1s, P1e, S1, S2, _ReadPos)
		var_set.append(_tmp_var)
		all_var_sets.append(var_set)

	
	all_var_editd = {}
	for j, var_set in enumerate(all_var_sets):
		var_set_editd = 0
		var_set_editd_wo_intron = 0
		var_set_copy  = []
		
		for _tmp_var in var_set:
			EDITVAL = calculate_edit_distance(_tmp_var, Junctions, REFPOS)
			var_set_editd += EDITVAL[0]
			var_set_editd_wo_intron += EDITVAL[1]
			var_set_copy.extend(EDITVAL[2])

		all_var_sets[j]  = var_set_copy
		all_var_editd[j] = (var_set_editd_wo_intron, var_set_editd)

	try:
		min_j = min(all_var_editd, key = all_var_editd.get)
		Final_Var_Set   = all_var_sets[min_j]
		Final_Var_EditD = all_var_editd[min_j][1]
		
	except:
		Final_Var_Set   = [var]
		Final_Var_EditD = edit_distance_lazy_calc(var)

	return Final_Var_Set, Final_Var_EditD



def merge_introns(VAR_LIST, SM):
    
    intron_coords = HTSeq.GenomicArrayOfSets(chroms = "auto", stranded = False)

    sorted_introns = []
    
    for POS, ALLELE_dict in VAR_LIST.items():
        if 'INTRON' in POS:
            for ALLELE in ALLELE_dict:
            	if SM in ALLELE_dict[ALLELE]:
	                sorted_introns.append((POS, ALLELE))

    sorted_introns.sort(key = lambda x : x[0][1]) 

    sorted_introns_norep = []

    for intron in sorted_introns:
        POS, ALLELE = intron
        (CHROM, START, var_group), (REF, ALT, END) = POS, ALLELE

        if END - START < 20:
            continue  	  

        similar_found = False
        minD = 6 + len(ALT)

        for j in range(1, 6):
            if j > len(sorted_introns_norep): 
                break

            prev_POS, prev_ALLELE = sorted_introns_norep[-j]

            pSTART = prev_POS[1]
            pEND   = prev_ALLELE[2]

            if abs(START - pSTART) <= minD and abs(END - pEND) <= minD:
                similar_found = True
                VAR_LIST[prev_POS][prev_ALLELE][SM]["READS"].update(VAR_LIST[POS][ALLELE][SM]["READS"])
                iGI = HTSeq.GenomicInterval(CHROM, pSTART, pEND)
                intron = sorted_introns_norep[-j]
                break
        
        if not similar_found:
            sorted_introns_norep.append(intron)
            iGI = HTSeq.GenomicInterval(CHROM, START, END)
        
        intron_coords[iGI] += intron
        

    for iGI, intron_set in intron_coords.steps():

        if iGI.length > 20 and len(intron_set) >= 2:
            ip_id = (iGI.chrom, iGI.start, "INTRONIC_PART")
            ip_allele = ('', '', iGI.end)

            VAR_LIST[ip_id][ip_allele] = {SM : {"READS" : {}, "FEAT" : {}, "IP_ALLELES": [], "IP_ALLELES_AC" : []} }
            VAR_LIST[ip_id][ip_allele][SM]['FEAT']['FILTER'] = "PASS"
            VAR_LIST[ip_id][ip_allele][SM]['FEAT']['QUAL']   = 10.0
            VAR_LIST[ip_id][ip_allele]['variant_type'] = "INTRONIC_PART"
            VAR_LIST[ip_id][ip_allele]['variant_id']   = str(iGI)
            
            intron_set = list(intron_set)

            for i, intron in enumerate(intron_set):
                i_reads = VAR_LIST[intron[0]][intron[1]][SM]["READS"]

                VAR_LIST[ip_id][ip_allele][SM]["IP_ALLELES"].append(intron) 

                AC = 0
                for read_name, read_allele in i_reads.items():
                    if read_allele:
                        VAR_LIST[ip_id][ip_allele][SM]["READS"][read_name] = i
                        AC += 1

                VAR_LIST[ip_id][ip_allele][SM]["IP_ALLELES_AC"].append(AC)


    return VAR_LIST



def count_tandem_repeat_around_variant(var, REFPOS):

	tmp_vals = []
	for repLen in range(1, 5):
		for (l1, l2) in [(var.g_start, var.g_start + repLen), (var.g_start - repLen, var.g_start)]:

			repSeq = REFPOS.get_sequence(l1, l2)
			
			repCount_up = 0
			start_pos = l1
			while start_pos - repLen > REFPOS.gi.start + 2 and REFPOS.get_sequence(start_pos - repLen, start_pos) == repSeq:
				repCount_up += 1
				start_pos -= repLen

			repCount_dn = 0
			end_pos = l1
			while end_pos + repLen < REFPOS.gi.end - 2 and REFPOS.get_sequence(end_pos, end_pos + repLen) == repSeq:
				repCount_dn += 1
				end_pos += repLen

			repCount = repCount_up + repCount_dn
			tmp_vals.append((start_pos, end_pos, repSeq, repCount))

	(repeat_pos, repeat_end, repSeq, REP_COUNT) = max(tmp_vals, key = lambda x : x[3])

	REP_REGION = HTSeq.GenomicInterval(REFPOS.gi.chrom, repeat_pos - 1, repeat_end + 1)

	return REP_COUNT, REP_REGION



def trim_var_with_var(var1, var2):
	if var1 is None or var2 is None:
		return None
	elif var2.REF not in var1.REF or var2.ALT not in var1.ALT:
		return None
	elif var2.REF == var1.REF and var2.ALT == var1.ALT:
		return None
	elif var2.L_ref() > var1.L_ref() or var2.L_alt() > var1.L_alt():
		return None
	elif var2.get_variant_type() == "NO_VARIANT":
		return None

	if var2.g_start == var1.g_start:
		Pfx_ref, Sfx_ref, REF_12, ref_2c = trim_sequences(var1.REF, var2.REF, 'prefix')
		Pfx_alt, Sfx_alt, ALT_12, alt_2c = trim_sequences(var1.ALT, var2.ALT, 'prefix')

	elif var2.g_end == var1.g_end:
		Pfx_ref, Sfx_ref, REF_12, ref_2c = trim_sequences(var1.REF, var2.REF, 'sufix')
		Pfx_alt, Sfx_alt, ALT_12, alt_2c = trim_sequences(var1.ALT, var2.ALT, 'sufix')
	else:
		return None

	if ref_2c or alt_2c:
		return None

	_ReadPos = {}
	for r, (r_start, r_end) in var1.RIPos.items():
		_ReadPos[r] = (r_start + Pfx_alt, r_end - Sfx_alt)

	variant_3 = variant(var1.g_start + Pfx_ref, var1.g_end - Sfx_ref, REF_12, ALT_12, _ReadPos)
	variant_3.trim(min_match = 1)

	return variant_3




def find_equivalent_indels(VAR_COLLECTION, REFPOS):

	### First trim
	for i, variant in enumerate(VAR_COLLECTION):
		VAR_COLLECTION[i].trim(min_match = 1)
		if variant.get_variant_type() == "NO_VARIANT":
			VAR_COLLECTION[i] = None
		

	i = 0
	while i < len(VAR_COLLECTION):

		variant_1 = VAR_COLLECTION[i]
		
		if variant_1 is None or variant_1.get_variant_type() != "COMPLEX":
			i += 1
			continue
		
		split_candidates = {}
		for variant_2 in VAR_COLLECTION:
			variant_3 = trim_var_with_var(variant_1, variant_2)
			
			if variant_3 is None:
				continue
			
			e2 = calculate_edit_distance(variant_2, [], REFPOS)
			e3 = calculate_edit_distance(variant_3, [], REFPOS)
			split_candidates[(variant_2, variant_3)] = e2[0] + e3[0]
		

		if split_candidates:

			(variant_2, variant_3) = min(split_candidates, key = split_candidates.get)
			if split_candidates[(variant_2, variant_3)]	< 100.1:
				VAR_COLLECTION.append(variant_3)
				_ReadPos = {}
				Pfx = variant_2.g_start - variant_1.g_start
				Sfx = variant_2.g_end - variant_1.g_end
				for RI, (r_start, r_end) in variant_1.RIPos.items():
					_ReadPos[RI] = (r_start + Pfx, r_end - Sfx)
				variant_2.RIPos.update(_ReadPos)
				VAR_COLLECTION[i] = None
		i += 1
	

	### Find repeats around variants and move it to consensus location (beggining of rep)
	for i, variant in enumerate(VAR_COLLECTION):
		if variant is None:
			continue
		elif variant.get_variant_type() == "INTRON":
			variant.rep_region = HTSeq.GenomicInterval(REFPOS.gi.chrom, variant.g_start - 1, variant.g_end + 1)
			continue
		elif variant.get_variant_type() == "SNV":
			REP_COUNT, REP_REGION = count_tandem_repeat_around_variant(variant, REFPOS)
			variant.rep_region = HTSeq.GenomicInterval(REFPOS.gi.chrom, variant.g_start, variant.g_end)
			variant.rep_count  = REP_COUNT
			continue
		else:
			REP_COUNT, REP_REGION = count_tandem_repeat_around_variant(variant, REFPOS)
			variant.rep_count  = REP_COUNT
			variant.rep_region = REP_REGION

	### Find equivalent indels (same haplotype)
	for i, variant_1 in enumerate(VAR_COLLECTION):
		if variant_1 is None:
			continue
		if variant_1.get_variant_type() in ("INSERTION", "DELETION"):
			variant_1.leftalign(REFPOS)
		for j, variant_2 in enumerate(VAR_COLLECTION):
			if i == j or variant_2 is None:
				continue
			if variant_1.same_variant(variant_2, REFPOS):
				if variant_2.get_variant_type() in ("INSERTION", "DELETION"):
					variant_2.leftalign(REFPOS)
				VAR_COLLECTION[i].RIPos.update(variant_2.RIPos)
				VAR_COLLECTION[j] = None

	### Revome variants outside of RC coordinates
	for i, variant in enumerate(VAR_COLLECTION):
		if variant is None:
			continue
		elif variant.get_variant_type() == "INTRON":
			continue
		else:
			if REFPOS.gi.start <= variant.g_start and variant.g_start <= REFPOS.gi.end:
				pass
			else:
				VAR_COLLECTION[i] = None

	### Final trimming
	for i, variant in enumerate(VAR_COLLECTION):
		if variant is None:
			continue
		elif variant.get_variant_type() == "NO_VARIANT":
			VAR_COLLECTION[i] = None
			continue
		VAR_COLLECTION[i].trim(min_match = 1)
		if variant.get_variant_type() in ("INSERTION", "DELETION"):
			VAR_COLLECTION[i].leftalign(REFPOS)
	
	return [v for v in VAR_COLLECTION if v is not None]

