#!/usr/bin/python
import os
import sys
import copy
import numpy as np ; np.seterr(all = "ignore")
import _pickle as pickle
import itertools
import datetime
import scipy.stats as stats
from collections import defaultdict, Counter
from time import time
import pyfaidx 
import HTSeq

from . import gqv_software_management as SOFTMAN
from . import gqv_dbg_utils
from . import gqv_vartool

import warnings

warnings.filterwarnings('ignore')


def assign_vars_to_reads(src_snk_pairs, SetOfReadIndexes, SetOfReadGroups, options, REFPOS):

    REF_INDEX = SetOfReadGroups[0]['Index']
    k  = options.kmer
    OV = k - 2

    SetOfVariants = list()
    memoization_vars = defaultdict(set)

    r_total, r_unmapped = 0, 0

    for Read_Index, bubble_list_tmp in src_snk_pairs.items():
        _match_src_snk_to_reads(Read_Index, REF_INDEX, bubble_list_tmp, SetOfReadIndexes, 
                                options, memoization_vars, REFPOS, ronda = 1)

    for Read_Index in range(len(SetOfReadIndexes)):

        try:
            bubble_list_tmp = src_snk_pairs.pop(Read_Index)
        except KeyError:
            bubble_list_tmp = {}

        _match_src_snk_to_reads(Read_Index, REF_INDEX, bubble_list_tmp, SetOfReadIndexes, 
                                options, memoization_vars, REFPOS, ronda = 2)
        
        Graph, source, target = gqv_dbg_utils.make_read_dbg(Read_Index, REF_INDEX,  
                                                            SetOfReadIndexes)

        # gqv_dbg_utils.draw_rid_graph(Graph, f'rid_graph.{Read_Index}.txt')

        SetOfReadIndexes[Read_Index].delete_events()

        edges = defaultdict(list)
        
        all_paths = gqv_dbg_utils.Dijkstra_find_end_to_end_paths(Read_Index, REF_INDEX, 
                                                                SetOfReadIndexes, Graph, 
                                                                source, target, edges, REFPOS)
        # Killer switch 
        if all_paths is None:
            SOFTMAN.print_time_stamp("WARNING: skipping RC {}. Too many recursions".format(REFPOS.gi))
            break

        RI_alignment, var_list = _find_best_read_path(Read_Index, REF_INDEX, all_paths, 
                                                        REFPOS, k, SetOfReadIndexes, 
                                                        SetOfReadGroups)
        r_total += SetOfReadIndexes[Read_Index].n

        if not RI_alignment and REF_INDEX != Read_Index:
            r_unmapped += SetOfReadIndexes[Read_Index].n 
        
        SetOfVariants.extend(var_list)        
        SetOfReadIndexes[Read_Index].add_mapping(RI_alignment)
        
    # print("Unmapped_rate {} {:.3f} n= {} {}".format(REFPOS.gi, r_unmapped/r_total, r_unmapped, r_total))
    return SetOfVariants



def _src_snk_comb(src_list, snk_list, src_seq, snk_seq, k):
    '''
    Src - Snk pairs in homopolymers/tandem repeats. 
    Only keep outer nodes
    '''

    eff_src_list = src_list[:]
    eff_snk_list = snk_list[:]

    if len(src_list) > 5:
        for i in range(len(src_list) - 1):
            if src_list[i] - src_list[i - 1] == 1 and src_list[i + 1] - src_list[i] == 1:
                eff_src_list[i] = 'NA'
        eff_src_list = [x for x in eff_src_list if x != 'NA']

    if len(snk_list) > 5:
        for i in range(1, len(snk_list) - 1):
            if snk_list[i] - snk_list[i - 1] == 1 and snk_list[i + 1] - snk_list[i] == 1:
                eff_snk_list[i] = 'NA'
        eff_snk_list = [x for x in eff_snk_list if x != 'NA']


    if 'N' in src_seq or 'N' in snk_seq:
        combs = []

    elif len(eff_snk_list) > 5 and len(eff_src_list) > 5:
        space1 = [eff_snk_list[i + 1] - eff_snk_list[i] for i in range(len(eff_snk_list) - 1)]
        space2 = [eff_src_list[i + 1] - eff_src_list[i] for i in range(len(eff_src_list) - 1)]
        
        if abs(np.mean(space1) - np.mean(space2)) < 1:
            combs = []
            for i in range(len(eff_snk_list[:-1])):
                snk_i = eff_snk_list[i + 1]
                if snk_i >= eff_src_list[0]:
                    combs.append((eff_src_list[0], snk_i))
        else:
            combs = []
            for s_i, s_j in itertools.product(eff_src_list, eff_snk_list):
                delta = s_j - s_i
                if delta >= 0 and delta < 2*k:
                    combs.append((s_i, s_j))
    else:
        combs = []
        for s_i, s_j in itertools.product(eff_src_list, eff_snk_list):
            if s_j >= s_i:
                combs.append((s_i, s_j))

    return combs



def _match_src_snk_to_reads(Read_Index, REF_Index, bubble_list_tmp, SetOfReadIndexes, 
        options, memoization_vars, REFPOS, ronda):

    k = options.kmer
    OV = k - 2
    Ref_Seq_path, Ref_node_index = SetOfReadIndexes[REF_Index].PathSeq(return_index=True)
    Ref_Seq_index = SetOfReadIndexes[REF_Index].PathIndex
    Ref_Long_Seq  = SetOfReadIndexes[REF_Index].LongSeq
    
    Rid_Seq_path, Rid_node_index = SetOfReadIndexes[Read_Index].PathSeq(return_index=True)
    Rid_Seq_index = SetOfReadIndexes[Read_Index].PathIndex
    Rid_Long_Seq  = SetOfReadIndexes[Read_Index].LongSeq
    
    memoization_read_haplotypes = defaultdict(list)


    if Rid_Long_Seq in Ref_Long_Seq:
        bubble_list_tmp = {}
    
    for (src_seq, snk_seq), SrcSnk_Path_Seqs in bubble_list_tmp.items():
        assert len(src_seq) == k - 1 and len(snk_seq) == k - 1
        
        ref_src_i = Ref_node_index[src_seq]
        ref_snk_i = Ref_node_index[snk_seq]
        rid_src_i = Rid_node_index[src_seq]
        rid_snk_i = Rid_node_index[snk_seq]
        ref_combs = _src_snk_comb(ref_src_i, ref_snk_i, src_seq, snk_seq, k)
        rid_combs = _src_snk_comb(rid_src_i, rid_snk_i, src_seq, snk_seq, k)

        for ref_i, ref_j in ref_combs:
            f1 = Ref_Seq_index[ref_i] + OV
            f2 = Ref_Seq_index[ref_j] + OV
            SrcSnk_Ref_Seq = Ref_Long_Seq[f1 : f2]
            comb_found = False
            for rid_i, rid_j in rid_combs:
                r1 = Rid_Seq_index[rid_i] + OV
                r2 = Rid_Seq_index[rid_j] + OV
                SrcSnk_Rid_Seq = Rid_Long_Seq[r1 : r2]

                if not SrcSnk_Rid_Seq in SrcSnk_Path_Seqs:
                    continue

                comb_found = True

                Ref_Long_Seq_Var = Ref_Long_Seq[:f1] + SrcSnk_Rid_Seq + Ref_Long_Seq[f2:]
                Rid_Long_Seq_Var = Rid_Long_Seq[:r1] + SrcSnk_Ref_Seq + Rid_Long_Seq[r2:]

                if Rid_Long_Seq_Var in memoization_read_haplotypes[Ref_Long_Seq_Var]: 
                    continue
                else:
                    memoization_read_haplotypes[Ref_Long_Seq_Var].append(Rid_Long_Seq_Var)

                # memoized via decorator 
                V_t = gqv_vartool.trim_sequences(SrcSnk_Rid_Seq, SrcSnk_Ref_Seq, 'sufix')

                Pfx, Sfx, d_Seq, f_Seq = V_t 
                # print(f"SRCSNK RI={Read_Index}", V_t)
                # print('srcsnk', f1 + Pfx, REFPOS.genome_pos(f1 + Pfx), len(d_Seq) - len(f_Seq), f_Seq, d_Seq)

                if len(d_Seq) - len(f_Seq) > 20:
                    continue

                memoization_vars[src_seq, snk_seq, SrcSnk_Ref_Seq, f1, f2].add(V_t)
                if ronda == 2:
                    RdP = gqv_dbg_utils.seq_interval_obj('RID', r1 + Pfx, r2 - Sfx)
                    RfP = gqv_dbg_utils.seq_interval_obj('REF', f1 + Pfx, f2 - Sfx)
                    SetOfReadIndexes[Read_Index].add_event(RfP, RdP, f_Seq, d_Seq)

            if ronda == 1:
                for SrcSnk_Path_Seq in SrcSnk_Path_Seqs:
                    # memoized via decorator
                    V_t = gqv_vartool.trim_sequences(SrcSnk_Path_Seq, SrcSnk_Ref_Seq, 'sufix')
                    memoization_vars[src_seq, snk_seq, SrcSnk_Ref_Seq, f1, f2].add(V_t)

    if ronda == 2:
        for memo_key in memoization_vars:
            	
            src_seq, snk_seq, SrcSnk_Ref_Seq, f1, f2 = memo_key
            ref_src_i = Ref_node_index[src_seq]
            ref_snk_i = Ref_node_index[snk_seq]
            rid_src_i = Rid_node_index[src_seq]
            rid_snk_i = Rid_node_index[snk_seq]

            vars_from_comb = memoization_vars[memo_key] 

            for rid_i in rid_src_i:
                if rid_i >= len(Rid_Seq_index) - 1:
                    continue
                r1_tmp = Rid_Seq_index[rid_i] + OV
                for Pfx, Sfx, d_Seq, f_Seq in vars_from_comb:

                    r1 = r1_tmp + Pfx
                    r2 = r1 + len(d_Seq)
                    # print(f"SRCSNK RI={Read_Index}", (Pfx, Sfx, d_Seq, f_Seq), 'rpos', r2, Sfx, len(Rid_Long_Seq), d_Seq, Rid_Long_Seq[r1 : r2] )
                    if r2 <= len(Rid_Long_Seq) and d_Seq == Rid_Long_Seq[r1 : r2] and r2 + Sfx >= len(Rid_Long_Seq) - 1:
                        
                        Ref_Long_Seq_Var = Ref_Long_Seq[:f1 + Pfx] + d_Seq + Ref_Long_Seq[f2 - Sfx:]
                        Rid_Long_Seq_Var = Rid_Long_Seq[:r1] + f_Seq + Rid_Long_Seq[r2:]
                        if Rid_Long_Seq_Var in memoization_read_haplotypes[Ref_Long_Seq_Var]: 
                            continue
                        else:
                            memoization_read_haplotypes[Ref_Long_Seq_Var].append(Rid_Long_Seq_Var)
                        RdP = gqv_dbg_utils.seq_interval_obj('RID', r1, r2)
                        RfP = gqv_dbg_utils.seq_interval_obj('REF', f1 + Pfx, f2 - Sfx)
                        SetOfReadIndexes[Read_Index].add_event(RfP, RdP, f_Seq, d_Seq)

            for rid_j in rid_snk_i:
                if rid_j <= 0:
                    continue
                r2_tmp = Rid_Seq_index[rid_j] + OV
                for Pfx, Sfx, d_Seq, f_Seq in vars_from_comb:
                    r2 = r2_tmp - Sfx
                    r1 = r2 - len(d_Seq)
                    if r1 >= 0 and d_Seq == Rid_Long_Seq[r1:r2] and r1 - Pfx - OV <= 1:
                        Ref_Long_Seq_Var = Ref_Long_Seq[:f1 + Pfx] + d_Seq + Ref_Long_Seq[f2 - Sfx:]
                        Rid_Long_Seq_Var = Rid_Long_Seq[:r1] + f_Seq + Rid_Long_Seq[r2:]
                        if Rid_Long_Seq_Var in memoization_read_haplotypes[Ref_Long_Seq_Var]: 
                            continue
                        else:
                            memoization_read_haplotypes[Ref_Long_Seq_Var].append(Rid_Long_Seq_Var)
                        RdP = gqv_dbg_utils.seq_interval_obj('RID', r1, r2)
                        RfP = gqv_dbg_utils.seq_interval_obj('REF', f1 + Pfx, f2 - Sfx)
                        SetOfReadIndexes[Read_Index].add_event(RfP, RdP, f_Seq, d_Seq)



def _find_best_read_path(Read_Index, REF_INDEX, paths, REFPOS, k, SetOfReadIndexes, SetOfReadGroups):  
    min_editd_1 = sys.maxsize
    min_editd_2 = 100.0
    min_n_vars  = 5
    min_map_overlap = sys.maxsize

    final_mapping  = list()
    final_var_list = dict()
    final_path_2   = None

    path_values = defaultdict(dict)

    memoization_alignments = dict()

    Junctions    = SetOfReadIndexes[Read_Index].junctions
    Rid_Long_Seq = SetOfReadIndexes[Read_Index].LongSeq
    Ref_Long_Seq = SetOfReadIndexes[REF_INDEX].LongSeq
    Rid_Long_Seq_ext = 'YY' + Rid_Long_Seq + 'YY'

    path_keys = list(paths.keys())
    path_keys.sort(key = lambda p: p.count('REF'))

    dijkstra_memory = defaultdict(lambda : {100.0 : ['nan']})


    for i, path_coord_tmp in enumerate(path_keys):      
        
        path_editd_1 = path_coord_tmp.count('REF') * 0.1
        path_editd_2 = path_editd_1
        path_vars    = []
        path_n_vars  = 0

        path_coord_split = path_coord_tmp.split(';') 

        last_node = (path_coord_split[-1], ) 
                
        for j, tmp_node in enumerate(path_coord_split):
            if tmp_node.startswith('RID'):
                

                if path_n_vars > min_n_vars:
                    break
                elif path_editd_2 > 30.1 and path_n_vars > 1:
                    break

                onward_path = tuple(path_coord_split[j : ]) 

                if path_editd_2 < min(dijkstra_memory[onward_path]) and path_editd_2 < min(dijkstra_memory[last_node]):
                    try: 
                        dijkstra_memory[onward_path][path_editd_2].append(i)
                    except:
                        dijkstra_memory[onward_path][path_editd_2] = [i]
                else:
                    break

                if j == 0 or j == len(path_coord_split) - 1:
                    continue

            prev_node = path_coord_split[j - 1]
            next_node = path_coord_split[j + 1]

            if tmp_node.startswith('REF'):
                Rid_Src = gqv_dbg_utils.seq_interval_obj(*prev_node.split('|'))
                Rid_Snk = gqv_dbg_utils.seq_interval_obj(*next_node.split('|'))
                Ref_Var = gqv_dbg_utils.seq_interval_obj(*tmp_node.split('|'))
                Ref_Seq = Ref_Var.get_seq(Ref_Long_Seq)
                Rid_Seq = Rid_Long_Seq_ext[Rid_Src.end : Rid_Snk.start]
                args = (Ref_Seq, Rid_Seq, Ref_Var.start, Ref_Var.end)

                RPos, REnd = Rid_Src.end - 2, Rid_Snk.start - 2

                VAR = gqv_vartool.convert_to_var(Ref_Seq, Rid_Seq, Ref_Var.start, Ref_Var.end, 
                                                 0, REnd - RPos, Read_Index, REFPOS)

                if str(VAR) in memoization_alignments:
                    VAROBJ_ED = memoization_alignments[str(VAR)]
                else:
                    VAROBJ_ED = gqv_vartool.calculate_edit_distance(VAR, Junctions, REFPOS)
                    memoization_alignments[str(VAR)] = VAROBJ_ED                

                VAROBJ_ED_copy = copy.deepcopy(VAROBJ_ED)
                for _var in VAROBJ_ED_copy[2]:
                    _var.shift_Rpos(RPos) 

                path_editd_1 += VAROBJ_ED_copy[0]
                path_editd_2 += VAROBJ_ED_copy[1]
                path_vars    += VAROBJ_ED_copy[2] 
                path_n_vars  += VAROBJ_ED_copy[3]

        path_values[i] = [path_editd_1, path_editd_2, path_vars, path_n_vars]

    if dijkstra_memory and path_values:

        final_path_editd = min(dijkstra_memory[last_node])

        for final_path_index in dijkstra_memory[last_node][final_path_editd]:

            if final_path_index != 'nan' and final_path_editd < 100:

                final_path_vars = path_values[final_path_index][2]     

                reference_position = paths[path_keys[final_path_index]] 

                mapping = _align_read_with_vars(Rid_Long_Seq, REFPOS, reference_position, final_path_vars)

                better_overlap, min_map_overlap = _mapping_overlap(SetOfReadIndexes, Read_Index, REF_INDEX,
                                                               mapping, SetOfReadGroups, min_map_overlap)

                if better_overlap:
                    final_var_list = path_values[final_path_index][2]
                    final_mapping  = mapping
                    final_path_2   = final_path_index 
                    final_editd2   = final_path_editd 

    return (final_mapping, final_var_list)


def _align_read_with_vars(Rid_Long_Seq, REFPOS, rc_pos, final_indel_list):
    mapping = []
    mapped_read_pos = REFPOS.genome_pos(max(rc_pos, 2)) 
    
    if "I" in Rid_Long_Seq or not final_indel_list:
        
        if Rid_Long_Seq.strip("XY") in REFPOS.Seq_ext:
            RC_pos = REFPOS.Seq_ext.index(Rid_Long_Seq.strip("XY"))
            prev_RID_end = 2
            prev_REF_end = mapped_read_pos
            mapping = []
            while "I" in Rid_Long_Seq[prev_RID_end : ]:
                i = Rid_Long_Seq.index("I", prev_RID_end)
                RI_len = i - prev_RID_end
                _map = (prev_RID_end, i, prev_REF_end, prev_REF_end + RI_len, "REF")
                mapping.append(_map)
                prev_REF_end = REFPOS.genome_pos(RC_pos + prev_RID_end + RI_len + 1)
                prev_RID_end = i + 1
            
            RI_len = len(Rid_Long_Seq) - prev_RID_end
            _map = (prev_RID_end, len(Rid_Long_Seq), prev_REF_end , prev_REF_end + RI_len, "REF")
            mapping.append(_map)

        return mapping
    
    prev_REF_end = mapped_read_pos
    prev_RID_end = 0

    for variant in sorted(final_indel_list, key = lambda x : x.g_start):
        RID_pos, RID_end = list(variant.RIPos.values())[0] 
        RID_block = Rid_Long_Seq[prev_RID_end : RID_pos].strip("XY")
        REF_block = REFPOS.get_sequence(prev_REF_end, variant.g_start) 
        assert RID_block == REF_block, "{} != {} RC = {} {}".format(RID_block, REF_block, REFPOS.RC_i, list(variant.RIPos.items()))
        mapping.append((prev_RID_end, RID_pos, prev_REF_end, variant.g_start, 'REF'))
        mapping.append((RID_pos, RID_end, variant.g_start, variant.g_end, 'ALT1'))
        prev_RID_end = RID_end
        prev_REF_end = variant.g_end

    RID_block = Rid_Long_Seq[prev_RID_end:].strip("XY")
    REF_block = REFPOS.get_sequence(prev_REF_end, prev_REF_end + len(RID_block))
    assert RID_block == REF_block, '{} != {}'.format(RID_block, REF_block)
    mapping.append((prev_RID_end, len(Rid_Long_Seq), 
                    prev_REF_end, prev_REF_end + len(RID_block), 'REF'))

    return mapping

    
def _mapping_overlap(SetOfReadIndexes, Read_Index, REF_INDEX, mapping, SetOfReadGroups, min_map_overlap):
    RI_map_blocks = [x[2:4] for x in mapping if x[4] == 'REF']
    
    overlap_ratios  = []
    overlap_weights = []
    overlap_diffs   = []

    def _find_intronic_parts(RG_map_blocks):
        cov = defaultdict(int)
        ccov = defaultdict(int)
        introns = []

        for s, e in list(RI_map_blocks) + list(RG_map_blocks):
            cov[s] += 1
            cov[e] -= 1

        cumsum = 0
        for coord in sorted(cov):
            cumsum += cov[coord]
            ccov[coord] = cumsum

        intron_started = False
        for coord in sorted(ccov):
            if ccov[coord] == 0 and not intron_started:
                intron_started = True
                intron_start = coord
            elif ccov[coord] > 0 and intron_started:
                if coord - intron_start > 20:
                    introns.append((intron_start, coord))
                intron_started = False

        return introns

    
    for RG in SetOfReadIndexes[Read_Index].reads:
        diffs = [np.nan]*len(SetOfReadGroups[RG]['Convert_Pos'])

        if RG == 0:
            overlap_diffs.append(diffs)
            continue

        RG_n   = len(SetOfReadGroups[RG]['reads'])
        RG_len = len(SetOfReadGroups[RG]['seq'])
        
        _introns = _find_intronic_parts(SetOfReadGroups[RG]['blocks'])

        for rg_i, rg_pos in enumerate(SetOfReadGroups[RG]['Convert_Pos']):

            if rg_pos >= 0 and rg_pos <= RG_len:
                ri_pos = SetOfReadIndexes[Read_Index].PathIndex[rg_i]
                genome_pos_reassembly = None

                for rid_s, rid_e, gen_s, gen_e, lab in mapping:
                    if lab == 'REF' and rid_s <= ri_pos and rid_e > ri_pos:
                        genome_pos_reassembly = gen_s + (ri_pos - rid_s)

                genome_pos_original = None
                tmp_rg_pos = 0
                for block_s, block_e in SetOfReadGroups[RG]['blocks']:
                    if block_e - block_s + tmp_rg_pos > rg_pos and genome_pos_original is None:
                        genome_pos_original = block_s + (rg_pos - tmp_rg_pos)
                    tmp_rg_pos += block_e - block_s

                if genome_pos_original is None or genome_pos_reassembly is None:
                    continue

                genome_diff = abs(genome_pos_original - genome_pos_reassembly)
                
                junction_diff = 0
                for jxn in _introns:
                    if genome_pos_reassembly <= jxn[0] and jxn[1] <= genome_pos_original:
                        junction_diff = jxn[1] - jxn[0]
                        break
                    elif genome_pos_original <= jxn[0] and jxn[1] <= genome_pos_reassembly:
                        junction_diff = jxn[1] - jxn[0]
                        break

                diffs[rg_i] = abs(genome_diff - junction_diff)

        overlap_diffs.append(diffs)
        
        if len(RI_map_blocks) == 1:
            overlap = 0.0 
        else:
            overlap = np.nansum(diffs)

        overlap_ratios.append(overlap)
        overlap_weights.append(RG_n)

    if not overlap_ratios:
        overlap_ratio_mean = np.nan
    else:
        overlap_ratio_mean = np.average(overlap_ratios, weights = overlap_weights)
    
    if overlap_ratio_mean <= min_map_overlap and overlap_ratio_mean < 500:
        better_overlap = True
        for RG, diffs in zip(SetOfReadIndexes[Read_Index].reads, overlap_diffs):
            SetOfReadGroups[RG]['Mapping_Diff'] = list(diffs)
        min_map_overlap = overlap_ratio_mean
    else:
        better_overlap = False

    return better_overlap, min_map_overlap



def process_var_list(chrom, AllOfVariants, chrom_VAR_LIST):

    new_AllOfVariants = HTSeq.GenomicArrayOfSets([chrom], stranded = False)

    AllOfVariants.sort(key = lambda x : x.g_start)

    n, m = 0, 0 

    filter_counter = Counter()

    for variant in AllOfVariants:        
        if variant.get_variant_type() == "INTRON" :
            continue

        n += 1
        variant_group = VariantGroup[variant.get_variant_type()]
        anchor_pos    = (chrom, variant.g_start + 1, variant_group)
        anchor_allele = (variant.REF, variant.ALT, variant.g_end)

        sm_filters = []
        if anchor_pos not in chrom_VAR_LIST:
            continue
        elif anchor_allele not in chrom_VAR_LIST[anchor_pos]:
            continue


        for feat, val in chrom_VAR_LIST[anchor_pos][anchor_allele].items():

            if isinstance(val, dict):
                sm_filter = val['FEAT']['FILTER']
                sm_filters.append(sm_filter)

        filter_count = sm_filters.count("PASS")
        filter_counter[filter_count] += 1

        if filter_count >= 3:    
            pass 
        else:
            continue
        m += 1
        var_gi = HTSeq.GenomicInterval(chrom, variant.g_start, variant.g_end)
        variant.RIPos = {}
        variant.RGPos = {}
        new_AllOfVariants[var_gi] += variant

    return new_AllOfVariants



def overlap_vars_and_reads_ref(RefVariants , SetOfReadGroups):
    
    for v_i, variant in enumerate(RefVariants):
        RC, AC, AAC, e = (0, 0, 0, 0)

        # We asume all these reads contain the REF allele
        Allele = "REF"

        for RG in SetOfReadGroups:

            found = False
            RG_pos = 0

            for block_s, block_e in SetOfReadGroups[RG]['blocks']:
                if variant.g_start >= block_s and variant.g_end <= block_e:
                    RG_pos += (variant.g_start - block_s) 
                    found = True
                    break
                else:
                    RG_pos += (block_e - block_s)      

            if not found:      
                continue

            RG_seq      = SetOfReadGroups[RG]['seq']
            RG_n        = len(SetOfReadGroups[RG]['reads'])
            RG_len      = len(RG_seq)
            RG_end      = RG_len - (RG_pos + variant.L_ref())

            if RG == 0:
                continue
            elif RG_pos < 0 or RG_pos > RG_len - 1:
                continue
            elif RG_end < 0 or RG_end > RG_len - 1:
                continue

            RC += RG_n

            for Read_Name in SetOfReadGroups[RG]['reads']:
                variant.RGPos[Read_Name] = Allele

        total = float(RC + AC + AAC + e)
        # print("VAR_ref", str(variant), (RC, AC, AAC, e))

    return RefVariants 





def overlap_vars_and_reads(SetOfVariants, SetOfReadIndexes, SetOfReadGroups, REFPOS):
    
    GA_coords = HTSeq.GenomicArrayOfSets([REFPOS.gi.chrom], stranded=False)
    for Read_Index in range(len(SetOfReadIndexes)):
        gStart = SetOfReadIndexes[Read_Index].gStart
        gEnd   = SetOfReadIndexes[Read_Index].gEnd
        if gStart is None or gEnd is None:
            continue
        rGI = HTSeq.GenomicInterval(REFPOS.gi.chrom, gStart, gEnd)
        GA_coords[rGI] += Read_Index


    SetOfVariants.sort(key = lambda x : x.g_start)

    for v_i, variant in enumerate(SetOfVariants):
        RC, AC, AAC, e = (0, 0, 0, 0)

        if variant.get_variant_type() == "INTRON":
            Overlapping_Read_Indexes = variant.RIPos.keys()
        else:
            var_pos = HTSeq.GenomicPosition(REFPOS.gi.chrom, variant.g_start - 1)
            var_end = HTSeq.GenomicPosition(REFPOS.gi.chrom, variant.g_end + 1)
            Overlapping_Read_Indexes = GA_coords[var_pos] | GA_coords[var_end]

        for Read_Index in Overlapping_Read_Indexes:

            if Read_Index in variant.RIPos:
                Allele = 'ALT'
                RI_pos, RI_end = variant.RIPos[Read_Index]
            else:
                RI_pos, RI_end, Allele = SetOfReadIndexes[Read_Index].find_position(variant.rep_region.start, 
                                                                                    variant.rep_region.end)
            if RI_pos is None or RI_end is None:
                continue

            _RI_pos_a, _RI_end_a = (None, None)
            for i, ri_i in enumerate(SetOfReadIndexes[Read_Index].PathIndex):
                if RI_pos >= ri_i:
                    _RI_pos_a, _RI_pos_b = i, RI_pos - ri_i
                if RI_end >= ri_i:
                    _RI_end_a, _RI_end_b = i, RI_end - ri_i


            if _RI_pos_a is None or _RI_end_a is None:
                continue
        
            for RG in SetOfReadIndexes[Read_Index].reads:

                RG_seq      = SetOfReadGroups[RG]['seq']
                RG_pos      = SetOfReadGroups[RG]['Convert_Pos'][_RI_pos_a] + _RI_pos_b
                RG_end      = SetOfReadGroups[RG]['Convert_Pos'][_RI_end_a] + _RI_end_b
                RG_diff_pos = SetOfReadGroups[RG]['Mapping_Diff'][_RI_pos_a] 
                RG_diff_end = SetOfReadGroups[RG]['Mapping_Diff'][_RI_end_a] 
                RG_n        = len(SetOfReadGroups[RG]['reads'])
                RG_len      = len(RG_seq)


                if RG == 0:
                    continue
                elif RG_pos < 0 or RG_pos > RG_len - 1:
                    continue
                elif RG_end < 0 or RG_end > RG_len - 1:
                    continue
                elif RG_diff_pos > 20 or RG_diff_end > 20:
                    continue

                is_error = False
                if Allele == 'ALT':
                    if RG_seq[RG_pos : RG_end] == variant.ALT:
                        AC += RG_n
                    else:
                        is_error = True

                elif Allele == 'REF':
                    rep_seq = REFPOS.get_sequence(variant.rep_region.start, variant.rep_region.end)
                    if rep_seq:
                        RC += RG_n
                    else:
                        is_error = True

                elif Allele == 'ALT1':
                    AAC += RG_n

                if is_error:
                    e += RG_n
                    continue

                for Read_Name in SetOfReadGroups[RG]['reads']:
                    variant.RGPos[Read_Name] = (Allele, RG_pos, RG_end, RG_diff_pos, RG_diff_end)

        total = float(RC + AC + AAC + e)
    del GA_coords

    return [v for v in SetOfVariants if v is not None]



def rough_stutter_noise(varLen, repLen):
    stutter_prob = 0.001

    if repLen >= 5:
        if abs(varLen) == 1:
            stutter_prob = 0.075
        elif abs(varLen) >= 2:
            stutter_prob = 0.035

    return stutter_prob



VariantGroup = {}
for variant_type in ['INSERTION', 'DELETION', 'COMPLEX']:
    VariantGroup[variant_type] = 'INDEL'
for variant_type in ['SNV', "INTRON"]:
    VariantGroup[variant_type] = variant_type



def feature_collection_ref(SetOfVariants, SetOfReadNames, options, SM, REFPOS):

    all2int = {'REF': 0, 'ALT': 1, 'ALT1': 0}
    VAR_LIST  = defaultdict(dict)

    for variant in SetOfVariants:

        variant_type  = variant.get_variant_type()
        variant_group = VariantGroup[variant_type]

        FEAT_ARRAY = list()
        VAR_READ_AND_ALLELES = dict()
        
        while variant.RGPos:
            Read_Name, Allele = variant.RGPos.popitem() 

            ALLELEINT = all2int[Allele]

            FEAT_ARRAY.append(ALLELEINT) 
            
            ReadPairName = ':'.join(Read_Name.split(':')[:-1])

            if options.run_mode != "Full":
                continue
            
            if ReadPairName in VAR_READ_AND_ALLELES:
                if VAR_READ_AND_ALLELES[ReadPairName] == ALLELEINT:
                    pass
                else:
                    VAR_READ_AND_ALLELES.pop(ReadPairName)
            else:
                VAR_READ_AND_ALLELES[ReadPairName] = ALLELEINT

        if not FEAT_ARRAY:
            continue

        ### Allele counts

        AC = FEAT_ARRAY.count(1)
        RC = FEAT_ARRAY.count(0)
        DP = len(FEAT_ARRAY)
        AB = round((AC + options.maxSeqErrorRate)/(DP + options.maxSeqErrorRate), 3)

        
        ### Base quality
        if variant_group == 'SNV':      
            if   (variant.REF, variant.ALT) in [('A', 'G')] and REFPOS.strand in ('+', '.'):
                variant_type = 'RNA_EDITING?'
            elif (variant.REF, variant.ALT) in [('T', 'C')] and REFPOS.strand in ('-', '.'):
                variant_type = 'RNA_EDITING?'

        ### Strand bias

        sttuter_error = rough_stutter_noise(variant.length(), variant.rep_count) 
        error = sttuter_error + options.maxSeqErrorRate

        pb = stats.binomtest(int(AC), int(DP), error).pvalue


        anchor_pos    = (REFPOS.gi.chrom, variant.g_start + 1, variant_group)
        anchor_allele = (variant.REF, variant.ALT, variant.g_end)

        Feature_Vector = {'FILTER': '',
                          'RCL': REFPOS.RC_i,
                          'pb': pb,
                          'dHAP': 0.0,
                          'ndHAP': 0.0,
                          'AC': AC,
                          'RC': RC,
                          'DP': DP,
                          'AB': AB,
                          'MI': np.nan,
                          'MI_n': 0,
                          'NVARS': 1,
                          'PcrDup': 0.0,
                          'OVER_PLOIDY' : False,
                          'n_PLOIDY' : 1,
                          'FISHER_STRAND': 0.01,
                          'READPOS_ALT_mean': np.nan,
                          'BASEQUAL_ALT_mean': np.nan,
                          'REP_COUNT' : variant.rep_count}

        VAR_Info   = {'STRAND'      : REFPOS.strand,
                      'REP_COUNT'   : variant.rep_count,
                      'REP_COORDS'  : variant.rep_region,
                      'INDEL_LEN'   : variant.length(),
                      'variant_type': variant_type,
                      'variant_id'  : REFPOS.gi.chrom + ":" + str(variant),
                      'error'       : error,
                      SM : {'READS' : VAR_READ_AND_ALLELES,
                            'FEAT'  : Feature_Vector} }

        VAR_LIST[anchor_pos][anchor_allele] = VAR_Info

    return VAR_LIST 



def feature_collection(SetOfVariants, SetOfReadNames, options, SM, REFPOS):

    all2int = {'REF': 0, 'ALT': 1, 'ALT1': 0}
    VAR_LIST  = defaultdict(dict)

    for variant in SetOfVariants:

        variant_type  = variant.get_variant_type()
        variant_group = VariantGroup[variant_type]

        FEAT_ARRAY = list()
        VAR_READ_AND_ALLELES = dict()
        
        while variant.RGPos:
            Read_Name, (Allele, RG_pos, RG_end, RG_diff_pos, RG_diff_end) = variant.RGPos.popitem() 

            Read_is_reverse = SetOfReadNames[Read_Name]['Dir']
            Read_Qualities  = SetOfReadNames[Read_Name]['Quals']
            Read_is_dup     = SetOfReadNames[Read_Name]['PcrDup']

            RL = len(Read_Qualities)

            QUAL_around = np.nan
        
            if Allele == "REF":
                QUAL_around = np.nan
            elif variant_type == 'SNV':
                QUAL_around = Read_Qualities[RG_pos]
            else:
                if Read_is_reverse:
                    QUAL_region = Read_Qualities[max(RG_pos - 5, 0) : RG_pos + 1]
                else:
                    RG_rep_end  = RG_pos + variant.rep_region.length
                    QUAL_region = Read_Qualities[RG_rep_end - 1 : RG_rep_end + 5]

                if QUAL_region:
                    QUAL_around = np.nanmean(QUAL_region) 
                else:
                    QUAL_around = np.nanmean(Read_Qualities) 
                
            

            READPOS   = min(RL - 1 - RG_end, RG_pos) 
            ALLELEINT = all2int[Allele]
            ALLELEDIR = ALLELEINT * 10 + int(Read_is_reverse)


            FEAT_ARRAY.append([ALLELEINT, QUAL_around, READPOS, ALLELEDIR, RG_diff_pos, RG_diff_end, Read_is_dup]) 
            
            if options.run_mode != "Full":
                continue
                
            ReadPairName = ':'.join(Read_Name.split(':')[:-1])
            
            if ReadPairName in VAR_READ_AND_ALLELES:
                if VAR_READ_AND_ALLELES[ReadPairName] == ALLELEINT:
                    pass
                else:
                    VAR_READ_AND_ALLELES.pop(ReadPairName)
            else:
                VAR_READ_AND_ALLELES[ReadPairName] = ALLELEINT

        if not FEAT_ARRAY:
            continue

        FEAT_ARRAY = np.array(FEAT_ARRAY)
        colNames = {'ALLELE': 0, 'DNQ': 1, 'READPOS': 2, 'STRAND': 3, 'MAP_DIFF_POS': 4, 'MAP_DIFF_END': 5, "DUP": 6}

        ### Allele counts
        alleles, counts = np.unique(FEAT_ARRAY[:, colNames['ALLELE']], return_counts = True)
        allele_counts = defaultdict(lambda : 0, zip(alleles, counts))


        AC = allele_counts[1]
        RC = allele_counts[0]
        DP = sum(allele_counts.values())
        AB = round((AC + options.maxSeqErrorRate)/(DP + options.maxSeqErrorRate), 3)

        ### Read position
        READPOS_ALT_mean = np.nanmean(FEAT_ARRAY[FEAT_ARRAY[:, 0] == 1, colNames['READPOS']])

        ### PCR-dup
        DUP_ALT_mean = np.nanmean(FEAT_ARRAY[FEAT_ARRAY[:, 0] == 1, colNames['DUP']])
        
        ### Base quality
        if variant_group == 'SNV':      
            BASEQUAL_ALT_mean = np.nanmedian(FEAT_ARRAY[FEAT_ARRAY[:, 0] == 1, colNames['DNQ']])
            if   (variant.REF, variant.ALT) in [('A', 'G')] and REFPOS.strand in ('+', '.'):
                variant_type = 'RNA_EDITING?'
            elif (variant.REF, variant.ALT) in [('T', 'C')] and REFPOS.strand in ('-', '.'):
                variant_type = 'RNA_EDITING?'
        
        elif variant_group == 'INDEL':
            BASEQUAL_ALT_mean = np.nanmedian(FEAT_ARRAY[FEAT_ARRAY[:, 0] == 1, colNames['DNQ']]) 
        else:            
            BASEQUAL_ALT_mean = np.nan

        ### Strand bias
        _strands, counts = np.unique(FEAT_ARRAY[:, colNames['STRAND']], return_counts = True)
        allele_directions = defaultdict(lambda : 0, zip(_strands, counts))
        FS, pFS = stats.fisher_exact([[allele_directions[0] + 1, allele_directions[1] + 1],
                                      [allele_directions[10] + 1, allele_directions[11] + 1]])


        sttuter_error = rough_stutter_noise(variant.length(), variant.rep_count) 
        error = sttuter_error + options.maxSeqErrorRate

        pb = stats.binomtest(int(AC), int(DP), error).pvalue


        anchor_pos    = (REFPOS.gi.chrom, variant.g_start + 1, variant_group)
        anchor_allele = (variant.REF, variant.ALT, variant.g_end)

        Feature_Vector = {'FILTER': '',
                          'RCL': REFPOS.RC_i,
                          'pb': pb,
                          'dHAP': 0.0,
                          'ndHAP': 0.0,
                          'AC': AC,
                          'RC': RC,
                          'DP': DP,
                          'AB': AB,
                          'MI': np.nan,
                          'MI_n': 0,
                          'NVARS': 1,
                          'PcrDup': DUP_ALT_mean,
                          'OVER_PLOIDY' : False,
                          'n_PLOIDY' : 1,

                          'FISHER_STRAND': FS,
                          
                          'READPOS_ALT_mean': READPOS_ALT_mean,
                          'BASEQUAL_ALT_mean': BASEQUAL_ALT_mean,
                          'REP_COUNT' : variant.rep_count}

        VAR_Info   = {'STRAND'    : REFPOS.strand,
                      'REP_COUNT' : variant.rep_count,
                      'REP_COORDS': variant.rep_region,
                      'INDEL_LEN' : variant.length(),
        	          'variant_type': variant_type,
                      'variant_id': REFPOS.gi.chrom + ":" + str(variant),
                      'error' : error,
                      SM : {'READS' : VAR_READ_AND_ALLELES,
                            'FEAT'  : Feature_Vector} }

        VAR_LIST[anchor_pos][anchor_allele] = VAR_Info

    return VAR_LIST 
 


def table_genotyping(VAR_LIST, SM, options):

    ReadPairNames = set()

    RC_variants = []
    for anchor_pos, anchor_alleles in VAR_LIST.items():
        for allele, sm_attr in anchor_alleles.items():
            if sm_attr["variant_type"] in ('COMPLEX', 'INSERTION', 'DELETION', 'SNV'):
                RC_variants.append([anchor_pos, allele])
            ReadPairNames |= set(sm_attr[SM]["READS"])

    ReadPairNames = list(ReadPairNames)

    RC_variants.sort(key = lambda h : h[0][1])
    n = len(RC_variants)
    m = len(ReadPairNames)

    if n >= 1 and m >= 1:

        rc_var_strings = ['{0[0][1]}:{0[1][0]}>{0[1][1]}'.format(v) for v in RC_variants]

        mat = np.zeros((m, n))
        
        POS_ALLELES = defaultdict(list)

        w = []

        for j, (POS, ALLELE) in enumerate(RC_variants):
            POS_ALLELES[POS].append((j, ALLELE))
            w.append(VAR_LIST[POS][ALLELE][SM]['FEAT']['BASEQUAL_ALT_mean'])

            for i, r in enumerate(ReadPairNames):
                var_reads = VAR_LIST[POS][ALLELE][SM]["READS"]
                if r in var_reads:
                    mat[i][j] = var_reads[r] + 1
                else:
                    mat[i][j] = 0

        mat = mat[~np.all(mat == 0, axis = 1), :]
        mat[mat == 0] = np.nan

        w = w/np.nanmax(w)
        w[np.isnan(w)] = np.nanmean(w)

        h_mean, h_std, haplotypes, clusters = _k_means_clustering(mat, options.ploidy, w)
        haplotypes = np.round(haplotypes, 0) 

        #### Haplotype difference
        hap_mat  = np.zeros(mat.shape)
        for hap_i, read_dict in clusters.items():
            for read_i in read_dict:
                hap_mat[read_i] = haplotypes[hap_i]
        
        diff_mat = hap_mat - mat

        # print(rc_var_strings)
        # print(haplotypes, h_mean, h_std)
        # print(w)
        # for i in range(mat.shape[0]):
        #     print(','.join(map(str, mat[i])))
    
        #### vars per read
        read_var_n = [0]*mat.shape[0]
        for i in range(mat.shape[0]):
            read_var_n[i] = np.nansum(mat[i] - 1.0)

        #### Remove over-ploidy alleles
        for POS, ALLELE_columns in POS_ALLELES.items():

            ACs = dict((ALLELE, np.nansum(hap_mat[:, j] - 1.0)) for j, ALLELE in ALLELE_columns)
            RCmat = hap_mat[:, [c[0] for c in ALLELE_columns]]
            # rows that have no ALTs 
            RCs = np.sum(~(RCmat - 1).any(axis = 1))
            ACs['REF'] = RCs 

            

            for j, ALLELE in ALLELE_columns:
                VAR_LIST[POS][ALLELE][SM]['FEAT']['n_PLOIDY'] = len(ACs) 

            for _ in range(min(options.ploidy, len(ACs))):
                max_allele = max(ACs, key = ACs.get)
                ACs.pop(max_allele)

            for ALLELE, AC in ACs.items():
                if ALLELE == "REF": 
                    continue
                VAR_LIST[POS][ALLELE][SM]['FEAT']['OVER_PLOIDY'] = True 
                

        #### Update Feats
        for j, (POS, ALLELE) in enumerate(RC_variants):
            hap_counts = hap_mat[ ~np.isnan(mat[:, j]), j] - 1.0

            if ALLELE not in VAR_LIST[POS]:
                continue
            VAR_LIST[POS][ALLELE][SM]['FEAT']['AC']    = np.nansum(hap_counts)
            VAR_LIST[POS][ALLELE][SM]['FEAT']['AB']    = np.nanmean(hap_counts)
            VAR_LIST[POS][ALLELE][SM]['FEAT']['DP']    = len(hap_counts)
            VAR_LIST[POS][ALLELE][SM]['FEAT']['RC']    = np.nansum(hap_counts == 0)
            VAR_LIST[POS][ALLELE][SM]['FEAT']['dHAP']  = np.nansum(diff_mat[:, j])
            VAR_LIST[POS][ALLELE][SM]['FEAT']['ndHAP'] = np.nanmean(diff_mat[:, j])
            VAR_LIST[POS][ALLELE][SM]['FEAT']['NVARS'] = np.nansum((mat[:, j] - 1.0) * read_var_n)/np.nansum(mat[:, j] - 1.0)


        del mat
        del diff_mat

    return VAR_LIST




def _k_means_clustering(mat, ploidy, w, iterations = 5, random_initialization = 3):
    m, n = mat.shape
    np.random.seed(1)
    j_rand = list(range(n))

    score_track = {}
   
    for r_iter in range(max(random_initialization, n)):
        
        centroids = np.zeros((ploidy, n))

        i_rand = np.random.randint(0, ploidy, n)
        
        centroids[i_rand, j_rand] = 1
        centroids = centroids + 1

        prev_centroids = centroids.copy()
        
        for j_iter in range(iterations):
            clusters = defaultdict(dict)

            for i in range(m):
                min_distance = 1e6
                for p_i, centroid in enumerate(centroids):
                    d = np.square((centroid - mat[i]) * w)
                    d = np.sqrt(np.nansum(d))
                    if d < min_distance:
                        min_distance = d 
                        best_centroid = p_i 

                clusters[best_centroid][i] = min_distance

            read_dist2hap = []

            for c_i, read_dict in clusters.items():
                tmp = np.nanmean(mat[list(read_dict), :], axis = 0)
                tmp[np.isnan(tmp)] = centroids[c_i][np.isnan(tmp)]
                centroids[c_i] = tmp

                read_dist2hap.extend(read_dict.values())
            
            if (centroids == prev_centroids).all():
                break
            else:
                prev_centroids = centroids.copy()

        score_track[r_iter] = [np.mean(read_dist2hap),
                             np.std(read_dist2hap),
                             centroids,
                             clusters]

    best_r_i = min(score_track, key=lambda h: score_track[h][0])
    return score_track[best_r_i]


def merge_SM_VAR_LIST( VAR_LIST):

    merged_VAR_LIST = defaultdict(dict)

    for POS, ALLELE_dict in VAR_LIST.items():
        # skip intronic parts
        if "INTRONIC_PART" in POS:
            continue

        for ALLELE, SM_dict in ALLELE_dict.items():     

            merged_VAR_LIST[POS][ALLELE] = {'merged_SM' : {"READS" : {}, "FEAT" : {}}}

            AC, DP, RC = 0, 0, 0

            fvars = []
            over_ploidy = []
            for SM, SM_attr in SM_dict.items():
                if not isinstance(SM_attr, dict):
                    merged_VAR_LIST[POS][ALLELE][SM] = SM_attr
                    continue
                # Adding samples
                merged_VAR_LIST[POS][ALLELE][SM] = {"FEAT" : copy.deepcopy(SM_attr['FEAT'])}
                # Pool reads
                merged_VAR_LIST[POS][ALLELE]['merged_SM']["READS"].update(SM_attr.pop("READS"))
                # Sum Feats
                AC += SM_attr['FEAT']['AC']
                DP += SM_attr['FEAT']['DP']
                RC += SM_attr['FEAT']['RC']

                fvars.append([SM_attr['FEAT']['READPOS_ALT_mean'], 
                              SM_attr['FEAT']['BASEQUAL_ALT_mean'], 
                              SM_attr['FEAT']['NVARS'], 
                              SM_attr['FEAT']['ndHAP'],
                              SM_attr['FEAT']['n_PLOIDY']])

                over_ploidy.append(SM_attr['FEAT']['OVER_PLOIDY'])

            merged_VAR_LIST[POS][ALLELE]['merged_SM']["FEAT"]["AC"] = AC
            merged_VAR_LIST[POS][ALLELE]['merged_SM']["FEAT"]["DP"] = DP
            merged_VAR_LIST[POS][ALLELE]['merged_SM']["FEAT"]["RC"] = RC
            merged_VAR_LIST[POS][ALLELE]['merged_SM']["FEAT"]["AB"] = (AC + 0.01)/(DP + 0.01)

            m_readpos, m_basequal, m_nvars, m_ndhap, m_n_ploidy = np.nanmean(np.array(fvars), axis = 0)

            merged_VAR_LIST[POS][ALLELE]['merged_SM']["FEAT"]['READPOS_ALT_mean']  = m_readpos
            merged_VAR_LIST[POS][ALLELE]['merged_SM']["FEAT"]['BASEQUAL_ALT_mean'] = m_basequal
            merged_VAR_LIST[POS][ALLELE]['merged_SM']["FEAT"]['NVARS']             = m_nvars
            merged_VAR_LIST[POS][ALLELE]['merged_SM']["FEAT"]['ndHAP']             = m_ndhap 
            merged_VAR_LIST[POS][ALLELE]['merged_SM']["FEAT"]['n_PLOIDY']          = m_n_ploidy

            merged_VAR_LIST[POS][ALLELE]['merged_SM']["FEAT"]['OVER_PLOIDY'] = any(over_ploidy)

            if DP == 0:
                pb = 1.0
            else:
                pb = stats.binomtest(int(AC), int(DP), merged_VAR_LIST[POS][ALLELE]['error']).pvalue

            merged_VAR_LIST[POS][ALLELE]['merged_SM']["FEAT"]['pb'] = pb
        
    return merged_VAR_LIST

            

def read_genome_fasta(gf):
    fasta = pyfaidx.Faidx(gf)
    return fasta


def write_feat_file(var_list, FeatFile):
    if os.path.exists(FeatFile):
        Feat_Handle = open(FeatFile, 'a')
        Title = True
    else:
        Feat_Handle = open(FeatFile, 'w')
        Title = False

    feat_list = []          

    for genomic_pos in sorted(var_list):
        for ALT, line in var_list[genomic_pos].items():

            flat_dict = {} 

            for feat, v in line.items():
                if not isinstance(v, dict):
                    flat_dict[feat] = v

            for sm, v in line.items():
                if isinstance(v, dict):
                                        
                    for feat, value in line[sm]['FEAT'].items():
                        flat_dict[feat] = value

                    feat_list  = sorted(list(flat_dict.keys()))

                    if not Title: # title is missing
                        title_line = [str(_feat) for _feat in feat_list] ; print(title_line)
                        Feat_Handle.write('\t'.join(["SM"] + title_line) + '\n')
                        Title = True

                    if abs(flat_dict['INDEL_LEN']) > 20:
                        continue

                    outline = [sm] + [str(flat_dict[_feat]) for _feat in feat_list] 
                    Feat_Handle.write('\t'.join(outline) + '\n')

    Feat_Handle.close()


def write_readcluster_file(ALL_READ_CLUSTERS, outfile):
    if os.path.exists(outfile):
        f = open(outfile, 'a')
    else:
        f = open(outfile, 'w')
        f.write("bam_file\tindex\tchrom\tstrand\tstart\tend\tmax_coverage\n")

    for (sm, bam), rc_list in ALL_READ_CLUSTERS.items():
        for rc in rc_list:
            outline = [sm] + list(rc)
            f.write("\t".join(map(str, outline)) + "\n")

    f.close()


def write_vcf_header(search_regions, genome_fasta, VcfFile, SM_names):

    genome_faidx = read_genome_fasta(genome_fasta)

    VCF_handle = open(VcfFile, 'w')
    VCF_handle.write('##fileformat=VCFv4.2\n')
    VCF_handle.write('##fileDate={}\n'.format(datetime.date.today()))
    VCF_handle.write('##CL=python {}">\n'.format(' '.join(sys.argv)))
    
    for chrom in set([search_region[0] for search_region in search_regions]):
        chrom_len = genome_faidx.index[chrom]['rlen']
        VCF_handle.write('##contig=<ID={},length={}>\n'.format(chrom, chrom_len))

    VCF_handle.write('##FILTER=<ID=PASS,Description="Variant passes all filters">\n')
    VCF_handle.write('##FILTER=<ID=LowQual,Description="Variant does not pass one or more filtering criteria">\n')
    VCF_handle.write('##INFO=<ID=varType,Number=1,Type=String,Description="Variant type">\n')
    VCF_handle.write('##INFO=<ID=varGroup,Number=1,Type=String,Description="Variant Group">\n')
    VCF_handle.write('##INFO=<ID=varLength,Number=1,Type=Integer,Description="Variant length">\n')
    VCF_handle.write('##INFO=<ID=varEnd,Number=1,Type=Integer,Description="End position of the variant">\n')
    VCF_handle.write('##INFO=<ID=TandemRep,Number=1,Type=Integer,Description="Count of repeat sequence adjacent to the variant">\n')
    VCF_handle.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of reads overlapping the variant position">\n')
    VCF_handle.write('##FORMAT=<ID=AC,Number=R,Type=Integer,Description="Number of reads containing the variant allele">\n')
    VCF_handle.write('##FORMAT=<ID=RC,Number=R,Type=Integer,Description="Number of reads containing the reference allele">\n')
    VCF_handle.write('##FORMAT=<ID=AB,Number=1,Type=Float,Description="Allelic balance of the variant">\n')
    VCF_handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Consensus genotype">\n')
    VCF_handle.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
    VCF_handle.write('##FORMAT=<ID=RCL,Number=1,Type=Integer,Description="Read Cluster id">\n')

    VCF_handle.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format('\t'.join(SM[1] for SM in SM_names)))



def write_vcf_file(var_list, VcfFile, SM_names):

    VCF_handle = open(VcfFile, 'a')

    for VAR in sorted(var_list):

        (CHROM, POS, VarGroup) = VAR

        if VarGroup in ("INTRON", "INTRONIC_PART"):
            continue

        for (REF, ALT, END), SAMPLE_dict in var_list[VAR].items():

            quals   = []
            filters = []
            for SM in SM_names + ["merged_SM"]:
                try:
                    quals.append(SAMPLE_dict[SM]["FEAT"]["QUAL"])
                    filters.append(SAMPLE_dict[SM]["FEAT"]["FILTER"])
                except:
                    quals.append(np.nan)
                    filters.append('')

            QUAL = np.nanmax(np.array(quals))

            if any(_filter == "PASS" for _filter in filters):
                FILTER = 'PASS'
            else:
                FILTER = 'LowQual'

            outline_basic = f'{CHROM}\t{POS}\t.\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t'

            outline_info  = f'varEnd={END};varGroup={VarGroup};'
            outline_info += 'varType={variant_type};varLength={INDEL_LEN};TandemRep={REP_COUNT}\t'.format(**SAMPLE_dict)

            outline_format = 'GT:GQ:DP:AC:RC:AB:RCL\t'
            outline_values =  ''
            
            for (bamf, SM) in SM_names:
                if SM in SAMPLE_dict:
                    outline_values += './.:0:{DP:.0f}:{AC:.0f}:{RC:.0f}:{AB:.2f}:{RCL:.0f}\t'.format(**SAMPLE_dict[SM]['FEAT'])
                else:
                    outline_values += './.:.:.:.:.:.:.\t'
            
            outline_values = outline_values.strip() 

            VCF_handle.write(outline_basic + outline_info + outline_format + outline_values + '\n')

    VCF_handle.close()


def write_mutinfo_header(outfile):
    with open(outfile, 'w') as f:
        f.write("SAMPLE\tVAR1_TYPE\tVAR2_TYPE\tVAR1\tQUAL\tMI(mean)\tVAR2(n)\tVAR2\tVAR2_MI\tVAR2(cov)\tSTRAND\n")
       

    
def write_mutinfo_file(mi_outlines, outfile):
    with open(outfile, 'a') as f:
        for mi_outline in mi_outlines:
            f.write("\t".join(map(str, mi_outline)) + "\n")


def write_introns_header(BedFile):
    BED_handle = open(BedFile, 'w')
    bed_line = f'CHROM\tIP_start\tIP_end\tSM\tn\tstrand\tintrons_string\n'
    BED_handle.write(bed_line)
    BED_handle.close()


def write_introns_bed(var_list, BedFile):

    BED_handle = open(BedFile, 'a')

    for VAR in sorted(var_list):

        (CHROM, IP_start, VarGroup) = VAR

        if VarGroup != "INTRONIC_PART":
            continue

        for (REF, ALT, IP_end), SAMPLE_dict in var_list[VAR].items():

            outline_values = ''

            for SM, SM_val in SAMPLE_dict.items():
                if not isinstance(SM_val, dict):
                    continue
                introns    = SM_val['IP_ALLELES']
                introns_ac = SM_val['IP_ALLELES_AC']

                n = len(introns)

                introns_string = ''
                for (intron, intron_ac) in zip(introns, introns_ac):

                    introns_string += ';{}-{}:{}'.format(intron[0][1], intron[1][2], intron_ac)
                
                bed_line = f'{CHROM}\t{IP_start}\t{IP_end}\t{SM}\t{n}\t.\t{introns_string}\n'
                BED_handle.write(bed_line)

    BED_handle.close()


