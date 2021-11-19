#!/usr/bin/python

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
import networkx as nx
import HTSeq

from . import gqv_software_management as gqv_SftMngr
from . import gqv_dbg_utils
from . import gqv_glm
from . import gqv_vartool




def Process_de_Bruijn_Graph(SetOfReadGroups, k, db, CUTOFFS):
    node_First  = SetOfReadGroups[0]['seq'][:k - 1]
    node_Last   = SetOfReadGroups[0]['seq'][-(k - 1):]
    node_Second = db.successors(node_First).__next__()
    
    xdb = nx.DiGraph()
    gqv_dbg_utils.dbg_compress_graph(db, xdb, k, node_First, node_Second)

    min_Allele_Ratio = min(CUTOFFS["P_maxSeqErrorRate"] * 2, 0.05)

    all_nodes = list(xdb.nodes)[:]

    gqv_dbg_utils.dbg_remove_alternative_ends(all_nodes, xdb, node_First, node_Last)

    Graph_Read_Tracker = defaultdict(dict)
    gqv_dbg_utils.dbg_get_read_paths(xdb, Graph_Read_Tracker, k)

    SetOfReadIndexes = gqv_dbg_utils.compress_read_paths(xdb, k, Graph_Read_Tracker, 
                                                            SetOfReadGroups)

    gqv_SftMngr.rm_nested_dict(Graph_Read_Tracker)

    return (xdb, SetOfReadIndexes)




def assign_vars_to_reads(src_snk_pairs, Ref_Index, SetOfReadIndexes, SetOfReadGroups, 
        k, REFPOS):
    
    OV = k - 2
    memoization_vars = defaultdict(set)
    memoization_trims = dict()
    SetOfVariants = list()
    
    for Read_Index, bubble_list_tmp in src_snk_pairs.items():
        _match_src_snk_to_reads(Read_Index, Ref_Index, bubble_list_tmp, SetOfReadIndexes, 
                                memoization_trims, memoization_vars, ronda = 1, RC_i = None)

    for Read_Index in range(len(SetOfReadIndexes)):

        try:
            bubble_list_tmp = src_snk_pairs.pop(Read_Index)
        except KeyError:
            bubble_list_tmp = {}

        _match_src_snk_to_reads(Read_Index, Ref_Index, bubble_list_tmp, SetOfReadIndexes, 
                                memoization_trims, memoization_vars, ronda = 2, RC_i = REFPOS.RC_i)
        
        Graph, source, target = gqv_dbg_utils.make_read_dbg(Read_Index, Ref_Index, 
                                                            SetOfReadIndexes, REFPOS.RC_i)

        SetOfReadIndexes[Read_Index].delete_events()
        edges = defaultdict(list)
        
        all_paths = gqv_dbg_utils.Dijkstra_find_end_to_end_paths(Read_Index, Ref_Index, 
                                                                SetOfReadIndexes, Graph, 
                                                                source, target, edges, REFPOS.RC_i)

        RI_alignment, var_list = _find_best_read_path(Read_Index, Ref_Index, all_paths, 
                                                        REFPOS, k, SetOfReadIndexes, 
                                                        SetOfReadGroups)
        # r_total += SetOfReadIndexes[Read_Index].n
        # if not RI_alignment:
        #     r_unmapped += SetOfReadIndexes[Read_Index].n

        SetOfVariants.extend(var_list)        
        SetOfReadIndexes[Read_Index].add_mapping(RI_alignment)
        
    mms = [memoization_vars, memoization_trims]
    map(gqv_SftMngr.rm_nested_dict, mms)

    return SetOfVariants



def _src_snk_comb(src_i_list, snk_i_list, src_seq, snk_seq, k):

    effective_src_i_list = src_i_list[:]
    effective_snk_i_list = snk_i_list[:]

    if len(src_i_list) > 5:
        for i in range(len(src_i_list) - 1):
            if src_i_list[i] - src_i_list[i - 1] == 1 and src_i_list[i + 1] - src_i_list[i] == 1:
                effective_src_i_list[i] = 'NA'
        effective_src_i_list = [x for x in effective_src_i_list if x != 'NA']

    if len(snk_i_list) > 5:
        for i in range(1, len(snk_i_list) - 1):
            if snk_i_list[i] - snk_i_list[i - 1] == 1 and snk_i_list[i + 1] - snk_i_list[i] == 1:
                effective_snk_i_list[i] = 'NA'
        effective_snk_i_list = [x for x in effective_snk_i_list if x != 'NA']


    if 'N' in src_seq or 'N' in snk_seq:
        combs = []

    elif len(effective_snk_i_list) > 5 and len(effective_src_i_list) > 5:
        space1 = [effective_snk_i_list[i + 1] - effective_snk_i_list[i] for i in range(len(effective_snk_i_list) - 1)]
        space2 = [effective_src_i_list[i + 1] - effective_src_i_list[i] for i in range(len(effective_src_i_list) - 1)]
        if abs(np.mean(space1) - np.mean(space2)) < 1:
            combs = []
            for i in range(len(effective_snk_i_list[:-1])):
                snk_i = effective_snk_i_list[i + 1]
                if snk_i >= effective_src_i_list[0]:
                    combs.append((effective_src_i_list[0], snk_i))
        else:
            combs = []
            for s_i, s_j in itertools.product(effective_src_i_list, effective_snk_i_list):
                delta = s_j - s_i
                if delta >= 0 and delta < 2*k:
                    combs.append((s_i, s_j))
    else:
        combs = []
        for s_i, s_j in itertools.product(effective_src_i_list, effective_snk_i_list):
            if s_j >= s_i:
                combs.append((s_i, s_j))

    return combs



def _match_src_snk_to_reads(Read_Index, REF_Index, bubble_list_tmp, SetOfReadIndexes, 
        memoization_trims, memoization_vars, ronda, RC_i):

    k = SetOfReadIndexes[REF_Index].k
    OV = k - 2
    Ref_Seq_path, Ref_node_index = SetOfReadIndexes[REF_Index].PathSeq(return_index=True)
    Ref_Seq_index = SetOfReadIndexes[REF_Index].PathIndex
    Ref_Long_Seq  = SetOfReadIndexes[REF_Index].LongSeq
    
    Rid_Seq_path, Rid_node_index = SetOfReadIndexes[Read_Index].PathSeq(return_index=True)
    Rid_Seq_index = SetOfReadIndexes[Read_Index].PathIndex
    Rid_Long_Seq  = SetOfReadIndexes[Read_Index].LongSeq
    
    memoization_read_haplotypes = []


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
                # print('0_srcsnk',Read_Index, SrcSnk_Rid_Seq , src_seq, snk_seq, SrcSnk_Rid_Seq in SrcSnk_Path_Seqs, SetOfReadIndexes[Read_Index].reads)
                if not SrcSnk_Rid_Seq in SrcSnk_Path_Seqs:
                    continue

                comb_found = True

                Ref_Long_Seq_Var = Ref_Long_Seq[:f1] + SrcSnk_Rid_Seq + Ref_Long_Seq[f2:]
                Rid_Long_Seq_Var = Rid_Long_Seq[:r1] + SrcSnk_Ref_Seq + Rid_Long_Seq[r2:]
                local_haplotype = (Ref_Long_Seq_Var, Rid_Long_Seq_Var)

                if local_haplotype in memoization_read_haplotypes: ## this stept akes long? 
                    continue
                else:
                    memoization_read_haplotypes.append(local_haplotype)

                try:
                    V_t = memoization_trims[SrcSnk_Rid_Seq, SrcSnk_Ref_Seq]
                except:
                    V_t = gqv_vartool.trim_sequences(SrcSnk_Rid_Seq, SrcSnk_Ref_Seq, first = 'sufix')
                    memoization_trims[SrcSnk_Rid_Seq, SrcSnk_Ref_Seq] = V_t

                Pfx, Sfx, d_Seq, f_Seq = V_t

                if len(d_Seq) - len(f_Seq) > 20:
                    continue
                memoization_vars[src_seq, snk_seq, SrcSnk_Ref_Seq, f1, f2].add(V_t)
                if ronda == 2:
                    RdP = gqv_dbg_utils.seq_interval_obj('RID', r1 + Pfx, r2 - Sfx)
                    RfP = gqv_dbg_utils.seq_interval_obj('REF', f1 + Pfx, f2 - Sfx)
                    SetOfReadIndexes[Read_Index].add_event(RfP, RdP, f_Seq, d_Seq)

            if ronda == 1:
                for SrcSnk_Path_Seq in SrcSnk_Path_Seqs:
                    try:
                        V_t = memoization_trims[SrcSnk_Path_Seq, SrcSnk_Ref_Seq]
                    except:
                        V_t = gqv_vartool.trim_sequences(SrcSnk_Path_Seq, SrcSnk_Ref_Seq, first = 'sufix')
                        memoization_trims[SrcSnk_Path_Seq, SrcSnk_Ref_Seq] = V_t

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
                    if r2 <= len(Rid_Long_Seq) and d_Seq == Rid_Long_Seq[r1 : r2] and r2 + Sfx >= len(Rid_Long_Seq) - 1:
                        
                        Ref_Long_Seq_Var = Ref_Long_Seq[:f1 + Pfx] + d_Seq + Ref_Long_Seq[f2 - Sfx:]
                        Rid_Long_Seq_Var = Rid_Long_Seq[:r1] + f_Seq + Rid_Long_Seq[r2:]
                        local_haplotype = (Ref_Long_Seq_Var, Rid_Long_Seq_Var)
                        if local_haplotype in memoization_read_haplotypes:
                            continue
                        else:
                            memoization_read_haplotypes.append(local_haplotype)
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
                        local_haplotype = (Ref_Long_Seq_Var, Rid_Long_Seq_Var)
                        if local_haplotype in memoization_read_haplotypes:
                            continue
                        else:
                            memoization_read_haplotypes.append(local_haplotype)
                        RdP = gqv_dbg_utils.seq_interval_obj('RID', r1, r2)
                        RfP = gqv_dbg_utils.seq_interval_obj('REF', f1 + Pfx, f2 - Sfx)
                        SetOfReadIndexes[Read_Index].add_event(RfP, RdP, f_Seq, d_Seq)


def _find_best_read_path(Read_Index, Ref_Index, paths, REFPOS, k, SetOfReadIndexes, SetOfReadGroups):  
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
    Ref_Long_Seq = SetOfReadIndexes[Ref_Index].LongSeq
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

                better_overlap, min_map_overlap = _mapping_overlap(SetOfReadIndexes, Read_Index, Ref_Index,
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
        
        if Rid_Long_Seq.strip("XY") in REFPOS.Seq:
            RC_pos = REFPOS.Seq.index(Rid_Long_Seq.strip("XY"))
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
        RID_pos, RID_end = list(variant.ReadPos.values())[0] 
        RID_block = Rid_Long_Seq[prev_RID_end : RID_pos].strip("XY")
        REF_block = REFPOS.get_sequence(prev_REF_end, variant.g_start) 
        assert RID_block == REF_block, "{} != {} RC = {} {}".format(RID_block, REF_block, REFPOS.RC_i, list(variant.ReadPos.items()))
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

    
def _mapping_overlap(SetOfReadIndexes, Read_Index, Ref_Index, mapping, SetOfReadGroups, min_map_overlap):
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
            overlap = 0.0 ###### why???????????????????
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





def overlap_vars_and_reads(SetOfVariants, SetOfReadIndexes, SetOfReadGroups, REFPOS):
    
    GA_coords = HTSeq.GenomicArrayOfSets([REFPOS.chrom], stranded=False)
    for Read_Index in range(len(SetOfReadIndexes)):
        gStart = SetOfReadIndexes[Read_Index].gStart
        gEnd   = SetOfReadIndexes[Read_Index].gEnd
        if gStart is None or gEnd is None:
            continue
        rGI = HTSeq.GenomicInterval(REFPOS.chrom, gStart, gEnd)
        GA_coords[rGI] += Read_Index


    SetOfVariants.sort(key = lambda x : x.g_start)

    for v_i, variant in enumerate(SetOfVariants):
        RC, AC, AAC, e = (0, 0, 0, 0)
        if variant.get_variant_type() == "INTRON":
            Overlapping_Read_Indexes = variant.ReadPos.keys()
        else:
            var_pos = HTSeq.GenomicPosition(REFPOS.chrom, variant.g_start - 1)
            var_end = HTSeq.GenomicPosition(REFPOS.chrom, variant.g_end + 1)
            Overlapping_Read_Indexes = GA_coords[var_pos] | GA_coords[var_end]

        for Read_Index in Overlapping_Read_Indexes:

            if Read_Index in variant.ReadPos:
                Allele = 'ALT'
                RI_pos, RI_end = variant.ReadPos[Read_Index]
            else:
                RI_pos, RI_end, Allele = SetOfReadIndexes[Read_Index].find_position(REFPOS, 
                                                                                    variant.rep_region.start, 
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

        if AC <= 1: #########################
            SetOfVariants[v_i] = None
        elif variant.get_variant_type() == 'SNV':
            if AC < 2 or AC <= e or AC/total < 0.05:
                SetOfVariants[v_i] = None
        
    del GA_coords

    return [v for v in SetOfVariants if v is not None]



def rough_stutter_noise(varLen, repLen):
    stutter_prob = 0.001

    if repLen >= 5:
        if abs(varLen) == 1:
            stutter_prob = 0.075
        elif abs(varLen) == 2:
            stutter_prob = 0.035

    return stutter_prob



VariantGroup = {}
for variant_type in ['INSERTION', 'DELETION', 'COMPLEX']:
    VariantGroup[variant_type] = 'INDEL'
for variant_type in ['SNV', 'INTRON']:
    VariantGroup[variant_type] = variant_type



def table_genotyping(SetOfVariants, SetOfReadNames, CUTOFFS, ploidy, REFPOS):

    all2int = {'REF': 0, 'ALT': 1, 'ALT1': 2}
    OUTPUT_LINES  = defaultdict(dict)
    ReadPairNames = set()
    New_SetOfVariants = defaultdict(dict)

    for variant in SetOfVariants:

        variant_type  = variant.get_variant_type()
        variant_group = VariantGroup[variant_type]

        FEAT_ARRAY = list()
        VAR_READ_AND_ALLELES = dict()
        
        t0 = time()

        while variant.RGPos:
            Read_Name, (Allele, RG_pos, RG_end, RG_diff_pos, RG_diff_end) = variant.RGPos.popitem() 

            Read_is_reverse = SetOfReadNames[Read_Name]['Dir']
            Read_Qualities  = SetOfReadNames[Read_Name]['Quals']
            RL = len(Read_Qualities)

            QUAL_around = np.nan
        
            if Allele == "REF":
                QUAL_around = np.nan
            elif variant_type == 'SNV':
                QUAL_around = Read_Qualities[RG_pos]
            else:
                QUAL_around = np.median(Read_Qualities[max(RG_pos - 7, 0) : RG_end + variant.rep_count + 7])
            
            READPOS = min(RL - 1 - RG_end, RG_pos)
            ALLELEINT = all2int[Allele]
            ALLELEDIR = ALLELEINT * 10 + int(Read_is_reverse)

            FEAT_ARRAY.append([ALLELEINT, QUAL_around, READPOS, ALLELEDIR, RG_diff_pos, RG_diff_end]) 
            
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

        ReadPairNames |= set(VAR_READ_AND_ALLELES)
        FEAT_ARRAY = np.array(FEAT_ARRAY)
        colNames = {'ALLELE': 0, 'DNQ': 1, 'READPOS': 2, 'STRAND': 3, 'MAP_DIFF_POS': 4, 'MAP_DIFF_END': 5}

        ### Allele counts
        alleles, counts = np.unique(FEAT_ARRAY[:, colNames['ALLELE']], return_counts = True)
        allele_counts = defaultdict(lambda : 0, zip(alleles, counts))


        AC = allele_counts[1]
        RC = allele_counts[0]
        DP = sum(allele_counts.values())
        AB = round((AC + CUTOFFS['P_maxSeqErrorRate'])/(DP + CUTOFFS['P_maxSeqErrorRate']), 3)

        ### MAP overlap
        MAP_OV1 = np.mean(FEAT_ARRAY[FEAT_ARRAY[:, 0] == 1, colNames['MAP_DIFF_POS']])
        MAP_OV2 = np.mean(FEAT_ARRAY[FEAT_ARRAY[:, 0] == 1, colNames['MAP_DIFF_END']])
        MAP_OV_mean = (MAP_OV1 + MAP_OV1)/2

        ### Read position
        READPOS_ALT_mean = np.median(FEAT_ARRAY[FEAT_ARRAY[:, 0] == 1, colNames['READPOS']])

        ### Base quality
        if variant_group == 'SNV':      
            BASEQUAL_ALT_mean = np.median(FEAT_ARRAY[FEAT_ARRAY[:, 0] == 1, colNames['DNQ']])
            if   (variant.REF, variant.ALT) in [('A', 'G')] and REFPOS.Strand in ('+', '.'):
                variant_type = 'RNA_EDITING?'
            elif (variant.REF, variant.ALT) in [('T', 'C')] and REFPOS.Strand in ('-', '.'):
                variant_type = 'RNA_EDITING?'
        
        elif variant_group == 'INDEL':
            BASEQUAL_ALT_mean = np.median(FEAT_ARRAY[FEAT_ARRAY[:, 0] == 1, colNames['DNQ']]) 
        else:            
            BASEQUAL_ALT_mean = np.nan

        ### Strand bias
        _strands, counts = np.unique(FEAT_ARRAY[:, colNames['STRAND']], return_counts = True)
        allele_directions = defaultdict(lambda : 0, zip(_strands, counts))
        FS, pFS = stats.fisher_exact([[allele_directions[0] + 1, allele_directions[1] + 1],
                                      [allele_directions[10] + 1, allele_directions[11] + 1]])

        sttuter_error = rough_stutter_noise(variant.length(), variant.rep_count) 
        error = sttuter_error + CUTOFFS['P_maxSeqErrorRate']
        pb = max(stats.binom_test(AC, DP, error), 1e-16)


        ### Filters
        if AC < CUTOFFS['AC'] or DP < CUTOFFS['DP'] or AB < CUTOFFS['AB']:
            FILTER = 'LowCov'
        elif BASEQUAL_ALT_mean < CUTOFFS['P_minBaseQual']:
            FILTER = 'LowCov'
        elif READPOS_ALT_mean < CUTOFFS['P_minReadPos']:
            FILTER = 'LowCov'
        else:
            FILTER = 'PASS'


        Feature_Vector = {'CHROM': REFPOS.chrom,
                         'POS': variant.g_start + 1,
                         'REF': variant.REF,
                         'ALT': variant.ALT, 
            	         'END': variant.g_end,
            	         'FILTER': FILTER,
                         'STRAND': REFPOS.Strand,
                         'REP_COUNT': variant.rep_count,
                         'REP_COORDS': variant.rep_region,
                         'RCL': REFPOS.RC_i,
                         'variant_group': variant_group,
            	         'variant_type': variant_type,
                         'variant_id': REFPOS.chrom + ":" + str(variant),
                         'pb': pb,
            	         'dHAP': np.nan,
            	         'ndHAP': np.nan,
            	         'AC': AC,
                         'RC': RC,
            	         'DP': DP,
            	         'AB': AB,
                         'MI': np.nan,
                         'MI_n': 0,
                         'RC_read_ov': REFPOS.Read_OV, 
            	         'INDEL_LEN': variant.length(),
                         'NVARS_ALT_mean': 1,
                         'FISHER_STRAND': FS,
            	         'READPOS_ALT_mean': READPOS_ALT_mean,
            	         'BASEQUAL_ALT_mean': BASEQUAL_ALT_mean,
            	         'MAP_OV_mean': MAP_OV_mean}

        anchor_pos    = (Feature_Vector["CHROM"], Feature_Vector["POS"], Feature_Vector["variant_group"])
        anchor_allele = (Feature_Vector["REF"], Feature_Vector["ALT"], Feature_Vector["END"])
        
        OUTPUT_LINES[anchor_pos][anchor_allele] = {"READS" : VAR_READ_AND_ALLELES, "FEAT" : Feature_Vector}
        
        

    ReadPairNames = list(ReadPairNames)

    RC_variants = []
    for anchor_pos, anchor_alleles in OUTPUT_LINES.items():
        for allele, attr in anchor_alleles.items():
            if attr["FEAT"]["variant_type"] in ('COMPLEX', 'INSERTION', 'DELETION', 'SNV'):
                RC_variants.append([anchor_pos, allele])

    RC_variants.sort(key = lambda h : h[0][1])
    n = len(RC_variants)
    m = len(ReadPairNames)


    if n >= 1 and m >= 1:
        mat = np.zeros((m, n))
        
        POS_ALLELES = defaultdict(list)

        for j, (POS, ALLELE) in enumerate(RC_variants):
            POS_ALLELES[POS].append((j, ALLELE))
            for i, r in enumerate(ReadPairNames):
                var_reads = OUTPUT_LINES[POS][ALLELE]["READS"]
                if r in var_reads:
                    mat[i][j] = 1 + int(var_reads[r] == 1)
                else:
                    mat[i][j] = 0

        mat = mat[~np.all(mat == 0, axis = 1), :]
        mat[mat == 0] = np.nan
        
        h_mean, h_std, haplotypes, clusters = _k_means_clustering(mat, ploidy)
        haplotypes = np.round(haplotypes, 0) #; print("haplotypes0", haplotypes[0], "\nhaplotypes1", haplotypes[1])

        #### Haplotype difference
        hap_mat  = np.zeros(mat.shape)
        for hap_i, read_dict in clusters.items():
            for read_i in read_dict:
                hap_mat[read_i] = haplotypes[hap_i]

        diff_mat = hap_mat - mat

        #### vars per read
        read_var_n = [0]*mat.shape[0]
        for i in range(mat.shape[0]):
            read_var_n[i] = np.nansum(mat[i] - 1.0)

        #### Remove over-ploidy alleles
        for POS, ALLELE_columns in POS_ALLELES.items():
            ACs = dict((ALLELE, np.nansum(hap_mat[:, j] - 1.0)) for j, ALLELE in ALLELE_columns)
            
            for _ in range(min(ploidy, len(ACs))):
                max_allele = max(ACs, key = ACs.get)
                ACs.pop(max_allele)

            for ALLELE, AC in ACs.items():
                OUTPUT_LINES[POS][ALLELE]['FEAT']['FILTER'] = 'LowCov' 

        #### Update Feats
        for j, (POS, ALLELE) in enumerate(RC_variants):
            hap_counts = hap_mat[ ~np.isnan(mat[:, j]), j] - 1.0
            OUTPUT_LINES[POS][ALLELE]['FEAT']['AC']    = np.nansum(hap_counts)
            OUTPUT_LINES[POS][ALLELE]['FEAT']['AB']    = np.nanmean(hap_counts)
            OUTPUT_LINES[POS][ALLELE]['FEAT']['DP']    = len(hap_counts)
            OUTPUT_LINES[POS][ALLELE]['FEAT']['RC']    = np.nansum(hap_counts == 0)
            OUTPUT_LINES[POS][ALLELE]['FEAT']['dHAP']  = np.nansum(diff_mat[:, j])
            OUTPUT_LINES[POS][ALLELE]['FEAT']['ndHAP'] = np.nanmean(diff_mat[:, j])
            OUTPUT_LINES[POS][ALLELE]['FEAT']['NVARS_ALT_mean'] = np.nansum((mat[:, j] - 1.0) * read_var_n)/np.nansum(mat[:, j] - 1.0)


        del mat
        del diff_mat

    return OUTPUT_LINES




def _k_means_clustering(mat, ploidy, iterations = 5, random_initialization = 10):
    m, n = mat.shape
    score_track = {}
    base_centroid_mean = np.nanmean(mat, axis=0)
    base_centroid_std = np.nanstd(mat, axis=0)
    for r_i in range(random_initialization):
        centroids = np.zeros((ploidy, n))
        for p_i in range(ploidy):
            random_element = np.random.random(n).reshape(1, n) - 0.5
            centroids[p_i] = base_centroid_mean + random_element * base_centroid_std

        prev_centroids = centroids.copy()
        for t_j in range(iterations):
            clusters = defaultdict(dict)
            read_dist2hap = []
            for i in range(m):
                distances = {}
                for p_i, centroid in enumerate(centroids):
                    d = np.sqrt(np.nansum(np.square(centroid - mat[i])))
                    distances[round(d, 2)] = p_i

                min_distance = min(distances)
                clusters[distances[min_distance]][i] = min_distance
                read_dist2hap.append(min_distance)

            for c_i, rows in clusters.items():
                tmp = np.nanmean(mat[list(rows), :], axis=0)
                tmp[np.isnan(tmp)] = centroids[c_i][np.isnan(tmp)]
                centroids[c_i] = tmp

            if (centroids == prev_centroids).all():
                break
            else:
                prev_centroids = centroids.copy()

        score_track[r_i] = [np.mean(read_dist2hap),
         np.std(read_dist2hap),
         centroids,
         clusters]

    best_r_i = min(score_track, key=lambda h: score_track[h][0])
    return score_track[best_r_i]



def read_genome_fasta(gf):
    fasta = pyfaidx.Faidx(gf)
    return fasta


def write_feat_file(var_list, FeatFile):
    Feat_Handle = open(FeatFile, 'w')
    Title = False
    for genomic_pos in sorted(var_list):
        for ALT, line in var_list[genomic_pos].items():

            if not Title:
                feat_list = list(line.keys())
                Feat_Handle.write('\t'.join(feat_list) + '\tLABEL\n')
                Title = True
            if abs(line['INDEL_LEN']) > 20:
                continue

            outline = [ line[feat] for feat in feat_list ] 
            Feat_Handle.write('\t'.join(map(str, outline)) + '\n')

    Feat_Handle.close()


def write_readcluster_file(ALL_READ_CLUSTERS, outfile):
    with open(outfile, 'w') as f:
        f.write("Index\tChrom\tStart\tEnd\tStrand\tMax_Coverage\n")
        for rc in ALL_READ_CLUSTERS:
            RC_Index, chrom, strand, RC_Start, RC_End, MaxCov = rc
            outline = [RC_Index, chrom, RC_Start, RC_End, strand, MaxCov]
            f.write("\t".join(map(str, outline)) + "\n")



def write_vcf_file(var_list, search_regions, genome_fasta, VcfFile, sample_name):
    genome_faidx = read_genome_fasta(genome_fasta)

    VCF_handle = open(VcfFile, 'w')
    VCF_handle.write('##fileformat=VCFv4.2\n')
    VCF_handle.write('##fileDate={}\n'.format(datetime.date.today()))
    VCF_handle.write('##CL=python {}">\n'.format(' '.join(sys.argv)))
    
    for search_region in search_regions:
        chrom = search_region[0]
        VCF_handle.write('##contig=<ID={},length={}>\n'.format(chrom, genome_faidx.index[chrom]['rlen']))

    VCF_handle.write('##FILTER=<ID=PASS,Description="Variant passes all filters">\n')
    VCF_handle.write('##FILTER=<ID=LowCov,Description="Variant does not pass one or more filtering criteria">\n')
    VCF_handle.write('##INFO=<ID=varType,Number=1,Type=String,Description="Variant type">\n')
    VCF_handle.write('##INFO=<ID=varLength,Number=1,Type=Integer,Description="Variant length">\n')
    VCF_handle.write('##INFO=<ID=varEnd,Number=1,Type=Integer,Description="End position of the variant">\n')
    VCF_handle.write('##INFO=<ID=TandemRep,Number=1,Type=Integer,Description="Count of repeat sequence adjacent to the variant">\n')
    VCF_handle.write('##INFO=<ID=MI,Number=1,Type=Float,Description="Mutual Information">\n')
    VCF_handle.write('##INFO=<ID=MI_n,Number=1,Type=Integer,Description="Count of other variants tested for MI">\n')
    VCF_handle.write('##INFO=<ID=dHAP,Number=1,Type=Float,Description="Distance from possible haplotypes">\n')
    VCF_handle.write('##INFO=<ID=ndHAP,Number=1,Type=Float,Description="Normalized dHap difference">\n')
    VCF_handle.write('##INFO=<ID=MAP_OV_mean,Number=1,Type=Float,Description="Mapping correction">\n') 
    VCF_handle.write('##INFO=<ID=NVARS_ALT_mean,Number=1,Type=Float,Description="Mean number of variants in ALT read">\n') 
    VCF_handle.write('##INFO=<ID=READPOS_ALT_mean,Number=1,Type=Float,Description="Mean read position of ALT reads">\n') 
    VCF_handle.write('##INFO=<ID=BASEQUAL_ALT_mean,Number=1,Type=Float,Description="Mean base quality of ALT reads">\n')
    VCF_handle.write('##INFO=<ID=RC_read_ov,Number=1,Type=Float,Description="Median overlap of reads within RC">\n')
    VCF_handle.write('##INFO=<ID=BinomProb,Number=1,Type=Float,Description="Binomial probability of AB">\n')
    VCF_handle.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of reads overlapping the variant position">\n')
    VCF_handle.write('##FORMAT=<ID=AC,Number=R,Type=Integer,Description="Number of reads containing the variant allele">\n')
    VCF_handle.write('##FORMAT=<ID=RC,Number=R,Type=Integer,Description="Number of reads containing the reference allele">\n')
    VCF_handle.write('##FORMAT=<ID=AB,Number=1,Type=Float,Description="Allelic balance of the variant">\n')
    VCF_handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Consensus genotype">\n')
    VCF_handle.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
    VCF_handle.write('##FORMAT=<ID=RCL,Number=1,Type=Integer,Description="Read Cluster id">\n')
    VCF_handle.write('##FORMAT=<ID=DQ,Number=1,Type=Float,Description="Change in Base Quality Score before vs After the Variant positions">\n')
    VCF_handle.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format(sample_name))
    
    for genomic_pos in sorted(var_list):
        for ALT, line_attr in var_list[genomic_pos].items():
            line = line_attr['FEAT']

            if line['variant_type'] == 'INTRON':
                continue

            outline_tmp_basic = '{CHROM}\t{POS}\t.\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t'.format(**line)

            outline_tmp_info  = 'varType={variant_type};varLength={INDEL_LEN:.0f};varEnd={END:.0f};TandemRep={REP_COUNT:.0f}'.format(**line)
            outline_tmp_info += ';MAP_OV_mean={MAP_OV_mean:.2f};dHAP={dHAP:.2f};ndHAP={ndHAP:.2f}'.format(**line)
            outline_tmp_info += ';READPOS_ALT_mean={READPOS_ALT_mean:.2f};BASEQUAL_ALT_mean={BASEQUAL_ALT_mean:.2f}'.format(**line)
            outline_tmp_info += ';MI={MI:.2f};MI_n={MI_n:.0f};BinomProb={pb:.2f};NVARS_ALT_mean={NVARS_ALT_mean:.2f};RC_read_ov={RC_read_ov:.2f}\t'.format(**line)    

            outline_tmp_format = 'GT:GQ:DP:AC:RC:AB:RCL\t'
            outline_tmp_values = './.:0:{DP:.0f}:{AC:.0f}:{RC:.0f}:{AB:.2f}:{RCL:.0f}\n'.format(**line)
            VCF_handle.write(outline_tmp_basic + outline_tmp_info + outline_tmp_format + outline_tmp_values)

    VCF_handle.close()



