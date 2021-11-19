#!/usr/bin/python
import sklearn.metrics
import HTSeq
import numpy as np ; np.seterr(all = "ignore")


def merge_introns(all_rcs_vars):
    
    intron_coords = HTSeq.GenomicArrayOfSets(chroms = "auto", stranded = False)

    sorted_introns = []
    
    for POS, ALLELE_dict in all_rcs_vars.items():
        if 'INTRON' in POS:
            for ALLELE, attr in ALLELE_dict.items():
                sorted_introns.append((POS, ALLELE))
    
    sorted_introns.sort(key = lambda x : x[0][1]) 

    sorted_introns_norep = []

    for intron in sorted_introns:
        POS, ALLELE = intron
        (chrom, start, var_group), (REF, ALT, end) = POS, ALLELE

        if end - start < 20:
            continue  	  

        similar_found = False
        minD = 6 + len(ALT)

        for j in range(1, 6):
            if j > len(sorted_introns_norep):
                break
            
            prev_POS, prev_ALLELE = sorted_introns_norep[-j]
            (prev_chrom, prev_start, prev_var_group), (prev_ref, prev_alt, prev_end) = prev_POS, prev_ALLELE

            if abs(start - prev_start) <= minD and abs(end - prev_end) <= minD:
                similar_found = True
                all_rcs_vars[prev_POS][prev_ALLELE]["READS"].update(all_rcs_vars[POS][ALLELE]["READS"])
                iGI = HTSeq.GenomicInterval(chrom, prev_start, prev_end)
                intron = sorted_introns_norep[-j]
                break
        
        if not similar_found:
            sorted_introns_norep.append(intron)
            iGI = HTSeq.GenomicInterval(chrom, start, end)
        
        intron_coords[iGI] += intron
        
    intronic_alleles = {}

    for iGI, intron_set in intron_coords.steps():
        if iGI.length > 20 and len(intron_set) >= 2:
            ip_id = (iGI.chrom, iGI.start, "INTRONIC_PART")
            ip_allele = ('', '', iGI.end)
            all_rcs_vars[ip_id][ip_allele] = {"READS" : {}, "FEAT" : {}} 

            all_rcs_vars[ip_id][ip_allele]["FEAT"]['variant_type'] = "INTRON"
            all_rcs_vars[ip_id][ip_allele]["FEAT"]['FILTER'] = "PASS"
            all_rcs_vars[ip_id][ip_allele]["FEAT"]['QUAL']   = 1e6
            all_rcs_vars[ip_id][ip_allele]["FEAT"]['AB']     = 0.5
            all_rcs_vars[ip_id][ip_allele]["FEAT"]['variant_id'] = str(iGI)
            
            intron_set = list(intron_set)

            for i, intron in enumerate(intron_set):
                i_reads = all_rcs_vars[intron[0]][intron[1]]["READS"]

                for read_name, read_allele in i_reads.items():
                    if read_allele:
                        all_rcs_vars[ip_id][ip_allele]["READS"][read_name] = i

            outline = ["[MI Log introns]", iGI.start,  iGI.end, len(intron_set), ",".join(["{0[1]}-{1[2]}".format(*_intron) for _intron in intron_set])]
            print("\t".join(map(str, outline)))

            intronic_alleles[str(iGI)] = intron_set
    return intronic_alleles




def mutual_information(all_rcs_vars, var_type_1, var_type_2, 
    mi_testable_distance = 2000, mi_testable_common_reads = 4, 
    mi_testable_mono = 2, intronic_alleles_dict = None):
    
    mi_info   = dict()


    var1_list = []
    
    for POS, ALLELE_dict in all_rcs_vars.items():
        if var_type_1 in POS:
            for ALLELE, attr in ALLELE_dict.items():
                if attr["FEAT"]["FILTER"] == "PASS" and attr["FEAT"]["QUAL"] >= 0 :
                    var1_list.append((POS, ALLELE))

    var2_list = []
    
    for POS, ALLELE_dict in all_rcs_vars.items():
        if var_type_2 in POS:
            for ALLELE, attr in ALLELE_dict.items():
                if attr["FEAT"]["FILTER"] == "PASS" and attr["FEAT"]["QUAL"] >= 0 :
                    var2_list.append((POS, ALLELE))

    var1_list.sort(key = lambda x : x[0][1])
    var2_list.sort(key = lambda x : x[0][1])


    for var1 in var1_list:
        var1_reads = all_rcs_vars[var1[0]][var1[1]]["READS"] 
        var1_id    = all_rcs_vars[var1[0]][var1[1]]["FEAT"]['variant_id']
        var1_pos   = var1[0][1]
        var1_end   = var1[1][2]
        
        MIs = []
        MI_n = []
        tag_vars = []

        for var2 in var2_list:
            var2_reads = all_rcs_vars[var2[0]][var2[1]]["READS"] 
            var2_id    = all_rcs_vars[var2[0]][var2[1]]["FEAT"]['variant_id']
            var2_pos   = var2[0][1]
            var2_end   = var2[1][2]
            

            if   var1_pos - var2_end > mi_testable_distance:
                continue
            elif var2_pos - var1_end > mi_testable_distance:
                break
           
            if var1_pos == var2_pos:
                continue
            
            common_reads = set(var1_reads) & set(var2_reads)
            lc = len(common_reads)

            if lc < mi_testable_common_reads:
                continue

            mat = np.zeros((lc, 2)).astype(int)
            for i, r in enumerate(common_reads):
                mat[i][0] = var1_reads[r]
                mat[i][1] = var2_reads[r]

            if intronic_alleles_dict is not None:
                umat, ucounts = np.unique(mat, axis = 0, return_counts = True)
                for u_i in range(umat.shape[0]):
                    jxn_coords = intronic_alleles_dict[var2_id][umat[u_i][1]]   
                             
                    jxn_s = jxn_coords[0][1]
                    jxn_e = jxn_coords[1][2]
                    jxn_mid = (jxn_e + jxn_s)//2
                    jxn_count = ucounts[u_i]
                    outline = ["[MI Log sashimi]", var1_id, var2_id, umat[u_i][0], "{}-{}".format(jxn_s, jxn_e), jxn_count]
                    print("\t".join(map(str, outline)))


            var1_alleles, var1_allele_counts = np.unique(mat[:, 0], return_counts = True)
            var2_alleles, var2_allele_counts = np.unique(mat[:, 1], return_counts = True)
            umat, ucounts = np.unique(mat, axis = 0, return_counts = True)

            s1, s2 = set(var1_reads), set(var2_reads)

            tag = "likely_heterozygous"

            for AC in var1_allele_counts:
                if AC > lc - mi_testable_mono or AC/lc > 0.9:
                    tag = "likely_homozygous"

            for AC in var2_allele_counts:
                if AC > lc - mi_testable_mono or AC/lc > 0.9:
                    tag = "likely_homozygous"

            if tag == "likely_homozygous":
                continue
           
            tag_vars.append((var2_id, var2_alleles))
            MI_score = sklearn.metrics.mutual_info_score(mat[:, 0], mat[:, 1])
            MIs.append(round(MI_score, 3))
            MI_n.append(lc)

        if MIs:
            n = len(MI_n)
            wAve_MI  = np.average(MIs, weights = MI_n)
            mi_info[var1] = [n, round(wAve_MI, 3), tag_vars, MIs, MI_n]
        else:
            mi_info[var1] = [0, np.nan, [], [], []]

    return mi_info



def mi_parse_variants(ALL_VARS_FEAT):

    printing_lines = []

    mi_info = mutual_information(ALL_VARS_FEAT, "SNV", "SNV")

    for var, (n, mi, tag_vars, mi_vars, mi_n) in mi_info.items():
        (POS, ALLELE) = var

        ALL_VARS_FEAT[POS][ALLELE]["FEAT"]["MI"]   = mi
        ALL_VARS_FEAT[POS][ALLELE]["FEAT"]["MI_n"] = n

        strand = ALL_VARS_FEAT[POS][ALLELE]["FEAT"]['STRAND']
        qual   = ALL_VARS_FEAT[POS][ALLELE]["FEAT"]['QUAL']
        var_id = ALL_VARS_FEAT[POS][ALLELE]["FEAT"]['variant_id']

        mi_vars = ["{:.3f}".format(_mi) for _mi in mi_vars]
        line = ["SNP", "SNP", var_id, qual, mi, n, ";".join([v[0] for v in tag_vars]), ";".join(mi_vars), ";".join(map(str, mi_n)), strand]
        printing_lines.append(line)



    
    mi_info = mutual_information(ALL_VARS_FEAT, "INDEL", "SNV")

    for var, (n, mi, tag_vars, mi_vars, mi_n) in mi_info.items():
        (POS, ALLELE) = var

        ALL_VARS_FEAT[POS][ALLELE]["FEAT"]["MI"]   = mi
        ALL_VARS_FEAT[POS][ALLELE]["FEAT"]["MI_n"] = n

        strand = ALL_VARS_FEAT[POS][ALLELE]["FEAT"]['STRAND']
        qual   = ALL_VARS_FEAT[POS][ALLELE]["FEAT"]['QUAL']
        var_id = ALL_VARS_FEAT[POS][ALLELE]["FEAT"]['variant_id']

        mi_vars = ["{:.3f}".format(_mi) for _mi in mi_vars]
        line = ["INDEL", "SNP", var_id, qual, mi, n, ";".join([v[0] for v in tag_vars]), ";".join(mi_vars), ";".join(map(str, mi_n)), strand]
        printing_lines.append(line)

    

    intronic_alleles_dict = merge_introns(ALL_VARS_FEAT) 

    for v_j, var1_type in enumerate(["INDEL", "SNV", "RNA_EDITING"]):

        mi_info = mutual_information(ALL_VARS_FEAT, var1_type, 
                    "INTRONIC_PART", intronic_alleles_dict = intronic_alleles_dict)

        for var, (n, mi, tag_vars, mi_vars, mi_n) in mi_info.items():
            (POS, ALLELE) = var

            ALL_VARS_FEAT[POS][ALLELE]["FEAT"]["MI"]   = mi
            ALL_VARS_FEAT[POS][ALLELE]["FEAT"]["MI_n"] = n

            strand = ALL_VARS_FEAT[POS][ALLELE]["FEAT"]['STRAND']
            qual   = ALL_VARS_FEAT[POS][ALLELE]["FEAT"]['QUAL']
            var_id = ALL_VARS_FEAT[POS][ALLELE]["FEAT"]['variant_id']
            
            for j, (ip, ip_alleles) in enumerate(tag_vars):
                _mi, lc = mi_vars[j], mi_n[j]
                ip_alleles_string = ["SP{0}:{1[1]}-{2[2]}".format(i + 1, intronic_alleles_dict[ip][i][0], intronic_alleles_dict[ip][i][1]) for i in ip_alleles]
                ip_alleles_string = ";".join(ip_alleles_string)
                line = [var1_type, "INTRONIC_PART", var_id, qual, _mi, 1, ip_alleles_string, _mi, lc, strand]
                printing_lines.append(line)


    return printing_lines



def write_mutinfo_file(mi_outlines, outfile):
    with open(outfile, 'w') as f:
        f.write("VAR1_TYPE\tVAR2_TYPE\tVAR1\tMI_mean\tMI_n\tTAG_VARS\tTAG_VARS_MIs\tTAG_VARS_LC\tSTRAND\n")
        for mi_outline in mi_outlines:
            f.write("\t".join(map(str, mi_outline)) + "\n")

    






