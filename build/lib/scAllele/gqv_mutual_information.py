#!/usr/bin/python
import sklearn.metrics
import HTSeq
import numpy as np ; np.seterr(all = "ignore")



def mutual_information(VAR_LIST, var_type_1, var_type_2, SM,
    mi_testable_distance = 12000, mi_testable_common_reads = 10, 
    mi_testable_mono = 0, mi_testable_ratio = 0.9, intronic_alleles_dict = None):
    
    mi_info   = dict()

    var1_list = []
    
    for POS, ALLELE_dict in VAR_LIST.items():
        if var_type_1 in POS:
            for ALLELE, attr in ALLELE_dict.items():
                if SM in attr and attr[SM]["FEAT"]["FILTER"] == "PASS" :
                    var1_list.append((POS, ALLELE))

    var2_list = []
    
    for POS, ALLELE_dict in VAR_LIST.items():
        if var_type_2 in POS:
            for ALLELE, attr in ALLELE_dict.items():
                if SM in attr and attr[SM]["FEAT"]["FILTER"] == "PASS" :
                    var2_list.append((POS, ALLELE))

    var1_list.sort(key = lambda x : x[0][1])
    var2_list.sort(key = lambda x : x[0][1])


    for var1 in var1_list:
        var1_reads = VAR_LIST[var1[0]][var1[1]][SM]["READS"] 
        var1_id    = VAR_LIST[var1[0]][var1[1]]['variant_id']
        var1_pos   = var1[0][1]
        var1_end   = var1[1][2]
        
        MIs = []
        MI_n = []
        tag_vars = []

        for var2 in var2_list:
            var2_reads = VAR_LIST[var2[0]][var2[1]][SM]["READS"] 
            var2_id    = VAR_LIST[var2[0]][var2[1]]['variant_id']
            var2_pos   = var2[0][1]
            var2_end   = var2[1][2]
                   
            common_reads = set(var1_reads) & set(var2_reads)
            lc = len(common_reads)


            if lc < mi_testable_common_reads:
                continue
            elif var1_pos - var2_end > mi_testable_distance:
                continue
            elif var2_pos - var1_end > mi_testable_distance:
                break
            elif var1_pos == var2_pos:
                continue

            mat = np.zeros((lc, 2)).astype(int)
            for i, r in enumerate(common_reads):
                mat[i][0] = var1_reads[r]
                mat[i][1] = var2_reads[r]

            if intronic_alleles_dict is not None:
                umat, ucounts = np.unique(mat, axis = 0, return_counts = True)
                for u_i in range(umat.shape[0]):
                    (jxn_s, jxn_e) = intronic_alleles_dict[var2_id][umat[u_i][1]]   
                    jxn_mid = (jxn_e + jxn_s)//2
                    jxn_count = ucounts[u_i]
                    # outline = ["[MI Log sashimi]", var1_id, var2_id, umat[u_i][0], "{}-{}".format(jxn_s, jxn_e), jxn_count]
                    # print("\t".join(map(str, outline)))


            var1_alleles, var1_allele_counts = np.unique(mat[:, 0], return_counts = True)
            var2_alleles, var2_allele_counts = np.unique(mat[:, 1], return_counts = True)

            tag = "likely_heterozygous"

            for AC in var1_allele_counts:
                if AC > lc - mi_testable_mono or AC/lc > mi_testable_ratio:
                    tag = "likely_homozygous"

            for AC in var2_allele_counts:
                if AC > lc - mi_testable_mono or AC/lc > mi_testable_ratio:
                    tag = "likely_homozygous"

            if tag == 'likely_homozygous':
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




def mi_parse_variants(VAR_LIST, SM, options):

    printing_lines = []

    mi_info = mutual_information(VAR_LIST, "SNV", "SNV", SM, 
                mi_testable_common_reads = options.link_min_count,
                mi_testable_mono  = 2, 
                mi_testable_ratio = 0.9)

    for var, (n, mi, tag_vars, mi_vars, mi_n) in mi_info.items():
        (POS, ALLELE) = var

        VAR_LIST[POS][ALLELE][SM]["FEAT"]["MI"]   = mi
        VAR_LIST[POS][ALLELE][SM]["FEAT"]["MI_n"] = n

        if n:
            strand = VAR_LIST[POS][ALLELE]['STRAND']
            qual   = VAR_LIST[POS][ALLELE][SM]["FEAT"]['QUAL']
            var_id = VAR_LIST[POS][ALLELE]['variant_id']

            if mi < 0.3 and n >= 2:
                VAR_LIST[POS][ALLELE]['variant_type'] = "RNA_EDITING"
        
            mi_vars = ["{:.3f}".format(_mi) for _mi in mi_vars]
            
            line = [SM, "SNV", "SNV", var_id, qual, mi, n, ";".join([v[0] for v in tag_vars]), 
                        ";".join(mi_vars), ";".join(map(str, mi_n)), strand]
            
            printing_lines.append(line)


    
    mi_info = mutual_information(VAR_LIST, "INDEL", "SNV", SM,
                mi_testable_common_reads = options.link_min_count,
                mi_testable_mono  = 2, 
                mi_testable_ratio = 0.9)

    for var, (n, mi, tag_vars, mi_vars, mi_n) in mi_info.items():
        (POS, ALLELE) = var

        VAR_LIST[POS][ALLELE][SM]["FEAT"]["MI"]   = mi
        VAR_LIST[POS][ALLELE][SM]["FEAT"]["MI_n"] = n

        if n:
            strand = VAR_LIST[POS][ALLELE]['STRAND']
            qual   = VAR_LIST[POS][ALLELE][SM]["FEAT"]['QUAL']
            var_id = VAR_LIST[POS][ALLELE]['variant_id']

            mi_vars = ["{:.3f}".format(_mi) for _mi in mi_vars]
            
            line = [SM, "INDEL", "SNP", var_id, qual, mi, n, ";".join([v[0] for v in tag_vars]), 
                        ";".join(mi_vars), ";".join(map(str, mi_n)), strand]
            
            printing_lines.append(line)

    
    IP_alleles_dict = {}

    for POS, ALLELE_dict in VAR_LIST.items():
        if "INTRONIC_PART" in POS:
            for ALLELE, attr in ALLELE_dict.items():
                if SM in attr:
                    ip_introns = [(intron[0][1], intron[1][2]) for intron in attr[SM]["IP_ALLELES"]]
                    ip_var_id  = attr["variant_id"]
                    IP_alleles_dict[ip_var_id] = ip_introns


    for var1_type in ["INDEL", "SNV", "RNA_EDITING"]:
        
        mi_info = mutual_information(VAR_LIST, var1_type, "INTRONIC_PART", SM, 
                                        mi_testable_common_reads = options.link_min_count,
                                        mi_testable_mono  = 0, 
                                        mi_testable_ratio = 1.0,
                                        intronic_alleles_dict = IP_alleles_dict)


        for var, (n, mi, tag_vars, mi_vars, mi_n) in mi_info.items():
            (POS, ALLELE) = var

            VAR_LIST[POS][ALLELE][SM]["FEAT"]["MI"]   = mi
            VAR_LIST[POS][ALLELE][SM]["FEAT"]["MI_n"] = n

            strand = VAR_LIST[POS][ALLELE]['STRAND']
            qual   = VAR_LIST[POS][ALLELE][SM]["FEAT"]['QUAL']
            var_id = VAR_LIST[POS][ALLELE]['variant_id']
            
            for j, (ip, ip_alleles) in enumerate(tag_vars):
                mi, lc = mi_vars[j], mi_n[j]
                
                ip_alleles_string = ["SP{0}:{1[0]}-{1[1]}".format(i + 1, IP_alleles_dict[ip][i]) \
                                        for i in ip_alleles]

                ip_alleles_string = ip + ";".join(ip_alleles_string)
                
                line = [SM, var1_type, "INTRONIC_PART", var_id, qual, mi, 1, ip_alleles_string, 
                         mi, lc, strand]
                
                printing_lines.append(line)

    return printing_lines









