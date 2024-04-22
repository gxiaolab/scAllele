#!/usr/bin/python
import time
import sys
import gc
import os
import re
import copy
import HTSeq
import collections
import multiprocessing 
import pickle as pickle
import numpy as np
from optparse import OptionParser, OptionGroup

import importlib.resources as pkg_resources

from . import glm
from . import gqv_utils
from . import gqv_dbg_utils
from . import gqv_bam_utils
from . import gqv_glm
from . import gqv_vartool
from . import gqv_software_management as SOFTMAN
from . import gqv_mutual_information  as MUTINFO



__version__ = "0.0.9.4"


def process_read_cluster_1(ARGs):
	
	RC_info, options, bam_file, SM = ARGs

	RC_Index, chrom, strand, RC_Start, RC_End, MaxCov = RC_info 

	genome_faidx = gqv_utils.read_genome_fasta(options.Genome)

	mp_id = multiprocessing.current_process().name
	
	### Obtain reads that overlap RC
	RC_vals = gqv_bam_utils.get_RC_reads(bam_file, RC_info, genome_faidx, SM, options)

	if RC_vals is None:
		return None
	else:
		SetOfReadNames, REFPOS, SetOfReadGroups = RC_vals
	
	### Construct de Bruijn graph
	db = gqv_dbg_utils.construct_de_Bruijn_graph(SetOfReadGroups, options) 
	
	### Process de Bruijn graph
	xdb, SetOfReadIndexes = gqv_dbg_utils.Process_de_Bruijn_Graph(SetOfReadGroups, options, db) 

	del db

	REF_INDEX = SetOfReadGroups[0]['Index'] 


	### DFT algorithm for src-snk pairs compressed
	Src_Snk_pairs = gqv_dbg_utils.find_Source_Target_pairs(xdb, REF_INDEX, SetOfReadIndexes, options)  
	del xdb
	
	if not Src_Snk_pairs:
		return None

	### Find variants for every read index
	SetOfVariants = gqv_utils.assign_vars_to_reads(Src_Snk_pairs, SetOfReadIndexes, 
													SetOfReadGroups, options, REFPOS) 
	
	SOFTMAN.rm_nested_dict(Src_Snk_pairs)
	
	### Process, filter and modify variants
	SetOfVariants = gqv_vartool.find_equivalent_indels(SetOfVariants, REFPOS) 

	### Obtian all overlapping read indexes for every variant
	SetOfVariants = gqv_utils.overlap_vars_and_reads(SetOfVariants, SetOfReadIndexes, 
										   SetOfReadGroups, REFPOS)
	SOFTMAN.rm_nested_dict(SetOfReadGroups)
	SOFTMAN.rm_nested_list(SetOfReadIndexes) 

	### Feature collection for every variant
	VAR_LIST = gqv_utils.feature_collection(SetOfVariants, SetOfReadNames, options, SM, REFPOS)

	### Genotyping table
	VAR_LIST = gqv_utils.table_genotyping(VAR_LIST, SM, options)	
	
	return bam_file, SM, RC_info, SetOfVariants, VAR_LIST




def process_read_cluster_2(ARGs):

	RC_info, options, bam_file, SM, SetOfVariants, AllOfVariants = ARGs

	RC_Index, chrom, strand, RC_Start, RC_End, MaxCov = RC_info 

	# ----- get variants -------
	RefVariants = []
	RC_GI = HTSeq.GenomicInterval(chrom, RC_Start, RC_End)
	for v_gi, variant_set in AllOfVariants[RC_GI].steps():
		for variant in variant_set:
			if variant not in SetOfVariants:
				RefVariants.append(variant)

	# --------------------------------------------

	if not RefVariants:
		return {}

	genome_faidx = gqv_utils.read_genome_fasta(options.Genome)

	mp_id = multiprocessing.current_process().name

	### Obtain reads that overlap RC
	RC_vals = gqv_bam_utils.get_RC_reads(bam_file, RC_info, genome_faidx, SM, options)

	if RC_vals is None:
		return {}
	else:
		SetOfReadNames, REFPOS, SetOfReadGroups = RC_vals

	### Obtian all overlapping read indexes for every variant
	SetOfVariants = gqv_utils.overlap_vars_and_reads_ref(RefVariants, SetOfReadGroups)

	SOFTMAN.rm_nested_dict(SetOfReadGroups)

	### Feature collection for every variant
	VAR_LIST = gqv_utils.feature_collection_ref(SetOfVariants, SetOfReadNames, options, SM, REFPOS)	

	return VAR_LIST




def parse_chroms(bam_files, options):
	region = options.search_region

	genome_chroms = gqv_utils.read_genome_fasta(options.Genome)

	bam_chroms = set()
	for bam_file in bam_files:
		bam_chroms |= set(gqv_bam_utils.get_chroms_from_bam(bam_file))

	common_chroms = set(genome_chroms.index) & bam_chroms
	
	Region_List = []
	
	r = re.match("([\d|\w]+):(\d+)\-(\d+)", region)

	if not common_chroms:
		raise ValueError("Invalid genome reference\n")

	elif region == "All":
		# all super-clusters
		Region_List = gqv_bam_utils.find_super_read_clusters(bam_files, options) 
	
	elif bool(r):
		chrom, fetch_start, fetch_end = r.groups()

		# all super-clusters that overlap this region
		if chrom not in common_chroms:
			raise ValueError(f"Region is not valid {region}\n")

		Region_List = gqv_bam_utils.find_super_read_clusters(bam_files, options, 
																chrom   = chrom, 
																f_start = int(fetch_start), 
																f_end   = int(fetch_end)) 
	
	else:
		# all super-clusters in this chrom
		if region in common_chroms: 		
			Region_List = gqv_bam_utils.find_super_read_clusters(bam_files, options, 
																	chrom = region)
		else:
			raise ValueError(f"Region is not valid {region}\n")

	return Region_List




def main():
	args = sys.argv

	description = """A variant caller and variant analysis tool for scRNA-seq data."""

	usage = """\n\tscAllele -b <file.bam> -g <genome.fa> -o <output prefix>"""

	strandedness_options = ['fr-firststrand', 'fr-secondstrand', 'fr-unstrand']


	parser = OptionParser(usage = usage, description = description)
	
	parser.add_option("-b", "--input-bam",   
		dest = "input_bam", 
		help = "[Required] Input bam file, (or comma-seprated list of bam files) sorted and indexed")
	parser.add_option("-o", "--output-vcf",  
		dest = "output_prefix", 
		help = "[Required] Prefix of the output files")
	parser.add_option("-g", "--genome-file", 
		dest = "Genome", 
		help = "[Required] Reference genome file (fasta format)")


	filtering_options = OptionGroup(parser, "Filtering Options")

	filtering_options.add_option("--AB", 
		dest = "minRatioVar", 
		type = "float",  
		default = 1e-2, 
		help = "Minimum allelic ratio for the variant allele. Default: 0.01")
	filtering_options.add_option("--AC", 
		dest = "minCountVar", 
		type = "int",    
		default = 2,    
		help = "Minimum read depth supporting the variant allele. Default: 2")
	filtering_options.add_option("--DP", 
		dest = "minCoverage", 
		type = "int",    
		default = 5,   
		help = "Minimum read depth at the position of the variant. Default: 5")
	filtering_options.add_option("--min-base_position", 
		dest = "minReadPos", 
		type = "int", 
		default = 10,
		help = "Minimum mean distance for the variant from the read ends. Default = 7")
	filtering_options.add_option("--min-base_quality", 
		dest = "minBaseQual", 
		type = "int", 
		default = 20,
		help = "Minimum mean base quality score to consider SNPs. Default = 20")

	parser.add_option_group(filtering_options)

	
	run_mode_options = OptionGroup(parser, "Run Mode Options")
	
	run_mode_options.add_option("--run_mode",  
		dest = "run_mode", 
		default = 'Full', 
		help = "Select <Variant_Caller> <Full> or <Training> mode. Default: Full")
	run_mode_options.add_option("--glm_clf_name",  
		dest = "glm_classifier_name", 
		default = None, 
		help = "Prefix of the GLM pickle objects with the GLM models")
	
	parser.add_option_group(run_mode_options)


	linkage_options = OptionGroup(parser, "Linkage Options")
	
	linkage_options.add_option("--link_min_count",  dest = "link_min_count", 
		default =  10,
		type = "int", 
		help = "Minimum number of common reads for linkage analysis. Default = 10")
	linkage_options.add_option("--link_min_mi",  dest = "link_min_mi", 
		default = 0.52, 
		type = "float",
		help = "Minimum mutual information for linkage analysis. Default = 0.52")
	
	parser.add_option_group(linkage_options)


	advanced_options = OptionGroup(parser, "Advanced Options")

	advanced_options.add_option("-c", "--region", 
		dest = "search_region",     
		default = 'All', 
		help = "Limit search to this region (chrom:start-end or chrom). Default: All")
	advanced_options.add_option("-n", "--nodes", 
		dest = "nodes", 
		type = "int", 
		default = 64, 
		help = "Number of threads for parallel execution. Default = 64")
	advanced_options.add_option("--min-map_quality", 
		dest = "minMapQ", 
		type = "int", 
		default = 40,
		help = "Minimum mapping quality to keep a read. Default = 40")
	advanced_options.add_option("--max_soft_clip", 
		dest = "maxSoftClip", 
		type = "int", 
		default = 5,
		help = "Maximum length of soft-clipping allow at each end of the read. Default = 5")
	advanced_options.add_option("--kmer-size", 
		dest = "kmer", 
		type = "int", 
		default = 15,
		help = "k-mer size for the de-Bruijn graph assembly. Default: 15")
	advanced_options.add_option("--strandedness", 
		dest = "strandedness", 
		default = "fr-unstrand",
		help = "Select from {}. Default: 'fr-unstrand'".format(strandedness_options))
	advanced_options.add_option("--maxSeqErrorRate", 
		dest = "maxSeqErrorRate", 
		type = "float", 
		default = 1e-2,
		help = "Maximum estimate of sequencing error rate. Default: 0.01")
	advanced_options.add_option("--Ploidy", 
		dest = "ploidy", 
		type = "int", 
		default = 2,
		help = "Maximum ploidy to be considered. Default: 2.")

	parser.add_option_group(advanced_options)

	(options, Args) = parser.parse_args()

	if len(sys.argv[1:]) < 3:
		parser.print_help()
		sys.exit(1)


	if not os.path.exists(options.Genome):
		parser.error('Genome file does not exist')

	if not options.output_prefix:   
		parser.error('Output file prefix not given')
	
	if options.strandedness not in strandedness_options:
		parser.error("\nStrandedness option: {} not valid. \
			Use the -h option to view usage\n".format(strandedness))

	if options.kmer > 20 or options.kmer < 10:
		SOFTMAN.print_time_stamp("WARNING: We recomend a k-mer size of 10-20")

	if options.run_mode in ("Full", "Variant_Caller"):
		if options.glm_classifier_name is not None:
			f_d = open(options.glm_classifier_name + '.DELETION.glm.pickle', 'rb')
			f_i = open(options.glm_classifier_name + '.INSERTION.glm.pickle', 'rb')
			f_s = open(options.glm_classifier_name + '.SNP.glm.pickle', 'rb')
		else:
			f_d = pkg_resources.open_binary(glm, 'GM12878_smartseq2.all.feature_matrix.tab.DELETION.glm.pickle')
			f_i = pkg_resources.open_binary(glm, 'GM12878_smartseq2.all.feature_matrix.tab.INSERTION.glm.pickle')
			f_s = pkg_resources.open_binary(glm, 'GM12878_smartseq2.all.feature_matrix.tab.SNP.glm.pickle')

		GLM_DEL = pickle.load(f_d)
		GLM_INS = pickle.load(f_i)
		GLM_SNP = pickle.load(f_s)

	if options.run_mode not in ("Variant_Caller", "Full", "Training"):
		parser.error("Run mode {} not valid\n".format(options.run_mode))


	

	sys.setrecursionlimit(int(1e5))

	SOFTMAN.print_time_stamp("COMMAND = '{}'".format(' '.join(args)))
	
	SOFTMAN.print_time_stamp("Running {} on {} mode".format(args[0], options.run_mode))


	bam_list = gqv_bam_utils.read_bam_input(options.input_bam)

	Search_Regions = parse_chroms(bam_list, options)

	SAMPLE_list = [(bam_file, gqv_bam_utils.get_sample_name(bam_file)) for bam_file in bam_list]


	if options.run_mode in ("Full", "Variant_Caller"):
		gqv_utils.write_vcf_header(Search_Regions, options.Genome, 
									f"{options.output_prefix}.vcf", SAMPLE_list)
		gqv_utils.write_introns_header(f"{options.output_prefix}.intronic_parts.bed")

	if options.run_mode in ("Full"):
		gqv_utils.write_mutinfo_header(f"{options.output_prefix}.mi_summary.tab")

	os.system(f"rm -f {options.output_prefix}.*")


	n = min(options.nodes, multiprocessing.cpu_count())
	pool = multiprocessing.Pool(n)


	for Search_Region in Search_Regions:	

		(chrom, fetch_start, fetch_end) = Search_Region
		L = fetch_end - fetch_start
		SOFTMAN.print_time_stamp("Processing contig {0[0]}:{0[1]}-{0[2]} Length = {1}".format(Search_Region, L))
			
		ALL_READ_CLUSTERS = collections.defaultdict(dict)
		chrom_VAR_LIST    = collections.defaultdict(dict)
		chrom_VAR_MUTINFO = []
		
		ARGs = []
		for (bam_file, SM) in SAMPLE_list:	
			RCs = gqv_bam_utils.find_read_clusters(bam_file, options, chrom, fetch_start, fetch_end)
			ARGs += [(rc, options, bam_file, SM) for rc in RCs]

		
		### ----- Pool processing -----

		raw_VAR_LIST = []

		for rc_output in pool.imap_unordered(process_read_cluster_1, ARGs): 
			if rc_output is None: 
				continue
			
			bam_file, SM, rc, SetOfVariants, rc_VAR_LIST = rc_output
			raw_VAR_LIST.extend(copy.deepcopy(SetOfVariants))
			SOFTMAN.merge_copy(chrom_VAR_LIST, rc_VAR_LIST)
			ALL_READ_CLUSTERS[(SM, bam_file)][tuple(rc)] = SetOfVariants				

		gc.collect()

		
		for (bam_file, SM) in SAMPLE_list:

			if options.run_mode in ("Full", "Variant_Caller"):
				chrom_VAR_LIST = gqv_vartool.merge_introns(chrom_VAR_LIST, SM)
				gqv_glm.update_vars(GLM_INS, GLM_DEL, GLM_SNP, chrom_VAR_LIST, SM, options)

			if options.run_mode in ("Full", ):
				mi_outlines = MUTINFO.mi_parse_variants(chrom_VAR_LIST, SM, options)
				chrom_VAR_MUTINFO.extend(mi_outlines)

		raw_VAR_LIST = gqv_utils.process_var_list(chrom, raw_VAR_LIST, chrom_VAR_LIST)

		
		ARGs = []
		for (SM, bam_file), RCs in ALL_READ_CLUSTERS.items():

			if not RCs: 
				continue
			
			for rc, SetOfVariants in RCs.items():
				ARGs.append((rc, options, bam_file, SM, SetOfVariants, raw_VAR_LIST))

		for rc_VAR_LIST in pool.imap_unordered(process_read_cluster_2, ARGs): 
			SOFTMAN.merge_copy(chrom_VAR_LIST, rc_VAR_LIST)


		
		if options.run_mode in ("Full", "Variant_Caller"):

			gqv_utils.write_introns_bed(chrom_VAR_LIST,  f"{options.output_prefix}.intronic_parts.bed")
			chrom_VAR_LIST = gqv_utils.merge_SM_VAR_LIST(chrom_VAR_LIST)
			chrom_VAR_LIST = gqv_vartool.merge_introns(chrom_VAR_LIST, 'merged_SM')
			gqv_glm.update_vars(GLM_INS, GLM_DEL, GLM_SNP, chrom_VAR_LIST, 'merged_SM', options)

		if options.run_mode in ("Full",):

			mi_outlines = MUTINFO.mi_parse_variants(chrom_VAR_LIST, 'merged_SM', options)
			chrom_VAR_MUTINFO.extend(mi_outlines)


		SOFTMAN.delete_reads(chrom_VAR_LIST)


		if options.run_mode in ("Full", "Variant_Caller"):
			gqv_utils.write_vcf_file(chrom_VAR_LIST,    f"{options.output_prefix}.vcf", SAMPLE_list)
			gqv_utils.write_introns_bed(chrom_VAR_LIST, f"{options.output_prefix}.intronic_parts.bed")

		if options.run_mode in ("Full",):
			gqv_utils.write_mutinfo_file(chrom_VAR_MUTINFO, f"{options.output_prefix}.mi_summary.tab")

		if options.run_mode == "Training" and chrom_VAR_LIST:
			gqv_utils.write_feat_file(chrom_VAR_LIST, f"{options.output_prefix}.feature_matrix.tab")


		gqv_utils.write_readcluster_file(ALL_READ_CLUSTERS, f"{options.output_prefix}.read_cluster_info.tab")


		SOFTMAN.rm_nested_dict(chrom_VAR_LIST)
		SOFTMAN.rm_nested_dict(ALL_READ_CLUSTERS)
		gc.collect()


	pool.close()
	pool.terminate()
	SOFTMAN.print_time_stamp("DONE!")



if __name__ == "__main__":
	main()





