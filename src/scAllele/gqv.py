#!/usr/bin/python

import sys
import multiprocessing 
import gc
import _pickle as pickle
import os
import re
import collections
from optparse import OptionParser, OptionGroup
from time import time

import importlib.resources as pkg_resources
from . import glm

from . import gqv_utils
from . import gqv_dbg_utils
from . import gqv_bam_utils
from . import gqv_glm
from . import gqv_software_management as gqv_SftMngr
from . import gqv_vartool
from . import gqv_mutual_information  as gqv_MutInfo





def process_read_cluster(args):

	RC_info, CUTOFFS, options, k = args

	RC_Index, chrom, strand, RC_Start, RC_End, MaxCov = RC_info

	if MaxCov < CUTOFFS["DP"]:
		return {}

	genome_faidx = gqv_utils.read_genome_fasta(options.Genome)

	mp_id = multiprocessing.current_process().name

	RC_vals = gqv_bam_utils.get_RC_reads(options.input_bam, RC_info, genome_faidx, 
										 options.strandedness, CUTOFFS)

	if RC_vals is None:
		return {}
	else:
		SetOfReadNames, REFPOS, SetOfReadGroups, gqv_reads_total, gqv_reads_lost = RC_vals

	
	### Construct de Bruijn graph
	db = gqv_dbg_utils.construct_de_Bruijn_graph(SetOfReadGroups, k)

	### Process de Bruijn graph
	xdb, SetOfReadIndexes = gqv_utils.Process_de_Bruijn_Graph(SetOfReadGroups, k, db, CUTOFFS) 

	del db

	REF_INDEX = SetOfReadGroups[0]['Index']

	### DFT algorithm for src-snk pairs compressed
	Src_Snk_pairs = gqv_dbg_utils.find_Source_Target_pairs(xdb, REF_INDEX, SetOfReadIndexes, k) 

	del xdb

	### Find variants for every read index
	SetOfVariants = gqv_utils.assign_vars_to_reads(Src_Snk_pairs, REF_INDEX, SetOfReadIndexes, 
													SetOfReadGroups, k, REFPOS) 
	gqv_SftMngr.rm_nested_dict(Src_Snk_pairs)
	
	### Process, filter and modify variants
	SetOfVariants = gqv_vartool.find_equivalent_indels(SetOfVariants, REFPOS) 

	### Obtian all overlapping read indexes for every variant
	SetOfVariants = gqv_utils.overlap_vars_and_reads(SetOfVariants, SetOfReadIndexes, 
													   SetOfReadGroups, REFPOS)

	gqv_SftMngr.rm_nested_dict(SetOfReadGroups)
	gqv_SftMngr.rm_nested_list(SetOfReadIndexes) 
	
	### Feature collection for every variant
	VAR_LIST = gqv_utils.table_genotyping(SetOfVariants, SetOfReadNames, 
											CUTOFFS, options.ploidy, REFPOS) 

	if RC_Index % 100 == 0:
		gc.collect()
	
	return VAR_LIST




def parse_chroms(fasta_file, bam_file, region):
	genome_chroms = gqv_utils.read_genome_fasta(fasta_file)
	bam_chroms    = gqv_bam_utils.get_chroms_from_bam(bam_file)

	common_chroms = set(genome_chroms.index) & set(bam_chroms)

	Region_List = []
	
	if not common_chroms:
		raise ValueError("There are no common chromosomes between the bam file and the genome reference. \
							Please check the Genome Assembly used in both files\n")
	elif region is None:
		Region_List = [(chrom, None, None) for chrom in common_chroms]
	elif bool(re.match("([\d|\w])+:(\d)+\-(\d)+", region)):
		chrom, fetch_start, fetch_end = re.split('[:,-]', region)
		if chrom not in common_chroms:
			raise ValueError("Region selection {} not valid/not present in bam and/or reference\n".format(region))
		Region_List = [(chrom, int(fetch_start), int(fetch_end))]
	elif os.path.exists(region):
		with open(region) as rf:
			for line in rf:
				try:
					chrom, start, end = line.split()
					start, end = int(start), int(end)
					Region_List.append((chrom, start, end))
				except:
					chrom = line.strip()
					Region_List.append((chrom, None, None))
				if chrom not in common_chroms:
					raise ValueError("Region selection {} not valid/not present in bam and/or reference\n".format(line))
	else:
		if region in common_chroms: 
			Region_List = [(region, None, None)]
		else:
			raise ValueError("Region selection {} not valid/not present in bam and/or reference\n".format(region))
	return Region_List



def update_vars(GLM_ins, GLM_del, GLM_snp, ALL_VARS_FEAT):
	for POS, ALT_dict in ALL_VARS_FEAT.items():
		for ALLELE, feat_dict in ALT_dict.items():
			if "INTRON" in POS or "INTRONIC_PART" in POS:
				ALL_VARS_FEAT[POS][ALLELE]["FEAT"]["QUAL"] = 1e6
			else:
				qual = gqv_glm.covert_to_Featurevector(feat_dict, GLM_ins, GLM_del, GLM_snp, pseudo_prob=1e-10)
				ALL_VARS_FEAT[POS][ALLELE]["FEAT"]["QUAL"] = qual




if __name__ == "__main__":
	main()


def main():
	args = sys.argv

	gqv_SftMngr.print_time_stamp("COMMAND = '{}'".format(' '.join(args)))

	description = """A variant caller and variant analysis for scRNA-seq data"""

	usage = """\n\tscAllele -b <file.bam> -g <genome.fa> -o <output prefix>"""

	strandedness_options = ['fr-firststrand', 'fr-secondstrand', 'fr-unstrand']


	parser = OptionParser(usage = usage, description = description)
	
	parser.add_option("-b", "--input-bam",   dest = "input_bam", 
		help = "[Required] Input bam file, sorted and indexed")
	parser.add_option("-o", "--output-vcf",  dest = "output_vcf", 
		help = "[Required] Prefix of the output files")
	parser.add_option("-g", "--genome-file", dest = "Genome", 
		help = "[Required] Reference Genome (Fasta format)")


	filtering_options = OptionGroup(parser, "Filtering Options")

	filtering_options.add_option("--AB", dest = "minRatioVar", type = "float",  
		default = 1e-2, 
		help = "Minimum allelic ratio for the variant allele. Default: 0.01")
	filtering_options.add_option("--AC", dest = "minCountVar", type = "int",    
		default = 2,    
		help = "Minimum read depth supporting the variant allele. Default: 2")
	filtering_options.add_option("--DP", dest = "minCoverage", type = "int",    
		default = 5,   
		help = "Minimum read depth at the position of the variant. Default: 10")
	filtering_options.add_option("--min-base_position", dest = "minReadPos", type = "int", 
		default = 10,
		help = "Minimum mean distance for the variant from the read ends. Default = 7")
	filtering_options.add_option("--min-base_quality", dest = "minQual", type = "int", 
		default = 20,
		help = "Minimum mean base quality score to consider SNPs. Default = 20")

	parser.add_option_group(filtering_options)

	
	run_mode_options = OptionGroup(parser, "Run Mode Options")
	
	run_mode_options.add_option("--run_mode",  dest = "run_mode", 
		default =  'Variant_Caller', 
		help = "Select <Variant_Caller> or <Training> mode. Default: Variant_Caller")
	run_mode_options.add_option("--glm_clf_name",  dest = "glm_classifier_name", 
		default = None, 
		help = "Prefix of the GLM pickle objects with the GLM models")
	
	parser.add_option_group(run_mode_options)


	advanced_options = OptionGroup(parser, "Advanced Options")

	advanced_options.add_option("-c", "--region", dest = "search_region",     
		default = None, 
		help = "Limit search to this region (chrom:start-end or chrom). Default: All")
	advanced_options.add_option("-n", "--nodes", dest = "nodes", type = "int", 
		default = 64, 
		help = "Number of threads for parallel execution. Default = 16")
	advanced_options.add_option("--min-map_quality", dest = "minMapQ", type = "int", 
		default = 40,
		help = "Minimum mapping quality to keep a read. Default = 40")
	advanced_options.add_option("--max_soft_clip", dest = "maxSoftClip", type = "int", 
		default = 5,
		help = "Maximum length of soft-clipping allow at each end of the read. Default = 5")
	advanced_options.add_option("--kmer-size", dest = "kmer", type = "int", 
		default = 15,
		help = "k-mer size for de-Bruijn graph assembly. Default: 15")
	advanced_options.add_option("--strandedness", dest = "strandedness", 
		default = "fr-unstrand",
		help = "Select from {}. Default: 'fr-unstrand'".format(strandedness_options))
	advanced_options.add_option("--maxSeqErrorRate", dest = "maxSeqErrorRate", type = "float", 
		default = 1e-2,
		help = "Maximum estimate of sequencing error rate. Default: 0.01")
	advanced_options.add_option("--Ploidy", dest = "ploidy", type = "int", 
		default = 2,
		help = "Maximum ploidy to be considered. Default: 2.")

	parser.add_option_group(advanced_options)

	(options, Args) = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		print(description)
		parser.print_usage()
		sys.exit()

	if not os.path.exists(options.input_bam):
		parser.error('Bam file does not exist')

	if not os.path.exists(options.Genome):
		parser.error('Genome file does not exist')

	if not options.output_vcf:   
		parser.error('Output file prefix not given')
	
	if options.strandedness not in strandedness_options:
		parser.error("\nStrandedness option: {} not valid. \
			Use the -h option to view usage\n".format(strandedness))

	if options.kmer > 20 or options.kmer < 10:
		gqv_SftMngr.print_time_stamp("WARNING: We recomend a k-mer size of 10-20")

	if options.run_mode == "Variant_Caller":
		if options.glm_classifier_name is not None:

			f_d = open(options.glm_classifier_name + '.DELETION.glm.pickle', 'rb')
			f_i = open(options.glm_classifier_name + '.INSERTION.glm.pickle', 'rb')
			f_s = open(options.glm_classifier_name + '.SNP.glm.pickle', 'rb')
		else:

			f_d = pkg_resources.open_binary(glm, 'smartseq2_gm12878_DS.txt.DELETION.glm.pickle')
			f_i = pkg_resources.open_binary(glm, 'smartseq2_gm12878_DS.txt.INSERTION.glm.pickle')
			f_s = pkg_resources.open_binary(glm, 'smartseq2_gm12878_DS.txt.SNP.glm.pickle')

		GLM_DEL = pickle.load(f_d)
		GLM_INS = pickle.load(f_i)
		GLM_SNP = pickle.load(f_s)

	if options.run_mode not in ("Variant_Caller", "Training"):
		parser.error("Run mode {} not valid\n".format(options.run_mode))



	k = options.kmer

	CUTOFFS = {	"P_minMapQual"  : options.minMapQ,
				"P_minBaseQual" : options.minQual,
				"P_minReadPos"  : options.minReadPos,
				"P_minSoftClipLen" : options.maxSoftClip,
				"P_maxSeqErrorRate" : options.maxSeqErrorRate,
				"AC" : options.minCountVar,
				"AB" : options.minRatioVar,
				"DP" : options.minCoverage}

	sys.setrecursionlimit(int(1e5))
	
	gqv_SftMngr.print_time_stamp("Running {} on {} mode".format(args[0], options.run_mode))

	Search_Regions = parse_chroms(options.Genome, options.input_bam, options.search_region)

	ALL_VARS_FEAT = collections.defaultdict(dict)
	ALL_VARS_MUTINFO = []
	ALL_READ_CLUSTERS = []

	for (chrom, fetch_start, fetch_end) in Search_Regions:
			
		RCs = gqv_bam_utils.find_read_clusters(options.input_bam, options.strandedness, 
												CUTOFFS, chrom, fetch_start, fetch_end)
		ALL_READ_CLUSTERS.extend(RCs)

		CHROM_VARS_FEAT = collections.defaultdict(dict)

		if RCs:
			gqv_SftMngr.print_time_stamp("\tContig {} : {} read clusters found.".format(chrom, len(RCs)))


		if options.nodes > 1:

			n = min(options.nodes, multiprocessing.cpu_count())

			gqv_SftMngr.print_time_stamp("\t\tRunning {} child processes".format(n))

			pool = multiprocessing.Pool(n)

			RCs = [(rc, CUTOFFS, options, k) for rc in RCs]

			gen_output = pool.imap_unordered(process_read_cluster, RCs)

			for key in gen_output: 
				var_list_rc = key
				gqv_SftMngr.merge_copy(CHROM_VARS_FEAT, var_list_rc)

			pool.terminate()
			
		elif options.nodes == 1: 
			for arg in RCs:
				var_list_rc = process_read_cluster(arg)
				gqv_SftMngr.merge_copy(CHROM_VARS_FEAT, var_list_rc)


		if options.run_mode == "Variant_Caller":

			gqv_SftMngr.print_time_stamp("\tCalculating variant scores")

			update_vars(GLM_INS, GLM_DEL, GLM_SNP, CHROM_VARS_FEAT)

			gqv_SftMngr.print_time_stamp("\tCalculating mutual information")

			mi_outlines = gqv_MutInfo.mi_parse_variants(CHROM_VARS_FEAT)

			ALL_VARS_MUTINFO.extend(mi_outlines)


		for POS, ALLELE_dict in CHROM_VARS_FEAT.items():
			for ALLELE in ALLELE_dict:
				CHROM_VARS_FEAT[POS][ALLELE].pop("READS")

		gqv_SftMngr.merge_copy(ALL_VARS_FEAT, CHROM_VARS_FEAT)

		gc.collect()
		
	gqv_SftMngr.print_time_stamp("\tWritting output to prefix {}...".format(options.output_vcf))


	
	if options.run_mode == "Variant_Caller":
		SM = gqv_bam_utils.get_sample_name(options.input_bam)
		gqv_utils.write_vcf_file(ALL_VARS_FEAT,	Search_Regions, options.Genome,
										options.output_vcf + ".vcf", SM)

		gqv_MutInfo.write_mutinfo_file(ALL_VARS_MUTINFO, 
										options.output_vcf + ".mi_summary.tab")

	elif options.run_mode == "Training":
		gqv_utils.write_feat_file(ALL_VARS_FEAT, 
										options.output_vcf + ".feature_matrix.tab")

	
	gqv_utils.write_readcluster_file(ALL_READ_CLUSTERS, 
										options.output_vcf + ".read_cluster_info.tab")

	gqv_SftMngr.print_time_stamp("DONE!")






