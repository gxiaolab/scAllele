#!/usr/bin/python

import numpy as np ; np.seterr(all = "ignore")
import pandas as pd
import matplotlib ; matplotlib.use('Agg')
import random
import pickle
import pyfaidx
import matplotlib.pyplot as plt

from glob import glob
from optparse import OptionParser, OptionGroup
from collections import Counter, defaultdict

from sklearn.metrics         import roc_curve, auc, f1_score
from sklearn.model_selection import ShuffleSplit, GridSearchCV
from sklearn.linear_model    import LogisticRegression
from sklearn.utils           import resample

import sys
import os
import re
import vcfpy





def update_vars(GLM_ins, GLM_del, GLM_snp, ALL_VARS_FEAT, SM, options):
	for POS, ALT_dict in ALL_VARS_FEAT.items():
		if "INTRON" in POS or "INTRONIC_PART" in POS:
			continue

		for ALLELE, var_dict in ALT_dict.items():
			if SM not in var_dict:
				continue
			else:
				qual = covert_to_Featurevector(var_dict, SM, GLM_ins, GLM_del, GLM_snp, options, )

			sm_feat = var_dict[SM]["FEAT"]
			
			ALL_VARS_FEAT[POS][ALLELE][SM]["FEAT"]["QUAL"] = qual

			### Filters
			if sm_feat['OVER_PLOIDY']:
				FILTER = 'Overploidy'
			elif sm_feat['AC'] < options.minCountVar:
				FILTER = 'LowAC'
			elif sm_feat['DP'] < options.minCoverage:
				FILTER = 'LowDP'
			elif sm_feat['AB'] < options.minRatioVar:
				FILTER = 'LowAB'
			elif sm_feat['BASEQUAL_ALT_mean'] < options.minBaseQual:
				FILTER = 'LowBaseQual'
			elif sm_feat['READPOS_ALT_mean']  < options.minReadPos:
				FILTER = 'LowReadPos'
			
			elif qual >= 3.01:
				FILTER = "PASS"
			else:
				FILTER = "LowQual"

	
			ALL_VARS_FEAT[POS][ALLELE][SM]["FEAT"]["FILTER"] = FILTER
		
			


def covert_to_Featurevector(var_dict, SM, model_INS, model_DEL, model_SNP, 
	options, pseudo_prob = 1e-10):

	sm_feat = var_dict[SM]["FEAT"]

	FV = [sm_feat["AB"], 
		  sm_feat['BASEQUAL_ALT_mean'], 
		  var_dict['REP_COUNT'], 
		  sm_feat["ndHAP"], 
		  sm_feat["NVARS"],
		  sm_feat["READPOS_ALT_mean"], 
		  np.log(sm_feat["pb"] + 1e-5 )]


	try:

		if var_dict["INDEL_LEN"] == 0:
			p_False, p_True = model_SNP.predict_proba([FV])[0]
		elif var_dict["INDEL_LEN"] > 0:
			p_False, p_True = model_INS.predict_proba([FV])[0]
		elif var_dict["INDEL_LEN"] < 0:
			p_False, p_True = model_DEL.predict_proba([FV])[0]

		QUAL = -10 * np.log10(p_False + pseudo_prob)

	except:
		QUAL = 0.0

	return round(QUAL, 2)



def variant_is_same(var_call, var_base, ref_genome, overhang = 25):
	base_chrom, base_pos, base_ref, base_alt = re.split('[:,>]', var_base)
	call_chrom, call_pos, call_ref, call_alt = re.split('[:,>]', var_call)
	
	base_pos, call_pos = int(base_pos), int(call_pos)
	
	if var_call == var_base:
		return True
	elif abs(base_pos - call_pos) > overhang:
		return False
	elif len(base_ref) * len(base_alt) * len(call_ref) * len(call_alt) == 1:
		return False
	elif (len(base_ref) - len(base_alt)) != (len(call_ref) - len(call_alt)):
		return False
	

	seq = ref_genome.fetch(call_chrom, base_pos - overhang, base_pos + overhang)
	seq = str(seq).upper()

	subseq = seq[overhang : overhang + len(base_ref)]

	assert subseq == base_ref.upper(), '{} != {}'.format(subseq, base_ref)
	
	ALT_BASE = seq[: overhang] + base_alt.upper() + seq[overhang + len(base_ref) : ]
	overhang = overhang + (call_pos - base_pos) 
	ALT_CALL = seq[: overhang] + call_alt.upper() + seq[overhang + len(call_ref) : ]

	if ALT_BASE == ALT_CALL:
		return True
	else:
		return False



def process_dataframe(gqv_training_file, baseline_vcf, ref_genome, overhang = 25):
	canonical_chroms = ["chr1","chr2", "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
						"chr11","chr12", "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
						"chr21","chr22", "chrX"]

	df = pd.read_csv(gqv_training_file, index_col = None, sep = '\t')

	print("\tOriginal: rows=", len(df))

	df = df[(abs(df.INDEL_LEN) <= 20 )] 
	print("\tFiltering by length: rows=", len(df))

	df = df[(df.AC >= 2) & ( df.AB >= 0.05)] 
	print("\tFiltering by Cov: rows=", len(df))
	
	df = df[(df.BASEQUAL_ALT_mean >= 20)]
	print("\tFiltering by BASEQUAL: rows=", len(df))

	df = df[(df.READPOS_ALT_mean >= 7)]
	print("\tFiltering by READPOS: rows=", len(df))

	df = df[(df.variant_type != "RNA_EDITING?")]
	print("\tFiltering by editing type rows=", len(df))

	df = df[(df.DP <= 10)]
	print("\tFiltering by coverage =", len(df))


	df = df.fillna(0)
	
	df['pb'] = np.log(df['pb'] + 1e-5)

	vcf_file_handle = vcfpy.Reader(filename = baseline_vcf)

	# vcf_file_handle.infos['datasetsmissingcall'] = VcfInfo('datasetsmissingcall', '.', 'String', "bar", 'bar', 'bar')
	# vcf_file_handle.infos['testable']            = VcfInfo('testable', '.', 'String', "bar", 'bar', 'foo')
	
	Y_alt = []
	var_ids = []

	for i, row in df.iterrows():
		CHROM, start, end, REF, ALT = re.split('[-,:,>]', row['variant_id'])
		start, end = int(start) + 1, int(end)

		var_id = f'{CHROM}:{start}:{REF}>{ALT}'

		found = False 
		label = False

		if CHROM in canonical_chroms:

			for record in vcf_file_handle.fetch(CHROM, start - overhang, end + overhang):
				for record_ALT in record.ALT:
					var_id_base = f"{CHROM}:{record.POS}:{record.REF}>{record_ALT.value}"
					label = variant_is_same(var_id, var_id_base, ref_genome)
					# print("checking", row["REP_COUNT"], row["variant_type"], var_id, var_id_base, label )
					if label: 
						found = True

				if found:
					break

		if i % 50000 == 0:
			sys.stderr.write("\tLine counter {}...\n".format(i))
			sys.stderr.flush()

		var_ids.append(var_id)
		Y_alt.append(int(label))


	df["LABEL"]  = Y_alt

	return df 




def glm_modelling(gqv_training_file, df):
	#

	Feats = ['AB', 'BASEQUAL_ALT_mean', "REP_COUNT", 'ndHAP', "NVARS", "READPOS_ALT_mean", "pb"]

	for VAR_TYPE, (l1, l2) in [("INSERTION", (1, 20)), ("DELETION", (-21, -1)), ("SNP", (0, 0))]:

		X = df[( df.INDEL_LEN >= l1 ) & ( df.INDEL_LEN <= l2) ] 
		X = X.reset_index(drop = True)
		Y = X.LABEL 
		X = X[Feats] 

		sys.stderr.write('\n\nVARIANT = {}\tLabel count: TP={} FP={}\n'.format(VAR_TYPE, sum(Y), len(Y) - sum(Y)))

		scores = list()
		fig = plt.figure()
		ax = fig.add_subplot(111)

		splits = ShuffleSplit(n_splits = 20, test_size = 0.33, random_state = 0)
		splits.get_n_splits(X, Y)

		LogRegClassifier = LogisticRegression(solver = 'liblinear', class_weight = 'balanced')

		for (train, test) in splits.split(X):
			X_train, Y_train = X.loc[train], Y[train]
			X_test,  Y_test  = X.loc[test],  Y[test]

			probas = LogRegClassifier.fit(X_train, Y_train).predict_proba(X_test)
			fpr, tpr, thresholds = roc_curve(Y_test, probas[:, 1])
			roc_auc = auc(fpr, tpr)

			Y_pred = LogRegClassifier.fit(X_train, Y_train).predict(X_test)
			f1_s = f1_score(Y_test, Y_pred)

			scores.append((f1_s, roc_auc))

			ax.plot(fpr, tpr, lw = 0.6)

		ax.plot([0, 1], [0, 1], transform = ax.transAxes, color='r', ls='--')
		ax.set_xlim([-0.05, 1.05])
		ax.set_ylim([-0.05, 1.05])
		ax.set_xlabel('False Positive Rate')
		ax.set_ylabel('True Positive Rate')

		mean_f1, mean_auroc = np.mean(scores, axis = 0)
		std_f1, std_auroc   = np.std(scores, axis = 0)

		Title1 = 'ROC GLM classifier for micro-indels (training set = {})\n'.format(gqv_training_file)
		Title2 = 'Mean F1 = {:.2f}, Mean AUROC = {:.2f}\n'.format(mean_f1, mean_auroc)
		Title3 = 'Std. Dev. F1 = {:.3f}, Std. Dev. AUROC = {:.3f}'.format(std_f1, std_auroc)

		sys.stderr.write('\n\t' + Title2 + '\t\t' + Title3 + '\n')

		ax.set_title(Title1 + Title2 + Title3)

		fig.savefig("{}.ROC_for_training_set.Var={}.pdf".format(gqv_training_file, VAR_TYPE))

		pickle_out_name = '{}.{}.glm.pickle'.format(gqv_training_file, VAR_TYPE)

		tmp_feat_vals = dict((k, round(v, 2)) for k, v  in zip(Feats, LogRegClassifier.coef_[0]))

		sys.stderr.write("\tcoefficients: {} intercept: {:.2f}\n".format(tmp_feat_vals, LogRegClassifier.intercept_[0]))
		sys.stderr.write("\n\tWriting output to {}\n".format(pickle_out_name))

		with open(pickle_out_name, 'wb') as fid:
			pickle.dump(LogRegClassifier, fid)

		ps = LogRegClassifier.predict_proba(df[Feats])
		p_False, p_True = ps[:,0], ps[:,1]
		df['QUAL_{}'.format(VAR_TYPE)] = -10 * np.log10(p_False + 1e-10)	
	
	df.to_csv(gqv_training_file + ".proc.csv", index = False, sep = "\t")



def main():
	args = sys.argv

	description  = """scAllele's training module:\n"""
	description += """With a feature matrix file (run scAllele in "Training" mode)\n"""
	description += """and a groud truth set of variants, train a new variant classifier.\n"""

	usage = """\n\tscAllele_train -i <Feat_Matrix.tab> -v <truth.vcf> -g <genome.fa> """

	parser = OptionParser(usage = usage, description = description)
	
	parser.add_option("-i", 
		dest = "featmat", 
		help = "[Required] scAllele's Feature Matrix file")
	parser.add_option("-v",     
		dest = "truth_vcf", 
		help = "[Required] Groundtruth set of variants (vcf format)")
	parser.add_option("-g",  
		dest = "Genome", 
		help = "[Required] Reference Genome (fasta format)")

	(options, Args) = parser.parse_args()

	if len(sys.argv) <= 4:
		parser.print_help()
		sys.exit(1)

	ref_file = pyfaidx.Faidx(options.Genome)

	df = process_dataframe(options.featmat, options.truth_vcf, ref_file)

	glm_modelling(options.featmat,  df)



if __name__ == '__main__':
	main()

