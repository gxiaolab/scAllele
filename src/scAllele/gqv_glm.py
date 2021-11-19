#!/usr/bin/python
import numpy as np ; np.seterr(all = "ignore")
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import random
import matplotlib.pyplot as plt

from sklearn.metrics         import roc_curve, auc, f1_score
from sklearn.model_selection import ShuffleSplit
from sklearn.linear_model    import LogisticRegression
from sklearn                 import svm
from collections import Counter, defaultdict
import _pickle as pickle
import sys
import os
import re
import vcf
from vcf.parser import _Info as VcfInfo, field_counts as vcf_field_counts

from glob import glob
import pyfaidx


def covert_to_Featurevector(var_dict, model_INS, model_DEL, model_SNP, pseudo_prob = 1e-10):

	var_dict = var_dict["FEAT"]

	if np.isnan(var_dict['BASEQUAL_ALT_mean']):
		var_dict['BASEQUAL_ALT_mean'] = 40.0

	if np.isnan(var_dict["ndHAP"]):
		var_dict["ndHAP"] = 0.0

	if np.isnan(var_dict["NVARS_ALT_mean"]):
		var_dict["NVARS_ALT_mean"] = 1.0

	if np.isnan(var_dict["pb"]):
		var_dict["pb"] = 1.0

	
	FV = [np.log(var_dict["AB"] + 0.005), 
		  abs(var_dict['BASEQUAL_ALT_mean']), 
		  var_dict['REP_COUNT'], 
		  var_dict["ndHAP"], 
		  var_dict["NVARS_ALT_mean"],
		  var_dict["READPOS_ALT_mean"], 
		  np.log(var_dict["pb"])]


	try:

		if var_dict["INDEL_LEN"] == 0:
			p_False, p_True = model_SNP.predict_proba([FV])[0]
		elif var_dict["INDEL_LEN"] > 0:
			p_False, p_True = model_INS.predict_proba([FV])[0]
		elif var_dict["INDEL_LEN"] < 0:
			p_False, p_True = model_DEL.predict_proba([FV])[0]

		QUAL = np.log((p_True + pseudo_prob)/(p_False + pseudo_prob))

	except:
		print("problematic Feat Vector", FV, var_dict["DP"], var_dict["AC"], var_dict["AB"])
		QUAL = 0.0

	
	return round(QUAL, 2)



def read_benchmark(baseline_vcf):
  
	print("Loading variants...")

	var_pickle_collection = '{}.pickle'.format(baseline_vcf)

	if not os.path.exists(var_pickle_collection):
		benchmarked_vars = defaultdict(dict)

		vcf_file_handle = vcf.Reader(filename = baseline_vcf)

		vcf_file_handle.infos['datasetsmissingcall'] = VcfInfo('datasetsmissingcall', '.', 'String', \
                        "Names of datasets that are missing a call or have an incorrect call at this location, and the high-confidence call is a variant",\
                        'bar', 'foo')

		vcf_file_handle.infos['testable'] = VcfInfo('testable', '.', 'String', "Whether this loci is testableor not", 'bar', 'foo')



		LC = 0 
		for record in vcf_file_handle:
			LC += 1
			if LC % 500000 == 0:
				sys.stderr.write("\tLine counter {}...\n".format(LC))

			for ALT in record.ALT:
				var = "{}:{}:{}>{}".format(record.CHROM, record.POS, record.REF, ALT)
				VarLen = len(ALT) - len(record.REF)
				POSBIN = record.POS//int(1e5)
				try:
					benchmarked_vars[(record.CHROM, VarLen)][POSBIN].append(var)
				except:
					benchmarked_vars[(record.CHROM, VarLen)][POSBIN] = [var]

		with open(var_pickle_collection, 'wb') as fid:
			pickle.dump(benchmarked_vars, fid)

	else:
		with open(var_pickle_collection, 'rb') as fid:
			benchmarked_vars = pickle.load(fid)

	print("done")
	return benchmarked_vars


def variant_is_same(var_call, var_base):
	base_chrom, base_pos, base_ref, base_alt = re.split('[:,>]', var_base)
	call_chrom, call_pos, call_ref, call_alt = re.split('[:,>]', var_call)
	
	base_pos, call_pos = int(base_pos), int(call_pos)
	
	if abs(base_pos - call_pos) > 50:
		return False
	elif var_call == var_base:
		return True

	overhang = 50

	seq = ref_file_handle.fetch(call_chrom, base_pos - overhang, base_pos + overhang)
	seq = str(seq).upper()

	subseq = seq[overhang : overhang + len(base_ref)]

	assert subseq == base_ref.upper(), '{} != {}'.format(subseq, base_ref)
	
	ALT_BASE = seq[: overhang] + base_alt.upper() + seq[overhang + len(base_ref) : ]
	overhang = overhang + (call_pos - base_pos) 
	ALT_CALL = seq[: overhang] + call_alt.upper() + seq[overhang + len(call_ref) : ]

	# print('checking', var_call, var_base, "\n", ALT_BASE, "\n", ALT_CALL, ALT_BASE == ALT_CALL)

	if ALT_BASE == ALT_CALL:
		return True
	else:
		return False


def glm_modelling(gqv_training_file, baseline_vcf):
	df = pd.read_csv(gqv_training_file, index_col = None, sep = '\t')

	print("\tOriginal: rows=", len(df))

	df = df[(abs(df.INDEL_LEN) <= 20 )] 
	print("\tFiltering by length: rows=", len(df))
	df = df[( df.FILTER == "PASS" )] 
	print("\tFiltering by FILTER: rows=", len(df))
	df = df[( df.AC >= 2 ) & ( df.AB >= 0.05)] 
	print("\tFiltering by Cov: rows=", len(df))
	df = df[~((df.REF == "A") & (df.ALT == "G"))]
	print("\tFiltering by A-to-G: rows=", len(df))
	df = df[~((df.REF == "T") & (df.ALT == "C"))]
	print("\tFiltering by T-to-C: rows=", len(df))

	df = df.fillna(0)

	df['AB'] = np.log(df['AB'] + 0.005)
	df['AC'] = np.log(df['AC'] + 0.005) 
	df['QUAL_ALT_mean'] = df['BASEQUAL_ALT_mean']
	df['MAP_OV_mean'] = df['MAP_OV_mean'] 
	df['ndHAP'] = df['ndHAP']
	df['pb'] = np.log(df['pb'].fillna(1))


	benchmarked_vars = read_benchmark(baseline_vcf) 
	
	Y_alt = []
	var_ids = []
	for i, row in df.iterrows():
		POSBIN = row['POS']//int(1e5)
		variant_id_call = "{}:{}:{}>{}".format(row['CHRAM'], row['POS'], row['REF'], row['ALT'])
		VarLen = row['INDEL_LEN']
		label = False
		if VarLen == 0:
			if POSBIN in benchmarked_vars[(row['CHRAM'], VarLen)]:
				if variant_id_call in benchmarked_vars[(row['CHRAM'], VarLen)][POSBIN]:
					label = True
		else:
			if POSBIN in benchmarked_vars[(row['CHRAM'], VarLen)]:
				if variant_id_call in benchmarked_vars[(row['CHRAM'], VarLen)][POSBIN]:
					label = True
				else:
					for variant_id_base in benchmarked_vars[(row['CHRAM'], VarLen)][POSBIN]:
						label = variant_is_same(variant_id_call, variant_id_base)
						if label:
							break
		# print(variant_id_call , label)
		var_ids.append(variant_id_call)
		Y_alt.append(label)

	df["LABEL"]  = Y_alt
	df['VAR_id'] = var_ids

	Feats = ['AB', 'QUAL_ALT_mean', "REP_COUNT", 'ndHAP', "NVARS_ALT_mean", "READPOS_ALT_mean", "pb"]

	#Feats = ['AB', 'QUAL_ALT_mean', "REP_COUNT", "pb"]



	for IN_or_DEL, (l1, l2) in [("INSERTION", (1, 20)), ("DELETION", (-21, -1)), ("SNP", (0, 0))]:

		X = df[( df.INDEL_LEN >= l1 ) & ( df.INDEL_LEN <= l2) ] 
		X = X.reset_index(drop = True)
		Y = X.LABEL 
		X = X[Feats] 
		
		sys.stderr.write('\n\nVARIANT = {}\tLabel count: TP={} FP={}\n'.format(IN_or_DEL, sum(Y), len(Y) - sum(Y)))
		
		pos_label = sum(Y)
		neg_label = len(Y) - pos_label
		tp_ratio = pos_label/len(Y)	

		if tp_ratio > 0.66 or tp_ratio  < 0.33 or len(Y) < 1000:
			sys.stderr.write('\tWARNING: The training data set is too small = {}\n'.format(len(Y)))
			sys.stderr.write('\tWARNING: The training data set ratio is weird = {:.3f}\n'.format(tp_ratio))

			if pos_label > neg_label:
				X_neg = X[~Y]
				X_pos = X[Y].sample(n = neg_label, random_state = 1)
				X = pd.concat([X_neg, X_pos], ignore_index=True)

				Y_neg = Y[~Y]
				Y_pos = Y[Y].sample(n = neg_label, random_state = 1) 
				Y = pd.concat([Y_neg, Y_pos], ignore_index=True)

				pos_label = sum(Y)
				neg_label = len(Y) - pos_label
				tp_ratio = pos_label/len(Y)	

				sys.stderr.write('\t\tupdate: The training data set ratio is weird = {:.3f}\n'.format(tp_ratio))

		scores = list()
		fig = plt.figure()
		ax = fig.add_subplot(111)

		splits = ShuffleSplit(n_splits = 20, test_size = 0.33, random_state = 0)
		splits.get_n_splits(X, Y)

		LogRegClassifier = LogisticRegression(solver = 'liblinear')

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
		std_f1, std_auroc = np.std(scores, axis = 0)

		Title1 = 'ROC GLM classifier for micro-indels (training set = {})\n'.format(gqv_training_file)
		Title2 = 'Mean F1 = {:.2f}, Mean AUROC = {:.2f}\n'.format(mean_f1, mean_auroc)
		Title3 = 'Std. Dev. F1 = {:.3f}, Std. Dev. AUROC = {:.3f}'.format(std_f1, std_auroc)

		sys.stderr.write('\n\t' + Title2 + '\t\t' + Title3 + '\n')

		ax.set_title(Title1 + Title2 + Title3)

		fig.savefig("ROC_for_training_set={}.Var={}.jpg".format(gqv_training_file, IN_or_DEL))

		pickle_out_name = '{}.{}.glm.pickle'.format(gqv_training_file, IN_or_DEL)

		tmp_feat_vals = dict((k, round(v, 2)) for k, v  in zip(Feats, LogRegClassifier.coef_[0]))
		sys.stderr.write("\tcoefficients: {} intercept: {:.2f}\n".format(tmp_feat_vals, LogRegClassifier.intercept_[0]))
		sys.stderr.write("\n\tWriting output to {}\n".format(pickle_out_name))

		with open(pickle_out_name, 'wb') as fid:
			pickle.dump(LogRegClassifier, fid)


		FV = df[Feats]
		ps = LogRegClassifier.predict_proba(FV)
		p_False, p_True = ps[:,0], ps[:,1]
		df['QUAL_{}'.format(IN_or_DEL)] = np.log((p_True + 1e-10)/(p_False + 1e-10))	
	

	df.to_csv(gqv_training_file + ".proc2.csv", index = False, sep = "\t")




if __name__ == '__main__':

	training_file, baseline_vcf, ref_file = sys.argv[1 : 4]

	ref_file_handle = pyfaidx.Faidx(ref_file)

	glm_modelling(training_file, baseline_vcf)




