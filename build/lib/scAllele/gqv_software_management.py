#!/usr/bin/python
import os
import sys
import time 
import copy 

try:
	import psutil
	def print_time_stamp(Message = ""):
		process = psutil.Process()
		Memory  = process.memory_info()[0]/float(1024**3)
		timeS   = time.strftime("[%a, %d %b %Y, %H:%M:%S]")
		sys.stderr.write("{}\t{:.3f}\t{}\n".format(timeS, Memory, Message))  
		sys.stdout.flush() 
		sys.stderr.flush() 

except:
	def print_time_stamp(Message = ""):
		timeS = time.strftime("[%a, %d %b %Y %H:%M:%S]")
		sys.stderr.write("{}\t{}\n".format(timeS, Message))


def delete_reads(d1):
	for key, val in list(d1.items()):
		if key == "READS":
			d1['READS'].clear()
			
		elif isinstance(val, dict): 
			delete_reads(val)



def merge_copy(d1, d2):
    for key in d2:
        if key in d1 and isinstance(d1[key], dict) and isinstance(d2[key], dict):
            merge_copy(d1[key], d2[key])
        elif key in d1 and isinstance(d1[key], set) and isinstance(d2[key], set):
            d1[key] |= d2[key]
        else:
            d1[key] = copy.deepcopy(d2[key])

            
def rm_nested_dict(my_dict):
	my_dict.clear()
	
	# for key, sub_dict in list(my_dict.items()):
	# 	try:
	# 		rm_nested_dict(sub_dict)
	# 	except:
	# 		pass
	# 	my_dict.pop(key)


def rm_nested_list(my_list):
	while my_list:
		del my_list[0]


def get_size(obj, seen=None):
	"""Recursively finds size of objects"""

	size = sys.getsizeof(obj)
	if seen is None:
		seen = set()
	obj_id = id(obj)
	if obj_id in seen:
		return 0

	seen.add(obj_id)
	if isinstance(obj, dict):
		size += sum([get_size(v, seen) for v in obj.values()])
		size += sum([get_size(k, seen) for k in obj.keys()])
	elif hasattr(obj, '__dict__'):
		size += get_size(obj.__dict__, seen)
	elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
		size += sum([get_size(i, seen) for i in obj])

	return size










