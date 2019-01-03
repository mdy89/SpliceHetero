import os
import sys
import commands
import numpy as np
import time
import math
import multiprocessing as mp
import itertools
import argparse
import random
from scipy.stats import entropy
import datetime

nn = "\n"
tt = "\t"
ss = "/"
cc = ","

def GET_TIME():
	return str(datetime.datetime.now()).split(".")[0]
	##End GET_TIME

def main():
	#[1]Argument Parsing
	print "STEP_1: Argument Parsing... [%s]"%GET_TIME()
	ArgumentParser = GenArgumentParser()

	#[2]Input Processing
	print "STEP_2: Input Processing... [%s]"%GET_TIME()
	InputHandler = GenInputHandler(ArgumentParser)

	#[3]sITH Calculation
	print "STEP_3: sITH Calculation... [%s]"%GET_TIME()
	UnitHandler = GenUnitHandler(ArgumentParser, InputHandler)

	#[4] 
	print "Completed [%s]"%GET_TIME()
	##End main

class GenArgumentParser(object):
	def __init__(self):
		#[1]
		self.Parser_obj = argparse.ArgumentParser()
		self.BuildParser_method()

		#[2]
		self.Handler_obj = self.Parser_obj.parse_args()
		self.CheckHandler_method()

		#[3]
		self.Out_dir = self.Handler_obj.out_directory_path; makePath(self.Out_dir)
		self.ArgumentExport_method()
		##End init
	
	def ArgumentExport_method(self):
		#[1]
		Argument_path = self.Out_dir + ss + "ARGUMENTS.txt"
		Argument_file = open(Argument_path, 'w')
		COMMAND_line = "CMD: %s"%makeLine(sys.argv, " ")
		Argument_file.write(COMMAND_line + nn)

		#[2]
		ARG_index = 1
		for ARG_NAME in sorted(self.Handler_obj.__dict__):
			ARG_VALUE = self.Handler_obj.__dict__[ARG_NAME]
			ARG_line = "[%s]%s=%s"%(ARG_index, ARG_NAME, ARG_VALUE); ARG_index += 1
			Argument_file.write(ARG_line + nn)
			##End for

		#[3]
		Argument_file.close()
		##End ArgumentExport_method
	
	def CheckHandler_method(self):
		#[1]
		self.Message_list = []
		
		#[2]
		self.Message_list.append("Usage: python SpliceHetero.py [arguments]")
		self.Message_list.append("\nEssemtial Arguments:")
		self.Message_list.append("'-csp', '--Case_Sample_path' [PATH]\t\tList of case-sample input (see web-page for more input format)")
		self.Message_list.append("'-rsp', '--Reference_Sample_path' [PATH]\t\tList of reference-sample input (see web-page for more input format)")
		self.Message_list.append("'-odp', '--Out_Dir_path' [PATH]\t\tDirectory path for output'")
		
		#[3]
		self.Message_list.append("\nOptional Arguments:")
		self.Message_list.append("'-prn', '--Process_number' [INT]\t\tNumber of processes to use (default: 1)")
		self.Message_list.append("'-slb', '--Stranded_Library_bool' [BOOL]\t\tIf it is stranded library (options: True/False, default: False)'")

		#[4]
		self.Message_line = makeLine(self.Message_list, nn)
		self.Essential_list = ["case_sample_path", "reference_sample_path", "out_directory_path"]
		for arg_name in self.Essential_list:
			if self.Handler_obj.__dict__[arg_name] == None:
				sys.exit(self.Message_line)
				##End if
			##End for
		##End CheckHandler_method
	
	def BuildParser_method(self):
		#[1]Essential
		self.Parser_obj.add_argument('-csp', '--Case_Sample_path', dest="case_sample_path",help = "")
		self.Parser_obj.add_argument('-rsp', '--Reference_Sample_path', dest="reference_sample_path",help = "")
		self.Parser_obj.add_argument('-odp', '--Out_Dir_path', dest="out_directory_path",help = "")

		#[2]Optional
		self.Parser_obj.add_argument('-prn', '--Process_number', dest="process_number_int",help = "", default=1)
		self.Parser_obj.add_argument('-slb', '--Stranded_Library_bool', dest="Stranded_Library_bool",help = "", default="False")
		##End BuildParser_method
	##End GenArgumentParser

class GenInputHandler(object):
	def __init__(self, ArgumentParser):
		#[1]
		self.Handler_obj = ArgumentParser.Handler_obj

		#[2]
		self.SampleRegister_method()
		##End init
	
	def SampleRegister_method(self):
		#[1]
		self.Case_list = []; self.Case_idx = 0
		self.Case_path = self.Handler_obj.case_sample_path
		for case_line in open(self.Case_path):
			sample_name, junction_path = case_line.split()[:2]
			sampleObj = GenSampleObject()
			sampleObj.Name_str = sample_name
			sampleObj.Junction_path = junction_path
			self.Case_list.append(sampleObj)
			##End for
		
		#[2]
		self.Ref_list = []; self.Ref_idx = 0
		self.Ref_path = self.Handler_obj.reference_sample_path
		for ref_line in open(self.Ref_path):
			sample_name, junction_path = ref_line.split()[:2]
			sampleObj = GenSampleObject()
			sampleObj.Name_str = sample_name
			sampleObj.Junction_path = junction_path
			self.Ref_list.append(sampleObj)
			##End for

		#[3]
		a, b = map(len, [self.Case_list, self.Ref_list])
		print " - CASE:%s, REF:%s samples are registered [%s]"%(a, b, GET_TIME())
		##End SampleRegister_method
	##End GenInputHandler

class GenUnitHandler(object):
	def __init__(self, ArgumentParser, InputHandler):
		#[1]
		self.Handler_obj = ArgumentParser.Handler_obj
		self.Ref_list = InputHandler.Ref_list
		self.Case_list = InputHandler.Case_list
		self.Sample_list = self.Ref_list + self.Case_list

		#[2]
		self.nProcess_int = int(self.Handler_obj.process_number_int)
		self.sLibrary_bool = self.Handler_obj.Stranded_Library_bool.upper() == "TRUE"
		self.Out_dir = self.Handler_obj.out_directory_path
		self.Unit_dir = self.Out_dir + ss + "tmp"
		makePath(self.Out_dir); makePath(self.Unit_dir)

		#[3]
		print "STEP_3_1: Process Managing... [%s]"%GET_TIME()
		self.HashExtract_method()

		#[4]
		self.UnitGenerate_method()

		#[5]
		self.UnitExport_method()

		#[6]
		self.SHD_Export_method()
		##End init
	
	def SHD_Export_method(self):
		#[1]
		self.Normalized_path = self.Out_dir + ss + "sITH_REGIONS.txt"
		self.Normalized_file = open(self.Normalized_path, 'r')
		self.SampleName_list = self.Normalized_file.readline().split()[4:]

		#[2]
		self.SHD_dic = {}
		for normalized_line in self.Normalized_file:
			sample_distance_list = map(float, normalized_line.split()[4:])
			for sample_index, sample_distance in zip(range(len(self.SampleName_list)), sample_distance_list):
				if sample_index not in self.SHD_dic:
					self.SHD_dic[sample_index] = []
					##End if
				self.SHD_dic[sample_index].append(sample_distance)
				##End for
			##End for

		#[3]
		self.SHD_path = self.Out_dir + ss + "sITH.txt"
		self.SHD_file = open(self.SHD_path, 'w')
		head_line = makeLine(["Sample", "sITH"], tt)
		self.SHD_file.write(head_line + nn)

		#[4]
		for sample_index in sorted(self.SHD_dic):
			sample_name = self.SampleName_list[sample_index]
			SHD_list = filter(lambda x:not math.isnan(x), self.SHD_dic[sample_index])
			if SHD_list:
				SHD_AVG = np.mean(SHD_list)
			else:
				SHD_AVG = float("nan")
				##End if-else
			SHD_line = makeLine([sample_name, SHD_AVG], tt)
			self.SHD_file.write(SHD_line + nn)
			##End for

		#[5]
		self.SHD_file.close()

		#[6]
		commands.getoutput("rm -r %s"%self.Unit_dir)
		##End SHD_Export_method
	
	def UnitExport_method(self):
		#[1]
		self.Count_path = self.Out_dir + ss + "COUNT.txt"
		self.Distance_path = self.Out_dir + ss + "sITH_REGIONS.txt"

		#[2]
		self.Count_file = open(self.Count_path, 'w')
		self.Distance_file = open(self.Distance_path, 'w')

		#[3]
		self.Hash_list = []; self.Hash_index = 0
		self.Hash_path = self.Unit_dir + ss + "HashList.txt"
		for hash_line in open(self.Hash_path):
			chr_name, hash_size = hash_line.split()
			hashObj = GenHashObject(); hashObj.Chr_name = chr_name
			hashObj.Job_name = "Hash.%s"%self.Hash_index; self.Hash_index += 1
			hashObj.Size_int = int(hash_size)
			self.Hash_list.append(hashObj)
			##End for

		#[4]
		self.Hash_index = 0
		for jobObj in sorted(self.Hash_list, key=lambda x:x.Size_int, reverse=True):
			jobObj.Index_int = self.Hash_index; self.Hash_index += 1
			jobObj.Count_path = self.Unit_dir + ss + "%s.COUNT_REGIONS.txt"%jobObj.Job_name
			jobObj.Distance_path = self.Unit_dir + ss + "%s.sITH_REGIONS.txt"%jobObj.Job_name
			self.JobExport_method(jobObj)
			##End for

		#[5]
		self.Count_file.close()
		self.Distance_file.close()
		##End UnitExport_method
		
	def JobExport_method(self, jobObj):
		#[1]
		jobObj.Count_file = open(jobObj.Count_path, 'r')
		jobObj.Distance_file = open(jobObj.Distance_path, 'r')

		#[2]
		if jobObj.Index_int == 0:
			self.Count_file.write(jobObj.Count_file.readline())
			self.Distance_file.write(jobObj.Distance_file.readline())
			head_line = makeLine(["Chromosome", "Shared_site", "Fixed_sites"], tt)
		else:
			jobObj.Count_file.readline()
			jobObj.Distance_file.readline()
			##End if
			
		#[3]
		self.Count_file.write(jobObj.Count_file.read())
		self.Distance_file.write(jobObj.Distance_file.read())

		#[4]
		jobObj.Count_file.close()
		jobObj.Distance_file.close()
		##End JobExport_method
	
	def UnitGenerate_method(self):
		#[1]
		self.Process_list = []
		self.Semaphore_obj = mp.Semaphore(self.nProcess_int)
		self.Hash_path = self.Unit_dir + ss + "HashList.txt"
		self.Hash_line_list = open(self.Hash_path).readlines()
		self.Hash_index_list = range(len(self.Hash_line_list))
		self.HASH_flag_list = map(lambda x:int(np.percentile(self.Hash_index_list,x)), range(0,101,10))

		#[2]
		self.Hash_index = 0
		for hash_line in open(self.Hash_path):
			chr_name, hash_size = hash_line.split()
			jobObj = GenHashObject(); jobObj.Chr_name = chr_name
			jobObj.Job_name = "Hash.%s"%self.Hash_index; self.Hash_index += 1
			self.Semaphore_obj.acquire(); argVect = tuple([jobObj])
			procObj = mp.Process(target=self.JobProcess_method, args=argVect)
			procObj.start(); self.Process_list.append(procObj)
			if (self.Hash_index-1) in self.HASH_flag_list:
				a, b = self.Hash_index, len(self.Hash_line_list)
				print " - %s/%s jobs are being processed... [%s]"%(a, b, GET_TIME())
				##End if
			##End for

		#[3]
		for procObj in self.Process_list:
			procObj.join()
			##End for
		##End UnitGenerate_method

	def JobProcess_method(self, jobObj):
		#[1]
		self.Chr_name = jobObj.Chr_name
		self.Job_name = jobObj.Job_name

		#[2]
		self.UnitExtract_method()

		#[3]
		self.CountExtract_method()

		#[4]
		self.DistanceExtract_method()

		#[5]
		self.Semaphore_obj.release()
		##End JobProcess_method
	
	def DistanceExtract_method(self):
		#[1]
		self.Count_path = self.Unit_dir + ss + "%s.COUNT_REGIONS.txt"%self.Job_name
		self.Count_file = open(self.Count_path, 'r')
		self.Count_file.readline()

		#[2]
		self.Distance_path = self.Unit_dir + ss + "%s.sITH_REGIONS.txt"%self.Job_name
		self.Distance_file = open(self.Distance_path, 'w')
		self.Name_list = map(lambda x:x.Name_str, self.Case_list)
		self.Head_line = makeLine(["Chromosome", "Shared_site", "Alternative_sites", "Strand"] + self.Name_list, tt)
		self.Distance_file.write(self.Head_line + nn)

		#[3]
		for count_line in self.Count_file:
			distance_line = self.Get_DistanceLine_method(count_line)
			self.Distance_file.write(distance_line + nn)
			##End for	

		#[4]
		self.Distance_file.close()
		##End DistanceExtract_method

	def Get_DistanceLine_method(self, count_line):
		#[1]
		chr_name, shared_site, alternative_sites, strand = count_line.split()[:4]
		count_line_list = count_line.split()[4:]
		count_vector_list = map(lambda x:map(float,x.split(cc)), count_line_list)

		#[2]
		profile_vector_list = []
		for count_vector in count_vector_list:
			#[2-1]
			count_sum = float(sum(count_vector))
			if count_sum == 0:
				profile_vector = [float("nan")] * len(count_vector)
			else:
				profile_vector = self.Get_ProfileVector_method(count_vector)
				##End if-else
			#[2-2]
			profile_vector_list.append(profile_vector)
			##End for
			
		#[3]
		avg_case_distance_dic = {}
		avg_case_distance_list = []
		for caseObj in self.Case_list:
			#[3-1]
			case_idx = caseObj.Index_int
			case_profile = profile_vector_list[case_idx]
			case_distance_list = []
			#[3-2]
			case_profile_key = tuple(case_profile)
			if case_profile_key in avg_case_distance_dic:
				avg_case_distance = avg_case_distance_dic[case_profile_key]
				avg_case_distance_list.append(avg_case_distance); continue
				##End if
			#[3-3]
			case_ref_distance_dic = {}
			for refObj in self.Ref_list:
				#[3-3-1]
				ref_idx = refObj.Index_int
				ref_profile = profile_vector_list[ref_idx]
				#[3-3-2]
				ref_profile_key = tuple(ref_profile)
				if ref_profile_key in case_ref_distance_dic:
					case_distance = case_ref_distance_dic[ref_profile_key]
					case_distance_list.append(case_distance); continue
					##End if
				#[3-3-3]
				case_distance = self.GetDistance_method(ref_profile, case_profile)	
				case_distance_list.append(case_distance)
				case_ref_distance_dic[ref_profile_key] = case_distance
				##End for
			#[3-4]
			case_distance_list = filter(lambda x:not math.isnan(x), case_distance_list)
			if case_distance_list:
				avg_case_distance = np.mean(case_distance_list)
			else:
				avg_case_distance = float("nan")
				##End if-else
			#[3-5]
			avg_case_distance_list.append(avg_case_distance)
			avg_case_distance_dic[case_profile_key] = avg_case_distance
			##End for

		#[4]
		distance_line = makeLine([chr_name, shared_site, alternative_sites, strand] + avg_case_distance_list, tt)
		return distance_line
		##End Get_ProfileLine_method

	def Get_ProfileVector_method(self, count_vector):
		#[1]
		count_sum = sum(count_vector); pseudo_count = count_sum/100.
		adjusted_count_vector = map(lambda x:x+pseudo_count, count_vector)
		adjusted_count_sum = sum(adjusted_count_vector)

		#[2]
		adjusted_profile_vector = map(lambda x:x/adjusted_count_sum, adjusted_count_vector)
		return adjusted_profile_vector
		##End Get_ProfileVector_method

	def GetDistance_method(self, profile_A, profile_B):
		#[1]
		if profile_A == profile_B:
			return 0.
			##End if

		#[2]
		profile_M = map(lambda x:(x[0]+x[1])/2, zip(profile_A, profile_B))

		#[3]
		profile_M_sum = sum(profile_M)
		if math.isnan(profile_M_sum):
			return float("nan")
			##End if

		#[4]
		distance_A = entropy(profile_A, profile_M, base=2)
		distance_B = entropy(profile_B, profile_M, base=2)

		#[5]
		return (distance_A + distance_B)/2
		##End GetDistance_method
	
	def CountExtract_method(self):
		#[1]
		self.CountRegister_method()

		#[2]
		self.CountExport_method()
		##End CountExtract_method
	
	def CountExport_method(self):
		#[1]
		self.Count_path = self.Unit_dir + ss + "%s.COUNT_REGIONS.txt"%self.Job_name
		self.Count_file = open(self.Count_path, 'w')
		self.Name_list = map(lambda x:x.Name_str, self.Sample_list)
		self.Head_line = makeLine(["Chromosome", "Shared_site[1-based]", "Alternative_sites[1-based]", "Strand"] + self.Name_list, tt)
		self.Count_file.write(self.Head_line + nn)

		#[2]
		case_name_set = set(map(lambda x:x.Name_str, self.Case_list))
		ref_name_set = set(map(lambda x:x.Name_str, self.Ref_list))
		for unitObj in sorted(self.Unit_list, key=lambda x:x.Fixed_site):
			#[2-1]
			sample_index_list = sorted(unitObj.Count_dic.keys())
			sample_name_set = set(map(lambda x:self.Sample_list[x].Name_str, sample_index_list))
			#[2-2]
			if not (sample_name_set & case_name_set):
				continue
				##End if	
			if not (sample_name_set & ref_name_set):
				continue
				##End if	
			#[2-3]
			unitObj.Var_sites = sorted(unitObj.Var_sites)
			count_line_list = []
			for sampleObj in self.Sample_list:
				var_count_list = []
				for var_site in unitObj.Var_sites:
					try: var_count = unitObj.Count_dic[sampleObj.Index_int][var_site]
					except: var_count = 0
					var_count_list.append(var_count)
					##End for
				var_count_line = makeLine(var_count_list, cc)
				count_line_list.append(var_count_line)
				##End for
			#[2-4]
			var_sites_line = makeLine(map(lambda x:x[0], unitObj.Var_sites), cc)
			unit_count_line = makeLine([unitObj.Chr_name, unitObj.Fixed_site[0], var_sites_line, unitObj.Strand_str] + count_line_list, tt)
			self.Count_file.write(unit_count_line + nn)
			##End for

		#[3]
		self.Count_file.close()
		##End CountExport_method
	
	def CountRegister_method(self):
		for raw_line in open(self.Raw_path):
			#[1]
			chr_name, junc_st, junc_ed, junc_str, sample_idx, junc_cnt = raw_line.split()
			junc_st, junc_ed, sample_idx = map(int, [junc_st, junc_ed, sample_idx])
			junc_cnt = float(junc_cnt)
			#[2]
			junc_st_key = int(junc_st), junc_str; junc_ed_key = int(junc_ed), junc_str
			#[3]
			if junc_st_key in self.Unit_dic:
				unitObj = self.Unit_dic[junc_st_key]
				if sample_idx not in unitObj.Count_dic:
					unitObj.Count_dic[sample_idx] = {}
					##End if
				unitObj.Count_dic[sample_idx][junc_ed_key] = junc_cnt
				##End if
			#[4]
			if junc_ed_key in self.Unit_dic:
				unitObj = self.Unit_dic[junc_ed_key]
				if sample_idx not in unitObj.Count_dic:
					unitObj.Count_dic[sample_idx] = {}
					##End if
				unitObj.Count_dic[sample_idx][junc_st_key] = junc_cnt
				##End if
			##End fo
		##End CountRegister_method
	
	def UnitExtract_method(self):
		#[1]
		self.Left_dic = {}; self.Right_dic = {}

		#[2]
		self.Raw_path = self.Unit_dir + ss + "%s.RAW.txt"%self.Job_name
		self.Raw_file = open(self.Raw_path,'w'); self.Sample_idx = 0
		for sampleObj in self.Sample_list:
			sampleObj.Index_int = self.Sample_idx; self.Sample_idx += 1
			self.SampleProcess_method(sampleObj)
			##End for

		#[3]
		self.Raw_file.close()
		self.Left_list = filter(lambda x:len(x.Var_sites)>=2, self.Left_dic.values())
		self.Right_list = filter(lambda x:len(x.Var_sites)>=2, self.Right_dic.values())
		self.Unit_list = self.Left_list + self.Right_list

		#[4]
		self.Unit_dic = {}
		for unitObj in self.Unit_list:
			self.Unit_dic[unitObj.Fixed_site] = unitObj
			##End for
		##End UnitExtract_method
	
	def SampleProcess_method(self, sampleObj):
		for bed_line in open(sampleObj.Junction_path):
			#[1]
			chr_name, bed_st, bed_ed, x, junc_cnt, junc_str = bed_line.split()
			junc_st = int(bed_st)+1; junc_ed = int(bed_ed); junc_cnt = float(junc_cnt)
			if not self.sLibrary_bool:
				junc_str = "_"
				##End if
			if chr_name != self.Chr_name:
				continue
				##End if
			if junc_cnt == 0:
				continue
				##End if
			#[2]
			junc_st_key = junc_st, junc_str
			junc_ed_key = junc_ed, junc_str
			#[3]
			if junc_st_key not in self.Left_dic:
				unitObj = GenUnitObject()
				unitObj.Chr_name = chr_name
				unitObj.Fixed_site = junc_st_key
				unitObj.Var_sites = set([])
				unitObj.Strand_str = junc_str
				unitObj.Count_dic = {}
				self.Left_dic[junc_st_key] = unitObj
				##End if
			unitObj = self.Left_dic[junc_st_key]
			unitObj.Var_sites.add(junc_ed_key)
			#[4]
			if junc_ed_key not in self.Right_dic:
				unitObj = GenUnitObject()
				unitObj.Chr_name = chr_name
				unitObj.Fixed_site = junc_ed_key
				unitObj.Var_sites = set([])
				unitObj.Strand_str = junc_str
				unitObj.Count_dic = {}
				self.Right_dic[junc_ed_key] = unitObj
				##End if
			unitObj = self.Right_dic[junc_ed_key]
			unitObj.Var_sites.add(junc_ed_key)
			#[5]
			raw_line = makeLine([chr_name, junc_st, junc_ed, junc_str, sampleObj.Index_int, junc_cnt], tt)
			self.Raw_file.write(raw_line + nn)
			##End for
		##End SampleProcess_method
	
	def HashExtract_method(self):
		#[1]
		self.Sample_list = self.Case_list + self.Ref_list

		#[2]
		random.seed(0)
		self.Random_cases = random.sample(self.Case_list, min(len(self.Case_list), 3))
		self.Random_refs = random.sample(self.Ref_list, min(len(self.Ref_list), 3))
		self.Random_list = self.Random_cases + self.Random_refs

		#[2]
		self.Hash_dic = {}
		for sampleObj in self.Random_list:
			self.Get_SampleHash_method(sampleObj)	
			##End for

		#[3]
		for hashObj in self.Hash_dic.values():
			hashObj.Size_int = hashObj.Range_list[1]-hashObj.Range_list[0]
			##End for

		#[4]
		self.Hash_path = self.Unit_dir + ss + "HashList.txt"
		self.Hash_file = open(self.Hash_path, 'w')
		for hashObj in sorted(self.Hash_dic.values(), key=lambda x:x.Size_int):
			hash_line = makeLine([hashObj.Chr_name, hashObj.Size_int], tt)
			self.Hash_file.write(hash_line + nn)
			##End for
	
		#[5]
		self.Hash_file.close()
		##End HashExtract_method
	
	def Get_SampleHash_method(self, sampleObj):
		for bed_line in open(sampleObj.Junction_path):
			chr_name, bed_st, bed_ed, x, x, junc_str = bed_line.split()[:6]
			if chr_name not in self.Hash_dic:
				hashObj = GenHashObject(); hashObj.Chr_name = chr_name
				hashObj.Reference_bit = False
				hashObj.Range_list = [float("inf"), -float("inf")]
				self.Hash_dic[chr_name] = hashObj
				##End if
			hashObj = self.Hash_dic[chr_name]
			if sampleObj in self.Ref_list:
				hashObj.Reference_bit = True
				##End if
			hashObj.Range_list[0] = min(hashObj.Range_list[0], int(bed_st))
			hashObj.Range_list[1] = max(hashObj.Range_list[1], int(bed_ed))
			##End for
		##End Get_SampleHash_method
	##End GenUnitHandler



class GenSampleObject(object):
	def __init__(self):
		pass
		##End init
	##End GenSampleObject

class GenCommObject(object):
	def __init__(self):
		pass
		##End init
	##End GenSampleObject

class GenHashObject(object):
	def __init__(self):
		pass
		##End init
	##End GenHashObject

class GenUnitObject(object):
	def __init__(self):
		pass
		##End init
	##End GenUnitObject

def makePath(dirPath):
	return commands.getoutput("mkdir %s"%dirPath)
	##End makePath

def makeLine(tokenList, sepToken):
	return sepToken.join(map(str, tokenList))	
	##End makeLine

if __name__ ==  "__main__" :
	main()
	sys.exit()
	##End if
