import sys
import math
import logging
import os
import subprocess
import logging

fraw_list = []#List of all .md files found (even the includes).
fmd_list = []#List of all config .md files that can be run (not the includes).
fmd_base_list = []#List of the names of the .md files that can be run.

dir_cwd = ""#Current working directory.
dir_openmd = ""#Absolute path for openmd
dir_base = ""#Directory where the script is run from.

FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT)

"""
Function sets up the dir_base and dir_openmd. If an openmd executable is not
found, script exits. Function looks for the openmd in the relative location
../build/bin but stores the openmd location as an absolute path.
@author Samuel Njoroge
"""
def setupDirectories():
	global dir_base, dir_openmd, dir_cwd
	logger = logging.getLogger("tcpserver")
	dir_base = os.getcwd()
	if(os.path.isfile("../build/bin/openmd")):
		os.chdir("../build/bin/")
		dir_openmd = os.getcwd()
		os.chdir(dir_base)
	else:
		logger.error("OpenMD : %s", "openmd executable not found at the expected location. Script Will Quit...")
		sys.exit()
	
	
"""
Function checks if the sample_file and validate_file (.md files) have the same
statusTime = interval time for the stats file.
@author Samuel Njoroge
@param string sample_file - .md file that is being run.
@param string validate_file - .md file the result is being compared to.
@return boolean
"""
def validate_md_time(sample_file, validate_file):
  sample_status_time = 0
  sample_sample_time = 0
  sample_run_time = 0
  validate_status_time = 0
  validate_sample_time = 0
  validate_run_time = 0
  logger = logging.getLogger("tcpserver")

  samplefh = open(sample_file, "r")
  validatefh = open(validate_file, "r")
  
  line = samplefh.readline()
  while line:
    if "statusTime" in line:
      arr = line.split('=')
      temp = arr[1]
      sample_status_time = float(temp.replace(';', ''))
    elif "sampleTime" in line:
      arr = line.split('=')
      temp = arr[1]
      sample_sample_time = float(temp.replace(';', ''))
    elif "runTime" in line:
      arr = line.split('=')
      temp = arr[1]
      sample_run_time = float(temp.replace(';', ''))

    line = samplefh.readline()

  line = validatefh.readline()
  while line:
    if "statusTime" in line:
      arr = line.split('=')
      temp = arr[1]
      validate_status_time = float(temp.replace(';', ''))
    elif "sampleTime" in line:
      arr = line.split('=')
      temp = arr[1]
      validate_sample_time = float(temp.replace(';', ''))
    elif "runTime" in line:
      arr = line.split('=')
      temp = arr[1]
      validate_run_time = float(temp.replace(';', ''))

    line = validatefh.readline()

  if (sample_status_time > 0) or (validate_status_time > 0):
    if sample_status_time == validate_status_time:
      return True

  if (sample_sample_time > 0) or (validate_sample_time > 0):
    if sample_sample_time == validate_sample_time:
      return True

  if (sample_run_time > 0) or (validate_run_time > 0):
    if sample_run_time == validate_run_time:
      return True

  logger.warning("MD File: %s", "Sample/Validation times do not match.")
  return False
  
"""
Function checks if an .md config file and not an include file.
@author Samuel Njoroge
@param string file - .md file name
@return boolean
"""
def file_is_md(filename):
	file_handle = open(filename)
	line = file_handle.readline()
	
	while line:
		if "<OpenMD version=" in line:
			return True
		line = file_handle.readline()
			
	return False
	
"""
Function compares two numbers.
@author Samuel Njoroge and ().
@param float one - first number to compare.
@param float two - second number to compare.
@param boolean ignore_sign - if sign is a factor in the comparison.
@return float diff - difference between the two numbers.
"""
def absDiff(one, two, ignore_sign = False):
	if ignore_sign:
		one, two = math.fabs(one), math.fabs(two)
	return max(one, two) - min(one, two)
	
"""
Function compares two files.
@author Samuel Njoroge and ().
@param string fExpected - name of the expected file.
@param string fNew - name of the new test file.
@param float epsilon - Precision of the comparison of the files.
@param boolean ignore_sign - if sign will be a factor in comparing the digits.
@return boolean
"""
def compare(fExpected, fNew, epsilon = 0.00001, ignore_sign=False):
	logger = logging.getLogger("tcpserver")
	fone = open(fExpected, 'r')
	ftwo = open(fNew, 'r')

	diffs = 0
	
	i = 0
	for lineone in fone:
		linetwo = ftwo.readline()

		elementsone = lineone.split()
		elementstwo = linetwo.split()

		i = i + 1

		lenone = len(elementsone)
		lentwo = len(elementstwo)

		if lenone != lentwo:
			diffs = diffs + 1
			logger.warning("Line: %d - %s", i, "no mach")
			return True
		else:
			for j in range(lenone):
				# used to ignore XYZ meta data such as # of frames and # of atoms
				try:
					feone = int(elementsone[j])
					fetwo = int(elementstwo[j])
					# these are ints -- skip this pair
					continue
				except ValueError:
					pass
				try:
					feone = float(elementsone[j])
					fetwo = float(elementstwo[j])

					fediff = absDiff(feone, fetwo, ignore_sign)

					if fediff > epsilon:
						diffs = diffs + 1
						print "Line " + str(i) + " do not match in the files."
						return True
				except ValueError:
					pass
	return False#diffs == 0
	
"""
Function scans the directory for .md config files which will be run for
testing.
@author Samuel Njoroge
"""
def scanForMdFiles():
	paths = []
	variable = "sam is here"
	for p in os.listdir("../samples/"):
		paths.append(os.path.abspath("../samples/" + p))
	
	while len(paths):
		p = paths.pop()
		if os.path.isfile(p):
			fraw_list.append(p)
		else:
			for np in os.listdir(p):
				paths.append(os.path.abspath(p + "/" + np))
				
	for f in fraw_list:
		if (".md" in f) and (not ".svn" in f) and file_is_md(f):
			fmd_list.append(f)
			temp = os.path.basename(f)
			temp = temp.replace(".md",'')
			fmd_base_list.append(temp)
		
"""
Function runs OpenMD on the files the compares them with the already run files (*_v.stat)
@author Samuel Njoroge
"""
def runMdFiles():
	logger = logging.getLogger("tcpserver")
	global dir_base, dir_openmd, dir_cwd
	output = []
	for x in range(0, len(fmd_list)):
		#subprocess.call(["export FORCE_PARAM_PATH=/Users/snjoroge/Documents/openmd/development/forceFields"])
		if "argon" in fmd_list[x]:
			logger.debug("Switching to Directory: %s", os.path.dirname(fmd_list[x]))
			os.chdir(os.path.dirname(fmd_list[x]))
			logger.debug("Running: %s", fmd_list[x])
			output = subprocess.call([dir_openmd + "/openmd", fmd_list[x]])
			if(os.path.exists(os.path.dirname(fmd_list[x]) + "/" + fmd_base_list[x] + ".stat")):
				#print "Renaming File: " + fmd_base_list[x] + ".stat - " + fmd_base_list[x] + "_v.stat"
				#subprocess.call(["cp", os.path.dirname(fmd_list[x]) + "/" + fmd_base_list[x] + ".stat", os.path.dirname(fmd_list[x]) + "/" + fmd_base_list[x] + "_v.stat"])
				logger.debug("Comparing: %s", "Comparing: " + fmd_base_list[x] + ".stat <=> " + fmd_base_list[x] + "_v.stat")
				if(compare(os.path.dirname(fmd_list[x]) + "/" + fmd_base_list[x] + ".stat", os.path.dirname(fmd_list[x]) + "/" + fmd_base_list[x] + "_v.stat")):
					logger.warning("Files: %s", "Files do not match.")
				else:
					logger.debug("Files Match")
		os.chdir(dir_base)

def cleanUp():
	print "Delete all files generated."
	for x in range(0, len(fmd_list)):
		if(os.path.exists(os.path.dirname(fmd_list[x]) + "/" + fmd_base_list[x] + ".eor")):
			print "DELETE:" + os.path.dirname(fmd_list[x]) + "/" + fmd_base_list[x] + ".eor"
			os.remove(os.path.dirname(fmd_list[x]) + "/" + fmd_base_list[x] + ".eor")
		if(os.path.exists(os.path.dirname(fmd_list[x]) + "/" + fmd_base_list[x] + ".stat")):
			print "DELETE:" + os.path.dirname(fmd_list[x]) + "/" + fmd_base_list[x] + ".stat"
			os.remove(os.path.dirname(fmd_list[x]) + "/" + fmd_base_list[x] + ".stat")
		if(os.path.exists(os.path.dirname(fmd_list[x]) + "/" + fmd_base_list[x] + ".dump")):
			print "DELETE:" + os.path.dirname(fmd_list[x]) + "/" + fmd_base_list[x] + ".dump"
			os.remove(os.path.dirname(fmd_list[x]) + "/" + fmd_base_list[x] + ".dump")
	
