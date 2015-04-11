#!/usr/bin/python2.7 -tt

import fileinput
import sys
import subprocess
import os
from numpy import *

statfile = open("statfile", "w")
for i in range(5, 20):
	o = open("output_verynew",'w')
	for line in open("B_bcc_001.omd"):
		line = line.replace("cutoffRadius = "+"%.1f"%i,"cutoffRadius = "+"%.1f"%(i+1))
		o.write(line)
	o.close()
	subprocess.call(['mv','output_verynew', 'B_bcc_001.omd'])
	subprocess.call(['../../../build/bin/openmd', 'B_bcc_001.omd'])
	f = open("B_bcc_001.stat")
	f.next()
	for line in f:
		statfile.write(line)
	f.close()
statfile.close()

data = loadtxt("statfile", double)
electrostaticPotential = data[:,8]
reciporcalSpacePotential = data[:,9]
total = electrostaticPotential + reciporcalSpacePotential

f= open("cutoffVsEnergy", 'w')
for i in range(6, 21):
	print >> f, i, total[i]
f.close()  

