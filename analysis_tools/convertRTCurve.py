#!/usr/bin/python

import os
import sys
import re

sourceDir = sys.argv[1]
rt_inputFiles = ['rt_curve_D1.txt', 'rt_curve_D2.txt', 'rt_curve_D3p.txt', 'rt_curve_D3m.txt']
t0_inputFiles = ['time_window_D1.txt', 'time_window_D2.txt', 'time_window_D3p.txt', 'time_window_D3m.txt']

map_detectorID = {'D1U': 1, 'D1Up': 2, 'D1X': 3, 'D1Xp': 4, 'D1V': 5, 'D1Vp': 6, 
                  'D2V': 7, 'D2Vp': 8, 'D2Xp': 9, 'D2X': 10, 'D2U': 11, 'D2Up': 12, 
                  'D3pVp': 13, 'D3pV': 14, 'D3pXp': 15, 'D3pX': 16, 'D3pUp': 17, 'D3pU': 18, 
                  'D3mVp': 19, 'D3mV': 20, 'D3mXp': 21, 'D3mX': 22, 'D3mUp': 23, 'D3mU': 24}
map_detectorName = {}
for (key, val) in map_detectorID.iteritems():
	map_detectorName[val] = key

# read t0
t0 = {}
for t0_input in t0_inputFiles:	
	fin = open(sourceDir+'/'+t0_input, 'r')
	for line in fin.readlines():
		data = line.strip().split()
		t0[data[0]] = float(data[1])
	fin.close()

# read rt curves
rt_curve = {} # key = detectorID, val = dictionary from tdcTime to 
for rt_input in rt_inputFiles:
	fin = open(sourceDir+'/'+rt_input, 'r')
	for line in fin.readlines():
		data = line.strip().split()

		detectorID = map_detectorID[data[0]]
		if not detectorID in rt_curve:
			rt_curve[detectorID] = {}

		rt_curve[detectorID][t0[data[0]] - float(data[4])] = float(data[5])
	fin.close()

# output to kTracker recognizable format
fout = open(sys.argv[2], 'w')
for detectorID in rt_curve.keys():
	rt_sorted = sorted(rt_curve[detectorID].iteritems(), key = lambda d: d[0])
	fout.write(str(detectorID)+'  '+str(len(rt_sorted))+'  '+str(rt_sorted[0][0])+' '+str(rt_sorted[-1][0])+' '+map_detectorName[detectorID]+'\n')
	for (tdcTime, driftDistance) in rt_sorted:
		fout.write('0   '+str(tdcTime)+'  '+str(driftDistance)+'\n')
fout.close()


