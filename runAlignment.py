#!/usr/bin/python

import os
import sys
import random
import time

def runCmd(cmd):
    print cmd
    os.system(cmd)

def nEventAll():
    n = 0
    for i in range(2167, 2169):
        output = os.popen('grep Process log_'+str(i)+' | wc -l')
	text = output.read().strip()
	n = n + int(text)

    print str(n)+'/157764 ... '
    return n
    

nCycle = int(sys.argv[1])
if len(sys.argv) == 3:
    offset = int(sys.argv[2])
else:
    offset = 1

for i in range(offset, nCycle+1):
    print 'Working on the '+str(i)+'th optimization cycle ... '

    runCmd('./update run_2167_raw.root run_2167_align_'+str(i)+'.root;./update run_2168_raw.root run_2168_align_'+str(i)+'.root');
    runCmd('./kFastTracking run_2167_align_'+str(i)+'.root'+' rec_2167_align_'+str(i)+'.root > log_2167 &')
    runCmd('./kFastTracking run_2168_align_'+str(i)+'.root'+' rec_2168_align_'+str(i)+'.root > log_2168 &')

    nMinutes = 0
    while nEventAll() != 157764:
        nMinutes = nMinutes+1
        print str(nMinutes)+' minutes passed and tracking is not finished, wait for another 1 minute ...'
	time.sleep(60)

    runCmd('hadd rec_align_'+str(i)+'.root rec_216[7,8]_align_'+str(i)+'.root')
    runCmd('./milleAlign rec_align_'+str(i)+'.root align_mille_'+str(i)+'.txt increament.log_'+str(i)+' > log_mille_'+str(i))
    runCmd('./makeRTProfile rec_align_'+str(i)+'.root calibration_'+str(i)+'.txt')
    runCmd('./hodoAlign rec_align_'+str(i)+'.root hodoAlign_temp.root alignment_hodo_'+str(i)+'.txt > log_hodo_'+str(i))
    runCmd('./propAlign rec_align_'+str(i)+'.root propAlign_temp.root alignment_prop_'+str(i)+'.txt > log_prop_'+str(i))
    runCmd('mv align_eval.root align_eval_'+str(i)+'.root')
    runCmd('cp align_mille_'+str(i)+'.txt align_mille.txt')
    runCmd('cp alignment_hodo_'+str(i)+'.txt alignment_hodo.txt')
    runCmd('cp alignment_prop_'+str(i)+'.txt alignment_prop.txt')
    runCmd('cp calibration_'+str(i)+'.txt calibration.txt')
