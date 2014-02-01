#!/usr/bin/python

import os
import sys
import random
import time

def runCmd(cmd):
    print cmd
    os.system(cmd)

def prepareConf(log_prev, conf):
    print 'Generating configuration file for millepede according to '+log_prev

    # read previous log and decide whether to turn on/off a detector alignment
    sigma = [[0.1, 0.005, 0.05] for i in range(24)]
    if os.path.isfile(log_prev):
        fin = open(log_prev, 'r')
        lines = fin.readlines()
        for index, line in enumerate(lines):
            delta = float(line.strip().split()[3])

            if abs(delta) < sigma[index][2]:
            	factor = abs(delta)/sigma[index][2]
            	for i in range(3):
            		sigma[index][i] = sigma[index][i]*factor
    
    # save the results
    fout = open(conf, 'w')
    for index, oneline in enumerate(sigma):
        fout.write('%d      %f       %f         %f\n' % (index+1, oneline[0], oneline[1], oneline[2]))
    fout.close()

## command line control
if( len(sys.argv) < 2 ):
  sys.exit("You need more args.  Script needs a help menu to say what they are.")
runID = sys.argv[1]
nCycle = int(sys.argv[2])
if len(sys.argv) == 4:
    offset = int(sys.argv[3])
else:
    offset = 1

## Key performance knob
nEvtMax = 100000
nJobs = 5

for i in range(offset, nCycle+1):
    print 'Working on the '+str(i)+'th optimization cycle ... '

    # define the file names and apply the current alignment parameters
    rawFile = 'run_'+runID+'_raw.root'
    alignFile = 'run_'+runID+'_align_'+str(i)+'.root'
    recFile_initial = 'rec_'+runID+'_align_'+str(i)

    runCmd('./update '+rawFile+' '+alignFile)

    # divide the task to nJobs jobs and submit them all to background
    nEvents_single = nEvtMax/nJobs
    for j in range(nJobs):
    	runCmd('./kFastTracking %s %s_%d.root %d %d > log_%s_%d &' % (alignFile, recFile_initial, j+1, j*nEvents_single, nEvents_single, runID, j+1))
    time.sleep(300)
    
    # check if all jobs are done running
    nMinutes = 4
    while int(os.popen('pgrep -u %s -g %d kFastTracking | wc -l' % (os.environ['USER'], os.getpgrp())).read().strip()) != 0:
        nMinutes = nMinutes+1
        sys.stdout.write('\r'+str(nMinutes)+' minutes passed and tracking is not finished, wait for another 1 minute ...')
        sys.stdout.flush()
        time.sleep(60)
    sys.stdout.write('\n')
    
    # combine the outputsxw
    runCmd('hadd '+recFile_initial+'.root '+recFile_initial+'_[1-'+str(nJobs)+'].root')

    # clean up space
    runCmd('rm log_'+str(runID)+'_[1-'+str(nJobs)+']')
    runCmd('rm '+recFile_initial+'_[1-'+str(nJobs)+'].root')
    runCmd('rm '+alignFile)

    # chamber alignment based on millepede
    prepareConf('increament.log_'+str(i-1), 'mille.conf')
    runCmd('./milleAlign '+recFile_initial+'.root align_mille_'+str(i)+'.txt increament.log_'+str(i)+' > log_mille_'+str(i))
    runCmd('mv align_eval.root align_eval_'+str(i)+'.root')
    runCmd('cp align_mille_'+str(i)+'.txt align_mille.txt')
    runCmd('mv mille.conf mille.conf_'+str(i))
    
    # hodoscope alignment
    runCmd('./hodoAlign '+recFile_initial+'.root alignment_hodo_'+str(i)+'.txt')
    runCmd('mv hodo_eval.root hodo_eval_'+str(i)+'.root')
    runCmd('cp alignment_hodo_'+str(i)+'.txt alignment_hodo.txt')

    # hodoscope alignment
    #runCmd('./propAlign '+recFile_initial+'.root alignment_prop_'+str(i)+'.txt')
    #runCmd('mv prop_eval.root prop_eval_'+str(i)+'.root')
    #runCmd('cp alignment_prop_'+str(i)+'.txt alignment_prop.txt')
    
    # chamber calibration
    #runCmd('./makeRTProfile '+recFile_initial+'.root calibration_'+str(i)+'.txt')
    #runCmd('mv cali_eval.root cali_eval_'+str(i)+'.root')
    #runCmd('cp calibration_'+str(i)+'.txt calibration.txt')
