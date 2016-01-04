#!/usr/bin/python

import os
import sys
import random
import time

def runCmd(cmd):
    print cmd
    os.system(cmd)

def checkOverflow(output):
    if not os.path.exists(output):
        return True

    for line in open(output).readlines():
        vals = [float(val) for val in line.strip().split()]
        for val in vals:
            if abs(val) > 1.:
                return False
    return True

def prepareConf(log_prev, conf, fudge = 1.):
    print 'Generating configuration file for millepede according to '+log_prev, 'with fudge factor of '+str(fudge)

    if os.path.isfile(conf):
    	return

    # read previous log and decide whether to turn on/off a detector alignment
    sigma = [[0.1, 0.005, 0.05] for i in range(24)]
    if os.path.isfile(log_prev):
        fin = open(log_prev, 'r')
        lines = fin.readlines()
        for index, line in enumerate(lines):
            delta = float(line.strip().split()[3])

            if abs(delta) < sigma[index][2]:
            	factor = fudge*abs(delta)/sigma[index][2]
            	for i in range(3):
            		sigma[index][i] = sigma[index][i]*factor

    # save the results
    fout = open(conf, 'w')
    for index, oneline in enumerate(sigma):
        fout.write('%d      %f       %f         %f\n' % (index+1, oneline[0], oneline[1], oneline[2]))
    fout.close()

## command line control
runID = sys.argv[1]
nCycle = int(sys.argv[2])
if len(sys.argv) > 3:
    offset = int(sys.argv[3])
else:
    offset = 1
skipTracking = (len(sys.argv) > 4)

## Key performance knob
nEvtMax = 600000
nJobs = 8

for i in range(offset, nCycle+1):
    print 'Working on the '+str(i)+'th optimization cycle ... '

    # define the file names and apply the current alignment parameters
    rawFile = 'run_'+runID+'_raw.root'
    alignFile = 'run_'+runID+'_align_'+str(i)+'.root'
    recFile_initial = 'rec_'+runID+'_align_'+str(i)

    if not skipTracking:
        runCmd('./update ac '+rawFile+' '+alignFile)

    # divide the task to nJobs jobs and submit them all to background
    nEvents_single = nEvtMax/nJobs
    for j in range(nJobs):
    	if not skipTracking:
            runCmd('./kFastTracking %s %s_%02d.root %d %d > log_%s_%d &' % (alignFile, recFile_initial, j+1, nEvents_single, j*nEvents_single, runID, j+1))

    # check if all jobs are done running
    nMinutes = 0
    while int(os.popen('pgrep -u %s -g %d kFastTracking | wc -l' % (os.environ['USER'], os.getpgrp())).read().strip()) != 0:
        nMinutes = nMinutes+1
        sys.stdout.write('\r'+str(nMinutes)+' minutes passed and tracking is not finished, wait for another 1 minute ...')
        sys.stdout.flush()
        time.sleep(60)
    sys.stdout.write('\n')

    # combine the outputs
    if not skipTracking:
        runCmd('hadd '+recFile_initial+'.root '+recFile_initial+'_*.root')

    # clean up space
    runCmd('rm log_'+str(runID)+'_*')
    runCmd('rm '+recFile_initial+'_*.root')
    runCmd('rm '+alignFile)

    # reset the skip tracking flag so it's only effective in the first cycle
    skipTracking = False

    # chamber alignment based on millepede
    nTry = 1
    while nTry < 20:
        prepareConf('increament.log_'+str(i-1), 'mille.conf', nTry)
        runCmd('./milleAlign '+recFile_initial+'.root align_mille_'+str(i)+'.txt increament.log_'+str(i)+' > log_mille_'+str(i))
        if checkOverflow('increament.log_'+str(i)):
            break
        runCmd('rm mille.conf')
        nTry = nTry + 1

    if not checkOverflow('increament.log_'+str(i)):
        sys.exit()
    runCmd('mv align_eval.root align_eval_'+str(i)+'.root')
    runCmd('cp align_mille_'+str(i)+'.txt alignment/run4/align_mille.txt')
    runCmd('mv mille.conf mille.conf_'+str(i)+'_'+str(nTry))

    # hodoscope alignment
    #runCmd('./hodoAlign '+recFile_initial+'.root alignment_hodo_'+str(i)+'.txt')
    #runCmd('mv hodo_eval.root hodo_eval_'+str(i)+'.root')
    #runCmd('cp alignment_hodo_'+str(i)+'.txt alignment_hodo.txt')

    # prop. tube alignment
    #for j in range(10):
    #	runCmd('./propAlign '+recFile_initial+'.root alignment_prop_temp.txt')
    #	runCmd('mv alignment_prop_temp.txt alignment_prop.txt')
    #runCmd('mv prop_eval.root prop_eval_'+str(i)+'.root')
    #runCmd('cp alignment_prop.txt alignment_prop_'+str(i)+'.txt')

    # chamber calibration
    #runCmd('./makeRTProfile '+recFile_initial+'.root calibration_'+str(i)+'.txt')
    #runCmd('mv cali_eval.root cali_eval_'+str(i)+'.root')
    #runCmd('cp calibration_'+str(i)+'.txt calibration.txt')
