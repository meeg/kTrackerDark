#!/usr/bin/python

import os
import sys
import random
import time

def runCmd(cmd):
    print cmd
    os.system(cmd)

def prepareConf(log_prev, conf):
    # read previous log and decide whether to turn on/off a detector alignment
    controls = []
    if os.path.isfile(log_prev):
        fin = open(log_prev, 'r')
        for line in fin.readlines():
            delta = float(line.strip().split()[3])

            if abs(delta) < 0.002:
                controls.append(0)
            else:
                controls.append(1)
    else:
        for detectorID in xrange(1, 25):
            controls.append(0)

    # adjust the controls for D1 and D2
    for index in xrange(0, 12, 2):
        if controls[index] == 1 or controls[index+1] == 1:
            controls[index] = 1
            controls[index+1] = 1

    # adjust the controls for D3p and D3m
    for index in xrange(12, 18):
        if controls[index] == 1:
            for detectorID in xrange(12, 18):
                controls[detectorID] = 1
            break

    for index in xrange(18, 24):
        if controls[index] == 1:
            for detectorID in xrange(18, 24):
                controls[detectorID] = 1
            break   
    
    # save the results
    fout = open(conf, 'w')
    for index, onoff in enumerate(controls):
        fout.write(str(index+1)+'  '+str(onoff)+'\n')
    fout.close()

runID = sys.argv[1]
nCycle = int(sys.argv[2])
if len(sys.argv) == 4:
    offset = int(sys.argv[3])
else:
    offset = 1

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
        runCmd('./kFastTracking '+alignFile+' '+recFile_initial+'_'+str(j+1)+'.root '+str(j*nEvents_single)+' '+str(nEvents_single)+' > log_'+runID+'_'+str(j+1)+' &')
    time.sleep(300)
    
    # check if all jobs are done running
    nMinutes = 4
    while int(os.popen('ps | grep kFastTracking | wc -l').read().strip()) != 0:
        nMinutes = nMinutes+1
        print str(nMinutes)+' minutes passed and tracking is not finished, wait for another 1 minute ...'
        time.sleep(60)
    
    # combine the outputs
    runCmd('hadd '+recFile_initial+'.root '+recFile_initial+'_[1-'+str(nJobs)+'].root')

    # clean up space
    runCmd('rm '+recFile_initial+'_[1-'+str(nJobs)+'].root')
    runCmd('rm '+alignFile)

    # chamber alignment based on millepede
    prepareConf('increament.log_'+str(i), 'mille.conf')
    runCmd('./milleAlign '+recFile_initial+'.root align_mille_'+str(i)+'.txt increament.log_'+str(i)+' > log_mille_'+str(i))
    runCmd('mv align_eval.root align_eval_'+str(i)+'.root')
    runCmd('cp align_mille_'+str(i)+'.txt align_mille.txt')
    runCmd('cp mille.conf mille.conf_'+str(i))
    
    # hodoscope alignment
    runCmd('./hodoAlign '+recFile_initial+'.root alignment_hodo_'+str(i)+'.txt')
    runCmd('mv hodo_eval.root hodo_eval_'+str(i)+'.root')
    runCmd('cp alignment_hodo_'+str(i)+'.txt alignment_hodo.txt')

    # hodoscope alignment
    runCmd('./propAlign '+recFile_initial+'.root alignment_prop_'+str(i)+'.txt')
    runCmd('mv prop_eval.root prop_eval_'+str(i)+'.root')
    runCmd('cp alignment_prop_'+str(i)+'.txt alignment_prop.txt')
    
    # chamber calibration
    runCmd('./makeRTProfile '+recFile_initial+'.root calibration_'+str(i)+'.txt')
    runCmd('mv cali_eval.root cali_eval_'+str(i)+'.root')
    runCmd('cp calibration_'+str(i)+'.txt calibration.txt')
