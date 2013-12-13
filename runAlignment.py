#!/usr/bin/python

import os
import sys
import random
import time

def runCmd(cmd):
    print cmd
    #os.system(cmd)

runID = sys.argv[1]
nCycle = int(sys.argv[2])
if len(sys.argv) == 4:
    offset = int(sys.argv[3])
else:
    offset = 1

nEvtMax = 60000
nJobs = 5

for i in range(offset, nCycle+1):
    print 'Working on the '+str(i)+'th optimization cycle ... '

    rawFile = 'run_'+runID+'_raw.root'
    alignFile = 'run'+runID+'_align_'+str(i)+'.root'
    recFile_initial = 'rec_'+runID+'_align_'+str(i)

    runCmd('./update '+rawFile+' '+alignFile)

    nEvents_single = nEvtMax/nJobs
    for j in range(nJobs):
        runCmd('./kFastTracking '+alignFile+' '+recFile_initial+'_'+str(j+1)+'.root '+str(j*nEvents_single)+' '+str(nEvents_single)+' > log_'+runID+'_'+str(j+1)+' &')

    time.sleep(300)
    
    nMinutes = 4
    while int(os.popen('ps | grep kFastTracking | wc -l').read().strip()) != 0:
        nMinutes = nMinutes+1
        print str(nMinutes)+' minutes passed and tracking is not finished, wait for another 1 minute ...'
	time.sleep(60)

    runCmd('hadd '+recFile_initial+'.root '+recFile_initial+'_[1-'+str(nJobs)+'].root')
    runCmd('./milleAlign '+recFile_initial+'.root align_mille_'+str(i)+'.txt increament.log_'+str(i)+' > log_mille_'+str(i))
    runCmd('mv align_eval.root align_eval_'+str(i)+'.root')
    runCmd('cp align_mille_'+str(i)+'.txt align_mille.txt')
