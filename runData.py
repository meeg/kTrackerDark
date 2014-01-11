#!/usr/bin/python

import os
import sys
import time

## Run one job on a given schema
def runCmd(cmd):
    print cmd
    #os.system(cmd)

## command line parser
exe = sys.argv[1]
runlist = sys.argv[2]
suffix1 = sys.argv[3]
suffix2 = sys.argv[4]

## Read in run list
fin = open(sys.argv[2], 'r')
schemas = [line.strip() for line in fin.readlines()]

## Decide how many jobs should be run at the same time
nJobsMax = 6
if len(sys.argv) > 5:
    nJobsMax = int(sys.argv[5])

## Check the active job list
nSubmitted = 0
nMinutes = 0.
while nSubmitted < len(schemas):
    nRunning = int(os.popen('pgrep '+exe+' | wc -l').read().strip())
    print(sys.argv[1]+': '+str(nMinutes)+' minutes passed, '+str(nSubmitted)+" submitted, "+str(nRunning)+' running ...' )
    for i in range(nRunning, nJobsMax):
        ## make the actual command
        cmd = './'+exe+' run_'+schemas[nSubmitted]+'_'+suffix1+'.root run_'+schemas[nSubmitted]+'_'+suffix2+'.root > log_'+exe+'_'+schemas[nSubmitted]+' &'
        runCmd(cmd)
	nSubmitted = nSubmitted + 1
        
	if nSubmitted >= len(schemas): break
    
    time.sleep(60)
    nMinutes = nMinutes + 1.

## Only when all jobs are finished should the script quit
while nRunning != 0:
    time.sleep(30)
    nMinutes = nMinutes + 0.5

    nRunning = int(os.popen('pgrep '+exe+' | wc -l').read().strip())
    print(exe+': '+str(nMinutes)+' minutes passed, '+str(nSubmitted)+" submitted, "+str(nRunning)+' running ...' )
