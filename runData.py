#!/usr/bin/python

import os
import sys
import time

## Run one job on a given schema
def runCmd(cmd):
    print cmd
    os.system(cmd)

## command line parser
exe = sys.argv[1]
suffix1 = sys.argv[2]
suffix2 = sys.argv[3]
runlist = 'runlist.txt'
if len(sys.argv) > 4:
    runlist = sys.argv[4]

## Decide how many jobs should be run at the same time
nJobsMax = 6
if len(sys.argv) > 5:
    nJobsMax = int(sys.argv[5])

## get the username
username = os.environ['USER']

## Read in run list
if os.path.isfile(runlist): 
    fin = open(runlist, 'r')
    schemas = [line.strip() for line in fin.readlines()]
    fin.close()
else:
    runlist = raw_input('Input the run list separated by space: ')
    schemas = [word.strip() for word in runlist.strip().split()]

## Check the active job list
nSubmitted = 0
nMinutes = 0.
while nSubmitted < len(schemas):
    # control the total number of running programs, no matter it's from this thread or nor
    nRunning = int(os.popen('pgrep -u %s %s | wc -l' % (username, exe)).read().strip())
    print(exe+': '+str(nMinutes)+' minutes passed, '+str(nSubmitted)+" submitted, "+str(nRunning)+' running ...' )
    for i in range(nRunning, nJobsMax):
    	## check if all the jobs are submitted
    	if nSubmitted >= len(schemas): break

        ## make the actual command
        inputFile = 'run_%s_%s.root' % (schemas[nSubmitted], suffix1)
        outputFile = 'run_%s_%s.root' % (schemas[nSubmitted], suffix2)
        logFile = 'log_%s_%s' % (schemas[nSubmitted], suffix2)
        cmd = './%s %s %s > %s &' % (exe, inputFile, outputFile, logFile)

        runCmd(cmd)
        nSubmitted = nSubmitted + 1
            
    time.sleep(60)
    nMinutes = nMinutes + 1.

## Only when all jobs are finished should the script quit, depanding on the jobs started by this script only
nRunning = int(os.popen('pgrep -u %s -g %d %s | wc -l' % (username, os.getpgrp(), exe)).read().strip())
while nRunning != 0:
    time.sleep(30)
    nMinutes = nMinutes + 0.5

    nRunning = int(os.popen('pgrep -u %s -g %d %s | wc -l' % (username, os.getpgrp(), exe)).read().strip())
    print(exe+': '+str(nMinutes)+' minutes passed, '+str(nSubmitted)+" submitted, "+str(nRunning)+' running ...' )
