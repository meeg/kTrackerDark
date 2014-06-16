#!/usr/bin/python

import os
import sys
import time
from optparse import OptionParser

## Run one job on a given schema
def runCmd(cmd):
    print cmd
    os.system('nice ' + cmd)

## command line parser
parser = OptionParser('Usage: %prog executable sources targets [options]')
parser.add_option('-l', '--list', type = 'string', dest = 'list', help = 'List of run IDs', default = '')
parser.add_option('-m', '--jobs', type = 'int', dest = 'nJobsMax', help = 'Maximum number of jobs running', default = 6)
parser.add_option('-n', '--notify', type = 'string', dest = 'notify', help = 'E-mail sent to notify the end of jobs', default = '')
parser.add_option('-o', '--output', type = 'string', dest = 'output', help = 'Output file name (i.e. call hadd at the end)', default = '')
parser.add_option('-s', '--suffix', type = 'string', dest = 'suffix', help = 'Additional arguments needed in commands', default = '')
(options, args) = parser.parse_args()

exe = args[0]
pattern1 = args[1]
pattern2 = args[2]
suffix = options.suffix
runlist = options.list
nJobsMax = options.nJobsMax

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
    print(exe+': '+str(nMinutes)+' minutes passed, '+str(nSubmitted)+'/'+str(len(schemas))+' submitted, '+str(nRunning)+' running ...' )
    for i in range(nRunning, nJobsMax):
        ## check if all jobs are submitted
    	if nSubmitted >= len(schemas): break

        ## make the actual command
        inputFile = pattern1.replace('?', schemas[nSubmitted])
        outputFile = pattern2.replace('?', schemas[nSubmitted])
        arguments = suffix.replace('?', schemas[nSubmitted])
        logFile = 'log_%s_%s' % (exe, schemas[nSubmitted])
        cmd = './%s %s %s %s > %s &' % (exe, inputFile, outputFile, arguments, logFile)

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
    print(exe+': '+str(nMinutes)+' minutes passed, '+str(nSubmitted)+'/'+str(len(schemas))+' submitted, '+str(nRunning)+' running ...' )

## run hadd to merge all final output to a single file
if '.root' in options.output:
	cmd = 'hadd -f ' + options.output
	for runID in schemas:
		cmd = cmd + ' ' + pattern2.replace('?', runID)
	runCmd(cmd)

## Send out notification if required
if '@' in options.notify:
	subject = '%s finished successfully on %d/%d jobs after %f minutes' % (exe, nSubmitted, len(schemas), nMinutes)
	content = str(sys.argv).strip('[]') + '\n' + str(schemas).strip('[]')
	os.system('echo "%s" | mail -s "%s" %s' % (content, subject, options.notify))