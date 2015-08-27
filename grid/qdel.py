#!/usr/bin/env python

import os
import sys
import re
import subprocess
from optparse import OptionParser

# parse all the commandline controls
parser = OptionParser('Usage: %prog [options]')
parser.add_option('-r', '--range', type = 'string', dest = 'range', help = 'range of job ID, e.g. 101-202', default = '')
parser.add_option('-p', '--pattern', type = 'string', dest = 'pattern', help = 'pattern of job name', default = '')
parser.add_option('-s', '--status', type = 'string', dest = 'status', help = 'running status of jobs', default = '')
(options, args) = parser.parse_args()

if len(sys.argv) < 2:
    parser.parse_args(['--help'])

jobIDmin = -1
jobIDmax = 999999
if options.range != '':
    jobIDmin, jobIDmax = options.range.split('-')

pattern = '[a-z]'
if options.pattern != '':
    pattern = options.pattern

status = options.status

output, err = subprocess.Popen('jobsub_q | grep liuk', stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True).communicate()
jobDetails = []
for line in output.split('\n'):
    vals = line.strip().split()
    jobID = int(re.findall(r'^(\d{7})', vals[0])[0])
    jobDetails.append((jobID, vals[0], vals[8], vals[5]))
jobDetails.sort(key = lambda x : x[0])

toBeKilled = []
for job in jobDetails:
    if job[0] < jobIDmin or job[0] > jobIDmax:
        continue
    if len(re.findall(pattern, job[2])) == 0:
        continue
    if status not in job[3]:
        continue

    toBeKilled.append(job[1])

for index, url in enumerate(toBeKilled):
    print '%d/%d' % (index+1, len(toBeKilled)),
    os.system('jobsub_rm --jobid=%s')
