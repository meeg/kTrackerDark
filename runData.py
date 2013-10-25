#!/usr/bin/python

import os
import sys
import time

## configure the job name, executable name and the final command submitted
exe_name = {}
command = {}

exe_name['update'] = 'update'
command['update']  = r"'./update data/'+schema+'.root data/'+schema+'_temp.root dsds &'"

exe_name['seed'] = 'kSeeder'
command['seed']  = r"'./kSeeder data/'+schema+'.root data/'+schema+'_seed.root > log_seed_'+schema+' &'"

exe_name['track'] = 'kTracker'
command['track']  = r"'./kTracker data/'+schema+'_seed.root data/'+schema+'_track.root > log_track_'+schema+' &'"

exe_name['vertex'] = 'kVertex'
command['vertex']  = r"'./kVertex data/'+schema+'_track.root data/'+schema+'_vertex.root > log_vertex_'+schema+' &'"

exe_name['vertexfast'] = 'kVertex_fast'
command['vertexfast']  = r"'./kVertex_fast data/'+schema+'_fast.root data/'+schema+'_vertexfast.root > log_vertexfast_'+schema+' &'"

exe_name['vertex_fast'] = 'kVertex'
command['vertex_fast']  = r"'./kVertex data/'+schema+'_fast.root data/'+schema+'_vertex_fast.root > log_vertex_fast_'+schema+' &'"

exe_name['fast'] = 'kFastTracking'
command['fast']  = r"'./kFastTracking data/'+schema+'.root data/'+schema+'_fast.root > log_fast_'+schema+' &'"

exe_name['read'] = 'sqlDataReader'
command['read']  = r"'./sqlDataReader '+schema+' data/'+schema+'.root &'" 

exe_name['test'] = 'sleep'
command['test']  = r"'sleep 61'"

## Run one job on a given schema
def runCmd(job, schema):
    schema = schema.rstrip()
    command_run = eval(command[job])
    print command
    os.system(command_run)

## Read in run list
fin = open(sys.argv[2], 'r')
schemas = fin.readlines()

## Decide how many jobs should be run at the same time
nJobsMax = 6
if len(sys.argv) > 3:
    nJobsMax = int(sys.argv[3])

## Check the active job list
nSubmitted = 0
while nSubmitted < len(schemas):
    nRunning = int(os.popen('ps | grep '+exe_name[sys.argv[1]]+' | wc -l').read().strip())
    print(sys.argv[1]+': '+str(nSubmitted)+" submitted, "+str(nRunning)+' running ...' )
    for i in range(nRunning, nJobsMax):
        runCmd(sys.argv[1], schemas[nSubmitted])
	nSubmitted = nSubmitted + 1
        
	if nSubmitted >= len(schemas): break

    time.sleep(60)

## Only when all jobs are finished should the script quit
while nRunning != 0:
    time.sleep(30)
    nRunning = int(os.popen('ps | grep '+exe_name[sys.argv[1]]+' | wc -l').read().strip())
    print(sys.argv[1]+': '+str(nSubmitted)+" submitted, "+str(nRunning)+' running ...' )
