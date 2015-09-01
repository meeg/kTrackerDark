#!/usr/bin/env python

import os
import sys
import time
import subprocess
from datetime import datetime
from optparse import OptionParser

import GridUtil

# parse all the commandline controls
parser = OptionParser('Usage: %prog [options]')
parser.add_option('-l', '--list', type = 'string', dest = 'list', help = 'List of run IDs', default = '')
parser.add_option('-c', '--config', type = 'string', dest = 'config', help = 'I/O configuration file for vertexing', default = '')
parser.add_option('-s', '--server', type = 'string', dest = 'server', help = 'MySQL server', default = 'seaquel.physics.illinois.edu')
parser.add_option('-p', '--port', type = 'int', dest = 'port', help = 'MySQL port', default = 3283)
parser.add_option('-o', '--output', type = 'string', dest = 'output', help = 'output schema pattern', default = '')
parser.add_option('-m', '--max', type = 'int', dest = 'nJobsMax', help = 'maximum number of uploading instance', default = 8)
parser.add_option('-e', '--errlog', type = 'string', dest = 'errlog', help = 'Failed command log', default = 'vertexingMonitor_err.log')
parser.add_option('-d', '--debug', action = 'store_true', dest = 'debug', help = 'Enable massive debugging output', default = False)
(options, args) = parser.parse_args()

if len(sys.argv) < 2:
    parser.parse_args(['--help'])

# initialize vertex configuration
conf = GridUtil.JobConfig(options.config)

# initialize the runlist
runIDs = [int(line.strip()) for line in open(options.list).readlines()]

# initiaize grid
GridUtil.gridInit()

# the uploader
uploader = os.path.join(os.getenv('KTRACKER_ROOT'), 'sqlResWriter')

# main part
finishedRuns = []
uploadedRuns = []
fout = open(GridUtil.workDir + '/' + options.errlog, 'w')
while len(uploadedRuns) < len(runIDs):
    # check the running status
    failedJobs = []
    for index, runID in enumerate(runIDs):
        # if it's already checked
        if runID in finishedRuns:
            continue

        targetFile = os.path.join(conf.outdir, 'vertex', GridUtil.version, GridUtil.getSubDir(runID), 'vertex_%06d_%s.root' % (runID, GridUtil.version))
        nTotalJobs, nFinishedJobs, failedOpts = GridUtil.getJobStatus(conf.outdir, 'vertex', runID)
        for opt in failedOpts:
            failedJobs.append(GridUtil.makeCommandFromOpts('vertex', opt, conf))
        if (not os.path.exists(targetFile)) or len(failedOpts) != 0:
            fout.write('%s: %06d %02d %02d %02d %s\n' % (datetime.now(), runID, nTotalJobs, nFinishedJobs, len(failedOpts), 'certain jobs failed'))
            continue
        if nTotalJobs == 0 or nTotalJobs != nFinishedJobs:
            continue

        finishedRuns.append(runID)
        uploadedRuns.append(runID)

    # resubmit failed jobs if needed
    GridUtil.submitAllJobs(failedJobs, 'vertexMonitor_err.log_'+GridUtil.getTimeStamp())

    # submit the uploader to background
    #nRunning = int(os.popen('pgrep %s | wc -l' % uploader.split('/')[-1]).read().strip())
    #toBeUploadedRuns = [runID for runID in finishedRuns if runID not in uploadedRuns]
    #nJobs = options.nJobsMax - nRunning
    #if nJobs > len(toBeUploadedRuns):
    #    nJobs = len(toBeUploadedRuns)

    #print '%s: %s uploader running, will submit %d more.' % (datetime.now(), nRunning, nJobs)
    #time.sleep(60)    # this is to ensure the copy from node to server is completed
    #for index in range(nJobs):
    #    runID = toBeUploadedRuns[index]
    #    sourceFile = os.path.join(conf.outdir, 'vertex', GridUtil.version, GridUtil.getSubDir(runID), 'vertex_%06d_%s.root' % (runID, GridUtil.version))
    #    targetSchema = options.output % runID
    #    uploadLog = os.path.join(GridUtil.workDir, 'log_upload_%06d' % runID)
    #    uploadErr = os.path.join(GridUtil.workDir, 'err_upload_%06d' % runID)

    #    cmd = '%s %s %s %s %d 1> %s 2> %s &' % (uploader, sourceFile, targetSchema, options.server, options.port, uploadLog, uploadErr)
    #    print cmd
    #    os.system(cmd)

    #    uploadedRuns.append(runID)

    # sleep for 1 minutes
    fout.flush()
    print '%s: %d/%d finished, %d uploaded. ' % (datetime.now(), len(finishedRuns), len(runIDs), len(uploadedRuns))
    time.sleep(600)

fout.close()
stopGridGuard()
