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
parser.add_option('-t', '--tconfig', type = 'string', dest = 'tconfig', help = 'I/O configuration file for tracking', default = '')
parser.add_option('-v', '--vconfig', type = 'string', dest = 'vconfig', help = 'I/O configuration file for vertexing', default = '')
parser.add_option('-e', '--errlog', type = 'string', dest = 'errlog', help = 'Failed command log', default = 'trackingMonitor_err.log')
parser.add_option('-d', '--debug', action = 'store_true', dest = 'debug', help = 'Enable massive debugging output', default = False)
(options, args) = parser.parse_args()

if len(sys.argv) < 2:
    parser.parse_args(['--help'])

# initialize general tracking configuration
tconf = GridUtil.JobConfig(options.tconfig)
vconf = GridUtil.JobConfig(options.vconfig)

# initialize the runlist
runIDs = [int(line.strip()) for line in open(options.list).readlines()]

# initiaize grid
GridUtil.gridInit()

# loop until all tracking jobs on list are done
nExist = 0
fout = open(GridUtil.workDir + '/' + options.errlog, 'w')
while nExist != len(runIDs):
    failedJobs = []
    vertexJobs = []
    nExist = 0
    for index, runID in enumerate(runIDs):
        targetFile = os.path.join(vconf.indir, 'track', GridUtil.version, GridUtil.getSubDir(runID), 'track_%06d_%s.root' % (runID, GridUtil.version))
        tempTargetFile = os.path.join(tconf.outdir, 'track', GridUtil.version, GridUtil.getSubDir(runID), 'track_%06d_%s.root' % (runID, GridUtil.version))
        vertexOpts = os.path.join(vconf.outdir, 'opts', GridUtil.version, GridUtil.getSubDir(runID), '%s_%06d_%s.opts' % (GridUtil.auxPrefix['vertex'], runID, GridUtil.version))

        #skip if this file has already been merged
        if os.path.exists(targetFile):
            nExist = nExist + 1
            if not os.path.exists(vertexOpts):
                vertexJobs.append(GridUtil.makeCommand('vertex', runID, vconf))
            continue

        #check the running status
        nTotalJobs, nFinishedJobs, failedOpts = GridUtil.getJobStatus(tconf.outdir, 'track', runID)
        if options.debug:
            print runID, nTotalJobs, nFinishedJobs, failedOpts
        for opt in failedOpts:
            failedJobs.append(GridUtil.makeCommandFromOpts('track', opt, tconf))
        if nTotalJobs != nFinishedJobs or len(failedOpts) != 0:
            if len(failedOpts) != 0:
                fout.write('%s: %06d %02d %02d %02d %s\n' % (datetime.now(), runID, nTotalJobs, nFinishedJobs, len(failedOpts), 'certain jobs failed'))
            continue

        # now proceed to merge
        mergeSuccessful = True
        sourceFile = os.path.join(tconf.outdir, 'track', GridUtil.version, GridUtil.getSubDir(runID), 'track_%06d_%s_*.root' % (runID, GridUtil.version))
        mergeOutput, mergeErr = subprocess.Popen('hadd %s %s' % (tempTargetFile, sourceFile), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True).communicate()

        for line in mergeErr.strip().split('\n'):
            if 'dictionary' not in line:
                mergeSuccessful = False

        if mergeSuccessful:
            print 'Run %06d finished and merged.' % runID
            if tempTargetFile == targetFile or GridUtil.runCommand('mv %s %s ' % (tempTargetFile, targetFile)):
                GridUtil.runCommand('rm ' + sourceFile)
            else:
                print 'Run %06d failed in moving to pnfs.' % runID
                fout.write('%s: %06d %02d %02d %02d %s\n' % (datetime.now(), runID, nTotalJobs, nFinishedJobs, len(failedOpts), 'moving to pnfs failed'))
            vertexJobs.append(GridUtil.makeCommand('vertex', runID, vconf))
        else:
            print 'Run %06d failed in merging.' % runID
            fout.write('%s: %06d %02d %02d %02d %s\n' % (datetime.now(), runID, nTotalJobs, nFinishedJobs, len(failedOpts), 'merging failed'))

    # submit the finished vertexing jobs
    GridUtil.submitAllJobs(vertexJobs, 'vertex_err.log_'+GridUtil.getTimeStamp())

    # submit all failed jobs
    GridUtil.submitAllJobs(failedJobs, 'track_err.log_'+GridUtil.getTimeStamp())

    # sleep for 1 minutes
    fout.flush()
    print '%s: %d/%d finished. ' % (datetime.now(), nExist, len(runIDs))
    time.sleep(600)

fout.close()
stopGridGuard()
