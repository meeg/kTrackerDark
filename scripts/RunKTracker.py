#!/usr/bin/env python
##################################################
# submit acceptance jobs given an input file of
#  points of x,y,z,p,theta
#------------
# Brian Tice - tice@anl.gov
# May 16, 2014
##################################################

import sys, os
from optparse import OptionParser, TitledHelpFormatter
from datetime import datetime

#make sure setup script has been sourced.  If it hasn't, this variable won't be known
ktracker_root = os.getenv("KTRACKER_ROOT")
if not ktracker_root:
  sys.exit( "ERROR: You must source the setup script for the kTracker package" )
seaquest_v = os.getenv("SEAQUEST_RELEASE")

timestamp = datetime.now().strftime( "%y%m%d_%H%M%S_%f" )

from CondorJobBuilder import CondorJobBuilder

# decode argument options
formatter = TitledHelpFormatter( max_help_position=22, width = 190 )
parser = OptionParser( formatter = formatter)

default_opts   = os.path.join( ktracker_root, "opts/default.opts" )
default_outdir = "/e906/data/users/%s/kTracker/%s" %( os.getenv("USER"), seaquest_v )

parser.add_option("--input", dest="input", help="Input file (required)", default=None )
parser.add_option("--opts", dest="template_opts", help="Options file to use as template (default=%default)", default=default_opts )
parser.add_option("--outdir", dest="outdir", help="Directory for output (default=%default)", default=default_outdir )
parser.add_option("--outtag", dest="outtag", help="Append this tag to output filenames", default="" )
parser.add_option("--n_events", dest="n_events", help="Number of events to process (default=all)", default=-1, type=int )
parser.add_option("--first_event", dest="first_event", help="First event to process (default=0)", default=0, type=int )
parser.add_option("--interactive", dest="interactive", help="Run interactively", action="store_true", default=False )

#opts for input file or runnumber
#opts for type of tracking

opts,args = parser.parse_args()

#make sure input file was given
if not opts.input:
  sys.exit( "ERROR: You must specify an input file with --input=<input>")

#use full path for input
opts.input = os.path.abspath(opts.input)
if not os.path.isfile(opts.input):
  sys.exit( "ERROR: Cannot find input file: %s", opts.input )

#separate grid/nogrid output areas for now
if opts.interactive:
  opts.outdir += "/nogrid"
else:
  opts.outdir += "/grid"

#outtag should start with _
if opts.outtag != "" and opts.outtag[0] != "_":
  opts.outtag = "_" + opts.outtag

#filenames
input_base = os.path.basename( opts.input )
output_base = input_base.replace( "run_","reco%s_" % opts.outtag )

#set output directories
reco_outdir = os.path.join( opts.outdir, "reco" )
if not os.path.exists( reco_outdir ):
  os.makedirs( reco_outdir )
  os.system( "chmod 01775 " + reco_outdir )

output = os.path.join( reco_outdir, output_base )

#setup logging
logdir = os.path.join( opts.outdir, "logfile")
logfile = os.path.join( logdir, "%s.log" % output_base )
if not os.path.exists( logdir ):
  os.makedirs(logdir)
  os.system( "chmod 01775 " + logdir )
if os.path.isfile( logfile ):
  os.unlink( logfile )



###
#create a wrapper script and opts file in /e906/app
wrapper_dir   = "/e906/app/users/%s/kTracker/gridwrapper/%s/" % ( os.getenv("USER"), seaquest_v )
wrapper_name  = os.path.join( wrapper_dir, "kFastTracking%s_%s.sh" % (opts.outtag, timestamp ) )
optsfile_name = os.path.join( wrapper_dir, "kFastTracking%s_%s.opts" % (opts.outtag, timestamp ) )

if not os.path.exists( wrapper_dir ):
  os.makedirs( os.path.dir(wrapper_dir) )

wrapper = open(wrapper_name, "w")
optsfile = open(optsfile_name, "w" )
os.system( "chmod 755 " + wrapper_name )
os.system( "chmod 755 " + optsfile_name )

#setup software
wrapper.write("source %s/setup_seaquest.sh \n" % os.getenv("SEAQUEST_SETUP_ROOT") )
wrapper.write("source %s/setup.sh \n" % os.getenv("KTRACKER_ROOT") )
if not opts.interactive:
  wrapper.write( "cd $_CONDOR_SCRATCH_DIR\n")

#form ktracker command
cmd = "%s/kFastTracking %s" % ( ktracker_root, optsfile_name )
wrapper.write(cmd + "\n")

wrapper.close()
os.system( "chmod 755 " + wrapper_name )

#write the opts file
for line in open( opts.template_opts ):
  line = line.strip()

  #don't touch spacing or comment lines
  if len(line) and line[0] != "#":

    #change template file as specified
    if "InputFile" in line:
      if opts.interactive:
        line = "InputFile %s" % opts.input
      else:
        line = "InputFile $CONDOR_DIR_INPUT/%s" % os.path.basename(opts.input)

    elif "OutputFile" in line:
      print "Changing output file"
      if opts.interactive:
        line = "OutputFile %s" % output
      else:
        line = "OutputFile $CONDOR_DIR_RECO/%s" % os.path.basename(output)

    elif opts.n_events>0 and "N_Events" in line:
      line = "N_Events %d" % opts.n_events

    elif "FirstEvent" in line:
      line = "FirstEvent %d" % opts.first_event

  optsfile.write( line + "\n" )
optsfile.close()

#########
#create the condor job
condor_job = CondorJobBuilder( executable = wrapper_name, logfile=logfile, interactive = opts.interactive, extra_args="-T" )
condor_job.addCondorDir("RECO", reco_outdir)
condor_job.addCondorFile( opts.input )

#want jobsub_prefix !
#condor_job.jobsub_prefix="kFastTracking"

#run/submit
condor_job.execute()
