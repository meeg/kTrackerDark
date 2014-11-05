#!/usr/bin/python

import os
import sys, optparse

# user options
usage = """usage: %prog [options]
Select only one option.  If no option is specified then do a complete make clean, then make
"""
parser = optparse.OptionParser(usage=usage)
parser.add_option( "-c", "--clean", dest="clean", default=False, action="store_true", help="Clean kTracker and TrackExtrapolator.")
parser.add_option( "-k", "--k-clean", dest="kclean", default=False, action="store_true", help="Clean kTracker.")

#get the args
(opts,args) = parser.parse_args()



#list of analysis tools that we usually want
analysisTools = ["sqlDataReader", "sqlResWriter"]

if len(sys.argv) == 1:
  os.chdir('TrackExtrapolator')
  os.system('make clean;make')
  os.chdir('..')
  os.system('make clean;make')
  for tool in analysisTools:
    if os.path.isfile(tool):
      os.unlink(tool)
    print "Compiling analysis tool:",tool
    os.system( "./compile analysis_tools/%s" % tool )
else:
  if opts.clean:
    os.chdir('TrackExtrapolator')
    os.system('make clean')
    os.chdir('..')
    os.system('make clean')
  elif opts.kclean:
    os.system('make clean')
