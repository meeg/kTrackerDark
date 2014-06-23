#!/usr/bin/python

import os
import sys

#list of analysis tools that we usually want
analysisTools = ["sqlDataReader"]

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
  if sys.argv[1] == 'c':
    os.chdir('TrackExtrapolator')
    os.system('make clean')
    os.chdir('..')
    os.system('make clean')
  elif sys.argv[1] == 'k':
    os.system('make clean')
