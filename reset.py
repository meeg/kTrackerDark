#!/usr/bin/python

import os
import sys

print "What are the options?"

if len(sys.argv) == 1:
  os.chdir("kTrackerServices")
  os.system('make clean;make')
  os.chdir('../TrackExtrapolator')
  os.system('make clean;make')
  os.chdir('..')
  os.system('make clean;make')
else:
  if sys.argv[1] == 'c':
    os.chdir("kTrackerServices")
    os.system("make clean")
    os.chdir('../TrackExtrapolator')
    os.system('make clean')
    os.chdir('..')
    os.system('make clean')
  elif sys.argv[1] == 'k':
    os.system('make clean')
  elif sys.argv[1] == "m":   #make all
    os.chdir("kTrackerServices")
    os.system("make")
    os.chdir('../TrackExtrapolator')
    os.system('make')
    os.chdir('..')
    os.system('make')
