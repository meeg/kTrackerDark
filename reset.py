#!/usr/bin/python

import os
import sys

if len(sys.argv) == 1:
  os.chdir('TrackExtrapolator')
  os.system('make clean;make')
  os.chdir('..')
  os.system('make clean;make')
else:
  if sys.argv[1] == 'c':
    os.chdir('TrackExtrapolator')
    os.system('make clean')
    os.chdir('..')
    os.system('make clean')
  elif sys.argv[1] == 'k':
    os.system('make clean')
