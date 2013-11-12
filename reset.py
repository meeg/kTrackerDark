#!/usr/bin/python

import os

os.chdir('TrackExtrapolator')
os.system('make clean;make')
os.chdir('..')
os.system('make clean;make')
