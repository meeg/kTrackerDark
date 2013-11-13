#!/usr/bin/python

import os
import sys

os.chdir('TrackExtrapolator')
os.system('make clean;make')
os.chdir('..')
os.system('make clean;make')
