#!/usr/bin/python

import sys
import fnmatch
import os
import shutil
import argparse
from argparse import RawTextHelpFormatter
import time
import math

  
currentdir = os.getcwd()

parser = argparse.ArgumentParser(description='\
Takes a file name and checks if the submatrices constructed during\
the Faugere-Lachartre implementation exist. If so, it combines them\
to an ABCD matrix.',
formatter_class=RawTextHelpFormatter)
parser.add_argument('-n', '--name', required=True,
    help='name of the benchmark')

args = parser.parse_args()

#old method
if (os.path.isfile(args.name+'-oldMethod-sub-A.pbm') and
    os.path.isfile(args.name+'-oldMethod-sub-B.pbm') and 
    os.path.isfile(args.name+'-oldMethod-sub-C.pbm') and
    os.path.isfile(args.name+'-oldMethod-sub-D.pbm')):
  os.system('pnmcat -lr '+args.name+'-oldMethod-sub-A.pbm '+\
      args.name+'-oldMethod-sub-B.pbm > sub-AB.pbm')
  os.system('pnmcat -lr '+args.name+'-oldMethod-sub-C.pbm '+\
      args.name+'-oldMethod-sub-D.pbm > sub-CD.pbm')
  os.system('pnmcat -tb sub-AB.pbm sub-CD.pbm >\
    '+args.name+'-oldMethdod-sub-ABCD.pbm')
  os.system('rm -f sub-AB.pbm sub-CD.pbm')

#new method
if (os.path.isfile(args.name+'-newMethod-sub-A.pbm') and
    os.path.isfile(args.name+'-newMethod-sub-B.pbm') and 
    os.path.isfile(args.name+'-newMethod-sub-C.pbm') and
    os.path.isfile(args.name+'-newMethod-sub-D.pbm')):
  os.system('pnmcat -lr '+args.name+'-newMethod-sub-A.pbm '+\
      args.name+'-newMethod-sub-B.pbm > sub-AB.pbm')
  os.system('pnmcat -lr '+args.name+'-newMethod-sub-C.pbm '+\
      args.name+'-newMethod-sub-D.pbm > sub-CD.pbm')
  os.system('pnmcat -tb sub-AB.pbm sub-CD.pbm >\
    '+args.name+'-newMethdod-sub-ABCD.pbm')
  os.system('rm -f sub-AB.pbm sub-CD.pbm')

#multiline method
if (os.path.isfile(args.name+'-multilineMethod-sub-A.pbm') and
    os.path.isfile(args.name+'-multilineMethod-sub-B.pbm') and 
    os.path.isfile(args.name+'-multilineMethod-sub-C.pbm') and
    os.path.isfile(args.name+'-multilineMethod-sub-D.pbm')):
  os.system('pnmcat -lr '+args.name+'-multilineMethod-sub-A.pbm '+\
      args.name+'-multilineMethod-sub-B.pbm > sub-AB.pbm')
  os.system('pnmcat -lr '+args.name+'-multilineMethod-sub-C.pbm '+\
      args.name+'-multilineMethod-sub-D.pbm > sub-CD.pbm')
  os.system('pnmcat -tb sub-AB.pbm sub-CD.pbm >\
    '+args.name+'-multilineMethdod-sub-ABCD.pbm')
  os.system('rm -f sub-AB.pbm sub-CD.pbm')
