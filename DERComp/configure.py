
#!/usr/bin/env python

import os
import sys
import shutil
import ROOT
import argparse

parser = argparse.ArgumentParser(description='Configure rqlib')
parser.add_argument('file', type=str, help='Sample LZap file to generate rqlib')
parser.add_argument('--force', action='store_true', help='Force rqlib generation')
args = parser.parse_args()

filename = args.file

if not os.path.isdir('lib'):
    file = ROOT.TFile(filename)
    file.MakeProject("lib","*","update+");
    file.Close()
elif args.force:
    shutil.rmtree('lib')
    file = ROOT.TFile(filename)
    file.MakeProject("lib","*","update+");
    file.Close()
else:
    print "rqlib already exists. doing nothing."

sys.exit()
