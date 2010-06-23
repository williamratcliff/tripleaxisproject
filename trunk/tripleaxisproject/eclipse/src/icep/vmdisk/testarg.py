import argparse


arg='myseq1 myseq2 -s -f -now '
parser = argparse.ArgumentParser(description='Frabble the foo and the bars')
parser.add_argument('-f', action='store_true', help='frabble the foos')
parser.add_argument('-s', action='store_true', help='frabble the foos')
parser.add_argument('-o', nargs=1, help='frabble the foos')
parser.add_argument('-now', action='store_true', help='frabble the foos')
parser.add_argument('scans', nargs='+', help='a bar to be frabbled')
myargs = parser.parse_args(arg.split())
isFile=myargs.f
isSeconds=myargs.s
isNow=myargs.now
scans=myargs.scans
print 'isFile',isFile
print 'isSeconds',isSeconds
print 'isNow',isNow
print 'scans',scans
print 'overhead',myargs.o