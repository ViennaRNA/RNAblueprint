#! /usr/bin/python

__author__ = "Stefan Hammer"
__copyright__ = "Copyright 2012, TBI"
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hammer@univie.ac.at"
__status__ = "creation"

import subprocess
from optparse import OptionParser

# This is the optionparser where possible options are defined
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="defines input FILE", metavar="FILE")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

# write command line options to variables
verbose = options.verbose
infile = options.filename
# check if infile is empty!

# function to run command line programs in the code
def check_output(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    output = process.communicate()
    retcode = process.poll()
    if retcode:
            raise subprocess.CalledProcessError(retcode, command, output=output[0])
    return output[0]

bin_subopt = "RNAsubopt"
bin_barriers = "../Barriers-install/bin/barriers"
bin_treekin = "../treekin-install/bin/treekin"

#command = bin_subopt + " -s -e 1 < %s" % infile + " | " 
#	+ bin_barriers + " -G RNA-noLP --bsize --max 1000 --rates" + " | " 
#	+ bin_treekin + " -m I --p0 4=1 --ratesfile=%s" % ("rates.out")

subopt = subprocess.Popen([bin_subopt + " -s -e 1 < %s" % infile], stdout=subprocess.PIPE)
barriers = subprocess.Popen([bin_barriers + " -G RNA-noLP --bsize --max 1000 --rates"], stdout=subprocess.PIPE)
#treekin = subprocess.Popen([bin_treekin + " -m I --p0 4=1 --ratesfile=%s" % "rates.out"], stdout=subprocess.PIPE)

subopt.stdout.close()
output = barriers.communicate()[0]


# run RNAsubopt and write output to variable
#subopt_command = " -s -e 1 < %s" % infile
#subopt_out =  check_output(subopt_command)
#print "output is: %s" % subopt_out

