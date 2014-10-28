# This file is a part of MolStat, which is distributed under the Creative
# Commons Attribution-NonCommercial 4.0 International Public License.
# MolStat (c) 2014, Northwestern University.

##
 # @file tests/simulator-dists.py
 # @brief Use the \"identity\" observable with a large number of trials to
 #    determine if it's working correctly.
 # 
 # @test Test suite for the simulator program.
 #
 # @author Matthew G.\ Reuter
 # @date October 2014

import subprocess
import os
import math

## @cond

# runtime variables
trials = 1000000
datfile = 'hist.dat'
bins = 10

# test 1 -- constant distribution
print "Constant distribution"
constval = 5.
process = subprocess.Popen('../molstat-simulator', \
	stdout=subprocess.PIPE, \
	stdin=subprocess.PIPE, \
	stderr=subprocess.PIPE)
output = process.communicate( \
'observable Identity ' + str(bins) + ' linear\n' \
'model IdentityModel\n' \
'	distribution parameter constant ' + str(constval) + '\n' \
'endmodel\n' \
'trials ' + str(trials) + '\n' \
'output ' + datfile)

# even though we requested 'bins' bins, the exception handler reduces it to 1
# because constant distribution always produces (in this case) 'constval'
# (i.e., there is a null data range)

# read in the histogram
hist = open(datfile, 'r')
j = 0 # line count
for bin in hist:
	# make sure that the output is correct. first, tokenize the string
	tokens = str.split(bin)
	assert(len(tokens) == 2)
	
	x = float(tokens[0])
	pdf = float(tokens[1])

	assert(math.fabs(x - constval) < 1.e-6)
	assert(math.fabs(pdf - float(trials)) < 1.e-6 * float(trials))
	j += 1

hist.close()
assert(j == 1)
# delete the output histogram file
os.remove(datfile)



# test 2 -- uniform distribution
print "\nUniform distribution"
minval = -2.
maxval = 2.
process = subprocess.Popen('../molstat-simulator', \
	stdout=subprocess.PIPE, \
	stdin=subprocess.PIPE, \
	stderr=subprocess.PIPE)
output = process.communicate( \
'observable Identity ' + str(bins) + ' linear\n' \
'model IdentityModel\n' \
'	distribution parameter uniform ' + str(minval) + ' ' + str(maxval) + '\n' \
'endmodel\n' \
'trials ' + str(trials) + '\n' \
'output ' + datfile)

# read in the histogram
hist = open(datfile, 'r')
j = 0 # line count
for bin in hist:
	# make sure that the output is correct. first, tokenize the string
	tokens = str.split(bin)
	assert(len(tokens) == 2)
	
	x = float(tokens[0])
	pdf = float(tokens[1])

	assert(math.fabs(x - (minval + (j + 0.5) * (maxval - minval)/bins)) < 1.e-2)
	print "Expected: " + str(float(trials)/bins) + ", Actual: " + str(pdf)
	assert(math.fabs(pdf - float(trials)/bins) < 1.e-2 * float(trials))
	j += 1

hist.close()
assert(j == bins)
# delete the output histogram file
os.remove(datfile)



# test 3 -- normal distribution
print "\nNormal distribution"
mean = 1.
stdev = 2.
bins = 17
process = subprocess.Popen('../molstat-simulator', \
	stdout=subprocess.PIPE, \
	stdin=subprocess.PIPE, \
	stderr=subprocess.PIPE)
output = process.communicate( \
'observable Identity ' + str(bins) + ' linear\n' \
'model IdentityModel\n' \
'	distribution parameter normal ' + str(mean) + ' ' + str(stdev) + '\n' \
'endmodel\n' \
'trials ' + str(trials) + '\n' \
'output ' + datfile)

# read in the histogram
hist = open(datfile, 'r')
x = []
pdf = []
for bin in hist:
	# make sure that the output is correct. first, tokenize the string
	tokens = str.split(bin)
	assert(len(tokens) == 2)
	
	x.append( float(tokens[0]) )
	pdf.append( float(tokens[1]) )

hist.close()
# delete the output histogram file
os.remove(datfile)
assert(len(x) == bins)

# make sure the distribution is correct
dx = 0.5 * (x[1] - x[0])
for j in range(len(x)):
	expected = trials * \
		(0.5 * (1. + math.erf((x[j] + dx - mean) / (math.sqrt(2.) * stdev))) \
		-0.5 * (1. + math.erf((x[j] - dx - mean) / (math.sqrt(2.) * stdev)))
	)

	print "Expected: " + str(expected) + ", Actual: " + str(pdf[j])
	assert(math.fabs(pdf[j] - expected) < 1.e-2 * float(trials))

## @endcond
