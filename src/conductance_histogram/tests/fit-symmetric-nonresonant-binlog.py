import subprocess

# test the log-binning aspect of the fitter program

process = subprocess.Popen('../fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate( \
'SymmetricNonresonant\n' \
'symmetric-nonresonant-binlog.dat\n' \
'bin log 10.')

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == 'Resid = 2.752118e+01\nc=5.8240e+01, d=9.6983e+00, norm=1.1800e+01\n')
