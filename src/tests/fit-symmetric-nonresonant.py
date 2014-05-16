import subprocess

process = subprocess.Popen('../fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate('SymmetricNonresonant\nsymmetric-nonresonant.dat\nnoprint')

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == 'Resid = 2.752039e+01\nc= 5.824e+01, d= 9.698e+00, norm= 1.180e+01\n')
