import subprocess

process = subprocess.Popen('../fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate('SymmetricResonant\nsymmetric-resonant.dat\nnoprint')

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == 'Resid = 4.310146e-01\ngamma= 9.946e+00, norm= 3.957e+00\n')