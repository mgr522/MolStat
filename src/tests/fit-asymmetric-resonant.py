import subprocess

process = subprocess.Popen('../fitter', stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
output = process.communicate('AsymmetricResonant\nasymmetric-resonant.dat\nnoprint')

# make sure no errors were reported
assert(output[1] == '')

# check the output string
assert(output[0] == 'Resid = 4.895159e+00\ngammaL=1.7085e+01, gammaR=1.1731e+01, r=2.5489e+01, norm=8.2639e+00\n')
