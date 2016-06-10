import subprocess
import os
import math
import numpy as np

# This script simulates a 2D conductance-distance histogram based on image charge effects
# and the symmetric one site channel model in MolStat.
# This simulation process works by simulating single traces using MolStat
# These single traces are then histogrammed

# To simulate a single trace, a file for the trace number is created
# Then, a single trial input is run in molstat for each distance point
# The output from MolStat will be a single conductance value attributed to each distance point
# The outputs are appended to the end of the single trace file

# We are assuming that the electrodes are modeled as infinite flat planes of gold atoms
# All interelectrode distances measured with respect to the electrode surface (nucleus)

z_e = .3 # interelctrode distance immediately after junction rupture
z_m = .5 # interelectrode distance when the molecule bridges the junction
z2 = z_m + 1 # interelectrode distance at the end of the trace
z_ip = .1 # distance from gold surface to image plane (THIS IS NOT INTERELECTRODE DISTANCE)
n = 300 # number of distance points
nt = 300 # number of traces
ke = 8.98e9 * 1.602e-19*1e9 # Coulomb's constant in units of nm*ev/e^2, where e is the charge of an electron
epsilon = [0 for i in range(n)] # initializing new site energy array

z = np.linspace(z_m, z2, num = n) # initialize distance vector with n points

# Computing image charge effects
IC = ke * (1/(4*np.fabs(z/2-z_ip)) + 1/(4*np.fabs(z/2-z_ip)))
#print IC

# Computing image charge renormalization relative to molecule in junction
epsilon_renorm = IC - IC[0]
z_eff = z - z[0] # distance 0 is set to the point where the molecule first attaches


for j in xrange(nt): 
    l = np.random.normal(.3,.075) # accounts for individual trace length variations
    nm = np.random.poisson(1) # allows for multiple molecules to bind in junction per trace
    for i in xrange(n):   
        if z_eff[i] <= l: # molecular regime
            # Open input and output files
            fin = open('/Users/benwu/Documents/molstat-1.3.1/src/MolecularInput.txt',"r")
            fout = open('/Users/benwu/Documents/molstat-1.3.1/src/input.txt',"w")
    
            for line in fin:
                # updates inputfile
                if "output" in line: # searchs input file for word "output"
                    data = line.split() # parses up line on which "output" is found
                    fout.write(line.replace(data[1], 'MolecularContribution'+str(j) + '.dat')) #replaces "output" with updated file output 
                elif "epsilon" in line: 
                    data = line.split() 
                    epsilon[i] = float(data[3])
                    fout.write(line.replace(data[3],str(epsilon[0] + epsilon_renorm[i]))) 
                elif "width" in line: 
                    data = line.split() 
                    fout.write(line.replace(data[3],str(z[i]-z_e)))
                elif "nm" in line:
                    data = line.split() # parses up line
                    fout.write(line.replace(data[3],str(nm)))
                else:
                    fout.write(line) # rewrites the rest of input file line into output file
            # close input and output files
            fin.close()
            fout.close()
            # runs molstat simulator with input.txt
            process = subprocess.Popen('../molstat-simulator', \
                    stdout=subprocess.PIPE, \
                    stdin=subprocess.PIPE, \
                    stderr=subprocess.PIPE)
            with open('/Users/benwu/Documents/molstat-1.3.1/src/input.txt', "r") as inputfile:
                text = inputfile.read()
            #print(text)
            output = process.communicate(text)
            #print(output)
    
        else: # tunneling regime
            fin = open('/Users/benwu/Documents/molstat-1.3.1/src/TunnelingInput.txt',"r")
            fout = open('/Users/benwu/Documents/molstat-1.3.1/src/input.txt',"w")
    
            for line in fin:
                if "output" in line:
                    data = line.split()
                    fout.write(line.replace(data[1], 'TunnelingContribution'+str(j) + '.dat'))
                elif "width" in line: # searching for epsilon
                    data = line.split() # parses up line
                    fout.write(line.replace(data[3],str(z[i]-z_e))) # writes the updated epsilon value into the output file
                else:
                    fout.write(line) # rewrites input file line into output file
            fin.close()
            fout.close()        
            process = subprocess.Popen('../molstat-simulator', \
                    stdout=subprocess.PIPE, \
                    stdin=subprocess.PIPE, \
                    stderr=subprocess.PIPE)
            with open('/Users/benwu/Documents/molstat-1.3.1/src/input.txt', "r") as inputfile:
                text = inputfile.read()
            #print(text)
            output = process.communicate(text)
            #print(output)
    
