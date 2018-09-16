#!/usr/bin/env python3

"""
coilprofile.py: A program that makes plots of the response of pickup coils in the VSM@HFML in Nijmegen, The Netherlands

This program is designed to work with Hung's data acquisition software and 
VSM control panel

When a magnetized sample in a VSM goes through a position loop, meaning that 
it will vibrate at a certain position for a set amount of time, stop, and then move on to the next position. A voltage will be induced in the pickup coils, this voltage will depend strongly on the position at which the sample is currently vibrating. When the voltage across the coil is measured during this
loop, there will be distinct peaks representing every position. 
This program will identify and plot these peaks versus the position at which 
they occur.
"""

__author__ = "Andrew Ammerlaan"
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Andrew Ammerlaan"
__email__ = "andrewammerlaan@riseup.net"
__status__ = "Production"

import sys
import os.path
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

## Config

# Set the paramaters that were used in the data acquisition
convert = 0.005376 # How many mm does the motor move 
# when it moves 1 motor position?
stepsize = 50 #0.269mm
extrapolstep = 5 # If fit=yes, how many steps should be extrapolated in plot
endstep = 8000 #43.011mm
beginstep = 3000 #16.129mm, it is assumed that the first peak is here
# if not, then the x-axis will not correspond correctly to the motor 
# coordinates

numstep = 100
stepdur = 20 #s
reacttime = 1 #s It takes a certain amount of time for the voltage to reach
# it's peak value, these data points should not be included in the peak.

disp = 200 #1.08mm
freq = 18.1 #Hz
deltafreq = 0.5 #Hz frequency must be between freq-deltafreq and freq+deltafreq for it to be a valid datapoint

minpeakV = 1e-4 #V signal must be minimally this high for it to be considered a data point

# The last peak can be problematic because the measurement does not stop 
# automaically, to remove the last point say 'yes' here, otherwise say 'no'
rm_lastpoint = 'yes'
rm_firstpoint = 'no' # Remove first point?
fit = 'yes' # Fit a sine through the data
theonetruecenter = 0 # If multiple possible coilcenters exist it is up to
# the user to choose which one is the correct one, e.g. if the first peak
# shoud be the coil center set this to 0
makecoordrel = 'yes' # Plot x coordinates in mm relative to coil center. 
# This option is only available if fit='yes'

# When measuring additional devices as well as the 2 SR830s,
# the index of the columns that correspond to the reference and signal 
# lock in's data should be inserted into their respective variable
# Note that time is always 0
# e.g: "Time(sec)	Dev2_freq	Dev2SR830_X	Dev2_Y	
# Dev3_freq	Dev3SR830_X	Dev3_Y"
sigfreqpos = 1	# Default: 1
sigXpos = 2	# Default: 2
sigYpos = 3	# Default: 3
reffreqpos = 4	# Default: 4
refXpos = 5	# Default: 5
refYpos = 6	# Default: 6

## Program

# Check that a data file is given as argument and set the filename parameter
if len(sys.argv) < 2:
	print("Please parse a data file as first argument")
	sys.exit("Error: Not enough arguments")

filename = sys.argv[1]

if os.path.isfile(filename) == False:
	print("Please parse a valid data file as first argument")
	sys.exit("Error: File not found")

# Open and read the file and save it's lines in the list 'file', then close it
with open(filename, mode='r') as fileopen:
	file = fileopen.readlines()

# Count the number of lines in the file
num_lines = len(file)

print("Loaded file:%s \nThis file has %i lines\n" %(filename,num_lines))

# Create arrays for all data columns
time = np.zeros(num_lines)
sigfreq = np.zeros(num_lines)
sigX = np.zeros(num_lines)
sigY = np.zeros(num_lines)
reffreq = np.zeros(num_lines)
refX = np.zeros(num_lines)
refY = np.zeros(num_lines)

a = 0 # Counting integer

for line in file:
    if line.startswith(('#','\s','\n','\t','Time(sec)')):
        # Skip the header and empty lines
        continue
    # For each line split the columns into dataline and move into data arrays
    dataline = line.split()
    time[a] = dataline[0]
    sigfreq[a] = dataline[sigfreqpos]
    sigX[a] = dataline[sigXpos]
    sigY[a] = dataline[sigYpos]
    reffreq[a] = dataline[reffreqpos]
    refX[a] = dataline[refXpos]
    refY[a] = dataline[refYpos]
    a = a + 1

# Make the amplitude relative to the reference amplitude, this becomes
# relavant near the beginning and end of the position loop where the 
# motor has to work harder against the spring
relsigX = sigX / refX
relsigY = sigY / refY

# Remove all extra zeros
time = np.trim_zeros(time)
sigfreq = np.trim_zeros(sigfreq)
sigX = np.trim_zeros(sigX)
sigY = np.trim_zeros(sigY)
reffreq = np.trim_zeros(reffreq)
refX = np.trim_zeros(refX)
refY = np.trim_zeros(refY)

relsigX = np.trim_zeros(relsigX)
relsigY = np.trim_zeros(relsigY)

# Introduce new arrays for the processed data
datapoint = np.zeros(time.size)
peaktime = np.zeros(time.size)
datastddev = np.zeros(time.size)
tmpindices = np.zeros(time.size)

timeofoneindex = time[1] - time[0] # Moving 1 index corresponds to this time
indexwidth = int(round(reacttime / timeofoneindex)) # How many indices 
# around a point should also be above minpeakV?
peakindexrange = int(round((stepdur - (2 * reacttime)) / timeofoneindex)) #
# How many indices is a peak wide?
peakbit = 0 # Bit to keep track of if a peak has been detected
k = 0 # Counting integer
j = 0 # Counting integer

for i in range(0, time.size - 1):
    # Only use data where the frequencies correspond to approximatly
    # the frequency that was input into the motor
    if freq - deltafreq <= sigfreq[i] <= freq + deltafreq and freq - deltafreq <= reffreq[i] <= freq + deltafreq:
        if abs(sigX[i]) >= minpeakV:
            # If this is true, a peak has been found, indices of peak are 
            # recorded in tmpindices
            tmpindices[j] = i
            j = j + 1
            peakbit = 1 # Found a peak, set bit to 1
        # When this is no longer true the peak has ended
        else:
            # Now process the peak, but do it only once (peakbit is zero 
            # again when this is done)
            if peakbit == 1:              
                tmpindices = np.trim_zeros(tmpindices) # Trim extra zeros
                
                if tmpindices.size <= indexwidth:
                    # The peak is to small, it must be a random fluctuation
                    tmpindices = np.zeros(time.size) # Reset tmp array
                    j = 0 # Reset counting integer
                    peakbit = 0 # Reset peakbit
                    print("A peak was found, however it lasted for less then %f seconds, it will be ignored. If it shouldn't be ignored decrease reacttime = %f s. It might be a good idea to increase minpeakV = %f V" %(reacttime, reacttime, minpeakV))
                    continue
                
                # The peak is not a random fluctuation, where does it begin
                # and where does it end?
                startindex = np.min(tmpindices)
                endindex = np.max(tmpindices)
                
                # Thus the middle of the peak has index:
                midindex = int(round((startindex + endindex) / 2))
                
                p = 0 # Counting integer
                peaksig = np.zeros(time.size) # New tmp array for the 
                # relevant data
                
                for t in range(midindex - peakindexrange , midindex + peakindexrange):
                    # Move relevant data into new array
                    peaksig[p] = relsigX[t]
                    p = p + 1
                
                # Trim zeros, save time index, average and stddev of 
                # peak data
                peaksig = np.trim_zeros(peaksig)
                peaktime[k] = time[midindex] 
                datapoint[k] = np.average(peaksig)
                datastddev[k] = np.std(peaksig)
                
                peakbit = 0 # Reset peakbit
                j = 0 # Reset counting integer
                tmpindices = np.zeros(time.size) # reset tmp array
                k = k +1

# Trim zeros of peak averages and peak times
datapoint = np.trim_zeros(datapoint)
peaktime = np.trim_zeros(peaktime)

# Last point can be problematic because the motor will stay on at the end
# of the position loop
if rm_lastpoint == 'yes' or rm_lastpoint == 'Yes':
    print("\nRemoving final point because rm_lastpoint='yes'")
    datapoint = np.delete(datapoint, -1)
    peaktime = np.delete(peaktime, -1)
    endstep = endstep - stepsize

if rm_firstpoint == 'yes' or rm_firstpoint =='Yes':
    print("\nRemoving first point because rm_firstpoint='yes'")
    datapoint = np.delete(datapoint, 0)
    peaktime = np.delete(peaktime, 0)
    beginstep = beginstep + stepsize

# New array for the erros (stddev)
dataerr = np.empty(datapoint.size)

# Use this for loop instead of trim_zeros because sometimes the last value is
# zero and then it would get removed as well
for l in range(0,datapoint.size-1):
    dataerr[l] = datastddev[l]


# Convert time to position, by calulating the 'speed' at which the peak
# positions change, it is assumed that the first peak occurs at beginstep
speed = stepsize / (peaktime[1] - peaktime[0])
offset = beginstep - (speed * peaktime[0])

motcoord = speed * peaktime + offset

coilcenterbit = 0 # By default, there is no coil center thus bit=0

# length in motcoord to extrapolated
motextrapol = extrapolstep * stepsize

# Fit if it is enabled in config
if fit == 'yes' or fit == 'Yes':
    
    def sine(x, f, amp, phase, offset):
        return amp * np.sin(x * f + phase) + offset

    # An educated guess to the fit paramaters
    guess_amp = (np.max(datapoint) - np.min(datapoint)) / 2
    guess_offset = np.max(datapoint) - guess_amp
    guess_freq = 1 / (abs(motcoord[np.argmax(datapoint)] - motcoord[np.argmin(datapoint)]) * 2)
    p0 = [guess_freq, guess_amp, 0, guess_offset]
    
    # Fit
    params, params_covar = curve_fit(sine, motcoord, datapoint, p0=p0)
    
    # Make points out of fit
    x = np.linspace(motcoord[0] - motextrapol, motcoord[-1] + motextrapol, num=1000)
    data_fit = sine(x, params[0], params[1], params[2], params[3])
    
    print("\nFit paramaters:")
    print("\nfrequency = \t%f \t covariance = \t" %(params[0]))
    print(params_covar[0])
    print("\namplitude = \t%f\t covariance = \t" %(params[1])) 
    print(params_covar[1])
    print("\nphase = \t%f\t covariance = \t" %(params[2]))
    print(params_covar[2])
    print("\noffset = \t%f\t covariance = \t" %(params[3]))
    print(params_covar[3])

    # Find indices where the sign between data_fit and offset changes 
    # (when the sine intersects with the offset)
    # These points (if they exist) must be the coil center
    coilcenterindex = np.argwhere(np.diff(np.sign(params[3] - data_fit)))
    if coilcenterindex.size < 1:
        print("\nCould not find coil center. This might be due to a incorrect fit, or the center is not in the data. Try increasing 'extrapolstep'")
        
    elif coilcenterindex.size > 1:
        allcoilcenter = x[coilcenterindex]
        print("\nFound more then 1 possible coil centers, possible candidates:")
        print(allcoilcenter)
        print("This might be due to a incorrect fit, or the data is extrapolated too far. Try decreasing 'extrapolstep'")
        coilcenter = allcoilcenter[theonetruecenter][0] # Let user decide
        # which center is the right one
        print("\nBecause 'theonetruecenter'= %i center-candidate number %i is choosen to be the correct one" %(theonetruecenter, theonetruecenter + 1))
        coilcenterbit = 1 # Found a coilcenter
        
    elif coilcenterindex.size == 1:
        coilcenter = x[coilcenterindex][0][0]
        print("\nFound coilcenter at motor coordinates: %f" %(coilcenter))
        coilcenterbit = 1 # Found a coilcenter
        
else:
    makecoordrel = 'no'
     
# Motor coords to mm
mmcoord = motcoord * convert
mmextrapol = motextrapol * convert
mmdisp = disp * convert
# Make the mm coords relative to coilcenter, if enabled
if makecoordrel == 'yes' or makecoordrel == 'Yes':
    if coilcenterbit == 1:
        coilcenter_in_mm = coilcenter * convert
        mmcoord = mmcoord - coilcenter_in_mm
        print("\n'makecoordrel' is enabeld, and a coil center has been found. Top x-axis coordinates will be relative to the coil center")

    else:
        print("\n'makecoordrel' is enabeld, but a coil center could not be found. Coordinates will not be relative")

# Plot
fig, ax1 = plt.subplots()

ax1.minorticks_on()
ax1.errorbar(motcoord, datapoint, xerr=disp, yerr=dataerr, linestyle='None')

ax1.set_ylabel("Relative peak voltage sigX / refX")
ax1.set_xlabel("Motor position")

ax1.set_ylim(np.min(datapoint) - minpeakV, np.max(datapoint) + minpeakV)
ax1.set_xlim(motcoord[0] - motextrapol, motcoord[-1] + motextrapol)

if fit == 'yes' or fit == 'Yes':
    ax1.plot(x, data_fit)

ax2 = ax1.twiny() # Second x-axis is twin to first

# Make second plot on the second axis, but make it transparent
# because otherwise the plot would be shown twice
ax2.minorticks_on()
ax2.errorbar(mmcoord, datapoint, xerr=mmdisp, yerr=dataerr, linestyle='None', color='white', alpha=0)

if makecoordrel == 'yes' or makecoordrel == 'Yes':
    ax2.set_xlabel("Sample position relative to the coil center (mm)")
else:
    ax2.set_xlabel("Sample position (mm)")
ax2.set_xlim(mmcoord[0] - mmextrapol, mmcoord[-1] + mmextrapol)

if (makecoordrel == 'yes' or makecoordrel == 'Yes') and coilcenterbit == 1:
    # If a center has been found, add lines to show it
    ax2.axvline(color='black', alpha=0.5)
    ax2.axhline(y=params[3] ,color='black', alpha=0.5)


ax1.grid(which='both', axis='both')

if len(sys.argv) == 3:
    # If there is a second argument, interpret it as output file
    fileout = sys.argv[2]
    if fileout.endswith(".tex"):
        # If the file ends in .tex, use tikz to save it
        from matplotlib2tikz import save as tikz_save
        tikz_save(fileout, figureheight = '\\figH', figurewidth = '\\figW')
        print("\nSaving plot as %s using matplotlib2tikz" %(fileout))
        sys.exit(0)
    else:
        # If not use regular save function
        plt.savefig(fileout)
        print("\nSaving plot as %s" %(fileout))
        sys.exit(0)

else:
    # No output argument, so just show the plot instead
    plt.show()
    print("\nShowing plot, to save the plot parse a output file as second argument. Valid outut files are files that are supported by matplotlib's savefig (e.g. .pdf, .png) or matplotlib2tikz (.tex)")
    sys.exit(0)
