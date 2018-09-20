#!/usr/bin/env python3

"""
coilprofile.py: A program that makes plots of the response of pickup coils
in the VSM@HFML in Nijmegen, The Netherlands

This program is designed to work with Hung's data acquisition software and
VSM control panel

When a magnetized sample in a VSM goes through a position loop, meaning that
it will vibrate at a certain position for a set amount of time, stop,
and then move on to the next position.
A voltage will be induced in the pickup coils, this voltage will depend
strongly on the position at which the sample is currently vibrating.
When the voltage across the coil is measured during this
loop, there will be distinct peaks representing every position.
This program will identify and plot these peaks versus the position at which
they occur.
"""

__author__ = "Andrew Ammerlaan"
__license__ = "GPLv3"
__version__ = "2.1"
__maintainer__ = "Andrew Ammerlaan"
__email__ = "andrewammerlaan@riseup.net"
__status__ = "Production"

import sys
import os.path
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

"""Config"""

convert = 5.376e-3  # How many mm does the motor move
# when it moves 1 motor position?
mmextrapol = 2  # (mm) How many mm to extrapolate/to add to xlim
yextrapol = 0.002  # How much to add to ylim

stepdur = 20  # (s) # How long is the VSM at 1 position?
reacttime = 5  # (s) It takes a certain amount of time for the voltage to reach
# it's peak value, these data points should not be included in the peak.

mmdisp = 1.08  # (mm) How much does the VSM vibrate at a position? (this is
# NOT the amplitude, but two times the amplitude)
freq = 18.1  # (Hz)
deltafreq = 0.3  # (Hz) frequency must be between freq-deltafreq and
# freq+deltafreq for it to be a valid datapoint

minpeakV = 0.4e-4  # (V) signal must be minimally this high for it to be
# considered a data point

# Potentio meter calibaration
# motor coordinates = a * potentiometer DC voltage + b
a = 8830.788615
b = 238.4372817

# The last peaks can be problematic because the measurement does not stop
# automaically at the end of the position loop,
# set this to remove a number of peaks from the end
rm_lastpoint = 0
rm_firstpoint = 0  # Amount of peaks to remove from beginning

fit = 'no'  # Fit theoretical coil profile to the data
guess_x0 = 15  # mm where is the coil center approximatly?
guess_r1 = 1.6  # mm approximate inner radius
guess_r2 = 2.8  # mm approximate outer radius
guess_L = 12  # mm approximate length
guess_N = 6120  # aproximate number of windings

makecoordrel = 'yes'  # Plot x coordinates in mm relative to coil center.
# This option is only available if fit='yes'
dividebyref = 'yes'  # Divide sigX by refX to make the data relative to the
# reference signal from recorded from the VSM motor

# When measuring additional devices as well as the 2 SR830s,
# the index of the columns that correspond to the reference and signal
# lock in's data should be inserted into their respective variable
# e.g: "Time(sec)	Dev2_freq	Dev2SR830_X	Dev2_Y	Dev3_freq
# Dev3SR830_X	Dev3_Y	Dev4_K2000(V)"
timepos = 0  # Default 0
sigfreqpos = 1  # Default: 1
sigXpos = 2  # Default: 2
sigYpos = 3  # Default: 3
reffreqpos = 4  # Default: 4
refXpos = 5  # Default: 5
refYpos = 6  # Default: 6
potmetDCpos = 7  # Default 7

timeinmillisec = 'yes'  # set to yes if time in data file is in milliseconds 
# instead of seconds


"""Read file"""

# Check that a data file is given as argument and set the filename parameter
if len(sys.argv) < 2:
    print("Please parse a data file as first argument")
    sys.exit("Error: Not enough arguments")

filename = sys.argv[1]

if os.path.isfile(filename) is False:
    print("Please parse a valid data file as first argument")
    sys.exit("Error: File not found")

# Open and read the file and save it's lines in the list 'file', then close it
with open(filename, mode='r') as fileopen:
    file = fileopen.readlines()

# Count the number of lines in the file
num_lines = len(file)

print("Loaded file: %s \nThis file has %i lines\n" % (filename, num_lines))


"""Sort data from file into arrays"""

# Create arrays for all data columns
time = np.zeros(num_lines)
sigfreq = np.zeros(num_lines)
sigX = np.zeros(num_lines)
sigY = np.zeros(num_lines)
reffreq = np.zeros(num_lines)
refX = np.zeros(num_lines)
refY = np.zeros(num_lines)
potmetDC = np.zeros(num_lines)

h = 0  # Counting integer

for line in file:
    # Some data files start there lines with a space or tab, this would cause
    # isdigit() to report false. This finds the first index that is not a space
    w = 0
    for q in range(0, len(line)-1):
        if line[q].isspace():
            w = w + 1
        else:
            break

    # Now check if the first non-space entry is a number, else it cannot be
    # a data line
    if line[w].isdigit():
        # For each line split the columns into dataline and move data to arrays
        dataline = line.split()
        time[h] = dataline[timepos]
        sigfreq[h] = dataline[sigfreqpos]
        sigX[h] = dataline[sigXpos]
        sigY[h] = dataline[sigYpos]
        reffreq[h] = dataline[reffreqpos]
        refX[h] = dataline[refXpos]
        refY[h] = dataline[refYpos]
        potmetDC[h] = dataline[potmetDCpos]
        h = h + 1

# Make the amplitude relative to the reference amplitude, this becomes
# relavant near the beginning and end of the position loop where the
# motor has to work harder against the spring
if dividebyref == 'yes' or dividebyref == 'Yes':
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
potmetDC = np.trim_zeros(potmetDC)

if dividebyref == 'yes' or dividebyref == 'Yes':
    relsigX[np.isneginf(relsigX)] = 0
    relsigY[np.isneginf(relsigY)] = 0
    relsigX[np.isposinf(relsigX)] = 0
    relsigY[np.isposinf(relsigY)] = 0
    relsigX = np.nan_to_num(relsigX)
    relsigY = np.nan_to_num(relsigY)
    relsigX = np.trim_zeros(relsigX)
    relsigY = np.trim_zeros(relsigY)


"""Find peaks"""

# Introduce new arrays for the processed data
datapoint = np.zeros(time.size)
peakpos = np.zeros(time.size)
datastddev = np.zeros(time.size)
tmpindices = np.zeros(time.size)

if timeinmillisec == 'yes' or timeinmillisec == 'Yes':
    time = time * 1e-3

timeofoneindex = np.average(np.diff(time))  # Moving 1 index corresponds
# to this much in time
indexwidth = int(round(reacttime / timeofoneindex))  # How many indices
# around a point should also be above minpeakV?
peakindexrange = int(round((stepdur - (2 * reacttime)) / timeofoneindex))
# How many indices is a peak wide?

peakbit = 0  # Bit to keep track of if a peak has been detected
k = 0  # Counting integer
j = 0  # Counting integer

# Find peaks
for i in range(0, time.size - 1):
    # Only use data where the frequencies correspond to approximatly
    # the frequency that was input into the motor
    if (freq - deltafreq <= sigfreq[i] <= freq + deltafreq and
       freq - deltafreq <= reffreq[i] <= freq + deltafreq and
       abs(sigX[i]) >= minpeakV):
        # If this is true, a peak has been found, indices of peak are
        # recorded in tmpindices
        tmpindices[j] = i
        j = j + 1
        peakbit = 1  # Found a peak, set bit to 1
    # When this is no longer true the peak has ended
    else:
        # Now process the peak, but do it only once (peakbit is zero
        # again when this is done)
        if peakbit == 1:
            tmpindices = np.trim_zeros(tmpindices)  # Trim extra zeros

            if tmpindices.size <= indexwidth:
                # The peak is too small, it must be a random fluctuation
                tmpindices = np.zeros(time.size)  # Reset tmp array
                j = 0  # Reset counting integer
                peakbit = 0  # Reset peakbit
                print("\nA peak was found, however it lasted for less",
                      "then %f seconds, it will be ignored." % (reacttime),
                      "If it shouldn't be ignored decrease",
                      "reacttime = %f s." % (reacttime),
                      "It might be a good idea to",
                      "increase minpeakV = %f V" % (minpeakV))
                continue

            # The peak is not a random fluctuation, where does it begin
            # and where does it end?
            startindex = np.min(tmpindices)
            endindex = np.max(tmpindices)

            # Thus the middle of the peak has index:
            midindex = int(round((startindex + endindex) / 2))

            p = 0  # Counting integer
            peaksig = np.zeros(time.size)  # New tmp array for the
            # relevant data

            for t in range(midindex - peakindexrange,
                           midindex + peakindexrange):
                # Move relevant data into new array
                if dividebyref == 'yes' or dividebyref == 'Yes':
                    peaksig[p] = relsigX[t]
                else:
                    peaksig[p] = sigX[t]
                p = p + 1

            # Trim zeros, save time index, average and stddev of
            # peak data
            peaksig = np.trim_zeros(peaksig)
            peakpos[k] = a * potmetDC[midindex] + b
            datapoint[k] = np.average(peaksig)
            datastddev[k] = np.std(peaksig)

            peakbit = 0  # Reset peakbit
            j = 0  # Reset counting integer
            tmpindices = np.zeros(time.size)  # reset tmp array
            k = k + 1


# Trim zeros of peak averages and peak times
datapoint = np.trim_zeros(datapoint)
datastddev = np.trim_zeros(datastddev)
peakpos = np.trim_zeros(peakpos)

# Remove begining points if necessary
for u in range(0, rm_firstpoint):
    datapoint = np.delete(datapoint, 0)
    datastddev = np.delete(datastddev, 0)
    peakpos = np.delete(peakpos, 0)
    print("\nRemoving peak %i" % (u),
          " because rm_firstpoint=%i" % (rm_firstpoint))


# Remove endpoints, if necessary
for t in range(peakpos.size - 1, (peakpos.size - 1) - rm_lastpoint, -1):
    datapoint = np.delete(datapoint, -1)
    datastddev = np.delete(datastddev, -1)
    peakpos = np.delete(peakpos, -1)
    print("\nRemoving peak %i" % (t),
          "because rm_lastpoint=%i" % (rm_lastpoint))


# Motor coords to millimeter and mm units from config to motor coords.
mmcoord = peakpos * convert
posextrapol = mmextrapol / convert
disp = mmdisp / convert


"""Fit"""

if fit == 'yes' or fit == 'Yes':

    def func(x, x0, A, r1, r2, L):
        dBdz = (
            np.log(r2 + np.sqrt(r2**2 + (x - x0 - L)**2)) -
            np.log(r1 + np.sqrt(r1**2 + (x - x0 - L)**2)) +
            np.log(r1 + np.sqrt(r1**2 + (x - x0)**2)) -
            np.log(r2 + np.sqrt(r2**2 + (x - x0)**2)) +
            (x - x0 - L)**2 * (
                          1 / (r2 * np.sqrt(r2**2 + (x - x0 - L)**2) +
                               r2**2 + (x - x0 - L)**2) -
                          1 / (r1 * np.sqrt(r1**2 + (x - x0 - L)**2) +
                               r1**2 + (x - x0 - L)**2)) +
            (x - x0)**2 * (
                    1 / (r1 * np.sqrt(r1**2 + (x - x0)**2) +
                         r1**2 + (x - x0)**2) -
                    1 / (r2 * np.sqrt(r2**2 + (x - x0)**2) +
                         r2**2 + (x - x0)**2))
                )
        V = A * dBdz
        return V

    # An educated guess to the fit paramaters
    if dividebyref == 'yes' or dividebyref == 'Yes':
        guess_amp = 1
    else:
        guess_amp = np.sqrt(2) * np.pi**2 * mmdisp * freq * guess_N * 1e-7
    p0 = [guess_x0, guess_amp, guess_r1, guess_r2, guess_L]

    # Fit
    params, params_covar = curve_fit(func, mmcoord, datapoint, p0=p0,
                                     sigma=datastddev, absolute_sigma=True)

    # Calculate stddev from covar
    params_err = np.sqrt(np.diag(params_covar))

    # Make points out of fit
    fit_pos = np.linspace(mmcoord[0] - mmextrapol,
                          mmcoord[-1] + mmextrapol, num=1000)
    data_fit = func(fit_pos, params[0], params[1], params[2],
                    params[3], params[4])

    print("\nFit paramaters:")
    print("x-offset = \t%f \u00B1 %f mm" % (params[0], params_err[0]))
    print("amplitude = \t%f \u00B1 %f" % (params[1], params_err[1]))
    print("inner radius = \t%f \u00B1 %f mm" % (params[2], params_err[2]))
    print("outer radius = \t%f \u00B1 %f mm" % (params[3], params_err[3]))
    print("length = \t%f \u00B1 %f mm" % (params[4], params_err[4]))

    # if enabled, move the x-axis such that the coilcenter is zero
    if (makecoordrel == 'yes' or makecoordrel == 'Yes'):
        coilcenterindex = np.argwhere(np.diff(np.sign(data_fit)))
        if coilcenterindex.size < 1:
            print("\nCan't find coil center, this is most likly due",
                  "to a incorrect fit. Try increasing 'deltafreq' or",
                  "'minpeakV'")
            makecoordrel = 'no'
        else:
            offset = fit_pos[coilcenterindex[0][0]]
            print("\nCoil center found at: %f =" % (offset / convert),
                  "%f mm, x-axis will be shifted accordingly" % (offset))
            mmcoord = mmcoord - offset
            fit_pos = np.linspace(mmcoord[0] - mmextrapol,
                                  mmcoord[-1] + mmextrapol, num=1000)
            data_fit = func(fit_pos, params[0] - offset, params[1], params[2],
                            params[3], params[4])
else:
    makecoordrel = 'no'


"""Plot"""

fig, ax1 = plt.subplots()

ax1.minorticks_on()
ax1.errorbar(mmcoord, datapoint, xerr=(mmdisp / 2), yerr=datastddev,
             linestyle='None', fmt='.', elinewidth=0.5, )

if ((makecoordrel == 'yes' or makecoordrel == 'Yes') and
   (fit == 'yes' or fit == 'Yes')):
    ax1.set_xlabel("Sample position relative to the coil center (mm)")
else:
    ax1.set_xlabel("Sample position (mm)")

if dividebyref == 'yes' or dividebyref == 'Yes':
    ax1.set_ylabel("Relative peak voltage sigX / refX")
else:
    ax1.set_ylabel("Peak signal voltage sigX (V)")

ax1.set_xlim(mmcoord[0] - mmextrapol, mmcoord[-1] + mmextrapol)

if fit == 'yes' or fit == 'Yes':
    ax1.plot(fit_pos, data_fit)
    ax1.set_ylim(min(np.min(datapoint) - yextrapol,
                     np.min(data_fit) - yextrapol),
                 max(np.max(datapoint) + yextrapol,
                     np.max(data_fit) + yextrapol))
else:
    ax1.set_ylim(np.min(datapoint) - yextrapol, np.max(datapoint) + yextrapol)

ax1.axhline(color='black', alpha=0.5)
if (makecoordrel == 'yes' or makecoordrel == 'Yes'):
    # If a center has been found, add lines to show it
    ax1.axvline(color='black', alpha=0.5)

ax2 = ax1.twiny()  # Second x-axis is twin to first

# Make second plot on the second axis, but make it transparent
# because otherwise the plot would be shown twice
ax2.minorticks_on()
ax2.errorbar(peakpos, datapoint, xerr=(disp / 2),
             yerr=datastddev, linestyle='None', color='white', alpha=0)

ax2.set_xlabel("Motor position")
ax2.set_xlim(peakpos[0] - posextrapol, peakpos[-1] + posextrapol)

ax1.grid(which='both', axis='both')
# Sometimes axis labels will fall of the plot without the tight_layout line
plt.tight_layout()


"""Save"""

if len(sys.argv) == 3:
    # If there is a second argument, interpret it as output file
    fileout = sys.argv[2]
    if fileout.endswith(".tex"):
        # If the file ends in .tex, use tikz to save it
        from matplotlib2tikz import save as tikz_save
        tikz_save(fileout, figureheight='\\figH', figurewidth='\\figW')
        print("\nSaving plot as %s using matplotlib2tikz" % (fileout))
        sys.exit(0)
    else:
        # If not use regular save function
        plt.savefig(fileout)
        print("\nSaving plot as %s" % (fileout))
        sys.exit(0)

else:
    # No output argument, so just show the plot instead
    plt.show()
    print("\nShowing plot, to save the plot parse a output file as",
          "second argument. Valid outut files are files that are supported",
          "by matplotlib's savefig (e.g. .pdf, .png)",
          "or matplotlib2tikz (.tex)")
    sys.exit(0)
