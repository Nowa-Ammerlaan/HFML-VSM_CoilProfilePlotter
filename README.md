# HFML-VSM pickupcoil profile plotter
A python program to make plots of the response of a VSM's pickup coils as a function of the position.

### Note:
I assume that you are already familiar with the concept of a VSM (Vibrating-sample magnetometer), if not then this is probably not what you are looking for.

## Long discription
To measure the response of a VSM's pickup coils, set it up to go through a 'position loop'. This means that the VSM will start vibrating for a set amount of time (stepdur) at a set begin position (beginstep), then stop. And then move (stepsize) and start vibrating there again. The vibrating motion will induce a (~90 degree phase shifted) AC voltage across the pickup coils, which will depend strongly on the position at which the VSM is currently vibrating. When this voltage is measured (e.g using 2 Lock-In amplifiers; 1 for reference and 1 for the VSM's signal) it can be plotted against the VSM's position. This program identifies each peak in the signal, calculates it's average value and standard deviation, and then plots is against the position. Optionally it can also fit a sine through the data, and identify the center of the coil (if it was within the range of the measured data). If a center has been found, it also has an option to shift the x-axis such that 0 is in the center of the coil.

## Dependencies
The versions indicated here are the versions that I am currently using, but there is no reason other versions shouldn't work to.
- python-3.6.6
- matplotlib-2.2.2
- numpy-1.14.5
- scipy-1.1.0
- Optionally: matplotlib2tikz-0.6.18 (when saving as .tex)

## How does it work?
Open the coilprofile.py file and edit the config options. See also the overview of options below.
Then run it like this:
```
./coilprofile.py inputfile
```
And it will process the data in `inputfile` and show you the plot using plt.show()
 
Or optionally use a second argument like this:
```
./coilprofile.py inputfile outputfile
```
And the program will use plt.savefig() to save the plot to `outputfile` instead.
Note that `outputfile` should have one of the filename extensions supported by matplotlib (e.g. `.pdf`, `.png`)
or the `.tex` extension, in which case it will use matplotlib2tikz instead.

## Can I use/adapt your program?
Of course you can, it's licensed under GPLv3.

## I have an improvemnt/comment
Awesome, please leave it in the bug reports and I will look at it ASAP.

## My plot looks like crap, what can I do?
Try increasing 'reacttime' and maybe play a bit with 'minpeakV' and 'deltafreq', see the overview of config options below.

## Help, it doesn't work
If you run into anything, and you can't work it out with the comments in the program or one of the options. You can always leave a bug report here and I'll see what I can do.

## Overview of options in the config (the config is just the top of the program)

**convert = _float_** 
How many millimeters does the VSM motor move when it moves 1 in motor coordinates.

**stepsize = _float_**
How far does the motor move between points (in motor coordinates)?

**extrapolstep = _int_**
If fit=yes, how many steps should be extrapolated in plot (in steps, not in motor coordinates)
I recommend you start with this set to 0 and then see how that looks first.

**endstep = _float_**
Where does the position loop end (in motor coordinates)?
This is used to remove points beyond endstep if 'remove_excess' is enabled

**beginstep = _float_**
Where does the position loop begin (in motor coordinates)?
It is assumed the first peak is going to be here, if there was no peak at the first position set this to wherever the first peak did occur.

**numstep = _int_**
How many steps where there in the position loop?
This is not used at the moment.

**stepdur = _float_**
How long does the VSM vibrate at one position (in seconds)? 
This is used in combination with reacttime to calculate how many points there are within a peak.

**reacttime = _float_**
The voltage will not jump instantly from background noise level to a peak, the data points that are above the threshold but not in the peak should be ignored. Set this according to how long it takes for the system to respond, too high is better then too low (in seconds).

**disp = _float_**
The displacement, or two times the amplitude of the VSM's vibration, this is used as the error in the position (in motor coordinates).

**freq = _float_**
The frequency of the VSM's vibration, this is compared with the frequency measured by the lock-in (in Hz).

**deltafreq = _float_** 
How close should the measured frequency be to the variable 'freq'?

**minpeakV = _float_** 
How high should the signal be before it is considered a peak?

**rm_lastpoint = _'yes'_**
Sometimes the last point is off, say yes here to remove it from the data
This happens if the VSM does not stop automatically after the position loop has been completed, but instead continues vibrating at that position until aborted by the user.
If there is more then 1 of such point, try remove_excess instead.

**rm_firstpoint = _'no'_**
Sometimes the first point is off, say yes here to remove it from the data.

**remove_excess = _'yes'_**
Remove points for which their motor coordinates are larger then stepsize, because these points must be recorded after the position loop has ended.

**fit = _'yes'_**
Fit a sine to the data, this only works well if at least both ends of the coil were within the range of the position loop. Otherwise there is no way for the fit to get a correct amplitude of frequency.

**theonetruecenter = _int_**
Because the fit function is sinusoidal it might have multiple points where it crosses it's rest position, but *there can be only one* coil center, here the user can select which one is the one true coil center (e.g. 0 means the first candidate is the right one).

**makecoordrel = _'yes'_**
Should the upper x-axis coordinates be shifted such that 0 corresponds to the coil center?
Note that this option is only available when fit=yes otherwise there is no coil center to shift to.

**dividebyref = _'yes'_**
When this is enabled, the VSM's signal is divided by the reference signal from the motor.
This should help in correcting any errors introduced due to position dependent friction.

**sigfreqpos = _int_**\
**sigXpos = _int_**\
**sigYpos = _int_**\
**reffreqpos = _int_**\
**refXpos = _int_**\
**refYpos = _int_**\
Sometimes additional devices are measured simultaneously with the lock-ins, here you can set in which columns the relevant data is stored. 
Note that time will always be in the zeroth column.
