# HFML-VSM pickupcoil profile plotter
A python program to make plots of the response of a VSM's pickup coils as a function of the position.

### Note:
I assume that you are already familiar with the concept of a VSM (Vibrating-sample magnetometer), if not then read the article in the Documents folder first.

## Long description
To measure the response of a VSM's pickup coils as a function of the position, set it up to go through a 'position loop'. This means that the VSM will start vibrating for a set amount of time (stepdur) at a set begin position, then stop. 
And then move a certain distance and start vibrating there again. 
The vibrating motion will induce a (~90 degree phase shifted) AC voltage across the pickup coils, which will depend strongly on the position at which the VSM is currently vibrating. 
When this voltage is measured (e.g using 2 Lock-In amplifiers; 1 for reference and 1 for the VSM's signal) it can be plotted against the VSM's position. 
This program identifies each peak in the signal, calculates it's average value and standard deviation, and then plots is against the position. 
Optionally it can also fit a theoretical curve through the data, and identify the center of the coil (if it was within the range of the measured data). 
If a center has been found, it also has an option to shift the x-axis such that 0 is in the center of the coil.

## Dependencies
The versions indicated here are the versions that I am currently using, but there is no reason other versions shouldn't work to.
- python-3.6.6
- matplotlib-2.2.2
- numpy-1.14.5
- scipy-1.1.0
- Optionally: matplotlib2tikz-0.6.18 (when saving as .tex)

## Other requirements
The program relies on the use of at least 1 lock-in to measure the VSM's signal, for best results use a second lock-in to record a reference signal from a potentiometer connected to the motor.
To determine the position of the data points correctly the program requires that the position be recorded simultaneously with the signal, this can be accomplished by using a 1Hz low pass filter on the reference signal and recording this voltage. 
This voltage should then be calibrated with the motor's position, the program has 2 config options for these calibration parameters (y=a*x+b).
Optionally it is possible to make this reference dimensionless by dividing it by the voltage across applied to the potentiometer.
Note that in this case the parameters 'a' and 'b' should be calculated using the dimensionless reference as well.

## How does it work?
Open the coilprofile.py file and edit the config options. See also the overview of options below.
Then run it like this:\
**Linux:**
```
./coilprofile.py inputfile
```
**Windows:**
```
C:\ProgramData\Anaconda3\python.exe D:\Coil-Profile-Plotter\HFML-VSM_CoilProfilePlotter\coilprofile.py inputfile
```
And it will process the data in `inputfile` and show you the plot using plt.show()
 
Or optionally use a second argument like this:\
**Linux:**
```
./coilprofile.py inputfile outputfile
```
**Windows:**
```
C:\ProgramData\Anaconda3\python.exe D:\Coil-Profile-Plotter\HFML-VSM_CoilProfilePlotter\coilprofile.py inputfile outputfile
```
And the program will use plt.savefig() to save the plot to `outputfile` instead.
Note that `outputfile` should have one of the filename extensions supported by matplotlib (e.g. `.pdf`, `.png`)
or the `.tex` extension, in which case it will use matplotlib2tikz.

## Can I use/adapt/copy your program?
Of course you can, it's licensed under GPLv3, just give credit when due and keep it open source.

## I have an improvement/comment
Awesome, please leave it in the bug reports and I will look at it ASAP.

## My plot looks like crap, what can I do?
Try increasing 'reacttime' and maybe play a bit with 'minpeakV' and 'deltafreq', see the overview of config options below.

## The scipy fit fails, what happened?
For a function as complex as this one, scipy's curve_fit requires very good starting estimates. 
Make sure that the 'guess_r1', 'guess_r2', 'guess_L' and 'guess_N' variables are set to very close estimates. 
Also provide a good guess for guess_x0, this can be done by first plotting with fit='no' and manually determining where approximately the first peak occurs.

## Help, it doesn't work
If you run into anything, and you can't work it out with the comments in the program or one of the options. You can always leave a bug report here and I'll see if I can help.

## Overview of options in the config (the config is just the top of the program)

**convert = _float_** 
How many millimeters does the VSM motor move when it moves 1 in motor coordinates.

**mmextrapol = _float_**
If fit=yes, how far should the fit be extrapolated (in mm)
I recommend you start with this set to 0 and then see how that looks first. And if it looks safe turn this up a bit. Note that this variable also controls the limits of the x-axis.

**yextrapol = _float_**
This parameter controls the limits of the y axis

**stepdur = _float_**
How long does the VSM vibrate at one position (in seconds)? 
This is used in combination with reacttime to calculate how many points there are within a peak.

**reacttime = _float_**
The voltage will not jump instantly from the background noise level to a peak, the data points that are above the threshold but not in the peak should be ignored. Set this according to how long it takes for the system to respond, too high is better then too low (in seconds).

**freq = _float_**
The frequency of the VSM's vibration, this is compared to the frequency measured by the lock-in (in Hz).

**deltafreq = _float_**
Any data points that are not between 'freq-deltafreq' and 'freq+deltafreq' are rejected. To accept all points just set this equal to 'freq' or very large.

**minpeakV = _float_**
How high should the signal be before it is considered a peak?

**a = _float_**\
**b = _float_**\
Motor_coordinates = a * potentiometer_DC_voltage + b

**rm_lastpoint = _int_**
Sometimes the last points are off, this integer controls the number of peaks to remove from the end.

**rm_firstpoint = _int_**
Same as 'rm_lastpoint' but for the first peaks

**fit = _'yes'_**
Fit the theoretical coil profile to the data, this works best if both ends of the coil were within the range of the position loop. If the fit is weird or fails, try adjusting the 'guess_' parameters.

**guess_x0 = _float_** Offset of the function, if unsure plot with fit='no' first, and determine the location of the first peak from the plot (in mm).\
**guess_amp = _float_** What was the VSM's amplitude (in mm).\
**guess_r1 = _float_** Approximate inner radius of the coil (in mm).\
**guess_r2 = _float_** Approximate outer radius of the coil (in mm).\
**guess_L = _float_** Approximate length of the coil (in mm).\
**guess_N = _float_** Approximate number of windings of the coil.

**makecoordrel = _'yes'_**
Should the upper x-axis coordinates be shifted such that 0 corresponds to the coil center?
Note that this option is only available when fit=yes otherwise there is no coil center to shift to.

**usetex = _'yes'_**
Use LaTeX to render all text, highly recommended when including plots in latex files.
Note that saving as .pgf or .tex gives best results in latex. 
This option has no effect when saving as .tex, as text won't be rendered until the document is compiled.

**plotheight = _float_**\
**plotwidth = _float_**\
Set the height and width of the plot (in mm).
Note that this is not used when saving plot as .tex, in this case the height and width need to be specified in the file that the plot is included into:
```
\newlength\figH
\newlength\figW
\setlength{\figH}{10cm}
\setlength{\figW}{\textwidth}
\input{filename.tex}
```
This has the advantage of being able to change the plots size without re-running this program, if however you have problems with .tex plots, try saving as .pgf instead. 

**dividebyref = _'yes'_**
When this is enabled, the VSM's signal is divided by the reference signal from the motor.
This should help in correcting any errors introduced due to position dependent friction.
Note that disabling this might cause fit to fail, it is recommended to keep this enabled at all times unless you have no reference data available.

**absval = _'yes'_**
If enabled the absolute value of the real and complex part is used instead of just the real part.
Note that the use of the lock-in's 'autophase' function is still required to correctly determine the sign.

**relref = _'yes'_**
This makes the reference dimensionless by dividing it by the voltage across the whole potentiometer.
This is useful when comparing measurements for which the potentiometer's voltage was different.
Note that if this is used the parameters 'a' and 'b' should be calculated using the dimensionless reference as well.

**timepos = _int_**\
**sigfreqpos = _int_**\
**sigXpos = _int_**\
**sigYpos = _int_**\
**reffreqpos = _int_**\
**refXpos = _int_**\
**refYpos = _int_**\
**potmetDCpos = _int_**\
**batDCpos = _int_**\
Here you can set in which columns of the data file the relevant data is stored.

**sig_corr = _float_**\
**ref_corr = _float_**\
The signal and reference are divided by these numbers to correct for any amplifications during measurement.

**timeinmillisec = _'yes'_** Set to yes if the time in the time-column is in millliseconds, if disabled time is assumed to be in seconds.
