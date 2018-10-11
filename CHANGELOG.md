### 1.0
- Initial commit

### 1.1
- Made the code more readable
- Added options 'remove_excess'
- Added option 'dividebyref'
- Lots of spelling errors removed from readme
- Combined the if statement that checks frequency with the if statement that compares sigX with minpeakV.

### 1.3
- Instead of using the first two points to calculate the time between steps, it now calculates the average of the difference between all steps.

### 2.0
- Replaced sine fit function by a correct fit for a coil.
- A separate column is now required for the motor position, assuming the first peak is at beginstep has proved unfeasible.
- Due to the above, a extra data column is now required with voltages recorded from the potentiometer.
- Removed now unnecessary peaktime variable.
- Switched motor postion x-axis and millimeter x-axis to prevent fit function overlaying the data.
- The visible data is now on the millimeter axis, and the motor coordinate axis now has the transparent axis. Doing this makes matplotlib automatically choose the right colors.
- Cleanup of unused variables
- Replaced extrapolstep with extrapol, because it is user-friendly to set the amount to extrapolate in millimeter.

### 2.1
- Instead of discarding a line if it starts with '#' '\n' etc, it now uses .isdigit(). For increased compatibility with data files that have header lines that do not start with a hashtag.

### 2.2
- NaN, inf, -inf are now set to 0 when dividebyref is enabled, this to remove probelmatic points where the reference is 0
- Added option to set time in milliseconds for datafiles that use this
- Added 'yextrapol' to control how limits of y axis

### 2.3
- Removed the errorbars on the x-axis, the amplitude of the VSM is not actually an error.

### 2.4
- Remove plt.tight_layout, it causes problems with minus signs when exporting as pgf or tex, it works with pdf etc. But pgf and tex are preferred formats for processing with latex further on.
- Now using latex default font by default, this can of course always be changed by the user to their preferred font.

### 2.5
- Added option to control figure width/height.
- Added option to render text with tex.

### 2.6
- Now also plots the stddev of the position measurement on the x-axis
