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
