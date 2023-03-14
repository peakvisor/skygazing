# skygazing
Lightweight header-only C++ library for astronomical calculations. Developed for use in [PeakVisor](https://peakvisor.com). 

![PeakVisor skygazing](https://raw.githubusercontent.com/peakvisor/skygazing/main/peakvisor-skygazing.jpg)

Based on [mourner/suncalc](https://github.com/mourner/suncalc), [AAPlus library](http://www.naughter.com/aa.html), [Fabiz/MeeusJs](https://github.com/Fabiz/MeeusJs) and The Book:

// [AA] “Astronomical Algorithms” 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.

To test against AAPlus library, use CMake argument -DSKYGAZING_AAPLUS_SRC="/path/to/aaplus_source/"
and build ${SKYGAZING_AAPLUS_SRC}/cmake-build-release/lib/libaaplus.dylib yourself. 

## Usage examples

For examples of usage, see printObservation method in [skygazing_test.h](test/skygazing_test.h)
