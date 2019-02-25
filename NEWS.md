
# Development GitHub Version (ongoing)
* UPDATEs in `modelTest()`:
  * Added an error message and a note in the help documentation warning users to not use calibration curves not supported by `uncalibrate()`.
  * NEW argument `gridclip` to add range constrain in the simulation envelope. When set to TRUE the sampling is executed within the observed range of dates.
  * UPDATE `timeRange` can now be left undefined, in which case the range of median calibrated observed dates is used. 
  * NEW argument `normalised` handles whether simulated dates should be normalised or not.
  * UPDATE: improved performance of `calsample` and `uncalsample` methods (update in `uncalibrate()`) `
  * Bug fix: users no longer need to specify `nsim` when `fitonly` is set to TRUE.
* UPDATE in `spd()`: aggregation matrix now spans beyond `timeRange` to avoid edge effects.
* UPDATE in `calibrate()`: the new agrument `F14C` enables calibratio in F14C space.
* UPDATE in `SPpermTest()`: function is deprecated and renamed `sptest()`
* NEW: `thinDates()` function to randomly select a maximum number of dates per site, bin or phase.
* NEW: `as.CalDates()` function converts dates in `BchronCalibratedDates` class ([Bchron](https://cran.r-project.org/web/packages/Bchron/index.html)) and `oxcAARCalibratedDatesList`class ([oxAAR](https://cran.r-project.org/web/packages/oxcAAR/index.html)) into CalDates` class (rcarbon).
* NEW: `sampleDates()` function for sampling random dates from calibrated dates or bins.
* NEW: `ckde()` function for generating composite kernel density estimates.
* NEW: Suite of functions for mapping the spatial intensity of a set of radiocarbon dates via Kernel Density Estimates. 
 * `stkde()` Map the spatio-temporal intensity of a set of radiocarbon dates across multiple years
 * `spkde()` Map the spatial intensity of a set of radiocarbon dates for a given focal year.
 *  Associated `plot()` function.
* Bug fix: fixed dislay error in some `plot()` functions when calendar is set to `BC/AD'

# Version 1.2.0 (1 October 2018)
* Bug fix in `SPpermTest()`: Functions was not working when not running when `raw=FALSE` and `ncores=1`.
* Bug fix in `summary()` for the output of `modelTest()`: The number of bins reported was incorrect.
* UPDATE: `binsense()` requires a smaller number of non-optional arguments and allows for binning based on median calibrated dates.
* UPDATE: `plot()` function for geometric growth rates now allows for BC/AD calendar display.  
* UPDATE: `binPrep()` can now group dates based on median calibrated dates.
* UPDATE: `modelTest()`can now test SPD generated from multiple calibration curves.
* UPDATE: `modelTest()` now uses randomised  back-calibration for the 'calsample' method.
* UPDATE: `uncalibrate()`'s random back-calibration now utilises both the user supplied lab error and the error of the calibration curve.
* NEW: `mixCurves()` function generates mixed terrestrial/marine calibration curves.
* Further minor updates in plot labels, help documentation, and vignette.

# Version 1.1.3 (26 June 2018)
* Bug fix in `permTest()`: p-values calculation algorithm was yielding p-values of 0.
* NEW: `summary()` functions for the outputs of `permTest()` and `modelTest()`.

# Version 1.1.2 (22 May 2018)
* Bug Fix in `calibrate()`: calibration of multiple dates with different calibration curves was based on the calibration curve of the first sample.
* Bug Fix in  `modelTest()`: back-calibration routine was ignoring  user-supplied calibration curve and was solely using "intcal13". 
* UPDATE: `modelTest()` now allows for two distinct procedures for generating random dates from fitted models.
* Minor corrections (e.g. typos) in the help documentation.

# Version 1.1.1 (28 April 2018)
* Added a Vignette.
* Added normally-distributed (non-14C) age in `calibrate()`.
* Fixed a bug in `SPpermTest()` generating opposite results (positive deviations were recorded as negative deviations) when `ncores` was larger than 1.  
* Fixed a minor bug in the implementation of the [North el al 2002] formula
* Several changes in the`plot.SpatialTest()`. 
* NEW function: `spd2gg()` to convert SPD curves into geometric growth rates for given temporal blocks and an associated plot function. 
* Minor updates in the documentation.

# Version 1.1.0 (12 March 2018)
* Improved performance of `modelTest()` function when running with multiple cores.
* NEW function: `p2pTest()` for comparing point to point differences in SPD (as in [Edinborough et al 2017](http://dx.doi.org/10.1073/pnas.1713012114)
)
* p-values of Monte-Carlo simulations are now all calculated using the formula in [North el al 2002](http://dx.doi.org/10.1086/341527)
* minor bug fixes


# Version 1.0.0 (9 Otober 2017)
First official CRAN release. 
