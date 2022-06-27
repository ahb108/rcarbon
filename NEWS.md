# Version 1.4.4 (In Development)
* Bug Fix: Fixed a bug in `spd()` reducing the `timeRange` when `datenormalised` was set to TRUE.
* UPDATE: `calibrate()` now issues a warning in case dates have calibrated probabilities outside the user defined `timeRange` value.

# Version 1.4.3 (17 February 2022)
* UPDATE: Reversed legend item order in `plot.stackCalSPD()` to match display
* UPDATE: `hpdi()` now has the option to return a character vector with HPDI ranges in either BP or BCAD.

# Version 1.4.2 (15 March 2021)
* Bug fix: Fixed a number of minor bugs and errors caused when the summed probability vector contains zero. Functions affected by this bugs were:
  * proportion plots in `plot.stackCalSPD()`
  * `summary()` and `plot()` functions for `SpdModelTest` class objects
* Bug Fix: Fixed issues with `summary()` function on `CalDates` class object with multiple dates (thanks to [nferebau](https://github.com/nfrerebeau) for the bug report and fix)
`plot.stackCalSPD()` proportion plot to handle instances with 0 summed probability in the time Range of analysis.
* UPDATE: `binMed()` function can now handle larger `caldates` objects.
* NEW: `poolDates()` function allows combination of 14C ages associated with the same event using Ward and Wilson 1978 method.
* Minor bug fixes.
* Minor typos and errors in the help documentation.
* Updated dependencies for `spatstat(>= 2.0-0)`. 


# Version 1.4.1 (6 October 2020)
* UPDATE: Improved performance of the `calibrate()` function (c.a 400% faster)
* UPDATE: `multiplot` allows for rescaled calibrated probability distribution (new argument `rescale`) for improved readability.
* UPDATE: The package vignette now includes examples from new functions such as `multiplot()` and `stackspd()`.
* UPDATE: `hdmi()` now computes the probability mass of each age bracket (thanks to [hanecakr](https://github.com/hanecakr) for the suggestion and sample code)
* NEW: `transformSPD()` allows for taphonomic (and other user-defined) corrections to SPDs.
* Bug fix: Fixed color matching and improved label placing  in the `multiplot` function.
* Bug fix: calCurve='normal' no longer causing error in `calibrate()`

# Version 1.4.0 (15 August 2020)
* IntCal20, ShCal20, and Marine20 curves added, with IntCal20 as default calibration curve for `calibrate()`.

# Version 1.3.3 (4 August 2020)
* Minor changes on NAMESPACE and test environment to ensure workable CRAN checks for all operating systems.

#  Verson 1.3.2 (25 July 2020)
* Bug fix : `F14C=TRUE` was producting an error message in `calibrate()` after refactoring in version 1.3.1. 
* Bug fix : fixed error in the handling of `resErrors` (Delta R error) for  marine reservoir effect.
* NEW: `stackspd()` creates a set of multiple SPDs as an object of class `stackCalSPD` with a dedicated plot function.
* NEW: `multiplot()` function for displaying multiple calibrated dates.
* NEW: `hpdi()` function for computing highest probability density intervals of calibrated dates.
* NEW: `combine()` function for concatenating multiple `calDates` class objects.
* UPDATE: `calibrate()`, `modelTest()`, and `sptest()` are now parallelised using the *doSNOW* package, enabling progress bar when running over multiple cores.

#  Version 1.3.1 (18 March 2020)
* `binPrep()` now accepts alternative clustering algorithms
* Refactoring of calibration related functions to increase readability
* Minor bug fixes

#  Version 1.3.0 (12 December 2019)
* UPDATEs in `modelTest()`:
  * Added an error message and a note in the help documentation warning users to not use calibration curves not supported by `uncalibrate()`.
  * NEW argument `gridclip` to add range constrain in the simulation envelope. When set to TRUE the sampling is executed within the observed range of dates.
  * UPDATE `timeRange` can now be left undefined, in which case the range of median calibrated observed dates is used. 
  * NEW argument `normalised` handles whether simulated dates should be normalised or not.
  * NEW arguments `backsight` and `changexpr` enabling the comparison of expected and observed rates of changes. 
  * UPDATE: improved performance of `calsample` and `uncalsample` methods (update in `uncalibrate()`) `
  * Bug fix: users no longer need to specify `nsim` when `fitonly` is set to TRUE.
* UPDATE in `permTest()`
  * NEW arguments `backsight` and `changexpr` enabling the comparison of expected and observed rates of changes. 
* UPDATE in `spd()`: aggregation matrix now spans beyond `timeRange` to avoid edge effects.
* UPDATE in `calibrate()`: the new argument `F14C` enables calibration in F14C space.
* UPDATE in `SPpermTest()`: 
  * function is deprecated and renamed `sptest()`
  * function now allows user defined expression for computing the rate of change
* UPDATE in `spd2gg()`:
  * function is deprecated and renamed `spd2rc()`
  * function now allows user defined expression for computing the rate of change
  * rates of changes can now be based on abutting time-blocks (e.g. 5600-5501 to 5500-5401) or a fixed interval across all years.
* NEW: `thinDates()` function to randomly select a maximum number of dates per site, bin or phase.
* NEW: `as.CalDates()` function converts dates in `BchronCalibratedDates` class ([Bchron](https://cran.r-project.org/package=Bchron)) and `oxcAARCalibratedDatesList`class ([oxAAR](https://cran.r-project.org/package=oxcAAR)) into CalDates` class (rcarbon).
* NEW: `sampleDates()` function for sampling random dates from calibrated dates or bins.
* NEW: `ckde()` function for generating composite kernel density estimates.
* NEW: Suite of functions for mapping the spatial intensity of a set of radiocarbon dates via Kernel Density Estimates. 
  * `stkde()` Map the spatio-temporal intensity of a set of radiocarbon dates across multiple years
  * `spkde()` Map the spatial intensity of a set of radiocarbon dates for a given focal year.
  *  Associated `plot()` function.
* NEW: `subset()` and `which.CalDates()` for subsetting and extracting indices of calibrated dates based on temporal intervals described by logical conditions (e.g. between 5500 and 4500 cal BP) and user defined probability mass.
* Bug fix: fixed display error in some `plot()` functions when calendar is set to `BC/AD'
* Bug fix: fixed error in the computation of the combined uncertainties in `mixCurves()`. 
* Further minor bug fixes for R 4.0.0 

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
* NEW function: `p2pTest()` for comparing point to point differences in SPD.
* p-values of Monte-Carlo simulations are now all calculated using a different formula.
* minor bug fixes


# Version 1.0.0 (9 October 2017)
First official CRAN release. 
