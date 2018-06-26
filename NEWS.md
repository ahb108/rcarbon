
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
