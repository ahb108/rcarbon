## Test environments
* Local ubuntu 18.04 install, R 4.0.2
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (via rhub)
* Ubuntu Linux 16.04 LTS, R-release, GCC (via rhub)
* Fedora Linux, R-devel, clang, gfortran (via rhub)


## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE

Uses the superseded packages: ‘doSNOW’, ‘snow’

  *Use of the 'doSNOW' package as opposed to the 'doParallel' package is required due to the support of the printed txtProgressBar in the 'doSNOW' package.*

