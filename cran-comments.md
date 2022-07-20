## Test environments
* Local Manjaro Linux 21.3.3 install, R 4.2.1
* Ubuntu Linux 20.04.1 LTS, R-release, GCC (via rhub)
* Fedora Linux, R-devel, clang, gfortran (via rhub)
* Windows Server 2022, R-devel, 64 bit (via rhub)


## R CMD check results

There were no ERRORs or WARNINGs.

There were 2 NOTEs:

* checking CRAN incoming feasibility ... NOTE

Uses the superseded packages: ‘doSNOW’, ‘snow’

  *Use of the 'doSNOW' package as opposed to the 'doParallel' package is required due to the support of the printed txtProgressBar in the 'doSNOW' package.*

* checking for detritus in the temp directory ... NOTE

Found the following files/directories:
    'lastMiKTeXException'

  * As noted in R-hub issue #503[https://github.com/r-hub/rhub/issues/503], this could be due to a bug/crash in MiKTeX and can likely be ignored.
