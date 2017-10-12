[![cran version](http://www.r-pkg.org/badges/version/rcarbon)](https://cran.rstudio.com/web/packages/rcarbon) 
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/rcarbon?)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/rcarbon?color=82b4e8)](https://github.com/metacran/cranlogs.app)
# rcarbon

This package enables the calibration and analysis of radiocarbon dates, often but not exclusively for the purposes of archaeological research. It includes functions not only for basic calibration, uncalibration, and plotting of one or more dates, but also a statistical framework for building demographic and related longitudinal inferences from aggregate radiocarbon date lists, including: Monte-Carlo simulation test (Timpson et al 2014), random mark permutation test (Crema et al 2016) and spatial permutation tests (Crema, Bevan, and Shennan 2017).

To install the stable version from CRAN:

```
install.packages('rcarbon')
```

To install the latest development version:

```
# "devtools" package required 
devtools::install_github('ahb108/rcarbon')
```

### References
Crema, E.R., Bevan, A., Shennan, S., 2017. [Spatio-temporal approaches to archaeological radiocarbon dates](https://doi.org/10.1016/j.jas.2017.09.007). Journal of Archaeological Science 87, 1–9. doi:10.1016/j.jas.2017.09.007

Crema, E.R., Habu, J., Kobayashi, K., Madella, M., 2016. [Summed Probability Distribution of 14 C Dates Suggests Regional Divergences in the Population Dynamics of the Jomon Period in Eastern Japan](https://doi.org/10.1371/journal.pone.0154809). PLOS ONE 11, e0154809. doi:10.1371/journal.pone.0154809

Timpson, A., Colledge, S., Crema, E., Edinborough, K., Kerig, T., Manning, K., Thomas, M.G., Shennan, S., 2014. [Reconstructing regional population fluctuations in the European Neolithic using radiocarbon dates: a new case-study using an improved method](https://doi.org/10.1016/j.jas.2014.08.011). Journal of Archaeological Science 52, 549–557. doi:10.1016/j.jas.2014.08.011


