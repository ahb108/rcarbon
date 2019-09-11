[![cran version](http://www.r-pkg.org/badges/version/rcarbon)](https://CRAN.R-project.org/package=rcarbon) 
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/rcarbon?)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/rcarbon?color=82b4e8)](https://github.com/metacran/cranlogs.app)
# rcarbon <img src="/logo/rcarbon_logo.png" align="right" />

This package enables the calibration and analysis of radiocarbon dates, often but not exclusively for the purposes of archaeological research. It includes functions not only for basic calibration, uncalibration, and plotting of one or more dates, but also a statistical framework for building demographic and related longitudinal inferences from aggregate radiocarbon date lists, including: Monte-Carlo simulation test (Timpson et al 2014), random mark permutation test (Crema et al 2016, Bevan et al 2017) point-to-point test (Edinborough et al 2017), and spatial permutation test (Crema, Bevan, and Shennan 2017).

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
Bevan, A. Colledge,S., Fuller, D., Fyfe, R., Shennan, S., Stevens, C. 2017.[Holocene fluctuations in human population demonstrate repeated links to food production and climate](https://doi.org/10.1073/pnas.1709190114). Proceedings of the National Academy of Sciences of the United States of America. doi:10.1073/pnas.1709190114  

Crema, E.R., Bevan, A., Shennan, S., 2017. [Spatio-temporal approaches to archaeological radiocarbon dates](https://doi.org/10.1016/j.jas.2017.09.007). Journal of Archaeological Science 87, 1–9. doi:10.1016/j.jas.2017.09.007

Crema, E.R., Habu, J., Kobayashi, K., Madella, M., 2016. [Summed Probability Distribution of 14 C Dates Suggests Regional Divergences in the Population Dynamics of the Jomon Period in Eastern Japan](https://doi.org/10.1371/journal.pone.0154809). PLOS ONE 11, e0154809. doi:10.1371/journal.pone.0154809

Edinborough, K., Porcic, M., Martindale, A., Brown, T.J., Supernant, K, Ames, K.M. 2017. [Radiocarbon test for demographic events in written and oral history](https://doi.org/10.1073/pnas.1713012114), Proceedings of the National Academy of Sciences of the United States of America, 114 (47), 12436–12441. doi:10.1073/pnas.1713012114

Timpson, A., Colledge, S., Crema, E., Edinborough, K., Kerig, T., Manning, K., Thomas, M.G., Shennan, S., 2014. [Reconstructing regional population fluctuations in the European Neolithic using radiocarbon dates: a new case-study using an improved method](https://doi.org/10.1016/j.jas.2014.08.011). Journal of Archaeological Science 52, 549–557. doi:10.1016/j.jas.2014.08.011

### Case Studies Using _rcarbon_

Aubán, J.B., Jiménez-Puerto, J., Ruiz, P.E., Gordo, S.P., 2018. [C14 y poblamiento en las comarcas centro-meridionales del País Valenciano (c. 7000-1500 cal BC)](https://www.raco.cat/index.php/RecerquesMuseuAlcoi/article/view/348376). Recerques del Museu d’Alcoi 0, 35-48–48.

Bevan, A. Colledge,S., Fuller, D., Fyfe, R., Shennan, S., Stevens, C. 2017.[Holocene fluctuations in human population demonstrate repeated links to food production and climate](https://doi.org/10.1073/pnas.1709190114). Proceedings of the National Academy of Sciences of the United States of America, 114 (49) E10524-E10531. doi:10.1073/pnas.1709190114  

Broughton,J.M., Weitzel,E.M. 2018. [Population reconstructions for humans and megafauna suggest mixed causes for North American Pleistocene extinctions](https://doi.org/10.1038/s41467-018-07897-1).Nature Communications, 9, 5441. doi:10.1038/s41467-018-07897-1

Brown, A. A., Crema, E. R. 2019. [Māori Population Growth in Pre-contact New Zealand: Regional Population Dynamics Inferred From Summed Probability Distributions of Radiocarbon Dates](https://doi.org/10.1080/15564894.2019.1605429). The Journal of Island and Coastal Archaeology, 1–19. 

Cheddadi, R., Palmisano, A., López-Sáez, J. A., Nourelbait, M., Zielhofer, C., Tabel, J., et al. 2019. [Human demography changes in Morocco and environmental imprint during the Holocene](https://doi.org/10.1177/0959683619826657) The Holocene, 0959683619826657. doi:10.1177/0959683619826657

Crema, E.R., Bevan, A., Shennan, S., 2017. [Spatio-temporal approaches to archaeological radiocarbon dates](https://doi.org/10.1016/j.jas.2017.09.007). Journal of Archaeological Science 87, 1–9. doi:10.1016/j.jas.2017.09.007

de Souza, J. G., Robinson, M., Maezumi, S. Y., Capriles, J., Hoggarth, J. A., Lombardo, U., et al. 2019. [Climate change and cultural resilience in late pre-Columbian Amazonia](https://doi.org/10.1038/s41559-019-0924-0). Nature Ecology & Evolution, 1. doi:10.1038/s41559-019-0924-0

Feeser, I., Dörfler, W., Kneisel, J., Hinz, M., & Dreibrodt, S. 2019. [Human impact and population dynamics in the Neolithic and Bronze Age: Multi-proxy evidence from north-western Central Europe](https://doi.org/10.1177/0959683619857223). The Holocene, 0959683619857223. doi:10.1177/0959683619857223

Freeman, J., Baggio, J.A., Robinson, E., Byers, D.A., Gayo, E., Finley, J.B., Meyer, J.A., Kelly, R.L., Anderies, J.M., 2018. [Synchronization of energy consumption by human societies throughout the Holocene](https://doi.org/10.1073/pnas.1802859115). Proceedings of the National Academy of Sciences of the United States of America, 115 (40) 9962-9967. doi:10.1073/pnas.1802859115

Fyfe, R.M., Woodbridge, J., Palmisano, A., Bevan, A., Shennan, S., Burjachs, F., Legarra Herrero, B., García Puchol, O., Carrión, J.-S., Revelles, J., Roberts, C.N., 2019. [Prehistoric palaeodemographics and regional land cover change in eastern Iberia](https://doi.org/10.1177/0959683619826643). The Holocene 0959683619826643. doi:10.1177/0959683619826643

Jørgensen, E. K., & Riede, F. 2019. [Convergent catastrophes and the termination of the Arctic Norwegian Stone Age: A multi-proxy assessment of the demographic and adaptive responses of mid-Holocene collectors to biophysical forcing.](https://doi.org/10.1177/0959683619862036) The Holocene, 0959683619862036. 

Liu, L., Chen, X., Wright, H., Xu, H., Li, Y., Chen, G., et al. 2019. [Rise and fall of complex societies in the Yiluo region, North China: The spatial and temporal changes.](https://doi.org/10.1016/j.quaint.2019.05.025) Quaternary International. 

MacInnes, D. 2019. [The impact of population dynamics on social complexity in Neolithic Orkney.](https://doi.org/10.1016/j.jasrep.2019.02.036) Journal of Archaeological Science: Reports, 24, 721–728. doi:10.1016/j.jasrep.2019.02.036

Nielsen, S.V., Persson, P., Solheim, S., 2019. [De-Neolithisation in southern Norway inferred from statistical modelling of radiocarbon dates](https://doi.org/10.1016/j.jaa.2018.11.004). Journal of Anthropological Archaeology 53, 82–91. doi:10.1016/j.jaa.2018.11.004

Palmisano, A., Bevan, A., Shennan, S. 2017. [Comparing archaeological proxies for long-term population patterns: An example from central Italy](https://doi.org/10.1016/j.jas.2017.10.001). Journal of Archaeological Science, 87, 59-72. doi:10.1016/j.jas.2017.10.001

Palmisano, A., Woodbridge, J., Roberts, N., Bevan, A., Fyfe, R., Shennan, S., et al. 2019. [Holocene landscape dynamics and long-term population trends in the Levant.](https://doi.org/10.1177/0959683619826642) The Holocene, 0959683619826642. doi:10.1177/0959683619826642

Riris, P. 2018 [Dates as data revisited: A statistical examination of the Peruvian preceramic radiocarbon record](https://doi.org/10.1016/j.jas.2018.06.008). Journal of Archaeological Science, 97, 67-76. doi:10.1016/j.jas.2018.06.008

Riris, P., Arroyo-Kalin, M. 2018. [Widespread Population Decline in South America Correlates with Mid-holocene Climate Change](https://doi.org/10.1038/s41598-019-43086-w) Scientific Reports, 9, Article number: 6850 (2019). doi:10.1038/s41598-019-43086-w.

Roberts, N., Woodbridge, J., Palmisano, A., Bevan, A., Fyfe, R., & Shennan, S. 2019. [Mediterranean landscape change during the Holocene: Synthesis, comparison and regional trends in population, land cover and climate](https://doi.org/10.1177/0959683619826697). The Holocene, 0959683619826697. 

Robinson, E., Nicholson, C., & Kelly, R. L. 2019. [The Importance of Spatial Data to Open-Access National Archaeological Databases and the Development of Paleodemography Research](https://doi.org/10.1017/aap.2019.29). Advances in Archaeological Practice, 1–14. doi:10.1017/aap.2019.29


Saag, L., Laneman, M., Varul, L., Malve, M., Valk, H., Razzak, M. A., et al. 2019. [The Arrival of Siberian Ancestry Connecting the Eastern Baltic to Uralic Speakers further East](https://doi.org/10.1016/j.cub.2019.04.026). Current Biology,

Schauer, P., Shennan, S., Bevan, A., Cook, G., Edinborough, K., Fyfe, R., et al. 2019. [Supply and demand in prehistory? Economics of Neolithic mining in northwest Europe](https://doi.org/10.1016/j.jaa.2019.03.001). Journal of Anthropological Archaeology, 54, 149–160. doi:10.1016/j.jaa.2019.03.001

Stoddart, S., Woodbridge, J., Palmisano, A., Mercuri, A. M., Mensing, S. A., Colombaroli, D., et al. 2019. [Tyrrhenian central Italy: Holocene population and landscape ecology.](https://doi.org/10.1177/0959683619826696) The Holocene, 0959683619826696. doi:10.1177/0959683619826696

Szpak, P., Savelle, J. M., Conolly, J., & Richards, M. P. 2019. [Variation in late holocene marine environments in the Canadian Arctic Archipelago: Evidence from ringed seal bone collagen stable isotope compositions](https://doi.org/10.1016/j.quascirev.2019.03.016). Quaternary Science Reviews, 211, 136–155. doi:10.1016/j.quascirev.2019.03.016

Tallavara, M., Personen, P. 2018. [Human ecodynamics in the north-west coast of Finland 10,000–2000 years ago](https://doi.org/10.1016/j.quaint.2018.06.032).Quaternary International,doi:/10.1016/j.quaint.2018.06.032

van der Bilt, W. G. M., Born, A., & Haaga, K. A. 2019. [Was Common Era glacier expansion in the Arctic Atlantic region triggered by unforced atmospheric cooling?](https://doi.org/10.1016/j.quascirev.2019.07.042) Quaternary Science Reviews, 105860. doi:10.1016/j.quascirev.2019.07.042

Weiberg, E., Bevan, A., Kouli, K., Katsianis, M., Woodbridge, J., Bonnier, A., et al. 2019. [Long-term trends of land use and demography in Greece: A comparative study](https://doi.org/10.1177/0959683619826641) The Holocene, 0959683619826641. doi:10.1177/0959683619826641

Woodbridge, J., Roberts, C. N., Palmisano, A., Bevan, A., Shennan, S., Fyfe, R., et al. 2019. [Pollen-inferred regional vegetation patterns and demographic change in Southern Anatolia through the Holocene.](https://doi.org/10.1177/0959683619826635) The Holocene, 0959683619826635. doi:10.1177/0959683619826635




