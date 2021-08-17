[![cran version](http://www.r-pkg.org/badges/version/rcarbon)](https://CRAN.R-project.org/package=rcarbon) 
[![development version](https://img.shields.io/badge/devel%20version-1.4.3-blue.svg)](https://github.com/ahb108/rcarbon)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/rcarbon?)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/rcarbon?color=82b4e8)](https://github.com/metacran/cranlogs.app)
[![R-CMD-check](https://github.com/ahb108/rcarbon/workflows/R-CMD-check/badge.svg)](https://github.com/ahb108/rcarbon/actions)
# rcarbon <img src="/logo/rcarbon_logo.png" align="right" />

This package enables the calibration and analysis of radiocarbon dates, often but not exclusively for the purposes of archaeological research. It includes functions not only for basic calibration, uncalibration, and plotting of one or more dates, but also a statistical framework for building demographic and related longitudinal inferences from aggregate radiocarbon date lists, including: Monte-Carlo simulation test (Timpson et al 2014), random mark permutation test (Crema et al 2016, Bevan et al 2017) point-to-point test (Edinborough et al 2017), and spatial permutation test (Crema, Bevan, and Shennan 2017). For more details see the associated paper (Crema and Bevan 2021).

To install the stable version from CRAN:

```
install.packages('rcarbon')
```

To install the latest development version:

```
# "devtools" package required 
devtools::install_github('ahb108/rcarbon')
```
### Citation
To cite _rcarbon_ in pubblications use:

Crema, E.R., Bevan, A. 2021 [Inference from Large Sets of Radiocarbon Dates: Software and Methods](https://doi.org/10.1017/RDC.2020.95) Radiocarbon 63(1),23-39  doi:10.1017/RDC.2020.95

### References
Bevan, A. Colledge,S., Fuller, D., Fyfe, R., Shennan, S., Stevens, C. 2017.[Holocene fluctuations in human population demonstrate repeated links to food production and climate](https://doi.org/10.1073/pnas.1709190114). Proceedings of the National Academy of Sciences of the United States of America. doi:10.1073/pnas.1709190114  

Crema, E.R., Bevan, A. 2021 [Inference from Large Sets of Radiocarbon Dates: Software and Methods](https://doi.org/10.1017/RDC.2020.95) Radiocarbon 63(1),23-39. doi:10.1017/RDC.2020.95

Crema, E.R., Bevan, A., Shennan, S., 2017. [Spatio-temporal approaches to archaeological radiocarbon dates](https://doi.org/10.1016/j.jas.2017.09.007). Journal of Archaeological Science 87, 1–9. doi:10.1016/j.jas.2017.09.007

Crema, E.R., Habu, J., Kobayashi, K., Madella, M., 2016. [Summed Probability Distribution of 14 C Dates Suggests Regional Divergences in the Population Dynamics of the Jomon Period in Eastern Japan](https://doi.org/10.1371/journal.pone.0154809). PLOS ONE 11, e0154809. doi:10.1371/journal.pone.0154809

Edinborough, K., Porcic, M., Martindale, A., Brown, T.J., Supernant, K, Ames, K.M. 2017. [Radiocarbon test for demographic events in written and oral history](https://doi.org/10.1073/pnas.1713012114), Proceedings of the National Academy of Sciences of the United States of America, 114 (47), 12436–12441. doi:10.1073/pnas.1713012114

Timpson, A., Colledge, S., Crema, E., Edinborough, K., Kerig, T., Manning, K., Thomas, M.G., Shennan, S., 2014. [Reconstructing regional population fluctuations in the European Neolithic using radiocarbon dates: a new case-study using an improved method](https://doi.org/10.1016/j.jas.2014.08.011). Journal of Archaeological Science 52, 549–557. doi:10.1016/j.jas.2014.08.011

### Case Studies Using _rcarbon_

Arroyo-Kalin, M., & Riris, P. (2021). [Did pre-Columbian populations of the Amazonian biome reach carrying capacity during the Late Holocene?](https://doi.org/10.1098/rstb.2019.0715) Philosophical Transactions of the Royal Society B: Biological Sciences, 376(1816), 20190715. 

Aubán, J.B., Jiménez-Puerto, J., Ruiz, P.E., Gordo, S.P., 2018. [C14 y poblamiento en las comarcas centro-meridionales del País Valenciano (c. 7000-1500 cal BC)](https://www.raco.cat/index.php/RecerquesMuseuAlcoi/article/view/348376). Recerques del Museu d’Alcoi 0, 35-48–48.

Barnett, R. L., Charman, D. J., Johns, C., Ward, S. L., Bevan, A., Bradley, S. L., et al. 2020. [Nonlinear landscape and cultural response to sea-level rise.](https://doi.org/10.1126/sciadv.abb637) Science Advances, 6(45), eabb6376. 6

Becker, F., Knitter, D., Nykamp, M., Schütt, B., 2020. [Meta-Analysis of Geomorphodynamics in the Western Lower Bakırçay Plain (Aegean Region, Turkey)](https://doi.org/10.3390/land9090338). Land 9, 338. 

Bergsvik, K. A., Darmark, K., Hjelle, K. L., Aksdal, J., & Åstveit, L. I. (2021). [Demographic developments in Stone Age coastal western Norway by proxy of radiocarbon dates, stray finds and palynological data](https://doi.org/10.1016/j.quascirev.2021.106898). Quaternary Science Reviews, 259, 106898. 

Bevan, A. Colledge,S., Fuller, D., Fyfe, R., Shennan, S., Stevens, C. 2017.[Holocene fluctuations in human population demonstrate repeated links to food production and climate](https://doi.org/10.1073/pnas.1709190114). Proceedings of the National Academy of Sciences of the United States of America, 114 (49) E10524-E10531. doi:10.1073/pnas.1709190114  

Bevan, A., & Crema, E. R. 2021. [Modifiable reporting unit problems and time series of long-term human activity](https://doi.org/10.1098/rstb.2019.0726). Philosophical Transactions of the Royal Society B: Biological Sciences, 376(1816), 20190726. 

Bird, D., Freeman, J., Robinson, E., Maughan, G., Finley, J.B., Lambert, P.M., Kelly, R.L., 2020. [A first empirical analysis of population stability in North America using radiocarbon records](https://doi.org/10.1177/0959683620919975). The Holocene. 

Broodbank, C., Lucarini, G. 2019. [The Dynamics of Mediterranean Africa, ca. 9600–1000 bc: An Interpretative Synthesis of Knowns and Unknowns](https://doi.org/10.1558/jma.40581), Journal of Mediterranean Archaeology, 32(2), 195-267.

Broughton,J.M., Weitzel,E.M. 2018. [Population reconstructions for humans and megafauna suggest mixed causes for North American Pleistocene extinctions](https://doi.org/10.1038/s41467-018-07897-1).Nature Communications, 9, 5441. doi:10.1038/s41467-018-07897-1

Brown, A. A., Crema, E. R. 2019. [Māori Population Growth in Pre-contact New Zealand: Regional Population Dynamics Inferred From Summed Probability Distributions of Radiocarbon Dates](https://doi.org/10.1080/15564894.2019.1605429). The Journal of Island and Coastal Archaeology, 1–19. 

Carter, V. A., Brunelle, A., Power, M. J., DeRose, R. J., Bekker, M. F., Hart, I., et al. (2021). [Legacies of Indigenous land use shaped past wildfire regimes in the Basin-Plateau Region, USA](https://doi.org/10.1038/s43247-021-00137-3). Communications Earth & Environment, 2(1), 1–9. 

Cascalheira, J., Alcaraz-Castaño, M., Alcolea-González, J., de Andrés-Herrero, M., Arrizabalaga, A., Aura Tortosa, J.E., Garcia-Ibaibarriaga, N., Iriarte-Chiapusso, M.-J., 2020. [Paleoenvironments and human adaptations during the Last Glacial Maximum in the Iberian Peninsula: A review.](https://doi.org/10.1016/j.quaint.2020.08.005) Quaternary International. 

Cheddadi, R., Palmisano, A., López-Sáez, J. A., Nourelbait, M., Zielhofer, C., Tabel, J., et al. 2019. [Human demography changes in Morocco and environmental imprint during the Holocene](https://doi.org/10.1177/0959683619826657) The Holocene, 0959683619826657. doi:10.1177/0959683619826657

Crema, E.R., Bevan, A. 2020 [Inference from Large Sets of Radiocarbon Dates: Software and Methods](https://doi.org/10.1017/RDC.2020.95) Radiocarbon, doi:10.1017/RDC.2020.95

Crema, E.R., Bevan, A., Shennan, S., 2017. [Spatio-temporal approaches to archaeological radiocarbon dates](https://doi.org/10.1016/j.jas.2017.09.007). Journal of Archaeological Science 87, 1–9. doi:10.1016/j.jas.2017.09.007

Crema, E. R., & Shoda, S. 2021. [A Bayesian approach for fitting and comparing demographic growth models of radiocarbon dates: A case study on the Jomon-Yayoi transition in Kyushu (Japan)](https://doi.org/10.1371/journal.pone.0251695). PLOS ONE, 16(5), e0251695. 

Crema, E.R. 2020 [Non-Stationarity and Local Spatial Analysis](https://doi.org/10.4324/9781351243858), In Gillings, M., Hacıgüzeller, P., Lock,G. (Eds.) Archaeological Spatial Analysis: A Methodological Guide, Taylor & Francis, London.

Crema, E.R., Kobayashi, K., 2020. [A multi-proxy inference of Jōmon population dynamics using bayesian phase models, residential data, and summed probability distribution of 14C dates](https://doi.org/10.1016/j.jas.2020.105136). Journal of Archaeological Science 117, 105136. 

Deforce K., Groenewoudt B., Haneca K. 2020. [2500 years of charcoal production in the Low Countries: The chronology and typology of charcoal kilns and their relation with early iron production.](https://doi.org/10.1016/j.quaint.2020.10.020). Quaternary International.

de Groot, B., 2020. [The impact of population fluctuations on the spatial spread of Neolithic ceramic traditions in West Anatolia and South-East Europe.](https://doi.org/10.1016/j.jaa.2019.101121) Journal of Anthropological Archaeology 57, 101121. doi:10.1016/j.jaa.2019.101121

de Souza, J. G., Robinson, M., Maezumi, S. Y., Capriles, J., Hoggarth, J. A., Lombardo, U., et al. 2019. [Climate change and cultural resilience in late pre-Columbian Amazonia](https://doi.org/10.1038/s41559-019-0924-0). Nature Ecology & Evolution, 1. doi:10.1038/s41559-019-0924-0

de Souza, P., Cartajena, I., Riquelme, R., Maldonado, A., Porras, M. E. de, Santander, B., et al. (2021). [Late Pleistocene–Early Holocene human settlement and environmental dynamics in the southern Atacama Desert highlands (24.0°S–24.5°S, Northern Chile)](https://doi.org/10.1002/gea.21849). Geoarchaeology.

DiNapoli, R. J., Crema, E. R., Lipo, C. P., Rieth, T. M., & Hunt, T. L. (2021). [Approximate Bayesian Computation of radiocarbon and paleoenvironmental record shows population resilience on Rapa Nui (Easter Island)](https://doi.org/10.1038/s41467-021-24252-z). Nature Communications, 12(1), 3939. 

Doering, B. N. 2021. [Subarctic landscape adaptations and paleodemography: A 14,000-year history of climate change and human settlement in central Alaska and Yukon](https://doi.org/10.1016/j.quascirev.2021.107139). Quaternary Science Reviews, 268, 107139. 

Edinborough, K., Martineau, R., Dufraisse, A., Shennan, S., Imbeaux, M., Dumontet, A., et al. (2021). [A Neolithic population model based on new radiocarbon dates from mining, funerary and population scaled activity in the Saint-Gond Marshes region of North East France](https://doi.org/10.1016/j.quaint.2021.03.001). Quaternary International. 

Ekholm, T. 2021. [Hunter-gatherer adaptions during the Early Holocene in Northern Sweden](https://doi.org/10.1177/0959683620961482). The Holocene, 31(1), 83–94. 

Feeser, I., Dörfler, W., Kneisel, J., Hinz, M., & Dreibrodt, S. 2019. [Human impact and population dynamics in the Neolithic and Bronze Age: Multi-proxy evidence from north-western Central Europe](https://doi.org/10.1177/0959683619857223). The Holocene, 0959683619857223. doi:10.1177/0959683619857223

Freeman, J., Baggio, J.A., Robinson, E., Byers, D.A., Gayo, E., Finley, J.B., Meyer, J.A., Kelly, R.L., Anderies, J.M., 2018. [Synchronization of energy consumption by human societies throughout the Holocene](https://doi.org/10.1073/pnas.1802859115). Proceedings of the National Academy of Sciences of the United States of America, 115 (40) 9962-9967. doi:10.1073/pnas.1802859115

Freeman, J., Hard, R. J., Mauldin, R. P., & Anderies, J. M. (2021). [Radiocarbon data may support a Malthus-Boserup model of hunter-gatherer population expansion](https://doi.org/10.1016/j.jaa.2021.101321). Journal of Anthropological Archaeology, 63, 101321. 

Fritz, C., Tosello, G., Fleury, G., Kasarhérou, E., Walter, P., Duranthon, F., et al. 2021. [First record of the sound produced by the oldest Upper Paleolithic seashell horn](https://doi.org/10.1126/sciadv.abe9510). Science Advances, 7(7), eabe9510. 

Fyfe, R.M., Woodbridge, J., Palmisano, A., Bevan, A., Shennan, S., Burjachs, F., Legarra Herrero, B., García Puchol, O., Carrión, J.-S., Revelles, J., Roberts, C.N., 2019. [Prehistoric palaeodemographics and regional land cover change in eastern Iberia](https://doi.org/10.1177/0959683619826643). The Holocene 0959683619826643. doi:10.1177/0959683619826643

Gil, A.F., Villalba, R., Franchetti, F.R., Otaola, C., Abbona, C.C., Peralta, E.A., Neme, G., 2020. [Between Foragers and Farmers: Climate Change and Human Strategies in Northwestern Patagonia](https://doi.org/10.3390/quat3020017). Quaternary 3, 17. 

Gjesfjeld, E., Etnier, M.A., Takase, K., Brown, W.A., Fitzhugh, B., 2020. [Biogeography and adaptation in the Kuril Islands, Northeast Asia.](https://doi.org/10.1080/00438243.2019.1715248) World Archaeology, 1–25. 

Gjesfjeld, E., Silvestro, D., Chang, J., Koch, B., Foster, J.G., Alfaro, M.E., 2020. [A quantitative workflow for modeling diversification in material culture.](https://doi.org/10.1371/journal.pone.0227579) PLOS ONE 15, e0227579. 

Gretzinger, J., Molak, M., Reiter, E., Pfrengle, S., Urban, C., Neukamm, J., et al. 2019. [Large-scale mitogenomic analysis of the phylogeography of the Late Pleistocene cave bear]( https://doi.org/10.1038/s41598-019-47073-z). Scientific Reports, 9. DOI: 10.1038/s41598-019-47073-z.

Harrison, S.P., Gaillard, M.-J., Stocker, B.D., Vander Linden, M., Klein Goldewijk, K., Boles, O., Braconnot, P., Dawson, A., Fluet-Chouinard, E., Kaplan, J.O., Kastner, T., Pausata, F.S.R., Robinson, E., Whitehouse, N.J., Madella, M., Morrison, K.D., 2020. [Development and testing scenarios for implementing land use and land cover changes during the Holocene in Earth system model experiments](https://doi.org/10.5194/gmd-13-805-2020). Geoscientific Model Development 13, 805–824. 

Jones, T. L., Coltrain, J. B., Jacobs, D. K., Porcasi, J., Brewer, S. C., Buckner, J. C., et al. (2021). [Causes and consequences of the late Holocene extinction of the marine flightless duck (Chendytes lawi) in the northeastern Pacific](https://doi.org/10.1016/j.quascirev.2021.106914). Quaternary Science Reviews, 260, 106914. 

Jørgensen, E. K., & Riede, F. 2019. [Convergent catastrophes and the termination of the Arctic Norwegian Stone Age: A multi-proxy assessment of the demographic and adaptive responses of mid-Holocene collectors to biophysical forcing.](https://doi.org/10.1177/0959683619862036) The Holocene, 0959683619862036. 

Jorgeson, I.A., Breslawski, R.P., Fisher, A.E., 2020. [Radiocarbon simulation fails to support the temporal synchroneity requirement of the Younger Dryas impact hypothesis](https://doi.org/10.1017/qua.2019.83). Quaternary Research 1–17. 

Kwak, S., Obata, H., Lee, G.-A., 2020. [Broad-spectrum foodways in southern coastal Korea in the Holocene: Isotopic and archaeobotanical signatures in Neolithic shell middens.](https://doi.org/10.1080/15564894.2020.1776427) The Journal of Island and Coastal Archaeology 0, 1–29. 

Lawrence, D., Palmisano, A., & Gruchy, M. W. de. 2021. [Collapse and continuity: A multi-proxy reconstruction of settlement organization and population trajectories in the Northern Fertile Crescent during the 4.2kya Rapid Climate Change event](https://doi.org/10.1371/journal.pone.0244871). PLOS ONE, 16(1), e0244871. 

Lima, M., Gayo, E.M., Latorre, C., Santoro, C.M., Estay, S.A., Cañellas-Boltà, N., Margalef, O., Giralt, S., Sáez, A., Pla-Rabes, S., Chr. Stenseth, N., 2020. [Ecology of the collapse of Rapa Nui society](https://doi.org/10.1098/rspb.2020.0662). Proceedings of the Royal Society B: Biological Sciences 287, 20200662. 

Liu, L., Chen, X., Wright, H., Xu, H., Li, Y., Chen, G., et al. 2019. [Rise and fall of complex societies in the Yiluo region, North China: The spatial and temporal changes.](https://doi.org/10.1016/j.quaint.2019.05.025) Quaternary International. 

López, J. M., Neme, G., & Gil, A. F. 2019. [Resource intensification and zooarchaeological record in the southern margins of pre-Hispanic Andean agriculture.](https://doi.org/10.1007/s12520-019-00857-w) Archaeological and Anthropological Sciences, 11(10), 5287–5300. 

Lucarini, G., Wilkinson, T., Crema, E.R., Palombini, A., Bevan, A., Broodbank, C., 2020. [The MedAfriCarbon Radiocarbon Database and Web Application. Archaeological Dynamics in Mediterranean Africa, ca. 9600–700 BC.](https://doi.org/10.5334/joad.60) Journal of Open Archaeology Data 8, 1. 

MacInnes, D. 2019. [The impact of population dynamics on social complexity in Neolithic Orkney.](https://doi.org/10.1016/j.jasrep.2019.02.036) Journal of Archaeological Science: Reports, 24, 721–728. doi:10.1016/j.jasrep.2019.02.036

Nielsen, S.V., Persson, P., Solheim, S., 2019. [De-Neolithisation in southern Norway inferred from statistical modelling of radiocarbon dates](https://doi.org/10.1016/j.jaa.2018.11.004). Journal of Anthropological Archaeology 53, 82–91. doi:10.1016/j.jaa.2018.11.004

Nielsen, S. V. (2021). [A Late Mesolithic Forager Dispersal Caused Pre-Agricultural Demographic Transition in Norway](https://doi.org/10.1111/ojoa.12218). Oxford Journal of Archaeology, 40(2), 153–175. 

Palmisano, A., Bevan, A., Shennan, S. 2017. [Comparing archaeological proxies for long-term population patterns: An example from central Italy](https://doi.org/10.1016/j.jas.2017.10.001). Journal of Archaeological Science, 87, 59-72. doi:10.1016/j.jas.2017.10.001

Palmisano, A., Lawrence, D., de Gruchy, M. W., Bevan, A., & Shennan, S. 2021. [Holocene regional population dynamics and climatic trends in the Near East: A first comparison using archaeo-demographic proxies](https://doi.org/10.1016/j.quascirev.2020.106739). Quaternary Science Reviews, 252, 106739. 

Palmisano, A., Woodbridge, J., Roberts, N., Bevan, A., Fyfe, R., Shennan, S., et al. 2019. [Holocene landscape dynamics and long-term population trends in the Levant.](https://doi.org/10.1177/0959683619826642) The Holocene, 0959683619826642. doi:10.1177/0959683619826642

Porčić, M., 2020. [Observations on the origin and demography of the Vinča culture](https://doi.org/10.1016/j.quaint.2020.04.012). Quaternary International. 

Porčić, M., Blagojević, T., Pendić, J., & Stefanović, S. 2021. [The Neolithic Demographic Transition in the Central Balkans: population dynamics reconstruction based on new radiocarbon evidence](https://doi.org/10.1098/rstb.2019.0712). Philosophical Transactions of the Royal Society B: Biological Sciences, 376(1816), 20190712. 

Prates, L., & Perez, S. I. (2021). [Late Pleistocene South American megafaunal extinctions associated with rise of Fishtail points and human population](https://doi.org/10.1038/s41467-021-22506-4). Nature Communications, 12(1), 2175. 

Prates, L., Politis, G.G., Perez, S.I., 2020. [Rapid radiation of humans in South America after the last glacial maximum: A radiocarbon-based study.](https://doi.org/10.1371/journal.pone.0236023) PLOS ONE 15, e0236023. 

Ren, X., Xu, J., Wang, H., Storozum, M., Lu, P., Mo, D., et al. 2021. [Holocene fluctuations in vegetation and human population demonstrate social resilience in the prehistory of the Central Plains of China](https://doi.org/10.1088/1748-9326/abdf0a). Environmental Research Letters. 

Riris, P. 2018 [Dates as data revisited: A statistical examination of the Peruvian preceramic radiocarbon record](https://doi.org/10.1016/j.jas.2018.06.008). Journal of Archaeological Science, 97, 67-76. doi:10.1016/j.jas.2018.06.008

Riris, P., Arroyo-Kalin, M. 2018. [Widespread Population Decline in South America Correlates with Mid-holocene Climate Change](https://doi.org/10.1038/s41598-019-43086-w) Scientific Reports, 9, Article number: 6850 (2019). doi:10.1038/s41598-019-43086-w.

Riris, P., & Silva, F. 2021. [Resolution and the detection of cultural dispersals: development and application of spatiotemporal methods in Lowland South America.](https://doi.org/10.1057/s41599-021-00717-w) Humanities and Social Sciences Communications, 8(1), 1–13. 

Roberts, N., Woodbridge, J., Palmisano, A., Bevan, A., Fyfe, R., & Shennan, S. 2019. [Mediterranean landscape change during the Holocene: Synthesis, comparison and regional trends in population, land cover and climate](https://doi.org/10.1177/0959683619826697). The Holocene, 0959683619826697. 

Robinson, E., Nicholson, C., & Kelly, R. L. 2019. [The Importance of Spatial Data to Open-Access National Archaeological Databases and the Development of Paleodemography Research](https://doi.org/10.1017/aap.2019.29). Advances in Archaeological Practice, 1–14. doi:10.1017/aap.2019.29

Robinson, E., Bocinsky, R. K., Bird, D., Freeman, J., & Kelly, R. L. (2021). [Dendrochronological dates confirm a Late Prehistoric population decline in the American Southwest derived from radiocarbon dates](https://doi.org/10.1098/rstb.2019.0718). Philosophical Transactions of the Royal Society B: Biological Sciences, 376(1816), 20190718. 

Roscoe, P., Sandweiss, D. H., & Robinson, E. 2021. [Population density and size facilitate interactive capacity and the rise of the state](https://doi.org/10.1098/rstb.2019.0725). Philosophical Transactions of the Royal Society B: Biological Sciences, 376(1816), 20190725. 

Saag, L., Laneman, M., Varul, L., Malve, M., Valk, H., Razzak, M. A., et al. 2019. [The Arrival of Siberian Ancestry Connecting the Eastern Baltic to Uralic Speakers further East](https://doi.org/10.1016/j.cub.2019.04.026). Current Biology,

Schauer, P., Shennan, S., Bevan, A., Cook, G., Edinborough, K., Fyfe, R., et al. 2019. [Supply and demand in prehistory? Economics of Neolithic mining in northwest Europe](https://doi.org/10.1016/j.jaa.2019.03.001). Journal of Anthropological Archaeology, 54, 149–160. doi:10.1016/j.jaa.2019.03.001

Schauer, P., Shennan, S., Bevan, A., Colledge, S., Edinborough, K., Kerig, T., & Pearson, M. P. (2021). [Cycles in Stone Mining and Copper Circulation in Europe 5500–2000 bc: A View from Space](https://doi.org/10.1017/eaa.2020.56). European Journal of Archaeology, 24(2), 204–225. 

Schauer, P., Bevan, A., Shennan, S., Edinborough, K., Kerig, T., Pearson, M.P., 2019. [British Neolithic Axehead Distributions and Their Implications](https://doi.org/10.1007/s10816-019-09438-6). Journal of Archaeological Method and Theory. doi:/10.1007/s10816-019-09438-6

Schirrmacher, J., Kneisel, J., Knitter, D., Hamer, W., Hinz, M., Schneider, R.R., Weinelt, M., 2020. [Spatial patterns of temperature, precipitation, and settlement dynamics on the Iberian Peninsula during the Chalcolithic and the Bronze Age.](https://doi.org/10.1016/j.quascirev.2020.106220) Quaternary Science Reviews 233, 106220. 

Schroeter, N., Mingram, J., Kalanke, J., Lauterbach, S., Tjallingii, R., Schwab, V. F., & Gleixner, G. (2021). [The Reservoir Age Effect Varies With the Mobilization of Pre-Aged Organic Carbon in a High-Altitude Central Asian Catchment](https://doi.org/10.3389/feart.2021.681931). Frontiers in Earth Science, 0. 

Seidensticker, D., Hubau, W., Verschuren, D., Fortes-Lima, C., Maret, P. de, Schlebusch, C. M., & Bostoen, K. 2021. [Population collapse in Congo rainforest from 400 CE urges reassessment of the Bantu Expansion](https://doi.org/10.1126/sciadv.abd8352). Science Advances, 7(7), eabd8352. 

Solheim, S., 2020. [Mesolithic coastal landscapes : Demography, settlement patterns and subsistence economy in southeastern Norway](https://doi.org/10.4324/9780203730942-4). In Schülke. A (Eds.) ,Coastal Landscapes of the Mesolithic, Taylor and Francis.

Stoddart, S., Woodbridge, J., Palmisano, A., Mercuri, A. M., Mensing, S. A., Colombaroli, D., et al. 2019. [Tyrrhenian central Italy: Holocene population and landscape ecology.](https://doi.org/10.1177/0959683619826696) The Holocene, 0959683619826696. doi:10.1177/0959683619826696

Szpak, P., Savelle, J. M., Conolly, J., & Richards, M. P. 2019. [Variation in late holocene marine environments in the Canadian Arctic Archipelago: Evidence from ringed seal bone collagen stable isotope compositions](https://doi.org/10.1016/j.quascirev.2019.03.016). Quaternary Science Reviews, 211, 136–155. doi:10.1016/j.quascirev.2019.03.016

Tallavara, M., Personen, P. 2018. [Human ecodynamics in the north-west coast of Finland 10,000–2000 years ago](https://doi.org/10.1016/j.quaint.2018.06.032).Quaternary International,doi:/10.1016/j.quaint.2018.06.032

Tallavaara, M., & Jørgensen, E. K. 2021. [Why are population growth rate estimates of past and present hunter–gatherers so different?](https://doi.org/10.1098/rstb.2019.0708) Philosophical Transactions of the Royal Society B: Biological Sciences, 376(1816), 20190708. 

Thompson, A. E., Feinman, G. M., Lemly, M., & Prufer, K. M. 2021. [Inequality, networks, and the financing of Classic Maya political power](https://doi.org/10.1016/j.jas.2021.105441). Journal of Archaeological Science, 133, 105441. 

Thomson, M.J., MacDonald, G.M., 2020. [Climate and growing season variability impacted the intensity and distribution of Fremont maize farmers during and after the Medieval Climate Anomaly based on a statistically downscaled climate model](https://doi.org/10.1088/1748-9326/aba57e). Environ. Res. Lett.

Van der Bilt, W. G. M., Born, A., & Haaga, K. A. 2019. [Was Common Era glacier expansion in the Arctic Atlantic region triggered by unforced atmospheric cooling?](https://doi.org/10.1016/j.quascirev.2019.07.042) Quaternary Science Reviews, 105860. doi:10.1016/j.quascirev.2019.07.042

Vander Linden, M., & Silva, F. 2021. [Dispersals as demographic processes: testing and describing the spread of the Neolithic in the Balkans](https://doi.org/10.1098/rstb.2020.0231). Philosophical Transactions of the Royal Society B: Biological Sciences, 376(1816), 20200231. 

Vogels, O., Fäder, E., Lenssen-Erz, T., 2020. [A matter of diversity? Identifying past hunter-gatherer aggregation camps through data driven analyses of rock art sites.](https://doi.org/10.1016/j.quaint.2020.05.057) Quaternary International. 

Weiberg, E., Bevan, A., Kouli, K., Katsianis, M., Woodbridge, J., Bonnier, A., et al. 2019. [Long-term trends of land use and demography in Greece: A comparative study](https://doi.org/10.1177/0959683619826641) The Holocene, 0959683619826641. doi:10.1177/0959683619826641

Weitzel, E.M., Codding, B.F., Carmody, S.B., Zeanah, D.W., 2020. [Food Production and Domestication Produced Both Cooperative and Competitive Social Dynamics in Eastern North America.](https://doi.org/10.1080/14614103.2020.1737394) Environmental Archaeology 0, 1–14. 

Wilson, K. M., & Hill, M. G. (2020). [Synthesis and assessment of the flat-headed peccary record in North America](https://doi.org/10.1016/j.quascirev.2020.106601). Quaternary Science Reviews, 248, 106601. 

Woodbridge, J., Roberts, C. N., Palmisano, A., Bevan, A., Shennan, S., Fyfe, R., et al. 2019. [Pollen-inferred regional vegetation patterns and demographic change in Southern Anatolia through the Holocene.](https://doi.org/10.1177/0959683619826635) The Holocene, 0959683619826635. doi:10.1177/0959683619826635




