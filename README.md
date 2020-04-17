# Time series analysis and modeling of the dynamics of Covid-19
## [Marco Arieli Herrera-Valdez](https://mahv13.wordpress.com)
### [Laboratorio de Fisiología de Sistemas](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=5&cad=rja&uact=8&ved=2ahUKEwi9p4KJidroAhUMi6wKHYrSBWcQFjAEegQIAhAB&url=https%3A%2F%2Fmarcoh48.wixsite.com%2Ffisiologiasistemasfc&usg=AOvVaw1RFgV1gOqxbpBJT3Bl6WEq), [Facultad de Ciencias](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=2ahUKEwjbiNnQrtvoAhUJA6wKHVI0BXMQFjAAegQIGRAD&url=http%3A%2F%2Fwww.fciencias.unam.mx%2F&usg=AOvVaw1dMRMU_F-IcpmaB1y1H4px), [UNAM](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=2ahUKEwivy6_irtvoAhUDaq0KHQVoCcAQFjAAegQIGhAD&url=https%3A%2F%2Fwww.unam.mx%2F&usg=AOvVaw0YWCGJ7FEpDwkcT3EYH-aM)

#### 1. Estimating the starting date, first wave peaks and relative contributions of subpopulations with different clinical profiles to the Covid-19 pandemic in Mexico
(Herrera-Nolasco CI, Herrera-Valdez MA, *In preparation*)

(a) [Estimating the starting time and peak of the Covid-19 epidemic in Mexico](tsam_Covid19_models/figures/Covid19_Mexico_InitialFit_Herrera-Valdez+Herrera-Nolasco_2020.png) 


#### 2. Fatality associated to reported cases during the first weeks of the Covid-19 pandemic

(Herrera-Nolasco CI, Herrera-McKiernan AJ, Herrera-McKiernan EA, O'Reilly-Regueiro E, Herrera-Valdez MA, *In preparation*)

(a) Delay and CFR analysis of case, death, and recovery reports from around the world.
Illustrating figures.

- [Delays between case reports, and between first death and first recovery with respect to the first report in all countries](tsam_Covid19_figures/tsam_Covid19_JHU_delaysAllCountries.png)

- [Cases vs deaths per 1000 habitants for some countries](tsam_Covid19_figures/tsam_Covid19_JHU_cases-deaths_x1000000_JHU.png)

- [CFRs for some countries from the first case in China](tsam_Covid19_figures/tsam_Covid19_cfr_JHU_fromFirstCaseInChina.png)

- [CFRs for some countries from the day the first case was reported locally](tsam_Covid19_figures/tsam_Covid19_JHU_cfr_fromFirstLocalCase.png)

- [CFRs for China taking into account local reports by provinces](tsam_Covid19_figures/tsam_Covid19_JHU_cfr_ProvincesChina_fromFirstLocalReport.png) Note the differences between the CFR from the total cases and the average CFR from the provinces.

- [CFRs from China, South Korea, and Italy, by age](tsam_Covid19_figures/tsam_Covid19_JHU_cfr+propDeathCases_ByAge_China+SKorea+Italy_OneFigure.png) Note the similarities between the weights of the CFRs by age from the different countries and consider their possible divergences in terms of the age structure of their populations, access to care, testing for Covid-19, etcetera.

**Inference of death cases for Mexico**
- [Estimation of death cases in Mexico for different age groups with a sub report factor of 1](tsam_Covid19_figures/tsam_Covid19_JHU_cfr+propDeathCasesByAgeTS_EstimatesMexico_subReportFactor1.png)

- [Estimation of death cases in Mexico for different age groups with a sub report factor of 10](tsam_Covid19_figures/tsam_Covid19_JHU_cfr+propDeathCasesByAgeTS_EstimatesMexico_subReportFactor10.png)

(b) Code
[JuPyTeR notebook](tsam_Covid19_cfrAnalysis/tsam_Covid19_JHU_cfr_Jan2020-.ipynb), its [html version](tsam_Covid19_cfrAnalysis/tsam_Covid19_JHU_cfr_Jan2020-.html), and the [Python  code](tsam_Covid19_cfrAnalysis/tsam_Covid19_JHU_cfr_Jan2020-.py). The [base code with general functions](tsam_Covid19_baseCode.py)

_(Created: MAHV, 20200301. Last modified: MAHV, 20200416)_
