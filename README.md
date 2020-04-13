# Time series analysis and modeling of the Covid-19 pandemic dynamics
## [Marco Arieli Herrera-Valdez](https://mahv13.wordpress.com)
### [Laboratorio de Fisiología de Sistemas](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=5&cad=rja&uact=8&ved=2ahUKEwi9p4KJidroAhUMi6wKHYrSBWcQFjAEegQIAhAB&url=https%3A%2F%2Fmarcoh48.wixsite.com%2Ffisiologiasistemasfc&usg=AOvVaw1RFgV1gOqxbpBJT3Bl6WEq), [Facultad de Ciencias](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=2ahUKEwjbiNnQrtvoAhUJA6wKHVI0BXMQFjAAegQIGRAD&url=http%3A%2F%2Fwww.fciencias.unam.mx%2F&usg=AOvVaw1dMRMU_F-IcpmaB1y1H4px), [UNAM](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=2ahUKEwivy6_irtvoAhUDaq0KHQVoCcAQFjAAegQIGhAD&url=https%3A%2F%2Fwww.unam.mx%2F&usg=AOvVaw0YWCGJ7FEpDwkcT3EYH-aM)
(Last modified: MAHV, Apri 12, 2020)


#### 1. Fatality associated to reported cases during the first weeks of the COVID-19 pandemic 

(Herrera-Nolasco CI, Herrera-McKiernan AJ, Herrera-McKiernan EA, O'Reilly-Regueiro E, Herrera-Valdez MA, *In preparation*)

(a) Illustrating figures.

- [Delays between case reports, and between first death and first recovery with respect to the first report in all countries](tsam_COVID19_JHU_delaysAllCountries.png)

- [Cases vs deaths per 1000 habitants for some countries](tsam_COVID19_JHU_cases-deaths_x1000000_JHU.png)

- [CFRs for some countries from the first case in China](tsam_COVID19_cfr_JHU_fromFirstCaseInChina.png)

- [CFRs for some countries from the day the first case was reported locally](tsam_COVID19_JHU_cfr_fromFirstLocalCase.png)

- [CFRs for China taking into account local reports by provinces](tsam_COVID19_JHU_cfr_ProvincesChina_fromFirstLocalReport.png) Note the differences between the CFR from the total cases and the average CFR from the provinces. 

- [CFRs from China, South Korea, and Italy, by age](tsam_COVID19_JHU_cfr+propDeathCases_ByAge_China+SKorea+Italy_OneFigure.png) Note the similarities between the weights of the CFRs by age from the different countries and consider their possible divergences in terms of the age structure of their populations, access to care, testing for COVID-19, etcetera. 

- [Estimation of death cases in Mexico for different age groups with a sub report factor of 1](tsam_COVID19_JHU_cfr+propDeathCasesByAgeTS_EstimatesMexico_subReportFactor1.png)

- [Estimation of death cases in Mexico for different age groups with a sub report factor of 10](tsam_COVID19_JHU_cfr+propDeathCasesByAgeTS_EstimatesMexico_subReportFactor10.png)

Delay and CFR analysis of case, death, and recovery reports from around the world. 

(b) Code 
[JuPyTeR notebook](tsam_COVID-19_JHU_cfr_Jan2020-.ipynb), its [html version](tsam_COVID-19_JHU_cfr_Jan2020-.html), and the [Python  code](tsam_COVID-19_JHU_cfr_Jan2020-.py). The [base code with general functions](tsam_COVID19_baseCode.py) 

