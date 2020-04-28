#dam_COVID19_cfrAnalysis_ReadMe.md

## Time series analysis from global and local data

### Fatality associated to reported cases during the first weeks of the Covid-19 pandemic

(Herrera-Nolasco CI, Herrera-McKiernan AJ, Herrera-McKiernan EA, O'Reilly-Regueiro E, Herrera-Valdez MA, *In preparation*)

(a) Delay and CFR analysis of case, death, and recovery reports from around the world.
Illustrating figures.

- [Delays between case reports, and between first death and first recovery with respect to the first report in all countries](../figures_COVID19_dataAnalysis/dam_COVID19_JHU_delays_caseDeaths.png)

- [Cases vs deaths per 1000 habitants for some countries](../figures_COVID19_dataAnalysis/dam_COVID19_JHU_cases-deaths_x1000000_JHU.png)

- [CFRs for some countries from the first case in China](../figures_COVID19_dataAnalysis/dam_COVID19_cfr_JHU_fromFirstCaseInChina.png)

- [CFRs for some countries from the day the first case was reported locally](../figures_COVID19_dataAnalysis/dam_COVID19_JHU_cfr_fromFirstLocalCase.png)

- [CFRs for China taking into account local reports by provinces](../figures_COVID19_dataAnalysis/dam_COVID19_JHU_cfr_ProvincesChina_fromFirstLocalReport.png) Note the differences between the CFR from the total cases and the average CFR from the provinces.

- [CFRs from China, South Korea, and Italy, by age](../figures_COVID19_dataAnalysis/dam_COVID19_JHU_cfr+propDeathCases_ByAge_China+SKorea+Italy_OneFigure.png) Note the similarities between the weights of the CFRs by age from the different countries and consider their possible divergences in terms of the age structure of their populations, access to care, testing for Covid-19, etcetera.

**Inference of death cases for Mexico**
- [Estimation of death cases in Mexico for different age groups with a sub report factor of 1](../figures_COVID19_dataAnalysis/dam_COVID19_JHU_cfr+propDeathCasesByAgeTS_EstimatesMexico_subReportFactor1.png)

- [Estimation of death cases in Mexico for different age groups with a sub report factor of 10](../figures_COVID19_dataAnalysis/dam_COVID19_JHU_cfr+propDeathCasesByAgeTS_EstimatesMexico_subReportFactor10.png)

(b) Code
[JuPyTeR notebook](./dam_COVID19_JHU_cfr_Jan2020-.ipynb), its [html version](./dam_COVID19_JHU_cfr_Jan2020-.html), and the [Python  code](./dam_COVID19_JHU_cfr_Jan2020-.py). The [base code with general functions](./dam_COVID19_baseCode.py)
