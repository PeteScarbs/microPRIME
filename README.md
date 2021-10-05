# microPRIME
microPRIME is a microsimulation model that forecasts myocardial infarction rates in England based on trends in risk factors and treatments. Full details about the microPRIME model and forecasted incidence, prevalence and event rates for England to 2035 are reported in a paper that is currently under review:

Scarborough P, Kaur A, Cobiac LJ. Forecast of myocardial infarction incidence, events and prevalence in England to 2035 using a microsimulation model with endogenous disease outcomes. Under review.

This repository hosts all of the R scripts that combine to run the microPRIME model. A full description of how the microPRIME model operates is provided in the supplementary material of the manuscript under review, and also in a word document in this repository: 
- microPRIME description.docx

The other files are R scripts for each of the modules that make up microPRIME. The underlying datasets that are used by microPRIME are not free for sharing. Please contact the lead author (peter.scarborough@ndph.ox.ac.uk) for more information.

The R scripts that make up microPRIME are as follows:
1. risk factors.R
2. demo.R
3. parameter draw.R
4. agent fill.R
5. microsim.R
6. emulator.R
