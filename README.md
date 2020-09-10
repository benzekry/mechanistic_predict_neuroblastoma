# Mechanistic modeling of metastasis for prediction of overall survival in neuroblastoma

This code generates the results and figures reported in the paper:  
Benzekry, S., Sentis, C., Coze, C., Tessonnier, L., & André, N. (2020). Descriptive and prognostic value of a computational model of metastasis in high-risk neuroblastoma. JCO: Clinical Cancer Informatics, in revision.

Specifically, the code is composed of python and matlab scripts and jupyter notebooks (in python and R) and performs:
  - Statistical analysis of prognosis factors of overall and progression-free survival (Kaplan-Meier, log-rank, Cox regression)
```
code/statistical_analysis.ipynb
```
with results exported in `statistical_analysis/`
  - Simulation of a mechanistic model of metastasis
```
code/main_simulate.m
```
  - Calibration of the model parameters from quantitative clinical data at diagnosis:  primary tumor size, lactate dehydrogenase (LDH) and SIOPEN score from nuclear imaging
```
code/mechanistic.ipynb
```
with results exported in `mechanistic/`
  - Assessment of the predictive power of patient-specific Cox regression-based models for overall survival
```
code/mechanistic.ipynb
```
and
```
code/cox_r_calibration_plot.ipynb
```
with results exported in `cox_regression` (clinical data alone) and `mechanistic/cox_regression` (clinical data + mathematical parameters)

The data used for the analysis is available in ``data/``

Previous preprint:
Descriptive and prognostic value of a computational model of metastasis in high-risk neuroblastoma
Sebastien Benzekry, Coline Sentis, Carole Coze, Laetitia Tessonnier, Nicolas Andre
medRxiv 2020.03.26.20042192; doi:https://doi.org/10.1101/2020.03.26.20042192
