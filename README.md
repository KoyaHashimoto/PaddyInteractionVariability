# R codes for Hashimoto et al. "Multifaceted effects of variable biotic interactions on population stability in complex interaction webs"

This repository assembles raw data csv files, markdown reports with source codes of R language, tables and figures for Hashimoto et al. "Multifaceted effects of variable biotic interactions on population stability in complex interaction webs". Readers may reproduce the results of the manuscript by running the source codes step-by-step.

### Contents of markdown reports

1. Data processing
2. Treatment effects on population density testetd by LMMs (including the results of Table S2)
3. Calculation of population sensitivity (including Fig. 2)
4. EDM
    1. Determine best embedding dimensions
    2. CCM (in 4_2_CCM directory) (including the results of Table S3)
        1. Test convergence in CCM analysis
        2. Test significance by surrogate analysis
        3. Determine optimal time lags by lagged-CCM
    3. S-map
        1. Estimation of per capita interaction effects by regularized S-map (including Fig. S3)
        2. Regularized S-map diagnostics
        3. Overview of S-map coefficient results (including Fig. S2 and S5)
5. Heterogeneity in interaction variability among different links (including Fig. 3)
6. Effects of interaction variability on population sensitivity (including Table 1 and Fig. 4 and S4)
7. Pesticide temporal dynamics (including Fig. S1)

### Source codes of regularized S-map (Cenci et al. 2019)

Regularized S-map, originally introduced by Cenci et al. 2019, was implemented here by using the custom codes from M. Ushio's github repository (https://github.com/ong8181/random-codes).