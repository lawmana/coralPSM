# coralPSM
Coral data analysis and proxy system modeling algorithms in MATLAB.

## Publication
 Lawman, A. E., Partin, J. W., Dee, S. G., Casadio, C. A., Di Nezio, P., & Quinn, T. M. (2020). Developing a coral proxy system model to compare coral and climate model estimates of changes in paleo‐ENSO variability. *Paleoceanography and Paleoclimatology*, *35*, e2019PA003836. https://doi.org/10.1029/2019PA003836.
 
 ## Overview
 This repository includes the following coral PSM algorithms:
1. **Variations in growth rate** - growthrate.m
2. **Analytical and calibration errors** - analytical_error.m
3. **Age model (depth to time transformation)** - psAgeModel.m

## Function Descriptions
### Variations in Annual Growth Rates
*See Section 2.4 in our paper for more details*

The function **growthrate.m** function perturbs the independent vector of a data set using an autoregressive (AR) model.

A coral’s growth rate may vary both within and between years. For example, a coral growing an average of 1.2 cm/year would achieve approximately monthly resolution if sampled in 1 mm increments. Although monthly resolution is targeted, one sample of coral powder may average 2-3 weeks of time when the coral is growing faster, or 5-6 weeks when the coral is growing slower. Due to variable growth rates, the net effect of equal sampling in the depth domain will lead to unequal sampling in the time domain. **growthrate.m** is used to assess how variations in coral growth impact the variance of a resulting geochemical time series when the coral is sampled at a fixed sampling resolution (e.g., 1 mm). Our study uses an AR(2) model in which the lag1 and lag2 coefficients are based on the measured annual growth rates for select *Porites* corals.

### Analytical Errors
*See Section 2.5.1 in our paper for more details*

The function **analytical_error.m** models analytical errors (e.g., laboratory analytical precision) as Gaussian white noise.

### Coral Age Modeling
*See section 2.5.2 in our paper for more details*

This function **psAgeModel.m** interpolates data to even sampling in the time domain using peaks and troughs in the data as chronological tie points. This function assumes that the largest signal in the input data is of constant frequency, but that this signal is distorted in time.

An application of this algorithm is to convert coral geochemical data from the depth to the time domain (coral age modeling) for paleoclimate studies. For example, corals are often sampled at approximately monthly resolution, with the annual cycle emerging as a dominant signal. 

The age model algorithm identifies the local minima/maxima (critical points) in the raw geochemical data (depth or sample number domain) and uses them as chronological tie points when interpolating the data to the target temporal resolution (e.g., monthly = 12 points-per-year). The local minima/maxima are assigned a calendar month based on knowledge of the climatology at the coral study site.

#### Example Files: https://github.com/lawmana/coralPSM/tree/master/Examples/AgeModel
1. **ageModelDemo.m** demonstrates how to use the function using generic sinusoidal data as the input
2. **coralAgeModel_VanuatuFossilCoral_SrCa.m** generates a relative age model using sample coral Sr/Ca data from Vanuatu in the southwest Pacific
