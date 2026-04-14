# **Rethinking Spatial Units for Malaria Prediction in Loreto, Peru: The Use of Watershed-Based Boundaries in the Amazon Rainforest**

## Description
*Malaria remains a major public health threat in the Peruvian Amazon, where Loreto accounts for more than 76% of national cases. Most spatial modelling studies rely on administrative boundaries (districts) that follow arbitrary political divisions, poorly reflecting the ecological dynamics that drive transmission.
This study evaluates whether watershed-based spatial units (hydrobasins) outperform traditional administrative boundaries for predicting Plasmodium vivax and Plasmodium falciparum cases in Loreto from 2009 to 2018. Four boundary types are compared: administrative districts and three hydrological levels (HydroSHEDS L6, L7, and ANA watersheds). Modelling combines Zero-Inflated Negative Binomial (ZINB) regression with three machine learning algorithms — Random Forest, Support Vector Machine, and XGBoost — using PCA-derived climatic, hydrological, and vegetation predictors from satellite imagery.*


---
 
## Key Results
 
Hydrobasin level 7 (78 units) consistently achieved the lowest prediction error across all models and both malaria species.
 
### P. falciparum — RMSE by boundary and model
 
| Boundary | XGBoost | Random Forest | SVM | ZINB |
|---|---:|---:|---:|---:|
| **Hydrobasin L7** | **15.79** | **17.21** | **19.42** | **24.09** |
| Districts | 17.10 | 17.32 | 26.26 | 36.99 |
| Hydrobasin ANA | 33.15 | 38.00 | 45.58 | 138.61 |
| Hydrobasin L6 | 43.06 | 39.96 | 56.80 | 90.67 |
 
### P. vivax — RMSE by boundary and model
 
| Boundary | XGBoost | Random Forest | SVM | ZINB |
|---|---:|---:|---:|---:|
| **Hydrobasin L7** | **38.44** | **40.98** | **46.24** | **77.08** |
| Districts | 50.76 | 51.84 | 76.27 | 97.34 |
| Hydrobasin ANA | 68.99 | 79.10 | 95.32 | 215.61 |
| Hydrobasin L6 | 118.93 | 122.48 | 148.68 | 229.05 |
 
> For *P. falciparum*, hydrobasin L7 reduced RMSE by 24–68% relative to district and ANA boundaries. For *P. vivax*, reductions ranged from 28–74%. Spatiotemporal cross-validation (2009–2014 train / 2015–2018 test) confirmed these patterns, with lowest temporal RMSE at hydrobasin L7 for both species (28.67 and 73.85 respectively).
 
---

## Data Availability
*Malaria surveillance data from the Peruvian Ministry of Health are not publicly available. Access to anonymised data may be requested subject to institutional approval. Contact: imtavh.innovalab@oficinas-upch.pe*
