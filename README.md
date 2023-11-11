# MSCN_missing
Model-based clustering with missing data using the MSCN

Reference: H. Tong and, C. Tortora. Missing values and directional outlier detection 
in model-based clustering. Journal of Classification, 2023. 
https://rdcu.be/dpYfD

Abstract
Robust model-based clustering tackles the task of uncovering heterogeneity in data
sets characterized by outliers. However, the use of many methods in this area becomes
severely limited in applications where partially observed records are common since their
existing frameworks often assume complete data only. Here, a mixture of multiple scaled
contaminated normal (MSCN) distributions is extended using the expectation-conditional
maximization (ECM) algorithm to accommodate data sets with values missing at random.
The newly proposed extension preserves the mixture’s capability to yield robust parameter
estimates and perform automatic outlier detection separately for each principal
component. In this fitting framework, the MSCN marginal density is approximated using
the inversion formula for the characteristic function. 


- mscnm contains the code
- Example_MSCN shows an example on a simulated data set

This material is based upon work supported by the National Science Foundation under
Grant No. 2209974 (Tortora).
