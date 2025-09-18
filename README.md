# Preoperative_fMRI_AscendingArousal
The following repository includes all relevant code for the analysis corresponding to the following paper "Preoperative Cholinergic Signatures Drive Segregated Brain Architecture in Postoperative Delirium" written by Dr. Natasha L Taylor et al.

### fMRI Analysis

All relevant for fMRI Analysis - note some of this code is overlaping with the following [repository](https://github.com/NatashaLTaylor/Preoperative_fMRI_FunctionalNetworks/tree/main/FC_Analysis) <br>


### Graph Theory Analysis
All graph theory analysis was conduced with Matlab V.2024b and can be found in [Graph Theory](/GraphTheory/) <br>
note - that relevant package is required for running graph theory calculations https://sites.google.com/site/bctnet/ <br>

**Cross-correlation** analysis to determine the relationship between peaks in phasic activity of the noradrenergic and cholinergic system can be found in [Graph Theory](/GraphTheory/ASS_Cross_Corr.m)


### Dynamic Functional Connectivity & LDA
Dynamic Functional connectivity was calculated from the following paper https://pubmed.ncbi.nlm.nih.gov/26231247/ and relevant code for calculating [MTD](https://github.com/macshine/coupling)<br>
Related code for calculating the Linear Discriminant Analysis to differentiate spatial-temporal brain dynamics unique to delirium and non-delirious [LDA for dFC](/dFC_LDA/).

### Logistic regression
All relevant code for logistic regression was done in R - see relevant packages required in [logistic regression code](/LogisticRegressions/).

### Figures
All generated figures can be foud in [Figures](/Figures/).



