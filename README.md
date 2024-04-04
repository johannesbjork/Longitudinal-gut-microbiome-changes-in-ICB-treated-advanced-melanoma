
# Index
** fit_pibble_model.R
   * The Logistic-Normal model (as implemented in Pibble, fido package) with three 3-way interactions


  * * R_vs_NR.R
    * Function for computing, from the fitted model, the marginal average for PFS>=12 (R) vs PFS<12 (NR) averaging across all levels of W1 (therapy regimen), W2 (colitis) and W3 (PPI-use), and plotting the posteriors for taxa whose 90% (can be adjusted) credible interval (CI) does not cover 0


  * * R_vs_NR_monotherapy.R
    * Function for computing and plotting the marginal average for PFS12>=12 (R) vs PFS<12 (NR) on monotherapy (W1=0), no colitis (W2=0) and no PPIs (W3=0)


  * * R_vs_NR_combitherapy.R
    * Function for computing and plotting the marginal average for PFS12>=12 (R) vs PFS<12 (NR) on combination therapy (W1=1), no colitis (W2=0) and no PPIs (W3=0)


  * * colitisyes_vs_colitisno.R
    * Function for computing and plotting the marginal average for colitis yes vs no, averaging over Z (PFS12), (therapy regimen) and W3 (PPI-use)

  * * survival_analysis.R 
    * Code for the survival analysis, including the calculation of the longitudinal balance

   
For the underlying calculations, please see the Supplementary Methods in Björk et al. 2024, Nature Medicine.

Data can be found on Open Science Framework: https://osf.io/rakdf/
