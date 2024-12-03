# Set up Pibble model

# interaction model
f <- reformulate(termlabels=c("visit*PFS12*combiIO", 
                              "visit*PFS12*colitis",
                              "visit*PFS12*ppi",
                              "centre", 
                              "time_since_IO", 
                              "antibiotics", 
                              "toxicity", 
                              "previous_therapy", 
                              "age_c", 
                              "sex", 
                              "bmi_c", 
                              "patientid"))

# Create the design matrix
X <- t(model.matrix(f, data=mdat))

# Because we skip the multinomial part of the Pibble model, we have to deal with the 0s
Y <- ftbl # it is already close compositions (i.e. sum constraint to 1)
Y <- as.matrix(zCompositions::cmultRepl(Y, z.delete = F, z.warning = F)) # multiplicative 0 replacement 
Y <- t(Y) # pibble assumes features on rows

N <- ncol(Y) # number of samples
D <- nrow(Y) # number of features (categories)

# Specify uninformative (default) priors 
upsilon <- D+3
Omega <- diag(D)
G <- cbind(diag(D-1), -1)
Xi <- (upsilon-D)*G%*%Omega%*%t(G)
Theta <- matrix(0, D-1, nrow(X))
Gamma <- diag(nrow(X))

priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)

# Fit the Pibble model
eta_init <- t(driver::alr(t(Y))) # The model implemented using the ALR transform as it is computationally simple and fast; the results of the model can be viewed as if any number of transforms had been used, i.e. we can switch between any.
eta_array <- array(eta_init, dim=c(nrow(eta_init), ncol(eta_init), 2000))
# uncollapsePibble can be fitted on proportions (normally, this is an intermediate step; if we had counts, we would start with the optimPibbleCollapsed function)
posterior <- uncollapsePibble(eta_array, priors$X, priors$Theta, priors$Gamma, priors$Xi, priors$upsilon, seed=2849)

# Attach dimnames
dimnames(posterior$Lambda)[[2]] <- rownames(X)
dimnames(posterior$Lambda)[[1]] <- rownames(Y)[-length(rownames(Y))]
dimnames(posterior$Sigma)[[1]] <- dimnames(posterior$Sigma)[[2]] <- rownames(Y)[-length(rownames(Y))]

posterior <- pibblefit(D=D,
                       N=N,
                       Q=nrow(X),
                       coord_system="alr",
                       iter=2000L,
                       alr_base=D,
                       Eta=eta_array,
                       Lambda=posterior$Lambda,
                       Sigma=posterior$Sigma,
                       Y=Y,
                       X=X,
                       names_categories=rownames(Y),
                       names_samples=colnames(Y),
                       names_covariates=rownames(X),
                       upsilon=priors$upsilon,
                       Theta=priors$Theta,
                       Gamma=priors$Gamma,
                       Xi=priors$Xi)

# Change to CLR transform
fit_species_clr <- to_clr(posterior)

# Attach dimnames
dimnames(fit_species_clr$Lambda)[[2]] <- rownames(X)
dimnames(fit_species_clr$Lambda)[[1]] <- rownames(Y)
dimnames(fit_species_clr$Sigma)[[1]] <- dimnames(fit_species_clr$Sigma)[[2]] <- rownames(Y)




