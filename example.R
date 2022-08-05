## Example of how the algorithm works
## Arguments can optionally change and optimize (criteria, learning rates, etc.)
## Assumption about number of pka values is needed (but can be any number and can be tried for different numbers)
## Depending on the randomization, with more pka values the algorithm is more prone to mistakes
## This happens if randomized pKa values are similar and not well covered in pH range or spectra are too similar
## This shouldn't be a problem since experiment design is meant to avoid this issues at most cases

## Import functions in written in 4 scripts
source('CF_and_derivatives.R')
source('CI_data_reduction.R')
source('spectra_sythesis_conc_matrix_new.R')
source('pKa_ML.R')

## set.seed is pseudorandom seed
set.seed(10)

## Randomly (uniform) create 50 pH points in range 3-13
pH <- runif(50,3,13)

## Create synthetic spectra (randomly as gaussians, more details in spectra_sythesis_conc_matrix_new.R)
## This spectra have 2 pka values, in range 4-12
## IF WANT TO TEST FOR NOISE, ERRORS, ETC. DO IT AFTER THIS LINE
example <- spectras(2, pH, pK=runif(2,4,12))

## Use PCA to reduce data to standardizes PCA scores that explain 99.9% of the variance
example1 <- reduce(example$spectra, crit = 0.999)

## Plot the spectra (assuming lambda range 200-800 nm)
matplot(y=t(example$spectra),x=200:800, typ="l")

## Save actual pKa values
pK=example$pK

## Run the algorithm, arguments are spectra, pH values, number of pKa values (2), number of iterations (1000),
## criterion for termination (10^-6) and alpha as learning rate for gradient descent (1)
try  <- ML_pka(example1, pH, 2, n_it=1000, crit=10^-6, alpha = 1)

## Print actual values
pK

## Print confidence intervals of determined pKa values; argument is output of ML_pka
try$pK
CI_ML(try)




#### rhodamine B - example
set.seed(1)

## Load data of spectra and pH values (change path if needed)
rb_x <- read.csv("RB.csv", header=FALSE, sep=";")
pH_rb=c(2.03,2.21,2.48,2.64,3.07,3.54,3.84,4.46,5.35,5.67) ## pH of experiment
rb_x <- t(rb_x) # Here transpose is needed

## Reduce spectra with PCA (99.9% ex. variance)
rb_x_pca <- reduce(rb_x, crit = 0.999)

## Platting spectra
matplot(y=t(rb_x),x=300:650, typ="l")

## Running the algorithm for 1 pKa
rb_try  <- ML_pka(rb_x_pca, pH_rb, 1,n_it=10000, crit=10^-9)

## Showing the results and confidence intervals
rb_try$pK
CI_ML(rb_try)

## Plotting the error (J) how it deceases with iterations
plot(rb_try$J~seq(1:length(rb_try$J)))
