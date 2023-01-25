## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  out.width = '70%', dpi = 300,
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(poolHelper)

## ----number of pools, message=FALSE, warning=FALSE, tidy=TRUE-----------------
# create a list with a single pool of 100 individuals
pools <- list(100)
# compute average absolute difference between allele frequencies
onePool <- maePool(nDip = 100, nloci = 1000, pools = pools, pError = 100, sError = 0.01, mCov = 100, vCov = 250, min.minor = 0)

# create a list with 10 pools, each with 10 individuals
pools <- list(rep(10, 10))
# compute average absolute difference between allele frequencies
tenPool <- maePool(nDip = 100, nloci = 1000, pools = pools, pError = 100, sError = 0.01, mCov = 100, vCov = 250, min.minor = 0)

# combine both  
final <- rbind(onePool, tenPool)
# convert the number of individuals in the pool to a factor
final$nPools <- as.factor(final$nPools)

# load the ggplot package
library(ggplot2)
# MAE value in the y-axis and the number of individuals in the pool in the x-axis
ggplot(final, aes(x = nPools, y = absError)) + 
  geom_boxplot() + theme_classic()

## ----coverage, message=FALSE, tidy=TRUE---------------------------------------
# create a vector with various mean coverages
mCov <- c(20, 50, 100)
# create a  vector with the variance of the coverage
vCov <- c(100, 250, 500)
# compute average absolute difference between allele frequencies
mydf <- maeFreqs(nDip = 100, nloci = 1000, pError = 100, sError = 0.01, mCov, vCov, min.minor = 0)

# convert the mean coverage into a factor
mydf$mean <- as.factor(mydf$mean)
# boxplot the MAE value in the y-axis and the coverage in the x-axis
ggplot(mydf, aes(x = mean, y = absError)) +
  geom_boxplot() + theme_classic()

## ----pool size, message=FALSE, tidy=TRUE--------------------------------------
# create a vector with various mean coverages
nDip <- c(10, 50, 100)

# compute average absolute difference between allele frequencies
mydf <- maeFreqs(nDip = nDip, nloci = 1000, pError = 100, sError = 0.01, mCov = 100, vCov = 250, min.minor = 0)

# convert the number of individuals into a factor
mydf$nDip <- as.factor(mydf$nDip)
# boxplot the MAE value in the y-axis and the coverage in the x-axis
ggplot(mydf, aes(x = nDip, y = absError)) +
  geom_boxplot() + theme_classic()

## ----combinations, fig.height=3, fig.width=5, message=FALSE, warning=FALSE----
# create a vector with various mean coverages
mCov <- c(20, 50, 100)
# create a  vector with the variance of the coverage
vCov <- c(100, 250, 500)
# create a vector with various pool errors 
pError <- c(5, 100, 250)

# compute average absolute difference between allele frequencies
mydf <- maeFreqs(nDip = 100, nloci = 1000, pError, sError = 0.01, mCov, vCov, min.minor = 0)

# convert the mean coverage into a factor
mydf$mean <- as.factor(mydf$mean)
# convert the pooling error to a factor
mydf$PoolError <- as.factor(mydf$PoolError)
  
# boxplot the MAE value in the y-axis and the pool error in the x-axis
# producing one boxplot for each of the different coverages
ggplot(mydf, aes(x = PoolError, y = absError, fill = mean)) +
  geom_boxplot() + theme_classic()

## ----reads one population-----------------------------------------------------
# simulate number of reads for one population
reads <- simulateCoverage(mean = 50, variance = 250, nSNPs = 100, nLoci = 1)
# display the structure of the reads object
str(reads)

## ----reads two populations----------------------------------------------------
# create a vector with the mean coverage of each population
mcov <- c(50, 100)
# create a vector with the variance of the coverage for each population
vcov <- c(250, 500)
# simulate number of reads for two populations
reads <- simulateCoverage(mean = mcov, variance = vcov, nSNPs = 100, nLoci = 1)
# display the structure of the reads object
str(reads)

## ----plot reads---------------------------------------------------------------
# create a vector with the mean coverage of each population
mcov <- c(50, 100)
# create a vector with the variance of the coverage for each population
vcov <- c(250, 500)
# simulate number of reads for two populations
reads <- simulateCoverage(mean = mcov, variance = vcov, nSNPs = 10000, nLoci = 1)
# plot the coverage of the first population
hist(reads[[1]][1,], col = rgb(0,0,1,1/4), xlim = c(0, 200), main = "", xlab = "")
# add the coverage of the second population
hist(reads[[1]][2,], col = rgb(1,0,0,1/4), add = TRUE)

## ----remove reads-------------------------------------------------------------
# check the minimum and maximum coverage before removal
x <- range(unlist(reads))
# remove sites with coverage below 25x and above 150x
reads <- remove_by_reads(nLoci = 1, reads = reads, minimum = 25, maximum = 150)
# display the structure of the reads object after removal
str(reads)
# check the minimum and maximum coverage after removal
range(unlist(reads))

## ----probability pool---------------------------------------------------------
# four pools with low sequencing error
poolProbs(nPools = 4, vector_np = c(10, 10, 10, 10), nSNPs = 6, pError = 5)
# four pools with high sequencing error
poolProbs(nPools = 4, vector_np = c(10, 10, 10, 10), nSNPs = 6, pError = 250)
# four pools but one is much larger
poolProbs(nPools = 4, vector_np = c(10, 100, 10, 10), nSNPs = 6, pError = 5)

## ----contribution pool low error----------------------------------------------
# simulate total coverage per site
reads <- unlist(simulateCoverage(mean = 100, variance = 250, nSNPs = 10, nLoci = 1))
# compute the proportion of contribution of each pool
probs <- poolProbs(nPools = 4, vector_np = rep(10, 4), nSNPs = 10, pError = 5)
# simulate the contribution in actual read numbers
pReads <- poolReads(nPools = 4, coverage = reads, probs = probs)
# output the number of reads per pool and per site 
pReads

## ----contribution pool high error---------------------------------------------
# simulate total coverage per site
reads <- unlist(simulateCoverage(mean = 100, variance = 250, nSNPs = 10, nLoci = 1))
# compute the proportion of contribution of each pool
probs <- poolProbs(nPools = 4, vector_np = rep(10, 4), nSNPs = 10, pError = 250)
# simulate the contribution in actual read numbers
pReads <- poolReads(nPools = 4, coverage = reads, probs = probs)
# output the number of reads per pool and per site 
pReads

## ----plot pool, tidy=TRUE-----------------------------------------------------
# simulate total coverage per site
reads <- simulateCoverage(mean = 100, variance = 250, nSNPs = 10000, nLoci = 1)
# unlist to create a vector with the coverage
reads <- unlist(reads)

# compute the proportion of contribution of each pool
probs <- poolProbs(nPools = 4, vector_np = rep(10, 4), nSNPs = 10000, pError = 5)
# simulate the contribution in actual read numbers
low.pReads <- poolReads(nPools = 4, coverage = reads, probs = probs)

# compute the proportion of contribution of each pool
probs <- poolProbs(nPools = 4, vector_np = rep(10, 4), nSNPs = 10000, pError = 250)
# simulate the contribution in actual read numbers
high.pReads <- poolReads(nPools = 4, coverage = reads, probs = probs)

# create the plot of the contribution with low pool error
h1 <- hist(unlist(low.pReads), plot = FALSE)
# create the plot of the contribution with high pool error
h2 <- hist(unlist(high.pReads), plot = FALSE)
# get the maximum x-value from the two plots
xmax <- max(h1[["breaks"]], h2[["breaks"]])
# and the maximum y-value
ymax <- max(h1[["counts"]], h2[["counts"]])
# set the color for the contribution computed with low pool error
col1 <- rgb(0,0,1,1/4)
# set the color for the contribution computed with high pool error
col2 <- rgb(1,0,0,1/4)

# plot the contribution computed with low pool error 
plot(h1, col = col1, xlim = c(0, xmax), ylim = c(0, ymax), main = "", xlab = "")
# add the plot of the contribution computed with high pool error
plot(h2, col = col2, add = TRUE)

## ----probability individual---------------------------------------------------
# compute the probability of contribution of each individual
indProbs(np = 10, nSNPs = 6, pError = 5)

## ----probability individual with higher error---------------------------------
# compute the probability of contribution of each individual
round(indProbs(np = 10, nSNPs = 5, pError = 150), digits = 5)

## ----individual reads low error-----------------------------------------------
# simulate total coverage per site
reads <- unlist(simulateCoverage(mean = 100, variance = 250, nSNPs = 12, nLoci = 1))
# compute the proportion of contribution of each pool
probs <- poolProbs(nPools = 4, vector_np = rep(10, 4), nSNPs = 12, pError = 5)
# simulate the contribution in actual read numbers
pReads <- poolReads(nPools = 4, coverage = reads, probs = probs)
# compute the proportion of contribution of each pool
probs <- indProbs(np = 10, nSNPs = 12, pError = 5)
# simulate the contribution in actual read numbers of each individual
indReads(np = 10, coverage = pReads[1,], probs = probs)

## ----individual reads high error----------------------------------------------
# simulate total coverage per site
reads <- unlist(simulateCoverage(mean = 100, variance = 250, nSNPs = 12, nLoci = 1))
# compute the proportion of contribution of each pool
probs <- poolProbs(nPools = 4, vector_np = rep(10, 4), nSNPs = 12, pError = 150)
# simulate the contribution in actual read numbers
pReads <- poolReads(nPools = 4, coverage = reads, probs = probs)
# compute the proportion of contribution of each pool
probs <- indProbs(np = 10, nSNPs = 12, pError = 150)
# simulate the contribution in actual read numbers of each individual
indReads(np = 10, coverage = pReads[1,], probs = probs)

## ----reference reads----------------------------------------------------------
# simulate total coverage per site
reads <- unlist(simulateCoverage(mean = 100, variance = 250, nSNPs = 12, nLoci = 1))
# compute the proportion of contribution of each pool
probs <- poolProbs(nPools = 4, vector_np = rep(10, 4), nSNPs = 12, pError = 5)
# simulate the contribution in actual read numbers
pReads <- poolReads(nPools = 4, coverage = reads, probs = probs)
# compute the proportion of contribution of each pool
probs <- indProbs(np = 10, nSNPs = 12, pError = 5)
# simulate the contribution in actual read numbers of each individual
iReads <- indReads(np = 10, coverage = pReads[1,], probs = probs)
# create fake genotypes - half the matrix is 0 and the other half is 2
geno <- rbind(matrix(0, nrow = 5, ncol = 12), matrix(2, nrow = 5, ncol = 12))
# simulate the number of reference reads
computeReference(genotypes = geno, indContribution = iReads, error = 0.001)

