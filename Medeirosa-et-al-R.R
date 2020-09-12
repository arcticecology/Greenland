# Load the required packages
# (vegan must be loaded after ade4 to avoid some conflicts)
library(ape)
library(spdep)
library(ade4)
library(adegraphics)
library(adespatial)
library(vegan)
library(AEM)


# Load additionnal functions
# (files must be in the working directory)
source("plot.links.R")
source("sr.value.R")
source("quickMEM.R")

#######################################################################################################
#CADbaffin
#######################################################################################################
# Import the data from CSV files
# (files must be in the working directory)
chiro <- read.csv("CADspeNewbaffin.csv")
chiro.env <- read.csv("CADenvNewbaffin3.csv")
chiro.xy <- read.csv("CADspaNewbaffin.csv")
####################################################################################################
#####################################################################################################

# Transform the data
chiro.h <- decostand (chiro, "hellinger")
chiro.xy.c <- scale(chiro.xy, center=TRUE, scale=FALSE)

####################################################################################################
####################################################################################################



######MEM analisys###### 
## 1. Construct the matrix of dbMEM variables
chiro.dbmem.tmp <- dbmem(chiro.xy, silent=FALSE)
chiro.dbmem <- as.data.frame(chiro.dbmem.tmp)
# Truncation distance used above:
(thr <- give.thresh(dist(chiro.xy)))

###################################################################################################
###################################################################################################
# 2. Test and forward selection of the environmental variables
# Forward selection of the environmental variables
chiro.env.rda <- rda(chiro.h ~., chiro.env)
(chiro.env.R2a <- RsquareAdj(chiro.env.rda)$adj.r.squared)
chiro.env.fwd <- 
  forward.sel(chiro.h, chiro.env,
              adjR2thresh = chiro.env.R2a, 
              nperm = 9999)
env.sign <- sort(chiro.env.fwd$order)
env.red <- chiro.env[ ,c(env.sign)]
colnames(env.red)

# 3. Test and forward selection of the dbMEM variables
# Run the global dbMEM analysis on the *detrended* chiro data
chiro.det.dbmem.rda <- rda(chiro.h ~., chiro.dbmem)
anova(chiro.det.dbmem.rda)
# Since the analysis is significant, compute the adjusted R2
# and run a forward selection of the dbMEM variables
(chiro.det.dbmem.R2a <- 
    RsquareAdj(chiro.det.dbmem.rda)$adj.r.squared)
(chiro.det.dbmem.fwd <- 
    forward.sel(chiro.h, 
                as.matrix(chiro.dbmem), 
                adjR2thresh = chiro.det.dbmem.R2a))
# Number of significant dbMEM
(nb.sig.dbmem <- nrow(chiro.det.dbmem.fwd))
# Identify the significant dbMEM sorted in increasing order
(dbmem.sign <- sort(chiro.det.dbmem.fwd$order))
# Write the significant dbMEM to a new object (reduced set)
dbmem.red <- chiro.dbmem[ ,c(dbmem.sign)]


## 5. Chiro - environment - trend - dbMEM variation partitioning
(chiro.varpart <- varpart(chiro.h, env.red, dbmem.red))
dev.new(title="Chiro - environment - dbMEM variation partitioning", width=12, height=6)
par(mfrow=c(1,2))
showvarparts(4) # To show the symbols of the fractions
plot(chiro.varpart, digits=2)
anova(
  rda(chiro.h, env.red, dbmem.red))
# Tests of the unique fractions [a] and [b]
# Fraction [a], pure environmental
anova(rda(chiro.h, env.red, cbind(dbmem.red)))
# Fraction [b], pure scale spatial
anova(rda(chiro.h, dbmem.red, cbind(env.red)))
###############################################################################################################################
###############################################################################################################################
#ICELAND
#######################################################################################################
# Import the data from CSV files
# (files must be in the working directory)
chiro <- read.csv("ICEspeNew.csv")
chiro.env <- read.csv("ICEenvNew4.csv")
chiro.xy <- read.csv("ICEspaNew2.csv")
####################################################################################################
#####################################################################################################

# Transform the data
chiro.h <- decostand (chiro, "hellinger")
chiro.xy.c <- scale(chiro.xy, center=TRUE, scale=FALSE)

####################################################################################################
####################################################################################################



######MEM analisys###### 
## 1. Construct the matrix of dbMEM variables
chiro.dbmem.tmp <- dbmem(chiro.xy, silent=FALSE)
chiro.dbmem <- as.data.frame(chiro.dbmem.tmp)
# Truncation distance used above:
(thr <- give.thresh(dist(chiro.xy)))

###################################################################################################
###################################################################################################
# 2. Test and forward selection of the environmental variables
# Forward selection of the environmental variables
chiro.env.rda <- rda(chiro.h ~., chiro.env)
(chiro.env.R2a <- RsquareAdj(chiro.env.rda)$adj.r.squared)
chiro.env.fwd <- 
  forward.sel(chiro.h, chiro.env,
              adjR2thresh = chiro.env.R2a, 
              nperm = 9999)
env.sign <- sort(chiro.env.fwd$order)
env.red <- chiro.env[ ,c(env.sign)]
colnames(env.red)

# 3. Test and forward selection of the dbMEM variables
# Run the global dbMEM analysis on the *detrended* chiro data
chiro.det.dbmem.rda <- rda(chiro.h ~., chiro.dbmem)
anova(chiro.det.dbmem.rda)
# Since the analysis is significant, compute the adjusted R2
# and run a forward selection of the dbMEM variables
(chiro.det.dbmem.R2a <- 
    RsquareAdj(chiro.det.dbmem.rda)$adj.r.squared)
(chiro.det.dbmem.fwd <- 
    forward.sel(chiro.h, 
                as.matrix(chiro.dbmem), 
                adjR2thresh = chiro.det.dbmem.R2a))
# Number of significant dbMEM
(nb.sig.dbmem <- nrow(chiro.det.dbmem.fwd))
# Identify the significant dbMEM sorted in increasing order
(dbmem.sign <- sort(chiro.det.dbmem.fwd$order))
# Write the significant dbMEM to a new object (reduced set)
dbmem.red <- chiro.dbmem[ ,c(dbmem.sign)]


## 5. Chiro - environment - trend - dbMEM variation partitioning
(chiro.varpart <- varpart(chiro.h, env.red, dbmem.red))
dev.new(title="Chiro - environment - dbMEM variation partitioning", width=12, height=6)
par(mfrow=c(1,2))
showvarparts(4) # To show the symbols of the fractions
plot(chiro.varpart, digits=2)
anova(
  rda(chiro.h, env.red, dbmem.red))
# Tests of the unique fractions [a] and [b]
# Fraction [a], pure environmental
anova(rda(chiro.h, env.red, cbind(dbmem.red)))
# Fraction [b], pure scale spatial
anova(rda(chiro.h, dbmem.red, cbind(env.red)))
###############################################################################################################################
###############################################################################################################################

######################################################################################################
#GREENLAND
#######################################################################################################
# Import the data from CSV files
# (files must be in the working directory)
chiro <- read.csv("GREENspeNew.csv")
chiro.env <- read.csv("GREENenvNew7.csv")
chiro.xy <- read.csv("GREENspaNew2.csv")
####################################################################################################
#####################################################################################################

# Transform the data
chiro.h <- decostand (chiro, "hellinger")
chiro.xy.c <- scale(chiro.xy, center=TRUE, scale=FALSE)

####################################################################################################
####################################################################################################



######MEM analisys###### 
## 1. Construct the matrix of dbMEM variables
chiro.dbmem.tmp <- dbmem(chiro.xy, silent=FALSE)
chiro.dbmem <- as.data.frame(chiro.dbmem.tmp)
# Truncation distance used above:
(thr <- give.thresh(dist(chiro.xy)))

###################################################################################################
###################################################################################################
# 2. Test and forward selection of the environmental variables
# Forward selection of the environmental variables
chiro.env.rda <- rda(chiro.h ~., chiro.env)
(chiro.env.R2a <- RsquareAdj(chiro.env.rda)$adj.r.squared)
chiro.env.fwd <- 
  forward.sel(chiro.h, chiro.env,
              adjR2thresh = chiro.env.R2a, 
              nperm = 9999)
env.sign <- sort(chiro.env.fwd$order)
env.red <- chiro.env[ ,c(env.sign)]
colnames(env.red)

# 3. Test and forward selection of the dbMEM variables
# Run the global dbMEM analysis on the *detrended* chiro data
chiro.det.dbmem.rda <- rda(chiro.h ~., chiro.dbmem)
anova(chiro.det.dbmem.rda)
# Since the analysis is significant, compute the adjusted R2
# and run a forward selection of the dbMEM variables
(chiro.det.dbmem.R2a <- 
    RsquareAdj(chiro.det.dbmem.rda)$adj.r.squared)
(chiro.det.dbmem.fwd <- 
    forward.sel(chiro.h, 
                as.matrix(chiro.dbmem), 
                adjR2thresh = chiro.det.dbmem.R2a))
# Number of significant dbMEM
(nb.sig.dbmem <- nrow(chiro.det.dbmem.fwd))
# Identify the significant dbMEM sorted in increasing order
(dbmem.sign <- sort(chiro.det.dbmem.fwd$order))
# Write the significant dbMEM to a new object (reduced set)
dbmem.red <- chiro.dbmem[ ,c(dbmem.sign)]


## 5. Chiro - environment - trend - dbMEM variation partitioning
(chiro.varpart <- varpart(chiro.h, env.red, dbmem.red))
dev.new(title="Chiro - environment - dbMEM variation partitioning", width=12, height=6)
par(mfrow=c(1,2))
showvarparts(4) # To show the symbols of the fractions
plot(chiro.varpart, digits=2)
anova(
  rda(chiro.h, env.red, dbmem.red))
# Tests of the unique fractions [a] and [b]
# Fraction [a], pure environmental
anova(rda(chiro.h, env.red, cbind(dbmem.red)))
# Fraction [b], pure scale spatial
anova(rda(chiro.h, dbmem.red, cbind(env.red)))
###############################################################################################################################
###############################################################################################################################
#NUUK
#######################################################################################################
# Import the data from CSV files
# (files must be in the working directory)
chiro <- read.csv("NUKspeNew.csv")
chiro.env <- read.csv("NukenvNew9.csv")
chiro.xy <- read.csv("NUKspaNew2.csv")
####################################################################################################
#####################################################################################################

# Transform the data
chiro.h <- decostand (chiro, "hellinger")
chiro.xy.c <- scale(chiro.xy, center=TRUE, scale=FALSE)

####################################################################################################
####################################################################################################



######MEM analisys###### 
## 1. Construct the matrix of dbMEM variables
chiro.dbmem.tmp <- dbmem(chiro.xy, silent=FALSE)
chiro.dbmem <- as.data.frame(chiro.dbmem.tmp)
# Truncation distance used above:
(thr <- give.thresh(dist(chiro.xy)))

###################################################################################################
###################################################################################################
# 2. Test and forward selection of the environmental variables
# Forward selection of the environmental variables
chiro.env.rda <- rda(chiro.h ~., chiro.env)
(chiro.env.R2a <- RsquareAdj(chiro.env.rda)$adj.r.squared)
chiro.env.fwd <- 
  forward.sel(chiro.h, chiro.env,
              adjR2thresh = chiro.env.R2a, 
              nperm = 9999)
env.sign <- sort(chiro.env.fwd$order)
env.red <- chiro.env[ ,c(env.sign)]
colnames(env.red)

# 3. Test and forward selection of the dbMEM variables
# Run the global dbMEM analysis on the *detrended* chiro data
chiro.det.dbmem.rda <- rda(chiro.h ~., chiro.dbmem)
anova(chiro.det.dbmem.rda)
# Since the analysis is significant, compute the adjusted R2
# and run a forward selection of the dbMEM variables
(chiro.det.dbmem.R2a <- 
    RsquareAdj(chiro.det.dbmem.rda)$adj.r.squared)
(chiro.det.dbmem.fwd <- 
    forward.sel(chiro.h, 
                as.matrix(chiro.dbmem), 
                adjR2thresh = chiro.det.dbmem.R2a))
# Number of significant dbMEM
(nb.sig.dbmem <- nrow(chiro.det.dbmem.fwd))
# Identify the significant dbMEM sorted in increasing order
(dbmem.sign <- sort(chiro.det.dbmem.fwd$order))
# Write the significant dbMEM to a new object (reduced set)
dbmem.red <- chiro.dbmem[ ,c(dbmem.sign)]


## 5. Chiro - environment - trend - dbMEM variation partitioning
(chiro.varpart <- varpart(chiro.h, env.red, chiro.dbmem))
dev.new(title="Chiro - environment - dbMEM variation partitioning", width=12, height=6)
par(mfrow=c(1,2))
showvarparts(4) # To show the symbols of the fractions
plot(chiro.varpart, digits=2)
anova(
  rda(chiro.h, env.red, chiro.dbmem))
# Tests of the unique fractions [a] and [b]
# Fraction [a], pure environmental
anova(rda(chiro.h, env.red, cbind(chiro.dbmem)))
# Fraction [b], pure scale spatial
anova(rda(chiro.h, chiro.dbmem, cbind(env.red)))
###############################################################################################################################
###############################################################################################################################