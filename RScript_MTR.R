# MURCIA TWIN REGISTER 
# EATING BEHAVIOURS AND BMI - multivairate analyses 
# Moritz Herle, Jose Morosoli Garcia 

#-------------------------------------------------------------------------------
#---------------- MULTI ACE Model 
#-------------------------------------------------------------------------------			


rm(list=ls())
ls()


install.packages("psych")
install.packages ('Hmisc')
install.packages ('OpenMx')

require(OpenMx)
require(Hmisc)
require(psych)

setwd("R:/XXX")

## Reading in data 
library("foreign")
#mxOption(NULL,"Default optimizer","SLSQP")
data <- read.table("TFEQ3.csv", sep=",", header= TRUE)
class(data)
dim (data)
describe (data)
head (data)
tail (data)

# View (data)
data$BMI_T2[data$BMI_T2==-9] <- NA

describe (data)

## regress EOE, cogR and UncE by age

data$rEOE <- residuals (lm (data$EE ~ data$Age_True, na.action="na.exclude"))
data$rCogR <- residuals (lm (data$CR ~ data$Age_True, na.action="na.exclude"))
data$rUncE <- residuals (lm (data$UE ~ data$Age_True, na.action="na.exclude"))
data$rBMI_T2 <- residuals (lm (data$BMI_T2 ~ data$Age_True, na.action="na.exclude"))

hist (data$rBMI_T2)
summary (data)
data$TrBMI <- (data$rBMI_T2 + 11)
hist (data$TrBMI)

data$TrEOE <- (data$rEOE+2)
hist (data$TrEOE)

data$TrCogR <- (data$rCogR +7)
hist (data$TrCogR )

data$TrUncE <- (data$rUncE +5)*2
hist (data$TrUncE )
describe (data)

# create sub dataset with variables for analyses

sub <- data [, c ("NIP1","NIP4", "Zyg","TrEOE", "TrCogR", "TrUncE", "TrBMI")]

describe (sub)

# reshape dataset
TwinData 		<- reshape(sub, idvar = c("NIP1", "Zyg"), timevar = "NIP4", direction = "wide")
describe (TwinData)
tail (TwinData)

#View (TwinData)

Data <- TwinData [, c("Zyg", "TrEOE.A","TrEOE.B", "TrCogR.A", "TrCogR.B", "TrUncE.A" ,"TrUncE.B", "TrBMI.A",  "TrBMI.B")]

names (Data) <- c("zyg", "eoe1", "eoe2", "cogr1", "cogr2", "unce1", "unce2", "bmi1", "bmi2")
head (Data)
describe (Data)

mxOption(NULL,"Default optimizer","SLSQP")


# Select variables for analysis
nv        <- 4       # number of traits

ntv       <- nv*2    # number of total variables

nlower	<- ntv*(ntv+1)/2 	# number of free elements in a lower matrix ntv*ntv
ncor	<-(nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

Vars      <- c("eoe", "cogr", "unce", "bmi")
selVars   <- c("eoe1", "cogr1", "unce1", "bmi1", "eoe2", "cogr2", "unce2","bmi2")


mzData <- subset(Data , zyg==1, selVars)
dzData <- subset(Data , zyg==2, selVars)  


# To create Labels for Lower Triangular Matrices
aLabs <- paste("a", do.call(c, sapply(seq(1, nv), function(x){ paste(x:nv, x,sep="") })), sep="")
cLabs <- paste("c", do.call(c, sapply(seq(1, nv), function(x){ paste(x:nv, x,sep="") })), sep="")
eLabs <- paste("e", do.call(c, sapply(seq(1, nv), function(x){ paste(x:nv, x,sep="") })), sep="")


mLabs	<- paste("m",1:nv,sep="")

# To create start values 
Stmean <-colMeans(mzData[,1:4],na.rm=TRUE)



# ----------------------------------------------------------------------------------------
# 1a Specify Fully Saturated Model (Phenotypic Cholesky Decomposition)
# ----------------------------------------------------------------------------------------

LabsMMz <- c("MZm11", "MZm21", "MZm31", "MZm41", "MZm12", "MZm22","MZm32", "MZm42")
LabsMDz <- c("DZm11", "DZm21", "DZm31", "DZm41", "DZm12", "DZm22","DZm32", "DZm42")

LabsSMz <- c("MZs11", "MZs21", "MZs31", "MZs41", "MZs12", "MZs22", "MZs32","MZs42")
LabsSDz <- c("DZs11", "DZs21", "DZs31", "DZs41", "DZs12", "DZs22", "DZs32","DZs42")

LabsCorMz	<- c("MZCor12", "MZCor13", "MZCor23", "MZCor14","MZCor24","MZCor34", "MZCor15", "MZCor25", "MZCor35", "MZCor45", "MZCor16","MZCor26","MZCor36","MZCor46","MZCor56","MZCor17","MZCor27","MZCor37","MZCor47", "MZCor57","MZCor67","MZCor18","MZCor28","MZCor38","MZCor48","MZCor58","MZCor68","MZCor78")
LabsCorDz	<- c("DZCor12", "DZCor13", "DZCor23", "DZCor14","DZCor24","DZCor34", "DZCor15", "DZCor25", "DZCor35", "DZCor45", "DZCor16","DZCor26","DZCor36","DZCor46","DZCor56","DZCor17","DZCor27","DZCor37","DZCor47", "DZCor57","DZCor67","DZCor18","DZCor28","DZCor38","DZCor48","DZCor58","DZCor68","DZCor78")



# Matrix & Algebra for expected means and covariances 
MZlow		<-mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=T, values=5, name="LowMZ" )
CovMZ		<-mxAlgebra( expression=LowMZ %*% t(LowMZ), name="ExpCovMZ" )
MeanMZ	<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=Stmean,labels= LabsMMz, name="ExpMeanMZ" )

DZlow		<-mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=T, values=5, name="LowDZ" )
CovDZ		<-mxAlgebra( expression=LowDZ %*% t(LowDZ), name="ExpCovDZ" )
MeanDZ	<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=Stmean,labels= LabsMDz,  name="ExpMeanDZ" )

# Data objects for Multiple Groups
dataMZ   <- mxData( observed=mzData, type="raw" )
dataDZ   <- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ    <- mxExpectationNormal( covariance="ExpCovMZ", means="ExpMeanMZ", dimnames=selVars)
objDZ    <- mxExpectationNormal( covariance="ExpCovDZ", means="ExpMeanDZ", dimnames=selVars)

fitFunction <- mxFitFunctionML()

# Combine Groups
modelMZ	<- mxModel( MZlow, CovMZ, MeanMZ, dataMZ, objMZ, fitFunction, name="MZ")
modelDZ	<- mxModel( DZlow, CovDZ, MeanDZ, dataDZ, objDZ, fitFunction, name="DZ")
minus2ll	<- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<- mxFitFunctionAlgebra( "m2LL" )
SatModel 	<- mxModel( "Sat", modelMZ, modelDZ, minus2ll, obj)

# -------------------------------------------------------------------------------------------------------------------------------
# 1a)	RUN Saturated Model (Cholesky Decomposition)

SatFit		<- mxRun(SatModel)
(SatSum		<- summary(SatFit))

# Generate SatModel Output
round(SatFit@output$estimate,4)



# -----------------------------------------------------------------------------------------------------------------------------------
# 1b)	Specify Multivariate Saturated Model (Gaussian Decomposition) 
# 	We will use this specification to fit a constrained model
# -----------------------------------------------------------------------------------------------------------------------------------
#svM	<-c(2,4,2,11,2,4,2,11)				# these are start values
svSD	<-c(2,2,1,5,2,2,1,5)
svMZ	<-c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5,0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
svDZ	<-c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5,0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

LabsMMz <- c("MZm11", "MZm21", "MZm31", "MZm41", "MZm12", "MZm22","MZm32", "MZm42")
LabsMDz <- c("DZm11", "DZm21", "DZm31", "DZm41", "DZm12", "DZm22","DZm32", "DZm42")

LabsSMz <- c("MZs11", "MZs21", "MZs31", "MZs41", "MZs12", "MZs22", "MZs32","MZs42")
LabsSDz <- c("DZs11", "DZs21", "DZs31", "DZs41", "DZs12", "DZs22", "DZs32","DZs42")

LabsCorMz	<- c("MZCor12", "MZCor13", "MZCor23", "MZCor14","MZCor24","MZCor34", "MZCor15", "MZCor25", "MZCor35", "MZCor45", "MZCor16","MZCor26","MZCor36","MZCor46","MZCor56","MZCor17","MZCor27","MZCor37","MZCor47", "MZCor57","MZCor67","MZCor18","MZCor28","MZCor38","MZCor48","MZCor58","MZCor68","MZCor78")
LabsCorDz	<- c("DZCor12", "DZCor13", "DZCor23", "DZCor14","DZCor24","DZCor34", "DZCor15", "DZCor25", "DZCor35", "DZCor45", "DZCor16","DZCor26","DZCor36","DZCor46","DZCor56","DZCor17","DZCor27","DZCor37","DZCor47", "DZCor57","DZCor67","DZCor18","DZCor28","DZCor38","DZCor48","DZCor58","DZCor68","DZCor78")


# Matrix & Algebra for expected means and covariances 
MeanMZ	<-mxMatrix( "Full", 1, ntv, free=T, values=Stmean, labels= LabsMMz, name="ExpMeanMZ" )
MZsd		<-mxMatrix( "Diag", ntv, ntv, free=T, values=svSD, labels =LabsSMz, name="sdMZ" )
Cormz		<-mxMatrix( "Stand", ntv, ntv, free=T, values=0.5, labels=LabsCorMz, name="MZCor") 
CovMZ		<-mxAlgebra( expression=sdMZ %*% MZCor %*% t(sdMZ), name="ExpCovMZ" )

MeanDZ	<-mxMatrix( "Full", 1, ntv, free=T, values=Stmean, labels=LabsMDz, name="ExpMeanDZ" )
DZsd		<-mxMatrix( "Diag", ntv, ntv, free=T, values=svSD, labels=LabsSDz, name="sdDZ" )
Cordz		<-mxMatrix( "Stand", ntv, ntv, free=T, values=0.5, labels=LabsCorDz,  name="DZCor")
CovDZ		<-mxAlgebra( expression=sdDZ %*% DZCor %*% t(sdDZ), name="ExpCovDZ" ) 

# Data objects for Multiple Groups
dataMZ   <- mxData( observed=mzData, type="raw" )
dataDZ   <- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ    <- mxExpectationNormal( covariance="ExpCovMZ", means="ExpMeanMZ", dimnames=selVars)
objDZ    <- mxExpectationNormal( covariance="ExpCovDZ", means="ExpMeanDZ", dimnames=selVars)

fitFunction <- mxFitFunctionML()

# Combine Groups
modelMZ	<- mxModel( MeanMZ, MZsd, Cormz, CovMZ, dataMZ, objMZ, fitFunction, name="MZ")
modelDZ	<- mxModel( MeanDZ, DZsd, Cordz, CovDZ, dataDZ, objDZ, fitFunction, name="DZ")
minus2ll	<- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<- mxFitFunctionAlgebra( "m2LL" )
Conf1		<- mxCI (c ('MZ.MZCor[5,1]', 'MZ.MZCor[6,2]', 'MZ.MZCor[7,3]', 'MZ.MZCor[8,4]') )
Conf2		<- mxCI (c ('DZ.DZCor[5,1]', 'DZ.DZCor[6,2]', 'DZ.DZCor[7,3]', 'DZ.DZCor[8,4]') )
SatGModel 	<- mxModel( "SatG", modelMZ, modelDZ, minus2ll, obj, Conf1, Conf2)

# -------------------------------------------------------------------------------------------------------------------------------
# 1b)	RUN Saturated Model (Gaussian Decomposition)

SatGFit	<- mxRun(SatGModel, intervals=T)
(SatGSum	<- summary(SatGFit))

# Generate SatGModel Output
round(SatGFit@output$estimate,4)

mxEval(MZ.ExpMeanMZ, SatGFit)
mxEval(MZ.ExpCovMZ, SatGFit)
mxEval(DZ.ExpMeanDZ, SatGFit)
mxEval(DZ.ExpCovDZ, SatGFit)


# --------------------------------------------------------------------------------------------
# 1b)	Specify & Run Sub1: Constrained Bivariate Model (Gaussian Decomposition):
#	Equal Means & Variances across Twin Order 
#	One overall set of Within-person cross-trait correlations
#	Symmetric xtwin-xtrait cor matrices in MZ and DZ group 
# --------------------------------------------------------------------------------------------

# To manipulate  the parameters in matrices we wish to change in the full model, we use the
# 'omxSetParameters' function with the original 'labels' and 'newlabels' to indicate the changes
# i.e. specifying the same label effectively constraints the parameters to be the same 

Sub1bModel	<- mxModel(SatGModel, name="Sub1" )

# means

Sub1bModel	<- omxSetParameters( Sub1bModel, free=T, values=Stmean, labels=c("MZm11", "MZm21", "MZm31", "MZm41", "MZm12", "MZm22","MZm32", "MZm42"),
                                newlabels=c("MZm1", "MZm2", "MZm3", "MZm4","MZm1", "MZm2", "MZm3", "MZm4"))

Sub1bModel	<- omxSetParameters( Sub1bModel, free=T, values=Stmean, labels=c("DZm11", "DZm21", "DZm31", "DZm41", "DZm12", "DZm22","DZm32", "DZm42"), 
                                newlabels=c("DZm1", "DZm2", "DZm3", "DZm4","DZm1", "DZm2", "DZm3", "DZm4"))

# SDs

Sub1bModel	<- omxSetParameters( Sub1bModel, free=T, values=svSD, labels=c("MZs11", "MZs21", "MZs31", "MZs41", "MZs12", "MZs22", "MZs32","MZs42"),
                                newlabels=c("MZs11", "MZs21", "MZs31", "MZs41","MZs11", "MZs21", "MZs31", "MZs41"))

Sub1bModel	<- omxSetParameters( Sub1bModel, free=T, values=svSD, labels=c("DZs11", "DZs21", "DZs31", "DZs41", "DZs12", "DZs22", "DZs32","DZs42"),
                                newlabels=c("DZs11", "DZs21", "DZs31", "DZs41","DZs11", "DZs21", "DZs31", "DZs41"))


# ------------------------------------------------------------------------------
# 2) RUN Sub1Model
Sub1bFit     <- mxRun( Sub1bModel, intervals=T )
(Sub1bSum     <- summary( Sub1bFit ))


mxCompare(SatGFit, Sub1bFit)


# --------------------------------------------------------------------------------------------
# 2)	Specify & Run Sub1: Constrained Bivariate Model (Gaussian Decomposition):
#	Equal Means & Variances across Twin Order & zyg group
#	One overall set of Within-person cross-trait correlations
#	Symmetric xtwin-xtrait cor matrices in MZ and DZ group 
# --------------------------------------------------------------------------------------------

# To manipulate  the parameters in matrices we wish to change in the full model, we use the
# 'omxSetParameters' function with the original 'labels' and 'newlabels' to indicate the changes
# i.e. specifying the same label effectively constraints the parameters to be the same 

Sub1Model	<- mxModel(SatGModel, name="Sub1" )

# means

Sub1Model	<- omxSetParameters( Sub1Model, free=T, values=Stmean, labels=c("MZm11", "MZm21", "MZm31", "MZm41", "MZm12", "MZm22","MZm32", "MZm42"),
                               newlabels=c("m1", "m2", "m3", "m4","m1", "m2", "m3", "m4"))

Sub1Model	<- omxSetParameters( Sub1Model, free=T, values=Stmean, labels=c("DZm11", "DZm21", "DZm31", "DZm41", "DZm12", "DZm22","DZm32", "DZm42"), 
                               newlabels=c("m1", "m2", "m3", "m4","m1", "m2", "m3", "m4"))

# SDs

Sub1Model	<- omxSetParameters( Sub1Model, free=T, values=svSD, labels=c("MZs11", "MZs21", "MZs31", "MZs41", "MZs12", "MZs22", "MZs32","MZs42"),
                               newlabels=c("s11", "s21", "s31", "s41","s11", "s21", "s31", "s41"))

Sub1Model	<- omxSetParameters( Sub1Model, free=T, values=svSD, labels=c("DZs11", "DZs21", "DZs31", "DZs41", "DZs12", "DZs22", "DZs32","DZs42"),
                               newlabels=c("s11", "s21", "s31", "s41","s11", "s21", "s31", "s41"))

# cor 

Sub1Model	<- omxSetParameters( Sub1Model, free=T, values=svMZ, labels=c("MZCor12", "MZCor13", "MZCor23", "MZCor14","MZCor24","MZCor34", "MZCor15", "MZCor25", "MZCor35", "MZCor45", "MZCor16","MZCor26","MZCor36","MZCor46","MZCor56","MZCor17","MZCor27","MZCor37","MZCor47", "MZCor57","MZCor67","MZCor18","MZCor28","MZCor38","MZCor48","MZCor58","MZCor68","MZCor78"),
                               newlabels=c("rph12","rph13", "rph23", "rph14", "rph24", "rph34", "MzEOE", "MZxtxt_12", "MZxtxt_13", "MZxtxt_14", "MZxtxt_12", "MzCR","MZxtxt_32", "MZxtxt_42","rph12","MZxtxt_13","MZxtxt_32","MzUn","MZxtxt_43","rph13","rph23", "MZxtxt_14", "MZxtxt_42","MZxtxt_43", "MzBmi", "rph14","rph24","rph34"))

Sub1Model	<- omxSetParameters( Sub1Model, free=T, values=svDZ, labels=c("DZCor12", "DZCor13", "DZCor23", "DZCor14","DZCor24","DZCor34", "DZCor15", "DZCor25", "DZCor35", "DZCor45", "DZCor16","DZCor26","DZCor36","DZCor46","DZCor56","DZCor17","DZCor27","DZCor37","DZCor47", "DZCor57","DZCor67","DZCor18","DZCor28","DZCor38","DZCor48","DZCor58","DZCor68","DZCor78"),
                               newlabels=c("rph12","rph13", "rph23", "rph14", "rph24", "rph34", "DzEOE", "DZxtxt_12", "DZxtxt_13", "DZxtxt_14", "DZxtxt_12", "DzCR","DZxtxt_32", "DZxtxt_42","rph12","DZxtxt_13","DZxtxt_32","DzUn","DZxtxt_43","rph13","rph23", "DZxtxt_14", "DZxtxt_42","DZxtxt_43", "DzBmi", "rph14","rph24","rph34"))

# ------------------------------------------------------------------------------
# 2) RUN Sub1Model
Sub1Fit      <- mxRun( Sub1Model, intervals=T )
(Sub1Sum     <- summary( Sub1Fit ))


mxCompare(SatGFit, Sub1bFit)
mxCompare(Sub1bFit, Sub1Fit)

mxCompare(SatFit, Sub1Fit)



#########################################
#### specify model 	ACE	
#########################################


# Matrix & Algebra for expected means vector
# Matrices declared to store a, c, and e Path Coefficients
pathA	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=.4, labels=aLabs, name="a" )
pathC	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=.1, labels=cLabs, name="c" )
pathE	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=.2, labels=eLabs, name="e" )

# Matrices generated to hold A, C, and E computed Variance Components
covA	<- mxAlgebra( expression=a %*% t(a), name="A" )
covC	<- mxAlgebra( expression=c %*% t(c), name="C" ) 
covE	<- mxAlgebra( expression=e %*% t(e), name="E" )
covP	<- mxAlgebra( expression=A+C+E, name="V" )
StA	<- mxAlgebra( expression=A/V, name="h2" )
StC	<- mxAlgebra( expression=C/V, name="c2" )
StE	<- mxAlgebra( expression=E/V, name="e2" )

# Algebra to compute Correlations
matI	<- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
Rph	<- mxAlgebra( expression=solve(sqrt(I*V)) %&% V , name="Phcor")
Rg	<- mxAlgebra( expression=solve(sqrt(I*A)) %&% A , name="Acor")
Rc	<- mxAlgebra( expression=solve(sqrt(I*C)) %&% C , name="Ccor")
Re	<- mxAlgebra( expression=solve(sqrt(I*E)) %&% E , name="Ecor")

# Algebra for expected Mean and Variance/Covariance Matrices in MZ & DZ twins
Mean	<- mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=Stmean, labels=c(mLabs,mLabs), name="ExpMean" )

covMZ	<- mxAlgebra( expression= rbind( cbind(A+C+E , A+C),
                                       cbind(A+C , A+C+E)),		name="ExpCovMZ" )
covDZ	<- mxAlgebra( expression= rbind( cbind(A+C+E       , 0.5%x%A+C),
                                       cbind(0.5%x%A+C , A+C+E)),	name="ExpCovDZ" )

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ	<- mxExpectationNormal( covariance="ExpCovMZ", means="ExpMean", dimnames=selVars )
objDZ	<- mxExpectationNormal( covariance="ExpCovDZ", means="ExpMean", dimnames=selVars )

fitFunction <- mxFitFunctionML()

# Combine Groups
pars	<- list( pathA, pathC, pathE, covA, covC, covE, covP, StA, StC, StE, matI, Rph, Rg, Rc, Re) 
modelMZ	<- mxModel( pars, covMZ, Mean, dataMZ, objMZ, fitFunction, name="MZ" )
modelDZ	<- mxModel( pars, covDZ, Mean, dataDZ, objDZ, fitFunction, name="DZ" )
minus2ll<- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj	<- mxFitFunctionAlgebra( "m2LL" )
conf1	<- mxCI (c ('MZ.h2[1,1]', 'MZ.h2[2,2]', 'MZ.h2[3,3]', 'MZ.h2[4,4]' ) )		# h2
conf2	<- mxCI (c ('MZ.c2[1,1]', 'MZ.c2[2,2]', 'MZ.c2[3,3]', 'MZ.c2[4,4]') )		# c2
conf3	<- mxCI (c ('MZ.e2[1,1]', 'MZ.e2[2,2]', 'MZ.e2[3,3]', 'MZ.e2[4,4]') )		# e2
conf4	<- mxCI (c ('MZ.Acor[2,1]','MZ.Acor[3,1]','MZ.Acor[4,1]','MZ.Acor[3,2]','MZ.Acor[4,2]','MZ.Acor[4,3]' )) #Rg
conf5	<- mxCI (c ('MZ.Ccor[2,1]','MZ.Ccor[3,1]','MZ.Ccor[4,1]','MZ.Ccor[3,2]','MZ.Ccor[4,2]','MZ.Ccor[4,3]' )) #Rc			
conf6	<- mxCI (c ('MZ.Ecor[2,1]','MZ.Ecor[3,1]','MZ.Ecor[4,1]','MZ.Ecor[3,2]','MZ.Ecor[4,2]','MZ.Ecor[4,3]' )) #Re
conf7	<- mxCI (c ('Phcor[2,1]', 'Phcor[3,1]', 'Phcor[4,1]', 'Phcor[3,2]', 'Phcor[4,2]', 'Phcor[4,3]')) #Rph
CholAceModel<- mxModel( "CholACE", pars, modelMZ, modelDZ, minus2ll, obj, conf1, conf2, conf3, conf4, conf5, conf6, conf7)

# -------------------------------------------------
# 	RUN full ACE MODEL
CholAceFit	<- mxTryHard(CholAceModel, intervals=T)
(CholAceSum	<- summary(CholAceFit))




mxCompare(SatFit , CholAceFit)





######################################################################################
# Explorative: Fit constrained model based on univariate estimates
# BMI: AE model 
# CR: CE model 
# EE, UE: ACE model 
# order of vars: EE, CR, UE, BMI




# Matrix & Algebra for expected means vector
# Matrices declared to store a, c, and e Path Coefficients
pathA	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, 
                   free=c(T, T, T, T, F, T, T, T, T, T), values=c(0.4, 0.4  , 0.4, 0.4, 0, 0.4, 0.4, 0.4, 0.4, 0.4), free=F, values=0), labels=c("a11", "a21", "a31", "a41", "a22", "a32", "a42", "a33", "a43", "a44"), name="a" )

pathC	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, 
                   free=c(T, T, T, T, T, T, T, T, T,F), 
                   values=c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, 0),free=F, values=0) ,  
                   labels=c("c11", "c21", "c31", "c41", "c22", "c32", "c42", "c33", "c43", "c44"), name="c" )

pathE	<- mxMatrix( type="Lower", nrow=nv, ncol=nv, 
                   free=c(T, T, T, T, T, T, T, T, T, T), 
                   values=c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,.02), free=F, values=0) ,
                  labels=c("e11", "e21", "e31", "e41", "e22", "e32", "e42", "e33", "e43", "e44"), name="e" )

# Matrices generated to hold A, C, and E computed Variance Components
covA	<- mxAlgebra( expression=a %*% t(a), name="A" )
covC	<- mxAlgebra( expression=c %*% t(c), name="C" ) 
covE	<- mxAlgebra( expression=e %*% t(e), name="E" )
covP	<- mxAlgebra( expression=A+C+E, name="V" )
StA	<- mxAlgebra( expression=A/V, name="h2" )
StC	<- mxAlgebra( expression=C/V, name="c2" ) 
StE	<- mxAlgebra( expression=E/V, name="e2" )

# Algebra to compute Correlations
matI	<- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
Rph	<- mxAlgebra( expression=solve(sqrt(I*V)) %&% V , name="Phcor")
Rg	<- mxAlgebra( expression=solve(sqrt(I*A)) %&% A , name="Acor")
Rc	<- mxAlgebra( expression=solve(sqrt(I*C)) %&% C , name="Ccor")
Re	<- mxAlgebra( expression=solve(sqrt(I*E)) %&% E , name="Ecor")

# Algebra for expected Mean and Variance/Covariance Matrices in MZ & DZ twins
Mean	<- mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=Stmean, labels=c(mLabs,mLabs), name="ExpMean" )

covMZ	<- mxAlgebra( expression= rbind( cbind(A+C+E , A+C),
                                       cbind(A+C , A+C+E)),		name="ExpCovMZ" )
covDZ	<- mxAlgebra( expression= rbind( cbind(A+C+E       , 0.5%x%A+C),
                                       cbind(0.5%x%A+C , A+C+E)),	name="ExpCovDZ" )

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ	<- mxExpectationNormal( covariance="ExpCovMZ", means="ExpMean", dimnames=selVars )
objDZ	<- mxExpectationNormal( covariance="ExpCovDZ", means="ExpMean", dimnames=selVars )

fitFunction <- mxFitFunctionML()

# Combine Groups
pars	<- list( pathA, pathC, pathE, covA, covC, covE, covP, StA, StC, StE, matI, Rph, Rg, Rc, Re) 
modelMZ	<- mxModel( pars, covMZ, Mean, dataMZ, objMZ, fitFunction, name="MZ" )
modelDZ	<- mxModel( pars, covDZ, Mean, dataDZ, objDZ, fitFunction, name="DZ" )
minus2ll<- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj	<- mxFitFunctionAlgebra( "m2LL" )
conf1	<- mxCI (c ('MZ.h2[1,1]', 'MZ.h2[3,3]', 'MZ.h2[4,4]' ) )		# h2
conf2	<- mxCI (c ('MZ.c2[1,1]', 'MZ.c2[2,2]', 'MZ.c2[3,3]' ) )		# c2
conf3	<- mxCI (c ('MZ.e2[1,1]', 'MZ.e2[2,2]', 'MZ.e2[3,3]', 'MZ.e2[4,4]') )		# e2
conf4	<- mxCI (c ('MZ.Acor[3,1]','MZ.Acor[4,1]','MZ.Acor[4,3]' )) #Rg
conf5	<- mxCI (c ('MZ.Ccor[2,1]','MZ.Ccor[3,1]','MZ.Ccor[3,2]' )) #Rc			
conf6	<- mxCI (c ('MZ.Ecor[2,1]','MZ.Ecor[3,1]','MZ.Ecor[4,1]','MZ.Ecor[3,2]','MZ.Ecor[4,2]','MZ.Ecor[4,3]' )) #Re
conf7	<- mxCI (c ('Phcor[2,1]', 'Phcor[3,1]', 'Phcor[4,1]', 'Phcor[3,2]', 'Phcor[4,2]', 'Phcor[4,3]')) #Rph
ace2model<- mxModel( "ace2", pars, modelMZ, modelDZ, minus2ll, obj, conf1, conf2, conf3, conf4, conf5, conf6, conf7)

# -------------------------------------------------
# 	RUN full ACE MODEL
ace2	<- mxTryHard(ace2model, intervals=T)
(ace2Sum	<- summary(ace2))



mxCompare(CholAceFit , ace2)













