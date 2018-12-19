library(simplerspec)
library(dplyr)
library(matrixStats)
library(bimixt)
library(car)
library(pls)
library(ranger)
library(resemble)
library(Cubist)
library(Rcpp)
library(doParallel)
library(hexbin)
library(RColorBrewer)

source("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/outlier/functions_calibTransfer.R")


###step 1: extract MIR spectra(opus format) from KSSL database for all soil properties using simplerspec
### see instruction for updating the simplerspec package -- since recent version of R requires to update Rcpp files

#list directory and files containing MIR in opus format
dirs <- list.dirs("/mnt/WHRC/sdangal/KSSL/MIRJune2018/MIR_Spectral_Library", full.names=TRUE)
all.files <- list.files(dirs, pattern= "*.0", recursive=TRUE,full.names=TRUE)

# read opus format files as list and gather spc
spc_list <- read_opus_univ(fnames = all.files, extract = c("spc"))
soilspec_tbl <- spc_list %>%
  # Gather list of spectra data into tibble data frame
  gather_spc()

#process to output spectra in desired format, average four replicates and save all spectra as a dataframe to build models
spc <- soilspec_tbl$spc
spc.df <- as.data.frame(matrix(unlist(spc), nrow=length(spc), byrow=T))
colnames(spc.df) <- colnames(spc[[1]])
spc.df <- data.frame(sample_id = soilspec_tbl$sample_id, spc.df)
spc.df$sample_id <- str_sub(spc.df$sample_id,1,str_length(spc.df$sample_id)-4)
spc.avg.df <- aggregate(.~sample_id, data = spc.df, FUN=mean, na.rm=TRUE)
spc.mat <- spc.avg.df[,1:2]
spc.mat$spc <- as.matrix(spc.avg.df[,2:ncol(spc.avg.df)])
spc.mat <- spc.mat[,-2]
spc.mat$sample_id <- as.numeric(spc.mat$sample_id)
save(spc.mat, file = "spc.mat.RData")

#### step 2: combine soil properties with MIR spectra to build models
##load csv files containing soil properties data with sample id
soilp <- read.csv("/whrc/sdangal/KSSL/Data/soilpropwithMIR.csv")
spc.mat$smp_id <- spc.mat$sample_id
merged.data <- merge(soilp, test, by = "smp_id")
save(merged.data, file = "merged.data.RData")


#### step 3: truncate spectra to 6000 - 600 cm-1 and remove CO2 sensitive region (2389-2268 cm-1)
#remove CO2 sensitive bands 2389-2268 
col.names <- colnames(merged.data$spc)
col.names <- as.numeric(substring(col.names,2))
min.index <- which(col.names <= 2389)[1]
max.index <- which(col.names <= 2268)[1]
merged.data$spc <- merged.data$spc[,-c(min.index:max.index)] 
##truncate the spectra from 6000 to 600 cm-1
merged.data$spc <- merged.data$spc[,-c(1:778)]

### step 4: perform baseline transformation and divide data by soil properties
merged.data$spc <- base_offset(merged.data$spc)
full.oc <- merged.data[!is.na(merged.data$OC),] ## repeat this step for other soil properties -- Al, Ca, CO3, pH, Fe, Ca, BD, OCD, Clay.
save(full.oc, file = "full.oc.RData")

### step 5: remove outliers by selecting using a standard deviaton threshold
fit.oc <- plsr(sqrt(OC)~spc, ncomp=20, data = full.oc, valid="CV", segments = 50)
pred.coc <- predict(fit.oc, newdata = full.oc$spc,ncomp=20)
sd.outlier <- optimum_sd_outlier(pred, resp, seq(0.1,3, by =0.02))
row.index <- outlier(pred, resp, sd.outlier[1])
full.oc <- full.oc[row.index,]


### step 6: perform kennard stone to separate data into 80% calibration and 20% validation sets
ks.OC <- kenStone(X = full.oc$spc, k = as.integer(0.8*nrow(full.oc)), metric = "mahal", pc = 10) ## repeat this step for other soil properties -- Al, Ca, CO3, pH, Fe, Ca, BD, OCD, Clay.
calib.oc <- full.oc[ks.OC$model, ]
valid.oc <- full.oc[ks.OC$test, ]
save(calib.oc, file = "calib.oc.RData")
save(valid.oc, file = "valid.oc.RData")


###step 7: build models using box-cox, log, squareroot and untransformed soil properties to test model performance with and without normal distribution
### only pH and OC carbon models were built for testing model performance using different transformations

sub.calib.oc$bc.oc <- for.trans(sub.calib.oc$oc)$var
lam <- for.trans(sub.calib.oc$OC)$lam
lam  ##need to keep track of lam when doing back transformation
sub.calib.ph$bc.ph <- for.trans(sub.calib.ph$ph)$var
lam <- for.trans(sub.calib.ph$ph)$lam
lam

#PLSR models  -- only untransformed OC shown
fit.plsr.untrans.oc <- plsr(oc~spc, ncomp=20, data = sub.calib.oc, valid="CV")
ncomp.onesigma <- selectNcomp(fit.plsr.untrans.oc, method = "onesigma", plot = TRUE,ylim = c(0, .2))
ncomp.onesigma
val.pred <- predict(fit.plsr.untrans.oc, newdata = sub.valid.oc$spc, ncomp=ncomp.onesigma)

##Random Forest models -- only untransformed OC shown
 X1 <- data.frame(sub.calib.oc$spc)
Y1 <- sub.calib.oc$oc
X1 <- X1[!is.na(Y1), ]
Y1 <- Y1[!is.na(Y1)]
tmp.calib <- cbind(Y1,X1)
fit.rf.untrans.oc <- ranger(Y1~., data = tmp.calib, quantreg=TRUE, keep.inbag=TRUE, num.trees=150)
val.pred <- predict(fit.rf.untrans.oc, data = sub.valid.oc$spc, type = "se")

## Spectrum based learner model  -- only untransformed OC shown
Xu <- sub.valid.oc$spc
Yu <- sub.valid.oc$oc 
Yr <- sub.calib.oc$oc
Xr <- sub.calib.oc$spc
Xu <- Xu[!is.na(Yu),]
Yu <- Yu[!is.na(Yu)]
Xr <- Xr[!is.na(Yr),]
Yr <- Yr[!is.na(Yr)]
ctrl <- mblControl(sm = 'pc', pcSelection = list('opc', 50),
                   valMethod = 'loc_crossval',center=TRUE,scale=FALSE,allowParallel=FALSE)
sbl.untrans.oc <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
                      mblCtrl = ctrl,
                      dissUsage = 'none',
                      k = seq(40, 100, by = 20),
                      method = 'pls', pls.c = 6)
obs <- sbl.untrans.oc$results$Nearest_neighbours_40$yu.obs
pred <- sbl.untrans.oc$results$Nearest_neighbours_40$pred
val.pred <- data.frame(obs, pred)

##Cubist Model -- only untransformed OC shown
resp <- sub.calib.oc$oc
sub.calib.oc <- sub.calib.oc[!is.na(resp),]
resp <- resp[!is.na(resp)]
cub.untrans.oc <- cubist(x=sub.calib.oc$spc, y = resp)
valid.pred <- predict(cub.untrans.oc, sub.valid.oc$spc)
val.pred <- data.frame(sub.valid.oc$oc, valid.pred)


###step 8: start building models using square root transformed soil properties
### only OC shown -- need to change the variable name to built models for all other soil properties

#load data
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/full_data_trun6000/sub.calib.oc.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/full_data_trun6000/sub.valid.oc.RData")
sub.calib.oc$oc <- sub.calib.oc$OC
sub.valid.oc$oc <- sub.valid.oc$OC

## PLSR model
fit.oc <- plsr(sqrt(oc)~spc, ncomp=20, data = sub.calib.oc, valid="CV")
save(fit.oc, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/process_pls/fit.oc.RData")
ncomp.onesigma <- selectNcomp(fit.oc, method = "onesigma", plot = TRUE,ylim = c(0, .2))
ncomp.onesigma
cal.pred <- predict(fit.oc, newdata = sub.calib.oc$spc, ncomp=ncomp.onesigma)
val.pred <- predict(fit.oc, newdata = sub.valid.oc$spc, ncomp=ncomp.onesigma)
pls.calib.pred.oc <- data.frame(sub.calib.oc$oc, cal.pred^2)
pls.valid.pred.oc <- data.frame(sub.valid.oc$oc, val.pred^2)
names(pls.calib.pred.oc) <- c("obs", "pred")
names(pls.valid.pred.oc) <- c("obs", "pred")
write.csv(pls.calib.pred.oc, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/process_pls/calib.pred.oc.csv")
write.csv(pls.valid.pred.oc, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/process_pls/valid.pred.oc.csv")

## Spectrum based learner model
Xu <- sub.valid.oc$spc
Yu <- sqrt(sub.valid.oc$oc) 
Yr <- sqrt(sub.calib.oc$oc)
Xr <- sub.calib.oc$spc
Xu <- Xu[!is.na(Yu),]
Yu <- Yu[!is.na(Yu)]
Xr <- Xr[!is.na(Yr),]
Yr <- Yr[!is.na(Yr)]
ctrl <- mblControl(sm = 'pc', pcSelection = list('opc', 50),
                   valMethod = 'loc_crossval',center=TRUE,scale=FALSE,allowParallel=FALSE)
sbl.sqrt.oc <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
                   mblCtrl = ctrl,
                   dissUsage = 'none',
                   k = seq(40, 100, by = 20),
                   method = 'pls', pls.c = 6)
save(sbl.sqrt.oc, file = "/home/sdangal/test/localreg/sbl.sqrt.oc.RData")

## Random Forest Model
X1 <- data.frame(sub.calib.oc$spc)
Y1 <- sqrt(sub.calib.oc$oc)
X1 <- X1[!is.na(Y1), ]
Y1 <- Y1[!is.na(Y1)]
tmp.calib <- cbind(Y1,X1)
plsr.sqrt.oc <- ranger(Y1~., data = tmp.calib, quantreg=TRUE, keep.inbag=TRUE, num.trees=150)
save(plsr.sqrt.oc, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/process_rf/plsr.sqrt.oc.RData")

pred <- predict(plsr.sqrt.oc, data = sub.calib.oc$spc, type = "se")
pred.new <- predict(plsr.sqrt.oc, data = sub.valid.oc$spc, type = "se")
calib.pred.oc <- data.frame(sub.calib.oc$oc,pred$predictions,pred$se) ##untransformed predictions
valid.pred.oc <- data.frame(sub.valid.oc$oc,pred.new$predictions,pred.new$se)
names(calib.pred.oc) <- c("obs", "pred","se")
names(valid.pred.oc) <- c("obs", "pred", "se")
write.csv(calib.pred.oc, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/process_rf/calib.pred.oc.csv")
write.csv(valid.pred.oc, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/process_rf/valid.pred.oc.csv")

## Cubist Model
resp <- sqrt(sub.calib.oc$oc)
sub.calib.oc <- sub.calib.oc[!is.na(resp),]
resp <- resp[!is.na(resp)]
cub.sqrt.oc <- cubist(x=sub.calib.oc$spc, y = resp)
save(cub.sqrt.oc, file = "/home/sdangal/test/localreg/cub.sqrt.oc.RData")
calib.pred <- predict(cub.sqrt.oc, sub.calib.oc$spc)
valid.pred <- predict(cub.sqrt.oc, sub.valid.oc$spc)
calib.pred.oc <- data.frame(resp^2, calib.pred^2)
valid.pred.oc <- data.frame(sub.valid.oc$oc, valid.pred^2)
write.csv(calib.pred.oc, file = "/home/sdangal/test/localreg/calib.pred.oc.csv")
write.csv(valid.pred.oc, file = "/home/sdangal/test/localreg/valid.pred.oc.csv")

## step 9: check the best model -- create files list of all model output in csv format to get quick summary
##set working directory to each model (PLSR, Cubist,RF, SBL) sub-directory -- only PLSR model shown here
setwd("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/process_pls")
#calib outputs
out.files <- list.files(getwd(), pattern = "calib", full.names=TRUE)
test <- lapply(out.files, function(x) read.csv(x))
var.names <- c("al", "bd", "ca", "cec", "clay", "co3", "fe", "OC", "ocden", "ph")
names(test) <- var.names
for(i in 1:length(out.files)){
  getsummary(test[[i]]$obs, test[[i]]$pred)
}
#valid outputs
out.files <- list.files(getwd(), pattern = "valid", full.names=TRUE)
test <- lapply(out.files, function(x) read.csv(x))
var.names <- c("al", "bd", "ca", "cec", "clay", "co3", "fe", "OC", "ocden", "ph")
names(test) <- var.names
for(i in 1:length(out.files)){
  getsummary(test[[i]]$obs, test[[i]]$pred)
}


## step 10: create figure 2 used in manuscript showing model performance 
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
#set wd to place where all validation output files are stored
out.files <- list.files(getwd(), pattern = "valid", full.names=TRUE)
nfile <- c("al", "bd", "ca", "cec", "clay", "co3", "fe", "oc", "ocden", "ph")
test <- lapply(out.files, function(x) read.csv(x))
var.names <- c("al", "bd", "ca", "cec", "clay", "co3", "fe", "OC", "ocden", "ph")
names(test) <- var.names

##repeat process below by changing the index from 2:10 to plot other variable performance
tiff(file = "plsr.al.tiff", width = 5400, height = 4000, units = "px", res = 800) 
hexbinplot(test[[1]]$obs~test[[1]]$pred, colramp=rf, main="", 
           ylab=paste("Observed Al (%wt)"), asp=1, 
           xlab=paste("Predicted Al (%wt)"), lwd=1, 
           lcex=8, inner=.4, cex.labels=1, xbins=50, 
           xlim = c(0,4.5), ylim = c(0,4.5),
           colorcut=c(0,0.005,0.01,0.03,0.07,0.15,0.25,0.5,1),panel=pfun.lm)
dev.off() 

## step 11: calculate u-deviation based on PLSR model
## Udeviation based on local model and random forest model are included in the model output
## set wd to directory where plsr model outputs are available
setwd("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/process_pls")
files <- list.files(getwd(), pattern = "\\.RData", full.names=TRUE)

##read calib data
cal.data <- list.files("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/full_data_trun6000", pattern ="calib", full.names=TRUE)

##read calib scores
temp <- lapply(files, function(x) mget(load(x)))
cal <- lapply(cal.data, function(x) mget(load(x)))
x.cal.scores <- lapply(1:length(temp), function(x){
  temp[[x]][[1]]$scores
})

##read loadings data
loadings <- lapply(1:length(temp), function(x){
  temp[[x]][[1]]$loadings
})

#read cal spc
x.cal.spc <- lapply(1:length(temp), function(x){
  cal[[x]][[1]]$spc
})

##Read valid data
val.data <- list.files("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/full_data_trun6000", pattern ="valid", full.names=TRUE)

#Read x.val.scores
val <- lapply(val.data, function(x) mget(load(x)))
x.val.scores <- lapply(1:length(temp), function(x){
  predict(temp[[x]][[1]], newdata = val[[x]][[1]]$spc, type = "scores")
})

#Read y.val.pred
y.val.pred <- lapply(1:length(temp), function(x){
  predict(temp[[x]][[1]], newdata = val[[x]][[1]]$spc)
})

## read valid spc
x.val.spc <- lapply(1:length(temp), function(x){
  val[[x]][[1]]$spc
}) 

##define pred variable names
#val.data
pred.var.names <- c("al", "bd", "ca", "cec", "clay", "co3", "fe", "OC", "ocden", "ph")

#get observed values of validation samples
val.obs <- lapply(1:length(temp), function(x){
  sqrt(val[[x]][[1]][pred.var.names[x]])
})

##this step is critical to ensure that there are no NA's values in the observations
val.pred.revised <- lapply(1:length(temp), function(x){
  y.val.pred[[x]][!is.na(val.obs[[x]]),,]
})

val.obs.revised <- lapply(1:length(temp), function(x){
  val.obs[[x]][!is.na(val.obs[[x]])]
})

#### start estimating prediction interval using PLSR model
#1 get XvalSampleVar and XValSampleTot
ResXValSamp <- lapply(1:length(temp), function(x){
  getResXValSamp(x.val.spc[[x]], x.cal.spc[[x]], x.val.scores[[x]], loadings[[x]])
})

ResXValTot <- lapply(1:length(temp), function(x){
  colMeans(ResXValSamp[[x]])
})

#2 Get ResYValVar -- residual variance for response variable in validation sets
#need observed and predicted values to calculate ResYValVar
ResYValVar <- lapply(1:length(temp), function(x){
  getResYValVar(val.obs.revised[[x]], val.pred.revised[[x]])
})


#3 Leverage
Hi <- lapply(1:length(temp), function(x){
  getLeverage(x.cal.scores[[x]], x.val.scores[[x]])
})

#4 Number of calibration object
nobj <- lapply(1:length(temp), function(x){
  dim(x.cal.spc[[x]])[1]
})

#5 compute Y-dev
ydev <- lapply(1:length(temp), function(x){
  getYdev(ResYValVar[[x]], ResXValSamp[[x]], ResXValTot[[x]], Hi[[x]], nobj[[x]])
})

names(ydev) <- pred.var.names
for(i in 1:length(ydev)){
  write.csv(ydev[i], file = paste0(names(ydev[i]), ".udev.csv"))
}


##Step 12: flagging samples that are untrustworthy
## see fratio_plsr, fratio_sbl and reld_rf.R scripts to assess the trustworthines of new prediction
