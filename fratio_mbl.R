
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/process_sbl/sbl.sqrt.oc.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/full_data_trun6000/sub.calib.oc.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/full_data_trun6000/sub.valid.oc.RData")
obs <- sqrt(sub.calib.oc$OC)
pred <- sqrt(sub.valid.oc$OC)
neigh <- sbl.sqrt.oc$neighbors.stat$samples.id[,1:40]
calib.spc <- sub.calib.oc$spc[!is.na(obs),]
valid.spc <- sub.valid.oc$spc[!is.na(pred),]

com.spc <- rbind(calib.spc, valid.spc) #calib: 1:33793, valid: 33794:42305

#get pca
sbl.pca <- prcomp(com.spc)
save(sbl.pca, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/sbl.pca.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/sbl.pca.RData")
calib.scores <- sbl.pca$x[1:33793,]
valid.scores <- sbl.pca$x[33794:42305,]
scores <- sbl.pca$x
loadings <- sbl.pca$rotation
spc.recons <- scores %*% t(loadings)
obs.spc <- scale(com.spc, scale = FALSE)

plot(spc.recons[1,], type = "l")
lines(obs.spc[1,], col= "red")

spc.res <- (spc.recons - obs.spc)^2
spc.res <- sqrt(rowSums(spc.res))
for(i in 1:length(spc.res)){
  sample.fratio <- (length(spc.res)-1) * spc.res^2/sum((spc.res[-i])^2)
}

ok <- pf(sample.fratio, 1, length(spc.res)) 
rows <- which(ok>0.99)
plot(scores)
points(scores[rows,], col="red", pch=16)

#search for samples in validation sets if any of the neighbours contains outliers
valueFound <- apply(neigh, 1, function(x){
  if(any(x %in% rows)){
    1
  }
  else {
    0
  }
})

temp1 <- which(valueFound==1) ##index for outliers
##plot Xu scores with outliers
plot(sbl.sqrt.oc$pcAnalysis$scores_Xu)
points(sbl.sqrt.oc$pcAnalysis$scores_Xu[temp1,], col="red", pch=16)

##extract untrustworthy samples and all prediction samples from sbl
obs <- sbl.sqrt.oc$results$Nearest_neighbours_40$yu.obs^2
pred <- sbl.sqrt.oc$results$Nearest_neighbours_40$pred^2
udev <- sbl.sqrt.oc$results$Nearest_neighbours_40$ydev/sqrt(pred) * pred

sbl.oc.nooutliers <- data.frame(obs[-temp1],pred[-temp1],udev[-temp1])
sbl.oc.outliers <- data.frame(obs[temp1], pred[temp1], udev[temp1])
names(sbl.oc.outliers) <- c("obs", "pred", "udev")
names(sbl.oc.nooutliers) <- c("obs", "pred", "udev")
write.csv(sbl.oc.outliers, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/sbl/sbl.oc.outliers.csv")
write.csv(sbl.oc.nooutliers, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/sbl/sbl.oc.nooutliers.csv")

##output for score plot
write.csv(valid.scores[-temp1,1:2], file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pcplot/sbl.scores.nooutliers.csv")
write.csv(valid.scores[temp1,1:2], file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pcplot/sbl.scores.outliers.csv")







####repeat it for bulk density ###############

load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/BD/sbl/sbl.sqrt.all.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/BD/sub.calib.all.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/BD/sub.valid.all.RData")

neigh <- sbl.sqrt.all$neighbors.stat$samples.id[,1:20]

obs <- sqrt(sub.calib.all$analyte21)
pred <- sqrt(sub.valid.all$analyte21)

calib.spc <- sub.calib.all$spc[!is.na(obs),]
valid.spc <- sub.valid.all$spc[!is.na(pred),]

com.spc <- rbind(calib.spc, valid.spc) #calib: 1:33793, valid: 33794:42305

#get pca
sbl.bd.pca <- prcomp(com.spc)
save(sbl.bd.pca, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/sbl.bd.pca.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/sbl.bd.pca.RData")

calib.scores <- sbl.bd.pca$x[1:13895,]
valid.scores <- sbl.bd.pca$x[13896:17201,]
scores <- sbl.bd.pca$x
loadings <- sbl.bd.pca$rotation
spc.recons <- scores %*% t(loadings)
obs.spc <- scale(com.spc, scale = FALSE)

plot(spc.recons[1,], type = "l")
lines(obs.spc[1,], col= "red")

spc.res <- (spc.recons - obs.spc)^2
spc.res <- sqrt(rowSums(spc.res))
for(i in 1:length(spc.res)){
  sample.fratio <- (length(spc.res)-1) * spc.res^2/sum((spc.res[-i])^2)
}

ok <- pf(sample.fratio, 1, length(spc.res)) 
rows <- which(ok>0.99)
plot(scores)
points(scores[rows,], col="red", pch=16)

#search for samples in validation sets if any of the neighbours contains outliers
valueFound <- apply(neigh, 1, function(x){
  if(any(x %in% rows)){
    1
  }
  else {
    0
  }
})

temp1 <- which(valueFound==1) ##index for outliers
##plot Xu scores with outliers
#plot(sbl.sqrt.all$pcAnalysis$scores_Xu)
#points(sbl.sqrt.all$pcAnalysis$scores_Xu[temp1,], col="red", pch=16)

write.csv(valid.scores[-temp1,1:2], file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pcplot/sbl.bd.scores.nooutliers.csv")
write.csv(valid.scores[temp1,1:2], file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pcplot/sbl.bd.scores.outliers.csv")




##extract untrustworthy samples and all prediction samples from sbl
obs <- sbl.sqrt.all$results$Nearest_neighbours_20$yu.obs^2
pred <- sbl.sqrt.all$results$Nearest_neighbours_20$pred^2
udev <- sbl.sqrt.all$results$Nearest_neighbours_20$ydev/sqrt(pred) * pred

sbl.bd.nooutliers <- data.frame(obs[-temp1],pred[-temp1],udev[-temp1])
sbl.bd.outliers <- data.frame(obs[temp1], pred[temp1], udev[temp1])
names(sbl.bd.outliers) <- c("obs", "pred", "udev")
names(sbl.bd.nooutliers) <- c("obs", "pred", "udev")
write.csv(sbl.bd.outliers, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/sbl/sbl.bd.outliers.csv")
write.csv(sbl.bd.nooutliers, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/sbl/sbl.bd.nooutliers.csv")
