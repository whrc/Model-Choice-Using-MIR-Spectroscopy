

###PLSR Fratio for OC############
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/process_pls/fit.oc.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/full_data_trun6000/sub.valid.oc.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/full_data_trun6000/sub.calib.oc.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/sbl.pca.RData")
obs <- sqrt(sub.calib.oc$OC)
pred <- sqrt(sub.valid.oc$OC)

calib.spc <- sub.calib.oc$spc[!is.na(obs),]
valid.spc <- sub.valid.oc$spc[!is.na(pred),]

com.spc <- rbind(calib.spc, valid.spc) #calib: 1:33793, valid: 33794:42305


all.scores <- sbl.pca$x
all.loads <- sbl.pca$rotation

#scores <- fit.all$scores
#loads <- fit.all$loadings
#weights <- fit.all$loading.weights
spc.recons.test1 <- all.scores %*% t(all.loads)
spc.recons.test <- spc.recons.test1

spc.orig.all <- scale(com.spc, scale = FALSE)
plot(spc.recons.test[1,], type = "l", ylim = c(-0.2,0.2))
lines(spc.orig.all[1,], col="red")

spc.res <- (spc.recons.test - spc.orig.all)^2
spc.res <- sqrt(rowSums(spc.res))



for(i in 1:length(spc.res)){
  sample.fratio <- (length(spc.res)-1) * spc.res^2/sum((spc.res[-i])^2)
}

ok <- pf(sample.fratio, 1, length(spc.res)) 
rows <- which(ok>0.99)
plot(all.scores, xlim = c(-25,30), ylim = c(-15,15))
points(all.scores[rows,1], all.scores[rows,2], col="blue", pch=16)
#points(fit.all$scores[rows.valid,1], fit.all$scores[rows.valid,2], col="red", pch=16)

#extract validation sets
valid.scores <- all.scores[33794:42305,]
valid.ok <- ok[33794:42305]
valid.rows <- which(valid.ok > 0.99)
plot(valid.scores)
points(valid.scores[valid.rows,], col="red", pch = 16)

write.csv(valid.scores[-valid.rows,1:2], file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pcplot/pls.scores.nooutliers.csv")
write.csv(valid.scores[valid.rows,1:2], file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pcplot/pls.scores.outliers.csv")

pls.valid <- read.csv("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/process_pls/valid.pred.oc.csv")
#remove nas to match with sbl prediction
tmp <- sqrt(sub.valid.oc$OC)
na.rows <-which(is.na(tmp)) 
pls.valid <- pls.valid[-na.rows,]
pls.oc.nooutliers <- pls.valid[-valid.rows,]
pls.oc.outliers <- pls.valid[valid.rows,]
write.csv(pls.oc.nooutliers, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pls/pls.oc.nooutliers.csv")
write.csv(pls.oc.outliers, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pls/pls.oc.outliers.csv")








####Bulk Density#####################
###################################################
###PLSR Fratio for OC############
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/BD/pls/fit.all.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/BD/sub.valid.all.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/BD/sub.calib.all.RData")

load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/sbl.bd.pca.RData")
all.scores <- sbl.bd.pca$x
all.loads <- sbl.bd.pca$rotation
com.all <- rbind(sub.calib.all$spc, sub.valid.all$spc)
#scores <- fit.all$scores
#loads <- fit.all$loadings
#weights <- fit.all$loading.weights
spc.recons.test1 <- all.scores %*% t(all.loads)
spc.recons.test <- spc.recons.test1

spc.orig.all <- scale(com.all, scale = FALSE)
plot(spc.recons.test[1,], type = "l", ylim = c(-0.2,0.2))
lines(spc.orig.all[1,], col="red")

spc.res <- (spc.recons.test - spc.orig.all)^2
spc.res <- sqrt(rowSums(spc.res))



for(i in 1:length(spc.res)){
  sample.fratio <- (length(spc.res)-1) * spc.res^2/sum((spc.res[-i])^2)
}

ok <- pf(sample.fratio, 1, length(spc.res)) 
rows <- which(ok>0.99)
plot(all.scores, xlim = c(-25,30), ylim = c(-15,15))
points(all.scores[rows,1], all.scores[rows,2], col="blue", pch=16)
#points(fit.all$scores[rows.valid,1], fit.all$scores[rows.valid,2], col="red", pch=16)

#extract validation sets
valid.scores <- all.scores[13896:17201,]
valid.ok <- ok[13896:17201]
valid.rows <- which(valid.ok > 0.99)
plot(valid.scores)
points(valid.scores[valid.rows,], col="red", pch = 16)

##write outliers file for making plot
write.csv(valid.scores[-valid.rows,1:2], file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pcplot/sbl.bd.scores.nooutliers.csv")
write.csv(valid.scores[valid.rows,1:2], file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pcplot/sbl.bd.scores.outliers.csv")



pls.valid <- read.csv("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/BD/pls/all.pred.csv")
pls.bd.nooutliers <- pls.valid[-valid.rows,]
pls.bd.outliers <- pls.valid[valid.rows,]
write.csv(pls.bd.nooutliers, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pls/pls.bd.nooutliers.csv")
write.csv(pls.bd.outliers, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pls/pls.bd.outliers.csv")

