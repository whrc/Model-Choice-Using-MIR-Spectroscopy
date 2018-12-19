

rf.oc <- read.csv("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/uncertainty/rf.pred.oc.csv")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/full_data_trun6000/sub.valid.oc.RData")
tmp <- sqrt(sub.valid.oc$OC)
na.rows <-which(is.na(tmp)) 
rf.oc <- rf.oc[-na.rows,]

rf.oc$reld <- rf.oc$udev/rf.oc$pred
quantile(rf.oc$udev/rf.oc$pred, 0.99)  ##0.6385
index.rf.oc <- which(rf.oc$reld > 0.6385)

rf.oc.nooutliers <- rf.oc[-index.rf.oc,]
rf.oc.outliers <- rf.oc[index.rf.oc,]

write.csv(rf.oc.nooutliers, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/rf/rf.oc.nooutliers.csv")
write.csv(rf.oc.outliers, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/rf/rf.oc.outliers.csv")

##output oulier index for plot
write.csv(valid.scores[-index.rf.oc,1:2], file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pcplot/rf.scores.nooutliers.csv")
write.csv(valid.scores[index.rf.oc,1:2], file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pcplot/rf.scores.outliers.csv")




##RF bd
rf.bd <- read.csv("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/BD/rf/valid.pred.all.csv")

rf.bd$reld <- rf.bd$se/rf.bd$pred
quantile(rf.bd$se/rf.bd$pred, 0.99)  ##0.2158
index.rf.bd <- which(rf.bd$reld > 0.2158)

rf.bd.nooutliers <- rf.bd[-index.rf.bd,]
rf.bd.outliers <- rf.bd[index.rf.bd,]

write.csv(rf.bd.nooutliers, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/rf/rf.bd.nooutliers.csv")
write.csv(rf.bd.outliers, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/rf/rf.bd.outliers.csv")

write.csv(valid.scores[-index.rf.bd,1:2], file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pcplot/rf.bd.scores.nooutliers.csv")
write.csv(valid.scores[index.rf.bd,1:2], file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/Fratio/pcplot/rf.bd.scores.outliers.csv")

