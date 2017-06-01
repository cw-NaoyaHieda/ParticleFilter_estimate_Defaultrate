old_d<-getwd()
setwd("C:/Users/naoya/Desktop/Graduage_study/2016_ParticleFilter_estimate_Defaultrate/作業フォルダ/csv")
data<-read.csv("tanomu.csv",header=FALSE)
summary(data)
colnames(data)<-c("全体","建設","卸売","不動産","小売飲食","その他サービス","製造")
head(data)
dt=seq(as.Date("2001-03-01"),as.Date("2016-10-01"),by="month")
data<-cbind(dt,data)
gsub2 <- function(x){
  gsub("%","",x)
}
data[,c(2:8)]<-apply(data[,c(2:8)],2,gsub2)
data[,c(2:8)]<-apply(data[,c(2:8)],2,as.numeric)
data[,c(2:8)]<-data[,c(2:8)]/100
data[is.na(data)]<-100
data<-data[data$全体!=100,]
setwd(old_d)
