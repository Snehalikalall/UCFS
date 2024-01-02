library('copula')
library(reshape2)
library(foreach)
library(doParallel)
library('prodlim')
library(Matrix)

dataorg=read.csv("kleinorg.csv", header=FALSE)
data=read.csv("kleinsemifea.csv", header=FALSE)
#Feature Selection
cl <- makeCluster(40)
registerDoParallel(cl)
nf=100
set.seed(1000)
datas2<-t(data)
n=nrow(datas2)
count=ncol(datas2)
start.time <- Sys.time()
# Feature by CBFS parallel 
fea<- matrix(0, nrow=1,ncol =nf)
fea[1,1]=1
theta=-0.5
fc=claytonCopula(theta,dim=2)
for (m in 2:nf)
{
  feas<-fea[1,(m-1)]
  parl<-foreach(j=1:count, .combine=c,.packages='copula') %dopar%
    {
      u<-pobs(cbind(datas2[,feas],datas2[,j]))
      a=pCopula(u,fc)
      res=mean(a)
    }
  result<-as.matrix(parl) 
  result[fea]<-1
  #feamid<-which(result==min(result[which(result>0)]))[1]
  feamid=which.min(result)
  fea[1,m]<-feamid
  print(m)
}

stopCluster(cl)
end.time <- Sys.time()  
time.taken <- end.time - start.time
time.taken  
registerDoSEQ()


#take main data to match feature, cells in row, genes in coloumn
#dataorg= as.matrix(read.csv("baron.csv", header=FALSE))
dataorg=as.matrix(t(dataorg))
datasemfea<-as.data.frame(data)
index875<-row.match(datasemfea,dataorg)
feabiase<-as.matrix(index875[as.matrix(fea)])
write.table(fea,file="ucfsfea.csv",sep=",",row.names = FALSE,col.names = FALSE)

