# testing project
rm(list = ls())
library(raster)
library(geosphere)
shp <- shapefile('C:/Users/joeja/Desktop/MATH/MATH522/Project/CANUS_Boundary_Files/CANUS_Watershed_Boundaries.shp')


#rotate a given shape
changeshape<-function(cur_x,cur_y,curloc_x,curloc_y,change,type){
  if(type=="rotate"){
    changed_x<-cur_x*cos(change)-cur_y*sin(change)
    changed_y<-cur_x*sin(change)+cur_y*cos(change)
    changed_loc_x<-curloc_x*cos(change)-curloc_y*sin(change)
    changed_loc_y<-curloc_x*sin(change)+curloc_y*cos(change)
  }
  if(type=="resize"){
    changed_x<-cur_x*change
    changed_y<-cur_y*change
    changed_loc_x<-curloc_x*change
    changed_loc_y<-curloc_y*change
  }
  if(type=="shear"){
    changemat<-diag(2)
    changemat[1,1]<-change^2+changemat[1,1]
    changemat[2,1]<-change+changemat[2,1]
    changemat[1,2]<-change+changemat[1,2]
    newmat<-changemat%*%matrix(c(cur_x,cur_y),nrow=2,byrow=T)
    newloc<-changemat%*%matrix(c(curloc_x,curloc_y),nrow=2,byrow=T)
    changed_x<-as.numeric(newmat[1,])
    changed_y<-as.numeric(newmat[2,])
    changed_loc_x<-newloc[1]
    changed_loc_y<-newloc[2]
  }
  if(type=="start"){
    newindex<- sample(length(cur_x),1)
    changed_x=cur_x
    changed_y=cur_y
    changed_loc_x=cur_x[newindex]
    changed_loc_y=cur_y[newindex]
  }
  if(type=="reflect"){
    changed_x= cur_x
    changed_y= -cur_y
    changed_loc_x=curloc_x
    changed_loc_y=curloc_y
  }
  
  return(list(changed_x=changed_x,changed_y=changed_y,changed_loc_x=changed_loc_x,changed_loc_y=changed_loc_y))
}

make_size_inv<- function(z,type){
  if(type=="sd"){
    oldsd<-norm(z-mean(z),type="2")
    newz<-z/oldsd
  }
  if(type=="f1"){
    zfft=fft(z)
    zfft<-zfft/abs(zfft[2])
    newz<-fft(zfft,inverse = T) * (1/length(zfft))
  }
  newz
}

make_rotate_inv<- function(z,type){
  if(type=="rot"){
    zfft=fft(z)
    c1_dir=zfft[2]/abs(zfft[2])
    cm1_dir=zfft[length(z)]/abs(zfft[length(z)])
    ang_c1=atan(Im(c1_dir)/Re(c1_dir))
    ang_cm1=atan(Im(cm1_dir)/Re(cm1_dir))
    zfft<-zfft * exp(sqrt(-1+0i) * -(ang_c1+ang_cm1)/2)
    newz<-fft(zfft,inverse = T) * (1/length(zfft))
  }
  if(type=="rot_start"){
    zfft=fft(z)
    c1_dir=zfft[2]/abs(zfft[2])
    cm1_dir=zfft[length(z)]/abs(zfft[length(z)])
    k=c(0:floor(length(z)/2),-floor(length(z)/2):-1)
    ang_c1=atan(Im(c1_dir)/Re(c1_dir))
    ang_cm1=atan(Im(cm1_dir)/Re(cm1_dir))
    zfft<-zfft * exp(sqrt(-1+0i) * -(ang_c1+ang_cm1)/2)
    zfft<-zfft * exp(sqrt(-1+0i) * k* -(ang_c1-ang_cm1)/2) #for starting point invariance
    newz<-fft(zfft,inverse = T) * (1/length(zfft))
  }
  newz
}

make_shear_inv<- function(z,type){
  library(fastICA)
  if(type=="ICA"){
    z=cbind(Re(z),Im(z))
    ICA<-fastICA(z,n.comp = 2)
    newz=ICA$S
    newz= complex(length.out = nrow(newz),real=newz[,1],imaginary=newz[,2])
  }
  newz
}

make_reflect_inv<- function(z, type){
  if(type=="sp"){
    newy= Im(z)
    if(Im(z)[1]> median(Im(z))) newy= - Im(z)
    newz= complex(length.out = length(z),real=Re(z),imaginary=newy)
  }
  newz
}


#plot roation examples
jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH522/Project/rotations.jpg"),width = 7, height = 7, units = 'in',res = 600)
par(mfrow=c(2,2),mai = c(.1, 0.1, 0.1, 0.1))
par(pty="s")
changes<-c(0,30,110,200)
s<-1
curshp<-shp@polygons[[s]]@Polygons[[1]]@coords
#curshp<-matrix(c(1,1,-1,-1,-1,1,-1,1),ncol=2)
curloc<-c(shp@data$Longitude[s],shp@data$Latitude[s])
for(i in 1:length(changes)){
  changedshape<-changeshape(cur_x=curshp[,1],cur_y = curshp[,2],curloc_x = curloc[1],curloc_y = curloc[2],change = changes[i],type = "rotate")
  plot(changedshape$changed_x,changedshape$changed_y,type="l",asp=1, xaxt='n',yaxt='n', ann=FALSE, frame.plot=F)
  points(changedshape$changed_loc_x,changedshape$changed_loc_y,col="red",pch=20,cex=3)
}
dev.off()


#plot size invariance examples
jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH522/Project/sizes.jpg"),width = 7, height = 7, units = 'in',res = 600)
par(mfrow=c(2,2),mai = c(.1, 0.1, 0.1, 0.1))
par(pty="s")
changes<-c(.25,.5,.85,1)
s<-1
curshp<-shp@polygons[[s]]@Polygons[[1]]@coords
#curshp<-matrix(c(1,1,-1,-1,-1,1,-1,1),ncol=2)
curloc<-c(shp@data$Longitude[s],shp@data$Latitude[s])
for(i in 1:length(changes)){
  changedshape<-changeshape(cur_x=curshp[,1],cur_y = curshp[,2],curloc_x = curloc[1],curloc_y = curloc[2],change = changes[i],type="resize")
  plot(changedshape$changed_x-mean(changedshape$changed_x),changedshape$changed_y-mean(changedshape$changed_y),type="l",xlim=c(-3,3),ylim=c(-1,1),
       asp=1, xaxt='n',yaxt='n', ann=FALSE, frame.plot=F)
  points(changedshape$changed_loc_x-mean(changedshape$changed_x),changedshape$changed_loc_y-mean(changedshape$changed_y),col="red",pch=20,cex=3)
}
dev.off()



#plot invariant shear examples
jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH522/Project/shears.jpg"),width = 7, height = 7, units = 'in',res = 600)
par(mfrow=c(2,2),mai = c(.1, 0.1, 0.1, 0.1))
par(pty="s")
changes<-c(0,.4,.7,1)
s<-1
curshp<-shp@polygons[[s]]@Polygons[[1]]@coords
#curshp<-matrix(c(1,1,-1,-1,-1,1,-1,1),ncol=2)
curloc<-c(shp@data$Longitude[s],shp@data$Latitude[s])
for(i in 1:length(changes)){
  changedshape<-changeshape(cur_x=curshp[,1],cur_y = curshp[,2],curloc_x = curloc[1],curloc_y = curloc[2],change = changes[i],type="shear")
  plot(changedshape$changed_x,changedshape$changed_y,type="l",
       asp=1, xaxt='n',yaxt='n', ann=FALSE, frame.plot=F)
  points(changedshape$changed_loc_x,changedshape$changed_loc_y,col="red",pch=20,cex=3)
}
dev.off()


#plot starting point invariance example
jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH522/Project/start.jpg"),width = 7, height = 7, units = 'in',res = 600)
par(mfrow=c(2,2),mai = c(.1, 0.1, 0.1, 0.1),pty="s")
changes<-c(0,.4,.7,1)
s<-1
curshp<-shp@polygons[[s]]@Polygons[[1]]@coords
#curshp<-matrix(c(1,1,-1,-1,-1,1,-1,1),ncol=2)
curloc<-c(shp@data$Longitude[s],shp@data$Latitude[s])
for(i in 1:length(changes)){
  changedshape<-changeshape(cur_x=curshp[,1],cur_y = curshp[,2],curloc_x = curloc[1],curloc_y = curloc[2],change = changes[i],type="start")
  plot(changedshape$changed_x,changedshape$changed_y,type="l",
       asp=1, xaxt='n',yaxt='n', ann=FALSE, frame.plot=F)
  points(changedshape$changed_loc_x,changedshape$changed_loc_y,col="red",pch=20,cex=3)
}
dev.off()


#ica demo

jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH522/Project/ICAdemo.jpg"),width = 5, height = 5, units = 'in',res = 600)

par(mfrow=c(2,2))
s1<-runif(10000,min=-sqrt(3),max=sqrt(3))
s2<-runif(10000,min=-sqrt(3),max=sqrt(3))

plot(s1,s2)

s=rbind(s1,s2)

A=matrix(c(.96,.28,-.28,.96),nrow=2)

x=A %*% s

plot(x[1,],x[2,])

#test for gaussian
s1<-rnorm(10000)
s2<-rnorm(10000)

plot(s1,s2)

s=rbind(s1,s2)

A=matrix(c(.96,.28,-.28,.96),nrow=2)

x=A %*% s

plot(x[1,],x[2,])

dev.off()




############################################################################################
########################### start hydrology application  ##################################
############################################################################################
rm(list = ls())
library(raster)
library(geosphere)
library(ranger)
library(xgboost)
shp <- shapefile('C:/Users/joeja/Desktop/MATH/MATH522/Project/CANUS_Boundary_Files/CANUS_Watershed_Boundaries.shp')

#import hongyi data
Hongyi_dataset <- read.csv("C:/Users/joeja/Desktop/MATH/MATH522/Project/Hongyi_dataset.csv")

#clean and organize data
Hongyi_dataset$Event<-NULL
Hongyi_dataset$a_event<-NULL
Hongyi_dataset$IQR_loga<-NULL
Hongyi_dataset$IQR_b<-NULL
Hongyi_dataset$loga_pointcloud<-NULL
Hongyi_dataset$a_pointcloud<-NULL
Hongyi_dataset$b_pointcloud<-NULL
Hongyi_dataset$Storage_event<-NULL
Hongyi_dataset$Storage_pointcloud<-NULL
Hongyi_dataset$Curvature<-NULL
Hongyi_dataset$Convexity<-NULL
Hongyi_dataset$X.15<-NULL
Hongyi_dataset$X.50<-NULL
Hongyi_dataset$X.85<-NULL


apply(Hongyi_dataset,2, function(x) sum(is.na(x)))

Hongyi_dataset$DrainageDensity<-NULL


#predicting log(a)

nfold=5
folds= sample(1:nfold,nrow(Hongyi_dataset),replace = T)
rfa<-rep(0,nfold)
rfb<-rep(0,nfold)

xgba<-rep(0,nfold)
xgbb<-rep(0,nfold)
for(fold in 1:nfold){
  train<-Hongyi_dataset[folds!=fold,]
  test<-Hongyi_dataset[folds==fold,]
  
  rfmod<-ranger(y=train$loga_event,x=train[,4:ncol(train)])
  preds<- predict(rfmod,test)$predictions
  rfa[fold]<-1-(sum((preds-test$loga_event)^2))/(sum((test$loga_event-mean(test$loga_event))^2))
  
  rfmod<-ranger(y=train$b_event,x=train[,4:ncol(train)])
  preds<- predict(rfmod,test)$predictions
  rfb[fold]<-1-(sum((preds-test$b_event)^2))/(sum((test$b_event-mean(test$b_event))^2))
  
  mod<-xgb.train(list(eta=.018,subsample=.839,max_depth=13,min_child_weight=2.06,colsample_bytree=.752,colsample_bylevel=.585)
                 ,xgb.DMatrix(as.matrix(train[,4:ncol(train)]),label=train$loga_event),nrounds=4168,verbose = 0)
  preds<-predict(mod,xgb.DMatrix(as.matrix(test[,4:ncol(train)])))
  xgba[fold]<-1-(sum((preds-test$loga_event)^2))/(sum((test$loga_event-mean(test$loga_event))^2))
  
  mod<-xgb.train(list(eta=.018,subsample=.839,max_depth=13,min_child_weight=2.06,colsample_bytree=.752,colsample_bylevel=.585)
                 ,xgb.DMatrix(as.matrix(train[,4:ncol(train)]),label=train$b_event),nrounds=4168,verbose = 0)
  preds<-predict(mod,xgb.DMatrix(as.matrix(test[,4:ncol(train)])))
  xgbb[fold]<-1-(sum((preds-test$b_event)^2))/(sum((test$b_event-mean(test$b_event))^2))
}
mean(rfa)
mean(rfb)
mean(xgba)
mean(xgbb)

rfmod<-ranger(y=Hongyi_dataset$loga_event,x=Hongyi_dataset[,4:ncol(Hongyi_dataset)])
rfmod$r.squared

#predicting b
rfmod<-ranger(y=Hongyi_dataset$b_event,x=Hongyi_dataset[,4:ncol(Hongyi_dataset)])
rfmod$r.squared


#plot a few random catchmnets
set.seed(1)
randcat<-sample(Hongyi_dataset$gridcode,16)


jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH522/Project/CatchmentShapesEX.jpg"),width = 7, height = 5, units = 'in',res = 600)
par(mfrow=c(4,4),mai = c(.1, 0, 0.1, 0),pty="s")
for(i in 1:length(randcat)){
  gc<-randcat[i]
  which_shp<-which(as.numeric(shp@data$gridcode)==gc)
  curshp<-shp@polygons[[which_shp]]@Polygons
  
  curshp<-curshp[[which.max(lapply(curshp, object.size))]]@coords
  #curloc<-c(shp@data$Longitude[which_shp],shp@data$Latitude[which_shp])
  curloc<-c(Hongyi_dataset$Longitude[Hongyi_dataset$gridcode==gc],Hongyi_dataset$Latitude[Hongyi_dataset$gridcode==gc])
  
  plot(curshp,type="l",
       asp=1, xaxt='n',yaxt='n', ann=FALSE, frame.plot=F)
  points(curloc[1],curloc[2],col="red",pch=20,cex=3)
}
dev.off()


shplist<-list()
for(i in 1:nrow(Hongyi_dataset)){
  gc<-Hongyi_dataset$gridcode[i]
  which_shp<-which(as.numeric(shp@data$gridcode)==gc)
  curshp<-shp@polygons[[which_shp]]@Polygons
  curshp<-curshp[[which.max(lapply(curshp, object.size))]]@coords
  curloc<-c(Hongyi_dataset$Longitude[Hongyi_dataset$gridcode==gc],Hongyi_dataset$Latitude[Hongyi_dataset$gridcode==gc])
  start<-which.min((curshp[,1]-curloc[1])^2+(curshp[,2]-curloc[2])^2)
  curshp<-curshp[c(start:nrow(curshp),1:(start-1)),]
  
  z=complex(length.out = nrow(curshp),real=curshp[,1],imaginary=curshp[,2])
  shplist[[i]]<-z
}

#test methods for size invariance
jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH522/Project/size_inv_sd.jpg"),width = 4, height = 7, units = 'in',res = 600)
par(mfrow=c(2,1),mai = c(0, 0, 0, 0),pty="s")

curshp=shplist[[1]]-mean(shplist[[1]])
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 1.2,type="resize")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
plot(newz,type="l",col="red",asp=1, xaxt='n',yaxt='n', frame.plot=F,ann=F,lty="dotdash")
lines(curshp,col="black",lty = "dashed")
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 0.7,type="resize")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
lines(newz,col="blue",lty = "dotted")

curshp=shplist[[1]]-mean(shplist[[1]])
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 1.2,type="resize")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
newz=make_size_inv(newz,"sd")
plot(newz,type="l",col="red",asp=1, xaxt='n',yaxt='n', ann=FALSE, frame.plot=F,lty="dotdash")
newz=make_size_inv(curshp,"sd")
lines(newz,col="black",lty = "dashed")
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 0.7,type="resize")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
newz=make_size_inv(newz,"sd")
lines(newz,col="blue",lty = "dotted")
dev.off()


jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH522/Project/size_inv_f1.jpg"),width = 4, height = 7, units = 'in',res = 600)
par(mfrow=c(2,1),mai = c(0, 0, 0, 0),pty="s")

curshp=shplist[[1]]-mean(shplist[[1]])
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 1.2,type="resize")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
plot(newz,type="l",col="red",asp=1, xaxt='n',yaxt='n', frame.plot=F,ann=F,lty="dotdash")
lines(curshp,col="black",lty = "dashed")
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 0.7,type="resize")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
lines(newz,col="blue",lty = "dotted")

curshp=shplist[[1]]-mean(shplist[[1]])
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 1.2,type="resize")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
newz=make_size_inv(z=newz,type="f1")
plot(newz,type="l",col="red",asp=1, xaxt='n',yaxt='n', ann=FALSE, frame.plot=F,lty="dotdash")
newz=make_size_inv(curshp,"f1")
lines(newz,col="black",lty = "dashed")
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 0.7,type="resize")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
newz=make_size_inv(newz,"f1")
lines(newz,col="blue",lty = "dotted")
dev.off()



#test methods for rotation invariance
jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH522/Project/rot_inv_can.jpg"),width = 4, height = 7, units = 'in',res = 600)
par(mfrow=c(2,1),mai = c(0, 0, 0, 0),pty="s")

curshp=shplist[[1]]-mean(shplist[[1]])
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 1.2,type="rotate")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
newz<-newz[c(700:length(newz),1:699)] #change starting point
plot(newz,type="l",col="red",asp=1, xaxt='n',yaxt='n', frame.plot=F,ann=F,lty="dotdash",ylim=c(-1,1),xlim=c(-1,1))
points(newz[1],col="red",pch=20,cex=3)
lines(curshp,col="black",lty = "dashed")
points(curshp[1],col="black",pch=20,cex=3)
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 4,type="rotate")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
newz<-newz[c(5000:length(newz),1:4999)] #change starting point
lines(newz,col="blue",lty = "dotted")
points(newz[1],col="blue",pch=20,cex=3)

curshp=shplist[[1]]-mean(shplist[[1]])
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 1.2,type="rotate")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
newz<-newz[c(700:length(newz),1:699)] #change starting point
newz=make_rotate_inv(z=newz,type="rot")
plot(newz,type="l",col="red",asp=1, xaxt='n',yaxt='n', ann=FALSE, frame.plot=F,lty="dotdash",ylim=c(-1,1),xlim=c(-1,1))
points(newz[1],col="red",pch=20,cex=3)
newz=make_rotate_inv(curshp,"rot")
lines(newz,col="black",lty = "dashed")
points(newz[1],col="black",pch=20,cex=3)
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 0.7,type="rotate")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
newz<-newz[c(5000:length(newz),1:4999)] #change starting point
newz=make_rotate_inv(newz,"rot")
lines(newz,col="blue",lty = "dotted")
points(newz[1],col="blue",pch=20,cex=3)
dev.off()


jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH522/Project/rot_start_inv_can.jpg"),width = 4, height = 7, units = 'in',res = 600)
par(mfrow=c(2,1),mai = c(0, 0, 0, 0),pty="s")

curshp=shplist[[1]]-mean(shplist[[1]])
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 1.2,type="rotate")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
newz<-newz[c(700:length(newz),1:699)] #change starting point
plot(newz,type="l",col="red",asp=1, xaxt='n',yaxt='n', frame.plot=F,ann=F,lty="dotdash",ylim=c(-1,1),xlim=c(-1,1))
points(newz[1],col="red",pch=20,cex=3)
lines(curshp,col="black",lty = "dashed")
points(curshp[1],col="black",pch=20,cex=3)
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 4,type="rotate")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
newz<-newz[c(5000:length(newz),1:4999)] #change starting point
lines(newz,col="blue",lty = "dotted")
points(newz[1],col="blue",pch=20,cex=3)

curshp=shplist[[1]]-mean(shplist[[1]])
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 1.2,type="rotate")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
newz<-newz[c(700:length(newz),1:699)] #change starting point
newz=make_rotate_inv(z=newz,type="rot_start")
plot(newz,type="l",col="red",asp=1, xaxt='n',yaxt='n', ann=FALSE, frame.plot=F,lty="dotdash",ylim=c(-1,1),xlim=c(-1,1))
points(newz[1],col="red",pch=20,cex=3)
newz=make_rotate_inv(curshp,"rot_start")
lines(newz,col="black",lty = "dashed")
points(newz[1],col="black",pch=20,cex=3)
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 0.7,type="rotate")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
newz<-newz[c(5000:length(newz),1:4999)] #change starting point
newz=make_rotate_inv(newz,"rot_start")
lines(newz,col="blue",lty = "dotted")
points(newz[1],col="blue",pch=20,cex=3)
dev.off()



#test methods for shear invariance
jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH522/Project/shear_inv_ICA.jpg"),width = 4, height = 7, units = 'in',res = 600)
par(mfrow=c(2,1),mai = c(0, 0, 0, 0),pty="s")

curshp=shplist[[1]]-mean(shplist[[1]])
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 1,type="shear")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
plot(newz,type="l",col="red",asp=1, xaxt='n',yaxt='n', frame.plot=F,ann=F,lty="dotdash")
lines(curshp,col="black",lty = "dashed")
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = -0.7,type="shear")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
lines(newz,col="blue",lty = "dotted")

curshp=shplist[[1]]-mean(shplist[[1]])
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 1,type="shear")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
newz=make_shear_inv(z=newz,type="ICA")
plot(newz,type="l",col="red",asp=1, xaxt='n',yaxt='n', ann=FALSE, frame.plot=F,lty="dotdash",ylim=c(-2,2),xlim=c(-2,2))
newz=make_shear_inv(curshp,"ICA")
lines(newz,col="black",lty = "dashed")
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = -0.7,type="shear")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
newz=make_shear_inv(newz,"ICA")
lines(newz,col="blue",lty = "dotted")
dev.off()



#test methods for reflection invariance
jpeg(paste0("C:/Users/joeja/Desktop/MATH/MATH522/Project/reflect_inv.jpg"),width = 7, height = 4, units = 'in',res = 600)
par(mfrow=c(1,2),mai = c(0, 0, 0, 0),pty="s")

curshp=shplist[[1]]-mean(shplist[[1]])
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 1,type="reflect")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
plot(newz,type="l",col="red",asp=1, xaxt='n',yaxt='n', frame.plot=F,ann=F,lty="dotdash")
points(newz[1],col="red",pch=20,cex=3)
lines(curshp,col="black",lty = "dashed")
points(curshp[1],col="black",pch=20,cex=3)

curshp=shplist[[1]]-mean(shplist[[1]])
changedshape<-changeshape(cur_x=Re(curshp),cur_y = Im(curshp),curloc_x = curshp[1],curloc_y = curshp[2],change = 1,type="reflect")
newz=complex(length.out = length(curshp),real=changedshape[[1]],imaginary=changedshape[[2]])
newz=make_reflect_inv(z=newz,type="sp")
plot(newz,type="l",col="red",asp=1, xaxt='n',yaxt='n', ann=FALSE, frame.plot=F,lty="dotdash")
points(newz[1],col="red",pch=20,cex=3)
newz=make_reflect_inv(curshp,"sp")
lines(newz,col="black",lty = "dashed")
points(newz[1],col="black",pch=20,cex=3)
dev.off()



#make everything equal size
what<-lapply(shplist,length)
bigshp<-max(as.numeric(what))
powersof2<-2^(1:30)
newlength<-powersof2[which(bigshp<powersof2)[1]]
shp_pad<-list()
for(i in 1:length(shplist)){
  curshp<-shplist[[i]]
  newshp<-rep(0,newlength)
  newshp[1:length(curshp)]<-curshp
  shp_pad[[i]]<-newshp
}


#get first numFD fourier coefficients for each shape
numFD<-10
shp_fd_cc<-matrix(0,nrow = length(shp_pad),ncol = numFD-1)
for(i in 1:nrow(shp_fd_cc)){
  curfft<- fft(shp_pad[[i]])
  shp_fd_cc[i,]<-abs(curfft[c(3:(1+numFD/2),length(curfft):(length(curfft)-(-1+numFD/2)))])/abs(curfft[2])
}

Hongyi_dataset_fd_cc<-cbind(Hongyi_dataset,shp_fd_cc)

rfmod<-ranger(y=Hongyi_dataset$loga_event,x=Hongyi_dataset[,4:(ncol(Hongyi_dataset))],importance = "permutation")
rfmod$r.squared
rfmod$variable.importance

rfmod<-ranger(y=Hongyi_dataset$loga_event,x=Hongyi_dataset_fd_cc[,4:(ncol(Hongyi_dataset_fd_cc))],importance = "permutation")
rfmod$r.squared
rfmod$variable.importance

rfmod<-ranger(y=Hongyi_dataset$b_event,x=Hongyi_dataset[,4:(ncol(Hongyi_dataset))],importance = "permutation")
rfmod$r.squared
rfmod$variable.importance

rfmod<-ranger(y=Hongyi_dataset$b_event,x=Hongyi_dataset_fd_cc[,4:(ncol(Hongyi_dataset_fd_cc))],importance = "permutation")
rfmod$r.squared
rfmod$variable.importance


#distance to outlet parameterization
shplist_dt<-list()
for(i in 1:nrow(Hongyi_dataset)){
  gc<-Hongyi_dataset$gridcode[i]
  which_shp<-which(as.numeric(shp@data$gridcode)==gc)
  curshp<-shp@polygons[[which_shp]]@Polygons
  curshp<-curshp[[which.max(lapply(curshp, object.size))]]@coords
  curloc<-c(Hongyi_dataset$Longitude[Hongyi_dataset$gridcode==gc],Hongyi_dataset$Latitude[Hongyi_dataset$gridcode==gc])
  start<-which.min((curshp[,1]-curloc[1])^2+(curshp[,2]-curloc[2])^2)
  curshp<-curshp[c(start:nrow(curshp),1:(start-1)),]
  
  d=rowSums(t(t(curshp)-curloc)^2)
  shplist_dt[[i]]<-d
}

what<-lapply(shplist_dt,length)
bigshp<-max(as.numeric(what))
powersof2<-2^(1:30)
newlength<-powersof2[which(bigshp<powersof2)[1]]
shp_pad_dt<-list()
for(i in 1:length(shplist_dt)){
  curshp<-shplist_dt[[i]]
  newshp<-rep(0,newlength)
  newshp[1:length(curshp)]<-curshp
  shp_pad_dt[[i]]<-newshp
}


numFD<-10
shp_fd_dt<-matrix(0,nrow = length(shp_pad_dt),ncol = numFD)
for(i in 1:nrow(shp_fd_dt)){
  curfft<- fft(shp_pad_dt[[i]])
  shp_fd_dt[i,]<-abs(curfft[2:(numFD+1)])/abs(curfft[1])
}

Hongyi_dataset_fd_dt<-cbind(Hongyi_dataset,shp_fd_dt)

rfmod<-ranger(y=Hongyi_dataset$loga_event,x=Hongyi_dataset_fd_dt[,4:(ncol(Hongyi_dataset_fd_dt))],importance = "permutation")
rfmod$r.squared
rfmod$variable.importance


rfmod<-ranger(y=Hongyi_dataset$b_event,x=Hongyi_dataset_fd_dt[,4:(ncol(Hongyi_dataset_fd_dt))],importance = "permutation")
rfmod$r.squared
rfmod$variable.importance







#calculate complexity
Hongyi_dataset$shp_complexity<-NA

poss_numfor<-seq(from=3,by=2,to=101)
error_max_norm<-matrix(0,nrow=nrow(Hongyi_dataset),ncol=length(poss_numfor))
error_2_norm<-matrix(0,nrow=nrow(Hongyi_dataset),ncol=length(poss_numfor))

for(i in 1:length(shp_pad)){
  curshp<-shplist[[i]]
  zfft=fft(curshp)
  
  for(j in 1:length(poss_numfor)){
    numFor<-poss_numfor[j]
    zfft_thresholded<-rep(0,length(zfft))
    to_keep<-c(seq(from=1,by=1,ceiling(numFor/2)),seq(from=length(zfft),by=-1,length.out = floor(numFor/2)))
    zfft_thresholded[to_keep]<-zfft[to_keep]
    estimate= fft(zfft_thresholded,inverse = T)/length(zfft_thresholded)
    error_max_norm[i,j]= max(abs(curshp- estimate))
    error_2_norm[i,j]= mean(abs(curshp- estimate)^2)
  }
  error_max_norm[i,]= error_max_norm[i,]/error_max_norm[i,1]
  error_2_norm[i,]=error_2_norm[i,]/error_2_norm[i,1]
  print(i)
}

for(i in 1:nrow(Hongyi_dataset)){
  x=log(poss_numfor)
  y=log(error_2_norm[i,])
  x<-x[y>-14]
  y<-y[y>-14]
  mod<-lm(y~x)
  Hongyi_dataset$shp_complexity[i]<-as.numeric(mod$coefficients[2])
}


rfmod<-ranger(y=Hongyi_dataset$loga_event,x=Hongyi_dataset[,4:(ncol(Hongyi_dataset))],importance = "permutation")
rfmod$r.squared
rfmod$variable.importance

#predicting b
rfmod<-ranger(y=Hongyi_dataset$b_event,x=Hongyi_dataset[,4:(ncol(Hongyi_dataset))],importance = "permutation")
rfmod$r.squared
rfmod$variable.importance



min_grid<-which.min(Hongyi_dataset$shp_complexity)
max_grid<-which.max(Hongyi_dataset$shp_complexity)

plot(shplist[[min_grid]],type="l")
plot(shplist[[max_grid]],type="l")

look<-which(Hongyi_dataset$shp_complexity< -2.5)
plot(shplist[[look[4]]],type="l")


i=look[1]
x=log(poss_numfor)
y=log(error_2_norm[i,])
x<-x[y>-14]
y<-y[y>-14]
plot(x,y)












##################################################################################################################
#########################3     old code for camels dataset       ##################################################
#################################################################################################################


poss_numfor<-seq(from=1,by=2,to=101)
error_max_norm<-matrix(0,nrow=nrow(dat),ncol=length(poss_numfor))
error_2_norm<-matrix(0,nrow=nrow(dat),ncol=length(poss_numfor))

for(i in 1:nrow(dat)){
  ref_point<-c(dat$gauge_lon[i],dat$gauge_lat[i])
  cur_shp<-shp[i,]
  cur_shp_boundary<-cur_shp@polygons[[1]]@Polygons
  cur_shp_boundary<-cur_shp_boundary[[which.max(lapply(cur_shp_boundary, object.size))]]@coords
  
  #get distance in x direction
  cur_shp_boundary_x<-cur_shp_boundary
  cur_shp_boundary_x[,2]<-ref_point[2]
  x_dist<-distHaversine(ref_point,cur_shp_boundary_x)*.001
  x_dist[cur_shp_boundary_x[,1]<ref_point[1]]<-x_dist[cur_shp_boundary_x[,1]<ref_point[1]]*-1
  
  #get distance in y direction
  cur_shp_boundary_y<-cur_shp_boundary
  cur_shp_boundary_y[,1]<-ref_point[1]
  y_dist<-distHaversine(ref_point,cur_shp_boundary_y)*.001
  y_dist[cur_shp_boundary_y[,2]<ref_point[2]]<-y_dist[cur_shp_boundary_y[,2]<ref_point[2]]*-1
  
  gauge_loc<-which.min(sqrt((x_dist)^2+(y_dist)^2))
  z=complex(length.out = length(x_dist),real=x_dist[c(gauge_loc:length(x_dist),1:(gauge_loc-1))],imaginary=y_dist[c(gauge_loc:length(y_dist),1:(gauge_loc-1))])
  
  zfft=fft(z)
  
  for(j in 1:length(poss_numfor)){
    numFor<-poss_numfor[j]
    zfft_thresholded<-rep(0,length(zfft))
    to_keep<-c(seq(from=1,by=1,ceiling(numFor/2)),seq(from=length(zfft),by=-1,length.out = floor(numFor/2)))
    zfft_thresholded[to_keep]<-zfft[to_keep]
    estimate= fft(zfft_thresholded,inverse = T)/length(zfft_thresholded)
    error_max_norm[i,j]= max(abs(z- estimate))
    error_2_norm[i,j]= mean(abs(z- estimate)^2)
  }
}

plot(log(poss_numfor),log(apply(error_max_norm, 2, mean)),type="l")
plot(log(poss_numfor),log(apply(error_2_norm, 2, mean)),type="l")
plot(estimate)
plot(z)



