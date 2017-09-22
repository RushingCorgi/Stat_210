#input the dataset#
rm(list=ls());
mydata <- read.csv(file="winequality-red.csv",head=TRUE,sep=";");

#plot the quality of sample wine#
Y<-mydata$quality;
sample.num<-length(Y);
N<-seq(1,by=1,length=sample.num);
plot(N,Y,main='The quality of the all the samples',xlab='serial number',ylab='wine quality level');

#compute the correlation matrix#
#use abbreviation for data header#
X.old<-mydata[,1:11];
names(X.old);
names(X.old)[names(X.old)=='fixed.acidity']='f.a.';
names(X.old)[names(X.old)=='volatile.acidity']='v.a.';
names(X.old)[names(X.old)=='citric.acid']='c.a.';
names(X.old)[names(X.old)=='residual.sugar']='r.s.';
names(X.old)[names(X.old)=='chlorides']='chlo.';
names(X.old)[names(X.old)=='free.sulfur.dioxide']='f.s.d.';
names(X.old)[names(X.old)=='total.sulfur.dioxide']='t.s.d';
names(X.old)[names(X.old)=='density']='dens.';
names(X.old)[names(X.old)=='sulphates']='sul.';
names(X.old)[names(X.old)=='alcohol']='alco.';
cor.matrix<-cor(X.old);
cor.matrix<-round(cor.matrix,2);
library(xtable);
xtable(cor.matrix);
# print(xtable(cor.matrix,
#              caption="String",
#              label="t:"),
#       type="latex",
#       file="cor.matrix",
#       table.placement="tp",
#       latex.environments=c("center", "footnotesize"))

#initial multiple linear regression#
X1<-X.old[,1];
X2<-X.old[,2];
X3<-X.old[,3];
X4<-X.old[,4];
X5<-X.old[,5];
X6<-X.old[,6];
X7<-X.old[,7];
X8<-X.old[,8];
X9<-X.old[,9];
X10<-X.old[,10];
X11<-X.old[,11];
fit0<-lm(Y~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11);
plot(fit0$residuals,main='The residuals of the initial fit',xlab='serial number',ylab='residuals');

#dichotomize "density"
plot(X8,main="The scatterplot of density",xlab="serial number",ylab="density");
abline(h=0.997,col='red',lwd=4);
D<-matrix(X8,nrow=1599);
count<-0;
for (val in X8) {
  count=count+1;
  if(val >=0.997) D[count,1]<-1
  else            D[count,1]<-0
}
#Set new dataset#
X.new=matrix(c(X3,D,X4,X5,X6,X9,X10,X11),ncol=8,nrow=1599);
X1<-X.new[,1];
X2<-X.new[,2];
X3<-X.new[,3];
X4<-X.new[,4];
X5<-X.new[,5];
X6<-X.new[,6];
X7<-X.new[,7];
X8<-X.new[,8];

#decimal base to binary base#
ten2two<-function(x){
  binary.x=matrix(0,8,1)
  i=8
  while(x>0&&i>0){
    q=floor(x/2)
    r=x-q*2
    binary.x[i]=r
    i=i-1
    x=q
  }
  return(binary.x)
}
#Compute all 256 AIC and BIC#
serial.matrix=matrix(0,256,8);
AIC=matrix(0,256,1);
BIC=matrix(0,256,1);
for(i in 1:256){
  i=i-1;
  c=ten2two(i);
  serial.matrix[i+1,]=c
  
  X=matrix(1,1599,1);
  count=0;
  for(k in 1:8){
    if(c[k]==1){
      X=cbind(X,X.new[,k]);
      count=count+1;
    }
  }
  n=1599;
  XtX = t(X)%*%X;
  XtXinv = solve(XtX);
  XtY = t(X)%*%Y;
  betahat = XtXinv%*%XtY;
  Yhat = X%*%betahat;
  R = Y - Yhat;
  SSE = sum(R^2);
  df = n-(count+1);
  sigmahat2=SSE/df;
  
  AIC[i+1,1]=log(sigmahat2)+(n+count+1)/(n-count-3);
  BIC[i+1,1]=log(sigmahat2)+(count+1)*log(n)/n;
}

#find the best model by AIC#
best.num1=which.min(AIC);
best.AIC=min(AIC);
c=serial.matrix[best.num1,];

#find the best model by BIC#
best.num2=which.min(BIC);
best.BIC=min(BIC);
c=serial.matrix[best.num2,];

#compute in the best model above#
X=matrix(1,1599,1);
for(k in 1:8){
  if(c[k]==1){
    X=cbind(X,X.new[,k]);
  }
}
XtX = t(X)%*%X;
XtXinv = solve(XtX);
XtY = t(X)%*%Y;
betahat = XtXinv%*%XtY;
Yhat = X%*%betahat;
R = Y - Yhat;

#Forward method#
X.forw=X.new;
X0=matrix(1,1599,1);
count=1;
temp=0;
max.iter=8;
j=1;
while(j<9){
  X0tX0 = t(X0)%*%X0;
  X0tX0inv = solve(X0tX0);
  X0tY = t(X0)%*%Y;
  gamma0hat = X0tX0inv%*%X0tY;
  Y0hat =  X0%*%gamma0hat;
  R0 = Y - Y0hat;
  SSE0 = sum(R0^2);
  df0 = n-count;
  i=0;
  while (i<max.iter){
    i=i+1;
    X=cbind(X0,X.forw[,i]);
    XtX = t(X)%*%X;
    XtXinv = solve(XtX);
    XtY = t(X)%*%Y;
    betahat = XtXinv%*%XtY;
    Yhat = X%*%betahat;
    R = Y - Yhat;
    SSE1 = sum(R^2);
    df1 = n - (count+1);
    sigmahat2=SSE1/df1;
    
    fobs = ((SSE0 - SSE1)/(df0 - df1))/(SSE1/df1);
    tau = qf(0.99, (df0-df1), df1)
    if(fobs>tau){
      temp=1;
      break;
    }
  }
  if(temp==0 && i==max.iter){
    break;
  }
  if(max.iter==1){
    break;
  }
  X0=X;
  if(max.iter==2){
    X.forw<-X.forw[,-i];
    X.forw<-matrix(X.forw,nrow=1599,ncol=1)
  }
  else{
    X.forw<-X.forw[,-i]; 
  }
  count=count+1;
  max.iter=max.iter-1;
  temp=0;
  j=j+1;
}


#backward method#
X0=matrix(1,1599,1);
X.back=cbind(X0,X.new);

count=8;
temp=0;
max.iter=9;
j=1;
while(j<9){
  X=X.back;
  XtX = t(X)%*%X;
  XtXinv = solve(XtX);
  XtY = t(X)%*%Y;
  betahat = XtXinv%*%XtY;
  Yhat = X%*%betahat;
  R = Y - Yhat;
  SSE1 = sum(R^2);
  df1 = n - (count+1);
  sigmahat2=SSE1/df1;
  
  i=1;
  while (i<max.iter){
    i=i+1;
    X0=X.back[,-i];
    X0tX0 = t(X0)%*%X0;
    X0tX0inv = solve(X0tX0);
    X0tY = t(X0)%*%Y;
    gamma0hat = X0tX0inv%*%X0tY;
    Y0hat =  X0%*%gamma0hat;
    R0 = Y - Y0hat;
    SSE0 = sum(R0^2);
    df0 = n-count;
    
    fobs = ((SSE0 - SSE1)/(df0 - df1))/(SSE1/df1);
    tau = qf(0.99, (df0-df1), df1)
    if(fobs<tau){
      temp=1;
      break;
    }
  }
  X.back=X0;
  if(temp==0 && i==max.iter){
    break;
  }
  if(max.iter==2){
    break;
  }
  count=count-1;
  max.iter=max.iter-1;
  temp=0;
  j=j+1;
}

#two way interaction#
X1=X.new[,1];
X2=X.new[,4];
X3=X.new[,6];
X4=X.new[,7];
X5=X.new[,8];
#1st case#
n=1599;
tw=X1*X2;
fit12=lm(Y~X1+X2+tw);
R1=fit12$residuals;
SSE1 = sum(R1^2);
df1 = n-4;
fit12.0=lm(Y~X1+X2);
R0=fit12.0$residuals;
SSE0=sum(R0^2);
df0= n-3;
fobs = ((SSE0 - SSE1)/(df0 - df1))/(SSE1/df1);
tau = qf(0.95, (df0-df1), df1)
if(fobs>tau){
  print("add beta12")
}

#2nd case#
n=1599;
tw=X1*X3;
fit13=lm(Y~X1+X3+tw);
R1=fit13$residuals;
SSE1 = sum(R1^2);
df1 = n-4;
fit13.0=lm(Y~X1+X3);
R0=fit13.0$residuals;
SSE0=sum(R0^2);
df0= n-3;
fobs = ((SSE0 - SSE1)/(df0 - df1))/(SSE1/df1);
tau = qf(0.95, (df0-df1), df1)
if(fobs>tau){
  print("add beta13")
}

#3rd case#
n=1599;
tw=X1*X4;
fit14=lm(Y~X1+X4+tw);
R1=fit14$residuals;
SSE1 = sum(R1^2);
df1 = n-4;
fit14.0=lm(Y~X1+X4);
R0=fit14.0$residuals;
SSE0=sum(R0^2);
df0= n-3;
fobs = ((SSE0 - SSE1)/(df0 - df1))/(SSE1/df1);
tau = qf(0.95, (df0-df1), df1)
if(fobs>tau){
  print("add beta14")
}

#4th case#
n=1599;
tw=X1*X5;
fit15=lm(Y~X1+X5+tw);
R1=fit15$residuals;
SSE1 = sum(R1^2);
df1 = n-4;
fit15.0=lm(Y~X1+X5);
R0=fit15.0$residuals;
SSE0=sum(R0^2);
df0= n-3;
fobs = ((SSE0 - SSE1)/(df0 - df1))/(SSE1/df1);
tau = qf(0.95, (df0-df1), df1)
if(fobs>tau){
  print("add beta15")
}

#5th case#
n=1599;
tw=X2*X3;
fit23=lm(Y~X2+X3+tw);
R1=fit23$residuals;
SSE1 = sum(R1^2);
df1 = n-4;
fit23.0=lm(Y~X2+X3);
R0=fit23.0$residuals;
SSE0=sum(R0^2);
df0= n-3;
fobs = ((SSE0 - SSE1)/(df0 - df1))/(SSE1/df1);
tau = qf(0.95, (df0-df1), df1)
if(fobs>tau){
  print("add beta23")
}

#6th case#
n=1599;
tw=X2*X4;
fit24=lm(Y~X2+X4+tw);
R1=fit24$residuals;
SSE1 = sum(R1^2);
df1 = n-4;
fit24.0=lm(Y~X2+X4);
R0=fit24.0$residuals;
SSE0=sum(R0^2);
df0= n-3;
fobs = ((SSE0 - SSE1)/(df0 - df1))/(SSE1/df1);
tau = qf(0.95, (df0-df1), df1)
if(fobs>tau){
  print("add beta24")
}

#7th case#
n=1599;
tw=X2*X5;
fit25=lm(Y~X2+X5+tw);
R1=fit25$residuals;
SSE1 = sum(R1^2);
df1 = n-4;
fit25.0=lm(Y~X2+X5);
R0=fit25.0$residuals;
SSE0=sum(R0^2);
df0= n-3;
fobs = ((SSE0 - SSE1)/(df0 - df1))/(SSE1/df1);
tau = qf(0.95, (df0-df1), df1)
if(fobs>tau){
  print("add beta25")
}

#8th case#
n=1599;
tw=X3*X4;
fit34=lm(Y~X3+X4+tw);
R1=fit34$residuals;
SSE1 = sum(R1^2);
df1 = n-4;
fit34.0=lm(Y~X3+X4);
R0=fit34.0$residuals;
SSE0=sum(R0^2);
df0= n-3;
fobs = ((SSE0 - SSE1)/(df0 - df1))/(SSE1/df1);
tau = qf(0.95, (df0-df1), df1)
if(fobs>tau){
  print("add beta34")
}

#9th case#
n=1599;
tw=X3*X5;
fit35=lm(Y~X3+X5+tw);
R1=fit35$residuals;
SSE1 = sum(R1^2);
df1 = n-4;
fit35.0=lm(Y~X3+X5);
R0=fit35.0$residuals;
SSE0=sum(R0^2);
df0= n-3;
fobs = ((SSE0 - SSE1)/(df0 - df1))/(SSE1/df1);
tau = qf(0.95, (df0-df1), df1)
if(fobs>tau){
  print("add beta35")
}

#10th case#
n=1599;
tw=X4*X5;
fit45=lm(Y~X4+X5+tw);
R1=fit45$residuals;
SSE1 = sum(R1^2);
df1 = n-4;
fit45.0=lm(Y~X4+X5);
R0=fit45.0$residuals;
SSE0=sum(R0^2);
df0= n-3;
fobs = ((SSE0 - SSE1)/(df0 - df1))/(SSE1/df1);
tau = qf(0.95, (df0-df1), df1)
if(fobs>tau){
  print("add beta45")
}

#full model with all the two way terms#
X14=X1*X4;
X23=X2*X3;
X24=X2*X4;
X34=X3*X4;
X35=X3*X5;
X45=X4*X5;
fit.twoway=lm(Y~X1+X2+X3+X4+X5+X14+X23+X24+X34+X35+X45);
fit.reduced=lm(Y~X1+X2+X3+X4+X5);
anova(fit.twoway,fit.reduced)
