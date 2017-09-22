datahere=matrix(ncol=4,scan('datamt'),byrow=T);
n=100;
N=datahere[,1];
x=datahere[,2];
U0=datahere[,3];
U1=datahere[,4];
y=U1-U0;
N1=N[1:50];
N2=N[51:100];
x1=x[1:50];
x2=x[51:100];
y1=y[1:50];
y2=y[51:100];
old.par<-par(mfrow=c(1,2))
plot(N1,x1,main='group 1',ylab="amydala x_i",xlab="i from 1 to 50")
plot(N2,x2,main='group 2',ylab="amydala x_i",xlab="i from 51 to 100")
par(old.par)

old.par<-par(mfrow=c(1,2))
plot(N1,y1,main='group 1',ylab="improved scorse",xlab="i from 1 to 50")
plot(N2,y2,main='group 2',ylab="improved scorse",xlab="i from 51 to 100")
par(old.par)

X1=c(rep(1,length=50),rep(0,length=50));
X2=c(x1,rep(0,length=50));
X3=c(rep(0,length=50),rep(1,length=50));
X4=c(rep(0,length=50),x2);
X=cbind(X1,X2,X3,X4);
XtX = t(X)%*%X;
XtXinv = solve(XtX);
XtY = t(X)%*%y;
betahat = XtXinv%*%XtY;
Yhat = X%*%betahat;
R = y - Yhat;
SSE1 = sum(R^2);
df1 = n - 4;
sigmahat2=SSE1/df1;

X0=cbind(rep(1,length=100),x);
X0tX0 = t(X0)%*%X0;
X0tX0inv = solve(X0tX0);
X0tY = t(X0)%*%y;
gamma0hat = X0tX0inv%*%X0tY;
Y0hat =  X0%*%gamma0hat;
R0 = y - Y0hat;
SSE0 = (t(R0))%*%R0;
df0 = n-2;

fobs = ((SSE0 - SSE1)/(df0 - df1))/(SSE1/df1)
fobs

tau = qf(0.95, (df0-df1), df1)
pvalue = 1 - pf(fobs, (df0-df1), df1)

c0<-c(1,0,-1,0);
c0=t(c0);
tc0=t(invc0);
c1<-c(0,1,0,-1);
c1=t(c1);
tc1=t(tc1);

tau2=qt(0.975,df1);
t0obs=c0%*%betahat/(sqrt(c0%*%XtXinv%*%tc0*sigmahat2));
t1obs=c1%*%betahat/(sqrt(c1%*%XtXinv%*%tc1*sigmahat2));