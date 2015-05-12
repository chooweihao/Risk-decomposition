a=(0:999)/1000;

#uniform
m_uni=2*(1-a);

#exponential
m_exp=a/a;

#pareto (mean=beta/(gamma-1))
beta=0.5; gamma=1.5;
m_par=beta/gamma/(1-a)^(1/gamma);

#weibull (mean=beta*Gamma(1+1/k))
k=2; beta=2/(pi^0.5);
m_wei=beta/k*(-log(1-a))^(1/k-1);


setEPS();
postscript("meancont.eps")
plot(a,m_uni,type="l",
xlab=expression(alpha),ylab=expression(m[alpha]),ylim=c(0,3),col="green");
lines(a,m_exp,col="black");
lines(a,m_par,col="blue");
lines(a,m_wei,col="red");
legend(0.2,3,c("Uniform","Exponential","Pareto","Weibull"),
lty=c(1,1,1,1), col=c("green","black","blue","red"));
dev.off();


#mean and risk densities

rr=(a-a^3)/(1-a); #risk ratio

setEPS();
postscript("unirisk.eps")
par(mar=c(5,5,1,1));
plot(a,m_uni,type="l",
xlab=expression(alpha),
ylab=expression(paste(list(m[alpha],r[alpha]))),cex.lab=2,ylim=c(0,5));
title("Uniform", line = -2, cex.main=2);
polygon(c(a,rev(a)),c(m_uni,rep(0,length(a))),col="blue");
polygon(c(a,rev(a)),c(m_uni*(1+rr),rev(m_uni)),col="red");
dev.off();

setEPS();
postscript("exprisk.eps")
par(mar=c(5,5,1,1));
plot(a,m_exp,type="l",
xlab=expression(alpha),
ylab=expression(paste(list(m[alpha],r[alpha]))),cex.lab=2,ylim=c(0,5));
title("Exponential", line = -2, cex.main=2);
polygon(c(a,rev(a)),c(m_exp,rep(0,length(a))),col="blue");
polygon(c(a,rev(a)),c(m_exp*(1+rr),rev(m_exp)),col="red");
dev.off();

setEPS();
postscript("parrisk.eps")
par(mar=c(5,5,1,1));
plot(a,m_par, type="l",
xlab=expression(alpha),
ylab=expression(paste(list(m[alpha],r[alpha]))),cex.lab=2,ylim=c(0,5));
title("Pareto", line = -2, cex.main=2);
polygon(c(a,rev(a)),c(m_par,rep(0,length(a))),col="blue");
polygon(c(a,rev(a)),c(m_par*(1+rr),rev(m_par)),col="red");
dev.off();

setEPS();
postscript("weirisk.eps")
par(mar=c(5,5,1,1));
plot(a,m_wei,type="l",
xlab=expression(alpha),
ylab=expression(paste(list(m[alpha],r[alpha]))),cex.lab=2,ylim=c(0,5));
title("Weibull", line = -2, cex.main=2);
polygon(c(a,rev(a)),c(m_wei,rep(0,length(a))),col="blue");
polygon(c(a,rev(a)),c(m_wei*(1+rr),rev(m_wei)),col="red");
dev.off();





#comparison of VaR and ES capital buffers

#weibull
a=(0:999)/1000;
k=2; beta=2/(pi^0.5);
m=beta/k*(-log(1-a))^(1/k-1);
m2=(a>0.85)*m;


setEPS();
postscript("shortfallcap.eps");
plot(a,m,type="l",xlab=expression(alpha),ylab=expression(m[alpha]),ylim=c(0,3),cex.lab=2);
polygon(c(a,rev(a)),c(m2,rep(0,length(a))),col="red");
dev.off();



x=(1:1000)/1000*4;
fx=k*x^(k-1)/(beta^k)*exp(-(x/beta)^k);
f2=(x>2.2)*fx;

setEPS();
postscript("varcap.eps");
plot(x,fx,type="l",xlab="x",ylab="f(x)",cex.lab=2);
polygon(c(x,rev(x)),c(f2,rep(0,length(x))),col="red");
dev.off();


