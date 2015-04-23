#market returns
mydata = read.csv("returns1.csv");

nasdaq=mydata[,4];
sandp=mydata[,5];
ftse=mydata[,6];
nasdaq=nasdaq[1:7284];
sandp=sandp[1:7284];
ftse=ftse[1:7284];
n=length(nasdaq);


#losses
x1=-100*nasdaq; u1=rank(x1)/n; x1d=density(x1); vp1=diff(sort(x1))*n;
x2=-100*sandp; u2=rank(x2)/n; x2d=density(x2); vp2=diff(sort(x2))*n;
x3=-100*ftse; u3=rank(x3)/n; x3d=density(x3); vp3=diff(sort(x3))*n;
xs=x1+x2+x3; us=rank(xs)/n; xsd=density(xs); vps=diff(sort(xs))*n;




#distribution plots
setEPS();
postscript("density.eps")
par(mar=c(5,5,1,1));
plot(x1d,col="red",xlab="x",ylab="f(x)",
xlim=c(-5,5),ylim=c(0,0.58),main=NA,cex=1.5,cex.lab=2)
lines(x2d,col="blue",cex=1.5,cex.lab=2)
lines(x3d,col="green",cex=1.5,cex.lab=2)
dev.off();


setEPS();
postscript("quantile.eps")
par(mar=c(5,5,1,1));
plot((1:n)/n,sort(x1),type="l",col="red",cex=1,ylim=c(-7,7),
xlab=expression(alpha),ylab=expression(V[alpha]),cex.lab=2);
lines((1:n)/n,sort(x2),lwd=1,col="blue")
lines((1:n)/n,sort(x3),,lwd=1,col="green")
dev.off();




#systematic risk densities
t=0.75;
phi1=(u1>=t)/(1-t);
phi2=(u2>=t)/(1-t);
phi3=(u3>=t)/(1-t);

r1=0; r2=0; r3=0;
m1=0; m2=0; m3=0;
for (i in 1:(n-1))
{
r1[i]=cov((u1>i/n),phi1)*vp1[i];
r2[i]=cov((u2>i/n),phi2)*vp2[i];
r3[i]=cov((u3>i/n),phi3)*vp3[i];

m1[i]=(1-i/n)*vp1[i];
m2[i]=(1-i/n)*vp2[i];
m3[i]=(1-i/n)*vp3[i];
}

r1z=movavg(r1,250); r2z=movavg(r2,250); r3z=movavg(r3,250); 
m1z=movavg(m1,250); m2z=movavg(m2,250); m3z=movavg(m3,250); 


setEPS();
postscript("mean.eps")
par(mar=c(5,5,1,1));
plot((1:length(m1z))/length(m1z),m1z,type="l",cex=1, ylim=c(0,4.2),
xlab=expression(alpha),ylab=expression(m[alpha]),cex.lab=2, col="red");
lines((1:length(m2z))/length(m2z),m2z,lty=1,lwd=1,col="blue");
lines((1:length(m3z))/length(m3z),m3z,lty=1,lwd=1,col="green");
dev.off();


setEPS();
postscript("risk.eps")
par(mar=c(5,5,1,1));
plot((1:length(r1z))/length(r1z),r1z,type="l",cex=1, ylim=c(0,4.2),
xlab=expression(alpha),ylab=expression(r[alpha]),cex.lab=2, col="red");
lines((1:length(r2z))/length(r2z),r2z,lty=1,lwd=1,col="blue");
lines((1:length(r3z))/length(r3z),r3z,lty=1,lwd=1,col="green");
dev.off();













#copula plots

setEPS();
postscript("nasdaqsp.eps")
par(mar=c(5,5,1,1));
plot(u1,u2,type="p",cex=0.01,
xlab=expression(u[1]),ylab=expression(u[2]),cex.lab=2);
dev.off();


setEPS();
postscript("spftse.eps")
par(mar=c(5,5,1,1));
plot(u2,u3,type="p",cex=0.01,
xlab=expression(u[2]),ylab=expression(u[3]),cex.lab=2);
dev.off();


setEPS();
postscript("nasdaqftse.eps")
par(mar=c(5,5,1,1));
plot(u1,u3,type="p",cex=0.01,
xlab=expression(u[1]),ylab=expression(u[3]),cex.lab=2);
dev.off();



#standalone risk densities
a=(1:(n-1))/n; t=0.75; Phi=(a>=t)*(a-t)/(1-t);
#r1=(a-Phi)*vp1; r2=(a-Phi)*vp2;

#systematic risk densities
rs=(a-Phi)*vps;
phis=(us>=t)/(1-t);
phi1=(u1>=t)/(1-t);
phi2=(u2>=t)/(1-t);
phi3=(u3>=t)/(1-t);

r1s=0; r2s=0; r3s=0; r1=0; r2=0; r3=0;
for (i in 1:(n-1))
{
r1s[i]=cov((u1>i/n),phis)*vp1[i]; 
r2s[i]=cov((u2>i/n),phis)*vp2[i];
r3s[i]=cov((u3>i/n),phis)*vp3[i];

r1[i]=cov((u1>i/n),phi1)*vp1[i];
r2[i]=cov((u2>i/n),phi2)*vp2[i];
r3[i]=cov((u3>i/n),phi3)*vp3[i];
}

r1z=movavg(r1,250); r1sz=movavg(r1s,250); 
r2z=movavg(r2,250); r2sz=movavg(r2s,250); 
r3z=movavg(r3,250); r3sz=movavg(r3s,250); 


setEPS();
postscript("nasdaqrisk.eps")
par(mar=c(5,5,1,1));
plot((1:length(r1z))/length(r1z),r1z,type="l",cex=1,
xlab=expression(alpha),ylab="Risk densities",cex.lab=2, col="red");
lines((1:length(r1z))/length(r1z),r1sz,lty=1,lwd=1,col="blue");
dev.off();


setEPS();
postscript("sprisk.eps")
par(mar=c(5,5,1,1));
plot((1:length(r2z))/length(r2z),r2z,type="l",cex=1,
xlab=expression(alpha),ylab="Risk densities",cex.lab=2, col="red");
lines((1:length(r2z))/length(r2z),r2sz,lty=1,lwd=1,col="blue");
dev.off();


setEPS();
postscript("ftserisk.eps")
par(mar=c(5,5,1,1));
plot((1:length(r3z))/length(r3z),r3z,type="l",cex=1,
xlab=expression(alpha),ylab="Risk densities",cex.lab=2, col="red");
lines((1:length(r3z))/length(r3z),r3sz,lty=1,lwd=1,col="blue");
dev.off();


setEPS();
postscript("risk1.eps")
par(mar=c(5,5,1,1));
plot((1:length(r1z))/length(r1z),r1z,type="l",cex=1, ylim=c(0,4.2),
xlab=expression(alpha),ylab="Standalone risk densities",cex.lab=2, col="red");
lines((1:length(r2z))/length(r2z),r2z,lty=1,lwd=1,col="blue");
lines((1:length(r3z))/length(r3z),r3z,lty=1,lwd=1,col="green");
dev.off();


setEPS();
postscript("risk2.eps")
par(mar=c(5,5,1,1));
plot((1:length(r1sz))/length(r1sz),r1sz,type="l",cex=1, ylim=c(0,4.2),
xlab=expression(alpha),ylab="Diversified risk densities",cex.lab=2, col="red");
lines((1:length(r2sz))/length(r2sz),r2sz,lty=1,lwd=1,col="blue");
lines((1:length(r3sz))/length(r3sz),r3sz,lty=1,lwd=1,col="green");
dev.off();




#overall risk
r1=cov(x1,phi1); r2=cov(x2,phi2); r3=cov(x3,phi3); rs=cov(xs,phis);
r1bar=cov(x1,phis); r2bar=cov(x2,phis); r3bar=cov(x3,phis); 

#overall risk after put option
q=0.95;
q1=quantile(x1,q); q2=quantile(x2,q); q3=quantile(x3,q); qs=quantile(xs,q); 
x1q=pmin(x1,q1); x2q=pmin(x2,q2); x3q=pmin(x3,q3); xsq=pmin(xs,qs); 
r1q=cov(x1q,phi1); r2q=cov(x2q,phi2); r3q=cov(x3q,phi3); rsq=cov(xsq,phis);
xsnew=x1q+x2q+x3q; usnew=rank(xsnew)/length(xsnew); phisnew=(usnew>=t)/(1-t);
r1qs=cov(x1q,phisnew); r2qs=cov(x2q,phisnew); r3qs=cov(x3q,phisnew);





#risk ratios
n=100; cutoff=matrix((1:n)/n);
ratio1=0; ratio2=0; ratio3=0;
for (i in 1:n)
{
k=cutoff[i]; ind1=(u1>k); ind2=(u2>k); ind3=(u3>k);
ratio1[i]=cov(ind1,phis)/cov(ind1,phi1);
ratio2[i]=cov(ind2,phis)/cov(ind2,phi2);
ratio3[i]=cov(ind3,phis)/cov(ind3,phi3);
}

setEPS();
postscript("nasdaqcop.eps")
par(mar=c(5,5,1,1));
plot(u1,us,type="p",cex=0.01,
xlab=expression(u[1]),ylab=expression(u["+"]),cex.lab=2);
lines(cutoff,ratio1,lty=1,lwd=5,col="red");
dev.off();

setEPS();
postscript("spcop.eps")
par(mar=c(5,5,1,1));
plot(u2,us,type="p",cex=0.01,
xlab=expression(u[2]),ylab=expression(u["+"]),cex.lab=2);
lines(cutoff,ratio2,lty=1,lwd=5,col="red");
dev.off();


setEPS();
postscript("ftsecop.eps")
par(mar=c(5,5,1,1));
plot(u3,us,type="p",cex=0.01,
xlab=expression(u[3]),ylab=expression(u["+"]),cex.lab=2);
lines(cutoff,ratio3,lty=1,lwd=5,col="red");
dev.off();


#aggregate densities
rs=(a-Phi)*vps; ms=(1-a)*vps;
rsz=movavg(rs,250); msz=movavg(ms,250); 


combined=cbind(xs,x1,x2,x3);
combineds=combined[order(xs),];
combinedd=diff(combineds);

xsp=combinedd[,1]; x1p=combinedd[,2]; x2p=combinedd[,3]; x3p=combinedd[,4];
xspz=movavg(xsp,250); x1pz=movavg(x1p,250); x2pz=movavg(x2p,250); x3pz=movavg(x3p,250);
p1z=x1pz/xspz; p2z=x2pz/xspz; p3z=x3pz/xspz; 
p1z=pmax(0,pmin(p1z,1)); p2z=pmax(0,pmin(p2z,1)); p3z=pmax(0,pmin(p3z,1)); tot=p1z+p2z+p3z;
p1z=p1z/psz; p2z=p2z/psz; p3z=p3z/psz; 
p1z=movavg(p1z,250); p2z=movavg(p2z,250); p3z=movavg(p3z,250);


t=length(rsz); s=length(p1z);
for (i in (s+1):t)
{
p1z[i]=(p1z[i-1]+p1z[i-2]+p1z[i-3])/3;
p2z[i]=(p2z[i-1]+p2z[i-2]+p2z[i-3])/3;
p3z[i]=(p3z[i-1]+p3z[i-2]+p3z[i-3])/3;
}

rsz1=rsz*p1z; msz1=msz*p1z; 
rsz2=rsz*p2z; msz2=msz*p2z;
rsz3=rsz*p3z; msz3=msz*p3z;





setEPS();
postscript("aggrisk.eps");
par(mar=c(5,5,1,1));
cutoff=(1:t)/t;
plot(cutoff,rsz1,type="l",xlab="Layer",ylab="Risk density",ylim=c(0,9),cex.lab=2,col="red");
lines(cutoff,rsz1+rsz2,col="blue")
lines(cutoff,rsz1+rsz2+rsz3,col="green")
dev.off();

setEPS();
postscript("aggmean.eps");
par(mar=c(5,5,1,1));
cutoff=(1:t)/t;
plot(cutoff,msz1,type="l",xlab="Layer",ylab="Mean density",ylim=c(0,9),cex.lab=2,col="red");
lines(cutoff,msz1+msz2,col="blue")
lines(cutoff,msz1+msz2+msz3,col="green")
dev.off();


setEPS();
postscript("aggrisk.eps");
par(mar=c(5,5,1,1));
cutoff=(1:t)/t;
plot(cutoff,rsz1,type="l",xlab="Layer",ylab="Risk density",ylim=c(0,9),cex.lab=2);
polygon(c(cutoff,rev(cutoff)),c(rsz1,rep(0,length(cutoff))),col="red");
polygon(c(cutoff,rev(cutoff)),c(rsz1+rsz2,rev(rsz1)),col="blue");
polygon(c(cutoff,rev(cutoff)),c(rsz1+rsz2+rsz3,rev(rsz1+rsz2)),col="green");
dev.off();


setEPS();
postscript("aggmean.eps");
par(mar=c(5,5,1,1));
cutoff=(1:t)/t;
plot(cutoff,msz1,type="l",xlab="Layer",ylab="Mean density",ylim=c(0,9),cex.lab=2);
polygon(c(cutoff,rev(cutoff)),c(msz1,rep(0,length(cutoff))),col="red");
polygon(c(cutoff,rev(cutoff)),c(msz1+msz2,rev(msz1)),col="blue");
polygon(c(cutoff,rev(cutoff)),c(msz1+msz2+msz3,rev(msz1+msz2)),col="green");
dev.off();






#moving average functions
movavg=function(series,lag)
{
smooth=0; z=length(series);

for (i in 1:lag)
{
smooth=smooth+series[(i):(z-lag+i)];
}
smooth=smooth/lag;

return(smooth)
} 
