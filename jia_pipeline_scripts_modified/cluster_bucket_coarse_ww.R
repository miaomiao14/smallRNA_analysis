plot_cluster_bucket <- function (input1,input2,numberofclusters) {
win=100;  ##???
name=paste(input1,'_',input2,'.pdf',sep="");
pdf(name,15,20,onefile=TRUE);  ##what does onefile parameter do?

for (i in 1:numberofclusters) {
#for (i in 1:15) {
file1=paste(input1,'.cluster',i,sep="");
file2=paste(input2,'.cluster',i,sep="");

if (file.exists(file1)==FALSE && file.exists(file2)==FALSE) { next;}

par(mfrow=c(4,3),cex=1,font.axis=2,font.lab=2,font.sub=2);

#lendis for reads and lendis for species
str1=paste(file1,'.lendis.reads',sep="");
str2=paste(file2,'.lendis.reads',sep="");

str3=paste(file1,'.lendis.species',sep="");
str4=paste(file2,'.lendis.species',sep="");

if(file.exists(str1)==TRUE && file.exists(str2)==TRUE)
{
L1=read.table(str1,header=FALSE);
M1=read.table(str2,header=FALSE);

ylim_plus=max(L1$V2,M1$V2)*1.2;ylim_minus=-max(L1$V3,M1$V3)*1.2;
reads=c("reads");
lendis.barplot(ylim_plus,ylim_minus,L1,reads);

L2=read.table(str3,header=FALSE);
M2=read.table(str4,header=FALSE);

ylim_plus=max(L2$V2,M2$V2)*1.2;ylim_minus=-max(L2$V3,M2$V3)*1.2;
species=c("species");
lendis.barplot(ylim_plus,ylim_minus,L2,species);

plot.new();

ylim_plus=max(L1$V2,M1$V2)*1.2;ylim_minus=-max(L1$V3,M1$V3)*1.2;
reads=c("reads");
lendis.barplot(ylim_plus,ylim_minus,M1,reads);

ylim_plus=max(L2$V2,M2$V2)*1.2;ylim_minus=-max(L2$V3,M2$V3)*1.2;
species=c("species");
lendis.barplot(ylim_plus,ylim_minus,M2,species);

plot.new();

}
else if (file.exists(str1)==TRUE)
{
    L1=read.table(str1,header=FALSE);
    ylim_plus=max(L1$V2)*1.2;ylim_minus=-max(L1$V3)*1.2;    
    reads=c("reads");
    lendis.barplot(ylim_plus,ylim_minus,L1,reads);
    
    L2=read.table(str3,header=FALSE);
    ylim_plus=max(L2$V2)*1.2;ylim_minus=-max(L2$V3)*1.2;
    species=c("species");
    lendis.barplot(ylim_plus,ylim_minus,L2,species);
    
    
    plot.new();
}
else if (file.exists(str2)==TRUE)
{
    M1=read.table(str2,header=FALSE);
    ylim_plus=max(M1$V2)*1.2;ylim_minus=-max(M1$V3)*1.2;
    reads=c("reads");
    lendis.barplot(ylim_plus,ylim_minus,M1,reads);

    M2=read.table(str4,header=FALSE);
    ylim_plus=max(M2$V2)*1.2;ylim_minus=-max(M2$V3)*1.2;
    species=c("species");
    lendis.barplot(ylim_plus,ylim_minus,M2,species);
    
    plot.new();    
}

if (file.exists(file1)==TRUE && file.exists(file2)==TRUE )
{ 

    #piRNA bucket 
    p=c();n=c();win=10;  ##step size for piRNA profile in each cluster
    s=c();z=c();
    
    F=read.table(file1,header=FALSE);
    G=read.table(file2,header=FALSE);
    
    pos=F$V2[F$V3=='+'];
    neg=F$V2[F$V3=='-'];
    
    G_pos=G$V2[G$V3=='+'];
    G_neg=G$V2[G$V3=='-'];
    
    for (i in 1:(max(F$V1)/win)){
    p[i]=sum(pos[(i*win-win+1):(i*win)]);
    n[i]=sum(neg[(i*win-win+1):(i*win)]);
    
    s[i]=sum(G_pos[(i*win-win+1):(i*win)]);
    z[i]=sum(G_neg[(i*win-win+1):(i*win)]);
    
    }

    x=(1:(max(F$V1)/win))*win-win+1;
    ylim_plus=max(p,s)*1.2;ylim_minus=-max(n,z)*1.2;
    plot(x,p,main='piRNA bucket',xlab='position',ylab='ppm',ylim=c(ylim_minus,ylim_plus),col=4,'l');
    lines(x,-n,col=2);
    
    #piRNA pp
    str=paste(file1,'.mapper2.pp',sep="");  
    R=read.table(str,header=FALSE);
    str=paste(file2,'.mapper2.pp',sep="");  
    R2=read.table(str,header=FALSE);
    
    str=paste(file1,'.pvalue',sep="");
    P=read.table(str,header=FALSE);
    str=paste(file2,'.pvalue',sep="");
    P2=read.table(str,header=FALSE);

    ylim_plus=max(R$V2,R2$V2)*1.2;ylim_minus=0;
    m=paste('piRNA pingpong  zscore=',sprintf("%.1f",P$V1),sep="");
    barplot(R$V2,names.arg=R$V1,main=m,col=4,ylab='pp raw score',ylim=c(ylim_minus,ylim_plus));

    #piRNA seqlogo
    str=paste(file1,'.seqlogo',sep="");
    S=read.table(str,header=FALSE);
    barplot(as.matrix(S),main='piRNA seqlogo',col=c("green","blue","yellow","red"));
#######################################################################################################
    #piRNA bucket
    ylim_plus=max(p,s)*1.2;ylim_minus=-max(n,z)*1.2;
    plot(x,s,main='piRNA bucket',xlab='position',ylab='ppm',ylim=c(ylim_minus,ylim_plus),col=4,'l');
    lines(x,-z,col=2);
    
    #piRNA pp
    ylim_plus=max(R$V2,R2$V2)*1.2;ylim_minus=0;
    m=paste('piRNA pingpong  zscore=',sprintf("%.1f",P2$V1),sep="");
    barplot(R2$V2,names.arg=R2$V1,main=m,col=4,ylab='pp raw score',ylim=c(ylim_minus,ylim_plus));
    
    #piRNA seqlogo
    str=paste(file2,'.seqlogo',sep="");
    Y=read.table(str,header=FALSE);
    barplot(as.matrix(Y),main='piRNA seqlogo',col=c("green","blue","yellow","red"));

}

else if (file.exists(file1)==TRUE )
{ 
#piRNA bucket 
p=c();n=c();win=10;  ##step size for piRNA profile in each cluster
F=read.table(file1,header=FALSE);
pos=F$V2[F$V3=='+']; neg=F$V2[F$V3=='-'];
for (i in 1:(max(F$V1)/win)){
p[i]=sum(pos[(i*win-win+1):(i*win)]);
n[i]=sum(neg[(i*win-win+1):(i*win)]);
}
x=(1:(max(F$V1)/win))*win-win+1;
ylim_plus=max(p)*1.2;ylim_minus=-max(n)*1.2;
plot(x,p,main='piRNA bucket',xlab='position',ylab='ppm',ylim=c(ylim_minus,ylim_plus),col=4,'l');
lines(x,-n,col=2);

#piRNA pp
str=paste(file1,'.mapper2.pp',sep="");  
R=read.table(str,header=FALSE);
str=paste(file1,'.pvalue',sep="");
P=read.table(str,header=FALSE);
ylim_plus=max(R$V2)*1.2;ylim_minus=0;
m=paste('piRNA pingpong  zscore=',sprintf("%.1f",P$V1),sep="");
barplot(R$V2,names.arg=R$V1,main=m,col=4,ylab='pp raw score',ylim=c(ylim_minus,ylim_plus));

#piRNA seqlogo
str=paste(file1,'.seqlogo',sep="");
S=read.table(str,header=FALSE);
barplot(as.matrix(S),main='piRNA seqlogo',col=c("green","blue","yellow","red"));


    plot.new();
    plot.new();
    plot.new();

}

else if (file.exists(file2)==TRUE )
{
    plot.new();
    plot.new();
    plot.new();

#piRNA bucket 
p=c();n=c();win=10;  ##step size for piRNA profile in each cluster
F=read.table(file2,header=FALSE);
pos=F$V2[F$V3=='+']; neg=F$V2[F$V3=='-'];
for (i in 1:(max(F$V1)/win)){
p[i]=sum(pos[(i*win-win+1):(i*win)]);
n[i]=sum(neg[(i*win-win+1):(i*win)]);
}
x=(1:(max(F$V1)/win))*win-win+1;
ylim_plus=max(p)*1.2;ylim_minus=-max(n)*1.2;
plot(x,p,main='piRNA bucket',xlab='position',ylab='ppm',ylim=c(ylim_minus,ylim_plus),col=4,'l');
lines(x,-n,col=2);

#piRNA pp
str=paste(file2,'.mapper2.pp',sep="");  
R=read.table(str,header=FALSE);
str=paste(file2,'.pvalue',sep="");
P=read.table(str,header=FALSE);
ylim_plus=max(R$V2)*1.2;ylim_minus=0;
m=paste('piRNA pingpong  zscore=',sprintf("%.1f",P$V1),sep="");
barplot(R$V2,names.arg=R$V1,main=m,col=4,ylab='pp raw score',ylim=c(ylim_minus,ylim_plus));

#piRNA seqlogo
str=paste(file2,'.seqlogo',sep="");
S=read.table(str,header=FALSE);
barplot(as.matrix(S),main='piRNA seqlogo',col=c("green","blue","yellow","red"));

}


#loc
str=paste(file1,'.loc',sep="");

if (file.exists(str)==TRUE )
{
    L=read.table(str,header=FALSE);
}
else
{
    str=paste(file2,'.loc',sep="");
    L=read.table(str,header=FALSE);
}
fi1=strsplit(file1, "\\.")[[1]][3];
cl1=strsplit(file1, "\\.")[[1]][11];
fi2=strsplit(file2, "\\.")[[1]][3];
cl2=strsplit(file2, "\\.")[[1]][11];
t=paste(fi1,'_',cl1,fi2,'_',cl2,' ',L$V1,sep="");
mtext(t,side=3,line=-1,outer=TRUE,cex=1.2);
}
dev.off();
}

lendis.barplot <- function (ymax,ymin,f,reads)
{
    L=f;ylim_plus=ymax;ylim_minus=ymin;
    sense_fration=sum(L$V2[L$V1>=23 & L$V1<=29])/(sum(L$V2[L$V1>=23 & L$V1<=29])+sum(L$V3[L$V1>=23 & L$V1<=29]));
    s=paste("plus fraction=",sprintf("%.2f",sense_fration),sep="");
    m=paste(reads," lendis",sep="");
    lab=paste('number of normalized ',reads,sep="");
    barplot(L$V2,main=m,col=4,ylab=lab,ylim=c(ylim_minus,ylim_plus),sub=s);
    barplot(-L$V3,names.arg=L$V1,col=2,xlab='length (nt)',ylim=c(ylim_minus,ylim_plus),add=TRUE);
}




