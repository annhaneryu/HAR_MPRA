############# HepG2 confirmation of data

setwd("/Users/haneryu/Dropbox/Hane_HARMPRAdata_Sean/data")
load("/Users/haneryu/Dropbox/Hane_HARMPRAdata_Sean/workspace.RData")

getColGrad <- function(a,b,x) {
	rgb(
		round((b[1]-a[1])*c(1:x)/x+a[1],2),
		round((b[2]-a[2])*c(1:x)/x+a[2],2),
		round((b[3]-a[3])*c(1:x)/x+a[3],2)
	)
}
whiteColorVector <- c(255,255,255)/255
redColorVector <- c(255,0,0)/255
cvWR <- getColGrad(whiteColorVector,redColorVector,100)

getHARdata <- function (path,d1f,d2f,d3f,r1f,r2f,r3f) {
	d1 <- read.delim(paste(path,d1f,sep="/"),header=F)
	d2 <- read.delim(paste(path,d2f,sep="/"),header=F)
	d3 <- read.delim(paste(path,d3f,sep="/"),header=F)
	r1 <- read.delim(paste(path,r1f,sep="/"),header=F)
	r2 <- read.delim(paste(path,r2f,sep="/"),header=F)
	r3 <- read.delim(paste(path,r3f,sep="/"),header=F)
	data.1 <- merge(r1,d1,by="V1")[,-3]
	data.2 <- merge(r2,d2,by="V1")[,-3]
	data.3 <- merge(r3,d3,by="V1")[,-3]
	colnames(data.1) <- c("barcodes","X","Y","name")
	colnames(data.2) <- c("barcodes","X","Y","name")
	colnames(data.3) <- c("barcodes","X","Y","name")
	data.1$category <- sapply(strsplit(as.character(data.1$name),":",fixed=T),FUN=function(x){ paste(x[1:(length(x)-1)],sep="",collapse=":") })
	data.1$ratio <- (data.1$X/(sum(data.1$X)*10^6))/(data.1$Y/(sum(data.1$Y)*10^6))
	data.2$category <- sapply(strsplit(as.character(data.2$name),":",fixed=T),FUN=function(x){ paste(x[1:(length(x)-1)],sep="",collapse=":") })
	data.2$ratio <- (data.2$X/(sum(data.2$X)*10^6))/(data.2$Y/(sum(data.2$Y)*10^6))
	data.3$category <- sapply(strsplit(as.character(data.3$name),":",fixed=T),FUN=function(x){ paste(x[1:(length(x)-1)],sep="",collapse=":") })
	data.3$ratio <- (data.3$X/(sum(data.3$X)*10^6))/(data.3$Y/(sum(data.3$Y)*10^6))
	res.1 <- as.data.frame(t(sapply(unique(data.1$category),FUN=function(x) { sel <- which(data.1$category == x); c((sum(data.1$X[sel])/length(sel))/sum(data.1$X)*10^6,(sum(data.1$Y[sel])/length(sel))/sum(data.1$Y)*10^6,length(sel)) } )))
	res.2 <- as.data.frame(t(sapply(unique(data.2$category),FUN=function(x) { sel <- which(data.2$category == x); c((sum(data.2$X[sel])/length(sel))/sum(data.2$X)*10^6,(sum(data.2$Y[sel])/length(sel))/sum(data.2$Y)*10^6,length(sel)) } )))
	res.3 <- as.data.frame(t(sapply(unique(data.3$category),FUN=function(x) { sel <- which(data.3$category == x); c((sum(data.3$X[sel])/length(sel))/sum(data.3$X)*10^6,(sum(data.3$Y[sel])/length(sel))/sum(data.3$Y)*10^6,length(sel)) } )))
	res.1$name <- rownames(res.1)
	res.2$name <- rownames(res.2)
	res.3$name <- rownames(res.3)
	res.1$ratio <- res.1$V1/res.1$V2
	res.2$ratio <- res.2$V1/res.2$V2
	res.3$ratio <- res.3$V1/res.3$V2
	v <- c()
	v$res1 <- res.1
	v$res2 <- res.2
	v$res3 <- res.3
	return(v)
}

gc1.l1 <- getHARdata(getwd(),"GC1-1-DNA.tsv","GC1-2-DNA.tsv","GC1-3-DNA.tsv","GC1-1-RNA.tsv","GC1-2-RNA.tsv","GC1-3-RNA.tsv")
wtc.gpc.l1 <- getHARdata(getwd(),"WTc-1-DNA.tsv","WTc-2-DNA.tsv","WTc-3-DNA.tsv","WTc-1-RNA.tsv","WTc-2-RNA.tsv","WTc-3-RNA.tsv")
wtc.gpc.l2 <- getHARdata(getwd(),"WTC-GPC-1-DNA.tsv","WTC-GPC-2-DNA.tsv","WTC-GPC-3-DNA.tsv","WTC-GPC-1-RNA.tsv","WTC-GPC-2-RNA.tsv","WTC-GPC-3-RNA.tsv")
pt2a.gpc.l1 <- getHARdata(getwd(),"Pt2A-1-DNA.tsv","Pt2A-2-DNA.tsv","Pt2A-3-DNA.tsv","Pt2A-1-RNA.tsv","Pt2A-2-RNA.tsv","Pt2A-3-RNA.tsv")
pt2a.gpc.l2 <- getHARdata(getwd(),"Pt2A-GPC-1-DNA.tsv","Pt2A-GPC-2-DNA.tsv","Pt2A-GPC-3-DNA.tsv","Pt2A-GPC-1-RNA.tsv","Pt2A-GPC-2-RNA.tsv","Pt2A-GPC-3-RNA.tsv")
pt5c.gpc.l1 <- getHARdata(getwd(),"Pt5C-GPC-1-DNA.tsv","Pt5C-GPC-2-DNA.tsv","Pt5C-GPC-3-DNA.tsv","Pt5C-GPC-1-RNA.tsv","Pt5C-GPC-2-RNA.tsv","Pt5C-GPC-3-RNA.tsv")
pt5c.npc.l1 <- getHARdata(getwd(),"Pt5C-NPC-1-DNA.tsv","Pt5C-NPC-2-DNA.tsv","Pt5C-NPC-3-DNA.tsv","Pt5C-NPC-1-RNA.tsv","Pt5C-NPC-2-RNA.tsv","Pt5C-NPC-3-RNA.tsv")
hs1.gpc.l1 <- getHARdata(getwd(),"HS1-11-GPC-1-DNA.tsv","HS1-11-GPC-2-DNA.tsv","HS1-11-GPC-3-DNA.tsv","HS1-11-GPC-1-RNA.tsv","HS1-11-GPC-2-RNA.tsv","HS1-11-GPC-3-RNA.tsv")
hs1.npc.l1 <- getHARdata(getwd(),"HS1-11-NPC-1-DNA.tsv","HS1-11-NPC-2-DNA.tsv","HS1-11-NPC-3-DNA.tsv","HS1-11-NPC-1-RNA.tsv","HS1-11-NPC-2-RNA.tsv","HS1-11-NPC-3-RNA.tsv")

########
# get full matrix of ratios
{
uid <- data.frame(sort(unique(c(wtc.gpc.l1$res1$name,wtc.gpc.l1$res2$name,wtc.gpc.l1$res3$name,wtc.gpc.l2$res1$name,wtc.gpc.l2$res2$name,wtc.gpc.l2$res3$name,pt2a.gpc.l1$res1$name,pt2a.gpc.l1$res2$name,pt2a.gpc.l1$res3$name,pt2a.gpc.l2$res1$name,pt2a.gpc.l2$res2$name,pt2a.gpc.l2$res3$name,gc1.l1$res1$name,gc1.l1$res2$name,gc1.l1$res3$name,hs1.gpc.l1$res1$name,hs1.gpc.l1$res2$name,hs1.gpc.l1$res3$name,pt5c.npc.l1$res1$name,pt5c.npc.l1$res2$name,pt5c.npc.l1$res3$name))))
uid$pos <- 0
uid$neg <- 0
colnames(uid)[1] <- "name"
uid$pos[grep("POSITIVE_CONTROL",uid$name)] <- 1
uid$neg[grep("NEGATIVE_CONTROL",uid$name)] <- 1
vista.pos <- "brain_vista_pos_ctrl10__cut=1_of_1__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED"
encode.neg <- "ENCODE_neg_ctrl42__cut=1_of_1__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED"
encode.pos <- "ENCODE_pos_ctrl5__cut=1_of_1__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED"
k27me3.neg <- "N1_K27me3_neg_ctrl42__cut=1_of_2__length=171__range_inclusive=1_to_171__VALIDATED_STRONG_CONTROL__no_modifications__NO_EcoRI_RESTRICTION_SITES_MODIFIED__NO_SbfI_RESTRICTION_SITES_MODIFIED"

x <-  match(uid$name,wtc.gpc.l1$res1$name)
w <- which(!is.na(x))
uid$wtc.gpc.l1.r1 <- rep(NA,dim(uid)[1])
uid$wtc.gpc.l1.r1[w] <- wtc.gpc.l1$res1$ratio[x[w]]
x <-  match(uid$name,wtc.gpc.l1$res2$name)
w <- which(!is.na(x))
uid$wtc.gpc.l1.r2 <- rep(NA,dim(uid)[1])
uid$wtc.gpc.l1.r2[w] <- wtc.gpc.l1$res2$ratio[x[w]]
x <-  match(uid$name,wtc.gpc.l1$res3$name)
w <- which(!is.na(x))
uid$wtc.gpc.l1.r3 <- rep(NA,dim(uid)[1])
uid$wtc.gpc.l1.r3[w] <- wtc.gpc.l1$res3$ratio[x[w]]

x <-  match(uid$name,wtc.gpc.l2$res1$name)
w <- which(!is.na(x))
uid$wtc.gpc.l2.r1 <- rep(NA,dim(uid)[1])
uid$wtc.gpc.l2.r1[w] <- wtc.gpc.l2$res1$ratio[x[w]]
x <-  match(uid$name,wtc.gpc.l2$res2$name)
w <- which(!is.na(x))
uid$wtc.gpc.l2.r2 <- rep(NA,dim(uid)[1])
uid$wtc.gpc.l2.r2[w] <- wtc.gpc.l2$res2$ratio[x[w]]
x <-  match(uid$name,wtc.gpc.l2$res3$name)
w <- which(!is.na(x))
uid$wtc.gpc.l2.r3 <- rep(NA,dim(uid)[1])
uid$wtc.gpc.l2.r3[w] <- wtc.gpc.l2$res3$ratio[x[w]]

x <-  match(uid$name,pt2a.gpc.l1$res1$name)
w <- which(!is.na(x))
uid$pt2.gpc.l1.r1 <- rep(NA,dim(uid)[1])
uid$pt2.gpc.l1.r1[w] <- pt2a.gpc.l1$res1$ratio[x[w]]
x <-  match(uid$name,pt2a.gpc.l1$res2$name)
w <- which(!is.na(x))
uid$pt2.gpc.l1.r2 <- rep(NA,dim(uid)[1])
uid$pt2.gpc.l1.r2[w] <- pt2a.gpc.l1$res2$ratio[x[w]]
x <-  match(uid$name,pt2a.gpc.l1$res3$name)
w <- which(!is.na(x))
uid$pt2.gpc.l1.r3 <- rep(NA,dim(uid)[1])
uid$pt2.gpc.l1.r3[w] <- pt2a.gpc.l1$res3$ratio[x[w]]

x <-  match(uid$name,pt2a.gpc.l2$res1$name)
w <- which(!is.na(x))
uid$pt2.gpc.l2.r1 <- rep(NA,dim(uid)[1])
uid$pt2.gpc.l2.r1[w] <- pt2a.gpc.l2$res1$ratio[x[w]]
x <-  match(uid$name,pt2a.gpc.l2$res2$name)
w <- which(!is.na(x))
uid$pt2.gpc.l2.r2 <- rep(NA,dim(uid)[1])
uid$pt2.gpc.l2.r2[w] <- pt2a.gpc.l2$res2$ratio[x[w]]
x <-  match(uid$name,pt2a.gpc.l2$res3$name)
w <- which(!is.na(x))
uid$pt2.gpc.l2.r3 <- rep(NA,dim(uid)[1])
uid$pt2.gpc.l2.r3[w] <- pt2a.gpc.l2$res3$ratio[x[w]]

x <-  match(uid$name,pt5c.npc.l1$res1$name)
w <- which(!is.na(x))
uid$pt5c.npc.l1.r1 <- rep(NA,dim(uid)[1])
uid$pt5c.npc.l1.r1[w] <- pt5c.npc.l1$res1$ratio[x[w]]
x <-  match(uid$name,pt5c.npc.l1$res2$name)
w <- which(!is.na(x))
uid$pt5c.npc.l1.r2 <- rep(NA,dim(uid)[1])
uid$pt5c.npc.l1.r2[w] <- pt5c.npc.l1$res2$ratio[x[w]]
x <-  match(uid$name,pt5c.npc.l1$res3$name)
w <- which(!is.na(x))
uid$pt5c.npc.l1.r3 <- rep(NA,dim(uid)[1])
uid$pt5c.npc.l1.r3[w] <- pt5c.npc.l1$res3$ratio[x[w]]

x <-  match(uid$name,hs1.gpc.l1$res1$name)
w <- which(!is.na(x))
uid$hs1.gpc.l1.r1 <- rep(NA,dim(uid)[1])
uid$hs1.gpc.l1.r1[w] <- hs1.gpc.l1$res1$ratio[x[w]]
x <-  match(uid$name,hs1.gpc.l1$res2$name)
w <- which(!is.na(x))
uid$hs1.gpc.l1.r2 <- rep(NA,dim(uid)[1])
uid$hs1.gpc.l1.r2[w] <- hs1.gpc.l1$res2$ratio[x[w]]
x <-  match(uid$name,hs1.gpc.l1$res3$name)
w <- which(!is.na(x))
uid$hs1.gpc.l1.r3 <- rep(NA,dim(uid)[1])
uid$hs1.gpc.l1.r3[w] <- hs1.gpc.l1$res3$ratio[x[w]]

x <-  match(uid$name,gc1.l1$res1$name)
w <- which(!is.na(x))
uid$gc1.r1 <- rep(NA,dim(uid)[1])
uid$gc1.r1[w] <- gc1.l1$res1$ratio[x[w]]
x <-  match(uid$name,gc1.l1$res2$name)
w <- which(!is.na(x))
uid$gc1.r2 <- rep(NA,dim(uid)[1])
uid$gc1.r2[w] <- gc1.l1$res2$ratio[x[w]]
x <-  match(uid$name,gc1.l1$res3$name)
w <- which(!is.na(x))
uid$gc1.r3 <- rep(NA,dim(uid)[1])
uid$gc1.r3[w] <- gc1.l1$res3$ratio[x[w]]

x <-  match(uid$name,hs1.npc.l1$res1$name)
w <- which(!is.na(x))
uid$hs1.npc.l1.r1 <- rep(NA,dim(uid)[1])
uid$hs1.npc.l1.r1[w] <- hs1.npc.l1$res1$ratio[x[w]]
x <-  match(uid$name,hs1.npc.l1$res2$name)
w <- which(!is.na(x))
uid$hs1.npc.l1.r2 <- rep(NA,dim(uid)[1])
uid$hs1.npc.l1.r2[w] <- hs1.npc.l1$res2$ratio[x[w]]
x <-  match(uid$name,hs1.npc.l1$res3$name)
w <- which(!is.na(x))
uid$hs1.npc.l1.r3 <- rep(NA,dim(uid)[1])
uid$hs1.npc.l1.r3[w] <- hs1.npc.l1$res3$ratio[x[w]]

x <-  match(uid$name,pt5c.gpc.l1$res1$name)
w <- which(!is.na(x))
uid$pt5c.gpc.l1.r1 <- rep(NA,dim(uid)[1])
uid$pt5c.gpc.l1.r1[w] <- pt5c.gpc.l1$res1$ratio[x[w]]
x <-  match(uid$name,pt5c.gpc.l1$res2$name)
w <- which(!is.na(x))
uid$pt5c.gpc.l1.r2 <- rep(NA,dim(uid)[1])
uid$pt5c.gpc.l1.r2[w] <- pt5c.gpc.l1$res2$ratio[x[w]]
x <-  match(uid$name,pt5c.gpc.l1$res3$name)
w <- which(!is.na(x))
uid$pt5c.gpc.l1.r3 <- rep(NA,dim(uid)[1])
uid$pt5c.gpc.l1.r3[w] <- pt5c.gpc.l1$res3$ratio[x[w]]
}

all <- uid
allm <- as.matrix(all[,4:ncol(all)])
allmi <- allm
allmi[which(!is.na(allmi))] <- 1
allmi[which(is.na(allmi))] <- 0
allmi.colors <- rep("white",dim(allmi)[1])
allmi.colors[which(info$perm==1)] <- "red"
allmi.colors[which(info$species!=0)] <- "blue"
heatmap(allmi,scale="none",col=c("white","midnightblue"),labRow=F,RowSideColors=allmi.colors)

aID <- rep("",nrow(all))
species <- rep(0,nrow(all))
cut <- rep(0,nrow(all))
pos <- rep(0,nrow(all))
neg <- rep(0,nrow(all))
perm <- rep(0,nrow(all))
info <- data.frame(aID,species,cut,pos,neg,perm)

info$perm[which(apply(all[,1:2],1,function(a) grep(as.character(a[1]),pattern="permutation"))==1)] <- 1
info$pos[which(apply(all[,1:2],1,function(a) grep(as.character(a[1]),pattern="POSITIVE_CONTROL"))==1)] <- 1
info$neg[which(apply(all[,1:2],1,function(a) grep(as.character(a[1]),pattern="NEGATIVE_CONTROL"))==1)] <- 1
info$species[which(apply(all[,1:2],1,function(a) grep(as.character(a[1]),pattern="species=human"))==1)] <- 1
info$species[which(apply(all[,1:2],1,function(a) grep(as.character(a[1]),pattern="species=chimp"))==1)] <- 2
info$cut[which(apply(all[,1:2],1,function(a) grep(as.character(a[1]),pattern="cut=1_"))==1)] <- 1
info$cut[which(apply(all[,1:2],1,function(a) grep(as.character(a[1]),pattern="cut=2_"))==1)] <- 2
info$cut[which(apply(all[,1:2],1,function(a) grep(as.character(a[1]),pattern="cut=3_"))==1)] <- 3
info$cut[which(apply(all[,1:2],1,function(a) grep(as.character(a[1]),pattern="cut=4_"))==1)] <- 4
info$cut[which(apply(all[,1:2],1,function(a) grep(as.character(a[1]),pattern="cut=5_"))==1)] <- 5
info$cut[which(apply(all[,1:2],1,function(a) grep(as.character(a[1]),pattern="cut=6_"))==1)] <- 6
info$aID <- apply(all,1,function(a) unlist(strsplit(as.character(a[1]),split="__"))[1] )

uaID <- unique(sort(info$aID))
hasSp <- 1:length(uaID)*0
for (i in 1:length(uaID)) {
	ws <- which(as.character(info$aID)==uaID[i])
	if (any(info$species[ws]==1) & any(info$species[ws]==2) ) {hasSp[i] <- 1}
}

uaID.f <- uaID[which(hasSp==1)]
sm <- matrix(0,length(uaID.f),length(c(paste("Hs_",colnames(allm),sep=""),paste("Pt_",colnames(allm),sep=""))))
colnames(sm) <- c(paste("Hs_",colnames(allm),sep=""),paste("Pt_",colnames(allm),sep=""))
rownames(sm) <- uaID.f
for (i in 1:length(uaID.f)) {
	whs <- which(as.character(info$aID)==uaID.f[i] & info$species==1 & info$perm==0)
	wpt <- which(as.character(info$aID)==uaID.f[i] & info$species==2 & info$perm==0)
	
	tv1 <- apply(rbind(allm[whs,],allm[whs,]),2,function(a) if (any(is.finite(a))) {max(a,na.rm=T)} else {NA})
	tv2 <- apply(rbind(allm[wpt,],allm[wpt,]),2,function(a) if (any(is.finite(a))) {max(a,na.rm=T)} else {NA})
	sm[i,] <- c(tv1,tv2)
}
write.table(sm,"2017jan15_allSppData.txt",sep="\t",quote=F)
sm0 <- sm
sm0[which(is.na(sm0))] <- 0
p <- prcomp(sm0)
pdf("2017jan15_pcaHeatmap.pdf")
#heatmap(p$rotation,col=cvBWO,scale="none",Rowv=NA,Colv=NA)
heatmap(p$rotation,scale="none",Rowv=NA,Colv=NA)
dev.off()
pdf("2017jan15_pcaBarplotSdev.pdf")
barplot(p$sdev,border=NA)
dev.off()

tm <- matrix(0,ncol(sm),ncol(sm))
for (i in 1:dim(sm)[2]) { for (j in 1:dim(sm)[2]) { tm[i,j] <- cor(sm[which(!is.na(sm[,i]) & !is.na(sm[,j])),i],sm[which(!is.na(sm[,i]) & !is.na(sm[,j])),j]) }}
#heatmap(cbind(-1,tm),scale="none",col=cvBWO,labRow=colnames(sm),labCol=c(-1,colnames(sm)))
heatmap(cbind(-1,tm),scale="none",labRow=colnames(sm),labCol=c(-1,colnames(sm)))
write.table(tm, "test.xls", sep="\t", col.names=T, row.names=F, quote=F)

###########################################################
# end of saved file
###########################################################

# examining site effects
numShuffles <- 100
perms <- sort(unique(info$aID[which(info$perm==1)]))
vcorperms <- matrix(0,numShuffles,length(perms))
colnames(vcorperms) <- perms
ctEffect <- vcorperms
ctIntEffects <- vcorperms
siteEffects <- vcorperms

testFrac <- 0.8

for (i in 1:length(perms)) {
	wp <- which(info$aID==perms[i] & info$perm==1)
	ntest <- round(testFrac*length(wp),0)
	tncol <- length(unlist(strsplit(as.character(uid$name[wp]),perl=T,split="_")))/length(wp)
	tn <- matrix(unlist(strsplit(as.character(uid$name[wp]),perl=T,split="_")),length(wp),tncol,byrow=T)[,c(1,4,6,8,17)]
	permid <- apply(tn,1,function(a) paste(a[1],a[2],collapse="_"))
	tv <- strsplit(tn[1,4],split=",")
	tn2 <- matrix(unlist(strsplit(tn[,4],split=",")),length(wp),length(tv[[1]]),byrow=T)
	tn2m <- apply(tn2,c(1,2),function(a) if (length(grep("h",a))>0) {1} else {0} )
	rownames(tn2m) <- permid
	colnames(tn2m) <- matrix(unlist(strsplit(tn2[1,],split=">")),length(tn2[1,]),2,byrow=T)[,1]
	nsites <- length(colnames(tn2m))
        sbp <- as.numeric(tn[1,5])
        spos <- as.numeric(matrix(unlist(strsplit(colnames(tn2m),split=":")),nsites,2,byrow=T)[,1])

        singles.p <- matrix(0,numShuffles,sbp)
        singles.f <- matrix(0,numShuffles,sbp)
        spp.p <- matrix(0,numShuffles,sbp)
        spp.f <- matrix(0,numShuffles,sbp)
        
  	for (j in 1:numShuffles) {	
		od <- c(1,sample(2:length(wp)))
		
		c1t <- cbind(allm[wp,10],0,tn2m)[od,][1:ntest,]
		c2t <- cbind(allm[wp,11],0,tn2m)[od,][1:ntest,]
		c3t <- cbind(allm[wp,12],0,tn2m)[od,][1:ntest,]
		x1t <- cbind(allm[wp,4],1,tn2m)[od,][1:ntest,]
		x2t <- cbind(allm[wp,5],1,tn2m)[od,][1:ntest,]
		x3t <- cbind(allm[wp,6],1,tn2m)[od,][1:ntest,]
	
		c1v <- cbind(allm[wp,10],0,tn2m)[od,][(ntest+1):length(wp),]
		c2v <- cbind(allm[wp,11],0,tn2m)[od,][(ntest+1):length(wp),]
		c3v <- cbind(allm[wp,12],0,tn2m)[od,][(ntest+1):length(wp),]
		x1v <- cbind(allm[wp,4],1,tn2m)[od,][(ntest+1):length(wp),]
		x2v <- cbind(allm[wp,5],1,tn2m)[od,][(ntest+1):length(wp),]
		x3v <- cbind(allm[wp,6],1,tn2m)[od,][(ntest+1):length(wp),]

        	#dft <- data.frame(rbind(x1t,x2t,x3t))
        	#dfv <- data.frame(rbind(x1v,x2v,x3v))
	      	dft <- data.frame(rbind(x1t,x2t,x3t))
          dfv <- data.frame(rbind(x1v,x2v,x3v))
        	dft <- data.frame(rbind(c1t,c2t,c3t,x1t,x2t,x3t))
        	dfv <- data.frame(rbind(c1v,c2v,c3v,x1v,x2v,x3v))
  
          #colnames(dft) <- c("enh",colnames(tn2m))
        	#colnames(dfv) <- c("enh",colnames(tn2m))      		
        	colnames(dft) <- c("enh",colnames(tn2m))
        	colnames(dfv) <- c("enh",colnames(tn2m)) 
        	colnames(dft) <- c("enh","spp",colnames(tn2m))
        	colnames(dfv) <- c("enh","spp",colnames(tn2m))
        	
        	gx <- glm(enh ~ .^2, data=dft)
        	csl <- length(coef(summary(gx))[,1])        	
        	vp <- predict.glm(gx,newdata=dfv)
        	
        	tvcorperm <- round(cor(vp,dfv$enh),2)
        	vcorperms[j,i] <- tvcorperm
        	
        	e10p.single <- -log10(coef(summary(gx))[,4])[2:csl][1:(nsites+1)]
        	effc.single <-        coef(summary(gx))[,1][2:csl][1:(nsites+1)]
        	e10p.spp    <- -log10(coef(summary(gx))[,4])[2:csl][(nsites+2):(2*nsites+1)]
        	effc.spp <- coef(summary(gx))[,1][2:csl][(nsites+2):(2*nsites+1)]
        	e10p.dbl <- -log10(coef(summary(gx))[,4])[2:csl][(2*nsites+2):(csl-1)]
        	effc.dbl <- coef(summary(gx))[,1][2:csl][(2*nsites+2):(csl-1)]
        	
        	singles.p[j,spos] <- e10p.single[2:length(e10p.single)]
        	singles.f[j,spos] <- effc.single[2:length(effc.single)]
        	f.max <- max(abs(singles.f))
        	f.norm <- round(100*(singles.f/f.max+1),0)
        	f.norm[which(f.norm<1)] <- 1
        	
        	spp.p[j,spos] <- e10p.spp[1:length(e10p.spp)]
        	spp.f[j,spos] <- effc.spp[1:length(effc.spp)]
        	f.max <- max(abs(spp.f))
        	f.norm <- round(100*(spp.f/f.max+1),0)
        	f.norm[which(f.norm<1)] <- 1
        	
        	ctEffect[j,i] <- abs(coef(summary(gx))[2,1])
        	ctIntEffects[j,i] <- sum(abs(effc.spp))
        	siteEffects[j,i] <- sum(abs(effc.single))
        	#doubles.p <- matrix(0,sbp,sbp)
        	#doubles.f <- matrix(0,sbp,sbp)
        	#dref <- matrix(as.numeric(unlist(strsplit(names(e10p.dbl),split="[`:]"))),length(e10p.dbl),length(unlist(strsplit(names(e10p.dbl),split="[`:]")))/length(e10p.dbl),byrow=T)[,c(2,6)]
        	#for (j in 1:length(e10p.dbl)) {
        	#	doubles.p[dref[j,1],dref[j,2]] <- e10p.dbl[j]
        	#	doubles.p[dref[j,2],dref[j,1]] <- e10p.dbl[j]
        	#	doubles.f[dref[j,1],dref[j,2]] <- effc.dbl[j]
        	#	doubles.f[dref[j,2],dref[j,1]] <- effc.dbl[j]
        	#}
        	#doubles.pn <- doubles.p
        	#doubles.pn[which(doubles.pn>5)] <- 5
        	#fd.max <- max(abs(doubles.f))
        	#pdf(paste("pairEffect_pValHeatmap_",perms[i],".pdf",sep=""),useDingbats=F)
        	#heatmap(doubles.pn[rev(spos),spos],Rowv=NA,Colv=NA,scale="none",col=cvWR,labRow=rev(spos),labCol=spos)
        	#dev.off()
        	#pdf(paste("pairEffect_effectHeatmap_",perms[i],".pdf",sep=""),useDingbats=F)
        	#heatmap(cbind(-fd.max,fd.max,doubles.f[rev(spos),spos]),Rowv=NA,Colv=NA,scale="none",col=cvBWO,labRow=rev(spos),labCol=c(-1,1,spos))
        	#dev.off()
        	
        	doubles.p <- matrix(0,sbp,sbp)
        	doubles.f <- matrix(0,sbp,sbp)
        	dref <- matrix(as.numeric(unlist(strsplit(names(e10p.dbl),split="[`:]"))),length(e10p.dbl),length(unlist(strsplit(names(e10p.dbl),split="[`:]")))/length(e10p.dbl),byrow=T)[,c(2,6)]
        	for (j in 1:length(e10p.dbl)) {
        	doubles.p[dref[j,1],dref[j,2]] <- e10p.dbl[j]
        	doubles.p[dref[j,2],dref[j,1]] <- e10p.dbl[j]
        	doubles.f[dref[j,1],dref[j,2]] <- effc.dbl[j]
        	doubles.f[dref[j,2],dref[j,1]] <- effc.dbl[j]
        	}
        	doubles.pn <- doubles.p
        	doubles.pn[which(doubles.pn>5)] <- 5
        	fd.max <- max(abs(doubles.f))
        	pdf(paste("pairEffect_pValHeatmap_",perms[i],".pdf",sep=""),useDingbats=F)
        	heatmap(doubles.pn[rev(spos),spos],Rowv=NA,Colv=NA,scale="none",col=cvWR,labRow=rev(spos),labCol=spos)
        	dev.off()
        	pdf(paste("pairEffect_effectHeatmap_",perms[i],".pdf",sep=""),useDingbats=F)
        	#heatmap(cbind(-fd.max,fd.max,doubles.f[rev(spos),spos]),Rowv=NA,Colv=NA,scale="none",col=cvBWO,labRow=rev(spos),labCol=c(-1,1,spos))
        	heatmap(cbind(-fd.max,fd.max,doubles.f[rev(spos),spos]),Rowv=NA,Colv=NA,scale="none", labRow=rev(spos),labCol=c(-1,1,spos))
        	dev.off()
 	}
 	
 	pdf(paste("effectHeatmap_",perms[i],".pdf",sep=""),useDingbats=F)
 	par(mfrow=c(2,2))
 	boxplot(singles.p,col="red",pch=19,axes=F,main="per site effect"); axis(1,at=spos,las=3);axis(2)
 	boxplot(spp.p,col="red",pch=19,axes=F,main="species effect"); axis(1,at=spos,las=3);axis(2) 	
 	bpcol = rep("light grey",sbp)
 	bpcol[which(apply(singles.f,2,mean)>0)] <- "orange"
 	bpcol[which(apply(singles.f,2,mean)<0)] <- "blue"
 	boxplot(singles.f,col=bpcol,pch=19,axes=F); axis(1,at=spos,las=3);axis(2)
 	bpcol = rep("light grey",sbp)
 	bpcol[which(apply(spp.f,2,mean)>0)] <- "orange"
 	bpcol[which(apply(spp.f,2,mean)<0)] <- "blue"
 	boxplot(spp.f,col=bpcol,pch=19,axes=F); axis(1,at=spos,las=3);axis(2)
 	dev.off()

}

pdf("2016nov25_allValidationCorrelations.pdf",useDingbats=F)
boxplot(vcorperms,col="grey",pch=19,las=3)
dev.off()
pdf("2016nov30_effectSizes.pdf",useDingbats=F,height=10.5,width=4)
par(mfrow=c(4,1))
boxplot(ctEffect,col="grey",pch=19,notch=T,las=3,ylim=c(0,0.3),outline=F,main="cell type effect")
boxplot(ctIntEffects,col="grey",pch=19,notch=T,las=3,ylim=c(0,0.3),outline=F,main="cell type : site interaction effects")
boxplot(siteEffects,col="grey",pch=19,notch=T,las=3,ylim=c(-0,0.3),outline=F,main="site effects")
boxplot(log2(siteEffects/ctIntEffects),col="grey",pch=19,notch=F,las=3,outline=F,main="log2(site/cellType) effects")
dev.off()

###### examine other celltype data
setwd("/Users/haneryu/Dropbox/Hane_HARMPRAdata_Sean/data")
load("/Users/haneryu/Dropbox/Hane_HARMPRAdata_Sean/dataworkspace.RData")

library(gcrma)
library(limma)

targets <- data.frame(as.factor(c(rep("HsSeq",27),rep("PtSeq",27))),as.factor(c( rep( c( rep("HsCell",6),rep("PtCell",9),rep("HsCell",3),rep("MmTestis",3),rep("HsCell",3),rep("PtCell",3)),2))),as.factor(rep(c(rep("GPC",12),rep("NPC",3),rep("GPC",3),rep("na",3),rep("NPC",3),rep("GPC",3)),2)))
colnames(targets) <- c("seqOrigin","cellType","time")
harGroup <- factor(paste(targets$cellType,targets$seqOrigin,sep="."))
harDesign <- model.matrix(~ 0 + harGroup)
colnames(harDesign) <- levels(harGroup)
fit  <- lmFit(sm,harDesign)

# need to revise below this line
cm <- makeContrasts(
	HsVsPtSeq.HsCell   = HsCell.HsSeq  - HsCell.PtSeq,
	HsVsPtSeq.PtCell   = PtCell.HsSeq  - PtCell.PtSeq,
	HsVsPtCell.HsSeq   = HsCell.HsSeq  - PtCell.PtSeq,
	HsVsPtCell.PtSeq   = HsCell.PtSeq  - PtCell.PtSeq,
	levels=harDesign)

fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)

tt.seq.HsCell <- topTable(fit2, coef="HsVsPtSeq.HsCell",number=nrow(smng))
tt.seq.PtCell <- topTable(fit2, coef="HsVsPtSeq.PtCell",number=nrow(smng))
tt.cell.HsSeq <- topTable(fit2, coef="HsVsPtCell.HsSeq",number=nrow(smng))
tt.cell.PtSeq <- topTable(fit2, coef="HsVsPtCell.PtSeq",number=nrow(smng))

write.table(tt.seq.HsCell,"2016dec23_seqDiffsHs.txt",sep="\t",quote=F)
write.table(tt.seq.PtCell,"2016dec23_seqDiffsPt.txt",sep="\t",quote=F)
write.table(tt.cell.HsSeq,"2016dec23_cellDiffsHs.txt",sep="\t",quote=F)
write.table(tt.cell.PtSeq,"2016dec23_cellDiffsPt.txt",sep="\t",quote=F)

allDEE <- matrix(0,nrow(sm),4)
rownames(allDEE) <- rownames(sm)
colnames(allDEE) <- c("seq.HsCell","seq.PtCell","cell.HsSeq","cell.PtSeq")
allDEE[,1] <- tt.seq.HsCell$adj.P.Val[match(rownames(allDEE),rownames(tt.seq.HsCell))]
allDEE[,2] <- tt.seq.PtCell$adj.P.Val[match(rownames(allDEE),rownames(tt.seq.PtCell))]
allDEE[,3] <- tt.cell.HsSeq$adj.P.Val[match(rownames(allDEE),rownames(tt.cell.HsSeq))]
allDEE[,4] <- tt.cell.PtSeq$adj.P.Val[match(rownames(allDEE),rownames(tt.cell.PtSeq))]
minP <- apply(allDEE,1,min)
wsig <- names(minP)[which(minP<0.01)]
sm.dee <- log2(sm[match(wsig,rownames(sm)),])
sm.dee <- sm.dee[which(apply(sm.dee,1,function(a) max(a)-min(a)) >= 0.4),]
sm.dee[which(sm.dee>0.5)] <- 0.5
sm.dee[which(sm.dee < -0.5)] <- -0.5

heatmap(sm.dee,col=cvBWO,scale="none")

