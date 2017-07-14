library("RColorBrewer")
library("ggplot2")
library("cowplot")

####Figure 1 A-D####
####Transcript length####

t<-read.csv("transcript.length",header=F,sep="\t")
row.names(t)<-t$V3
t$Length <-t$V2-t$V1+1
t.c<-t[grep("^NC",row.names(t), invert=TRUE),]
t.nc<-t[grep("^NC",row.names(t)),]
codes<-c(rep("Coding", dim( t.c)[1]), rep("Non-coding", dim( t.nc)[1]))
df<-data.frame("Length"=c(t.c$Length,t.nc$Length), "Gene_type"=codes)
g1<-ggplot(df,aes(x=Length, fill=Gene_type))+geom_density(alpha=0.75)+scale_fill_brewer(palette="Accent")+ylab("Density")+xlab("Transcript length [bp]")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+coord_cartesian(xlim=c(0, 10000))+theme_classic()+theme(legend.title=element_blank())+theme(text = element_text(size=14))

ggsave("1A.pdf", width = 5, height = 4)

####Trasncript number####
require(plyr)
require(dplyr)
require(ggplot2)
t<-read.csv("transcript.number",header=F,sep="\t")
row.names(t)<-t$V2
t.c<-t[grep("^NC",row.names(t), invert=TRUE),]
t.nc<-t[grep("^NC",row.names(t)),]
codes<-c(rep("Coding", dim( t.c)[1]), rep("Non-coding", dim( t.nc)[1]))
df<-data.frame("Transcripts"=c(t.c$V1,t.nc$V1), "Gene_type"=codes)
df$Transcripts[df$Transcripts > 10] <- 10
df2<-rbind(data.frame(count(df,Gene_type,Transcripts)), c("Non-coding", 5, NA),  c("Non-coding", 6, NA), c("Non-coding", 7, NA), c("Non-coding", 8, NA), c("Non-coding", 9, NA), c("Non-coding", 10, NA))
df2$n2<-as.numeric(df2$n)
positions <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
g2<-ggplot(df2,aes(x=Transcripts,y=n2,  fill=Gene_type))+geom_bar(position=position_dodge(), width=0.7, colour="black", stat="identity")+scale_fill_brewer(palette="Accent")+ylab("Count")+xlab("Transcript number")+theme_classic()+scale_x_discrete(limits = positions)+theme(legend.title=element_blank())+theme(text = element_text(size=14))

ggsave("1B.pdf", width = 5, height = 4)

####Number of exons####

t<-read.csv("exon.number",header=F,sep="\t")
row.names(t)<-t$V2
t.c<-t[grep("^NC",row.names(t), invert=TRUE),]
t.nc<-t[grep("^NC",row.names(t)),]
codes<-c(rep("Coding", dim( t.c)[1]), rep("Non-coding", dim( t.nc)[1]))
df<-data.frame("Exons"=c(t.c$V1,t.nc$V1), "Gene_type"=codes)
df$Exons[df$Exons > 10] <- 10
df2<-rbind(data.frame(count(df,Gene_type,Exons)), c("Non-coding", 5, NA),  c("Non-coding", 6, NA), c("Non-coding", 7, NA), c("Non-coding", 8, NA), c("Non-coding", 9, NA), c("Non-coding", 10, NA))
df2$n2<-as.numeric(df2$n)
positions <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
g3<-ggplot(df2,aes(x=Exons,y=n2,  fill=Gene_type))+geom_bar(position=position_dodge(), width=0.7, colour="black", stat="identity")+scale_fill_brewer(palette="Accent")+ylab("Count")+xlab("Exon number")+theme_classic()+scale_x_discrete(limits = positions)+theme(legend.title=element_blank())+theme(text = element_text(size=14))

ggsave("1C.pdf", width = 5, height = 4)

#####FPKM density plot#######################

t<-read.csv("all.genes.fpkm",header=T, row.names=1,sep="\t")
c<-as.vector(as.matrix(t[grep("^NC",row.names(t), invert=TRUE),]))
nc<-as.vector(as.matrix(t[grep("^NC",row.names(t)),]))
codes<-c(rep("Coding", length(c)), rep("Non-coding", length(nc)))
df<-data.frame("FPKM"=c(log2(c),log2(nc)), "Gene_type"=codes)
g4<-ggplot(df,aes(x=FPKM, fill=Gene_type)) + geom_density(alpha=0.75)+scale_fill_brewer(palette="Accent")+ylab("Density")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme_classic()+theme(legend.title=element_blank())+theme(text = element_text(size=14))

ggsave("1D.pdf", width = 5, height = 4)

####Figure 1E####
####Circos#######

library("OmicCircos")
add.alpha <- function(col, alpha=1){
   if(missing(col))
     stop("Please provide a vector of colours.")
   apply(sapply(col, col2rgb)/255, 2,
                      function(x)
                        rgb(x[1], x[2], x[3], alpha=alpha))
}
options(stringsAsFactors=FALSE)
seg.frame<-read.csv("chr.length",header=F,sep="\t")
names(seg.frame)<-c("seg.name","seg.Start","seg.End")
seg.frame$the.v<-rep(NA,20)
seg.frame$NO<-rep(NA,20)
seg.name<-seg.frame$seg.name
colors<-rainbow(20,alpha=0.5)
db<-segAnglePo(seg.frame, seg=seg.name)
t<-read.csv("circos.table",sep="\t",header=F)
names(t)<-c("CHR","START","END","BIN","CODING","TOTAL","nonTE","TE")
coding<-data.frame(seg.name=t$CHR, seg.po=t$START+250000, name1=t$CODING)
all<-data.frame(seg.name=t$CHR, seg.po=t$START+250000, name1=t$TOTAL)
nte<-data.frame(seg.name=t$CHR, seg.po=t$START+250000, name1=t$nonTE)
te<-data.frame(seg.name=t$CHR, seg.po=t$START+250000, name1=t$TE)
t2<-read.csv("centromere.locations",sep="\t",header=F)
cent<-data.frame(seg.name=t2$V1, start=t2$V2, end=t2$V3)
pdf("1E.pdf")
par(mar=c(2,2,2,2))
plot(c(-200,1000),c(-200,1000),type="n",axes=FALSE,xlab="",ylab="",main="")
circos(R=400,cir=db,type="chr",col="grey",print.chr.lab=TRUE,W=4)
circos(R=405, W=2, cir=db, mapping=cent, type="arc2", B=FALSE, col="black", lwd=2, cutoff=0)
c=brewer.pal(12, 'Paired')
circos(R=360, W=40, cir=db, mapping=coding, col.v=3, type="l", B=FALSE, col=c[1], lwd=0.5)
circos(R=320, W=40, cir=db, mapping=all, col.v=3, type="l", B=FALSE, col=c[2], lwd=0.5)
circos(R=280, W=40, cir=db, mapping=nte, col.v=3, type="l", B=FALSE, col=c[[3]], lwd=0.5)
circos(R=240, W=40, cir=db, mapping=te, col.v=3, type="l", B=FALSE, col=c[[4]], lwd=0.5)
dev.off()

####Figure 1F####
####TE fraction####

t<-read.csv("TE.fraction",sep="\t",header=T)
ggplot(t,aes(x=TE, y=Fraction, fill=Location))+geom_bar(position=position_dodge(), width=0.7, colour="black", stat="identity")+scale_fill_brewer(palette="Accent")+ylab("Fraction")+xlab("Transposable Element")+theme_classic()+theme(legend.title=element_blank())+
coord_flip()+theme(text = element_text(size=14))

ggsave("1F.pdf", width = 5, height = 4)

####Fig 2B####
####simulations####
t<-read.csv("simulations.synteny", sep="\t",header=F)
names(t)<-c("Type","Count","Significance")
ggplot(t,aes(x=factor(Type), y=Count, color=factor(Type), shape=Significance)) + geom_point(alpha=0.75, size=4)+scale_color_brewer(palette="Accent")+ylab("Loci count")+xlab("Control dataset")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme_classic()+theme(legend.title=element_blank())+ expand_limits(y = 1000)+geom_hline(yintercept = 1434)+annotate("text", 3.0, 1416, label = "Biological dataset = 1434")+scale_shape_manual(values=c(19, 17))+theme(text = element_text(size=14))

ggsave("2B.pdf", width = 5, height = 4)

####Fig 2C####
####simulations####
t<-read.csv("simulations.homeology",sep="\t",header=F)
names(t)<-c("Type","Count","Significance")
ggplot(t,aes(x=factor(Type), y=Count, color=factor(Type), shape=Significance)) + geom_point(alpha=0.75, size=4)+scale_color_brewer(palette="Accent")+ylab("Loci count")+xlab("Control dataset")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme_classic()+theme(legend.title=element_blank())+ expand_limits(y = 1000)+geom_hline(yintercept = 1826)+annotate("text", 3.0, 1800, label = "Biological dataset = 1826")+scale_shape_manual(values=c(17,19))+theme(text = element_text(size=14))

ggsave("2C.pdf", width = 5, height = 4)

####Fig 2D####
####Ks density plot####
t<-read.csv("ks.conserved.random",sep="\t",header=F)
names(t)<-c("R","Type")
ggplot(t,aes(x=R, fill=Type)) + geom_density(alpha=0.75)+scale_fill_brewer(palette="Accent")+ylab("Density")+xlab("Ks value")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme_classic()+theme(legend.title=element_blank())+coord_cartesian(xlim=c(0, 2))+theme(text = element_text(size=14))

ggsave("2D.pdf", width = 5, height = 4)

####Fig 2E####
####Violin plot####
t<-read.csv("homeologous.ori.r",sep="\t",header=F)
names(t)<-c("R","Type")
ggplot(t,aes(x=Type,y=R))+geom_violin(aes(fill=Type),alpha=0.75)+geom_boxplot(aes(fill=Type),alpha=0.75, width=0.1)+scale_fill_brewer(palette="Accent")+ylab("Correlation coefficient")+xlab("")+theme_classic()+theme(legend.title=element_blank(), axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(text = element_text(size=14))

ggsave("2E.pdf", width = 5, height = 4)

####Fig3A####
####Samples with expression#########
t<-read.csv("all.genes.fpkm",header=T, row.names=1,sep="\t")
nonZero<-apply(t, 1, function(c)sum(c>0.1))
t$Expression<-nonZero
t.c<-t[grep("^NC",row.names(t), invert=TRUE),]
t.nc<-t[grep("^NC",row.names(t)),]
codes<-c(rep("Coding", dim(t.c)[1]), rep("Non-coding", dim(t.nc)[1]))
df<-data.frame("Samples_with_expression"=c(t.c$Expression,t.nc$Expression), "Gene_type"=codes)
ggplot(df,aes(x=Samples_with_expression, fill=Gene_type)) + geom_density(alpha=0.75)+scale_fill_brewer(palette="Accent")+ylab("Density")+xlab("Samples with expression")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme_classic()+theme(legend.title=element_blank())+theme(text = element_text(size=14))

ggsave("3A.pdf", width = 5, height = 4)

####Fig3B####
####Counts plot#####
t<-read.csv("sample.expression",header=F,sep="\t")
names(t)<-c("Sample", "Count","Type")
ggplot(t,aes(x=Sample,y=Count,  fill=Type))+geom_bar(position=position_dodge(), width=0.7, colour="black", stat="identity")+scale_fill_brewer(palette="Accent")+ylab("Count")+xlab("Sample")+theme_classic()+theme(legend.title=element_blank())+coord_flip()+theme(text = element_text(size=14))

ggsave("3B.pdf", width = 5, height = 5)

####Fig3C####
####Counts plot2####
t<-read.csv("sample.merged.expression",header=F,sep="\t")
names(t)<-c("Sample", "Count","Type")
ggplot(t,aes(x=Sample,y=Count,  fill=Type))+geom_bar(position=position_dodge(), width=0.7, colour="black", stat="identity")+scale_fill_brewer(palette="Accent")+ylab("Count")+xlab("Sample")+theme_classic()+theme(legend.title=element_blank())+coord_flip()+theme(text = element_text(size=14))
ggsave("3C.pdf", width = 5, height = 4)

####Figure 3D####
####Heatmep####
#https://support.bioconductor.org/p/76250/
require(gtools)
require(RColorBrewer)
require(gplots)
library(pheatmap)
t<-read.csv("all.genes.fpkm",header=T, row.names=1,sep="\t")
difE<-apply(t, 1, function(c)sum(c>0.1))
t$DiffE<-difE
t.s<-subset(t, DiffE>0)
t.s$DiffE<-NULL
tl<-log1p(as.matrix(t.s))
tl.nc<-tl[grep("^NC",row.names(tl)),]
tls.nc<-1-cor(tl.nc,method="spearman")
mat<-tls.nc
o<-rownames(mat)
hc = hclust(as.dist(mat), method="complete")
mat = mat[hc$order, hc$order]
mat[lower.tri(mat)] = NA
diag(mat) <- NA
mat = mat[o, o]
pdf("3D.pdf", height=6, width=8)
pheatmap(mat, show_colnames = T, cluster_col = hc, cluster_row = hc, treeheight_col=0)
dev.off()


####Fig 3E####
####Violin plot####
t<-read.csv("specificity.tau", sep="\t",header=F)
names(t)<-c("R","Type","Place")
ggplot(t,aes(x=Place,y=R))+geom_violin(aes(fill=Type),alpha=0.75)+geom_boxplot(aes(fill=Type),alpha=0.75, width=0.1)+scale_fill_brewer(palette="Accent")+ylab("Correlation coefficient")+xlab("")+theme_classic()+theme(legend.title=element_blank(), axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(text = element_text(size=14))

ggsave("3E.pdf", width = 8, height = 5)

####Figure 4B and 4C####
####EXPRESSION####

t<-read.csv("GWAS.linc.ncgenes.fixed.fpkm",sep="\t",header=T)
t.s<-subset(t,GENE=="NC_GMAXST00018683")
ch2<-data.frame(SAMPLE=t.s$SAMPLE, FPKM=t.s$FPKM)
cp<-c("SEEDLING", "YOUNG_LEAF", "SAM","FLOWER_BUD1","FLOWER_BUD2","FLOWER_BUD3","FLOWER_BUD4","FLOWER1","FLOWER2","FLOWER3","FLOWER4")
cp2<-c(0,0,0.632526, 2.49547, 0.713063, 0, 1.97213, 3.82082, 1.41218, 0.352595, 0.912167)
ch<-data.frame(SAMPLE=cp,FPKM=cp2)
g1<-ggplot(data=ch, aes(x=SAMPLE, y=FPKM, group=1)) + geom_area(size=1.0, color="plum4", fill="plum4", alpha=0.5) + theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank())+ ggtitle("Chickpea lincRNA - Ca_linc_1449")+scale_x_discrete(expand=c(0.05,0.05))
g2<-ggplot(data=ch2, aes(x=SAMPLE, y=FPKM, group=1)) + geom_area(size=1.0, color="steelblue", fill="steelblue", alpha=0.5) + theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank())+ ggtitle("Soybean lincRNA - NC_GMAXST00018683")+scale_x_discrete(expand=c(0.05,0.05))
dev.off()
pdf("4B.pdf", height=6, width=8)
plot_grid(g2, g1, nrow=2, ncol=1)
dev.off()

p1<-c(0.0,0.0,0.359037492756,0.723781185486,0.136196655642,0.0349411497815,0.384902953211,0.0,0.0,0.0,0.0,29.1768408019,5.99656328175,0.392811269718,2.88122427581,0.1862844619,1.50786654119,0.0245239082223,0.0209491769573,0.193062458728,0.0,0.0197321674107,0.0,0.0,0.0,1.00701166474,10.0982992345,0.0,0.0,0.0,0.0,0.0,0.576311021003,0.494857718187,0.530754192615,0.983681162337,1.26807824449)
p2<-c("ROOT","LEAF3","SEED4","SEED3","POD_SEED3","POD_SEED2","POD_SEED1","STEM1","COTYLEDON1","LEAFBUD1","STEM2","FLOWER2","FLOWER3","FLOWER5","FLOWER4","POD1","SEED1","POD2","POD3","SEED2","LEAF1","LEAFBUD2","COTYLEDON2","LEAF2","LEAFBUD3","SHOOT_MERISTEM","FLOWER1","SEED5","LEAF_SD0","LEAF_SD1","LEAF_SD2","LEAF_SD3","SAM_SD0","SAM_SD1","SAM_SD2","SAM_SD3","SAM_SD4")
pf<-data.frame(SAMPLE=p2, FPKM=p1)
g3<-ggplot(data=pf, aes(x=SAMPLE, y=FPKM, group=1)) + geom_area(size=1.0, color="darkgreen", fill="darkgreen", alpha=0.5) + theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank())+ ggtitle("Soybean gene downstream of lincRNA - AGL6, rho=0.48")+scale_x_discrete(expand=c(0.05,0.05))

p1<-c(5.16267663573,2.66399499629,2.71984271395,2.36117122087,4.49864441777,4.63622400833,4.96950877402,5.6084154633,2.37606321208,8.94356841468,6.77180303802,8.39368087879,6.00860942944,11.6316891071,8.04675014791,3.67176204564,4.32686466366,3.96018133685,3.87971805957,4.4450919087,1.7403630026,7.39779390664,3.38775265246,4.84710767657,11.3235315503,10.9845393577,11.9847304047,1.04286209623,1.50018612503,0.717131888751,0.846707475277,1.16394792358,6.3452814983,8.89450579817,6.408468214,7.66002988792,5.91267054429)
p2<-c("ROOT","LEAF3","SEED4","SEED3","POD_SEED3","POD_SEED2","POD_SEED1","STEM1","COTYLEDON1","LEAFBUD1","STEM2","FLOWER2","FLOWER3","FLOWER5","FLOWER4","POD1","SEED1","POD2","POD3","SEED2","LEAF1","LEAFBUD2","COTYLEDON2","LEAF2","LEAFBUD3","SHOOT_MERISTEM","FLOWER1","SEED5","LEAF_SD0","LEAF_SD1","LEAF_SD2","LEAF_SD3","SAM_SD0","SAM_SD1","SAM_SD2","SAM_SD3","SAM_SD4")
pf<-data.frame(SAMPLE=p2, FPKM=p1)
g4<-ggplot(data=pf, aes(x=SAMPLE, y=FPKM, group=1)) + geom_area(size=1.0, color="darkgreen", fill="darkgreen", alpha=0.5) + theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank())+ ggtitle("Soybean gene upstream from lincRNA - THO complex subunit 7A, rho=0.5")+scale_x_discrete(expand=c(0.05,0.05))
pdf("4C.pdf", height=6, width=8)
plot_grid(g3, g4, nrow=2, ncol=1)
dev.off()

####Supplemetary####
####Fig S2####
t<-read.csv("centromere.fpkm",sep="\t",header=T)
g1<-ggplot(data=t, aes(x=SAMPLE, y=log1p(FPKM), group=GENE)) + geom_line(size=0.75, color="steelblue",alpha=0.2) + theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_x_discrete(expand=c(0,0))+theme(legend.title=element_blank())

t<-read.csv("centromere.control.fpkm",sep="\t",header=T)
g2<-ggplot(data=t, aes(x=SAMPLE, y=log1p(FPKM), group=GENE)) + geom_line(size=0.75, color="steelblue",alpha=0.2) + theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_x_discrete(expand=c(0,0))+theme(legend.title=element_blank())
pdf("S2.pdf", height=6, width=8)
plot_grid(g1, g2, labels=c("A", "B"), ncol = 1, nrow = 2)
dev.off()

####Fig S3A####
####Correlation coefficients####

t<-read.csv("correlations.0.5percent",sep="\t",header=T)
t2<-read.csv("correlations",sep="\t",header=T)
codes<-c(rep("All", dim(t)[1]), rep("Significant", dim(t2)[1]))
df<-data.frame("R"=c(t$R,t2$R), "Gene_type"=codes)
g1<-ggplot(df,aes(x=R, fill=Gene_type)) + geom_density(alpha=0.75)+scale_fill_brewer(palette="Accent")+ylab("Density")+xlab("Correlation coefficient")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme_classic()+theme(legend.title=element_blank())

####Fig S3B####
#####GO plot####
t<-read.csv("go.slim",header=F,sep="\t")
names(t)<-c("Desc", "GOtype","Gene_type", "Count")
t.s<-subset(t,Gene_type=="Non-Coding")
positions <-t.s$Desc
g2<-ggplot(t.s,aes(x=Desc,y=Count))+geom_bar(position=position_dodge(), width=0.7, colour="black", stat="identity", fill="plum4", alpha=0.75)+scale_fill_brewer(palette="Accent")+ylab("Count")+xlab("GO term")+theme_classic()+scale_x_discrete(limits = positions)+theme(legend.title=element_blank())+coord_flip()
pdf("S3.pdf", height=6, width=14)
plot_grid(g1, g2, labels=c("A", "B"), ncol = 2, nrow = 1)
dev.off()

####Fig S5####

t<-read.csv("GWAS.linc.ncgenes.fixed.fpkm",sep="\t",header=T)
sub<-"Days-to-flowering"
s<-subset(t, TYPE==sub)
g1<-ggplot(data=s, aes(x=SAMPLE, y=FPKM, group=GENE, colour=GENE)) + geom_line(size=1.0) + theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ggtitle(sub)+scale_colour_brewer(palette="Paired")

sub<-"Flowering-to-maturity"
s<-subset(t, TYPE==sub)
g2<-ggplot(data=s, aes(x=SAMPLE, y=FPKM, group=GENE, colour=GENE)) + geom_line(size=1.0) + theme_bw()+scale_colour_brewer(palette="Paired")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ggtitle(sub)

sub<-"Days-to-maturity"
s<-subset(t, TYPE==sub)
g3<-ggplot(data=s, aes(x=SAMPLE, y=FPKM, group=GENE, colour=GENE)) + geom_line(size=1.0) + theme_bw()+scale_colour_brewer(palette="Paired")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ggtitle(sub)

sub<-"Seed-weight"
s<-subset(t, TYPE==sub)
g4<-ggplot(data=s, aes(x=SAMPLE, y=FPKM, group=GENE, colour=GENE)) + geom_line(size=1.0) + theme_bw()+scale_colour_brewer(palette="Paired")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ggtitle(sub)

sub<-"Plant-height"
s<-subset(t, TYPE==sub)
g5<-ggplot(data=s, aes(x=SAMPLE, y=FPKM, group=GENE, colour=GENE)) + geom_line(size=1.0) + theme_bw()+scale_colour_brewer(palette="Paired")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ggtitle(sub)

sub<-"Seed-size"
s<-subset(t, TYPE==sub)
g6<-ggplot(data=s, aes(x=SAMPLE, y=FPKM, group=GENE, colour=GENE)) + geom_line(size=1.0) + theme_bw()+scale_colour_brewer(palette="Paired")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ggtitle(sub)

pdf("S5.pdf", height=12, width=20)
plot_grid(g1, g2, g3, g4, g5, g6, ncol = 2, nrow = 3)
dev.off()

