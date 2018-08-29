# Correlation of 6mer E-scores: AncBcd-HD VS Dro_Bcd-HD
library(dplyr)
library(ggplot2)
##################################################
# import 8mer E-score table
working_dir = "."
AB1 <- read.table(file.path(working_dir,"AB_a1_block3_8mers_11111111.txt"), header=T, sep="\t")
head(AB1)
AB2 <- read.table(file.path(working_dir,"AB_a2_block4_8mers_11111111.txt"), header=T, sep="\t")
head(AB2)

dro_8mer <- read.table(file.path(working_dir,"kf9_8mer_escores_forR_onlygoodNEW.txt"), header=T, sep="\t", skip=1)
head(dro_8mer)

# merge AncBcd and Bcd E-score data
t8mer = merge(dro_8mer, AB1 %>% select(-Median,-Z.score), by.x=c("X8mer","X8mer.rc"),by.y=c("X8.mer","X8.mer.1"), all=T)
t8mer = merge(t8mer, AB2 %>% select(-Median,-Z.score), by.x=c("X8mer","X8mer.rc"),by.y=c("X8.mer","X8.mer.1"), all=T)

colnames(t8mer)[7:8] = c("AncBcdHD_rep1", "AncBcdHD_rep2")
str(t8mer)

# quick plot on raw E-score
with(t8mer, {
  plot(New_GST.BcdHD_v14_181_R1, AncBcdHD_rep1, pch=16, cex=0.2)
  abline(a=0,b=1,col="red")
  
  plot(density(New_GST.BcdHD_v14_181_R1),pch=16)
  lines(density(AncBcdHD_rep1),col="red",pch=16)
  legend("topright",c("BcdHD","AncBcdHD"), col=c("black","red"),pch=16)
  
  plot(AncBcdHD_rep1, AncBcdHD_rep2, pch=16, cex=0.2)
  abline(a=0,b=1,col="red")
})

# normalize E-score
t8mer = t8mer %>% mutate(BcdHD_rep1.nn = (New_GST.BcdHD_v14_181_R1 - mean(New_GST.BcdHD_v14_181_R1))/sd(New_GST.BcdHD_v14_181_R1),
                         BcdHD_rep2.nn = (New_GST.BcdHD_v14_181_R2 - mean(New_GST.BcdHD_v14_181_R2))/sd(New_GST.BcdHD_v14_181_R2),
                         AncBcdHD_rep1.nn = (AncBcdHD_rep1 - mean(AncBcdHD_rep1))/sd(AncBcdHD_rep1),
                         AncBcdHD_rep2.nn = (AncBcdHD_rep2 - mean(AncBcdHD_rep2))/sd(AncBcdHD_rep2))
head(t8mer)

with(t8mer, {
  plot(BcdHD_rep2.nn, AncBcdHD_rep1.nn, pch=16, cex=0.2)
  abline(a=0,b=1,col="red")
  
  plot(density(BcdHD_rep2.nn),pch=16)
  lines(density(AncBcdHD_rep1.nn),col="red",pch=16)
  legend("topright",c("BcdHD","AncBcdHD"), col=c("black","red"),pch=16)
})

########################
load("kmer.6.RData")
head(kmer.6)
kmer.6[1,]
dim(kmer.6)

# extract 6mer E-scores by averaging E-scores of 8mers containing a 6mer
t6mer.mean <- apply(kmer.6, 1, function(motif) {
  kmer <- motif[1]
  idx <- unique(c(grep(kmer, t8mer$X8mer), grep(kmer, t8mer$X8mer.rc)))
  AB.6mer = mean(c(t8mer$AncBcdHD_rep1[idx], t8mer$AncBcdHD_rep2[idx]),na.rm=T)
  AB.nn.6mer = mean(c(t8mer$AncBcdHD_rep1.nn[idx], t8mer$AncBcdHD_rep2.nn[idx]),na.rm=T)
  Bcd.6mer = mean(c(t8mer$New_GST.BcdHD_v14_181_R1[idx],t8mer$New_GST.BcdHD_v14_181_R2[idx]),na.rm=T)
  Bcd.nn.6mer = mean(c(t8mer$BcdHD_rep1.nn[idx],t8mer$BcdHD_rep2.nn[idx]),na.rm=T)
  c(AB.6mer, AB.nn.6mer, Bcd.6mer, Bcd.nn.6mer)
})
dim(t6mer.mean)

t6mer = cbind(kmer.6, t(t6mer.mean))
str(t6mer)
colnames(t6mer) = c("X6mer", "X6mer.rc","AncBcdHD", "AncBcdHD.nn", "BcdHD", "BcdHD.nn")

with(t6mer, {
  print(cor(AncBcdHD, BcdHD))
  print(cor(AncBcdHD.nn, BcdHD.nn))
})

# quick plots
with(t6mer, {
  plot(BcdHD, AncBcdHD, pch=16, cex=0.2)
  abline(a=0,b=1,col="red")
  
  plot(density(BcdHD))
  lines(density(AncBcdHD),col="red")
  legend("topright",c("BcdHD","AncBcdHD"), col=c("black","red"),pch=16)
  
  plot(BcdHD.nn, AncBcdHD.nn, pch=16, cex=0.2)
  abline(a=0,b=1,col="red")
  
  plot(density(BcdHD.nn))
  lines(density(AncBcdHD.nn),col="red")
  legend("topright",c("BcdHD.norm","AncBcdHD.norm"), col=c("black","red"),pch=16)
})

AB_Bcd_plot = ggplot(t6mer, aes(x=AncBcdHD,y=BcdHD)) + geom_point(size=1.5) +
  scale_x_continuous(breaks=seq(-.5, .5, by=0.25), limits=c(-.35, .6), expand=c(0,0)) +
  scale_y_continuous(breaks=seq(-.5, .5, by=0.25), limits=c(-.35, .6), expand=c(0,0)) +
  stat_smooth(method="lm", se=FALSE, fill="black", size=1.25) +
  theme_bw() +
  theme(axis.line=element_line(colour="black")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.background = element_blank()) +
  theme(axis.title.y=element_text(size=12)) +
  ylab("E-score by Dro_BcdHD") +
  xlab("E-score by AncBcdHD")
ggsave(file.path(working_dir,"BcdHD_AncBcdHD_6mer_Escore_scatterplot.eps"),width=5.5,height=5)


##
