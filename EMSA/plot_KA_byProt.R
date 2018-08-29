## make box plot for summary KD values
library(ggplot2)
library(scales)

summary_Kd <- read.table("EMSA_summary_KD.txt", sep="\t", header=T)
head(summary_Kd)

KA = 1/summary_Kd$KD*10^9
KA.se = KA*summary_Kd$s.e./summary_Kd$KD
summary_Kd = cbind(summary_Kd, KA, KA.se)
summary_Kd

write.table(summary_Kd, "EMSA_summary_KA.txt", sep="\t")

#########################################
# generate multiple proteins on one plot
#########################################

# Kd_sum <- subset(summary_Kd, HD == "AncZenHD")
# Kd_sum <- subset(summary_Kd, HD== "AncBcdHD" | HD=="Q50K" 
#                          | HD=="K50Q" | HD== "AncZenHD")
# Kd_sum <- subset(summary_Kd, HD== "AncBcdHD" | HD== "AncZenHD")
# Kd_sum <- subset(summary_Kd, HD== "AncZenHD" | HD== "Q50K")
# Kd_sum <- subset(summary_Kd, HD=="I13L" | HD=="K50Q" | HD=="T42E" | HD=="A43R"
#                  | HD=="I58K" | HD=="R2K" | HD=="R55K" | HD=="T7A" | HD=="R54M"
#                  | HD=="A28R" | HD=="L31R" | HD=="AncBcdHD" | HD=="AncZenHD"
#                  | HD=="Q50K-M54R" | HD=="Q50K" | HD=="M54R")
# Kd_sum <- subset(summary_Kd, HD=="I13L" | HD=="K50Q" | HD=="T42E" | HD=="A43R"
#                  | HD=="I58K" | HD=="R2K" | HD=="R55K" | HD=="T7A" | HD=="R54M"
#                  | HD=="A28R" | HD=="L31R" | HD=="AncBcdHD")
Kd_sum <- subset(summary_Kd, HD=="AncZenHD" | HD=="Q50K-M54R" | HD=="Q50K" | HD=="M54R")

Kd_sum <- subset(Kd_sum, Motif == "BM" | Motif == "ZM")
Kd_sum

# orderlist <- c("AncZenHD", "Q50K", "K50Q", "AncBcdHD")
# orderlist <- c("AncZenHD", "AncBcdHD")
# orderlist <- c("AncZenHD", "Q50K")
# orderlist <- c("AncBcdHD", "R2K", "T7A", "I13L", "A28R", "L31R", "T42E", "A43R", "K50Q",
#                "R54M", "R55K", "I58K", "Q50K", "M54R", "Q50K-M54R", "AncZenHD")
# orderlist <- c("AncBcdHD", "R2K", "T7A", "I13L", "A28R", "L31R", "T42E", "A43R", 
#                "R54M", "R55K", "I58K", "Q50K", "Q50K-M54R", "K50Q", "M54R", "AncZenHD")
# orderlist <- c("AncBcdHD", "R2K", "T7A", "I13L", "A28R", "L31R", "T42E", "A43R", 
#                "K50Q", "R54M", "R55K", "I58K")
orderlist <- c("AncZenHD", "Q50K", "M54R", "Q50K-M54R")


Kd_sum <- transform(Kd_sum, HD = factor(HD, levels = orderlist))

# orderlist <- c("BCTL", "BM5", "BM2", "BM4", "BM", "BM6", "ZM")
# Kd_sum <- transform(Kd_sum, Motif = factor(Motif, levels = orderlist))

# g <- ggplot(data = Kd_sum, aes(x = Motif, y = KA, fill=HD)) +
#   geom_bar(width=0.5, position=position_dodge(width=.5), stat="identity") +
#   scale_fill_manual(values=c("#000000", "#999999")) +
#   ylab(bquote(K[A] ~ M^{-1})) + xlab("") +
#   scale_y_log10(breaks = 10^(6:9),
#                 labels = trans_format("log10", math_format(10^.x)), expand=c(0, 0)) +
#   coord_cartesian(ylim=c(1e6,1e9)) +
#   annotation_logticks(sides="l")
# g

g <- ggplot(data = Kd_sum, aes(x = factor(HD), y = KA, colour = Motif)) +
  geom_bar(aes(fill = Motif), width=0.5, position=position_dodge(width=.5), stat="identity") +
  ylab(bquote(K[A] ~ M^{-1})) + xlab("") +
  scale_y_log10(breaks = 10^(6:10),
                labels = trans_format("log10", math_format(10^.x)), expand=c(0, 0)) +
  coord_cartesian(ylim=c(1e6,1e10)) +
  annotation_logticks(sides="l")
g

gg_Kd <- g +
  theme_bw() +
  theme(axis.line=element_line(colour="black")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank())
gg_Kd

gg_Kd <- gg_Kd + theme(axis.text.x=element_text(size=12)) + 
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title.x=element_text(size=16)) +
  theme(axis.title.y=element_text(size=16)) +
  theme(legend.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg_se <- geom_errorbar(aes(ymin=KA-KA.se, ymax=KA+KA.se), width=.1,
                       position=position_dodge(width=.5))

gg_plot <- gg_Kd + gg_se
gg_plot

ggsave(gg_plot, file = "barPlot_KA_twoAnc.eps", width=4, height=4)
ggsave(gg_plot, file = "barPlot_KA_allBack.eps", width=8, height=4)
ggsave(gg_plot, file = "barPlot_KA_allForward.eps", width=4, height=4)
ggsave(gg_plot, file = "barPlot_KA_7motifs.eps", width=6, height=4)

#########################################
# generate multiple plots on one graph
#########################################
library(gridExtra)

summary_Kd <- subset(summary_Kd, Motif == "BM" | Motif == "ZM")
summary_Kd

HDs <- c("AncZenHD", "Q50K", "K50Q", "AncBcdHD")
HD.plots <- list()
for (prot in HDs) {
  Kd_sum <- subset(summary_Kd, HD==prot)
  
  gg_Kd <- ggplot(data = Kd_sum, aes(x = factor(HD), y = KA, colour = Motif)) +
    geom_bar(aes(fill = Motif), width=0.5, position=position_dodge(width=.5), stat="identity") +
    theme_bw() +   ylab(bquote(K[A] ~ M^{-1})) + xlab("") +
    scale_y_log10(breaks = 10^(6:9),
                  labels = trans_format("log10", math_format(10^.x)), expand=c(0, 0)) +
    coord_cartesian(ylim=c(1e6,3e9)) +
    annotation_logticks(sides="l") +
    theme(axis.line=element_line(colour="black")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank())
  
  gg_Kd <- gg_Kd + theme(axis.text.x=element_text(size=13)) + 
    theme(axis.text.y=element_text(size=13)) +
    theme(axis.title.x=element_text(size=16)) +
    theme(axis.title.y=element_text(size=16)) +
    theme(legend.text=element_text(size=14))
  gg_se <- geom_errorbar(aes(ymin=KA-s.e.KA, ymax=KA+s.e.KA, colour=Motif), width=.1,
                         position=position_dodge(width=.5))
  
  gg_plot <- gg_Kd + gg_se
  HD.plots[[prot]] = gg_plot
}

setEPS()
postscript(file="EMSA_barplot_4Anc_KA.eps", width=10, height=4, pointsize = 10)
grid.arrange(HD.plots[[1]], HD.plots[[2]], HD.plots[[3]], HD.plots[[4]], ncol=4)
dev.off()









