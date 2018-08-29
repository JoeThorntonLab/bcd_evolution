## kon-koff plot
library(ggplot2)
library(MASS)
library(scales)

data <- read.table("SPR_fourAnc_summary.txt", sep="\t", header=T)
data
summary_motifs <- data
summary_motifs <- data[1:8, 1:8]
summary_motifs

# summary_motifs <- summary_motifs[summary_motifs$HD == "AncZenHD" | summary_motifs$HD == "AncBcdHD",]
# summary_motifs <- summary_motifs[summary_motifs$HD != "AncBcdHD-K50q",]
# summary_motifs <- summary_motifs[summary_motifs$HD != "AncZenHD-q50K",]
# summary_motifs

# summary_motifs <- summary_motifs[summary_motifs$Motif == "ZM",]
summary_motifs <- summary_motifs[summary_motifs$Motif == "BM",]
summary_motifs

group1 <- c(1, 1, 2, 3)
summary_motifs <- cbind(summary_motifs, group1)
summary_motifs

# gg <- ggplot(summary_motifs, aes(x=koff, y=kon, colour=factor(group1), shape=HD))
gg <- ggplot(summary_motifs, aes(x=koff, y=kon, colour=Motif, shape=HD))
B_gg <- gg + geom_point(size = 6) +
#  scale_shape_manual(values=c(15, 17)) +
  scale_shape_manual(values=c(15, 10, 17, 9)) +
  scale_y_log10(limits=c(1e4, 1e7), breaks = 10^(4:7),
                labels = trans_format("log10", math_format(10^.x)), expand=c(0, 0)) +
  scale_x_log10(limits=c(.0001, 0.1), breaks = 10^(-4:-1),
                labels = trans_format("log10", math_format(10^.x)), expand=c(0, 0)) +
  ylab("kon (1/Ms)") + xlab("koff (1/s)") +
  theme_bw() +
  theme(axis.line=element_line(colour="black")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black"),
        plot.background = element_blank()) +
  annotation_logticks(sides="trbl")

B_gg

B_gg <- B_gg + theme(axis.text.x=element_text(size=12)) + 
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title.x=element_text(size=18)) +
  theme(axis.title.y=element_text(size=18)) +
  theme(title=element_text(size=16)) +
  theme(legend.text=element_text(size=12))
B_gg

# create a layer to plot the standard error of the means
B_gg_se_y <- geom_errorbar(aes(ymin=kon-s.e..kon, ymax=kon+s.e..kon),
                         data=summary_motifs, width=.2)
B_gg_se_x <- geom_errorbarh(aes(xmin=koff-s.e..koff, xmax=koff+s.e..koff),
                           data=summary_motifs, height=.1)
B_gg + B_gg_se_y + B_gg_se_x

# add an iso-affinity line to the plot
#B_gg_line <- geom_abline(intercept = .01, slope = 1e-8)  
#complete_gg <- B_gg + B_gg_se_y + B_gg_se_x + B_gg_line
#complete_gg

# Save png to disk
ggsave("kon_koff_AB_AZ.eps", width = 6.5, height = 5)
ggsave("kon_koff_AB_q50K_AZ.eps", width = 6.5, height = 5)
ggsave("kon_koff_AB_K50q_AZ.eps", width = 6.5, height = 5)

ggsave("kon_koff_ZM_14vs15compare.eps", width = 6.5, height = 5)
ggsave("kon_koff_BM_14vs15compare.eps", width = 6.5, height = 5)


rm(summary_motifs)
rm(gg)


