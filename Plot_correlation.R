library("ggplot2")
library("scales")

#save data to data frame
#ct_negrar <- c(35, 24.1, 27.0, 31.8, 29.2, 26.5, 33.4, 40.2, 35.2, 36.8, 37.8)
ct_LCM_N <- c(33.2, 27.2, 30.6, 33.9, 31.6, 28.1, 34.3, 33.4, 30.6, 38.4, 40.5)
ct_LCM_ORF1ab <- c(29.4, 33.4, 27.0, 32.5, 28.4, 26.4, 36.5, NA, 27, 38.5, NA)
num_reads <- c(235406, 169969, 87604, 42324, 41062, 35351, 4337, 2265, 2028, 30, 6)
data <- data.frame(ct = c(ct_LCM_N, ct_LCM_ORF1ab), num_reads_rep = rep(num_reads, 2), Gene = c(rep("N", length(ct_LCM_N)), rep("ORF1ab", length(ct_LCM_ORF1ab))))

#scatterplot
ggplot(data, aes(x = ct, y = num_reads_rep, color = Gene)) + 
  geom_point(size = 3) +
  geom_smooth(method = lm, aes(fill = Gene)) +
  scale_x_continuous(limits = c(25, 40)) +
  scale_y_continuous(trans=log10_trans(),
  breaks = trans_breaks("log10", function(x) 10^x),
  labels = trans_format("log10", math_format(10^.x))) +
  labs(x="Ct", y = "Num. reads") +
  theme(legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16))

ggsave(filename = "Ct_numReads.png", device = "png")

#check for normality distribution (H0: data follow a normal distribution)
shapiro.test(log10(num_reads))
shapiro.test(ct_LCM_N)
shapiro.test(ct_LCM_ORF1ab)

#evaluate strength of correlation using pearson method
cor_method <- "pearson"

correlation_N <- cor.test(ct_LCM_N, log10(num_reads), method = cor_method, alternative = "two.sided")
R2_N <- correlation_N$estimate^2

correlation_ORF1ab <- cor.test(ct_LCM_ORF1ab, log10(num_reads), method = cor_method, alternative = "two.sided")
R2_ORF1ab <- correlation_ORF1ab$estimate^2
