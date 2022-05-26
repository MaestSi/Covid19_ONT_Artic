library("ggplot2")
library("scales")
library("pheatmap")

#save data to data frame
#ct_negrar <- c(35, 24.1, 27.0, 31.8, 29.2, 26.5, 33.4, 40.2, 35.2, 36.8, 37.8)
ct_LCM_N <- c(33.2, 27.2, 30.6, 33.9, 31.6, 28.1, 34.3, 33.4, 30.6, 38.4, 40.5)
ct_LCM_ORF1ab <- c(29.4, 33.4, 27.0, 32.5, 28.4, 26.4, 36.5, NA, 27, 38.5, NA)
num_reads <- c(235406, 169969, 87604, 42324, 41062, 35351, 4337, 2265, 2028, 30, 6)
ct_PC <- c(11.3, 11.3, 17.3, 17.3, 23.7, 23.7, 30.8, 30.8)
num_copies_PC <- c(10**9, 10**9, 10**7, 10**7, 10**5, 10**5, 10**3, 10**3)
num_reads_PC <- c(1799, 1943, 44, 22, 9, 5, 8, 7)
#data <- data.frame(ct = c(ct_LCM_N, ct_LCM_ORF1ab), num_reads_rep = rep(num_reads, 2), Gene = c(rep("N", length(ct_LCM_N)), rep("ORF1ab", length(ct_LCM_ORF1ab))))

data <- data.frame(ct = ct_LCM_N, num_reads_rep = num_reads, Gene = rep("N", length(ct_LCM_N)))
data_PC <- data.frame(ct = ct_PC, num_reads_rep = num_reads_PC, Amplicon = rep("Ampl. 96", length(num_copies_PC)))

#scatterplot
ggplot(data, aes(x = ct, y = num_reads_rep, color = Gene)) + 
  geom_point(size = 3) +
  geom_smooth(method = lm, aes(fill = Gene)) +
  scale_x_continuous(limits = c(25, 42)) +
  scale_y_continuous(trans=log10_trans(),
  breaks = trans_breaks("log10", function(x) 10^x),
  labels = trans_format("log10", math_format(10^.x))) +
  labs(x="RT-qPCR (Ct)", y = "Sequencing (Reads number)") +
  theme(legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16))

ggsave(filename = "Ct_numReads.png", device = "png")

#scatterplot PC
ggplot(data_PC, aes(x = ct, y = num_reads_rep, color = Amplicon)) + 
  geom_point(size = 3) +
  geom_smooth(method = lm, aes(fill = Amplicon)) +
  scale_x_continuous(limits = c(10, 32)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  labs(x="RT-qPCR (Ct)", y = "Sequencing (Reads number)") +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave(filename = "Ct_numReads_PC.png", device = "png")

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

#check for normality distribution (H0: data follow a normal distribution)
shapiro.test(log10(num_reads_PC))
shapiro.test(ct_PC)

#evaluate strength of correlation using pearson method
cor_method <- "pearson"

correlation_PC <- cor.test(ct_PC, log10(num_reads_PC), method = cor_method, alternative = "two.sided")
R2_PC <- correlation_PC$estimate^2

#populate matrix with presence/absence/not genotypability of variants
variants <- matrix(data = 0, nrow = 26, ncol = 10)
rownames(variants) <- c("L18F", "T20N", "P26S", "Q52R", "del21765:6", "D80A", "D138Y", "del21991:3", "R190S", "D215G", "K417N", "N439K", "E484K", "N501Y", "A570D", "D614G", "H655Y", "Q677H", "P681H", "A701V", "T716I", "S982A", "T1027I", "Y73C", "E92K", "D3L")
colnames(variants) <- c("sample_128", "sample_80", "sample_326", "sample_18", "sample_123", "sample_241", "sample_41", "sample_331", "sample_325", "sample_282")

#excluded "sample_172", since no variant was genotypable
#variants[c("L18F", "T20N", "P26S", "Q52R", "del21765:6", "D80A", "D138Y", "del21991:3", "R190S", "D215G", "K417N", "N439K", "E484K", "N501Y", "A570D", "D614G", "H655Y", "Q677H", "P681H", "A701V", "T716I", "S982A", "T1027I", "Y73C", "E92K", "D3L"), "sample_172"] <- NA

variants[c("L18F", "T20N", "P26S", "Q52R", "del21765:6", "D80A", "D138Y", "del21991:3", "R190S", "D215G", "K417N", "N439K", "E484K", "N501Y", "A570D", "D614G", "S982A", "T1027I", "Y73C", "E92K", "D3L"), "sample_282"] <- NA

variants[c("K417N", "D614G"), "sample_326"] <- 1

variants[c("H655Y", "Q677H", "P681H", "A701V", "T716I"), "sample_241"] <- NA
variants["D614G", "sample_241"] <- 1

variants["D614G", "sample_123"] <- 1

variants[c("D614G", "del21991:3"), "sample_80"] <- 1

variants[c("del21765:6", "del21991:3", "N501Y", "A570D", "D614G", "P681H", "T716I", "S982A", "Y73C", "D3L"), "sample_128"] <- 1

variants[c("L18F", "T20N", "P26S", "Q52R", "del21765:6", "D80A", "D138Y", "del21991:3", "R190S", "D215G", "N439K", "E484K", "N501Y", "A570D", "D614G", "Y73C", "E92K", "D3L"), "sample_41"] <- NA
variants["Q677H", "sample_41"] <- 1

variants[c("L18F", "T20N", "P26S", "R190S", "D215G", "K417N", "A570D", "D614G", "S982A", "T1027I"), "sample_331"] <- NA

variants["D614G", "sample_18"] <- 1

variants[c("L18F", "T20N", "P26S", "Q52R", "del21765:6", "D80A", "D138Y", "del21991:3", "R190S", "D215G", "K417N", "N439K", "E484K", "N501Y", "A570D", "D614G", "H655Y", "Q677H", "P681H", "A701V", "T716I", "S982A", "T1027I"), "sample_325"] <- NA

#plot heatmap
png("STArS_heatmap.png", width = 600, height = 600)
pheatmap(variants, cluster_rows = FALSE, cluster_cols = FALSE, cellwidth = 20, cellheight = 20, legend = FALSE, na_col = "grey", show_rownames = TRUE, show_colnames = TRUE)
dev.off()