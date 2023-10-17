# install.packages("doBy")
# install.packages("reshape2")
# install.packages("patchwork")
# install.packages("ggpubr")
library(reshape2)
library(doBy)
library(ggplot2)
library(patchwork)



plant_all_stat1 <- read.csv("/Users/yan/Desktop/eurecom_code/motif-search/script/draw-fig/6cov.csv", sep=";")
plant_all_stat2 <- read.csv("/Users/yan/Desktop/eurecom_code/motif-search/script/draw-fig/13cov.csv", sep=";")
plant_all_stat3 <- read.csv("/Users/yan/Desktop/eurecom_code/motif-search/script/draw-fig/20cov.csv", sep=";")
plant_all_stat4 <- read.csv("/Users/yan/Desktop/eurecom_code/motif-search/script/draw-fig/27cov.csv", sep=";")
plant_all_stat5 <- read.csv("/Users/yan/Desktop/eurecom_code/motif-search/script/draw-fig/34cov.csv", sep=";")


#首先将根长数据转化为负值，便于作图
plant_all_stat1[which(plant_all_stat1$variable == 'FP'), c('value.mean', 'value.sd')] <- plant_all_stat1[which(plant_all_stat1$variable == 'FP'), c('value.mean', 'value.sd')] * -1
plant_all_stat2[which(plant_all_stat2$variable == 'FP'), c('value.mean', 'value.sd')] <- plant_all_stat2[which(plant_all_stat2$variable == 'FP'), c('value.mean', 'value.sd')] * -1
plant_all_stat3[which(plant_all_stat3$variable == 'FP'), c('value.mean', 'value.sd')] <- plant_all_stat3[which(plant_all_stat3$variable == 'FP'), c('value.mean', 'value.sd')] * -1
plant_all_stat4[which(plant_all_stat4$variable == 'FP'), c('value.mean', 'value.sd')] <- plant_all_stat4[which(plant_all_stat4$variable == 'FP'), c('value.mean', 'value.sd')] * -1
plant_all_stat5[which(plant_all_stat5$variable == 'FP'), c('value.mean', 'value.sd')] <- plant_all_stat5[which(plant_all_stat5$variable == 'FP'), c('value.mean', 'value.sd')] * -1
 
#ggplot2 作图
require(scales)

p1 <- ggplot(plant_all_stat1, aes(days, value.mean, fill = plant)) +
geom_col(position = position_dodge(width = 0.5), width = 0.5, size = 0.3, colour = 'black') +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(), axis.text.x = element_text(angle=90)) +
labs(x = '6 cov', y = 'Nb. of oligos') +
geom_hline(yintercept = 0, size = 0.3) +
scale_y_continuous(breaks = seq(-64, 256, 32), labels = as.character(abs(seq(-64, 256, 32))), limits = c(-64, 256)) +
annotate('text', label = 'TP', 1, 250) +
annotate('text', label = 'FP', 1, -50)
 
p2 <- ggplot(plant_all_stat2, aes(days, value.mean, fill = plant)) +
geom_col(position = position_dodge(width = 0.5), width = 0.5, size = 0.3, colour = 'black') +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(), axis.text.x = element_text(angle=90)) +
labs(x = '13 cov', y = 'Nb. of oligos') +
geom_hline(yintercept = 0, size = 0.3) +
scale_y_continuous(breaks = seq(-64, 256, 32), labels = as.character(abs(seq(-64, 256, 32))), limits = c(-64, 256)) 
# annotate('text', label = 'TP', 1, 250) +
# annotate('text', label = 'FP', 1, -50)

p3 <- ggplot(plant_all_stat3, aes(days, value.mean, fill = plant)) +
geom_col(position = position_dodge(width = 0.5), width = 0.5, size = 0.3, colour = 'black') +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(), axis.text.x = element_text(angle=90)) +
labs(x = '20 cov', y = 'Nb. of oligos') +
geom_hline(yintercept = 0, size = 0.3) +
scale_y_continuous(breaks = seq(-64, 256, 32), labels = as.character(abs(seq(-64, 256, 32))), limits = c(-64, 256)) 
# annotate('text', label = 'TP', 1, 250) +
# annotate('text', label = 'FP', 1, -50)

p4 <- ggplot(plant_all_stat4, aes(days, value.mean, fill = plant)) +
geom_col(position = position_dodge(width = 0.5), width = 0.5, size = 0.3, colour = 'black') +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(), axis.text.x = element_text(angle=90)) +
labs(x = '27 cov', y = 'Nb. of oligos') +
geom_hline(yintercept = 0, size = 0.3) +
scale_y_continuous(breaks = seq(-128, 256, 32), labels = as.character(abs(seq(-128, 256, 32))), limits = c(-128, 256))
# annotate('text', label = 'TP', 1, 250) +
# annotate('text', label = 'FP', 1, -100)

p5 <- ggplot(plant_all_stat5, aes(days, value.mean, fill = plant)) +
geom_col(position = position_dodge(width = 0.5), width = 0.5, size = 0.3, colour = 'black') +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(), axis.text.x = element_text(angle=90)) +
labs(x = '34 cov', y = 'Nb. of oligos') +
geom_hline(yintercept = 0, size = 0.3) +
scale_y_continuous(breaks = seq(-160, 256, 32), labels = as.character(abs(seq(-160, 256, 32))), limits = c(-160, 256))
# annotate('text', label = 'TP', 1, 260) +
# annotate('text', label = 'FP', 1, -150)
# # annotate("text", x = 12.5, y = 3.5, label = "Arbitrary text") +
# #     coord_cartesian(ylim = c(4, 8), clip = "off")

# p <- p1|p2|p3+p4+p5
# ggsave('/Users/yan/Desktop/eurecom_code/data/motif/boa-expand.pdf', p, width = 15, height = 3)

p <- p1 + p2 + p3 + p4+p5+ guide_area() + plot_layout(ncol = 3, guides = "collect")
ggsave('/Users/yan/Desktop/eurecom_code/data/motif/boa-expand.pdf', p, width = 9, height = 6)

# geom_text(aes(label=value.mean), position=position_dodge(width=0.9), vjust=-0.25) +
# scale_x_continuous(trans = log2_trans(),
#     breaks = trans_breaks("log2", function(x) 2^x)(c(64, 1073741824)),
#     labels = trans_format("log2", math_format(2^.x)),
#     limits = c(64, 1073741824)
#     ) + 
