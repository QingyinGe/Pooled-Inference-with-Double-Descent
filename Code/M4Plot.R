library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyverse)
source("Code/func.R")

ylim_upper = 10
ylim_lower = 0.1


theme_slides = theme(text=element_text(size=15),
                     legend.position = "top")


load("Result/MonthlyFrom1To68.Rdata")
ratio_all = data.frame(ratio_mat)

load("Result/DailyFrom1To2.Rdata")
ratio_all = data.frame(rbind(ratio_all,
                             ratio_mat))

load("Result/HourlyFrom1To1.Rdata")
ratio_all = data.frame(rbind(ratio_all,
                             ratio_mat))

separate_methods = c("S:EW", "S:Linear", "S:Cor", "S:Var")

colnames(ratio_all)[1:7] = c("P:Linear", "P:PCR80", "P:PCR70", separate_methods)

df4plot_long = pivot_longer(ratio_all, cols = `P:Linear`:`S:Var`,
                            values_to = "Ratio", names_to = "Method")
df4plot_long$Method = factor(df4plot_long$Method, 
                             levels = c("P:Linear", "P:PCR80", "P:PCR70", "S:EW", "S:Var", "S:Cor", "S:Linear"))

# Figure 2
plot_monthly1 = ggplot(df4plot_long%>% filter(Method == "S:EW", dataset == "Monthly1"), aes(x = Method, y = Ratio)) +
  geom_boxplot(outlier.alpha = 0.05, width = 0.3) + coord_cartesian(ylim = c(ylim_lower,ylim_upper)) + geom_hline(yintercept = 1, col = 'red') + 
  theme_bw() + ylab("Loss Ratio (EW/Ours)") + xlab("") + theme_slides + #theme(text=element_text(size=25), legend.position = "top") + 
  scale_x_discrete(breaks = NULL) + 
  # scale_y_continuous(breaks = 0:7)  + # Removes x-axis labels and ticks
  scale_y_log10() + 
  annotation_logticks(sides = "l")

print(plot_monthly1)

pdf("Figures/Figure2.pdf")
print(plot_monthly1)
dev.off()

# Figure 4
plot1 = ggplot(df4plot_long %>% filter(Method %in% separate_methods & dataset != 'Hourly1'), aes(x = Method, y = Ratio)) +
  geom_boxplot(outlier.alpha = 0.05) + coord_cartesian(ylim = c(ylim_lower,ylim_upper)) + geom_hline(yintercept = 1, col = 'red') + 
  theme_bw() + theme_slides + ylab("Loss Ratio (Benchmarks/P:Ridgeless)") + xlab("Benchmarks") + 
  # scale_y_continuous(breaks = 0:7) + #ggtitle("Overall") 
  scale_y_log10() + 
  annotation_logticks(sides = "l")

print(plot1)

pdf("Figures/Figure4.pdf")
print(plot1)
dev.off()

# Figure EC.?
plot1 = ggplot(df4plot_long %>% filter(Method %in% separate_methods & dataset == 'Hourly1'), aes(x = Method, y = Ratio)) +
  geom_boxplot(outlier.alpha = 0.05) + coord_cartesian(ylim = c(ylim_lower,ylim_upper)) + geom_hline(yintercept = 1, col = 'red') + 
  theme_bw() + theme_slides + ylab("Loss Ratio (Benchmarks/P:Ridgeless)") + xlab("Benchmarks") + 
  # scale_y_continuous(breaks = 0:7) + #ggtitle("Overall") 
  scale_y_log10() + 
  annotation_logticks(sides = "l")

print(plot1)

pdf("Figures/Hourly1.pdf")
print(plot1)
dev.off()

# Figure 5
df4plot_top8 = df4plot_long %>% filter(dataset %in% c(paste0("Monthly",1:7), "Daily1"))

df4plot_top8$dataset = factor(df4plot_top8$dataset,
                              levels = c("Monthly1", "Monthly2", "Monthly3", "Monthly4", "Daily1", "Monthly5", "Monthly6", "Monthly7"),
                              labels = c("(a) Monthly1, #Pooled = 6782", "(b) Monthly2, #Pooled = 6258",
                                         "(c) Monthly3, #Pooled = 2643", "(d) Monthly4, #Pooled = 1517",
                                         "(e) Daily1, #Pooled = 1506", "(f) Monthly5, #Pooled = 1285",
                                         "(g) Monthly6, #Pooled = 1225", "(h) Monthly7, #Pooled = 1048"))

plot_top8 = ggplot(df4plot_top8 %>% filter(Method %in% separate_methods), aes(x = Method, y = Ratio)) +
  geom_boxplot(outlier.alpha = 0.2) + coord_cartesian(ylim = c(ylim_lower,ylim_upper)) + geom_hline(yintercept = 1, col = 'red') + 
  theme_bw() + theme(legend.position = "top") + theme_slides + ylab("Loss Ratio (Benchmarks/P:Ridgeless)") + xlab("Benchmarks") + 
  # scale_y_continuous(breaks = 0:7)
  scale_y_log10() + 
  annotation_logticks(sides = "l")

print(plot_top8 + facet_wrap(~ dataset, nrow = 4))

pdf("Figures/Figure5.pdf")
print(plot_top8 + facet_wrap(~ dataset, nrow = 4))
dev.off()

# Figure EC1: P:Linear vs P:Ridgeless (Overall, all datasets)
plot_linear_overall = ggplot(df4plot_long %>% filter(Method %in% c("P:Linear", "P:PCR80", "P:PCR70") & dataset != 'Hourly1'),
                             aes(x = Method, y = Ratio)) +
  geom_boxplot(outlier.alpha = 0.05, width = 0.3) + coord_cartesian(ylim = c(ylim_lower,ylim_upper)) + geom_hline(yintercept = 1, col = 'red') +
  theme_bw() + ylab("Loss Ratio (Pooled Methods/P:Ridgeless)") + theme_slides +
  scale_x_discrete() +
  scale_y_log10() + 
  annotation_logticks(sides = "l")

print(plot_linear_overall)

pdf("Figures/FigureEC1.pdf")
print(plot_linear_overall) 
dev.off()
