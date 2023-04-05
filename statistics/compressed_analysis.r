#install.packages('dplyr', repos='http://cran.us.r-project.org')
#install.packages('reshape2', repos='http://cran.us.r-project.org')
#install.packages('RColorBrewer', repos='http://cran.us.r-project.org')
#install.packages('ggplot2', repos='http://cran.us.r-project.org')
#install.packages('gridExtra', repos='http://cran.us.r-project.org')
#install.packages('stringr', repos='http://cran.us.r-project.org')

library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(data.table)
library(gridExtra)

print("load data")
rawData <- read.table("result_compressed.csv", header=TRUE, sep = ";")

print(paste("Read", length(rawData), " cases"))

my_binwidth <- 4
histogramFull <- ggplot(rawData, aes(x=PIRCHE_II)) +
  geom_histogram(binwidth=my_binwidth, colour="black", fill="white") +
  geom_density(aes(y = ..density.. * (nrow(rawData) * my_binwidth)), col=2, alpha=.2, fill="#FF6666") +
  ylab("count") +
  geom_vline(aes(xintercept=mean(PIRCHE_II)), color="green", linetype="dashed", size=1) +
  ggtitle("PIRCHE-II sum histogram") +
  theme_light()

histogramA <- ggplot(rawData, aes(x=PIRCHE_II_A)) +
  geom_histogram(binwidth=my_binwidth, colour="black", fill="white") +
  geom_density(aes(y = ..density.. * (nrow(rawData) * my_binwidth)), col=2, alpha=.2, fill="#FF6666") +
  ylab("count") +
  geom_vline(aes(xintercept=mean(PIRCHE_II_A)), color="green", linetype="dashed", size=1) +
  ggtitle("PIRCHE-DRB1-A histogram") +
  theme_light()

histogramB <- ggplot(rawData, aes(x=PIRCHE_II_B)) +
  geom_histogram(binwidth=my_binwidth, colour="black", fill="white") +
  geom_density(aes(y = ..density.. * (nrow(rawData) * my_binwidth)), col=2, alpha=.2, fill="#FF6666") +
  ylab("count") +
  geom_vline(aes(xintercept=mean(PIRCHE_II_B)), color="green", linetype="dashed", size=1) +
  ggtitle("PIRCHE-DRB1-B histogram")+
  theme_light()

histogramC <- ggplot(rawData, aes(x=PIRCHE_II_C)) +
  geom_histogram(binwidth=my_binwidth, colour="black", fill="white") +
  geom_density(aes(y = ..density.. * (nrow(rawData) * my_binwidth)), col=2, alpha=.2, fill="#FF6666") +
  ylab("count") +
  geom_vline(aes(xintercept=mean(PIRCHE_II_C)), color="green", linetype="dashed", size=1) +
  ggtitle("PIRCHE-DRB1-C histogram")+
  theme_light()

histogramDRB1 <- ggplot(rawData, aes(x=PIRCHE_II_DRB1)) +
  geom_histogram(binwidth=my_binwidth, colour="black", fill="white") +
  geom_density(aes(y = ..density.. * (nrow(rawData) * my_binwidth)), col=2, alpha=.2, fill="#FF6666") +
  ylab("count") +
  geom_vline(aes(xintercept=mean(PIRCHE_II_DRB1)), color="green", linetype="dashed", size=1) +
  ggtitle("PIRCHE-DRB1-DRB1 histogram")+
  theme_light()

histogramDQB1 <- ggplot(rawData, aes(x=PIRCHE_II_DQB1)) +
  geom_histogram(binwidth=my_binwidth, colour="black", fill="white") +
  geom_density(aes(y = ..density.. * (nrow(rawData) * my_binwidth)), col=2, alpha=.2, fill="#FF6666") +
  ylab("count") +
  geom_vline(aes(xintercept=mean(PIRCHE_II_DQB1)), color="green", linetype="dashed", size=1) +
  ggtitle("PIRCHE-DRB1-DQB1 histogram")+
  theme_light()


grob <- arrangeGrob(histogramFull, histogramA, histogramB, histogramC, histogramDRB1, histogramDQB1,
                    ncol = 3)
ggsave("histograms_pirche.png", plot = grob, dpi = 300, width = 30, height = 20, units = "cm", limitsize = FALSE)
