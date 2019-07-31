#install.packages('dplyr', repos='http://cran.us.r-project.org')
#install.packages('reshape2', repos='http://cran.us.r-project.org')
#install.packages('RColorBrewer', repos='http://cran.us.r-project.org')
#install.packages('ggplot2', repos='http://cran.us.r-project.org')
#install.packages('gridExtra', repos='http://cran.us.r-project.org')
#install.packages('stringr', repos='http://cran.us.r-project.org')

library(dplyr)
library(reshape2)
library(ggplot2)
library(data.table)

print("load data")
rawData <- fread("pirche_result.csv", header=TRUE, skip=1, sep = ",", fill=TRUE, data.table=FALSE)

print("filter")
# only remain the donor rows, as the results are in these rows
rawData <- rawData[!is.na(rawData$PIRCHE_II),]

# rename columns to avoid trouble down the road
names(rawData)[names(rawData) == "A*"] <- c("A*_1", "A*_2")
names(rawData)[names(rawData) == "B*"] <- c("B*_1", "B*_2")
names(rawData)[names(rawData) == "C*"] <- c("C*_1", "C*_2")
names(rawData)[names(rawData) == "DRB1*"] <- c("DRB1*_1", "DRB1*_2")
names(rawData)[names(rawData) == "DRB3*"] <- c("DRB3*_1", "DRB3*_2")
names(rawData)[names(rawData) == "DRB4*"] <- c("DRB4*_1", "DRB4*_2")
names(rawData)[names(rawData) == "DRB5*"] <- c("DRB5*_1", "DRB5*_2")
names(rawData)[names(rawData) == "DQA1*"] <- c("DQA1*_1", "DQA1*_2")
names(rawData)[names(rawData) == "DQB1*"] <- c("DQB1*_1", "DQB1*_2")
names(rawData)[names(rawData) == "DPA1*"] <- c("DPA1*_1", "DPA1*_2")
names(rawData)[names(rawData) == "DPB1*"] <- c("DPB1*_1", "DPB1*_2")
names(rawData)[names(rawData) == "MICA*"] <- c("MICA*_1", "MICA*_2")
names(rawData)[names(rawData) == "MICB*"] <- c("MICB*_1", "MICB*_2")

# sum up the peptides in peptide detail columns
rawData$DRB1_Presents_A_Epitopes[is.na(rawData$DRB1_Presents_A_Epitopes)]<-""
rawData$DRB1_Presents_A_Epitopes_count <- str_count(rawData$DRB1_Presents_A_Epitopes, "\\(")
rawData$DRB1_Presents_B_Epitopes[is.na(rawData$DRB1_Presents_B_Epitopes)]<-""
rawData$DRB1_Presents_B_Epitopes_count <- str_count(rawData$DRB1_Presents_B_Epitopes, "\\(")
rawData$DRB1_Presents_C_Epitopes[is.na(rawData$DRB1_Presents_C_Epitopes)]<-""
rawData$DRB1_Presents_C_Epitopes_count <- str_count(rawData$DRB1_Presents_C_Epitopes, "\\(")
rawData$DRB1_Presents_DRB1_Epitopes[is.na(rawData$DRB1_Presents_DRB1_Epitopes)]<-""
rawData$DRB1_Presents_DRB1_Epitopes_count <- str_count(rawData$DRB1_Presents_DRB1_Epitopes, "\\(")
rawData$DRB1_Presents_DRB3_Epitopes[is.na(rawData$DRB1_Presents_DRB3_Epitopes)]<-""
rawData$DRB1_Presents_DRB3_Epitopes_count <- str_count(rawData$DRB1_Presents_DRB3_Epitopes, "\\(")
rawData$DRB1_Presents_DRB4_Epitopes[is.na(rawData$DRB1_Presents_DRB4_Epitopes)]<-""
rawData$DRB1_Presents_DRB4_Epitopes_count <- str_count(rawData$DRB1_Presents_DRB4_Epitopes, "\\(")
rawData$DRB1_Presents_DRB5_Epitopes[is.na(rawData$DRB1_Presents_DRB5_Epitopes)]<-""
rawData$DRB1_Presents_DRB5_Epitopes_count <- str_count(rawData$DRB1_Presents_DRB5_Epitopes, "\\(")
rawData$DRB1_Presents_DRB345_Epitopes_count <- apply(dplyr::select(rawData, DRB1_Presents_DRB3_Epitopes_count, DRB1_Presents_DRB4_Epitopes_count, DRB1_Presents_DRB5_Epitopes_count), 1, sum)
rawData$DRB1_Presents_DQA1_Epitopes[is.na(rawData$DRB1_Presents_DQA1_Epitopes)]<-""
rawData$DRB1_Presents_DQA1_Epitopes_count <- str_count(rawData$DRB1_Presents_DQA1_Epitopes, "\\(")
rawData$DRB1_Presents_DQB1_Epitopes[is.na(rawData$DRB1_Presents_DQB1_Epitopes)]<-""
rawData$DRB1_Presents_DQB1_Epitopes_count <- str_count(rawData$DRB1_Presents_DQB1_Epitopes, "\\(")
rawData$DRB1_Presents_DPA1_Epitopes[is.na(rawData$DRB1_Presents_DPA1_Epitopes)]<-""
rawData$DRB1_Presents_DPA1_Epitopes_count <- str_count(rawData$DRB1_Presents_DPA1_Epitopes, "\\(")
rawData$DRB1_Presents_DPB1_Epitopes[is.na(rawData$DRB1_Presents_DPB1_Epitopes)]<-""
rawData$DRB1_Presents_DPB1_Epitopes_count <- str_count(rawData$DRB1_Presents_DPB1_Epitopes, "\\(")

rawData$DRB3_Presents_A_Epitopes[is.na(rawData$DRB3_Presents_A_Epitopes)]<-""
rawData$DRB3_Presents_A_Epitopes_count <- str_count(rawData$DRB3_Presents_A_Epitopes, "\\(")
rawData$DRB3_Presents_B_Epitopes[is.na(rawData$DRB3_Presents_B_Epitopes)]<-""
rawData$DRB3_Presents_B_Epitopes_count <- str_count(rawData$DRB3_Presents_B_Epitopes, "\\(")
rawData$DRB3_Presents_C_Epitopes[is.na(rawData$DRB3_Presents_C_Epitopes)]<-""
rawData$DRB3_Presents_C_Epitopes_count <- str_count(rawData$DRB3_Presents_C_Epitopes, "\\(")
rawData$DRB3_Presents_DRB1_Epitopes[is.na(rawData$DRB3_Presents_DRB1_Epitopes)]<-""
rawData$DRB3_Presents_DRB1_Epitopes_count <- str_count(rawData$DRB3_Presents_DRB1_Epitopes, "\\(")
rawData$DRB3_Presents_DRB3_Epitopes[is.na(rawData$DRB3_Presents_DRB3_Epitopes)]<-""
rawData$DRB3_Presents_DRB3_Epitopes_count <- str_count(rawData$DRB3_Presents_DRB3_Epitopes, "\\(")
rawData$DRB3_Presents_DRB4_Epitopes[is.na(rawData$DRB3_Presents_DRB4_Epitopes)]<-""
rawData$DRB3_Presents_DRB4_Epitopes_count <- str_count(rawData$DRB3_Presents_DRB4_Epitopes, "\\(")
rawData$DRB3_Presents_DRB5_Epitopes[is.na(rawData$DRB3_Presents_DRB5_Epitopes)]<-""
rawData$DRB3_Presents_DRB5_Epitopes_count <- str_count(rawData$DRB3_Presents_DRB5_Epitopes, "\\(")
rawData$DRB3_Presents_DRB345_Epitopes_count <- apply(dplyr::select(rawData, DRB3_Presents_DRB3_Epitopes_count, DRB3_Presents_DRB4_Epitopes_count, DRB3_Presents_DRB5_Epitopes_count), 1, sum)
rawData$DRB3_Presents_DQA1_Epitopes[is.na(rawData$DRB3_Presents_DQA1_Epitopes)]<-""
rawData$DRB3_Presents_DQA1_Epitopes_count <- str_count(rawData$DRB3_Presents_DQA1_Epitopes, "\\(")
rawData$DRB3_Presents_DQB1_Epitopes[is.na(rawData$DRB3_Presents_DQB1_Epitopes)]<-""
rawData$DRB3_Presents_DQB1_Epitopes_count <- str_count(rawData$DRB3_Presents_DQB1_Epitopes, "\\(")
rawData$DRB3_Presents_DPA1_Epitopes[is.na(rawData$DRB3_Presents_DPA1_Epitopes)]<-""
rawData$DRB3_Presents_DPA1_Epitopes_count <- str_count(rawData$DRB3_Presents_DPA1_Epitopes, "\\(")
rawData$DRB3_Presents_DPB1_Epitopes[is.na(rawData$DRB3_Presents_DPB1_Epitopes)]<-""
rawData$DRB3_Presents_DPB1_Epitopes_count <- str_count(rawData$DRB3_Presents_DPB1_Epitopes, "\\(")

rawData$DRB4_Presents_A_Epitopes[is.na(rawData$DRB4_Presents_A_Epitopes)]<-""
rawData$DRB4_Presents_A_Epitopes_count <- str_count(rawData$DRB4_Presents_A_Epitopes, "\\(")
rawData$DRB4_Presents_B_Epitopes[is.na(rawData$DRB4_Presents_B_Epitopes)]<-""
rawData$DRB4_Presents_B_Epitopes_count <- str_count(rawData$DRB4_Presents_B_Epitopes, "\\(")
rawData$DRB4_Presents_C_Epitopes[is.na(rawData$DRB4_Presents_C_Epitopes)]<-""
rawData$DRB4_Presents_C_Epitopes_count <- str_count(rawData$DRB4_Presents_C_Epitopes, "\\(")
rawData$DRB4_Presents_DRB1_Epitopes[is.na(rawData$DRB4_Presents_DRB1_Epitopes)]<-""
rawData$DRB4_Presents_DRB1_Epitopes_count <- str_count(rawData$DRB4_Presents_DRB1_Epitopes, "\\(")
rawData$DRB4_Presents_DRB3_Epitopes[is.na(rawData$DRB4_Presents_DRB3_Epitopes)]<-""
rawData$DRB4_Presents_DRB3_Epitopes_count <- str_count(rawData$DRB4_Presents_DRB3_Epitopes, "\\(")
rawData$DRB4_Presents_DRB4_Epitopes[is.na(rawData$DRB4_Presents_DRB4_Epitopes)]<-""
rawData$DRB4_Presents_DRB4_Epitopes_count <- str_count(rawData$DRB4_Presents_DRB4_Epitopes, "\\(")
rawData$DRB4_Presents_DRB5_Epitopes[is.na(rawData$DRB4_Presents_DRB5_Epitopes)]<-""
rawData$DRB4_Presents_DRB5_Epitopes_count <- str_count(rawData$DRB4_Presents_DRB5_Epitopes, "\\(")
rawData$DRB4_Presents_DRB345_Epitopes_count <- apply(dplyr::select(rawData, DRB4_Presents_DRB3_Epitopes_count, DRB4_Presents_DRB4_Epitopes_count, DRB4_Presents_DRB5_Epitopes_count), 1, sum)
rawData$DRB4_Presents_DQA1_Epitopes[is.na(rawData$DRB4_Presents_DQA1_Epitopes)]<-""
rawData$DRB4_Presents_DQA1_Epitopes_count <- str_count(rawData$DRB4_Presents_DQA1_Epitopes, "\\(")
rawData$DRB4_Presents_DQB1_Epitopes[is.na(rawData$DRB4_Presents_DQB1_Epitopes)]<-""
rawData$DRB4_Presents_DQB1_Epitopes_count <- str_count(rawData$DRB4_Presents_DQB1_Epitopes, "\\(")
rawData$DRB4_Presents_DPA1_Epitopes[is.na(rawData$DRB4_Presents_DPA1_Epitopes)]<-""
rawData$DRB4_Presents_DPA1_Epitopes_count <- str_count(rawData$DRB4_Presents_DPA1_Epitopes, "\\(")
rawData$DRB4_Presents_DPB1_Epitopes[is.na(rawData$DRB4_Presents_DPB1_Epitopes)]<-""
rawData$DRB4_Presents_DPB1_Epitopes_count <- str_count(rawData$DRB4_Presents_DPB1_Epitopes, "\\(")

rawData$DRB5_Presents_A_Epitopes[is.na(rawData$DRB5_Presents_A_Epitopes)]<-""
rawData$DRB5_Presents_A_Epitopes_count <- str_count(rawData$DRB5_Presents_A_Epitopes, "\\(")
rawData$DRB5_Presents_B_Epitopes[is.na(rawData$DRB5_Presents_B_Epitopes)]<-""
rawData$DRB5_Presents_B_Epitopes_count <- str_count(rawData$DRB5_Presents_B_Epitopes, "\\(")
rawData$DRB5_Presents_C_Epitopes[is.na(rawData$DRB5_Presents_C_Epitopes)]<-""
rawData$DRB5_Presents_C_Epitopes_count <- str_count(rawData$DRB5_Presents_C_Epitopes, "\\(")
rawData$DRB5_Presents_DRB1_Epitopes[is.na(rawData$DRB5_Presents_DRB1_Epitopes)]<-""
rawData$DRB5_Presents_DRB1_Epitopes_count <- str_count(rawData$DRB5_Presents_DRB1_Epitopes, "\\(")
rawData$DRB5_Presents_DRB3_Epitopes[is.na(rawData$DRB5_Presents_DRB3_Epitopes)]<-""
rawData$DRB5_Presents_DRB3_Epitopes_count <- str_count(rawData$DRB5_Presents_DRB3_Epitopes, "\\(")
rawData$DRB5_Presents_DRB4_Epitopes[is.na(rawData$DRB5_Presents_DRB4_Epitopes)]<-""
rawData$DRB5_Presents_DRB4_Epitopes_count <- str_count(rawData$DRB5_Presents_DRB4_Epitopes, "\\(")
rawData$DRB5_Presents_DRB5_Epitopes[is.na(rawData$DRB5_Presents_DRB5_Epitopes)]<-""
rawData$DRB5_Presents_DRB5_Epitopes_count <- str_count(rawData$DRB5_Presents_DRB5_Epitopes, "\\(")
rawData$DRB5_Presents_DRB345_Epitopes_count <- apply(dplyr::select(rawData, DRB5_Presents_DRB3_Epitopes_count, DRB5_Presents_DRB4_Epitopes_count, DRB5_Presents_DRB5_Epitopes_count), 1, sum)
rawData$DRB5_Presents_DQA1_Epitopes[is.na(rawData$DRB5_Presents_DQA1_Epitopes)]<-""
rawData$DRB5_Presents_DQA1_Epitopes_count <- str_count(rawData$DRB5_Presents_DQA1_Epitopes, "\\(")
rawData$DRB5_Presents_DQB1_Epitopes[is.na(rawData$DRB5_Presents_DQB1_Epitopes)]<-""
rawData$DRB5_Presents_DQB1_Epitopes_count <- str_count(rawData$DRB5_Presents_DQB1_Epitopes, "\\(")
rawData$DRB5_Presents_DPA1_Epitopes[is.na(rawData$DRB5_Presents_DPA1_Epitopes)]<-""
rawData$DRB5_Presents_DPA1_Epitopes_count <- str_count(rawData$DRB5_Presents_DPA1_Epitopes, "\\(")
rawData$DRB5_Presents_DPB1_Epitopes[is.na(rawData$DRB5_Presents_DPB1_Epitopes)]<-""
rawData$DRB5_Presents_DPB1_Epitopes_count <- str_count(rawData$DRB5_Presents_DPB1_Epitopes, "\\(")

rawData$DRB345_Presents_A_Epitopes_count <- apply(dplyr::select(rawData, DRB3_Presents_A_Epitopes_count, DRB4_Presents_A_Epitopes_count, DRB5_Presents_A_Epitopes_count), 1, sum)
rawData$DRB345_Presents_B_Epitopes_count <- apply(dplyr::select(rawData, DRB3_Presents_B_Epitopes_count, DRB4_Presents_B_Epitopes_count, DRB5_Presents_B_Epitopes_count), 1, sum)
rawData$DRB345_Presents_C_Epitopes_count <- apply(dplyr::select(rawData, DRB3_Presents_C_Epitopes_count, DRB4_Presents_C_Epitopes_count, DRB5_Presents_C_Epitopes_count), 1, sum)
rawData$DRB345_Presents_DRB1_Epitopes_count <- apply(dplyr::select(rawData, DRB3_Presents_DRB1_Epitopes_count, DRB4_Presents_DRB1_Epitopes_count, DRB5_Presents_DRB1_Epitopes_count), 1, sum)
rawData$DRB345_Presents_DRB3_Epitopes_count <- apply(dplyr::select(rawData, DRB3_Presents_DRB3_Epitopes_count, DRB4_Presents_DRB3_Epitopes_count, DRB5_Presents_DRB3_Epitopes_count), 1, sum)
rawData$DRB345_Presents_DRB4_Epitopes_count <- apply(dplyr::select(rawData, DRB3_Presents_DRB4_Epitopes_count, DRB4_Presents_DRB4_Epitopes_count, DRB5_Presents_DRB4_Epitopes_count), 1, sum)
rawData$DRB345_Presents_DRB5_Epitopes_count <- apply(dplyr::select(rawData, DRB3_Presents_DRB3_Epitopes_count, DRB4_Presents_DRB5_Epitopes_count, DRB5_Presents_DRB5_Epitopes_count), 1, sum)
rawData$DRB345_Presents_DRB345_Epitopes_count <- apply(dplyr::select(rawData, DRB345_Presents_DRB3_Epitopes_count, DRB345_Presents_DRB4_Epitopes_count, DRB345_Presents_DRB5_Epitopes_count), 1, sum)
rawData$DRB345_Presents_DQA1_Epitopes_count <- apply(dplyr::select(rawData, DRB3_Presents_DQA1_Epitopes_count, DRB4_Presents_DQA1_Epitopes_count, DRB5_Presents_DQA1_Epitopes_count), 1, sum)
rawData$DRB345_Presents_DQB1_Epitopes_count <- apply(dplyr::select(rawData, DRB3_Presents_DQB1_Epitopes_count, DRB4_Presents_DQB1_Epitopes_count, DRB5_Presents_DQB1_Epitopes_count), 1, sum)
rawData$DRB345_Presents_DPA1_Epitopes_count <- apply(dplyr::select(rawData, DRB3_Presents_DPA1_Epitopes_count, DRB4_Presents_DPA1_Epitopes_count, DRB5_Presents_DPA1_Epitopes_count), 1, sum)
rawData$DRB345_Presents_DPB1_Epitopes_count <- apply(dplyr::select(rawData, DRB3_Presents_DPB1_Epitopes_count, DRB4_Presents_DPB1_Epitopes_count, DRB5_Presents_DPB1_Epitopes_count), 1, sum)

rawData$DQA1_DQB1_Presents_A_Epitopes[is.na(rawData$DQA1_DQB1_Presents_A_Epitopes)]<-""
rawData$DQA1_DQB1_Presents_A_Epitopes_count <- str_count(rawData$DQA1_DQB1_Presents_A_Epitopes, "\\(")
rawData$DQA1_DQB1_Presents_B_Epitopes[is.na(rawData$DQA1_DQB1_Presents_B_Epitopes)]<-""
rawData$DQA1_DQB1_Presents_B_Epitopes_count <- str_count(rawData$DQA1_DQB1_Presents_B_Epitopes, "\\(")
rawData$DQA1_DQB1_Presents_C_Epitopes[is.na(rawData$DQA1_DQB1_Presents_C_Epitopes)]<-""
rawData$DQA1_DQB1_Presents_C_Epitopes_count <- str_count(rawData$DQA1_DQB1_Presents_C_Epitopes, "\\(")
rawData$DQA1_DQB1_Presents_DRB1_Epitopes[is.na(rawData$DQA1_DQB1_Presents_DRB1_Epitopes)]<-""
rawData$DQA1_DQB1_Presents_DRB1_Epitopes_count <- str_count(rawData$DQA1_DQB1_Presents_DRB1_Epitopes, "\\(")
rawData$DQA1_DQB1_Presents_DRB3_Epitopes[is.na(rawData$DQA1_DQB1_Presents_DRB3_Epitopes)]<-""
rawData$DQA1_DQB1_Presents_DRB3_Epitopes_count <- str_count(rawData$DQA1_DQB1_Presents_DRB3_Epitopes, "\\(")
rawData$DQA1_DQB1_Presents_DRB4_Epitopes[is.na(rawData$DQA1_DQB1_Presents_DRB4_Epitopes)]<-""
rawData$DQA1_DQB1_Presents_DRB4_Epitopes_count <- str_count(rawData$DQA1_DQB1_Presents_DRB4_Epitopes, "\\(")
rawData$DQA1_DQB1_Presents_DRB5_Epitopes[is.na(rawData$DQA1_DQB1_Presents_DRB5_Epitopes)]<-""
rawData$DQA1_DQB1_Presents_DRB5_Epitopes_count <- str_count(rawData$DQA1_DQB1_Presents_DRB5_Epitopes, "\\(")
rawData$DQA1_DQB1_Presents_DRB345_Epitopes_count <- apply(dplyr::select(rawData, DQA1_DQB1_Presents_DRB3_Epitopes_count, DQA1_DQB1_Presents_DRB4_Epitopes_count, DQA1_DQB1_Presents_DRB5_Epitopes_count), 1, sum)
rawData$DQA1_DQB1_Presents_DQA1_Epitopes[is.na(rawData$DQA1_DQB1_Presents_DQA1_Epitopes)]<-""
rawData$DQA1_DQB1_Presents_DQA1_Epitopes_count <- str_count(rawData$DQA1_DQB1_Presents_DQA1_Epitopes, "\\(")
rawData$DQA1_DQB1_Presents_DQB1_Epitopes[is.na(rawData$DQA1_DQB1_Presents_DQB1_Epitopes)]<-""
rawData$DQA1_DQB1_Presents_DQB1_Epitopes_count <- str_count(rawData$DQA1_DQB1_Presents_DQB1_Epitopes, "\\(")
rawData$DQA1_DQB1_Presents_DPA1_Epitopes[is.na(rawData$DQA1_DQB1_Presents_DPA1_Epitopes)]<-""
rawData$DQA1_DQB1_Presents_DPA1_Epitopes_count <- str_count(rawData$DQA1_DQB1_Presents_DPA1_Epitopes, "\\(")
rawData$DQA1_DQB1_Presents_DPB1_Epitopes[is.na(rawData$DQA1_DQB1_Presents_DPB1_Epitopes)]<-""
rawData$DQA1_DQB1_Presents_DPB1_Epitopes_count <- str_count(rawData$DQA1_DQB1_Presents_DPB1_Epitopes, "\\(")

rawData$DPA1_DPB1_Presents_A_Epitopes[is.na(rawData$DPA1_DPB1_Presents_A_Epitopes)]<-""
rawData$DPA1_DPB1_Presents_A_Epitopes_count <- str_count(rawData$DPA1_DPB1_Presents_A_Epitopes, "\\(")
rawData$DPA1_DPB1_Presents_B_Epitopes[is.na(rawData$DPA1_DPB1_Presents_B_Epitopes)]<-""
rawData$DPA1_DPB1_Presents_B_Epitopes_count <- str_count(rawData$DPA1_DPB1_Presents_B_Epitopes, "\\(")
rawData$DPA1_DPB1_Presents_C_Epitopes[is.na(rawData$DPA1_DPB1_Presents_C_Epitopes)]<-""
rawData$DPA1_DPB1_Presents_C_Epitopes_count <- str_count(rawData$DPA1_DPB1_Presents_C_Epitopes, "\\(")
rawData$DPA1_DPB1_Presents_DRB1_Epitopes[is.na(rawData$DPA1_DPB1_Presents_DRB1_Epitopes)]<-""
rawData$DPA1_DPB1_Presents_DRB1_Epitopes_count <- str_count(rawData$DPA1_DPB1_Presents_DRB1_Epitopes, "\\(")
rawData$DPA1_DPB1_Presents_DRB3_Epitopes[is.na(rawData$DPA1_DPB1_Presents_DRB3_Epitopes)]<-""
rawData$DPA1_DPB1_Presents_DRB3_Epitopes_count <- str_count(rawData$DPA1_DPB1_Presents_DRB3_Epitopes, "\\(")
rawData$DPA1_DPB1_Presents_DRB4_Epitopes[is.na(rawData$DPA1_DPB1_Presents_DRB4_Epitopes)]<-""
rawData$DPA1_DPB1_Presents_DRB4_Epitopes_count <- str_count(rawData$DPA1_DPB1_Presents_DRB4_Epitopes, "\\(")
rawData$DPA1_DPB1_Presents_DRB5_Epitopes[is.na(rawData$DPA1_DPB1_Presents_DRB5_Epitopes)]<-""
rawData$DPA1_DPB1_Presents_DRB5_Epitopes_count <- str_count(rawData$DPA1_DPB1_Presents_DRB5_Epitopes, "\\(")
rawData$DPA1_DPB1_Presents_DRB345_Epitopes_count <- apply(dplyr::select(rawData, DPA1_DPB1_Presents_DRB3_Epitopes_count, DPA1_DPB1_Presents_DRB4_Epitopes_count, DPA1_DPB1_Presents_DRB5_Epitopes_count), 1, sum)
rawData$DPA1_DPB1_Presents_DQA1_Epitopes[is.na(rawData$DPA1_DPB1_Presents_DQA1_Epitopes)]<-""
rawData$DPA1_DPB1_Presents_DQA1_Epitopes_count <- str_count(rawData$DPA1_DPB1_Presents_DQA1_Epitopes, "\\(")
rawData$DPA1_DPB1_Presents_DQB1_Epitopes[is.na(rawData$DPA1_DPB1_Presents_DQB1_Epitopes)]<-""
rawData$DPA1_DPB1_Presents_DQB1_Epitopes_count <- str_count(rawData$DPA1_DPB1_Presents_DQB1_Epitopes, "\\(")
rawData$DPA1_DPB1_Presents_DPA1_Epitopes[is.na(rawData$DPA1_DPB1_Presents_DPA1_Epitopes)]<-""
rawData$DPA1_DPB1_Presents_DPA1_Epitopes_count <- str_count(rawData$DPA1_DPB1_Presents_DPA1_Epitopes, "\\(")
rawData$DPA1_DPB1_Presents_DPB1_Epitopes[is.na(rawData$DPA1_DPB1_Presents_DPB1_Epitopes)]<-""
rawData$DPA1_DPB1_Presents_DPB1_Epitopes_count <- str_count(rawData$DPA1_DPB1_Presents_DPB1_Epitopes, "\\(")

# some basic stats on the data
print("Summary")
print(summary(dplyr::select(rawData, "Patient/Donor_ID",
A_Originated_Epitopes_count, B_Originated_Epitopes_count, C_Originated_Epitopes_count,
DRB1_Originated_Epitopes_count, DRB3_Originated_Epitopes_count, DRB4_Originated_Epitopes_count, DRB5_Originated_Epitopes_count,
DQA1_Originated_Epitopes_count, DQB1_Originated_Epitopes_count,
DPA1_Originated_Epitopes_count, DPB1_Originated_Epitopes_count,
PIRCHE_I, PIRCHE_II)))

print("Dot plot")
dotPlot <- ggplot(data = rawData, aes(x=1, y=PIRCHE_II)) +
stat_summary(fun.y="median", color="red", geom="point", shape=21, size=5, position=position_dodge(width=0.75)) +
stat_summary(fun.y="mean", color="green", geom="point", shape=23, size=5, position=position_dodge(width=0.75)) +
geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, binwidth=2) +
theme_light()

print(dotPlot)