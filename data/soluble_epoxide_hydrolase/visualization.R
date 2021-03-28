#!/usr/bin/env Rscript

# PLIP ANALYZER - VISUALIZATION
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

library(ggplot2)

df = read.csv2("results/result_sorted.csv", sep = ";", dec = ".")

ggplot(df, aes(x=reorder(Interaction, -Frequency), y=Frequency)) +
  geom_bar(stat="identity", color="black", fill="dodgerblue3", width = 0.7) +
  ggtitle("Number of Interaction Occurences in Soluble Epoxide Hydrolase") +
  xlab("Interaction") +
  ylab("Number of Occurences") +
  geom_text(aes(label=Frequency), vjust=-0.3, size=3.0) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

df = df[which(df["Frequency"]>9),]

ggplot(df, aes(x=reorder(Interaction, -Frequency), y=Frequency)) +
  geom_bar(stat="identity", color="black", fill="dodgerblue3", width = 0.7) +
  ggtitle("Number of Interaction Occurences in Soluble Epoxide Hydrolase") +
  xlab("Interaction") +
  ylab("Number of Occurences") +
  geom_text(aes(label=Frequency), vjust=-0.3, size=3.0) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

df = read.csv2("results/result_hoh_sorted.csv", sep = ";", dec = ".")

ggplot(df, aes(x=reorder(Interaction, -Frequency), y=Frequency)) +
  geom_bar(stat="identity", color="black", fill="dodgerblue3", width = 0.7) +
  ggtitle("Number of Interaction Occurences in Soluble Epoxide Hydrolase") +
  xlab("Interaction") +
  ylab("Number of Occurences") +
  geom_text(aes(label=Frequency), vjust=-0.3, size=3.0) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

df = df[which(df["Frequency"]>9),]

ggplot(df, aes(x=reorder(Interaction, -Frequency), y=Frequency)) +
  geom_bar(stat="identity", color="black", fill="dodgerblue3", width = 0.7) +
  ggtitle("Number of Interaction Occurences in Soluble Epoxide Hydrolase") +
  xlab("Interaction") +
  ylab("Number of Occurences") +
  geom_text(aes(label=Frequency), vjust=-0.3, size=3.0) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
