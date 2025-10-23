rm(list=ls())
library(readxl)
library(ggplot2)
library(tidyr)

# Set working directory and file path
setwd("//ad.monash.edu/home/User042/unig0001/Documents/Quasi-Periodic GP/Journal_Codes/supplementary_material/simulation_studies")
file_path <- "p100.xlsx"

# Read data from the sheet 'w'
w <- read_excel(file_path, sheet = "w")

data_long <- pivot_longer(w, cols = everything(), names_to = "k", values_to = "Value")
data_long$k <- factor(data_long$k, levels = colnames(w))
y_values <- as.numeric(w[1002,])

line_values <- data.frame(
  k = factor(colnames(w), levels = colnames(w)),
  y_value = y_values,
  x_min = as.numeric(factor(colnames(w), levels = colnames(w))) - 0.4,  # Left end of line within the boxplot
  x_max = as.numeric(factor(colnames(w), levels = colnames(w))) + 0.4   # Right end of line within the boxplot
)
ggplot(data_long, aes(x = k, y = Value)) +
  geom_boxplot(fill = "grey", color = "black") +
  geom_segment(data = line_values, aes(x = x_min, xend = x_max, y = y_value, yend = y_value), color = "red", size=1.2,linetype="dashed") +
  labs(
       x = "Sample Size (n)",
       y = "Bootstrap Standard Errors") +
  theme_minimal()+
  theme(
    axis.title = element_text(size = 16),   # axis labels
    axis.text = element_text(size = 14)     # axis tick labels
  )


th <- read_excel(file_path, sheet = "th")

data_long <- pivot_longer(th, cols = everything(), names_to = "k", values_to = "Value")
data_long$k <- factor(data_long$k, levels = colnames(th))
y_values <- as.numeric(th[1002,])

line_values <- data.frame(
  k = factor(colnames(th), levels = colnames(th)),
  y_value = y_values,
  x_min = as.numeric(factor(colnames(th), levels = colnames(th))) - 0.4,  # Left end of line within the boxplot
  x_max = as.numeric(factor(colnames(th), levels = colnames(th))) + 0.4   # Right end of line within the boxplot
)
ggplot(data_long, aes(x = k, y = Value)) +
  geom_boxplot(fill = "grey", color = "black") +
  geom_segment(data = line_values, aes(x = x_min, xend = x_max, y = y_value, yend = y_value), color = "red", size=1.2,linetype="dashed") +
  labs(
       x = "Sample Size (n)",
       y = "Bootstrap Standard Errors") +
  theme_minimal()+
  theme(
    axis.title = element_text(size = 16),   # axis labels
    axis.text = element_text(size = 14)     # axis tick labels
  )


sig2 <- read_excel(file_path, sheet = "sig2")

data_long <- pivot_longer(sig2, cols = everything(), names_to = "k", values_to = "Value")
data_long$k <- factor(data_long$k, levels = colnames(sig2))
y_values <- as.numeric(sig2[1002,])

line_values <- data.frame(
  k = factor(colnames(sig2), levels = colnames(sig2)),
  y_value = y_values,
  x_min = as.numeric(factor(colnames(sig2), levels = colnames(sig2))) - 0.4,  # Left end of line within the boxplot
  x_max = as.numeric(factor(colnames(sig2), levels = colnames(sig2))) + 0.4   # Right end of line within the boxplot
)
ggplot(data_long, aes(x = k, y = Value)) +
  geom_boxplot(fill = "grey", color = "black") +
  geom_segment(data = line_values, aes(x = x_min, xend = x_max, y = y_value, yend = y_value), color = "red", size=1.2,linetype="dashed") +
  labs(
       x = "Sample Size (n)",
       y = "Bootstrap Standard Errors") +
  theme_minimal()+
  theme(
    axis.title = element_text(size = 16),   # axis labels
    axis.text = element_text(size = 14)     # axis tick labels
  )

