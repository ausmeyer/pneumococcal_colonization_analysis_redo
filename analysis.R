## Clear workspace and set working directory
rm(list = ls())
setwd('~/Google Drive/Documents/PostPostDoc/pneumococcal_colonization_analysis_redo/')

## Library Imports
libraries.call <- c("dplyr", "ggplot2", "readxl", "reshape")
lapply(libraries.call, require, character.only = TRUE)

## Utility Functions
generate.matrix <- function(df, cutoff, variable, response) {
  t.positive <- sum(df[[variable]] > cutoff & as.numeric(as.character(df[[response]]) == 1))
  f.positive <- sum(df[[variable]] > cutoff & as.numeric(as.character(df[[response]]) == 0))
  f.negative <- sum(df[[variable]] < cutoff & as.numeric(as.character(df[[response]]) == 1))
  t.negative <- sum(df[[variable]] < cutoff & as.numeric(as.character(df[[response]]) == 0))
  return(data.frame(t.p = t.positive, f.p = f.positive, f.n = f.negative, t.n = t.negative))
}

calculate.test.stats <- function(df) {
  sensitivity <- df$t.p / (df$t.p + df$f.n)
  specificity <- df$t.n / (df$t.n + df$f.p)
  
  return(data.frame(sens = sensitivity, spec = specificity))
}

make.test.stats.df <- function(df, which.variable) {
  df.series <- seq(0, max(df$variable), by = 0.01)
  df.stats <- data.frame()
  invisible(sapply(df.series, function(x) df.stats <<- rbind(df.stats, (calculate.test.stats(generate.matrix(df, cutoff = x, variable = "variable", response = which.variable))))))
  return(df.stats)
}

## Function to left censor the data for multiple plots
## @agm Requires ordering the data after censoring to obtain desired results
left.censor <- function(df) {
  df <- data.frame(x = df$x, y = df$y)
  new.df <- data.frame()
  lapply(1:(nrow(df)-1), function(x) new.df <<- rbind(new.df, data.frame(x = df$x[x+1], y = df$y[x])))
  df <- rbind(df, new.df)
}

import.data <- function(variable){
  df <- read_excel('Database_BMJ_Open_25Jul2014.xlsx', col_types = c('text', 'text', 'numeric', 
                                                                          'numeric', 'numeric', 'text', 
                                                                          'text', 'numeric', 'numeric', 
                                                                          'numeric', 'numeric', 'numeric', 
                                                                          'numeric', 'numeric', 'text', 
                                                                          'text', 'text', 'text', 
                                                                          'text'))
  df[, 1:17]
  df <- data.frame(age = df$age, variable = df[[variable]], bacter = df$bacteremia, pneumo = df$Pneumococcal_diagnosis)
  df <- na.omit(df)
}

make.fit <- function(df) {
  require('glm')
  
  fit <- glm(factor(bacter) ~ variable, data = df, family = "binomial")
  predicted.bacteremia <- predict(fit, type="response", se = TRUE)
  
  fit <- glm(factor(pneumo) ~ variable, data = df, family = "binomial")
  predicted.pneumo <- predict(fit, type="response", se = TRUE)
  
  return(data.frame(bacter = predicted.bacteremia, pneumo = predicted.pneumo))
}

## Plotting Functions
make.roc.plot <- function(df, file.name) {
  df <- data.frame(x = 1 - df$spec, y = df$sens)
  df <- left.censor(df)
  df <- df[order(df$y), ]
  
  p <- ggplot() + 
    geom_point(data = df, aes(x = x, y = y), size = 1) +
    geom_line(data = df, aes(x = x, y = y)) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("1 - Specificity") +
    ylab("Sensitivity") +
    theme_bw()
  ggsave(plot = p, filename = paste("ROC_", file.name, sep = ''), height = 5, width = 6)
  return(p)
}

make.logistic.plot <- function(df, file.name) {
  fit.df <- make.fit(df)
  new.df <- data.frame(variable = df$variable, bacter = df$bacter, prob.bacter = fit.df$bacter.fit, pneumo = df$pneumo, prob.pneumo = fit.df$pneumo.fit)
  
  p <- ggplot() + 
    geom_point(data = new.df, aes(x = variable, y = as.numeric(as.character(pneumo)), color = "Diagnosis Status")) +
    geom_point(data = new.df, aes(x = variable, y = as.numeric(as.character(prob.pneumo)), color = 'Logistic Fit')) +
    geom_line(data = new.df, aes(x = variable, y = as.numeric(as.character(prob.pneumo)), color = 'Logistic Fit')) +
    xlab("LytA Concentration (units)") +
    ylab("Pneumococcal pneumonia diagnosis") +
    scale_color_discrete(name = "") +
    theme_bw()
  ggsave(plot = p, filename = paste("Logistic_", file.name, sep = ''), height = 5, width = 7)
  return(p)
}

lytA.filename <- 'lytA.pdf'
raw.lytA <- import.data(variable = 'lytA_NP')
make.logistic.plot(df = raw.lytA, file.name = lytA.filename)
roc.data.lytA <- make.test.stats.df(df = raw.lytA, which.variable = "pneumo")
make.roc.plot(df = roc.data.lytA, file.name = lytA.filename)

pct.filename <- 'pct.pdf'
raw.pct <- import.data(variable = 'PCT')
make.logistic.plot(df = raw.pct, file.name = pct.filename)
roc.data.pct <- make.test.stats.df(df = raw.pct, which.variable = "pneumo")
make.roc.plot(df = roc.data.pct, file.name = pct.filename)

# 
# reform <- data.frame(lytA = lytA.series, sens = test.stats$sens, spec = test.stats$spec)
# reform.melted <- melt(reform, id = c("lytA"))
# colnames(reform.melted) <- c('x', 'var', 'y')
# 
# reform.sens <- left.censor(filter(reform.melted, var == 'sens'))
# reform.sens <- cbind(reform.sens, var = rep('sens', nrow(reform.sens)))
# reform.sens <- reform.sens[rev(order(reform.sens$y)), ]
# 
# reform.spec <- left.censor(filter(reform.melted, var == 'spec'))
# reform.spec <- cbind(reform.spec, var = rep('spec', nrow(reform.spec)))
# reform.spec <- reform.spec[order(reform.spec$y), ]
# reform.censored <- rbind(reform.sens, reform.spec)
# 
# p <- ggplot() + 
#   geom_point(data = reform.censored, aes(x = x, y = y, color = var)) +
#   geom_line(data = reform.censored, aes(x = x, y = y, color = var)) +
#   theme_bw()
# show(p)
# 
# pretest <- sum(raw.lytA$pneumo == '1')/nrow(raw.lytA)
# reform <- reform[reform$sens != 0 & reform$spec != 0, ]
# 
# p <- ggplot() + 
#   geom_point(data = reform, aes(x = lytA, y = sens / (1 - spec), color = 'LR Positive')) +
#   geom_point(data = reform, aes(x = lytA, y = (1 - sens) / spec, color = 'LR Negative')) +
#   geom_smooth(data = reform, aes(x = lytA, y = sens / (1 - spec), color = 'LR Positive')) +
#   geom_smooth(data = reform, aes(x = lytA, y = (1 - sens) / spec, color = 'LR Negative')) +
#   geom_vline(xintercept = 8) + 
#   ylim(0,5) + 
#   theme_bw()
# show(p)
# 
# p <- ggplot() + 
#   geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'Pretest Probability')) +
#   geom_point(data = reform, aes(x = lytA, y = (pretest * (sens / (1 - spec)) / (1 - pretest)) / ((pretest * (sens / (1 - spec)) / (1 - pretest)) + 1), color = 'Positive Probability')) +
#   geom_point(data = reform, aes(x = lytA, y = (pretest * ((1 - sens) / spec) / (1 - pretest)) / ((pretest * ((1 - sens) / spec) / (1 - pretest)) + 1), color = 'Negative Probability')) +
#   geom_smooth(data = reform, aes(x = lytA, y = (pretest * (sens / (1 - spec)) / (1 - pretest)) / ((pretest * (sens / (1 - spec)) / (1 - pretest)) + 1), color = 'Positive Probability')) +
#   geom_smooth(data = reform, aes(x = lytA, y = (pretest * ((1 - sens) / spec) / (1 - pretest)) / ((pretest * ((1 - sens) / spec) / (1 - pretest)) + 1), color = 'Negative Probability')) +
#   geom_vline(xintercept = 8) + 
#   ylim(0,1) + 
#   theme_bw()
# show(p)
# 
# #reform$sens <- rep(1.0, nrow(reform))
# p <- ggplot() + 
#   #geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'Pretest Probability')) +
#   geom_point(data = reform, aes(x = lytA, y = (pretest * (sens / (1 - spec)) / (1 - pretest)) / ((pretest * (sens / (1 - spec)) / (1 - pretest)) + 1) - (pretest * ((1 - sens) / spec) / (1 - pretest)) / ((pretest * ((1 - sens) / spec) / (1 - pretest)) + 1), color = 'Probability Difference')) +
#   geom_smooth(data = reform, aes(x = lytA, y = (pretest * (sens / (1 - spec)) / (1 - pretest)) / ((pretest * (sens / (1 - spec)) / (1 - pretest)) + 1) - (pretest * ((1 - sens) / spec) / (1 - pretest)) / ((pretest * ((1 - sens) / spec) / (1 - pretest)) + 1), color = 'Probability Difference')) +
#   geom_vline(xintercept = 8) + 
#   ylim(0,1) + 
#   theme_bw()
# show(p)