## Clear workspace and set working directory
rm(list = ls())
setwd('~/Google Drive/Documents/PostPostDoc/pneumococcal_colonization_analysis_redo/')

## Library Imports
libraries.call <- c("dplyr", 'tidyr', "ggplot2", "readxl", "reshape")
lapply(libraries.call, require, character.only = TRUE)

## Utility Functions
calc.pretest <- function(df, variable) {
  pretest <- sum(df[[variable]] == '1') / nrow(df)
}

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
  df.series <- seq(0, max(df$variable), length = 1000)
  df.stats <- data.frame()
  invisible(sapply(df.series, function(x) df.stats <<- rbind(df.stats, (calculate.test.stats(generate.matrix(df, cutoff = x, variable = "variable", response = which.variable))))))
  df.stats <- cbind(x.value = df.series, df.stats)
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
    geom_line(data = df, aes(x = x, y = y), size = 1.5) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("1 - specificity") +
    ylab("sensitivity") +
    theme_bw()
  
  ggsave(plot = p, filename = paste("ROC_", file.name, sep = ''), height = 5, width = 6)
  return(p)
}

make.logistic.plot <- function(df, file.name) {
  fit.df <- make.fit(df)
  new.df <- data.frame(variable = df$variable, bacter = df$bacter, prob.bacter = fit.df$bacter.fit, pneumo = df$pneumo, prob.pneumo = fit.df$pneumo.fit)
  
  p <- ggplot() + 
    geom_point(data = new.df, aes(x = variable, y = as.numeric(as.character(pneumo)), color = "Diagnosis Status")) +
    geom_line(data = new.df, aes(x = variable, y = as.numeric(as.character(prob.pneumo)), color = 'Logistic Fit'), size = 1.5) +
    xlab("substrate concentration (units)") +
    ylab("pneumococcal pneumonia diagnosis") +
    scale_color_discrete(name = "") +
    theme_bw()
  
  ggsave(plot = p, filename = paste("Logistic_", file.name, sep = ''), height = 5, width = 7)
  return(p)
}

make.sens.spec.plot <- function(df, file.name) {
  df.reformed <- melt(df, id = c("x.value"))
  colnames(df.reformed) <- c('x', 'var', 'y')
  
  reformed.sens <- left.censor(filter(df.reformed, var == 'sens'))
  reformed.sens <- cbind(reformed.sens, var = rep('sens', nrow(reformed.sens)))
  reformed.sens <- reformed.sens[rev(order(reformed.sens$y)), ]
  
  reformed.spec <- left.censor(filter(df.reformed, var == 'spec'))
  reformed.spec <- cbind(reformed.spec, var = rep('spec', nrow(reformed.spec)))
  reformed.spec <- reformed.spec[order(reformed.spec$y), ]
  reformed.censored <- rbind(reformed.sens, reformed.spec)
  
  p <- ggplot() + 
    geom_line(data = reformed.censored, aes(x = x, y = y, color = var), size = 1.5) +
    xlab("substrate concentration (units)") +
    ylab("") +
    scale_color_discrete(name = "", labels = c(sens = "sensitivity", spec = "specificity")) +
    theme_bw()
  
  ggsave(plot = p, filename = paste("SensSpec_", file.name, sep = ''), height = 5, width = 7)
  return(list(plot = p, df.reformed = df.reformed))
}

make.likelihoodratio.plot <- function(df, file.name) {
  df <- df[df$sens != 0 & df$spec != 0, ]
  
  p <- ggplot() + 
    geom_segment(aes(x = 0, y = 1, xend = max(df$x.value), yend = 1), linetype = 2) +
    geom_point(data = df, aes(x = x.value, y = sens / (1 - spec), color = 'LR positive'), size = 0.5) +
    geom_point(data = df, aes(x = x.value, y = (1 - sens) / spec, color = 'LR negative'), size = 0.5) +
    geom_smooth(data = df, aes(x = x.value, y = sens / (1 - spec), color = 'LR positive'), size = 0.5, se = F) +
    geom_smooth(data = df, aes(x = x.value, y = (1 - sens) / spec, color = 'LR negative'), size = 0.5, se = F) +
    ylim(0,5) + 
    xlab("substrate concentration (units)") +
    ylab("likelihood ratio") +
    scale_color_discrete(name = "") +
    theme_bw()
  
  ggsave(plot = p, filename = paste("LR_", file.name, sep = ''), height = 5, width = 7)
  return(list(plot = p))
}

make.probability.plots <- function(df.raw, df, file.name) {
  require(cowplot)
  pretest <- calc.pretest(df.raw, "pneumo")
  
  df <- df[df$sens != 0 & df$spec != 0, ]
  
  p1 <- ggplot() + 
    #geom_vline(xintercept = 2) +
    geom_segment(aes(x = 0, y = pretest, xend = max(df$x.value), yend = pretest, color = 'a'), linetype = 2) +
    geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'a')) +
    geom_point(data = df, aes(x = x.value, y = (pretest * (sens / (1 - spec)) / (1 - pretest)) / ((pretest * (sens / (1 - spec)) / (1 - pretest)) + 1), color = "b"), size = 0.5) +
    geom_point(data = df, aes(x = x.value, y = (pretest * ((1 - sens) / spec) / (1 - pretest)) / ((pretest * ((1 - sens) / spec) / (1 - pretest)) + 1), color = "c"), size = 0.5) +
    geom_smooth(data = df, aes(x = x.value, y = (pretest * (sens / (1 - spec)) / (1 - pretest)) / ((pretest * (sens / (1 - spec)) / (1 - pretest)) + 1), color = "b"), size = 0.5, se = F) +
    geom_smooth(data = df, aes(x = x.value, y = (pretest * ((1 - sens) / spec) / (1 - pretest)) / ((pretest * ((1 - sens) / spec) / (1 - pretest)) + 1), color = "c"), size = 0.5, se = F) +
    ylim(0,1) + 
    xlab("substrate concentration (units)") +
    ylab("probability of pneumococcal pneumonia") +
    scale_color_discrete(name = "", labels = c(a = 'pretest probability', b = 'posttest probability positive', c = 'posttest probability negative')) +
    theme_bw()
  
  p2 <- ggplot() + 
    #geom_vline(xintercept = 2) + 
    geom_point(data = df, aes(x = x.value, 
                              y = (pretest * (sens / (1 - spec)) / (1 - pretest)) / ((pretest * (sens / (1 - spec)) / (1 - pretest)) + 1) - (pretest * ((1 - sens) / spec) / (1 - pretest)) / ((pretest * ((1 - sens) / spec) / (1 - pretest)) + 1), 
                              color = 'a'), 
               size = 0.5) +
    geom_smooth(data = df, aes(x = x.value, 
                               y = (pretest * (sens / (1 - spec)) / (1 - pretest)) / ((pretest * (sens / (1 - spec)) / (1 - pretest)) + 1) - (pretest * ((1 - sens) / spec) / (1 - pretest)) / ((pretest * ((1 - sens) / spec) / (1 - pretest)) + 1), 
                               color = 'a'), 
                size = 0.5,
                se = F) +
    ylim(0,1) + 
    xlab("substrate concentration (units)") +
    ylab("posttest probability difference") +
    scale_color_discrete(name = "", labels = c(a = 'probability difference')) +
    theme_bw()
  
  p <- plot_grid(p1, p2, labels = c('A', 'B'))
  
  ggsave(plot = p, filename = paste("Prob_", file.name, sep = ''), height = 5, width = 15)
  return(list(plot = p))
}

lytA.filename <- 'lytA.pdf'
raw.lytA <- import.data(variable = 'lytA_NP')
make.logistic.plot(df = raw.lytA, file.name = lytA.filename)
roc.data.lytA <- make.test.stats.df(df = raw.lytA, which.variable = "pneumo")
make.roc.plot(df = roc.data.lytA, file.name = lytA.filename)
p.sens.spec <- make.sens.spec.plot(df = roc.data.lytA, file.name = lytA.filename)
make.likelihoodratio.plot(df = roc.data.lytA, file.name = lytA.filename)
make.probability.plots(df.raw = raw.lytA, df = roc.data.lytA, file.name = lytA.filename)

pct.filename <- 'pct.pdf'
raw.pct <- import.data(variable = 'PCT')
make.logistic.plot(df = raw.pct, file.name = pct.filename)
roc.data.pct <- make.test.stats.df(df = raw.pct, which.variable = "pneumo")
make.roc.plot(df = roc.data.pct, file.name = pct.filename)
p.sens.spec <- make.sens.spec.plot(roc.data.pct, file.name = pct.filename)
make.likelihoodratio.plot(roc.data.pct, pct.filename)
make.probability.plots(df.raw = raw.pct, df = roc.data.pct, file.name = pct.filename)


