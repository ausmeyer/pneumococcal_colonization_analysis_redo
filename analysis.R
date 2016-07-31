## Clear workspace and set working directory
rm(list = ls())
setwd('~/Google Drive/Documents/PostPostDoc/pneumococcal_colonization_analysis_redo/')

## Library Imports
libraries.call <- c("dplyr", 'tidyr', "ggplot2", "readxl", "reshape", 'cowplot')
lapply(libraries.call, require, character.only = TRUE)

## Utility Functions
calc.pretest <- function(df, variable) {
  x <- sum(df[[variable]] == '1') / nrow(df)
}

calc.prob.difference <- function(df, pretest) {
  pretest.odds <- pretest / (1 - pretest)
  lr.pos <- (df$sens / (1 - df$spec))
  lr.neg <- ((1 - df$sens) / df$spec)
  pos.odds <- pretest.odds * lr.pos
  neg.odds <- pretest.odds * lr.neg
  pos.prob <- pos.odds / (1 + pos.odds)
  neg.prob <- neg.odds / (1 + neg.odds)
  
  return(data.frame(lr.pos = lr.pos, lr.neg = lr.neg, pos.prob = pos.prob, neg.prob = neg.prob))
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
  sigma2.pos <- (1 / df$t.p) - (1 / (df$t.p + df$f.n)) + (1 / df$f.p) - (1 / (df$f.p + df$t.n))
  sigma2.neg <- (1 / df$f.n) - (1 / (df$t.p + df$f.n)) + (1 / df$t.n) - (1 / (df$f.p + df$t.n))
  
  return(data.frame(sens = sensitivity, spec = specificity, sigma2.pos = sigma2.pos, sigma2.neg = sigma2.neg))
}

make.test.stats.df <- function(df, which.variable) {
  #df.series <- seq(0, max(df$variable), length = 1000)
  df.series <- df$variable[order(df$variable)]
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
  bacter.fit <- glm(factor(bacter) ~ variable, data = df, family = "binomial")
  predicted.bacteremia <- predict(bacter.fit, type="response", se = TRUE)
  
  pneumo.fit <- glm(factor(pneumo) ~ variable, data = df, family = "binomial")
  predicted.pneumo <- predict(pneumo.fit, type="response", se = TRUE)
  
  a <<- pneumo.fit
  
  return(list(predicted = data.frame(bacter = predicted.bacteremia, pneumo = predicted.pneumo), bacter.fit = bacter.fit, pneumo.fit = pneumo.fit))
}

## Plotting Functions
make.roc.plot <- function(df, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df, pretest)
  #df <- df[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  prob.diff <- calc.prob.difference(df, pretest)
  
  df.censored <- data.frame(x = 1 - df$spec, y = df$sens)
  df.censored <- left.censor(df.censored)
  df.censored <- df.censored[order(df.censored$y), ]
  
  p1 <- ggplot() + 
    geom_line(data = df.censored, aes(x = x, y = y, color = 'line'), size = 1) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("false positive rate (1 - sp)") +
    ylab("true positive rate (sn)") +
    xlim(0, 1) +
    ylim(0, 1) +
    scale_color_discrete(name = "") +
    theme_bw() +
    theme(legend.position = "none")
  
  p2 <- ggplot() + 
    geom_line(data = df, 
              aes(x = 1 - spec, y = sens + spec - 1, color = 'Youden index'), size = 1) +
    geom_segment(data = df, 
                 aes(x = min((1 - spec)[(sens + spec - 1) == max(sens + spec - 1)]), 
                     y = 0, 
                     xend = min((1 - spec)[(sens + spec - 1) == max(sens + spec - 1)]),
                     yend =  max(sens + spec - 1), 
                     color = 'optimal cutoff'),
                 size = 0.75) +
    geom_point(data = df,
               aes(x = min((1 - spec)[(sens + spec - 1) == max(sens + spec - 1)]), 
                   y = max(sens + spec - 1), 
                   color = 'optimal cutoff')) +
    geom_hline(yintercept = 0) +
    xlab("false positive rate (1 - sp)") +
    ylab("Youden index (sn + sp - 1)") +
    xlim(0, 1) +
    ylim(-1, 1) +
    theme_bw() +
    theme(legend.position = c(0.5, 0.27), legend.title=element_blank())
  
  print(median((df$x.value)[(df$sens + df$spec - 1) == max(df$sens + df$spec - 1)]))
  
  p3 <- ggplot() + 
    geom_line(data = df, aes(x = x.value, y = df$sens + df$spec - 1, color = 'Youden index'), size = 1) +
    geom_segment(data = df, 
                 aes(x = min((x.value)[(sens + spec - 1) == max(sens + spec - 1)]), 
                     y = 0, 
                     xend = min((x.value)[(sens + spec - 1) == max(sens + spec - 1)]),
                     yend =  max(sens + spec - 1), 
                     color = 'optimal cutoff'),
                 size = 0.75) +
    geom_point(data = df, 
               aes(x = min((x.value)[(sens + spec - 1) == max(sens + spec - 1)]),
                   y =  max(sens + spec - 1),
                   color = 'optimal cutoff')
    ) +
    geom_hline(yintercept = 0) +
    xlab(xaxis) +
    ylab("Youden index (sn + sp - 1)") +
    ylim(-1, 1) +
    theme_bw() +
    theme(legend.position = c(0.5, 0.27), legend.title=element_blank())
  
  return(list(plot1 = p1, plot2 = p2, plot3 = p3))
}

make.logistic.plot <- function(df, df.roc, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df.roc, pretest)
  df.roc <- df.roc[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  prob.diff <- calc.prob.difference(df.roc, pretest)
  
  df <- df[df$variable <= max(df.roc$x.value), ]
  
  fit.df <- (make.fit(df))[["predicted"]]
  new.df <- data.frame(variable = df$variable, bacter = df$bacter, prob.bacter = fit.df$bacter.fit, pneumo = df$pneumo, prob.pneumo = fit.df$pneumo.fit)
  
  p <- ggplot() + 
    geom_point(data = new.df, aes(x = variable, y = as.numeric(as.character(pneumo)), color = "diagnostic status")) +
    geom_line(data = new.df, aes(x = variable, y = as.numeric(as.character(prob.pneumo)), color = 'logistic fit'), size = 1) +
    xlab(xaxis) +
    ylab("pneumococcal diagnosis") +
    scale_color_discrete(name = "") +
    theme_bw()
  
  return(list(plot = p))
}

make.sens.spec.plot <- function(df, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df, pretest)
  #df <- df[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  prob.diff <- calc.prob.difference(df, pretest)
  
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
    geom_line(data = reformed.censored, aes(x = x, y = y, color = var), size = 1) +
    xlab(xaxis) +
    ylab("") +
    scale_color_discrete(name = "", labels = c(sens = "sensitivity", spec = "specificity")) +
    theme_bw()
  
  return(list(plot = p, df.reformed = df.reformed))
}

make.likelihoodratio.plot <- function(df, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df, pretest)
  df <- df[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  prob.diff <- calc.prob.difference(df, pretest)
  
  alpha <- 0.05
  max.ribbon <- 15
  
  lower.pos <- prob.diff$lr.pos * exp(-qnorm(1 - (alpha / 2)) * sqrt(df$sigma2.pos))
  upper.pos <- prob.diff$lr.pos * exp(qnorm(1 - (alpha / 2)) * sqrt(df$sigma2.pos)) 
  
  lower.neg <- prob.diff$lr.neg * exp(-qnorm(1 - (alpha / 2)) * sqrt(df$sigma2.neg))
  upper.neg <- prob.diff$lr.neg * exp(qnorm(1 - (alpha / 2)) * sqrt(df$sigma2.neg))
  
  upper.pos[upper.pos > max.ribbon] <- max.ribbon
  upper.neg[upper.neg > max.ribbon] <- max.ribbon
  
  df <- cbind(df, lr.pos = prob.diff$lr.pos, lr.neg = prob.diff$lr.neg, lower.pos, upper.pos, lower.neg, lower.pos)
  
  p1 <- ggplot() + 
    geom_ribbon(data = df, aes(x = x.value, ymin = lower.pos, ymax=upper.pos), alpha = 0.2) +
    geom_ribbon(data = df, aes(x = x.value, ymin = lower.neg, ymax=upper.neg), alpha = 0.2) +
    geom_segment(aes(x = 0, y = 1, xend = max(df$x.value), yend = 1), linetype = 2) +
    geom_point(data = df, aes(x = x.value, y = lr.pos, color = 'LR positive'), size = 0.5) +
    geom_point(data = df, aes(x = x.value, y = lr.neg, color = 'LR negative'), size = 0.5) +
    geom_smooth(data = df, aes(x = x.value, y = lr.pos, color = 'LR positive'), size = 0.5, se = F) +
    geom_smooth(data = df, aes(x = x.value, y = lr.neg, color = 'LR negative'), size = 0.5, se = F) +
    ylim(0, max.ribbon) + 
    xlab(xaxis) +
    ylab("likelihood ratio") +
    scale_color_discrete(name = "") +
    theme_bw() +
    theme(legend.position = c(0.35, 0.85), legend.title=element_blank())
  
  p2 <- ggplot() + 
    geom_segment(aes(x = 0, y = 1, xend = max(df$x.value), yend = 1), linetype = 2) +
    geom_point(data = df, aes(x = x.value, y = lr.pos, color = 'LR positive'), size = 0.5) +
    geom_point(data = df, aes(x = x.value, y = lr.neg, color = 'LR negative'), size = 0.5) +
    geom_smooth(data = df, aes(x = x.value, y = lr.pos, color = 'LR positive'), size = 0.5, se = F) +
    geom_smooth(data = df, aes(x = x.value, y = lr.neg, color = 'LR negative'), size = 0.5, se = F) +
    ylim(0, 10) +
    xlab(xaxis) +
    ylab("likelihood ratio") +
    scale_color_discrete(name = "") +
    theme_bw() +
    theme(legend.position = c(0.35, 0.85), legend.title=element_blank())
  
  return(list(plot1 = p1, plot2 = p2))
}

make.continuous.likelihoodratio.plot <- function(df, df.raw, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df, pretest)
  df <- df[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  prob.diff <- calc.prob.difference(df, pretest)
  
  df.raw <- df.raw[df.raw$variable <= max(df$x.value), ]
  fit <- make.fit(df.raw)
  
  generate.continuous.lr <- function(tmp.df, tmp.fit, tmp.pretest) {
    x.1 <- (log(tmp.pretest / (1 - tmp.pretest)) - (tmp.fit[['pneumo.fit']])[["coefficients"]][[1]]) / (tmp.fit[['pneumo.fit']])[["coefficients"]][[2]]
    tmp.lr.continuous <- exp((tmp.fit[['pneumo.fit']])[["coefficients"]][[2]] * (tmp.df$variable - x.1))
    tmp.df <- cbind(tmp.df, tmp = tmp.lr.continuous)
    colnames(tmp.df)[colnames(tmp.df) == 'tmp'] <- paste('p_', tmp.pretest, sep = '')
    return(tmp.df)
  }
  
  df.raw.tmp <- df.raw
  df.raw.tmp <- generate.continuous.lr(df.raw.tmp, fit, pretest)
  df.raw.tmp <- df.raw.tmp[, c(-1, -3, -4)]
  colnames(df.raw.tmp) <- c('x.value', 'y.value')
  
  p1 <- ggplot() + 
    geom_segment(data = df.raw.tmp, aes(x = 0, y = 1, xend = max(x.value), yend = 1), linetype = 2) +
    geom_point(data = df.raw.tmp, aes(x = x.value, y = y.value), size = 0.5) +
    ylim(0, 10) +
    xlab(xaxis) +
    ylab("likelihood ratio") +
    scale_color_discrete(name = "") +
    theme_bw() +
    theme(legend.position = "none")
  
  df.raw.tmp <- df.raw
  lapply(seq(0.01, 0.99, length = 10), function(x) df.raw.tmp <<- generate.continuous.lr(df.raw.tmp, fit, x))
  
  df.raw.tmp <- df.raw.tmp[, c(-1, -3, -4)]
  colnames(df.raw.tmp)[1] <- 'x.value'
  df.raw.tmp <- melt(df.raw.tmp, id = c('x.value'))

  p2 <- ggplot() + 
    geom_segment(data = df.raw.tmp, aes(x = 0, y = 1, xend = max(x.value), yend = 1), linetype = 2) +
    geom_point(data = df.raw.tmp, aes(x = x.value, y = value, color = variable), size = 0.5) +
    ylim(0, 10) +
    xlab(xaxis) +
    ylab("likelihood ratio") +
    scale_color_discrete(name = "") +
    theme_bw() +
    theme(legend.position = "none")
  
  return(list(plot1 = p1, plot2 = p2))
}

make.probability.plots <- function(df.raw, df, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df, pretest)
  df <- df[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  prob.diff <- calc.prob.difference(df, pretest)
  
  p1 <- ggplot() + 
    #geom_vline(xintercept = 2) +
    geom_segment(aes(x = 0, y = pretest, xend = max(df$x.value), yend = pretest, color = 'a'), linetype = 2) +
    geom_point(data = df, aes(x = x.value, y = prob.diff$pos.prob, color = "b"), size = 0.5) +
    geom_point(data = df, aes(x = x.value, y = prob.diff$neg.prob, color = "c"), size = 0.5) +
    geom_smooth(data = df, aes(x = x.value, y = prob.diff$pos.prob, color = "b"), size = 0.5, se = F) +
    geom_smooth(data = df, aes(x = x.value, y = prob.diff$neg.prob, color = "c"), size = 0.5, se = F) +
    geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'a')) +
    geom_point(data = data.frame(x = max(df$x.value), y = pretest), aes(x = x, y = y, color = 'a')) +
    ylim(0, 1) + 
    xlab(xaxis) +
    ylab("probability pneumococcal") +
    scale_color_discrete(name = "", labels = c(a = 'pretest probability', b = 'posttest prob. positive', c = 'posttest prob. negative')) +
    theme_bw()
  
  print(min(df$x.value[(prob.diff$pos.prob - prob.diff$neg.prob) == max(prob.diff$pos.prob - prob.diff$neg.prob)]))
  print(max(prob.diff$pos.prob - prob.diff$neg.prob))
  
  p2 <- ggplot() + 
    geom_point(data = df, aes(x = x.value, 
                              y = prob.diff$pos.prob - prob.diff$neg.prob, 
                              color = 'b'), 
               size = 0.5) +
    geom_smooth(data = df, aes(x = x.value, 
                               y = prob.diff$pos.prob - prob.diff$neg.prob, 
                               color = 'b'), 
                size = 0.5,
                se = F) +
    geom_segment(data = df, 
                 aes(x = min(x.value[(prob.diff$pos.prob - prob.diff$neg.prob) == max(prob.diff$pos.prob - prob.diff$neg.prob)]), 
                     y = 0, 
                     xend = min(x.value[(prob.diff$pos.prob - prob.diff$neg.prob) == max(prob.diff$pos.prob - prob.diff$neg.prob)]),
                     yend =  max((prob.diff$pos.prob - prob.diff$neg.prob)), 
                     color = 'a'),
                 size = 0.75) +
    geom_point(data = df, 
               aes(x = min(x.value[(prob.diff$pos.prob - prob.diff$neg.prob) == max(prob.diff$pos.prob - prob.diff$neg.prob)]),
                   y =  max((prob.diff$pos.prob - prob.diff$neg.prob)),
                   color = 'a')) +
    ylim(0, 1) + 
    xlab(xaxis) +
    ylab("posttest probability diff") +
    scale_color_discrete(name = "", labels = c(a = 'optimal cutoff', b = 'probability difference')) +
    theme_bw()
  
  return(list(plot1 = p1, plot2 = p2, df = df, prob.diff = prob.diff$pos.prob - prob.diff$neg.prob))
}

make.combined.probability.plots <- function(df.raw, df, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df, pretest)
  df <- df[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  prob.diff <- calc.prob.difference(df, pretest)

  combined.probability <- (pretest / (1 - pretest)) * prob.diff$lr.pos * prob.diff$lr.neg / ((pretest / (1 - pretest)) * prob.diff$lr.pos * prob.diff$lr.neg + 1)
  
  df <- cbind(df, combined.probability)
    
  p1 <- ggplot() + 
    geom_segment(aes(x = 0, y = pretest, xend = max(df$x.value), yend = pretest, color = 'a'), linetype = 2) +
    geom_point(data = df, aes(x = x.value, y = combined.probability, color = "b"), size = 0.5) +
    #geom_smooth(data = df, aes(x = x.value, y = combined.probability, color = "b"), size = 0.5, span = 1) +
    geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'a')) +
    geom_point(data = data.frame(x = max(df$x.value), y = pretest), aes(x = x, y = y, color = 'a')) +
    ylim(0, 1) + 
    xlab(xaxis) +
    ylab("probability pneumococcal") +
    scale_color_discrete(name = "", labels = c(a = 'pretest probability', b = 'combined posttest probability')) +
    theme_bw() +
    theme(legend.position = c(0.35, 0.85), legend.title=element_blank())
  
  df.raw <- df.raw[df.raw$variable <= max(df$x.value), ]
  fit <- make.fit(df.raw)
  
  generate.continuous.lr <- function(tmp.df, tmp.fit, tmp.pretest) {
    x.1 <- (log(tmp.pretest / (1 - tmp.pretest)) - (tmp.fit[['pneumo.fit']])[["coefficients"]][[1]]) / (tmp.fit[['pneumo.fit']])[["coefficients"]][[2]]
    tmp.lr.continuous <- exp((tmp.fit[['pneumo.fit']])[["coefficients"]][[2]] * (tmp.df$variable - x.1))
    tmp.posttest <- (tmp.pretest / (1 - tmp.pretest)) * tmp.lr.continuous / (1 + (tmp.pretest / (1 - tmp.pretest)) * tmp.lr.continuous)
    tmp.df <- cbind(tmp.df, tmp = tmp.posttest)
    colnames(tmp.df)[colnames(tmp.df) == 'tmp'] <- paste('p_', tmp.pretest, sep = '')
    return(tmp.df)
  }
  
  df.raw.tmp <- df.raw
  df.raw.tmp <- generate.continuous.lr(df.raw.tmp, fit, pretest)
  df.raw.tmp <- df.raw.tmp[, c(-1, -3, -4)]
  colnames(df.raw.tmp) <- c('x.value', 'y.value')
  
  p2 <- ggplot() + 
    geom_segment(aes(x = 0, y = pretest, xend = max(df.raw.tmp$x.value), yend = pretest, color = 'a'), linetype = 2) +
    geom_point(data = df.raw.tmp, aes(x = x.value, y = y.value, color = "b"), size = 0.5) +
    geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'a')) +
    geom_point(data = data.frame(x = max(df.raw.tmp$x.value), y = pretest), aes(x = x, y = y, color = 'a')) +
    ylim(0, 1) + 
    xlab(xaxis) +
    ylab("probability pneumococcal") +
    scale_color_discrete(name = "", labels = c(a = 'pretest probability', b = 'continuous posttest probability')) +
    theme_bw() +
    theme(legend.position = c(0.35, 0.85), legend.title=element_blank())
  
  return(list(plot1 = p1, plot2 = p2))
}

raw.pneumo <- import.data(variable = "Pneumococcal_diagnosis")
pretest.probability <- calc.pretest(raw.pneumo, 'pneumo')

## lytA setup
lytA.filename <- 'lytA.pdf'
lytA.xaxis <- 'lytA density (log10 copies/mL)'
raw.lytA <- import.data(variable = 'lytA_NP')

## Plot lytA
roc.data.lytA <- make.test.stats.df(df = raw.lytA, which.variable = "pneumo")

p.log.lytA <- make.logistic.plot(df = raw.lytA, df.roc = roc.data.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)
p.roc.lytA <- make.roc.plot(df = roc.data.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)

p.sens.spec.lytA <- make.sens.spec.plot(df = roc.data.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)
p.lr.lytA <- make.likelihoodratio.plot(df = roc.data.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)
p.continuous.lr.lytA <- make.continuous.likelihoodratio.plot(df = roc.data.lytA, df.raw = raw.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)
p.prob.lytA <- make.probability.plots(df.raw = raw.lytA, df = roc.data.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)
p.combined.prob.lytA <- make.combined.probability.plots(df.raw = raw.lytA, df = roc.data.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)

print(mean(p.prob.lytA$prob.diff[(p.prob.lytA[["df"]])[["x.value"]] < 6]))

## pct setup
pct.filename <- 'pct.pdf'
pct.xaxis <- 'procalcitonin concentration (ng/mL)'
raw.pct <- import.data(variable = 'PCT')

## Plot pct
roc.data.pct <- make.test.stats.df(df = raw.pct, which.variable = "pneumo")

p.log.pct <- make.logistic.plot(df = raw.pct, df.roc = roc.data.pct, file.name = pct.filename, xaxis = pct.xaxis, pretest = pretest.probability)
p.roc.pct <- make.roc.plot(df = roc.data.pct, file.name = pct.filename, xaxis = pct.xaxis, pretest = pretest.probability)

p.sens.spec.pct <- make.sens.spec.plot(roc.data.pct, file.name = pct.filename, xaxis = pct.xaxis, pretest = pretest.probability)
p.lr.pct <- make.likelihoodratio.plot(roc.data.pct, pct.filename, xaxis = pct.xaxis, pretest = pretest.probability)
p.continuous.lr.pct <- make.continuous.likelihoodratio.plot(df = roc.data.pct, df.raw = raw.pct, file.name = pct.filename, xaxis = pct.xaxis, pretest = pretest.probability)
p.prob.pct <- make.probability.plots(df.raw = raw.pct, df = roc.data.pct, file.name = pct.filename, xaxis = pct.xaxis, pretest = pretest.probability)
p.combined.prob.pct <- make.combined.probability.plots(df.raw = raw.pct, df = roc.data.pct, file.name = pct.filename, xaxis = pct.xaxis, pretest = pretest.probability)

print(mean(p.prob.pct$prob.diff[2 < (p.prob.pct[["df"]])[["x.value"]] & (p.prob.pct[["df"]])[["x.value"]] < 40]))

## crp setup
crp.filename <- 'crp.pdf'
crp.xaxis <- 'c-reactive protein concentration (mg/L)'
raw.crp <- import.data(variable = 'CRP')

## Plot crp
roc.data.crp <- make.test.stats.df(df = raw.crp, which.variable = "pneumo")

p.log.crp <- make.logistic.plot(df = raw.crp, df.roc = roc.data.crp, file.name = crp.filename, xaxis = crp.xaxis, pretest = pretest.probability)
p.roc.crp <- make.roc.plot(df = roc.data.crp, file.name = crp.filename, xaxis = crp.xaxis, pretest = pretest.probability)

p.sens.spec.crp <- make.sens.spec.plot(roc.data.crp, file.name = crp.filename, xaxis = crp.xaxis, pretest = pretest.probability)
p.lr.crp <- make.likelihoodratio.plot(roc.data.crp, crp.filename, xaxis = crp.xaxis, pretest = pretest.probability)
p.continuous.lr.crp <- make.continuous.likelihoodratio.plot(df = roc.data.crp, df.raw = raw.crp, file.name = crp.filename, xaxis = crp.xaxis, pretest = pretest.probability)
p.prob.crp <- make.probability.plots(df.raw = raw.crp, df = roc.data.crp, file.name = crp.filename, xaxis = crp.xaxis, pretest = pretest.probability)
p.combined.prob.crp <- make.combined.probability.plots(df.raw = raw.crp, df = roc.data.crp, file.name = crp.filename, xaxis = crp.xaxis, pretest = pretest.probability)

print(mean(p.prob.crp$prob.diff[100 < (p.prob.crp[["df"]])[["x.value"]]]))

p.log <- plot_grid(p.log.crp$plot, p.log.pct$plot, p.log.lytA$plot, labels = c('A', 'B', 'C'), ncol = 3)
ggsave(plot = p.log, filename = "logistic_fits.pdf", height = 3, width = 15)

p.roc <- plot_grid(p.roc.crp$plot1, p.roc.pct$plot1, p.roc.lytA$plot1, 
                   p.roc.crp$plot2, p.roc.pct$plot2, p.roc.lytA$plot2, 
                   p.roc.crp$plot3, p.roc.pct$plot3, p.roc.lytA$plot3, 
                   labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'), 
                   ncol = 3,
                   scale = 0.95)
ggsave(plot = p.roc, filename = 'rocs.pdf', height = 9.5, width = 10.5)

p.sens.spec <- plot_grid(p.sens.spec.crp$plot, p.sens.spec.pct$plot, p.sens.spec.lytA$plot, labels = c('A', 'B', 'C'), ncol = 3)
ggsave(plot = p.sens.spec, filename = "sensitivities_specificities.pdf", height = 3, width = 15)

p.lr <- plot_grid(p.lr.crp$plot1, p.lr.pct$plot1, p.lr.lytA$plot1, 
                  p.lr.crp$plot2,  p.lr.pct$plot2,p.lr.lytA$plot2, 
                  labels = c('A', 'B', 'C', 'D', 'E', 'F'), ncol = 3)
ggsave(plot = p.lr, filename = "likelihood_ratios.pdf", height = 7.5, width = 12)

p.continuous.lr <- plot_grid(p.continuous.lr.crp$plot1, p.continuous.lr.pct$plot1, p.continuous.lr.lytA$plot1, 
                             p.continuous.lr.crp$plot2, p.continuous.lr.pct$plot2, p.continuous.lr.lytA$plot2,
                             labels = c('A', 'B', 'C', 'D', 'E', 'F'), ncol = 3)
ggsave(plot = p.continuous.lr, filename = "continuous_likelihood_ratios.pdf", height = 7.5, width = 12)

p.prob <- plot_grid(p.prob.crp$plot1, p.prob.pct$plot1, p.prob.lytA$plot1, 
                    p.prob.crp$plot2, p.prob.pct$plot2, p.prob.lytA$plot2, 
                    labels = c('A', 'B', 'C', 'D', 'E', 'F'), ncol = 3)
ggsave(plot = p.prob, filename = "probabilities.pdf", height = 6, width = 16)

p.combined.prob <- plot_grid(p.combined.prob.crp$plot1, p.combined.prob.pct$plot1, p.combined.prob.lytA$plot1, 
                             p.combined.prob.crp$plot2, p.combined.prob.pct$plot2, p.combined.prob.lytA$plot2,
                             labels = c('A', 'B', 'C', 'D', 'E', 'F'), ncol = 3)
ggsave(plot = p.combined.prob, filename = "combined_probabilities.pdf", height = 7.5, width = 12)
