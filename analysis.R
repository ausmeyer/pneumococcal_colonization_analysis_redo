## Clear workspace and set working directory
rm(list = ls())
setwd('~/Google Drive/Documents/PostPostDoc/pneumococcal_diagnosis_project/pneumococcal_colonization_analysis_redo/')

## Library Imports
libraries.call <- c("dplyr", 'tidyr', "ggplot2", "readr", "reshape", 'cowplot')
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
  
  smoothed.pos.prob <- as.data.frame(predict(loess(data = data.frame(x = df$x.value, y = pos.prob), y ~ x, span = 0.5), data.frame(x = df$x.value)))[,1]
  smoothed.neg.prob <- as.data.frame(predict(loess(data = data.frame(x = df$x.value, y = neg.prob), y ~ x, span = 0.5), data.frame(x = df$x.value)))[,1]
  
  return(data.frame(lr.pos = lr.pos, lr.neg = lr.neg, pos.prob = pos.prob, neg.prob = neg.prob, smoothed.pos.prob = smoothed.pos.prob, smoothed.neg.prob = smoothed.neg.prob))
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

  return(data.frame(sens = sensitivity, spec = specificity, sigma2.pos = sigma2.pos, sigma2.neg = sigma2.neg, t.p = df$t.p, t.n = df$t.n, f.p = df$f.p, f.n = df$f.n))
}

calculate.empiric.lr <- function(df) {
  return(df$t.p * df$t.n / (df$t.p^2 + 2*df$f.n*df$t.p + df$f.n^2) / df$f.p / df$f.n)
}

make.test.stats.df <- function(df, which.variable) {
  df.series <- df$variable[order(df$variable)]
  df.stats <- data.frame()
  invisible(sapply(df.series, function(x) df.stats <<- rbind(df.stats, (calculate.test.stats(generate.matrix(df, cutoff = x, variable = "variable", response = which.variable))))))
  df.stats <- cbind(x.value = df.series, df.stats)
  
  #Adjust minimum counts here
  df.stats <- df.stats[df.stats$t.p >= 0 & df.stats$t.n >= 0 & df.stats$f.p >= 0 & df.stats$f.n >= 0, ]
  
  return(df.stats)
}

make.empiric.lr.df <- function(df, which.variable) {
  df.series <- df$variable[order(df$variable)]
  df.stats <- data.frame()
  invisible(sapply(df.series, function(x) df.stats <<- rbind(df.stats, (calculate.empiric.lr(generate.matrix(df, cutoff = x, variable = "variable", response = which.variable))))))
  df.stats <- cbind(x.value = df.series, df.stats)
  colnames(df.stats)[2] <- 'y.value'
  
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

##############################
# End of accessory functions #
##############################

import.data <- function(variable, variable2) {
  df <- read_csv('data.csv', 
                 col_types = cols(
                   study_ID = col_character(),
                   Gender = col_character(),
                   age = col_integer(),
                   Bartlett_Score = col_integer(),
                   CD4_count = col_integer(),
                   Bactrim = col_integer(),
                   HAART = col_integer(),
                   CURB65 = col_integer(),
                   MR_proANP = col_double(),
                   MR_proADM = col_double(),
                   MR_proADM = col_double(),
                   PCT = col_double(),
                   CRP = col_double(),
                   lytA_NP = col_double(),
                   bacteremia = col_double(),
                   Pneumococcal_diagnosis = col_integer(),
                   Pneumococcal_diagnosis_expanded = col_integer()
                 ))
  
  if(missing(variable2)) {
    df <- df[order(df[[variable]]), ]
    
    df <- data.frame(age = df$age, variable = df[[variable]], bacter = df$bacteremia, pneumo = df$Pneumococcal_diagnosis)
    df <- na.omit(df)
    
    df <- df[(nrow(df) * 0.025):(nrow(df) * 0.975), ]
  }
  else {
    tmp.df <- data.frame(age = df$age, variable1 = df[[variable]], variable2 = df[[variable2]], bacter = df$bacteremia, pneumo = df$Pneumococcal_diagnosis)
    tmp.df <- na.omit(tmp.df)
    df <- data.frame(age = tmp.df$age, variable = tmp.df$variable1/tmp.df$variable2, bacter = tmp.df$bacter, pneumo = tmp.df$pneumo)
  }
  return(df)
}

make.fit <- function(df) {
  bacter.fit <- glm(factor(bacter) ~ variable, data = df, family = "binomial")
  predicted.bacteremia <- predict(bacter.fit, type="response", se = TRUE)
  
  pneumo.fit <- glm(factor(pneumo) ~ variable, data = df, family = "binomial")
  predicted.pneumo <- predict(pneumo.fit, type="response", se = TRUE)
  
  return(list(predicted = data.frame(bacter = predicted.bacteremia, pneumo = predicted.pneumo), bacter.fit = bacter.fit, pneumo.fit = pneumo.fit))
}

## Plotting Functions
make.roc.plot <- function(df.roc, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df.roc, pretest)

  df.censored <- data.frame(x = 1 - df.roc$spec, y = df.roc$sens)
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
    geom_line(data = df.roc, 
              aes(x = 1 - spec, y = sens + spec - 1, color = 'Youden index'), size = 1) +
    geom_segment(data = df.roc, 
                 aes(x = min((1 - spec)[(sens + spec - 1) == max(sens + spec - 1)]), 
                     y = 0, 
                     xend = min((1 - spec)[(sens + spec - 1) == max(sens + spec - 1)]),
                     yend =  max(sens + spec - 1), 
                     color = 'optimal cutoff'),
                 size = 0.75) +
    geom_point(data = df.roc,
               aes(x = min((1 - spec)[(sens + spec - 1) == max(sens + spec - 1)]), 
                   y = max(sens + spec - 1), 
                   color = 'optimal cutoff')) +
    geom_hline(yintercept = 0) +
    xlab("false positive rate (1 - sp)") +
    ylab("Youden index (sn + sp - 1)") +
    xlim(0, 1) +
    ylim(-1, 1) +
    theme_bw() +
    theme(legend.position = c(0.5, 0.27), 
          legend.title=element_blank(),
          legend.background = element_rect(fill = alpha('white', 0.0)))
  
  print(paste("Optimal false positive rate: ", round((1 - median((df.roc$spec)[(df.roc$sens + df.roc$spec - 1) == max(df.roc$sens + df.roc$spec - 1)])), digits = 2)), sep = '')
  print(paste("Optimal concentration: ", median((df.roc$x.value)[(df.roc$sens + df.roc$spec - 1) == max(df.roc$sens + df.roc$spec - 1)])), sep = '')
  
  p3 <- ggplot() + 
    geom_line(data = df.roc, aes(x = x.value, y = df.roc$sens + df.roc$spec - 1, color = 'Youden index'), size = 1) +
    geom_segment(data = df.roc, 
                 aes(x = min((x.value)[(sens + spec - 1) == max(sens + spec - 1)]), 
                     y = 0, 
                     xend = min((x.value)[(sens + spec - 1) == max(sens + spec - 1)]),
                     yend =  max(sens + spec - 1), 
                     color = 'optimal cutoff'),
                 size = 0.75) +
    geom_point(data = df.roc, 
               aes(x = min((x.value)[(sens + spec - 1) == max(sens + spec - 1)]),
                   y =  max(sens + spec - 1),
                   color = 'optimal cutoff')
    ) +
    geom_hline(yintercept = 0) +
    xlab(xaxis) +
    ylab("Youden index (sn + sp - 1)") +
    ylim(-1, 1) +
    theme_bw() +
    theme(legend.position = c(0.5, 0.27), 
          legend.title=element_blank(),
          legend.background = element_rect(fill = alpha('white', 0.0)))
  
  return(list(plot1 = p1, plot2 = p2, plot3 = p3))
}

make.logistic.plot <- function(df, df.roc, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df.roc, pretest)
  prob.diff.tmp <- prob.diff[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  df.roc <- df.roc[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  prob.diff <- prob.diff.tmp
  
  df <- df[df$variable <= max(df.roc$x.value), ]
  
  fit.df <- (make.fit(df))[["predicted"]]
  new.df <- data.frame(variable = df$variable, bacter = df$bacter, prob.bacter = fit.df$bacter.fit, pneumo = df$pneumo, prob.pneumo = fit.df$pneumo.fit)
  
  #print(summary(make.fit(df)[["pneumo.fit"]]))
  #print(nrow(df))
  
  p <- ggplot() + 
    geom_point(data = new.df, aes(x = variable, y = as.numeric(as.character(pneumo)), color = "diagnostic status")) +
    geom_line(data = new.df, aes(x = variable, y = as.numeric(as.character(prob.pneumo)), color = 'logistic fit'), size = 1) +
    xlab(xaxis) +
    ylab("pneumococcal diagnosis") +
    scale_color_discrete(name = "") +
    theme_bw()
  
  return(list(plot = p))
}

make.sens.spec.plot <- function(df.roc, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df.roc, pretest)
  
  df.reformed <- melt(df.roc, id = c("x.value"))
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

make.likelihoodratio.plot <- function(df.roc, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df.roc, pretest)
  prob.diff.tmp <- prob.diff[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  df.roc <- df.roc[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  prob.diff <- prob.diff.tmp
  
  alpha <- 0.05
  max.ribbon <- 15
  
  lower.pos <- prob.diff$lr.pos * exp(-qnorm(1 - (alpha / 2)) * sqrt(df.roc$sigma2.pos))
  upper.pos <- prob.diff$lr.pos * exp(qnorm(1 - (alpha / 2)) * sqrt(df.roc$sigma2.pos)) 
  
  lower.neg <- prob.diff$lr.neg * exp(-qnorm(1 - (alpha / 2)) * sqrt(df.roc$sigma2.neg))
  upper.neg <- prob.diff$lr.neg * exp(qnorm(1 - (alpha / 2)) * sqrt(df.roc$sigma2.neg))
  
  upper.pos[upper.pos > max.ribbon] <- max.ribbon
  upper.neg[upper.neg > max.ribbon] <- max.ribbon
  
  df.roc <- cbind(df.roc, lr.pos = prob.diff$lr.pos, lr.neg = prob.diff$lr.neg, lower.pos, upper.pos, lower.neg, lower.pos)
  
  p1 <- ggplot() + 
    geom_ribbon(data = df.roc, aes(x = x.value, ymin = lower.pos, ymax=upper.pos), alpha = 0.2) +
    geom_ribbon(data = df.roc, aes(x = x.value, ymin = lower.neg, ymax=upper.neg), alpha = 0.2) +
    geom_segment(aes(x = 0, y = 1, xend = max(df.roc$x.value), yend = 1), linetype = 2) +
    geom_point(data = df.roc, aes(x = x.value, y = lr.pos, color = 'LR positive'), size = 0.5) +
    geom_point(data = df.roc, aes(x = x.value, y = lr.neg, color = 'LR negative'), size = 0.5) +
    geom_smooth(data = df.roc, aes(x = x.value, y = lr.pos, color = 'LR positive'), span = 0.5, size = 0.5, se = F) +
    geom_smooth(data = df.roc, aes(x = x.value, y = lr.neg, color = 'LR negative'), span = 0.5, size = 0.5, se = F) +
    ylim(0, max.ribbon) + 
    xlab(xaxis) +
    ylab("likelihood ratio") +
    scale_color_discrete(name = "") +
    theme_bw() +
    theme(legend.position = c(0.35, 0.85), 
          legend.title=element_blank(),
          legend.background = element_rect(fill = alpha('white', 0.0)))
  
  p2 <- ggplot() + 
    geom_segment(aes(x = 0, y = 1, xend = max(df.roc$x.value), yend = 1), linetype = 2) +
    geom_point(data = df.roc, aes(x = x.value, y = lr.pos, color = 'LR positive'), size = 0.5) +
    geom_point(data = df.roc, aes(x = x.value, y = lr.neg, color = 'LR negative'), size = 0.5) +
    geom_smooth(data = df.roc, aes(x = x.value, y = lr.pos, color = 'LR positive'), span = 0.5, size = 0.5, se = F) +
    geom_smooth(data = df.roc, aes(x = x.value, y = lr.neg, color = 'LR negative'), span = 0.5, size = 0.5, se = F) +
    ylim(0, 10) +
    xlab(xaxis) +
    ylab("likelihood ratio") +
    scale_color_discrete(name = "") +
    theme_bw() +
    theme(legend.position = c(0.35, 0.85), 
          legend.title=element_blank(),
          legend.background = element_rect(fill = alpha('white', 0.0)))
  
  return(list(plot1 = p1, plot2 = p2))
}

make.continuous.likelihoodratio.plot <- function(df.roc, df.raw, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df.roc, pretest)
  prob.diff.tmp <- prob.diff[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  df.roc <- df.roc[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  prob.diff <- prob.diff.tmp
  
  df.raw <- df.raw[df.raw$variable <= max(df.roc$x.value), ]
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
  lapply(seq(0.05, 0.95, length = 18), function(x) df.raw.tmp <<- generate.continuous.lr(df.raw.tmp, fit, x))
  
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

make.probability.plots <- function(df.raw, df.roc, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df.roc, pretest)
  prob.diff.tmp <- prob.diff[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  df.roc <- df.roc[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  prob.diff <- prob.diff.tmp
  
  p1 <- ggplot() + 
    #geom_vline(xintercept = 2) +
    geom_segment(aes(x = 0, y = pretest, xend = max(df.roc$x.value), yend = pretest, color = 'a'), linetype = 2) +
    geom_point(data = df.roc, aes(x = x.value, y = prob.diff$pos.prob, color = "b"), size = 0.5) +
    geom_point(data = df.roc, aes(x = x.value, y = prob.diff$neg.prob, color = "c"), size = 0.5) +
    geom_smooth(data = df.roc, aes(x = x.value, y = prob.diff$pos.prob, color = "b"), span = 0.5, size = 0.5, se = F) +
    geom_smooth(data = df.roc, aes(x = x.value, y = prob.diff$neg.prob, color = "c"), span = 0.5, size = 0.5, se = F) +
    geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'a')) +
    geom_point(data = data.frame(x = max(df.roc$x.value), y = pretest), aes(x = x, y = y, color = 'a')) +
    ylim(0, 1) + 
    xlab(xaxis) +
    ylab("probability pneumococcal") +
    scale_color_discrete(name = "", labels = c(a = 'pretest probability', b = 'posttest prob. positive', c = 'posttest prob. negative')) +
    theme_bw() +
    theme(legend.background = element_rect(fill=alpha('white', 0.0)))
  
  print(paste("Optimal Bayesian value: ", min(df.roc$x.value[(prob.diff$smoothed.pos.prob - prob.diff$smoothed.neg.prob) == max(prob.diff$smoothed.pos.prob - prob.diff$smoothed.neg.prob)])), sep = '')
  #print(max(prob.diff$pos.prob - prob.diff$neg.prob))
  
  p2 <- ggplot() + 
    geom_point(data = df.roc, aes(x = x.value, 
                              y = prob.diff$pos.prob - prob.diff$neg.prob, 
                              color = 'b'), 
               size = 0.5) +
    geom_smooth(data = df.roc, aes(x = x.value, 
                               y = prob.diff$pos.prob - prob.diff$neg.prob, 
                               color = 'b'), 
                span = 0.5, 
                size = 0.5,
                se = F) +
    geom_segment(data = df.roc, 
                 aes(x = min(x.value[(prob.diff$smoothed.pos.prob - prob.diff$smoothed.neg.prob) == max(prob.diff$smoothed.pos.prob - prob.diff$smoothed.neg.prob)]), 
                     y = 0, 
                     xend = min(x.value[(prob.diff$smoothed.pos.prob - prob.diff$smoothed.neg.prob) == max(prob.diff$smoothed.pos.prob - prob.diff$smoothed.neg.prob)]),
                     yend =  max((prob.diff$smoothed.pos.prob - prob.diff$smoothed.neg.prob)), 
                     color = 'a'),
                 size = 0.75) +
    geom_point(data = df.roc, 
               aes(x = min(x.value[(prob.diff$smoothed.pos.prob - prob.diff$smoothed.neg.prob) == max(prob.diff$smoothed.pos.prob - prob.diff$smoothed.neg.prob)]),
                   y =  max((prob.diff$smoothed.pos.prob - prob.diff$smoothed.neg.prob)),
                   color = 'a')) +
    ylim(0, 1) + 
    xlab(xaxis) +
    ylab("posttest probability diff") +
    scale_color_discrete(name = "", labels = c(a = 'optimal cutoff', b = 'probability difference')) +
    theme_bw()
  
  return(list(plot1 = p1, plot2 = p2, df.roc = df.roc, prob.diff = prob.diff$pos.prob - prob.diff$neg.prob))
}

make.combined.probability.plots <- function(df.raw, df.roc, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df.roc, pretest)
  prob.diff.tmp <- prob.diff[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  df.roc <- df.roc[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  prob.diff <- prob.diff.tmp
  
  df.raw <- df.raw[df.raw$variable <= max(df.roc$x.value), ]
  fit <- make.fit(df.raw)
  
  combined.probability <- (pretest / (1 - pretest)) * prob.diff$lr.pos * prob.diff$lr.neg / 
    ((pretest / (1 - pretest)) * prob.diff$lr.pos * prob.diff$lr.neg + 1)
  
  df.roc <- cbind(df.roc, combined.probability)
  
  p1 <- ggplot() + 
    geom_segment(aes(x = 0, y = pretest, xend = max(df.roc$x.value), yend = pretest, color = 'a'), linetype = 2) +
    geom_point(data = df.roc, aes(x = x.value, y = combined.probability, color = "b"), size = 0.5) +
    geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'a')) +
    geom_point(data = data.frame(x = max(df.roc$x.value), y = pretest), aes(x = x, y = y, color = 'a')) +
    ylim(0, 1) + 
    xlab(xaxis) +
    ylab("probability pneumococcal") +
    scale_color_discrete(name = "", labels = c(a = 'pretest probability', b = 'combined posttest')) +
    theme_bw() +
    theme(legend.position = c(0.29, 0.86), 
          legend.title=element_blank(),
          legend.background = element_rect(fill = alpha('white', 0.0)))
  
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
    scale_color_discrete(name = "", labels = c(a = 'pretest probability', b = 'continuous posttest')) +
    theme_bw() +
    theme(legend.position = c(0.29, 0.86), 
          legend.title=element_blank(),
          legend.background = element_rect(fill = alpha('white', 0.0)))
  
  p3 <- ggplot() + 
    geom_segment(aes(x = 0, y = pretest, xend = max(df.roc$x.value), yend = pretest, color = 'a'), linetype = 2) +
    geom_point(data = df.roc, aes(x = x.value, y = combined.probability, color = "b"), size = 0.5) +
    geom_point(data = df.raw.tmp, aes(x = x.value, y = y.value, color = "c"), size = 0.5) +
    geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'a')) +
    geom_point(data = data.frame(x = max(df.roc$x.value), y = pretest), aes(x = x, y = y, color = 'a')) +
    ylim(0, 1) + 
    xlab(xaxis) +
    ylab("probability pneumococcal") +
    scale_color_discrete(name = "", labels = c(a = 'pretest', b = 'combined', c = 'continuous')) +
    theme_bw() +
    theme(legend.position = c(0.2, 0.82), 
          legend.title=element_blank(),
          legend.background = element_rect(fill = alpha('white', 0.0)))
  
  return(list(plot1 = p1, plot2 = p2, plot3 = p3))
}

## Avoid global variable by running a main function
run.main.analysis <- function() {
  raw.pneumo <- import.data(variable = "Pneumococcal_diagnosis")
  pretest.probability <- calc.pretest(raw.pneumo, 'pneumo')
  
  ## lytA setup
  lytA.filename <- 'lytA.pdf'
  lytA.xaxis <- 'lytA density (log10 copies/mL)'
  raw.lytA <- import.data(variable = 'lytA_NP')
  
  ## Plot lytA
  roc.data.lytA <- make.test.stats.df(df = raw.lytA, which.variable = "pneumo")
  
  p.log.lytA <- make.logistic.plot(df = raw.lytA, df.roc = roc.data.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)
  p.roc.lytA <- make.roc.plot(df.roc = roc.data.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)
  
  p.sens.spec.lytA <- make.sens.spec.plot(df.roc = roc.data.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)
  p.lr.lytA <- make.likelihoodratio.plot(df.roc = roc.data.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)
  p.continuous.lr.lytA <- make.continuous.likelihoodratio.plot(df.roc = roc.data.lytA, df.raw = raw.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)
  p.prob.lytA <- make.probability.plots(df.raw = raw.lytA, df.roc = roc.data.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)
  p.combined.prob.lytA <- make.combined.probability.plots(df.raw = raw.lytA, df.roc = roc.data.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)
  
  #print(mean(p.prob.lytA$prob.diff[(p.prob.lytA[["df"]])[["x.value"]] < 6]))
  
  ## pct setup
  pct.filename <- 'pct.pdf'
  pct.xaxis <- 'procalcitonin concentration (ng/mL)'
  raw.pct <- import.data(variable = 'PCT')
  
  ## Plot pct
  roc.data.pct <- make.test.stats.df(df = raw.pct, which.variable = "pneumo")
  
  p.log.pct <- make.logistic.plot(df = raw.pct, df.roc = roc.data.pct, file.name = pct.filename, xaxis = pct.xaxis, pretest = pretest.probability)
  p.roc.pct <- make.roc.plot(df.roc = roc.data.pct, file.name = pct.filename, xaxis = pct.xaxis, pretest = pretest.probability)
  
  p.sens.spec.pct <- make.sens.spec.plot(roc.data.pct, file.name = pct.filename, xaxis = pct.xaxis, pretest = pretest.probability)
  p.lr.pct <- make.likelihoodratio.plot(df.roc = roc.data.pct, pct.filename, xaxis = pct.xaxis, pretest = pretest.probability)
  p.continuous.lr.pct <- make.continuous.likelihoodratio.plot(df.roc = roc.data.pct, df.raw = raw.pct, file.name = pct.filename, xaxis = pct.xaxis, pretest = pretest.probability)
  p.prob.pct <- make.probability.plots(df.raw = raw.pct, df.roc = roc.data.pct, file.name = pct.filename, xaxis = pct.xaxis, pretest = pretest.probability)
  p.combined.prob.pct <- make.combined.probability.plots(df.raw = raw.pct, df.roc = roc.data.pct, file.name = pct.filename, xaxis = pct.xaxis, pretest = pretest.probability)
  
  #print(mean(p.prob.pct$prob.diff[2 < (p.prob.pct[["df"]])[["x.value"]] & (p.prob.pct[["df"]])[["x.value"]] < 40]))
  
  ## crp setup
  crp.filename <- 'crp.pdf'
  crp.xaxis <- 'c-reactive protein concentration (mg/dL)'
  raw.crp <- import.data(variable = 'CRP')
  raw.crp$variable <- raw.crp$variable/10
  
  ## Plot crp
  roc.data.crp <- make.test.stats.df(df = raw.crp, which.variable = "pneumo")
  
  p.log.crp <- make.logistic.plot(df = raw.crp, df.roc = roc.data.crp, file.name = crp.filename, xaxis = crp.xaxis, pretest = pretest.probability)
  p.roc.crp <- make.roc.plot(df.roc = roc.data.crp, file.name = crp.filename, xaxis = crp.xaxis, pretest = pretest.probability)
  
  p.sens.spec.crp <- make.sens.spec.plot(roc.data.crp, file.name = crp.filename, xaxis = crp.xaxis, pretest = pretest.probability)
  p.lr.crp <- make.likelihoodratio.plot(df.roc = roc.data.crp, crp.filename, xaxis = crp.xaxis, pretest = pretest.probability)
  p.continuous.lr.crp <- make.continuous.likelihoodratio.plot(df.roc = roc.data.crp, df.raw = raw.crp, file.name = crp.filename, xaxis = crp.xaxis, pretest = pretest.probability)
  p.prob.crp <- make.probability.plots(df.raw = raw.crp, df.roc = roc.data.crp, file.name = crp.filename, xaxis = crp.xaxis, pretest = pretest.probability)
  p.combined.prob.crp <- make.combined.probability.plots(df.raw = raw.crp, df.roc = roc.data.crp, file.name = crp.filename, xaxis = crp.xaxis, pretest = pretest.probability)
  
  #print(mean(p.prob.crp$prob.diff[100 < (p.prob.crp[["df"]])[["x.value"]]]))
  
  p.log <- plot_grid(p.log.crp$plot, p.log.pct$plot, p.log.lytA$plot, labels = c('A', 'B', 'C'), ncol = 3)
  ggsave(plot = p.log, filename = "logistic_fits.pdf", height = 3, width = 15)
  
  p.roc <- plot_grid(p.roc.crp$plot1, p.roc.pct$plot1, p.roc.lytA$plot1, 
                     p.roc.crp$plot2, p.roc.pct$plot2, p.roc.lytA$plot2, 
                     p.roc.crp$plot3, p.roc.pct$plot3, p.roc.lytA$plot3, 
                     labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'), 
                     ncol = 3,
                     scale = 0.95)
  ggsave(plot = p.roc, filename = 'rocs.pdf', height = 10.5, width = 11.5)
  
  p.sens.spec <- plot_grid(p.sens.spec.crp$plot, p.sens.spec.pct$plot, p.sens.spec.lytA$plot, labels = c('A', 'B', 'C'), ncol = 3)
  ggsave(plot = p.sens.spec, filename = "sensitivities_specificities.pdf", height = 3, width = 15)
  
  p.lr <- plot_grid(p.lr.crp$plot1, p.lr.pct$plot1, p.lr.lytA$plot1, 
                    p.lr.crp$plot2,  p.lr.pct$plot2,p.lr.lytA$plot2, 
                    labels = c('A', 'B', 'C', 'D', 'E', 'F'), ncol = 3)
  ggsave(plot = p.lr, filename = "likelihood_ratios.pdf", height = 7.5, width = 12)
  
  p.continuous.lr <- plot_grid(p.continuous.lr.crp$plot1, p.continuous.lr.pct$plot1, p.continuous.lr.lytA$plot1, 
                               labels = c('A', 'B', 'C'), ncol = 3)
  ggsave(plot = p.continuous.lr, filename = "continuous_likelihood_ratios.pdf", height = 3, width = 10)
  
  p.prob <- plot_grid(p.prob.crp$plot1, p.prob.pct$plot1, p.prob.lytA$plot1, 
                      p.prob.crp$plot2, p.prob.pct$plot2, p.prob.lytA$plot2, 
                      labels = c('A', 'B', 'C', 'D', 'E', 'F'), ncol = 3)
  ggsave(plot = p.prob, filename = "probabilities.pdf", height = 6, width = 16)
  
  p.combined.prob <- plot_grid(#p.combined.prob.crp$plot1, p.combined.prob.pct$plot1, p.combined.prob.lytA$plot1, 
                               #p.combined.prob.crp$plot2, p.combined.prob.pct$plot2, p.combined.prob.lytA$plot2,
                               p.combined.prob.crp$plot3, p.combined.prob.pct$plot3, p.combined.prob.lytA$plot3,
                               labels = c('A', 'B', 'C'), ncol = 3)
  ggsave(plot = p.combined.prob, filename = "combined_probabilities.pdf", height = 3.5, width = 11.5)
}

run.two.variable.analysis <- function() {
  raw.pneumo <- import.data(variable = "Pneumococcal_diagnosis")
  pretest.probability <- calc.pretest(raw.pneumo, 'pneumo')
  
  ## lytA setup
  lytA.pct.filename <- 'lytA_pct.pdf'
  lytA.pct.xaxis <- 'lytA/PCT'
  raw.lytA.pct <- import.data(variable = 'PCT', variable2 = 'CRP')
  
  ## Plot lytA
  roc.data.lytA.pct <- make.test.stats.df(df = raw.lytA.pct, which.variable = "pneumo")
  
  p.log.lytA.pct <- make.logistic.plot(df = raw.lytA.pct, df.roc = roc.data.lytA.pct, file.name = lytA.pct.filename, xaxis = lytA.pct.xaxis, pretest = pretest.probability)
  p.roc.lytA.pct <- make.roc.plot(df = roc.data.lytA.pct, file.name = lytA.pct.filename, xaxis = lytA.pct.xaxis, pretest = pretest.probability)
  
  p.sens.spec.lytA.pct <- make.sens.spec.plot(df = roc.data.lytA.pct, file.name = lytA.pct.filename, xaxis = lytA.pct.xaxis, pretest = pretest.probability)
  p.lr.lytA.pct <- make.likelihoodratio.plot(df.roc = roc.data.lytA.pct, file.name = lytA.pct.filename, xaxis = lytA.pct.xaxis, pretest = pretest.probability)
  p.continuous.lr.lytA.pct <- make.continuous.likelihoodratio.plot(df = roc.data.lytA.pct, df.raw = raw.lytA.pct, file.name = lytA.pct.filename, xaxis = lytA.pct.xaxis, pretest = pretest.probability)
  p.prob.lytA.pct <- make.probability.plots(df.raw = raw.lytA.pct, df = roc.data.lytA.pct, file.name = lytA.pct.filename, xaxis = lytA.pct.xaxis, pretest = pretest.probability)
  p.combined.prob.lytA.pct <- make.combined.probability.plots(df.raw = raw.lytA.pct, df = roc.data.lytA.pct, file.name = lytA.pct.filename, xaxis = lytA.pct.xaxis, pretest = pretest.probability)
  
  show(p.roc.lytA.pct)
}


run.main.analysis()