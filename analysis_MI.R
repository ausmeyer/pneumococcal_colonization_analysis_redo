## Clear workspace and set working directory
rm(list = ls())
setwd('~/Desktop/')

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
  
  return(data.frame(lr.pos = lr.pos, lr.neg = lr.neg, pos.prob = pos.prob, neg.prob = neg.prob))
}

generate.matrix <- function(df, cutoff, covariate, response) {
  t.positive <- sum(df[[covariate]] > cutoff & as.numeric(as.character(df[[response]]) == 1))
  f.positive <- sum(df[[covariate]] > cutoff & as.numeric(as.character(df[[response]]) == 0))
  f.negative <- sum(df[[covariate]] < cutoff & as.numeric(as.character(df[[response]]) == 1))
  t.negative <- sum(df[[covariate]] < cutoff & as.numeric(as.character(df[[response]]) == 0))
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

make.test.stats.df <- function(df, covariate, response) {
  df.series <- df[[covariate]][order(df[[covariate]])]
  df.stats <- data.frame()
  invisible(sapply(df.series, function(x) df.stats <<- rbind(df.stats, (calculate.test.stats(generate.matrix(df, cutoff = x, covariate = covariate, response = response))))))
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

import.data <- function(covariate) {
  df <- read_csv('MI_data.csv', 
                 col_types = cols(
                   age.group = col_character(),
                   sex = col_character(),
                   sbp = col_integer(),
                   dbp = col_integer(),
                   bahsTnT = col_double(),
                   fuhsTnT = col_double(),
                   mortality.ac = col_integer(),
                   MACE = col_integer(),
                   coronary.event = col_double(),
                   cva = col_double()
                 ))
  
  df <- cbind(df, diff.TnT = df$fuhsTnT - df$bahsTnT)
  df <- na.omit(df)
  df <- df[order(df[[covariate]]), ]
  df <- df[(nrow(df) * 0.025):(nrow(df) * 0.975), ]
  
  return(df)
}

make.fit <- function(df, covariate.variable, response.variable) {
  fit <- glm(factor(df[[response.variable]]) ~ df[[covariate.variable]], family = "binomial")

  predicted.response <- predict(fit, type = "response", se = TRUE)

  return(list(predicted = data.frame(response.variable.predicted = predicted.response), fit = fit))
}

make.logistic.plot <- function(df, file.name, xaxis, covariate.variable, response.variable) {

  fit.df <- make.fit(df = df, covariate.variable = covariate.variable, response.variable = response.variable)[['predicted']]

  new.df <- data.frame(variable = df[[covariate.variable]], response = df[[response.variable]], prob.response = fit.df$response.variable.predicted.fit)

  p <- ggplot() + 
    geom_point(data = new.df, aes(x = variable, y = as.numeric(as.character(response)), color = "diagnostic status")) +
    geom_line(data = new.df, aes(x = variable, y = as.numeric(as.character(prob.response)), color = 'logistic fit'), size = 1) +
    xlab(xaxis) +
    ylab(paste(response.variable, 'probability', sep = ' ')) +
    scale_color_discrete(name = "") +
    theme_bw()
  
  return(list(plot = p))
}

make.combined.probability.plots <- function(df.raw, df, file.name, xaxis, yaxis, pretest, covariate.variable, response.variable) {
  prob.diff <- calc.prob.difference(df, pretest)
  df <- df[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  prob.diff <- calc.prob.difference(df, pretest)
  df.raw <- df.raw[df.raw[[covariate.variable]] <= max(df$x.value), ]
  
  fit <- make.fit(df = df.raw, covariate.variable = covariate.variable, response.variable = response.variable)
  
  combined.probability <- (pretest / (1 - pretest)) * prob.diff$lr.pos * prob.diff$lr.neg / 
    ((pretest / (1 - pretest)) * prob.diff$lr.pos * prob.diff$lr.neg + 1)
  
  df <- cbind(df, combined.probability)

  p1 <- ggplot() + 
    geom_segment(aes(x = 0, y = pretest, xend = max(df$x.value), yend = pretest, color = 'a'), linetype = 2) +
    geom_point(data = df, aes(x = x.value, y = combined.probability, color = "b"), size = 0.5) +
    geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'a')) +
    geom_point(data = data.frame(x = max(df$x.value), y = pretest), aes(x = x, y = y, color = 'a')) +
    ylim(0, 1) + 
    xlab(xaxis) +
    ylab(yaxis) +
    scale_color_discrete(name = "", labels = c(a = 'pretest probability', b = 'combined posttest')) +
    theme_bw() +
    theme(legend.position = c(0.29, 0.86), 
          legend.title=element_blank(),
          legend.background = element_rect(fill = alpha('white', 0.0)))
  
  generate.continuous.lr <- function(tmp.df, tmp.fit, tmp.pretest) {
    x.1 <- (log(tmp.pretest / (1 - tmp.pretest)) - (tmp.fit[['fit']])[["coefficients"]][[1]]) / (tmp.fit[['fit']])[["coefficients"]][[2]]
    tmp.lr.continuous <- exp((tmp.fit[['fit']])[["coefficients"]][[2]] * (tmp.df[[covariate.variable]] - x.1))
    tmp.posttest <- (tmp.pretest / (1 - tmp.pretest)) * tmp.lr.continuous / (1 + (tmp.pretest / (1 - tmp.pretest)) * tmp.lr.continuous)
    tmp.df <- cbind(tmp.df, tmp = tmp.posttest)
    colnames(tmp.df)[colnames(tmp.df) == 'tmp'] <- 'posttest'
    return(tmp.df)
  }
  
  df.raw.tmp <- df.raw
  df.raw.tmp <- generate.continuous.lr(df.raw.tmp, fit, pretest)

  p2 <- ggplot() + 
    geom_segment(aes(x = 0, y = pretest, xend = max(df.raw.tmp[[covariate.variable]]), yend = pretest, color = 'a'), linetype = 2) +
    geom_point(data = df.raw.tmp, aes_string(x = covariate.variable, y = "posttest", color = "'b'"), size = 0.5) +
    geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'a')) +
    geom_point(data = data.frame(x = max(df.raw.tmp[[covariate.variable]]), y = pretest), aes(x = x, y = y, color = 'a')) +
    ylim(0, 1) + 
    xlab(xaxis) +
    ylab(yaxis) +
    scale_color_discrete(name = "", labels = c(a = 'pretest probability', b = 'continuous posttest')) +
    theme_bw()
    theme(legend.position = c(0.29, 0.86), 
          legend.title=element_blank(),
          legend.background = element_rect(fill = alpha('white', 0.0)))
  
  smoothed <- with(df, smooth.spline(x.value, combined.probability, df = 5))
  
  p3 <- ggplot() + 
    geom_segment(aes(x = 0, y = pretest, xend = max(df.raw.tmp[[covariate.variable]]), yend = pretest, color = 'a'), linetype = 2) +
    geom_point(data = df, aes(x = x.value, y = combined.probability, color = 'b'), size = 0.5) +
    geom_smooth(data = df, aes(x = x.value, y = combined.probability, color = 'b'), size = 0.5) +
    geom_point(data = df.raw.tmp, aes_string(x = covariate.variable, y = "posttest", color = "'c'"), size = 0.5) +
    geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'a')) +
    geom_point(data = data.frame(x = max(df.raw.tmp[[covariate.variable]]), y = pretest), aes(x = x, y = y, color = 'a')) +
    ylim(0, 1) + 
    xlab(xaxis) +
    ylab(yaxis) +
    scale_color_discrete(name = "", labels = c(a = 'pretest', b = 'combined', c = 'continuous')) +
    theme_bw() +
    theme(legend.position = c(0.2, 0.82), 
          legend.title=element_blank(),
          legend.background = element_rect(fill = alpha('white', 0.0)))
  
  show(p3)
  return(list(plot1 = p1, plot2 = p2, plot3 = p3))
}

## Avoid global variable by running a main function
run.main.analysis <- function() {
  response <- 'MACE'
  covariate <- 'bahsTnT'
  
  raw.df <- import.data(covariate)
  pretest.probability <- calc.pretest(raw.df, response)
  roc.data.df <- make.test.stats.df(df = raw.df, covariate = covariate, response = response)
  print(pretest.probability)
  
  p.coronary.event <- make.logistic.plot(df = raw.df, file.name = 'coronary_event_logistic.pdf', 
                                         xaxis = paste(covariate, '(pg/mL)'), 
                                         covariate.variable = covariate, 
                                         response.variable = response)
  
  show(p.coronary.event$plot)
  
  p.combined.prob.event <- make.combined.probability.plots(df.raw = raw.df, df = roc.data.df, 
                                                           file.name = 'coronary_event_combined.pdf', 
                                                           xaxis = paste(covariate, '(pg/mL)'), 
                                                           yaxis = paste('probability', response),
                                                           covariate.variable = covariate, 
                                                           response.variable = response,
                                                           pretest = pretest.probability)
  
  ggsave(plot = p.combined.prob.event$plot3, filename = "combined_probabilities.pdf", height = 3.5, width = 4)
}

run.main.analysis()