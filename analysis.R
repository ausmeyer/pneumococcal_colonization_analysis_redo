# Clear workspace and set working directory
rm(list = ls())
setwd('~/Google Drive/Documents/PostPostDoc/pneumococcal_colonization_analysis_redo/')

#library imports
libraries.call <- c("dplyr", "ggplot2", "readxl", "reshape", "epitools", 'glm')
lapply(libraries.call, require, character.only = TRUE)

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

left.censor <- function(df) {
  df <- data.frame(x = df$x, y = df$y)
  new.df <- data.frame()
  lapply(1:(nrow(df)-1), function(x) new.df <<- rbind(new.df, data.frame(x = df$x[x+1], y = df$y[x])))
  df <- rbind(df, new.df)
}
  
raw.dat <- read_excel('Database_BMJ_Open_25Jul2014.xlsx', col_types = c('text', 'text', 'numeric', 'numeric', 'numeric', 'text', 'text', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'text', 'text', 'text', 'text', 'text'))
raw.dat[, 1:17]
raw.dat <- data.frame(age = raw.dat$age, lytA = raw.dat$lytA_NP, bacter = raw.dat$bacteremia, pneumo = raw.dat$Pneumococcal_diagnosis)
raw.dat <- na.omit(raw.dat)

fit <- glm(factor(bacter) ~ lytA, data = raw.dat, family = "binomial")
predicted.bacteremia <- predict(fit, type="response", se = TRUE)

fit <- glm(factor(pneumo) ~ lytA, data = raw.dat, family = "binomial")
predicted.pneumo <- predict(fit, type="response", se = TRUE)

new.dat <- data.frame(lytA = raw.dat$lytA, bacter = raw.dat$bacter, prob.bacter = predicted.bacteremia$fit, pneumo = raw.dat$pneumo, prob.pneumo = predicted.pneumo$fit)

p.lytA <- ggplot() + 
  geom_point(data = new.dat, aes(x = lytA, y = as.numeric(as.character(pneumo)))) +
  geom_point(data = new.dat, aes(x = lytA, y = as.numeric(as.character(prob.pneumo)), color = 'red')) +
  theme_bw()
show(p.lytA)

#p.procal <- ggplot() + 
#  geom_point(data = new.dat, aes(x = pct, y = as.numeric(as.character(pneumo)))) +
#  geom_point(data = new.dat, aes(x = pct, y = as.numeric(as.character(prob.pneumo)), color = 'red')) +
#  theme_bw()
#show(p.procal)

lytA.series <- seq(0, max(raw.dat$lytA), by = 0.1)
test.stats <- data.frame()
lapply(lytA.series, function(x) test.stats <<- rbind(test.stats, (calculate.test.stats(generate.matrix(raw.dat, cutoff = x, variable = "lytA", response = "pneumo")))))

roc.dat <- data.frame(x = 1 - test.stats$spec, y = test.stats$sens)
roc.dat <- left.censor(roc.dat)
roc.dat <- roc.dat[order(roc.dat$y), ]

p <- ggplot() + 
  geom_point(data = roc.dat, aes(x = x, y = y)) +
  geom_line(data = roc.dat, aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw()
show(p)

reform <- data.frame(lytA = lytA.series, sens = test.stats$sens, spec = test.stats$spec)
reform.melted <- melt(reform, id = c("lytA"))
colnames(reform.melted) <- c('x', 'var', 'y')

reform.sens <- left.censor(filter(reform.melted, var == 'sens'))
reform.sens <- cbind(reform.sens, var = rep('sens', nrow(reform.sens)))
reform.sens <- reform.sens[rev(order(reform.sens$y)), ]

reform.spec <- left.censor(filter(reform.melted, var == 'spec'))
reform.spec <- cbind(reform.spec, var = rep('spec', nrow(reform.spec)))
reform.spec <- reform.spec[order(reform.spec$y), ]
reform.censored <- rbind(reform.sens, reform.spec)

p <- ggplot() + 
  geom_point(data = reform.censored, aes(x = x, y = y, color = var)) +
  geom_line(data = reform.censored, aes(x = x, y = y, color = var)) +
  theme_bw()
show(p)

pretest <- sum(raw.dat$pneumo == '1')/nrow(raw.dat)
reform <- reform[reform$sens != 0 & reform$spec != 0, ]

p <- ggplot() + 
  geom_point(data = reform, aes(x = lytA, y = sens / (1 - spec), color = 'LR Positive')) +
  geom_point(data = reform, aes(x = lytA, y = (1 - sens) / spec, color = 'LR Negative')) +
  geom_smooth(data = reform, aes(x = lytA, y = sens / (1 - spec), color = 'LR Positive')) +
  geom_smooth(data = reform, aes(x = lytA, y = (1 - sens) / spec, color = 'LR Negative')) +
  geom_vline(xintercept = 8) + 
  ylim(0,5) + 
  theme_bw()
show(p)

p <- ggplot() + 
  geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'Pretest Probability')) +
  geom_point(data = reform, aes(x = lytA, y = (pretest * (sens / (1 - spec)) / (1 - pretest)) / ((pretest * (sens / (1 - spec)) / (1 - pretest)) + 1), color = 'Positive Probability')) +
  geom_point(data = reform, aes(x = lytA, y = (pretest * ((1 - sens) / spec) / (1 - pretest)) / ((pretest * ((1 - sens) / spec) / (1 - pretest)) + 1), color = 'Negative Probability')) +
  geom_smooth(data = reform, aes(x = lytA, y = (pretest * (sens / (1 - spec)) / (1 - pretest)) / ((pretest * (sens / (1 - spec)) / (1 - pretest)) + 1), color = 'Positive Probability')) +
  geom_smooth(data = reform, aes(x = lytA, y = (pretest * ((1 - sens) / spec) / (1 - pretest)) / ((pretest * ((1 - sens) / spec) / (1 - pretest)) + 1), color = 'Negative Probability')) +
  geom_vline(xintercept = 8) + 
  ylim(0,1) + 
  theme_bw()
show(p)

#reform$sens <- rep(1.0, nrow(reform))
p <- ggplot() + 
  #geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'Pretest Probability')) +
  geom_point(data = reform, aes(x = lytA, y = (pretest * (sens / (1 - spec)) / (1 - pretest)) / ((pretest * (sens / (1 - spec)) / (1 - pretest)) + 1) - (pretest * ((1 - sens) / spec) / (1 - pretest)) / ((pretest * ((1 - sens) / spec) / (1 - pretest)) + 1), color = 'Probability Difference')) +
  geom_smooth(data = reform, aes(x = lytA, y = (pretest * (sens / (1 - spec)) / (1 - pretest)) / ((pretest * (sens / (1 - spec)) / (1 - pretest)) + 1) - (pretest * ((1 - sens) / spec) / (1 - pretest)) / ((pretest * ((1 - sens) / spec) / (1 - pretest)) + 1), color = 'Probability Difference')) +
  geom_vline(xintercept = 8) + 
  ylim(0,1) + 
  theme_bw()
show(p)