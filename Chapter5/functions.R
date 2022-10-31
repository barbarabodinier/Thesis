RunLinear <- function(covars, predictor, outcome, conf, multiple_x = TRUE) {
  if (multiple_x) {
    N <- ncol(predictor)
  } else {
    N <- ncol(outcome)
  }
  summary <- matrix(NA, ncol = 3, nrow = N)
  summary_formatted <- matrix(NA, ncol = 3, nrow = N)
  for (k in 1:N) {
    if (multiple_x) {
      mydatauniv <- cbind(
        outcome = covars[, outcome],
        covars[, conf, drop = FALSE],
        x = predictor[, k]
      )
    } else {
      mydatauniv <- cbind(
        outcome = outcome[, k],
        covars[, conf, drop = FALSE],
        x = covars[, predictor]
      )
    }
    mydatauniv <- na.exclude(mydatauniv)
    f <- paste0("outcome~x+", paste(conf, collapse = "+"))
    mymodel <- lm(as.formula(f), data = mydatauniv)
    summary[k, ] <- c(nrow(mydatauniv), coef(mymodel)["x"], summary(mymodel)$coefficients["x", 4])
    summary_formatted[k, ] <- c(nrow(mydatauniv), formatC(coef(mymodel)["x"], format = "f", digits = 2),
      pval = formatC(summary(mymodel)$coefficients["x", 4], format = "e", digits = 2)
    )
  }
  colnames(summary) <- colnames(summary_formatted) <- c("nobs", "coef", "pval")
  if (multiple_x) {
    rownames(summary) <- rownames(summary_formatted) <- colnames(predictor)
  } else {
    rownames(summary) <- rownames(summary_formatted) <- colnames(outcome)
  }
  summary <- as.data.frame(summary)
  summary_formatted <- as.data.frame(summary_formatted)

  return(list(summary = summary, formatted = summary_formatted))
}


RunLogistic <- function(covars, predictors, outcome, conf) {
  summary <- summary_formatted <- matrix(NA, ncol = 3, nrow = ncol(predictors))
  for (k in 1:ncol(predictors)) {
    # print(k)
    mydatauniv <- cbind(covars[, c(outcome, conf)], x = predictors[, k])
    mydatauniv <- na.exclude(mydatauniv)
    mydatauniv$outcome <- as.numeric(eval(parse(text = paste0("mydatauniv$", outcome))))
    mydatauniv$x <- scale(mydatauniv$x)
    f <- paste0("outcome~x+", paste(conf, collapse = "+"))
    mymodel <- glm(as.formula(f), data = mydatauniv, family = "binomial")
    f0 <- paste0("outcome~", paste(conf, collapse = "+"))
    mymodel0 <- glm(as.formula(f0), data = mydatauniv, family = "binomial")
    myanova <- anova(mymodel0, mymodel, test = "LRT")
    summary[k, ] <- c(nrow(mydatauniv), exp(coef(mymodel)["x"]),
      pval = myanova$`Pr(>Chi)`[2]
    )
    summary_formatted[k, ] <- c(nrow(mydatauniv), formatC(exp(coef(mymodel)["x"]), format = "f", digits = 2),
      pval = formatC(myanova$`Pr(>Chi)`[2], format = "e", digits = 2)
    )
  }
  colnames(summary) <- colnames(summary_formatted) <- c("nobs", "coef", "pval")
  rownames(summary) <- rownames(summary_formatted) <- colnames(predictors)
  summary <- as.data.frame(summary)
  summary_formatted <- as.data.frame(summary_formatted)
  return(list(summary = summary, formatted = summary_formatted))
}


SubtypeResampling <- function(data, tau, strata, ...) {
  s <- NULL
  for (z in unique(strata)) {
    s <- c(s, sample(which(strata == z), size = tau * sum(strata == z)))
    s <- c(s, sample(which(strata == z), size = tau * sum(strata == z)))
  }
  return(s)
}
