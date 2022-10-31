SilhouetteScore <- function(dist, hclust = NULL) {
  score <- rep(NA, 20)
  if (!is.null(hclust)) {
    for (k in 2:20) {
      myclusters <- cutree(hclust, k = k)
      mysilhouette <- silhouette(x = myclusters, dist = dist)
      score[k] <- mean(mysilhouette[, 3])
    }
  } else {
    for (k in 2:20) {
      set.seed(1)
      mykmeans <- stats::kmeans(x = simul$data, centers = k)
      myclusters <- mykmeans$cluster
      mysilhouette <- silhouette(x = myclusters, dist = dist)
      score[k] <- mean(mysilhouette[, 3])
    }
  }
  return(score)
}


hclusCut <- function(x, k, d.meth = "euclidean", ...) {
  return(list(cluster = cutree(hclust(dist(x, method = d.meth), ...), k = k)))
}


kmeansCut <- function(x, k, d.meth = "euclidean", ...) {
  set.seed(1)
  return(list(cluster = stats::kmeans(x, centers = k)$cluster))
}


GapStatistic <- function(xdata, nc_max = 20, iters = 25, method = "hclust") {
  if (method == "hclust") {
    out <- clusGap(
      x = as.matrix(xdata),
      FUNcluster = hclusCut,
      K.max = nc_max, B = iters,
      verbose = FALSE
    )
  }
  if (method == "kmeans") {
    out <- clusGap(
      x = as.matrix(xdata),
      FUNcluster = kmeansCut,
      K.max = nc_max, B = iters,
      verbose = FALSE
    )
  }
  return(as.data.frame(out$Tab))
}


AreaUnderCDF <- function(M, thr_list = seq(0.01, 1, by = 0.01)) {
  n <- nrow(M)
  x <- M[upper.tri(M)]
  CDF <- rep(NA, length(thr_list))
  for (k in 1:length(thr_list)) {
    CDF[k] <- sum(x <= thr_list[k]) / (n * (n - 1) / 2)
  }

  area <- 0
  for (k in 2:length(CDF)) {
    area <- area + (thr_list[k] - thr_list[k - 1]) * CDF[k]
  }

  return(list(CDF = CDF, area = area))
}


DeltaArea <- function(areas) {
  delta <- areas[1]
  for (k in 2:(length(areas))) {
    delta <- c(delta, (areas[k] - areas[k - 1]) / areas[k - 1])
  }
  names(delta) <- names(areas)
  return(delta)
}


DeltaAreaCDF <- function(stability, thr_list = seq(0.01, 1, by = 0.01)) {
  areas <- rep(NA, dim(stability$coprop)[3])
  for (k in 2:dim(stability$coprop)[3]) {
    areas[k] <- AreaUnderCDF(M = stability$coprop[, , k])$area
  }
  delta <- DeltaArea(areas[-1])
  delta <- c(NA, delta)
  names(delta) <- stability$nc[, 1]
  return(delta)
}


PAC <- function(stab, x1 = 0.1, x2 = 0.9) {
  maxK <- max(stab$nc)
  Kvec <- 2:maxK
  score <- rep(NA, length(Kvec))
  for (i in Kvec) {
    M <- stab$coprop[, , i]
    Fn <- ecdf(M[lower.tri(M)])
    score[i - 1] <- Fn(x2) - Fn(x1)
  }
  score <- c(NA, score)
  names(score) <- stab$nc
  return(score)
}


MonteCarloScore <- function(x, stab, iters = 25, objective = "entropy") {
  # Running M3C for reference distribution
  out <- M3C(
    mydata = t(x),
    iters = iters,
    clusteralg = "hc",
    maxK = max(stab$nc),
    pItem = stab$params$tau,
    repsreal = 2,
    seed = 1,
    objective = objective,
    removeplots = TRUE,
    silent = TRUE
  )

  # Calculation of PAC scores
  if (objective == "PAC") {
    real <- data.frame(K = stab$nc, PAC_REAL = PAC(stab))
    real <- real[-1, ]
    rownames(real) <- 1:nrow(real)
  }

  # Calculation of entropy scores
  if (objective == "entropy") {
    entropies <- rep(NA, dim(stab$coprop)[3])
    for (i in 2:dim(stab$coprop)[3]) {
      entropies[i] <- M3C:::entropy(stab$coprop[, , i])
    }
    real <- data.frame(K = stab$nc, PAC_REAL = entropies)
    real <- real[-1, ]
    rownames(real) <- 1:nrow(real)
  }

  # Storing the reference
  ls <- out$refpacscores

  # Calculating reference mean
  real$PAC_REF <- colMeans(ls)

  # Checking PAC values
  ptemp <- real$PAC_REAL
  ptemp[ptemp == 0] <- 0.0001
  pacreal <- ptemp

  # Calculating RCSI score
  diffM <- sweep(log(ls), 2, log(pacreal))
  real$RCSI <- colMeans(diffM)
  real$RCSI_SE <- (apply(diffM, 2, sd)) / sqrt(nrow(ls))

  # Calculating p-value
  pvals <- vapply(seq_len(ncol(ls)), function(i) {
    distribution <- as.numeric(ls[, i])
    ((length(distribution[distribution < real$PAC_REAL[i]])) + 1) / (iters + 1) # (b+1)/(m+1)=pval
  }, numeric(1))
  real$MONTECARLO_P <- pvals

  if (objective == "PAC") {
    variance <- apply(ls, 2, var)
    pvals2 <- vapply(seq_len(nrow(real)), function(i) {
      mean <- real$PAC_REF[i]
      var <- variance[[i]]
      realpac <- real$PAC_REAL[i]
      params2 <- M3C:::estBetaParams(mu = mean, var = var)
      pbeta(realpac, params2[[1]], params2[[2]])
    }, numeric(1))
    real$BETA_P <- pvals2
    real$P_SCORE <- -log10(real$BETA_P)
  } else if (objective == "entropy") {
    variance <- apply(ls, 2, sd)
    pvals2 <- vapply(seq_len(nrow(real)), function(i) {
      mean <- real$PAC_REF[i]
      var <- variance[[i]]
      realpac <- real$PAC_REAL[i]
      pnorm(realpac, mean = mean, sd = var)
    }, numeric(1))
    real$NORM_P <- pvals2
    real$P_SCORE <- -log10(real$NORM_P)
    colnames(real)[2:3] <- c("ENTROPY_REAL", "ENTROPY_REF")
  }

  real <- rbind(c(1, rep(NA, ncol(real) - 1)), real)
  return(real)
}


AllPerf <- function(stab) {
  perf <- NULL
  for (i in 1:dim(stab$coprop)[3]) {
    perf <- rbind(
      perf,
      ClusteringPerformance(
        theta = Clusters(stab, argmax_id = i),
        theta_star = simul
      )
    )
  }
  perf <- cbind(nc = stab$nc, perf)
  return(perf)
}


ScatterPerf <- function(x, perf, xaxt = "s", xlab = "", ylab = "ARI", col = "navy") {
  id <- ManualArgmaxId(x)
  mycolours <- rep(col, nrow(perf))
  mycolours[id] <- "darkred"
  plot(x, perf$ari,
    panel.first = c(
      abline(h = perf$ari[id], col = "darkred", lty = 3),
      abline(v = x[id], col = "darkred", lty = 3)
    ),
    # panel.first=c(abline(h=axisTicks(range(perf$ari), log=FALSE), lty=3, col="grey"),
    #               abline(v=axisTicks(range(delta, na.rm=TRUE), log=FALSE), lty=3, col="grey")),
    pch = 19, cex = 3,
    las = 1,
    col = colorspace::lighten(mycolours, amount = 0.8),
    xaxt = xaxt,
    cex.lab = 1.5,
    xlab = xlab,
    ylab = ylab
  )
  text(x, perf$ari, labels = perf$nc, col = colorspace::darken(mycolours, amount = 0.2))
}


ManualCalibPlot <- function(y, x = NULL, xaxt = "s", xlab = "Number of clusters", ylab = "", col = "navy") {
  if (is.null(x)) {
    x <- 1:length(y)
  }
  id <- ManualArgmaxId(y)
  mycolours <- rep(col, length(y))
  mycolours[id] <- "darkred"
  plot(x, y,
    panel.first = c(
      abline(h = y[id], col = "darkred", lty = 3),
      abline(v = x[id], col = "darkred", lty = 3)
    ),
    # panel.first=c(abline(h=axisTicks(range(perf$ari), log=FALSE), lty=3, col="grey"),
    #               abline(v=axisTicks(range(delta, na.rm=TRUE), log=FALSE), lty=3, col="grey")),
    pch = 19, cex = 3,
    las = 1,
    col = colorspace::lighten(mycolours, amount = 0.8),
    xaxt = xaxt,
    cex.lab = 1.5,
    xlab = xlab,
    ylab = ylab
  )
  text(x, y, labels = x, col = colorspace::darken(mycolours, amount = 0.2))
}


ManualArgmaxId <- function(x, digits = 10) {
  return(max(which(round(x, digits = digits) == max(round(x, digits = digits), na.rm = TRUE))))
}


CalibrationCurve <- function(stability,
                             bty = "o", xlab = "", ylab = "",
                             col = c("navy", "forestgreen", "tomato"),
                             legend = TRUE, ncol = 1) {
  y <- stability$Sc
  x <- stability$nc
  z <- round(stability$Lambda, digits = 5)

  mycolours <- colorRampPalette(col)(length(unique(z)))
  names(mycolours) <- unique(z)
  plot(NA,
    xlim = c(0, max(stability$nc)), ylim = c(0, 1),
    xlab = xlab, ylab = ylab,
    las = 1, cex.lab = 1.5, bty = bty
  )
  for (lambda in unique(z)) {
    ids <- which(z == lambda)
    points(x[ids], y[ids], pch = 18, col = mycolours[as.character(lambda)])
    lines(x[ids], y[ids], lty = 1, lwd = 0.5, col = mycolours[as.character(lambda)])
  }
  abline(v = Argmax(stability)[1], lty = 2, col = "darkred")
  if (legend) {
    if (length(unique(stability$Q)) == 1) {
      legend("topright",
        legend = unique(formatC(stability$Lambda, format = "f", digits = 2)),
        pch = 15, col = mycolours, bty = "n", title = expression(lambda), ncol = ncol
      )
    } else {
      legend("topright",
        legend = paste0(
          unique(formatC(stability$Lambda[, 1], format = "f", digits = 2)),
          " (", unique(stability$Q[, 1]), ")"
        ),
        pch = 15, col = mycolours, bty = "n", title = expression(lambda), ncol = ncol
      )
    }
  }
}
