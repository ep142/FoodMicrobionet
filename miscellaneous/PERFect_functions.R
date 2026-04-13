# PERFect functions v 0.2 14/12/2024

# repository at https://github.com/cxquy91/PERFect
# library(devtools)
# install_github("cxquy91/PERFect")
# original author Ekaterina Smirnova <ekaterina.smirnova@vcuhealth.org>
# Smirnova E, Cao Q (2022). PERFect: Permutation filtration for microbiome data. R package version 1.12.0, https://github.com/cxquy91/PERFect.

# X is OTU table, where taxa are columns and samples are rows of the table. 
# It should be a in dataframe format with columns corresponding to taxa names.
# install/load packages 
.scran_packages <- c("sn", "parallel", "Matrix", "fitdistrplus", "zoo", "psych")
.sbioc_packages <- c("BiocManager", "phyloseq")
.attached_pckgs <- .packages()

# install and load Bioconductor packages if necessary
if(!all(.sbioc_packages %in% .attached_pckgs)){
  .sinst <- .sbioc_packages %in% installed.packages()
  if(any(!.sinst)) {
    if(!.sinst[1]) {
      install.packages("BiocManager")
      .sinst <- .sbioc_packages %in% installed.packages()
    }
    if(any(!.sinst[2:length(.sinst)])) {
      BiocManager::install(.sbioc_packages[!.sinst], ask = F)
    }
  }
  sapply(.sbioc_packages, require, character.only = TRUE)
}
if(!("sn" %in% installed.packages())) {
  install.packages("sn")
}
require(sn)

DiffFiltLoss <- function (X, Order_Ind, Plot = TRUE, Taxa_Names = NULL) 
{
  if (!is(X, "matrix")) {
    X <- as.matrix(X)
  }
  if (!is(Order_Ind, "integer")) 
    stop("Order_Ind argument must be a vector of integer values")
  if (!is(Plot, "logical")) 
    stop("Plot argument must be a logical value")
  p <- dim(X)[2]
  DFL <- rep(1, p - 1)
  X <- as.matrix(X)
  Netw <- t(X) %*% X
  for (j in seq_len(p - 1)) {
    DFL[j] <- DiffFiltLoss_j(Order_Ind, Netw, j)
  }
  DFL <- DFL/sum(Netw * Netw)
  if (Plot == TRUE) {
    df <- data.frame(Order_Ind[-1], rep(seq_len(length(DFL))), 
                     DFL)
    names(df)[2] <- "x"
    Lab <- seq_len(length(Order_Ind[-1]))
    df <- cbind(Lab, df)
    p_FL <- ggplot(df) + geom_line(aes(x = Lab, y = DFL, 
                                       group = 1), colour = "dodgerblue3") + theme(panel.background = element_rect(fill = "white"), 
                                                                                   panel.grid.major = element_line(colour = "grey90"), 
                                                                                   axis.text.x = element_text(size = 10, colour = "black", 
                                                                                                              angle = 90, hjust = 1)) + ggtitle("") + ylab("Differences in Filtering Loss") + 
      xlab("Taxa") + xlim(0, max(df$Lab))
  }
  if (!is.null(Taxa_Names)) {
    names(DFL) <- Taxa_Names[-1]
  }
  return(list(DFL = DFL, p_FL = p_FL))
}

FiltLoss <- function (X, Order = "NP", Order.user = NULL, type = "Cumu", 
                      Plot = TRUE) 
{
  p <- dim(X)[2]
  Norm_Ratio <- rep(1, p)
  if (!is(X, "matrix")) {
    X <- as.matrix(X)
  }
  if (!(Order %in% c("NP", "NC", "NCw"))) 
    stop("Order argument can only be \"NP\", \"NC\", or \"NCw\" ")
  if (!(type %in% c("Cumu", "Ind"))) 
    stop("type argument can only be \"Cumu\" or \"Ind\" ")
  if (!is(Plot, "logical")) 
    stop("Plot argument must be a logical value")
  if (is.null(Order.user)) {
    if (Order == "NP") {
      Order.vec <- NP_Order(X)
    }
    if (Order == "NC") {
      Order.vec <- NC_Order(X)
    }
    if (Order == "NCw") {
      Order.vec <- NCw_Order(X)
    }
  }
  else {
    Order.vec <- Order.user
  }
  X <- X[, Order.vec]
  Order_Ind <- seq_len(length(Order.vec))
  Netw <- t(X) %*% X
  for (i in seq_len(p)) {
    if (type == "Cumu") {
      Ind <- Order_Ind[-seq_len(i)]
    }
    else {
      Ind <- Order_Ind[-i]
    }
    Netw_R <- Netw[Ind, Ind]
    Norm_Ratio[i] <- sum(Netw_R * Netw_R)
  }
  FL <- 1 - Norm_Ratio/sum(Netw * Netw)
  if (Plot == TRUE) {
    df <- data.frame(Order.vec, rep(seq_len(length(FL))), 
                     FL)
    names(df)[2] <- "x"
    Lab <- seq_len(length(Order.vec))
    df <- cbind(Lab, df)
    p_FL <- ggplot(df) + geom_line(aes(x = Lab, y = FL, group = 1), 
                                   colour = "dodgerblue3") + theme(panel.background = element_rect(fill = "white"), 
                                                                   panel.grid.major = element_line(colour = "grey90"), 
                                                                   axis.text.x = element_text(size = 10, colour = "black", 
                                                                                              angle = 90, hjust = 1)) + ggtitle("") + ylab("Filtering Loss") + 
      xlab("Taxa") + xlim(0, max(df$Lab))
  }
  names(FL) <- colnames(X)
  return(list(FL = FL, p_FL = p_FL))
}

FL_J <- function (X, J) 
{
  if (!is(X, "matrix")) {
    X <- as.matrix(X)
  }
  if (!is(J, "character")) 
    stop("J argument must be a character vector containing names of taxa to be removed")
  Ind <- which(colnames(X) %in% J)
  X_R <- X[, -Ind]
  Netw <- t(as.matrix(X)) %*% as.matrix(X)
  Netw_R <- t(as.matrix(X_R)) %*% as.matrix(X_R)
  FL <- 1 - (sum(Netw_R * Netw_R)/sum(Netw * Netw))
  return(FL)
}



PERFect_perm <- function (X, infocol = NULL, Order = "NP", Order.user = NULL, 
                          normalize = "counts", algorithm = "fast", center = FALSE, 
                          quant = c(0.1, 0.25, 0.5), distr = "sn", alpha = 0.1, rollmean = TRUE, 
                          direction = "left", pvals_sim = NULL, k = 10000, nbins = 30, 
                          hist = TRUE, col = "red", fill = "green", hist_fill = 0.2, 
                          linecol = "blue") 
{
  info <- NULL
  if (!is.null(infocol)) {
    info <- X[, infocol]
    X <- X[, -infocol]
  }
  if (!is(X, "matrix")) {
    X <- as.matrix(X)
  }
  if (!(Order %in% c("NP", "pvals", "NC", "NCw"))) 
    stop("Order argument can only be \"NP\", \"pvals\", \"NC\", or \"NCw\" ")
  if (!(normalize %in% c("counts", "prop", "pres"))) 
    stop("normalize argument can only be \"counts\", \"prop\", or \"pres\" ")
  if (!is(center, "logical")) 
    stop("center argument must be a logical value")
  if (!is.vector(quant)) 
    stop("quant argument must be a vector")
  if (!(distr %in% c("sn", "norm", "t", "cauchy"))) 
    stop("normalize argument can only be \"sn\", \"norm\", \"t\", or \"cauchy\" ")
  if (!is.numeric(alpha)) 
    stop("alpha argument must be a numerical value")
  if (!is(pvals_sim, "NULL") & length(pvals_sim$pvals) == 0) 
    stop("pvals_sim object must be a result from simultaneous PERFect with taxa abundance ordering")
  if (is.null(Order.user)) {
    if (Order == "NP") {
      Order.vec <- NP_Order(X)
    }
    if (Order == "pvals") {
      Order.vec <- pvals_Order(X, pvals_sim)
    }
    if (Order == "NC") {
      Order.vec <- NC_Order(X)
    }
    if (Order == "NCw") {
      Order.vec <- NCw_Order(X)
    }
  }
  else {
    Order.vec <- Order.user
  }
  X <- X[, Order.vec]
  nzero.otu <- apply(X, 2, Matrix::nnzero) != 0
  X <- X[, nzero.otu]
  p <- dim(X)[2]
  Order.vec <- Order.vec[nzero.otu]
  X.orig <- X
  if (normalize == "prop") {
    X <- X/apply(X, 1, sum)
  }
  else if (normalize == "pres") {
    X[X != 0] <- 1
  }
  if (center) {
    X <- apply(X, 2, function(x) {
      x - mean(x)
    })
  }
  if (algorithm == "fast") {
    n <- sum_n(dim(X)[2])$idx - 1
    pvals <- rep(NA, p - 1)
    Order_Ind <- rep(seq_len(length(Order.vec)))
    DFL <- DiffFiltLoss(X = X, Order_Ind, Plot = TRUE, Taxa_Names = Order.vec)
    FL <- FiltLoss(X = X, Order.user = Order.vec, type = "Ind", 
                   Plot = TRUE)$FL
    names(pvals) <- names(DFL$DFL)
    dfl_distr <- sampl_distr(X = X, k = k)
    lfl <- lapply(dfl_distr, function(x) log(x[!x == 0]))
    check <- c()
    no_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(no_cores)
    if (distr == "norm") {
      parallel::clusterEvalQ(cl, {
        library(fitdistrplus)
      })
      lfl_main <- lfl[n]
      if (length(quant) > 2) {
        quant <- quant[(length(quant) - 1):length(quant)]
        print("Warning: more than 2 quantile values are given. \nLargest 2 quantiles are used.")
      }
      if (length(quant) < 2) {
        stop("At least two quantile values must be specified.")
      }
      parallel::clusterExport(cl, c("distr", "quant"), 
                              envir = environment())
      fit <- parallel::parLapply(cl, lfl_main, function(x) fitdistrplus::qmedist(x, 
                                                                                 distr, probs = quant))
      est <- lapply(fit, function(x) x$estimate)
      for (i in seq_len(length(n))) {
        pvals[n[i]] <- pnorm(q = log(DFL$DFL[n[i]]), 
                             mean = est[[i]][1], sd = est[[i]][2], lower.tail = FALSE, 
                             log.p = FALSE)
        if (pvals[n[i]] < alpha) {
          stop_idx <- n[i - 1]
          temp <- c((n[i - 1] + 1):(n[i] - 1), n[i] + 
                      1, n[i] + 2)
          lfl_temp <- lfl[temp]
          fit <- parallel::parLapply(cl, lfl_temp, function(x) fitdistrplus::qmedist(x, 
                                                                                     distr, probs = quant))
          est <- lapply(fit, function(x) x$estimate)
          pvals[temp] <- mapply(function(dat1, dat2) {
            vapply(dat1, function(z) 1 - pnorm(q = dat1, 
                                               mean = dat2[1], sd = dat2[2], lower.tail = FALSE, 
                                               log.p = FALSE), numeric(1))
          }, dat1 = log(DFL$DFL[temp]), est)
          break
        }
      }
      potential <- which(pvals < alpha)[1]
      for (i in round((match(stop_idx, n) * 9/10)):match(stop_idx, 
                                                         n)) {
        check <- c(check, names(which(DFL$DFL[n[i]:n[(i + 
                                                        1)]] > max(DFL$DFL[n[i + 1]], DFL$DFL[n[i]]) & 
                                        FL[n[i]:n[(i + 1)]] > max(FL[n[i + 1]], FL[n[i]]))))
      }
      check <- c(check, names(which(DFL$DFL[seq_len(potential - 
                                                      1)] > DFL$DFL[potential] & FL[seq_len(potential - 
                                                                                              1)] > FL[potential])))
      check <- check[is.na(pvals[which(names(pvals) %in% 
                                         check)])]
      if (length(check) != 0) {
        check_idx <- which(names(pvals) %in% check)
        lfl_check <- lfl[check_idx]
        fit <- parallel::parLapply(cl, lfl_check, function(x) fitdistrplus::qmedist(x, 
                                                                                    distr, probs = quant))
        est <- lapply(fit, function(x) x$estimate)
        pvals[check_idx] <- mapply(function(dat1, dat2) {
          vapply(dat1, function(z) 1 - pnorm(q = dat1, 
                                             mean = dat2[1], sd = dat2[2], lower.tail = FALSE, 
                                             log.p = FALSE), numeric(1))
        }, dat1 = log(DFL$DFL[check_idx]), est)
      }
    }
    if (distr == "sn") {
      parallel::clusterEvalQ(cl, {
        library(fitdistrplus)
        library(sn)
      })
      lp <- list(xi = mean(lfl[[1]]), omega = sd(lfl[[1]]), 
                 alpha = 1.5)
      lfl_main <- lfl[n]
      parallel::clusterExport(cl, c("distr", "quant", "lp"), 
                              envir = environment())
      suppressWarnings(fit <- parallel::parLapply(cl, lfl_main, 
                                                  function(x) fitdistrplus::qmedist(x, distr, probs = quant, 
                                                                                    start = lp)))
      est <- lapply(fit, function(x) x$estimate)
      for (i in seq_len(length(n))) {
        pvals[n[i]] <- 1 - psn(x = log(DFL$DFL[n[i]]), 
                               xi = est[[i]][1], omega = est[[i]][2], alpha = est[[i]][3])
        if (pvals[n[i]] < alpha) {
          stop_idx <- n[i - 1]
          temp <- c((n[i - 1] + 1):(n[i] - 1), n[i] + 
                      1, n[i] + 2)
          lfl_temp <- lfl[temp]
          suppressWarnings(fit <- parallel::parLapply(cl, 
                                                      lfl_temp, function(x) fitdistrplus::qmedist(x, 
                                                                                                  distr, probs = quant, start = lp)))
          est <- lapply(fit, function(x) x$estimate)
          pvals[temp] <- mapply(function(dat1, dat2) {
            vapply(dat1, function(z) 1 - psn(x = dat1, 
                                             xi = dat2[1], omega = dat2[2], alpha = dat2[3]), 
                   numeric(1))
          }, dat1 = log(DFL$DFL[temp]), est)
          break
        }
      }
      potential <- which(pvals < alpha)[1]
      for (i in round((match(stop_idx, n) * 9/10)):match(stop_idx, 
                                                         n)) {
        check <- c(check, names(which(DFL$DFL[n[i]:n[(i + 
                                                        1)]] > max(DFL$DFL[n[i + 1]], DFL$DFL[n[i]]) & 
                                        FL[n[i]:n[(i + 1)]] > max(FL[n[i + 1]], FL[n[i]]))))
      }
      check <- c(check, names(which(DFL$DFL[seq_len(potential - 
                                                      1)] > DFL$DFL[potential] & FL[seq_len(potential - 
                                                                                              1)] > FL[potential])))
      check <- check[is.na(pvals[which(names(pvals) %in% 
                                         check)])]
      if (length(check) != 0) {
        check_idx <- which(names(pvals) %in% check)
        lfl_check <- lfl[check_idx]
        suppressWarnings(fit <- parallel::parLapply(cl, 
                                                    lfl_check, function(x) fitdistrplus::qmedist(x, 
                                                                                                 distr, probs = quant, start = lp)))
        est <- lapply(fit, function(x) x$estimate)
        pvals[check_idx] <- mapply(function(dat1, dat2) {
          vapply(dat1, function(z) 1 - psn(x = dat1, 
                                           xi = dat2[1], omega = dat2[2], alpha = dat2[3]), 
                 numeric(1))
        }, dat1 = log(DFL$DFL[check_idx]), est)
      }
    }
    parallel::stopCluster(cl)
    if (rollmean) {
      non_na_ind <- which(!is.na(pvals))
      pvals_avg <- pvals
      pvals_avg[non_na_ind] <- zoo::rollapply(pvals[non_na_ind], 
                                              width = 3, mean, align = direction, fill = NA, 
                                              na.rm = TRUE)
      pvals_avg[non_na_ind][is.na(pvals_avg[non_na_ind])] <- pvals[non_na_ind][is.na(pvals_avg[non_na_ind])]
      names(pvals_avg) <- names(pvals)[(length(pvals) - 
                                          length(pvals_avg) + 1):length(pvals)]
    }
    else {
      pvals_avg <- pvals
    }
    Ind <- which(pvals_avg <= alpha & pvals_avg > 0)
    if (length(Ind) == 0) {
      Ind <- which(pvals <= alpha & pvals > 0)
      warning("Rolling average greatly modifies this result. Try one of these suggestions: \n\n              1. Increase the number of permutations \n\n              2. Set rollingmean to FALSE ")
    }
    if (length(Ind != 0)) {
      Ind <- min(Ind)
    }
    else {
      Ind <- dim(X)[2] - 1
      warning("No taxa are significant at a specified alpha level.")
    }
    filtX <- X.orig[, -seq_len(Ind)]
    return(list(filtX = filtX, info = info, fit = NULL, hist = NULL, 
                est = NULL, dfl_distr = NULL, pvals = pvals_avg))
  }
  else if (algorithm == "full") {
    pvals <- rep(0, p - 1)
    hist_list <- lapply(seq_len(p - 1), function(x) NULL)
    Order_Ind <- rep(seq_len(length(Order.vec)))
    DFL <- DiffFiltLoss(X = X, Order_Ind, Plot = TRUE, Taxa_Names = Order.vec)
    names(pvals) <- names(DFL$DFL)
    dfl_distr <- sampl_distr(X = X, k = k)
    lfl <- lapply(dfl_distr, function(x) log(x[!x == 0]))
    if (distr == "norm") {
      if (length(quant) > 2) {
        quant <- quant[(length(quant) - 1):length(quant)]
        print("Warning: more than 2 quantile values are given. \nLargest 2 quantiles are used.")
      }
      if (length(quant) < 2) {
        stop("At least two quantile values must be specified.")
      }
      fit <- lapply(lfl, function(x) fitdistrplus::qmedist(x, 
                                                           distr, probs = quant))
      est <- lapply(fit, function(x) x$estimate)
      pvals <- mapply(function(dat1, dat2) {
        vapply(dat1, function(z) 1 - pnorm(q = dat1, 
                                           mean = dat2[1], sd = dat2[2], lower.tail = FALSE, 
                                           log.p = FALSE), numeric(1))
      }, dat1 = log(DFL$DFL), est)
    }
    if (distr == "sn") {
      lp <- list(xi = mean(lfl[[1]]), omega = sd(lfl[[1]]), 
                 alpha = 1.5)
      suppressWarnings(fit <- lapply(lfl, function(x) fitdistrplus::qmedist(x, 
                                                                            distr, probs = quant, start = lp)))
      est <- lapply(fit, function(x) x$estimate)
      pvals <- mapply(function(dat1, dat2) {
        vapply(dat1, function(z) 1 - psn(x = dat1, xi = dat2[1], 
                                         omega = dat2[2], alpha = dat2[3]), numeric(1))
      }, dat1 = log(DFL$DFL), est)
    }
    if (hist == TRUE) {
      lfl <- lapply(lfl, function(x) data.frame(x))
      for (i in seq_len(p - 1)) {
        lfl <- data.frame(log(dfl_distr[[i]][!dfl_distr[[i]] == 
                                               0]))
        if (length(dfl_distr[[i]][dfl_distr[[i]] == 0]) > 
            0) {
          print(paste("taxon", i, "number of zeroes = ", 
                      length(dfl_distr[[i]][dfl_distr[[i]] == 0])))
        }
        names(lfl) <- c("DFL")
        x <- "DFL"
        ord_map <- aes_string(x = x)
        hist <- ggplot(lfl, ord_map) + geom_histogram(bins = nbins, 
                                                      aes(y = ..density..), col = col, fill = fill, 
                                                      alpha = hist_fill) + theme(panel.background = element_rect(fill = "white"), 
                                                                                 panel.grid.major = element_line(colour = "grey90"), 
                                                                                 axis.text.x = element_text(size = 10)) + ggtitle("") + 
          xlab("log differences in filtering loss") + 
          ylab("Density")
        if (distr == "norm") {
          hist <- hist + stat_function(fun = dnorm, args = list(mean = est[[i]][1], 
                                                                sd = est[[i]][2]), colour = linecol)
          hist_list[[i]] <- hist
        }
        if (distr == "sn") {
          hist <- hist + stat_function(fun = dsn, args = list(xi = est[[i]][1], 
                                                              omega = est[[i]][2], alpha = est[[i]][3]), 
                                       colour = linecol)
          hist_list[[i]] <- hist
        }
      }
    }
    if (rollmean) {
      pvals_avg <- zoo::rollapply(pvals, width = 3, mean, 
                                  align = direction, fill = NA, na.rm = TRUE)
      pvals_avg[is.na(pvals_avg)] <- pvals[is.na(pvals_avg)]
      names(pvals_avg) <- names(pvals)
    }
    else {
      pvals_avg <- pvals
    }
    Ind <- which(pvals_avg <= alpha)
    if (length(Ind) == 0) {
      Ind <- which(pvals <= alpha & pvals > 0)
    }
    if (length(Ind != 0)) {
      Ind <- min(Ind)
    }
    else {
      Ind <- dim(X)[2] - 1
      warning("no taxa are significant at a specified alpha level")
    }
    filtX <- X.orig[, -seq_len(Ind)]
    return(list(filtX = filtX, info = info, fit = fit, hist = hist_list, 
                est = est, dfl_distr = dfl_distr, pvals = pvals_avg))
  }
}

PERFect_perm_reorder<- function (X, Order = "NP", Order.user = NULL, res_perm, normalize = "counts", 
                                 center = FALSE, alpha = 0.1, distr = "sn", rollmean = TRUE, 
                                 direction = "left", pvals_sim = NULL) 
{
  if (!is(X, "matrix")) {
    X <- as.matrix(X)
  }
  if (!(Order %in% c("NP", "pvals", "NC", "NCw"))) 
    stop("Order argument can only be \"NP\", \"pvals\", \"NC\", or \"NCw\" ")
  if (!is(res_perm, "NULL") & length(res_perm$pvals) == 0) 
    stop("res_perm argument must be the output from the function PERFect_perm()")
  if (!(normalize %in% c("counts", "prop", "pres"))) 
    stop("normalize argument can only be \"counts\", \"prop\", or \"pres\" ")
  if (!is(center, "logical")) 
    stop("center argument must be a logical value")
  if (!(distr %in% c("sn", "norm", "t", "cauchy"))) 
    stop("normalize argument can only be \"sn\", \"norm\", \"t\", or \"cauchy\" ")
  if (!is.numeric(alpha)) 
    stop("alpha argument must be a numerical value")
  if (Order == "NP") {
    Order.vec <- NP_Order(X)
  }
  if (Order == "pvals") {
    Order.vec <- pvals_Order(X, pvals_sim)
  }
  if (Order == "NC") {
    Order.vec <- NC_Order(X)
  }
  if (Order == "NCw") {
    Order.vec <- NCw_Order(X)
  }
  else if (!is.null(Order.user)) {
    Order.vec = Order.user
  }
  X <- X[, Order.vec]
  X.orig <- X
  nzero.otu <- apply(X, 2, Matrix::nnzero) != 0
  X <- X[, nzero.otu]
  p <- dim(X)[2]
  Order.vec <- Order.vec[nzero.otu]
  if (normalize == "prop") {
    X <- X/apply(X, 1, sum)
  }
  else if (normalize == "pres") {
    X[X != 0] <- 1
  }
  if (center) {
    X <- apply(X, 2, function(x) {
      x - mean(x)
    })
  }
  Order_Ind <- rep(seq_len(length(Order.vec)))
  DFL <- DiffFiltLoss(X = X, Order_Ind, Plot = TRUE, Taxa_Names = Order.vec)
  pvals <- rep(0, length(DFL$DFL))
  names(pvals) <- names(DFL$DFL)
  if (distr == "sn") {
    pvals <- mapply(function(dat1, dat2) {
      vapply(dat1, function(z) 1 - psn(x = dat1, xi = dat2[1], 
                                       omega = dat2[2], alpha = dat2[3]), numeric(1))
    }, dat1 = log(DFL$DFL), res_perm$est)
  }
  if (distr == "norm") {
    pvals <- mapply(function(dat1, dat2) {
      vapply(dat1, function(z) 1 - pnorm(q = dat1, mean = dat2[1], 
                                         sd = dat2[2], lower.tail = FALSE, log.p = FALSE), 
             numeric(1))
    }, dat1 = log(DFL$DFL), res_perm$est)
  }
  if (rollmean) {
    pvals_avg <- zoo::rollmean(pvals, k = 3, align = direction, 
                               fill = NA)
  }
  else {
    pvals_avg <- pvals
  }
  pvals_avg[is.na(pvals_avg)] <- pvals[is.na(pvals_avg)]
  Ind <- which(pvals_avg <= alpha)
  if (length(Ind != 0)) {
    Ind <- min(Ind)
  }
  else {
    Ind <- dim(X)[2] - 1
    warning("no taxa are significant at a specified alpha level")
  }
  res_perm$filtX <- X.orig[, -seq_len(Ind)]
  res_perm$pvals <- pvals_avg
  return(res_perm)
}

PERFect_sim <- function (X, infocol = NULL, Order = "NP", Order.user = NULL, 
                         normalize = "counts", center = FALSE, quant = c(0.1, 0.25, 
                                                                         0.5), distr = "sn", alpha = 0.1, rollmean = TRUE, direction = "left", 
                         pvals_sim = NULL, nbins = 30, col = "red", fill = "green", 
                         hist_fill = 0.2, linecol = "blue") 
{
  pDFL <- NULL
  phist <- NULL
  info <- NULL
  if (!is.null(infocol)) {
    info <- X[, infocol]
    X <- X[, -infocol]
  }
  if (!is(X, "matrix")) {
    X <- as.matrix(X)
  }
  if (!(Order %in% c("NP", "pvals", "NC", "NCw"))) 
    stop("Order argument can only be \"NP\", \"pvals\", \"NC\", or \"NCw\" ")
  if (!(normalize %in% c("counts", "prop", "pres"))) 
    stop("normalize argument can only be \"counts\", \"prop\", or \"pres\" ")
  if (!is(center, "logical")) 
    stop("center argument must be a logical value")
  if (!is.vector(quant)) 
    stop("quant argument must be a vector")
  if (!(distr %in% c("sn", "norm", "t", "cauchy"))) 
    stop("normalize argument can only be \"sn\", \"norm\", \"t\", or \"cauchy\" ")
  if (!is.numeric(alpha)) 
    stop("alpha argument must be a numerical value")
  if (!is(pvals_sim, "NULL") & length(pvals_sim$pvals) == 0) 
    stop("pvals_sim object must be a result from simultaneous PERFect with taxa abundance ordering")
  if (is.null(Order.user)) {
    if (Order == "NP") {
      Order.vec <- NP_Order(X)
    }
    if (Order == "pvals") {
      Order.vec <- pvals_Order(X, pvals_sim)
    }
    if (Order == "NC") {
      Order.vec <- NC_Order(X)
    }
    if (Order == "NCw") {
      Order.vec <- NCw_Order(X)
    }
  }
  else {
    Order.vec <- Order.user
  }
  X <- X[, Order.vec]
  nzero.otu <- apply(X, 2, Matrix::nnzero) != 0
  X <- X[, nzero.otu]
  p <- dim(X)[2]
  Order.vec <- Order.vec[nzero.otu]
  X.orig <- X
  if (normalize == "prop") {
    X <- X/apply(X, 1, sum)
  }
  else if (normalize == "pres") {
    X[X != 0] <- 1
  }
  if (center) {
    X <- apply(X, 2, function(x) {
      x - mean(x)
    })
  }
  Order_Ind <- rep(seq_len(length(Order.vec)))
  DFL <- DiffFiltLoss(X = X, Order_Ind, Plot = TRUE, Taxa_Names = Order.vec)
  Taxa <- Order.vec[-length(Order.vec)]
  lfl <- data.frame(Taxa, log(DFL$DFL))
  names(lfl) <- c("Taxa", "DFL")
  hist <- ggplot(data = lfl, aes(lfl$DFL)) + geom_histogram(bins = nbins, 
                                                            aes(y = ..density..), col = col, fill = fill, alpha = hist_fill) + 
    theme(panel.background = element_rect(fill = "white"), 
          panel.grid.major = element_line(colour = "grey90"), 
          axis.text.x = element_text(size = 10)) + ggtitle("") + 
    xlab("log differences in filtering loss") + ylab("Density")
  if (distr == "norm") {
    if (length(quant) > 2) {
      quant <- quant[(length(quant) - 1):length(quant)]
      print("Warning: more than 2 quantile values are given. \nLargest 2 quantiles are used.")
    }
    if (length(quant) < 2) {
      stop("At least two quantile values must be specified.")
    }
    fit <- fitdistrplus::qmedist(lfl$DFL, distr, probs = quant)
    est <- fit$estimate
    hist <- hist + stat_function(fun = dnorm, args = list(mean = est[1], 
                                                          sd = est[2]), colour = linecol)
    pvals <- pnorm(q = lfl$DFL, mean = est[1], sd = est[2], 
                   lower.tail = FALSE, log.p = FALSE)
  }
  if (distr == "t") {
    if (length(quant) > 2) {
      quant <- quant[(length(quant) - 1):length(quant)]
      print("Warning: more than 2 quantile values are given. \nLargest 2  quantiles are used.")
    }
    if (length(quant) < 2) {
      stop("At least 2 quantile value must be specified.")
    }
    fit <- fitdistrplus::qmedist(lfl$DFL, distr, probs = quant, 
                                 start = list(df = 2, ncp = mean(lfl$DFL)))
    est <- fit$estimate
    hist <- hist + stat_function(fun = dt, args = list(df = est[1], 
                                                       ncp = est[2]), colour = linecol)
    pvals <- pt(q = lfl$DFL, df = est[1], ncp = est[2], lower.tail = FALSE, 
                log.p = FALSE)
  }
  if (distr == "cauchy") {
    if (length(quant) > 2) {
      quant <- quant[(length(quant) - 1):length(quant)]
      print("Warning: more than 2 quantile values are given. \nLargest 2 quantiles are used.")
    }
    if (length(quant) < 2) {
      stop("At least 2 quantile value must be specified.")
    }
    fit <- fitdistrplus::qmedist(lfl$DFL, distr, probs = quant)
    est <- fit$estimate
    hist <- hist + stat_function(fun = dcauchy, args = list(location = est[1], 
                                                            scale = est[2]), colour = linecol)
    pvals <- pcauchy(q = lfl$DFL, location = est[1], scale = est[2], 
                     lower.tail = FALSE, log.p = FALSE)
  }
  if (distr == "sn") {
    if (length(quant) > 3) {
      quant <- quant[(length(quant) - 2):length(quant)]
      print("Warning: more than 3 quantile values are given. \nLargest 3 quantiles are used.")
    }
    if (length(quant) < 3) {
      stop("At least 3 quantile values must be specified.")
    }
    lp <- list(xi = mean(lfl$DFL), omega = sd(lfl$DFL), alpha = 1.5)
    suppressWarnings(fit <- fitdistrplus::qmedist(lfl$DFL, 
                                                  distr, probs = quant, start = lp))
    est <- fit$estimate
    hist <- hist + stat_function(fun = dsn, args = list(xi = est[1], 
                                                        omega = est[2], alpha = est[3]), colour = linecol)
    pvals <- 1 - psn(x = lfl$DFL, xi = est[1], omega = est[2], 
                     alpha = est[3])
  }
  names(pvals) <- names(DFL$DFL)
  if (rollmean) {
    pvals_avg <- zoo::rollmean(pvals, k = 3, align = direction, 
                               fill = NA)
  }
  else {
    pvals_avg <- pvals
  }
  pvals_avg[is.na(pvals_avg)] <- pvals[is.na(pvals_avg)]
  Ind <- which(pvals_avg <= alpha)
  if (length(Ind != 0)) {
    Ind <- min(Ind)
  }
  else {
    Ind <- dim(X)[2] - 1
    warning("no taxa are significant at a specified alpha level")
  }
  filtX <- X.orig[, -seq_len(Ind)]
  return(list(filtX = filtX, info = info, pvals = round(pvals_avg, 
                                                        5), DFL = DFL$DFL, fit = fit, hist = hist, est = est, 
              pDFL = DFL$p + ylab("Difference in Filtering Loss")))
}

pvals_Order <- function (Counts, res_sim) 
{
  vec <- names(res_sim$pvals)
  Order <- c(setdiff(colnames(Counts), vec), vec)
  Order_pvals <- Order[c(1, sort.int(res_sim$pvals, index.return = TRUE, 
                                     decreasing = TRUE)$ix + 1)]
  return(Order_pvals)
}

pvals_Plots <- function (PERFect, X, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha = 0.1) 
{
  if (length(PERFect$pvals) == 0) {
    stop("PERFect object must be a result from PERFect_sim() or PERFect_perm()")
  }
  if (!is(X, "matrix")) {
    X <- as.matrix(X)
  }
  if (!is.vector(quantiles)) 
    stop("quantiles argument must be a vector")
  if (!is.numeric(alpha)) 
    stop("alpha argument must be a numerical value")
  if (!(0 %in% quantiles)) {
    quantiles <- c(0, quantiles)
  }
  if (!(1 %in% quantiles)) {
    quantiles <- c(1, quantiles)
  }
  quantiles <- sort(quantiles)
  pvals <- PERFect$pvals
  Order_pvals <- names(pvals)
  taxa_filt <- Order_pvals[!(Order_pvals %in% colnames(PERFect$filtX))]
  taxa <- max(which(taxa_filt %in% Order_pvals))
  res_FLu <- FiltLoss(X = X[, Order_pvals], Order.user = Order_pvals, 
                      type = "Ind", Plot = TRUE)$FL
  FLu_vals <- res_FLu
  breaks <- quantile(res_FLu, quantiles)
  seq1 <- round(100 * quantiles)[-length(quantiles)]
  seq2 <- round(100 * quantiles)[-1]
  labs <- paste(seq1, seq2, sep = "-")
  labs <- paste("[", labs, ")%", sep = "")
  FLu_vals <- cut(res_FLu, breaks, right = FALSE, labels = labs, 
                  include.lowest = TRUE)
  Ind <- which(names(res_FLu) %in% names(pvals))
  df <- data.frame(seq_len(length(pvals)), pvals, res_FLu[Ind], 
                   FLu_vals[Ind])
  names(df) <- c("Taxa", "p_value", "FLu", "Quantiles")
  p_pvals <- ggplot(df) + geom_point(aes(x = Taxa, y = p_value, 
                                         color = Quantiles)) + ggtitle("Permutation PERFect p-values") + 
    theme(panel.background = element_rect(fill = "white"), 
          panel.grid.major = element_line(colour = "grey90"), 
          axis.text.x = element_text(size = 10, colour = "black", 
                                     angle = 90, hjust = 1)) + guides(color = guide_legend(title = "FLu Quantiles"))
  p_pvals <- p_pvals + geom_hline(yintercept = alpha, color = "red", 
                                  linetype = "dashed")
  p_pvals <- p_pvals + ggtitle("") + geom_vline(xintercept = taxa, 
                                                color = "purple", linetype = "dashed")
  return(list(data = df, plot = p_pvals))
}

TraditR1 <- function (X, thresh = 5) 
{
  if (!is(X, "matrix")) {
    X <- as.matrix(X)
  }
  if (!is.numeric(thresh)) 
    stop("thresh argument must be a numerical value")
  TaxaAll <- colnames(X)
  NNzero <- apply(X, 2, Matrix::nnzero)
  filtX <- X[, NNzero >= thresh]
  return(filtX)
}

TraditR2 <- function (X, Ab_min = 0.001) 
{
  if (!is(X, "matrix")) {
    X <- as.matrix(X)
  }
  if (!is.numeric(Ab_min)) 
    stop("Ab_min argument must be a numerical value")
  if (!(all(apply(X, 1, sum) == 1))) {
    X <- sweep(X, MARGIN = 1, apply(X, 1, sum), "/")
  }
  n <- dim(X)[1]
  Abund <- apply(X, 2, max)
  Nsamp <- apply(X, 2, Matrix::nnzero)
  Psamp <- Nsamp/n
  if (is.null(names(Abund))) {
    names(Abund) <- seq_len(Abund)
  }
  Taxa <- names(Abund[Abund > Ab_min])
  Abund <- Abund[Taxa]
  Nsamp <- Nsamp[Taxa]
  Psamp <- Psamp[Taxa]
  C1 <- names(Abund[Abund > 0.01 & Nsamp >= 1])
  C2 <- names(Abund[Abund > 0.001 & Psamp >= 0.02])
  C3 <- names(Abund[Psamp >= 0.05])
  Taxa <- unique(c(C1, C2, C3))
  filtX <- X[, Taxa]
  return(filtX)
}


NCw_Order <- function (X){
  Counts <- as.matrix(X)
  Netw <- t(Counts) %*% Counts
  Netw2 <- Netw
  diag(Netw2) <- 0
  NC_val <- sort(apply(Netw2, 2, Matrix::nnzero), decreasing = FALSE)
  NC <- names(NC_val)
  n <- dim(Counts)[1]
  n_Tilde <- apply(Counts[, NC], 2, Matrix::nnzero)
  NCW_val <- sort((n_Tilde * NC_val)/n, decreasing = FALSE)
  NCW <- names(NCW_val)
}

NC_Order <- function (X) 
{
  Counts <- as.matrix(X)
  Netw <- t(Counts) %*% Counts
  Netw2 <- Netw
  diag(Netw2) <- 0
  NC_val <- sort(apply(Netw2, 2, Matrix::nnzero), decreasing = FALSE)
  NC <- names(NC_val)
  return(NC)
}

NP_Order <- function (X) 
{
  NP <- names(sort(apply(X, 2, Matrix::nnzero)))
  return(NP)
}

sum_n <- function(tot){
  p <- (-1+sqrt(1+8*tot))/2
  idx <- ceiling(cumsum(c(p:1)))
  return(list("p" = p, "idx" = idx))
}

DiffFiltLoss_j <- function(Perm_Order,Netw, j){
  J_jp1 <-  Perm_Order[seq_len(j+1)]
  #calculate corresponding norm ratios this is the value of the matrix ||
  DFL <-  (2*t(Netw[J_jp1[max(NROW(J_jp1))],-J_jp1])%*%Netw[J_jp1[max(NROW(J_jp1))],-J_jp1]+
             Netw[J_jp1[max(NROW(J_jp1))], J_jp1[max(NROW(J_jp1))]]^2)
  
  return(DFL)
}

Perm_j_s <- function(j, Netw, k,p, p2 = NULL){
  #Netw = X'X - data matrix with taxa in columns and samples in rows
  #k - number of permutations used
  #create a list of k*p arrangements of orders
  if(is.null(p2)){p2 <- p}
  labl <- lapply(seq_len(k),function(x) NULL)
  labl <- lapply(labl,function(x)  sample(seq_len(p),p2))
  FL_j <- vapply(labl,DiffFiltLoss_j, numeric(1), Netw = Netw, j=j)
  return(FL_j)
}

sampl_distr <- function(X, k){
  p <- dim(X)[2]
  Netw <- t(X)%*%X
  # full_norm <- psych::tr(t(Netw)%*%Netw)#this is full norm value
  full_norm <- sum(Netw*Netw)
  #For each taxon j, create a distribution of its DFL's by permuting the labels
  res_all <- lapply(seq_len(p-1),function(x) x)
  
  # Calculate the number of cores
  no_cores <- parallel::detectCores()-1
  # Initiate cluster, start parrallel processing
  cl <- parallel::makeCluster(no_cores)
  #load variables for each core
  parallel::clusterExport(cl,c("DiffFiltLoss_j","Perm_j_s","Netw","k","p"),envir=environment())
  #parallel apply
  FL_j <- parallel::parLapply(cl, res_all, function(x) Perm_j_s(j = x, Netw = Netw, k = k, p = p, p2 = x + 1))
  #FL_j <- lapply(res_all, function(x) Perm_j_s(j = x, Netw =Netw, k=k, p =p, p2 = x+1))
  # End the parallel processing
  parallel::stopCluster(cl)
  
  #divide by the full matrix norm values
  res_pres <- lapply(FL_j, function(x) {x/full_norm})
  return(res_pres)
}

