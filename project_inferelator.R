start_inferalator <- function (e, ...)
{
  ratios <- e$get.cluster.matrix()
  ks <- 1:e$k.clust
  data <- ratios
  
  if (!is.null(envMap)) {
    data <- rbind(ratios, t(as.matrix(envMap)))
    predictors <- c(predictors, colnames(envMap))
  }
  
  if (!is.null(predictors)) 
    predictors <- predictors[predictors %in% rownames(data)]
  
  out <- run_inferelator(ks, data, colMap, predictors, clusterStack = e$clusterStack, gene.prefix = e$genome.info$gene.prefix, ...)
  invisible(out)
}


cv.glmnet <-
  function (x, y, lambda, K = 10, cv.reps = 10, trace = FALSE, 
            plot.it = TRUE, se = TRUE, weights = NA, ...) 
  {
    all.folds <- do.call("c", lapply(1:cv.reps, function(i) cv.folds(length(y), 
                                                                     K)))
    apply.func <- get.apply.func()
    if (is.na(weights[1])) 
      weights <- rep(1, length(y))
    residmat <- do.call(cbind, apply.func(1:length(all.folds), 
                                          function(i) {
                                            omit <- all.folds[[i]]
                                            fit <- my.glmnet(x[-omit, ], y[-omit], lambda = lambda, 
                                                             weights = weights[-omit], ...)
                                            fit <- predict(fit, x[omit, , drop = FALSE], ...)
                                            if (length(omit) == 1) {
                                              fit <- matrix(fit, nrow = 1)
                                            }
                                            apply((y[omit] - fit)^2, 2, mean)
                                          }))
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    object <- list(cv = cv, cv.error = cv.error, lambda = lambda, 
                   fraction = log(lambda))
    if (plot.it) {
      plot(lambda, cv, type = "b", ylim = range(cv, cv + cv.error, 
                                                cv - cv.error), ...)
      if (se) 
        error.bars(lambda, cv + cv.error, cv - cv.error, 
                   width = 1/length(lambda))
    }
    invisible(object)
  }
fill.in.time.series <-
  function (cols, col.map, ...) {
    out.cols <- cols
    all.cols <- rownames(col.map)
    firsts <- which(col.map$is1stLast == "f")
    lasts <- which(col.map$is1stLast == "l")
    #1:length(firsts) = 200; 
    term <- lapply(1:length(firsts), function(i) all.cols[firsts[i]:lasts[i]])
    cols.frac <- length(cols)/length(all.cols)
    for (j in 1:length(term)) {
      frac <- mean(term[[j]] %in% cols) #0
      is.in <- which(term[[j]] %in% out.cols)
      expand <- unique(sort(c(is.in, is.in - 1, is.in + 
                                1)))
      expand <- expand[expand > 0 & expand <= length(term[[j]])]
      out.cols <- unique(c(out.cols, term[[j]][expand]))
    }
    out.cols
  }


get.apply.func <-
  function (plot = F) 
    lapply



get.cluster.predictors <-
  function (cluster.rows, cluster.profile, predictor.mat, predictor.mat.ands,
            group_count, env.names, r_threshold = NA, aic.filter.cutoff = NA, 
            quiet = F, force.include = env.names, ...)
  {
    group_count <- group_count[names(group_count) %in% rownames(predictor.mat)]
    force.include <- force.include[force.include %in% names(group_count)]
    predictor.mat.env <- predictor.mat[force.include, , drop = F]
    possibly.regulates <- numeric()
    
    high.cor <- apply(predictor.mat[names(group_count), ], 1, cor, 
                      cluster.profile, use = "pairwise")

      predictor.mat <- unique(rbind(predictor.mat, predictor.mat.env, 
                                    predictor.mat.ands))

    list(possibly.regulates = possibly.regulates, predictors = rownames(predictor.mat))
  }


get.input.matrix <-
  function (profile, predictor.mat, conds.use = "ALL", col.map = NULL, 
            tau = 10, ratio.cutoff = 3, quiet = F) 
  {
    out.tmp <- profile
    in.tmp <- predictor.mat
    if (!is.null(col.map) && !is.na(tau) && tau > 0) {
      conds <- colnames(predictor.mat)
      cm <- col.map[conds, ]
      good.i <- ((cm$isTs == TRUE) & (cm$is1stLast %in% c("m", 
                                                          "l"))) | (cm$isTs == FALSE & cm$is1stLast == "e")
      curs <- as.character(cm$condName[good.i])
      prevs <- as.character(cm$prevCol[good.i])
      prevs[is.na(prevs)] <- curs[is.na(prevs)]
      del.term <- as.numeric(as.character(cm$delta.t[good.i]))
      del.term[del.term < 1] <- 1
      tmp <- curs %in% names(out.tmp) & prevs %in% names(out.tmp) & 
        prevs %in% colnames(predictor.mat)
      prevs <- prevs[tmp]
      curs <- curs[tmp]
      del.term <- del.term[tmp]
      out.tmp <- ((tau/del.term) * (out.tmp[curs] - out.tmp[prevs])) + 
        out.tmp[prevs]
      in.tmp <- predictor.mat[, prevs]
      colnames(in.tmp) <- names(out.tmp)
    }
    out.tmp[out.tmp > ratio.cutoff] <- ratio.cutoff
    out.tmp[out.tmp < -ratio.cutoff] <- -ratio.cutoff
    in.tmp[in.tmp > ratio.cutoff] <- ratio.cutoff
    in.tmp[in.tmp < -ratio.cutoff] <- -ratio.cutoff
    if (conds.use[1] == "ALL") 
      conds.use <- names(out.tmp)
    df.tmp <- t(in.tmp[, names(out.tmp) %in% conds.use])
    df.tmp <- df.tmp[, !is.na(apply(df.tmp, 2, var, use = "pair")) & 
                       apply(df.tmp, 2, var, use = "pair") > 0.01]
    output <- as.numeric(out.tmp[names(out.tmp) %in% conds.use])
    names(output) <- names(out.tmp)[names(out.tmp) %in% conds.use]
    df.tmp[is.na(df.tmp)] <- 0
    list(inp = df.tmp, outp = output)
  }
get.predictor.matrices <-
  function (predictors, data, gene.prefix = "DVU", preclust.k = NA, 
            funcs = NULL, quiet = F, ...) 
  {
    predictors <- predictors[predictors %in% rownames(data)]
    predictors.genetic <- grep(paste("^", gene.prefix, sep = ""), 
                               predictors, ignore.case = T, value = T)
    predictors.genetic <- predictors.genetic[!grepl("knockout", 
                                                    predictors.genetic, ignore.case = T)]
    env.names <- setdiff(predictors, predictors.genetic)
    if (!is.na(preclust.k) && preclust.k != 0 && preclust.k < 
        length(predictors)) {
      tmp <- preclust.tfs.kmeans(data = data, tfs = predictors.genetic, 
                                 clust.count = preclust.k, ...)
      predictor.mat <- tmp$result
      group_count <- tmp$group_count
      rm(tmp)
    }
    else {
      predictor.mat <- data[predictors, ]
      rownames(predictor.mat) <- predictors
      group_count <- as.list(predictors)
      names(group_count) <- rownames(predictor.mat)
    }
    if (length(env.names) > 0 && any(!env.names %in% rownames(predictor.mat))) 
      predictor.mat <- rbind(predictor.mat, data[env.names[!env.names %in% 
                                                             rownames(predictor.mat)], ])
    tmp <- unique(predictor.mat)
    if (any(!rownames(tmp) %in% rownames(predictor.mat))) 
      cat("Predictor", rownames(tmp)[!rownames(tmp) %in% rownames(predictor.mat)], 
          "is not unique. Removing.\n")
    predictor.mat <- tmp
    rm(tmp)
    predictor.mat.ands <- NULL
    if (!is.na(funcs) && length(funcs) > 0) {
      predictor.mat.ands <- join_predictors(predictor.mat, 
                                                     funcs = funcs, ...)
    }
    list(predictor.mat = predictor.mat, predictor.mat.ands = unique(predictor.mat.ands), 
         genetic.names = unique(predictors.genetic), group_count = group_count, 
         env.names = unique(env.names))
  }

inferelate.one.cluster <-
  function (cluster, predictors, data, col.map = NULL, conds.use = c("clust", 
                                                                     "ALL")[1], quiet = F, plot = T, shrink.opt = c("glmnet", "lars")[1], 
            predictor.mats = NULL, weighted = T, weights = NA, 
            ...) {
    predictor.mats <<- get.predictor.matrices(predictors, 
                                              data, quiet = quiet, ...)
    cluster.rows <- cluster$rows
    if (exists("col.map") && !is.null(col.map)) 
      cluster.conds <- fill.in.time.series(cluster$cols, 
                                           col.map, fill.all.term.frac = 1.25, remove.all.term.frac = 0.1, 
                                           fill.term.gap.size = 1)
    else cluster.conds <- cluster$cols
    cluster.rows <- cluster.rows[cluster.rows %in% rownames(data)]
    cluster.conds <- cluster.conds[cluster.conds %in% colnames(data)]
    cluster.profile <- apply(data[cluster.rows, , drop = F], 
                             2, mean, na.rm = T)
    if (any(is.na(cluster.profile))) 
      cluster.conds <- cluster.conds[-which(is.na(cluster.profile))]
    cluster.weights <- NA
    cluster.vars <- apply(data[cluster.rows, , drop = F], 
                          2, var, na.rm = T)
    cluster.vars <- cluster.vars/(abs(cluster.profile) + 
                                    0.05)          #Co-efficient of variance
    cluster.vars[is.na(cluster.vars) | cluster.vars == 0] <- 1
    cluster.weights <- 1/cluster.vars
    cluster.weights[is.na(cluster.weights) | is.infinite(cluster.weights)] <- 1
    
    #Getting variable list
    tmp <- get.cluster.predictors(cluster.rows, cluster.profile[cluster.conds, 
                                                                drop = F], predictor.mats$predictor.mat[, cluster.conds, 
                                                                                                        drop = F], predictor.mats$predictor.mat.ands[, cluster.conds, 
                                                                                                                                                     drop = F], predictor.mats$group_count, predictor.mats$env.names, 
                                  quiet = quiet, ...)
    possibly.regulates <- tmp$possibly.regulates
    predictor.mat <- rbind(predictor.mats$predictor.mat, predictor.mats$predictor.mat.ands)
    predictor.mat <- predictor.mat[rownames(predictor.mat) %in% 
                                     tmp$predictors, , drop = F]
    predictor.mat <- mean.variance.normalize(predictor.mat, filter = NA)
    
    if (!quiet) 
      cat("Inferelating on biclust #", cluster$k, "using", 
          length(cluster.conds), "conditions and", nrow(predictor.mat), 
          "predictors:\n", paste(rownames(predictor.mat), sep = ", "), 
          "\n")
    
    if (shrink.opt == "glmnet") {
      coeffs <- inferelator.enet(cluster.profile, predictor.mat, 
                                 cluster.conds, col.map = col.map, quiet = quiet, 
                                 weights = cluster.weights, ...)
    }
    else if (shrink.opt == "lars") {
      coeffs <- inferelator(cluster.profile, predictor.mat, 
                            cluster.conds, col.map = col.map, quiet = quiet, 
                            ...)
    }

    singleresult <- coeffs$coeffs
    coeffs.boot <- coeffs$coeffs.boot
    coef.quantiles <- coeffs$coef.quantiles
    all.inputs <- coeffs$all.inputs
    coeffs$coeffs <- coeffs$coeffs.boot <- coeffs$coef.quantiles <- coeffs$all.inputs <- NULL
    coeffs$main <- paste("Bicluster", cluster$k, cluster$nrows, 
                         "genes")
    if (length(singleresult) > 0) {
      conds <- c(cluster.conds, colnames(data)[!colnames(data) %in% 
                                                 cluster.conds])
      coeffs$predictor.mat <- predictor.mat[names(singleresult), 
                                            conds, drop = F]
      coeffs$colors <- c("red", ifelse(singleresult > 0, "#ffaaaa", 
                                       "#aaffaa"))
      return(list(k = cluster$k, coeffs = singleresult, possibly.regulates = possibly.regulates, 
                  cluster.conds = cluster.conds, coeffs.boot = coeffs.boot, 
                  coef.quantiles = coef.quantiles, all.inputs = all.inputs, 
                  plot.info = coeffs))
    }
    else {
      return(list(k = cluster$k, coeffs = numeric(), possibly.regulates = possibly.regulates, 
                  cluster.conds = cluster.conds, coeffs.boot = coeffs.boot, 
                  coef.quantiles = coef.quantiles, all.inputs = all.inputs, 
                  plot.info = coeffs))
    }
  }

inferelator <-
  function (profile, predictor.mat, conds.use, col.map = NULL, 
            tau = 10, ratio.cutoff = 3, coef.cutoff = 0.02, cv.k = 10, 
            min_choice = "min+2se", boot_n = 1, boot.opt = c("resample", 
                                                            "cv"), rescale.coeffs = T, quiet = T, max.coeffs = NA, 
            min.coeffs = NA, ...) 
  {
    if (min_choice == "min") 
      min_choice <- "min+0se"
    tmp <- get.input.matrix(profile, predictor.mat, conds.use, 
                            col.map = col.map, tau = tau, ratio.cutoff = ratio.cutoff, 
                            quiet = quiet)
    df.tmp <- tmp$inp
    output <- tmp$outp
    rm(tmp)
    apply.func <- get.apply.func()
    out.coe <- apply.func(1:boot_n, function(boot) {
      cols <- 1:length(output)
      if (boot > 1 && boot.opt == "resample") 
        cols <- sample(cols, replace = T)
      lars.obj <- try(lars(df.tmp[cols, ], output[cols], type = "lasso", 
                           trace = F), silent = quiet)
      if (class(lars.obj) == "try-error") {
        tries <- 1
        while (tries <= 20 && class(lars.obj) == "try-error") {
          lars.obj <- try(lars(df.tmp[cols, ], output[cols], 
                               type = "lasso", trace = F), silent = quiet)
          tries <- tries + 1
        }
      }
      if (class(lars.obj) == "try-error") 
        return(numeric())
      cv.lars.obj <- try(cv.lars(df.tmp[cols, ], output[cols], 
                                 K = cv.k, type = "lasso", plot.it = F, trace = F), 
                         silent = quiet)
      if (class(cv.lars.obj) == "try-error") {
        tries <- 1
        while (tries <= 20 && class(cv.lars.obj) == "try-error") {
          cv.lars.obj <- try(cv.lars(df.tmp[cols, ], output[cols], 
                                     K = cv.k, type = "lasso", plot.it = F, trace = F), 
                             silent = quiet)
          tries <- tries + 1
        }
      }
      if (class(cv.lars.obj) == "try-error") 
        return(numeric())
      min.i <- which.min(cv.lars.obj$cv)
      min.err <- cv.lars.obj$cv.error[min.i]
      if (grepl("+", min_choice[1], fixed = T)) {
        se <- as.numeric(gsub("min+", "", gsub("se", "", 
                                               min_choice[1])))
        best.s <- min(which(cv.lars.obj$cv <= min(cv.lars.obj$cv) + 
                              se * min.err))
      }
      else best.s <- which.min(cv.lars.obj$cv)
      orig.coeffs <- coeffs <- coef.lars(lars.obj, s = cv.lars.obj$fraction[best.s], 
                                         mode = "fraction")
      sorted <- names(sort(abs(coeffs), decreasing = T))
      coeffs <- coeffs[abs(coeffs) >= coef.cutoff]
      if (!is.na(min.coeffs) && length(coeffs) < min.coeffs) 
        coeffs <- orig.coeffs[sorted[1:min.coeffs]]
      if (!is.na(max.coeffs) && length(coeffs) > max.coeffs) 
        coeffs <- orig.coeffs[sorted[1:max.coeffs]]
      if (!quiet) 
        cat(boot, min_choice[1], min.i, min.err, best.s, cv.lars.obj$cv[best.s], 
            cv.lars.obj$fraction[best.s], length(coeffs), 
            "\n")
      if (rescale.coeffs && length(coeffs) > 0) {
        ins <- df.tmp[, names(coeffs), drop = F]
        coeffs.s <- coef(lm(output[cols] ~ ins[cols, ] - 
                              1))
        names(coeffs.s) <- names(coeffs)
        coeffs <- coeffs.s[abs(coeffs.s) >= coef.cutoff]
      }
      if (boot == 1) {
        out <- list(coeffs = coeffs, lars.obj = lars.obj, 
                    cv.lars.obj = cv.lars.obj, best.s = best.s, se = se, 
                    min.err = min.err)
        return(out)
      }
      coeffs
    })
    lars.obj <- out.coe[[1]]$lars.obj
    cv.lars.obj <- out.coe[[1]]$cv.lars.obj
    best.s <- out.coe[[1]]$best.s
    se <- out.coe[[1]]$se
    min.err <- out.coe[[1]]$min.err
    out.coe[[1]] <- out.coe[[1]]$coeffs
    coeffs <- out.coe[[1]]
    coeffs <- coeffs[order(abs(coeffs), decreasing = T)]
    coef.quantiles <- NULL
    if (boot_n > 1) {
      tmp <- unlist(out.coe)
      tmp2 <- table(names(tmp))
      coef.quantiles <- t(sapply(names(tmp2), function(i) {
        tmp3 <- tmp[names(tmp) == i]
        tmp3 <- c(tmp3, rep(0, boot_n - length(tmp3)))
        c(n = sum(names(tmp) == i)/boot_n, quantile(abs(tmp3), 
                                                    prob = c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95)) * 
            sign(mean(tmp3[tmp3 != 0], na.rm = T)))
      }))
      coef.quantiles <- coef.quantiles[!apply(coef.quantiles, 
                                              1, function(i) all(i[-1] == 0)), ]
    }
    return(list(coeffs = coeffs, coeffs.boot = out.coe, coef.quantiles = coef.quantiles, 
                lars.obj = lars.obj, cv.lars.obj = cv.lars.obj, min_choice = min_choice, 
                best.s = best.s, se = se, min.err = min.err, all.inputs = rownames(df.tmp)))
  }

run_inferelator_enet <- function (profile, predictor.mat, conds.use, col.map = NULL, 
                                  tau = 10, ratio.cutoff = 3, coef.cutoff = 0.02, cv.k = 10, 
                                  min_choice = "min+2se", boot_n = 1, boot.opt = c("resample", 
                                                                                  "cv"), rescale.coeffs = T, quiet = T, alpha = 0.9, weights = NA, 
                                  penalties = NA, max.coeffs = NA, min.coeffs = NA, ...) 
{
  
}

inferelator.enet <-
  function (profile, predictor.mat, conds.use, col.map = NULL, 
            tau = 10, ratio.cutoff = 3, coef.cutoff = 0.02, cv.k = 10, 
            min_choice = "min+2se", boot_n = 1, boot.opt = c("resample", 
                                                            "cv"), rescale.coeffs = T, quiet = T, alpha = 0.9, weights = NA, 
            penalties = NA, max.coeffs = NA, min.coeffs = NA, ...) 
  {
    if (min_choice == "min") {
      min_choice <- "min+0se"
    }
    tmp <- get.input.matrix(profile, predictor.mat, conds.use, 
                            col.map = col.map, tau = tau, ratio.cutoff = ratio.cutoff, 
                            quiet = quiet)
    df.tmp <- tmp$inp[!is.na(tmp$outp), ]
    output <- tmp$outp[!is.na(tmp$outp)]
    rm(tmp)
    orig.penalties <- penalties
    in.penalties <- rep(1, ncol(df.tmp))
    names(in.penalties) <- colnames(df.tmp)
    if (exists("penalties") && !is.na(penalties) && any(names(penalties) %in% 
                                                        names(in.penalties))) {
      penalties <- penalties[names(penalties) %in% names(in.penalties)]
      in.penalties[names(penalties)] <- penalties
    }
    else {
      orig.penalties <- in.penalties
    }
    if (any(orig.penalties != 1)) {
      warning("PENALTIES were set to non-1!!!")
      tmp <- names(which(orig.penalties != 1))
      for (t in tmp) {
        g <- grep(t, names(in.penalties), val = T, fixed = T)
        g <- grep("~~", g, val = T, fixed = T)
        if (length(g) <= 0) 
          next
        gg <- strsplit(g, "~~", fixed = T)
        for (i in 1:length(g)) {
          p1 <- orig.penalties[gg[[i]][1]]
          if (is.na(p1)) 
            p1 <- 1
          p2 <- orig.penalties[gg[[i]][2]]
          if (is.na(p2)) 
            p2 <- 1
          in.penalties[g[i]] <- mean(c(p1, p2), na.rm = T)
        }
      }
      rm(orig.penalties, g, gg, p1, p2, t, tmp)
    }
    if (is.na(weights[1])) 
      weights <- rep(1, length(output))
    else weights <- weights[names(output)]
    names(weights) <- names(output)
    apply.func <- get.apply.func()
    if (!quiet) 
      cat("Alpha =", alpha, "\n")
    out.coe <- apply.func(1:boot_n, function(boot) {
      cols <- 1:length(output)
      if (boot > 1 && boot.opt == "resample") 
        cols <- sample(cols, replace = T)
      glmnet.obj <- my.glmnet(df.tmp[cols, ], output[cols], 
                              penalty.factor = in.penalties, weights = weights[cols], 
                              alpha = if (alpha == "min_choice") 
                                0
                              else alpha, ...)
      if ("try-error" %in% class(glmnet.obj)) {
        tries <- 1
        while (tries <= 20 && "try-error" %in% class(glmnet.obj)) {
          glmnet.obj <- try(my.glmnet(df.tmp[cols, ], output[cols], 
                                      penalty.factor = in.penalties[cols], weights = weights[cols], 
                                      alpha = if (alpha == "min_choice") 
                                        0
                                      else alpha, ...), silent = quiet)
          tries <- tries + 1
        }
      }
      if ("try-error" %in% class(glmnet.obj)) 
        return(numeric())
      cv.glmnet.obj <- try(cv.glmnet(df.tmp[cols, ], output[cols], 
                                     lambda = glmnet.obj$lambda, K = cv.k, trace = F, 
                                     penalty.factor = in.penalties, weights = weights[cols], 
                                     alpha = if (alpha == "min_choice") 
                                       1
                                     else alpha, plot.it = F), silent = quiet)
      if ("try-error" %in% class(cv.glmnet.obj)) {
        tries <- 1
        while (tries <= 20 && "try-error" %in% class(cv.glmnet.obj)) {
          cv.glmnet.obj <- try(cv.glmnet(df.tmp[cols, ], 
                                         output[cols], lambda = glmnet.obj$lambda, K = cv.k, 
                                         trace = F, penalty.factor = in.penalties, weights = weights[cols], 
                                         alpha = if (alpha == "min_choice") 
                                           1
                                         else alpha, plot.it = F), silent = quiet)
          tries <- tries + 1
        }
      }
      if ("try-error" %in% class(cv.glmnet.obj)) 
        return(numeric())
      cv.glmnet.obj$alpha <- alpha
      min.i <- which.min(cv.glmnet.obj$cv)
      min.err <- cv.glmnet.obj$cv.error[min.i]
      se <- 1
      if (grepl("+", min_choice[1], fixed = T)) {
        se <- as.numeric(gsub("min+", "", gsub("se", "", 
                                               min_choice[1])))
        best.s <- min(which(cv.glmnet.obj$cv <= min(cv.glmnet.obj$cv) + 
                              se * min.err))
      }
      else best.s <- which.min(cv.glmnet.obj$cv)
      coeffs.tmp <- as.matrix(coef(glmnet.obj, s = glmnet.obj$lambda[best.s]))
      coeffs <- coeffs.tmp[coeffs.tmp != 0, ]
      names(coeffs) <- rownames(coeffs.tmp)[coeffs.tmp != 0]
      orig.coeffs <- coeffs <- coeffs[names(coeffs) != "(Intercept)"]
      sorted <- names(sort(abs(coeffs), decreasing = T))
      coeffs <- coeffs[abs(coeffs) >= coef.cutoff]
      if (!is.na(min.coeffs) && length(coeffs) < min.coeffs && 
          length(orig.coeffs) >= min.coeffs) 
        coeffs <- orig.coeffs[sorted[1:min.coeffs]]
      if (!is.na(max.coeffs) && length(coeffs) > max.coeffs) 
        coeffs <- orig.coeffs[sorted[1:max.coeffs]]
      if (!quiet) 
        cat(boot, min_choice[1], min.i, min.err, best.s, cv.glmnet.obj$cv[best.s], 
            glmnet.obj$lambda[best.s], length(coeffs), "\n")
      if (rescale.coeffs && length(coeffs) > 1) {
        ins <- df.tmp[, names(coeffs), drop = F]
        glmnet.obj2 <- my.glmnet(ins[cols, , drop = F], output[cols], 
                                 alpha = alpha, penalty.factor = in.penalties, 
                                 weights = weights[cols], ...)
        coeffs.s <- t(as.matrix(coef(glmnet.obj2)))
        coeffs <- coeffs.s[nrow(coeffs.s), ]
        coeffs <- coeffs[abs(coeffs) >= coef.cutoff]
        coeffs <- coeffs[names(coeffs) != "(Intercept)"]
      }
      if (boot == 1) {
        out <- list(coeffs = coeffs, lars.obj = glmnet.obj, 
                    cv.lars.obj = cv.glmnet.obj, best.s = best.s, 
                    se = se, min.err = min.err)
        return(out)
      }
      coeffs
    })
    lars.obj <- out.coe[[1]]$lars.obj
    cv.lars.obj <- out.coe[[1]]$cv.lars.obj
    best.s <- out.coe[[1]]$best.s
    se <- out.coe[[1]]$se
    min.err <- out.coe[[1]]$min.err
    out.coe[[1]] <- out.coe[[1]]$coeffs
    coeffs <- out.coe[[1]]
    coeffs <- coeffs[order(abs(coeffs), decreasing = T)]
    coef.quantiles <- NULL
    if (boot_n > 1) {
      tmp <- unlist(out.coe)
      tmp2 <- sort(table(names(tmp)), decreasing = T)
      coef.quantiles <- t(sapply(names(tmp2), function(i) {
        tmp3 <- tmp[names(tmp) == i]
        tmp3 <- c(tmp3, rep(0, boot_n - length(tmp3)))
        c(n = sum(names(tmp) == i)/boot_n, quantile(abs(tmp3), 
                                                    prob = c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95)) * 
            sign(mean(tmp3[tmp3 != 0], na.rm = T)))
      }))
      coef.quantiles <- coef.quantiles[!apply(coef.quantiles, 
                                              1, function(i) all(i[-1] == 0)), ]

    }
    return(list(coeffs = coeffs, coeffs.boot = out.coe, coef.quantiles = coef.quantiles, 
                lars.obj = lars.obj, cv.lars.obj = cv.lars.obj, min_choice = min_choice, 
                best.s = best.s, se = se, min.err = min.err, all.inputs = rownames(df.tmp)))
  }
load.egrin.data <-
  function (path = ".", ...) 
  {
    load(paste(path, "data_orig_EGRIN/egrin_newcode_workspace.RData", 
               sep = "/"))
    load(paste(path, "data_orig_EGRIN/env_map_egrin.RData", sep = "/"))
    relevant.env <- c("oxygen", "illumination", "Fe", "Cu", "Co", 
                      "Mn", "Zn", "Ni", "gamma", "uv")
    env.map <- env.map[relevant.env, ]
    load(paste(path, "data_orig_EGRIN/col_map_egrin.RData", sep = "/"))
    col.map <- NULL
    for (i in 1:(length(colMap) - 1)) {
      col.map <- rbind(col.map, as.data.frame(colMap[[i]]))
    }
    rownames(col.map) <- names(colMap)[1:(length(colMap) - 1)]
    colnames(col.map)[colnames(col.map) == "del.t"] <- "delta.t"
    pc <- as.character(col.map$prevCol)
    pc[is.na(pc)] <- as.character(col.map$condName[is.na(pc)])
    col.map$prevCol <- as.factor(pc)
    col.map$delta.t[is.na(col.map$delta.t)] <- 9999
    predictors <- c(readLines(paste(path, "data/halo/halo_tfs.txt", 
                                    sep = "/")), rownames(env.map))
    data <- rbind(ratios.egrin, env.map)
    load(paste(path, "data_orig_EGRIN/egrin_coeffs.RData", sep = "/"))
    invisible(list(col.map = col.map, env.map = env.map, predictors = predictors, 
                   data = data, clusterStack.egrin = clusterStack.egrin))
  }
join_predictors <-
  function (predictor.mat, predictors = rownames(predictor.mat), 
            funcs = "min", filter_by_r = 0.7, ...) 
  {
    if (is.null(funcs) || is.na(funcs)) 
      return(predictor.mat)
    result <- NULL
    tmp <- t(combn(predictors, 2))
    tmp <- tmp[tmp[, 1] != tmp[, 2], ]
    apply.func <- get.apply.func()
    for (func in funcs) {
      tmp2 <- do.call(rbind, apply.func(1:nrow(tmp), function(i) apply(predictor.mat[tmp[i, 
                                                                                         ], ], 2, func, na.rm = T)))
      tmp2[is.infinite(tmp2)] <- NA
      rownames(tmp2) <- paste(tmp[, 1], tmp[, 2], rep(func, 
                                                      nrow(tmp)), sep = combine.symbol)
      result <- rbind(result, tmp2)
      rm(tmp2)
    }
    cat("Combined", funcs, "predictor matrix is", nrow(result), 
        "x", ncol(result), "\n")
    out <- result
    if (!is.na(filter_by_r) && filter_by_r > 0 && filter_by_r < 1) {
      apply.func <- get.apply.func()
      all.cors <- cor(t(predictor.mat), t(result), use = "pairwise")
      tmp <- apply.func(1:nrow(result), function(i) {
        nm <- strsplit(rownames(result)[i], combine.symbol, 
                       fixed = T)[[1]]
        tmp.out <- NULL
        ttmp <- all.cors[nm[1:2], i]
        ttmp[is.na(ttmp)] <- 0
        if (!any(ttmp > filter_by_r)) 
          tmp.out <- rownames(result)[i]
        tmp.out
      })
      tmp <- do.call("c", tmp)
      out <- result[tmp, ]
      cat("Filtered for cor <=", filter_by_r, ", combined predictor matrix is now ", 
          nrow(out), "x", ncol(out), "\n")
    }
    attr(out, "filter_by_r") <- filter_by_r
    return(out)
  }
mean.variance.normalize <-
  function (input.matrix, filter = 0.04) 
  {
    if (!is.na(filter)) {
      which.good <- which(apply(input.matrix, 1, function(i) mean(!is.na(i), 
                                                                  na.rm = T)) >= filter)
      input.matrix <- input.matrix[which.good, ]
    }
    means <- apply(input.matrix, 1, mean, na.rm = T)
    sds <- apply(input.matrix, 1, sd, na.rm = T)
    sds[sds == 0 | is.na(sds)] <- 1
    input.matrix <- apply(input.matrix, 2, "-", means)
    input.matrix <- apply(input.matrix, 2, "/", sds)
    return(input.matrix)
  }
my.glmnet <-
  function (x, y, family = c("gaussian", "binomial", "poisson", 
                             "multinomial", "cox")[1], weights, offset = NULL, alpha = 1, 
            nlambda = 100, lambda.min = 1e-06, lambda = NULL, standardize = TRUE, 
            intercept = TRUE, thresh = 1e-04, dfmax = ncol(x) + 1, pmax = min(dfmax * 
                                                                                1.2, ncol(x)), exclude, penalty.factor = rep(1, ncol(x)), 
            lower.limits = -Inf, upper.limits = Inf, maxit = 100, type.gaussian = ifelse(ncol(x) < 
                                                                                           500, "covariance", "naive"), type.logistic = c("Newton", 
                                                                                                                                          "modified.Newton"), standardize.response = FALSE, type.multinomial = c("ungrouped", 
                                                                                                                                                                                                                 "grouped"), ...) 
    glmnet(x, y, family, weights, offset, alpha, nlambda, lambda.min, 
           lambda, standardize, intercept, thresh, dfmax, pmax, exclude, 
           penalty.factor, lower.limits, upper.limits, maxit, type.gaussian, 
           type.logistic, standardize.response, type.multinomial)
plot.cluster.coeffs <-
  function (coefs, scale = 1, cex = 0.5, ...) 
  {
    require(igraph0)
    network <- data.frame()
    comb.cnt <- 1
    node.types <- character()
    for (coe in coefs) {
      if (length(coe$coeffs) <= 0) {
        network <- rbind(network, data.frame(n1 = sprintf("bic%s", 
                                                          coe$k), n2 = sprintf("bic%s", coe$k), weight = NA, 
                                             mode = "-"))
      }
      else {
        for (i in 1:length(coe$coeffs)) {
          n <- strsplit(names(coe$coeffs)[i], combine.symbol, 
                        fixed = T)[[1]]
          if (length(n) == 1) {
            network <- rbind(network, data.frame(n1 = n, 
                                                 n2 = sprintf("bic%s", coe$k), weight = coe$coeffs[i], 
                                                 mode = ">"))
          }
          else {
            n2 <- paste("AND", comb.cnt, sep = "")
            network <- rbind(network, data.frame(n1 = n2, 
                                                 n2 = sprintf("bic%s", coe$k), weight = coe$coeffs[i], 
                                                 mode = ">"))
            network <- rbind(network, data.frame(n1 = n[1], 
                                                 n2 = n2, weight = 0, mode = "-"))
            network <- rbind(network, data.frame(n1 = n[2], 
                                                 n2 = n2, weight = 0, mode = "-"))
            comb.cnt <- comb.cnt + 1
          }
        }
      }
      if (!is.null(coe$possibly.regulates) && length(coe$possibly.regulates) > 
          0) {
        for (i in 1:length(coe$possibly.regulates)) {
          network <- rbind(network, data.frame(n1 = names(coe$possibly.regulates)[i], 
                                               n2 = sprintf("bic%s", coe$k), weight = 0, mode = "*"))
        }
      }
    }
  }
plot_objects <-
  function (coeffs, ...) 
  {
    layout(matrix(c(1, 1, 1, 2, 2, 2, 4, 4, 1, 1, 1, 2, 2, 2, 
                    4, 4, 3, 3, 3, 3, 5, 5, 5, 5, 3, 3, 3, 3, 5, 5, 5, 5, 
                    3, 3, 3, 3, 5, 5, 5, 5), nrow = 5, ncol = 8, byrow = T))
    my.plotCVLars <- function(cv.lars.object, se = TRUE, ...) {
      attach(cv.lars.object)
      plot(fraction, cv, type = "b", ylim = range(cv, cv + 
                                                    cv.error, cv - cv.error), ...)
      if (se) 
        error.bars(fraction, cv + cv.error, cv - cv.error, 
                   width = 1/length(fraction))
      detach(cv.lars.object)
      invisible()
    }
    if (!is.null(coeffs$boot_n)) 
      boot_n <- 1
    else if (!is.null(coeffs[[1]]$boot_n)) {
      boot_n <- coeffs[[1]]$boot_n
      coeb <- coeffs
      coeffs <- coeffs[[1]]
    }
    pi <- coeffs$plot.info
    require(lars)
    require(glmnet)
    if ("glmnet" %in% class(pi$lars.obj)) 
      plot(pi$lars.obj, "lambda")
    else plot(pi$lars.obj)
    lines(rep(pi$cv.lars.obj$fraction[pi$best.s], 2), c(-999, 
                                                        999), col = 2, lty = 2, lwd = 3)
    my.plotCVLars(pi$cv.lars.obj, se = TRUE, main = class(pi$lars.obj)[1])
    if ("glmnet" %in% class(pi$lars.obj)) 
      legend("bottomleft", pi$min_choice)
    else legend("topright", pi$min_choice)
    lines(rep(pi$cv.lars.obj$fraction[pi$best.s], 2), c(-999, 
                                                        999), col = 2, lty = 2, lwd = 3)
    if (grepl("+", pi$min_choice, fixed = T)) {
      lines(rep(pi$cv.lars.obj$fraction[which.min(pi$cv.lars.obj$cv)], 
                2), rep(min(pi$cv.lars.obj$cv), 2) + c(0, pi$se * 
                                                         pi$min.err), col = 2, lty = 2, lwd = 1)
      lines(c(pi$cv.lars.obj$fraction[pi$best.s], pi$cv.lars.obj$fraction[which.min(pi$cv.lars.obj$cv)]), 
            rep(min(pi$cv.lars.obj$cv), 2) + pi$se * pi$min.err, 
            col = 2, lty = 2, lwd = 1)
    }
    conds <- pi$clust.conds.plot
    if (length(coeffs$coeffs) > 0) {
      matplot(t(rbind(coeffs$observed[conds], pi$predictor.mat[, 
                                                               conds])), col = pi$colors, ylab = "Normalized expression", 
              xlab = "Conditions", type = "l", main = pi$main)
      legend("bottomright", c("biclust", names(coeffs$coeffs)), 
             col = pi$colors, lty = 1, cex = 0.5)
      lines(pi$cluster.profile, col = "red")
    }
    else {
      plot(coeffs$observed[conds], col = "red", ylab = "Normalized expression", 
           xlab = "Conditions", type = "l", main = pi$main)
      legend("bottomright", "biclust", col = "red", lty = 1, 
             cex = 0.5)
    }
    if (!is.null(coeffs$pred.term) && nrow(coeffs$pred.term) > 1) {
      matlines(t(apply(coeffs$pred.ss[, pi$clust.conds.plot], 
                       2, quantile, prob = c(0.1, 0.9), na.rm = T)), col = rep("lightblue", 
                                                                               2), lty = 1, lwd = 3)
      matlines(t(apply(coeffs$pred.term[, pi$clust.conds.plot], 
                       2, quantile, prob = c(0.1, 0.9), na.rm = T)), col = rep("gray", 
                                                                               2), lty = 1, lwd = 3)
    }
    if (boot_n > 1) {
      coeb <- coeb[sapply(coeb, length) > 0]
      pred.term <- t(sapply(coeb, "[[", "pred.term"))
      tmp <- t(apply(pred.term, 2, quantile, prob = c(0.05, 0.5, 
                                                    0.95)))
      rownames(tmp) <- colnames(coeb[[1]]$pred.term)
      matlines(tmp[pi$clust.conds.plot, ], typ = "l", lty = 1, 
               col = c("gray", "red", "gray"), lwd = 3)
    }
    lines(coeffs$pred.ss[1, pi$clust.conds.plot], col = "blue")
    lines(coeffs$pred.term[1, pi$clust.conds.plot], col = "black")
    lines(rep(pi$n.conds, 2), c(-999, 999), col = "gray", lty = 2, 
          lwd = 3)
    legend("bottomleft", c("pred.ss", "pred.term"), col = c("blue", 
                                                          "black"), lty = 1, cex = 0.5)
    lines(coeffs$observed[pi$clust.conds.plot], col = "red")
    out.net <- try(plot.cluster.coeffs(list(coeffs)))
    if (class(out.net) == "try-error") 
      plot(1:10)
    invisible(out.net)
  }
plot.coeff.stats <-
  function (coeffs) 
  {
    par(mfrow = c(2, 2))
    hist(sapply(coeffs, function(i) length(i$coeffs)), breaks = 20, 
         xlab = "Number of coeffs per bicluster")
    hist(unlist(lapply(coeffs, function(i) i$coeffs)), breaks = 20, 
         xlab = "Coefficient values")
    legend("topleft", as.character(sum(unlist(lapply(coeffs, 
                                                     function(i) i$coeffs < 0)))), text.col = "green", cex = 0.7)
    legend("topright", as.character(sum(unlist(lapply(coeffs, 
                                                      function(i) i$coeffs > 0)))), text.col = "red", cex = 0.7)
    hist(sapply(coeffs, function(i) i$rmsd["term"]), breaks = 20, 
         xlim = c(0, 1), xlab = "RMSD, In")
    hist(sapply(coeffs, function(i) i$rmsd["term.out"]), breaks = 20, 
         xlim = c(0, 1), xlab = "RMDS, Out")
  }
preclust.tfs.kmeans <-
  function (data, tfs, clust.count, n.iter = 200, n.start = 25, 
            seed = 31337, r_threshold = 0.85, ...) 
  {
    if (!is.na(seed)) 
      set.seed(seed)
    tf.matrix <- data[tfs[tfs %in% rownames(data)], ]
    data.c <- kmeans(tf.matrix, clust.count, iter.max = n.iter, 
                     nstart = n.start)
    result <- data.c$centers
    rownames(result) <- paste("TFGROUP", 1:clust.count, sep = "")
    group_count <- lapply(1:length(data.c$size), function(i) names(which(data.c$cluster == 
                                                                         i)))
    names(group_count) <- rownames(result)
    cat("Preclustered with k-means, predictor matrix is", nrow(result), 
        "x", ncol(result), "\n")
    tmp <- apply(result, 1, function(i) apply(data[tfs, ], 1, 
                                              cor, i))
    if (any(tmp > r_threshold)) {
      high.cors <- apply(tmp, 2, function(i) which(i > r_threshold))
      to.be.added <- lapply(names(group_count), function(i) names(high.cors[[i]])[!names(high.cors[[i]]) %in% 
                                                                                  group_count[[i]]])
      for (i in 1:length(group_count)) group_count[[i]] <- unique(c(group_count[[i]], 
                                                                to.be.added[[i]]))
      result <- t(sapply(group_count, function(i) apply(data[i, 
                                                           , drop = F], 2, mean, na.rm = T)))
    }
    for (i in 1:length(group_count)) if (length(group_count[[i]]) == 
                                       1) 
      names(group_count)[i] <- rownames(result)[i] <- group_count[[i]][1]
    return(list(result = result, group_count = group_count))
  }
predictelate <-
  function (cluster.rows, coeffs, ratios, predictor.mats = NULL, 
            group_count = NULL, col.map = NULL, tau = 10, max.coeffs = length(coeffs), 
            ...) 
  {
    if (length(coeffs) <= 0) {
      out <- ratios[1, ] * 0
      return(out)
    }
    if (max.coeffs < length(coeffs)) 
      coeffs <- sort(coeffs, decreasing = T)[1:max.coeffs]
    coeff.names <- unique(unlist(strsplit(names(coeffs), combine.symbol, 
                                          fixed = T)))
    if (is.null(predictor.mats) && !is.null(group_count) && any(coeff.names %in% 
                                                              names(group_count))) {
      tfgroup.ratios <- t(sapply(group_count[which(names(group_count) %in% 
                                                   coeff.names)], function(i) apply(ratios[i, , drop = F], 
                                                                                    2, mean)))
      rownames(tfgroup.ratios) <- names(group_count)[names(group_count) %in% 
                                                     coeff.names]
      ratios <- rbind(ratios, tfgroup.ratios)
    }
    else {
      ratios <- rbind(ratios, predictor.mats$predictor.mat[, 
                                                           colnames(ratios)], predictor.mats$predictor.mat.ands[, 
                                                                                                                colnames(ratios)])
    }
    out.ss <- 0
    for (j in 1:length(coeffs)) {
      if (coeffs[j] == 0) 
        next
      nodes <- unlist(strsplit(names(coeffs)[j], combine.symbol, 
                               fixed = T))
      if (length(nodes) == 1) {
        if (!nodes[1] %in% rownames(ratios)) 
          next
        tmp <- ratios[nodes, ]
      }
      else if (length(nodes) == 3) {
        if (names(coeffs)[j] %in% rownames(ratios)) {
          tmp <- ratios[names(coeffs)[j], ]
        }
        else {
          if (!all(nodes[1:2] %in% rownames(ratios))) 
            next
          tmp <- apply(ratios[c(nodes[1], nodes[2]), ], 
                       2, FUN = nodes[3], na.rm = T)
        }
      }
      tmp[is.na(tmp)] <- 0
      out.ss <- out.ss + tmp * coeffs[j]
    }
    out <- out.ss
    if (!is.null(col.map) && !is.na(tau) && tau > 0) {
      conds <- colnames(ratios)
      prevs <- as.character(col.map[conds, "prevCol"])
      del.term <- as.numeric(as.character(col.map[conds, "delta.t"]))
      del.term[del.term < 1] <- 1
      cluster.prof <- apply(ratios[cluster.rows, , drop = F], 
                            2, mean, na.rm = T)
      tmp1 <- cluster.prof[prevs]
      tmp1[is.na(tmp1)] <- 0
      tmp2 <- out.ss[conds]
      tmp2[is.na(tmp2)] <- 0
      out.term <- (tau * tmp1 + del.term * tmp2)/(tau + del.term)
      names(out.term) <- conds
      out <- out.term
    }
    out[out > 3] <- 3
    out[out < -3] <- -3
    out
  }

run_inferelator <-
  function (ks, data, col.map, predictors, clusterStack, tau = 10, 
            plot = T, coeffs = NULL, group_count = Inf, boot_n = 1, 
            boot.opt = c("resample.lars", "resample.rows", "resample", "lars")[1], ...) 
  {
    in.args <- c(mget(names(formals()), env = as.environment(-1)), 
                 sapply(as.list(substitute({
                   ...
                 })[-1]), deparse))
    data <- mean.variance.normalize(data, filter = 0.04)
    predictors <- predictors[predictors %in% rownames(data)]
    apply.func <- get.apply.func(plot)
    
    out <- apply.func(ks, function(i) {
      cluster <- clusterStack[[i]]
      k <- cluster$k
      apply.func <- get.apply.func()
      apply.func <- lapply
      out.k <- apply.func(1:boot_n, function(boot) {
        clust <- cluster
        if (length(clust$cols) <= 2) 
          return(NULL)
        coeffs <- inferelate.one.cluster(clust, predictors, 
                                         data, predictor.mats = predictor.mats, tau = tau, 
                                         col.map = col.map, boot_n = 1, boot.opt = "resample", 
                                         quiet = boot_n > 1, ...)
        clust.rows <- clust$rows[clust$rows %in% rownames(data)]
        clust.conds <- sort(coeffs$cluster.conds)
        clust.conds <- clust.conds[clust.conds %in% colnames(data)]
        observed <- apply(data[clust.rows, , drop = F], 2, 
                          mean, na.rm = T)
        apply.func <- get.apply.func()
        pred.ss <- do.call(rbind, apply.func(coeffs$coeffs.boot, 
                                             function(b) predictelate(clust.rows, b, data, 
                                                                      predictor.mats = predictor.mats, tau = tau, 
                                                                      ...)))
        pred.term <- do.call(rbind, apply.func(coeffs$coeffs.boot, 
                                             function(b) predictelate(clust.rows, b, data, 
                                                                      predictor.mats = predictor.mats, tau = tau, 
                                                                      col.map = col.map, ...)))
        if (is.null(pred.ss)) 
          pred.ss <- t(observed * 0)
        if (is.null(pred.term)) 
          pred.term <- t(observed * 0)
        vars <- apply(data[clust.rows, , drop = F], 2, 
                      var, na.rm = T)
        vars <- vars/(abs(observed) + 0.05)
        vars[is.na(vars) | vars == 0] <- 1
        weights <- 1/vars
        weights <- weights/sum(weights) * length(weights)
        rmsd.ss <- sqrt(weighted.mean((pred.ss[nrow(pred.ss), 
                                               ] - observed)[clust.conds]^2, weights[clust.conds], 
                                      na.rm = T))
        rmsd.term <- sqrt(weighted.mean((pred.term[nrow(pred.term), 
                                               ] - observed)[clust.conds]^2, weights[clust.conds], 
                                      na.rm = T))
        not.clust.conds <- colnames(data)[!colnames(data) %in% 
                                            clust.conds]
        rmsd.term.out <- sqrt(weighted.mean((pred.term[nrow(pred.term), 
                                                   ] - observed)[not.clust.conds]^2, weights[not.clust.conds], 
                                          na.rm = T))
        coeffs$plot.info$main <- paste("Bicluster", cluster$k, 
                                       cluster$nrows, "genes")
        coeffs$plot.info$clust.conds.plot <- c(clust.conds, 
                                               sort(colnames(data)[!colnames(data) %in% clust.conds]))
        coeffs$plot.info$n.conds <- length(clust.conds)
        #cat(k, tau, rmsd.ss, rmsd.term, rmsd.term.out, "\n")
        coeffs$pred.ss <- pred.ss
        coeffs$pred.term <- pred.term
        coeffs$rmsd <- c(ss = rmsd.ss, term = rmsd.term, term.out = rmsd.term.out)
        coeffs$observed <- observed
        coeffs$boot_n <- boot_n
        coeffs$boot.opt <- boot.opt
        attr(coeffs, "class") <- "coeff.obj"
        if (boot > 1) 
          coeffs$plot.info <- NULL
        coeffs
      })
      names(out.k) <- paste(k, 1:boot_n, sep = ".")
      try(plot_objects(out.k, ...))
      
      attr(out.k, "class") <- "coeff.obj"
      out.k
    })
    out <- do.call("c", out)
    attr(out, "CALL") <- match.call(expand.dots = T)
    attr(out, "class") <- "coeff.obj"
    invisible(out)
  }


