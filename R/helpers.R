getmode <- function(v) {
  tbl <- table(v)
  if (all(tbl == 1)) {
    median(v)
  } else {
    as.numeric(names(which.max(tbl)))
  }
}

## CVXR 1.8 compatibility: get variable value (works with CVXR 1.0.x and 1.8.x)
.cvxr_get_value <- function(result, v) {
  if (exists("value", mode = "function", envir = asNamespace("CVXR"))) {
    tryCatch(
      as.vector(CVXR::value(v)),
      error = function(e) {
        if (is.function(result$getValue)) result$getValue(v) else stop(e)
      }
    )
  } else if (is.function(result$getValue)) {
    as.vector(result$getValue(v))
  } else {
    stop("Cannot extract variable value from CVXR result")
  }
}

## CVXR 1.8 compatibility: get status (works with CVXR 1.0.x and 1.8.x)
.cvxr_get_status <- function(result, prob) {
  if (!is.null(result$status)) {
    result$status
  } else if (exists("status", mode = "function", envir = asNamespace("CVXR"))) {
    CVXR::status(prob)
  } else {
    stop("Cannot get status from CVXR result")
  }
}

## Treat both "optimal" and "optimal_inaccurate" as success (CVXR 1.8)
.cvxr_is_optimal <- function(status) {
  status %in% c("optimal", "optimal_inaccurate")
}

relevant.funs <- function(intercept = TRUE, model = c("linear", "logistic", "logistic_alter")) {
  model <- match.arg(model)

  if (model == "linear") {
    train.fun <- function(X, y, lambda = NULL) {
      if (is.null(lambda)) lambda <- "CV.min"
      p <- ncol(X)
      htheta <- if (lambda == "CV.min") {
        outLas <- cv.glmnet(X, y,
          family = "gaussian", alpha = 1,
          intercept = intercept, standardize = T
        )
        as.vector(coef(outLas, s = outLas$lambda.min))
      } else if (lambda == "CV") {
        outLas <- cv.glmnet(X, y,
          family = "gaussian", alpha = 1,
          intercept = intercept, standardize = T
        )
        as.vector(coef(outLas, s = outLas$lambda.1se))
      } else {
        outLas <- glmnet(X, y,
          family = "gaussian", alpha = 1,
          intercept = intercept, standardize = T
        )
        as.vector(coef(outLas, s = lambda))
      }
      if (intercept == FALSE) htheta <- htheta[2:(p + 1)]

      return(list(lasso.est = htheta))
    }

    cond_var.fun <- function(pred, y = NULL, sparsity = NULL) {
      n <- length(y)
      sigma.sq <- sum((y - pred)^2) / max(0.7 * n, n - sparsity)
      return(rep(sigma.sq, n))
    }
  } else {
    train.fun <- function(X, y, lambda = NULL) {
      if (is.null(lambda)) lambda <- "CV.min"
      p <- ncol(X)
      htheta <- if (lambda == "CV.min") {
        outLas <- cv.glmnet(X, y,
          family = "binomial", alpha = 1,
          intercept = intercept, standardize = T
        )
        as.vector(coef(outLas, s = outLas$lambda.min))
      } else if (lambda == "CV") {
        outLas <- cv.glmnet(X, y,
          family = "binomial", alpha = 1,
          intercept = intercept, standardize = T
        )
        as.vector(coef(outLas, s = outLas$lambda.1se))
      } else {
        outLas <- glmnet(X, y,
          family = "binomial", alpha = 1,
          intercept = intercept, standardize = T
        )
        as.vector(coef(outLas, s = lambda))
      }
      if (intercept == FALSE) htheta <- htheta[2:(p + 1)]

      return(list(lasso.est = htheta))
    }

    cond_var.fun <- function(pred, y = NULL, sparsity = NULL) {
      cond_var <- pred * (1 - pred)
      return(cond_var)
    }
  }

  if (model == "linear") {
    pred.fun <- function(x) x
    deriv.fun <- function(x) rep(1, length(x))
    weight.fun <- function(x) rep(1, length(x))
  } else if (model == "logistic") {
    pred.fun <- function(x) exp(x) / (1 + exp(x))
    deriv.fun <- function(x) exp(x) / (1 + exp(x))^2
    weight.fun <- function(x) (1 + exp(x))^2 / exp(x)
  } else if (model == "logistic_alter") {
    pred.fun <- function(x) exp(x) / (1 + exp(x))
    deriv.fun <- function(x) exp(x) / (1 + exp(x))^2
    weight.fun <- function(x) rep(1, length(x))
  }

  return(list(
    train.fun = train.fun,
    pred.fun = pred.fun,
    deriv.fun = deriv.fun,
    weight.fun = weight.fun,
    cond_var.fun = cond_var.fun
  ))
}

stars.pval <- function(p.value) {
  unclass(symnum(p.value,
    corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")
  ))
}

Direction_searchtuning <- function(X, loading, weight, deriv, resol = 1.5, maxiter = 10) {
  p <- ncol(X)
  n <- nrow(X)
  mu <- sqrt(2.01 * log(p) / n)
  opt.sol <- rep(0, p + 1)
  loading.norm <- sqrt(sum(loading^2))
  H <- cbind(loading / loading.norm, diag(1, p))

  ## 1st iteration to decide whether increase mu or decrease mu
  iter <- 1
  v <- Variable(p + 1)
  adj.XH <- sqrt(weight) * sqrt(deriv) * (X %*% H)
  obj <- 1 / 4 * sum_squares(adj.XH %*% v) / n + sum((loading / loading.norm) * (H %*% v)) + mu * sum(abs(v))
  # obj = 1/4*sum(((X%*%H%*%v)^2)*weight*deriv)/n + sum((loading/loading.norm)*(H%*%v)) + mu*sum(abs(v))
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  status <- .cvxr_get_status(result, prob)
  if (.cvxr_is_optimal(status)) {
    incr <- -1
    v_opt <- .cvxr_get_value(result, v)
  } else {
    incr <- 1
  }

  ## while loop to find the best mu (the smallest mu satisfying optimal status)
  while (iter <= maxiter) {
    laststatus <- status
    mu <- mu * (resol^incr)
    obj <- 1 / 4 * sum_squares(adj.XH %*% v) / n + sum((loading / loading.norm) * (H %*% v)) + mu * sum(abs(v))
    prob <- Problem(Minimize(obj))
    result <- solve(prob)
    status <- .cvxr_get_status(result, prob)
    if (incr == -1) {
      if (.cvxr_is_optimal(status)) {
        v_opt <- .cvxr_get_value(result, v)
        iter <- iter + 1
        next
      } else {
        step <- iter - 1
        break
      }
    }
    if (incr == 1) {
      if (!.cvxr_is_optimal(status)) {
        iter <- iter + 1
        next
      } else {
        step <- iter
        v_opt <- .cvxr_get_value(result, v)
        break
      }
    }
  }
  if (iter > maxiter) step <- maxiter

  direction <- -(1 / 2) * (v_opt[-1] + v_opt[1] * loading / loading.norm)
  return(list(
    proj = direction,
    step = step,
    incr = incr,
    laststatus = laststatus,
    curstatus = status,
    mu = mu
  ))
}

Direction_fixedtuning <- function(X, loading, weight, deriv, mu = NULL, resol = 1.5, step = 3, incr = -1) {
  p <- ncol(X)
  n <- nrow(X)
  if (is.null(mu)) {
    mu <- sqrt(2.01 * log(p) / n)
    mu <- mu * resol^{
      incr * step
    }
  }
  loading.norm <- sqrt(sum(loading^2))

  H <- cbind(loading / loading.norm, diag(1, p))
  v <- Variable(p + 1)
  adj.XH <- sqrt(weight) * sqrt(deriv) * (X %*% H)
  obj <- 1 / 4 * sum_squares(adj.XH %*% v) / n + sum((loading / loading.norm) * (H %*% v)) + mu * sum(abs(v))
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  opt.sol <- .cvxr_get_value(result, v)
  status <- .cvxr_get_status(result, prob)
  direction <- (-1) / 2 * (opt.sol[-1] + opt.sol[1] * loading / loading.norm)
  return(list(
    proj = direction,
    status = status,
    mu = mu
  ))
}

compute_direction <- function(loading, X, weight, deriv, mu = NULL, verbose = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  loading.norm <- sqrt(sum(loading^2))

  if (loading.norm <= 1e-5) {
    ## loading.norm too small, direction is set as 0
    message("loading norm too small, set proj direction as 0's \n")
    direction <- rep(0, length(loading))
  } else {
    if (n >= 6 * p) {
      ## low-dimensional case
      temp <- sqrt(weight * deriv) * X
      Sigma.hat <- t(temp) %*% temp / n
      direction <- solve(Sigma.hat) %*% loading / loading.norm
    } else {
      ## CVXR sometimes break down accidentally, we catch the error and offer an alternative
      direction <- NULL
      tryCatch(
        expr = {
          if (is.null(mu)) {
            ## when mu is not specified
            if (n >= 0.9 * p) {
              step.vec <- incr.vec <- rep(NA, 3)
              for (t in 1:3) {
                index.sel <- sample(1:n, size = round(0.9 * p), replace = FALSE)
                Direction.Est.temp <- Direction_searchtuning(X[index.sel, , drop = F], loading, weight = weight[index.sel], deriv = deriv[index.sel])
                step.vec[t] <- Direction.Est.temp$step
                incr.vec[t] <- Direction.Est.temp$incr
              }
              step <- getmode(step.vec)
              incr <- getmode(incr.vec)
              Direction.Est <- Direction_fixedtuning(X, loading, weight = weight, deriv = deriv, step = step, incr = incr)
              while (!.cvxr_is_optimal(Direction.Est$status)) {
                step <- step + incr
                Direction.Est <- Direction_fixedtuning(X, loading, weight = weight, deriv = deriv, step = step, incr = incr)
              }
              if (verbose) {
                cat(paste0(
                  "The projection direction is identified at mu = ", round(Direction.Est$mu, 6),
                  "at step =", step, "\n"
                ))
              }
            } else {
              Direction.Est <- Direction_searchtuning(X, loading, weight = weight, deriv = deriv)
              if (verbose) {
                cat(paste0(
                  "The projection direction is identified at mu = ", round(Direction.Est$mu, 6),
                  "at step =", Direction.Est$step, "\n"
                ))
              }
            }
          } else {
            ## when mu is specified
            Direction.Est <- Direction_fixedtuning(X, loading, weight = weight, deriv = deriv, mu = mu)
            while (!.cvxr_is_optimal(Direction.Est$status)) {
              mu <- mu * 1.5
              Direction.Est <- Direction_fixedtuning(X, loading, weight = weight, deriv = deriv, mu = mu)
              if (verbose) cat(paste0("The projection direction is identified at mu = ", round(Direction.Est$mu, 6), "\n"))
            }
          }
          direction <- Direction.Est$proj
        },
        warning = function(w) {
          message("Caught a warning using CVXR!")
          print(w)
        },
        error = function(e) {
          message("Caught an error using CVXR! Alternative method is applied for proj direction.")
          print(e)
        }
      )
      ## Use diagonal fallback when CVXR fails or direction was not set
      if (is.null(direction)) {
        temp <- sqrt(weight * deriv) * X
        Sigma.hat <- t(temp) %*% temp / n
        diag_sigma <- diag(Sigma.hat)
        diag_sigma[diag_sigma < 1e-10] <- 1e-10
        Sigma.hat.inv <- diag(1 / diag_sigma)
        direction <- as.vector(Sigma.hat.inv %*% loading / loading.norm)
      }
    }
  }
  return(direction)
}
