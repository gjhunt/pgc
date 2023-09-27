nloglik <- function(T_seq, N_seq, n_seq, alpha, beta) {
    tryCatch({
        lls <- -sapply(1:length(T_seq), function(i) suppressWarnings({
            extraDistr::dbetapr(T_seq[i], shape1 = n_seq[i], shape2 = alpha * N_seq[i],
                scale = beta, log = TRUE)
        }))
    }, error = function(e) {
        lls <- Inf
    })
    if (any(!is.finite(lls)))
        lls <- Inf

    return(sum(lls))
}

nloglik_phi <- function(T_seq, N_seq, n_seq, phi, beta) {
    return(nloglik(T_seq, N_seq, n_seq, alpha = phi * beta, beta = beta))
}

#' Estimate poisson parameters phi = alpha/beta and beta from sequence of times (T_seq), number of centers (N_seq), and number of recruits (n_seq), across a collection of clinical trials.
#'
#' Estimation procedure calculates maximum likelihood estimates via coordinate descent on phi and beta. By default each iteration of coordinate descent uses NLOPTR non-linear optimization package using the Constrained Optimization BY Linear Approximation (COBYLA) method. 
#' 
#' @export
#' @param T_seq Sequence of recruitment times across trials.
#' @param N_seq Sequence of number of centers across trials.
#' @param n_seq Sequence of number of recruits across trials.
#' @param phi0 Initial estimate for phi for coordinate descent. If not specified, phi0 is randomly chosen uniformly between phi_bnds.
#' @param beta0 Initial estimate for beta for coordinate descent. If not specified, beta0 is randomly chosen uniformly between beta_bnds.
#' @param phi_bnds Bounds for initalization of phi0 and optimization bounds. If not specified, defaults to 0 to 1. 
#' @param beta_bnds Bounds for initalization of beta0 and optimization bounds. If not specified, defaults to 0 to 1. 
#' @param maxiter Max number of coordinate descent iterations.
#' @param eps_stop Stopping criterion for coordinate descent. Stops when absolute value of difference between either estimate of phi or beta falls below eps_stop. 
#' @param verbose Boolean of whether or not to print out coordinate descent optimization results.
#' @param opts Options for NLOPTR optimzation algorithm. Defaults to list(algorithm = 'NLOPT_LN_COBYLA', xtol_rel = 1e-07, maxeval = 100).
#' @param tritter If verbose=TRUE, how often to print out optimization results. Default is every 100 iterations. 
#' @return A list of output:
#' \itemize{
#' \item{est:} Vector of MLE estimates for phi=alpha/beta and beta. 
#' \item{ll:} value of negative of log-likelihood at MLEs. 
#' }
get_mle <- function(T_seq, N_seq, n_seq, phi0 = NULL, beta0 = NULL, phi_bnds = NULL,
    beta_bnds = NULL, maxiter = 1000, eps_stop = 1e-07, verbose = TRUE, opts = NULL,
    tritter = floor(maxiter/100)) {

    # initial randomized alpha and beta starting values
    if (is.null(beta_bnds))
        beta_bnds <- c(0, 1)

    if (is.null(phi_bnds))
        phi_bnds <- c(0, 1)

    if (is.null(beta0))
        beta0 <- stats::runif(1, beta_bnds[1], beta_bnds[2])

    if (is.null(phi0))
        phi0 <- stats::runif(1, phi_bnds[1], phi_bnds[2])

    phi_estimate <- phi0
    beta_estimate <- beta0

    # initialize optimzation params
    if (is.null(opts))
        opts <- list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-07, maxeval = 100)

    if (maxiter > 0) {
        i <- 0
        tol_not_reached <- TRUE
        while ((i <= maxiter) && tol_not_reached) {
            # coordinate-wise iteration
            phi_old <- phi_estimate
            beta_old <- beta_estimate

            # estimate phi given beta
            loss <- function(x) nloglik_phi(T_seq = T_seq, N_seq = N_seq, n_seq = n_seq,
                phi = x, beta = beta_estimate)
            output <- nloptr::nloptr(x0 = phi_estimate, eval_f = loss, lb = phi_bnds[1],
                ub = phi_bnds[2], opts = opts)
            phi_estimate <- output$solution  #update phi

            # estimate beta given phi
            loss <- function(x) nloglik_phi(T_seq = T_seq, N_seq = N_seq, n_seq = n_seq,
                phi = phi_estimate, beta = x)
            output <- nloptr::nloptr(x0 = beta_estimate, eval_f = loss, lb = beta_bnds[1],
                ub = beta_bnds[2], opts = opts)
            beta_estimate <- output$solution  #update beta

            # print out progress
            if (verbose && (i%%tritter == 0)) {
                tmp_ll <- nloglik(T_seq, N_seq, n_seq, phi_estimate, beta_estimate)
                cat(paste0("iter...", i, ": phi: ", phi_estimate, " beta: ", beta_estimate,
                  " ll: ", tmp_ll, "\n"))
                utils::flush.console()
            }

            i <- i + 1
            tol_not_reached <- (abs(phi_estimate - phi_old) > eps_stop) || (abs(beta_estimate -
                beta_old) > eps_stop)
        }
    }

    ll_final <- nloglik(T_seq, N_seq, n_seq, phi_estimate, beta_estimate)

    return(list(est = c(phi_estimate, beta_estimate), ll = ll_final))
}

#' Estimate poisson parameters phi = alpha/beta and beta from sequence of times (T), number of centers (N), and number of recruits (n), across a collection of clinical trials.
#'
#' Estimation procedure calculates maximum likelihood estimates via coordinate descent on phi and beta with some number of random restarts and chooses the estimates with the highest log-likelihood. 
#' 
#' @export
#' @param T_seq Sequence of recruitment times across trials.
#' @param N_seq Sequence of number of centers across trials.
#' @param n_seq Sequence of number of recruits across trials.
#' @param maxiter Max number of coordinate descent iterations.
#' @param N_restart Number of restarts for initial randomization of coordinate descent variables.
#' @param phi0 Initial estimate for phi for coordinate descent. If not specified, phi0 is randomly chosen uniformly between phi_bnds.
#' @param beta0 Initial estimate for beta for coordinate descent. If not specified, beta0 is randomly chosen uniformly between beta_bnds.
#' @param phi_bnds Bounds for initalization of phi0 and optimization bounds. If not specified, defaults to 0 to 1. 
#' @param beta_bnds Bounds for initalization of beta0 and optimization bounds. If not specified, defaults to 0 to 1. 
#' @param verbose Boolean of whether or not to print out coordinate descent optimization results.
#' @return A list of output:
#' \itemize{
#' \item{phi_hat:} MLE of phi=alpha/beta
#' \item{beta_hat:} MLE of beta
#' \item{ll:} value of negative of log-likelihood at MLEs. 
#' }
estimate <- function(T_seq, N_seq, n_seq, maxiter = 1000, N_restart = 3, phi0 = NULL,
    beta0 = NULL, phi_bnds = NULL, beta_bnds = NULL, verbose = FALSE) {

    opts <- lapply(1:N_restart, function(i) {
        if (verbose) {
            cat("Try : ", i, "\n")
            utils::flush.console()
        }
        get_mle(T_seq = T_seq, N_seq = N_seq, n_seq = n_seq, phi0 = phi0, beta0 = beta0,
            phi_bnds = phi_bnds, beta_bnds = beta_bnds, maxiter = maxiter, verbose = verbose)
    })

    besti <- which.min(sapply(opts, "[[", "ll"))
    bestest <- lapply(opts, "[[", "est")[[besti]]
    phi_estimate <- bestest[1]
    beta_estimate <- bestest[2]

    return(list(phi_hat = phi_estimate, beta_hat = beta_estimate, ll = opts[[besti]]$ll))
}


