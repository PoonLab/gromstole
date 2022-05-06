#' Transform a vector to the interior of the feasible region to be used as inits for constrOptim
#' 
#' @param x A vector of numbers
#' 
#' @return A vector with the same length of x such that the sum of the values is larger than 0 but less than 1 and all values are between 0 and 1
to_feasible <- function(x){
    x <- x - min(x, na.rm = TRUE) + 0.0001
    x <- 0.99 * x / sum(x, na.rm = TRUE)
    x
}


#' Initialize a vector of proportions for the variants of concern, prioritizing current (March 2022) most probable VoCs.
#' 
#' @param varmat The variant matrix to be used in the study; only used for the rownames (which should be VOCs in the expected format)
#' 
#' @return varmat a vector with the same length as the number of rows of varmat, such that the values sum to less than one and each value is between 0 and 1
rho_initializer <- function(varmat) {
    rho_init <- rep(10, length = nrow(varmat))
    # Omicron
    probables1 <- c("BA.1$", "BA.2$", "^B.1.1.529$")
    for (i in seq_along(probables1)) {
        var_in_names <- grepl(probables1[i], rownames(varmat))
        if(any(var_in_names)) {
            rho_init[which(var_in_names)] <- 50
        }
    }

    # Delta
    probables2 <- grepl("AY*", rownames(varmat))
    rho_init[which(probables2)] <- 25

    # Ensure it's not on the boundary
    0.99 * rho_init / sum(rho_init)
}


#' Estimate the proportions of VOCs using Constrained OPTimization
#' 
#' @param coco A data frame containing counts, coverage, and mutation names
#' @param varmat The variant matrix to be used in the study. The rownames must be the VoCs and the colnames must be the mutation names (in the same format as the mutation names in `coco`)
#' 
#' @return A list containing the results (\code{$par}) as well as convergence information from \code{constrOptim}.
#' 
#' \describe{
#'      \item{par}{The estimated proportions of the variants of concern (\eqn{rho})}
#'      \item{value}{The estimated optimum value of the objective function. Not often useful.}
#'      \item{counts}{The number of times the function was calculated (i.e. total iterations). Has nothing to do with the counts column in \code{coco}.}
#'      \item{message}{If the algorithm didn't converge, some explanation as to why.}
#'      \item{outer.iterations}{The number of calls to \code{optim}}
#'      \item{barrier.value}{value of the barrier function at the optimum}
#'      \item{init_method}{The method used to initialize \eqn{rho}. Either Prior_Assumption (assumes mostly Omicron, then Delta, then fills out the others to make it sum to 1), "Uniform" for equal values of each variant, or "Nuclear" for a routine that tries 100 random perturbations of the previous values.}
#' }
#' 
#' @export
#' 
#' @details The estimates are found by minimizing the squared difference between the frequency of each mutation and the prediction of a binomial model where the proportion is equal to the sum of rho times the relevant column of varmat and the size parameter is equal to the coverage. 
#' 
#' The algorithm will first try a prior guess based on the current (March 2022) most common VOCs, then will try a uniform proportion, then (the nuclear option) will try 100 random perturbations until it works. Fails gracefully, with list elements indicating the convergence status and the initialization of rho. 
#' 
#' This function currently does not return any estimate of variance for the proportions and should not be trusted beyond a quick check. In practice it produces results equivalent to \code{coda_binom()} (even when it doesn't converge), but always check the residuals.
#' 
#' @seealso \code{\link[stats]{constrOptim}}
#' 
#' @examples
#' varmat = simulate_varmat() # default values (Omicron)
#' coco = simulate_coco(varmat, rel_counts = c(100, 200, 300)) # expect 1/6, 2/6, and 3/6
#' res = copt_binom(coco, varmat)
#' res$par
provoc_optim <- function(coco, varmat) {
    #print(coco)
    muts <- coco$mutation
    cou2 <- coco$count
    cov2 <- coco$coverage
    vari2 <- varmat
    rho_init <- rho_initializer(vari2)

    objective <- function(rho, count, varmat, coverage) {
        prob <- as.numeric(rho %*% varmat)
        prob[coverage == 0] <- 0
        prob[prob == 0 & count != 0] <- 0.000001
        -sum(stats::dbinom(x = count, size = coverage, prob = prob,
            log = TRUE))
    }

    # Constraints will kill me --------------------------------
    # sum(p) < 1 => sum(p) - 1 < 0 => -sum(p) + 1 > 0
    u_sum1 <- rep(-1, length(rho_init))
    c_sum1 <- -1

    # sum(p) > 0
    u_sum0 <- rep(1, length(rho_init))
    c_sum0 <- 0

    # p_i > 0 => 1p_i + 0p_j > 0
    u_p0 <- diag(length(rho_init))
    c_p0 <- rep(0, length(rho_init))

    # p_i < 1 => -p_i > -1 => -p_i + 1 > 0
    u_p1 <- diag(length(rho_init))
    c_p1 <- rep(-1, length(rho_init))

    ui <- rbind(u_sum1, u_sum0, u_p0, u_p1)
    ci <- c(c_sum1, c_sum0, c_p0, c_p1)


    res <- stats::constrOptim(rho_init,
        f = objective, grad = NULL,
        ui = ui, ci = ci,
        count = cou2, coverage = cov2, varmat = vari2,
        control = list(maxit = 100000))
    res$init_method <- "Prior_Assumption"

    if(res$convergence) { # if not converged,
        # Try a different initialization
        rho_init <- 0.99*rep(1/nrow(vari2), nrow(vari2))
        res <- stats::constrOptim(rho_init,
            f = objective, grad = NULL,
            ui = ui, ci = ci,
            count = cou2, coverage = cov2, varmat = vari2,
            control = list(maxit = 100000))
        res$init_method <- "Uniform"
    }

    bestres <- res
    if(res$convergence) { 
        print("Trying the nuclear option for constrOptim.")
        print("This has never actually worked before.")
        print("Godspeed.")
        # Uniform inititialization
        i <- 0
        converged <- FALSE
        while(i < 50 & !converged) {
            i <- i + 1
            if(!i %% 10) print(paste0("Attempt ", i, " of 50."))

            # Add noise to previous iteration
            rho_init <- res$par +
                stats::rnorm(length(rho_init), 0, 0.05)
            # Constrain to interior of feasible region
            rho_init <- to_feasible(rho_init)

            res <- stats::constrOptim(rho_init,
                f = objective, grad = NULL,
                ui = ui, ci = ci,
                count = cou2, coverage = cov2, varmat = vari2,
                control = list(maxit = 10000))
            }
            res$init_method <- "Nuclear"

            # Only take the best results out of all nuclear tries
            if(res$value < bestres$value) {
                bestres <- res
            }

            converged <- !res$convergence
    }


    res_df <- data.frame(rho = bestres$par,
            ci_low = NA,
            ci_high = NA, 
            variant = rownames(varmat))
    convergence <- ifelse(bestres$convergence == 0, TRUE, bestres$convergence)
    convergence_note <- paste("Optim results: ", 
        bestres$convergence, 
        "; Initialization: ", bestres$init_method, sep = '')
    return(list(res_df = res_df, convergence = convergence, convergence_note = convergence_note))
}
