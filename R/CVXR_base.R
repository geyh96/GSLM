# Fit_L2 = function(y,K,lambda, print_loss){
#     N=length(y)
#     y=matrix(y,N,1)
    
#     Kmatrix = K@.Data
#     # print(dim(Kmatrix))
#     # print(N)
    
#     DSL_alpha <- ginv(Kmatrix + N*lambda * diag(rep(1,N)),0) %*% y
#     if(print_loss){
#         print("FittingRSS")
#         print(mean((y - Kmatrix %*% DSL_alpha)^2))
#     }
#     return(DSL_alpha)
# }

#' The L2 loss function with RKHS ridge penalty
#' 
#' @param y is the response variable.
#' @param K is the Kernel Matrix generated from the package kernlab. K@.Data is a matrix contains the data.
#' @param lambda is the tunable parameter of RKHS ridge penalty.
#' @param print_loss a logic value to decide whether print loss or not
CVXR_L2 <- function(y, K,  lambda, print_loss=FALSE) {
    Kmatrix <- K@.Data
    y <- matrix(y, length(y), 1)
    n <- length(y)

    CVXR_beta <- CVXR::Variable(n)
    obj <- (1 / n * sum((y - Kmatrix %*% CVXR_beta)^2) 
            + lambda * CVXR::quad_form(CVXR_beta, Kmatrix))
    problem <- CVXR::Problem(CVXR::Minimize(obj))
    result <- solve(problem,
                    solver = "OSQP", 
                    feastol = 1e-7, 
                    reltol_inacc = 1e-5)
    alpha <- result$getValue(CVXR_beta)

    Rss <- mean((y - Kmatrix %*% alpha)^2)
    alpha_min <- alpha
    Rss_min <- Rss
    #######################################################################################
    if ((any(is.na( alpha))) | (result$status != "optimal")) {
        print("isnaalpha")
        result <- solve(problem, 
                        solver = "ECOS", 
                        feastol = 1e-7, 
                        reltol_inacc = 1e-5)
        alpha <- result$getValue(CVXR_beta)
        Rss <- mean((y - Kmatrix %*% alpha)^2)

        if (Rss_min > Rss) {
            alpha_min <- alpha
            Rss_min <- min(Rss, Rss_min)

        }
    }
    #######################################################################################
    if ((any(is.na( alpha))) | (result$status != "optimal")) {
        print("isnaalpha2")
        result <- solve(problem, 
                        solver = "SCS", 
                        num_iter = 10000)
        alpha <- result$getValue(CVXR_beta)
        Rss <- mean((y - Kmatrix %*% alpha)^2)
        if (Rss_min > Rss) {
            alpha_min <- alpha
            Rss_min <- min(Rss, Rss_min)
        }
    }
    #######################################################################################
    if ((any(is.na( alpha))) | (result$status != "optimal")) {
        print("isnaalpha2")
        result <- solve(problem)
        alpha <- result$getValue(CVXR_beta)
        Rss <- mean((y - Kmatrix %*% alpha)^2)
        if (Rss_min > Rss) {
            alpha_min <- alpha
        }
        Rss_min <- min(Rss, Rss_min)
    }
    
    if(print_loss){
        print("FittingRSS")
        print(mean((y - Kmatrix %*% alpha_min)^2))
    }
    return(alpha_min)
}


#' The Huber loss function with RKHS ridge penalty
#' 
#' @param y is the response variable.
#' @param K is the Kernel Matrix generated from the package kernlab. K@.Data is a matrix contains the data.
#' @param lambda is the tunable parameter of RKHS ridge penalty.
#' @param print_loss a logic value to decide whether print loss or not
CVXR_Huber <- function(y, K, lambda, print_loss=FALSE) {
    Kmatrix <- K@.Data
    n <- length(y)
    y <- matrix(y, n, 1)
    CVXR_beta <- CVXR::Variable(n)

    L2_alpha <- ginv(Kmatrix + n * diag(rep(1,n)),0) %*% y
    huber_M <- 1.345 * 1.483 * median(abs(y - Kmatrix %*% L2_alpha))

    obj <- (1 / (n) * sum(CVXR::huber((y - Kmatrix %*% CVXR_beta), huber_M))
            + lambda* CVXR::quad_form(CVXR_beta, Kmatrix))

    problem <- CVXR::Problem(CVXR::Minimize(obj))
    result <- solve(problem,
                    solver = "OSQP", 
                    feastol = 1e-7, 
                    reltol_inacc = 1e-5)
    alpha <- result$getValue(CVXR_beta)

    Rss <- mean((y - Kmatrix %*% alpha)^2)
    alpha_min <- alpha
    Rss_min <- Rss
    #######################################################################################
    if (any(is.na( alpha)) | (result$status != "optimal")) {
        print("isnaalpha")
        result <- solve(problem, 
                        solver = "ECOS", 
                        feastol = 1e-7, 
                        reltol_inacc = 1e-5)
        alpha <- result$getValue(CVXR_beta)
        Rss <- mean((y - Kmatrix %*% alpha)^2)

        if (Rss_min > Rss) {
            alpha_min <- alpha
            Rss_min <- min(Rss, Rss_min)

        }
    }
    #######################################################################################
    if (any(is.na( alpha)) | (result$status != "optimal")) {
        print("isnaalpha2")
        result <- solve(problem, 
                        solver = "SCS", 
                        num_iter = 10000)
        alpha <- result$getValue(CVXR_beta)
        Rss <- mean((y - Kmatrix %*% alpha)^2)
        if (Rss_min > Rss) {
            alpha_min <- alpha
            Rss_min <- min(Rss, Rss_min)
        }
    }
    #######################################################################################
    if (any(is.na( alpha)) | (result$status != "optimal")) {
        print("isnaalpha2")
        result <- solve(problem)
        alpha <- result$getValue(CVXR_beta)
        Rss <- mean((y - Kmatrix %*% alpha)^2)
        if (Rss_min > Rss) {
            alpha_min <- alpha
        }
        Rss_min <- min(Rss, Rss_min)
    }

    if(print_loss){
        print("FittingRSS")
        print(mean((y - Kmatrix %*% alpha_min)^2))
    }
    return(alpha_min)
}


#' The Quantile loss function with RKHS ridge penalty
#' 
#' @param y is the response variable.
#' @param K is the Kernel Matrix generated from the package kernlab. K@.Data is a matrix contains the data.
#' @param lambda is the tunable parameter of RKHS ridge penalty.
#' @param tau is the quantile level we choose. tau = 0.5 by default.
#' @param print_loss a logic value to decide whether print loss or not
CVXR_Quantile <- function(y, K,  lambda,tau = 0.5, print_loss=FALSE) {
    Kmatrix <- K@.Data
    y <- matrix(y, length(y), 1)
    n <- length(y)

    CVXR_beta <- CVXR::Variable(n)
    obj <- (1 / n * sum(abs(y - Kmatrix %*% CVXR_beta)) 
            + lambda * CVXR::quad_form(CVXR_beta, Kmatrix))
    problem <- CVXR::Problem(CVXR::Minimize(obj))
    result <- solve(problem,
                    solver = "OSQP", 
                    feastol = 1e-7, 
                    reltol_inacc = 1e-5)
    alpha <- result$getValue(CVXR_beta)

    Rss <- mean((y - Kmatrix %*% alpha)^2)
    alpha_min <- alpha
    Rss_min <- Rss
    #######################################################################################
    if (any(is.na( alpha)) | (result$status != "optimal")) {
        print("isnaalpha")
        result <- solve(problem, 
                        solver = "ECOS", 
                        feastol = 1e-7, 
                        reltol_inacc = 1e-5)
        alpha <- result$getValue(CVXR_beta)
        Rss <- mean((y - Kmatrix %*% alpha)^2)

        if (Rss_min > Rss) {
            alpha_min <- alpha
            Rss_min <- min(Rss, Rss_min)

        }
    }
    #######################################################################################
    if (any(is.na( alpha)) | (result$status != "optimal")) {
        print("isnaalpha2")
        result <- solve(problem, 
                        solver = "SCS", 
                        num_iter = 10000)
        alpha <- result$getValue(CVXR_beta)
        Rss <- mean((y - Kmatrix %*% alpha)^2)
        if (Rss_min > Rss) {
            alpha_min <- alpha
            Rss_min <- min(Rss, Rss_min)
        }
    }
    #######################################################################################
    if (any(is.na( alpha)) | (result$status != "optimal")) {
        print("isnaalpha2")
        result <- solve(problem)
        alpha <- result$getValue(CVXR_beta)
        Rss <- mean((y - Kmatrix %*% alpha)^2)
        if (Rss_min > Rss) {
            alpha_min <- alpha
        }
        Rss_min <- min(Rss, Rss_min)
    }
    
    if(print_loss){
        print("FittingRSS")
        print(mean((y - Kmatrix %*% alpha_min)^2))
    }
    return(alpha_min)
}

#' The Logistic loss function with RKHS ridge penalty
#' 
#' @param y is the response variable.
#' @param K is the Kernel Matrix generated from the package kernlab. K@.Data is a matrix contains the data.
#' @param lambda is the tunable parameter of RKHS ridge penalty.
#' @param print_loss a logic value to decide whether print loss or not
CVXR_Logistic <- function(y, K, lambda, print_loss=FALSE) {

    y <- matrix(y, length(y), 1)
    n <- length(y)
    CVXR_beta <- CVXR::Variable(n)
    Kmatrix <- K@.Data

    obj <- 1 / n * (-1 * sum(CVXR::logistic((-1) * Kmatrix[(y <= 0), ] %*% CVXR_beta)) 
                    - sum(CVXR::logistic(Kmatrix[(y == 1), ] %*% CVXR_beta)))
                    - lambda * CVXR::quad_form(CVXR_beta, Kmatrix)
    
    problem <- CVXR::Problem(CVXR::Maximize(obj))
    result <- solve(problem, 
                    solver = "ECOS", 
                    feastol = 1e-7, 
                    reltol_inacc = 1e-5)
    alpha <- result$getValue(CVXR_beta)
    Rss <- result$value
    if(is.na(Rss)){
        Rss=1e4
    }
    alpha_min <- alpha
    Rss_min <- Rss
    ############################################################################

    ##############################################################################
    if (any(is.na( alpha)) | (result$status != "optimal") ) {
        result <- solve(problem, 
                        solver = "SCS", 
                        num_iter = 10000)
        alpha <- result$getValue(CVXR_beta)
        Rss <- result$value
        if(is.na(Rss)){
            Rss=1e4
        }
        if (Rss_min > Rss) {
            alpha_min <- alpha
        }
        Rss_min <- min(Rss, Rss_min)
    }
    #################################################################################
    if (any(is.na( alpha)) | (result$status != "optimal")) {
        result <- solve(problem)
        alpha <- result$getValue(CVXR_beta)
        Rss <- result$value
        if(is.na(Rss)){
            Rss=1e4
        }
        if (Rss_min > Rss) {
            alpha_min <- alpha
        }
        Rss_min <- min(Rss, Rss_min)
    }
    if(print_loss){
        print("FittingRSS")
        print(Rss_min)
    }
    return(alpha_min)
}

#' The SVM-Classification loss function with RKHS ridge penalty
#' 
#' @param y is the response variable.
#' @param K is the Kernel Matrix generated from the package kernlab. K@.Data is a matrix contains the data.
#' @param lambda is the tunable parameter of RKHS ridge penalty.
#' @param print_loss a logic value to decide whether print loss or not
CVXR_SVC <- function(y, K, lambda, print_loss=FALSE) {
    n <- length(y)
    Kmatrix <- K@.Data
    ysvc <- (y == 1) * 1 - (y == 0)
    ysvc <- matrix(ysvc, n, 1)

    CVXR_beta <- CVXR::Variable(n)
    CVXR_v <- CVXR::Variable(1)
    
    obj <- (1 / n * sum(((CVXR::pos(1 - ysvc * (Kmatrix %*% CVXR_beta - CVXR_v))))) 
            + lambda * CVXR::quad_form(CVXR_beta, Kmatrix))

    problem <- CVXR::Problem(CVXR::Minimize(obj))
   result <- solve(problem, 
                    solver = "OSQP", 
                    feastol = 1e-7, 
                    reltol_inacc = 1e-5)
    alpha <- result$getValue(CVXR_beta)
    Rss <- result$value
    alpha_min <- alpha
    Rss_min <- Rss
    ############################################################################
    if (any(is.na( alpha)) | (result$status != "optimal") ) {
        print("isnaalpha")
        result <- solve(problem, 
                        solver = "ECOS", 
                        feastol = 1e-7, 
                        reltol_inacc = 1e-5)
        alpha <- result$getValue(CVXR_beta)
        Rss <- result$value

        if (Rss_min > Rss) {
            alpha_min <- alpha
        }
        Rss_min <- min(Rss, Rss_min)
    }
    ##############################################################################
    if (any(is.na( alpha)) | (result$status != "optimal") ) {
        print("isnaalpha2")
        result <- solve(problem, 
                        solver = "SCS", 
                        num_iter = 10000)
        alpha <- result$getValue(CVXR_beta)
        Rss <- result$value
        if (Rss_min > Rss) {
            alpha_min <- alpha
        }
        Rss_min <- min(Rss, Rss_min)
    }
    #################################################################################
    if (any(is.na( alpha)) | (result$status != "optimal")) {
        print("isnaalpha2")
        result <- solve(problem)
        alpha <- result$getValue(CVXR_beta)
        Rss <- result$value
        if (Rss_min > Rss) {
            alpha_min <- alpha
        }
        Rss_min <- min(Rss, Rss_min)
    }
    if(print_loss){
        print("FittingRSS")
        print(Rss_min)
    }
    return(alpha_min)
}


CVXR_SVR <- function(y, K, lambda, print_loss=FALSE) {
    n <- length(y)
    y <- matrix(y, n, 1)
    Kmatrix <- K@.Data
    CVXR_beta <- CVXR::Variable(n)
    CVXR_v <- CVXR::Variable(1)


    obj <- (1 / n * sum(((CVXR::pos(y - (Kmatrix %*% CVXR_beta - CVXR_v))))) + lambda * CVXR::quad_form(CVXR_beta, Kmatrix))

    problem <- CVXR::Problem(CVXR::Minimize(obj))
    result <- solve(problem, 
                    solver = "OSQP", 
                    feastol = 1e-7, 
                    reltol_inacc = 1e-5)
    alpha <- result$getValue(CVXR_beta)
    Rss <- mean((y - Kmatrix %*% alpha)^2)
    alpha_min <- alpha
    Rss_min <- Rss
    if (any(is.na( alpha)) | (result$status != "optimal")) {
        print("isnaalpha")
        result <- solve(problem, 
                        solver = "ECOS", 
                        feastol = 1e-7, 
                        reltol_inacc = 1e-5)
        alpha <- result$getValue(CVXR_beta)
        Rss <- mean((y - Kmatrix %*% alpha)^2)

        if (Rss_min > Rss) {
            alpha_min <- alpha
        }
        Rss_min <- min(Rss, Rss_min)
    }
    if (any(is.na( alpha)) | (result$status != "optimal")) {
        print("isnaalpha2")
        result <- solve(problem, 
                        solver = "SCS", 
                        num_iter = 10000)
        alpha <- result$getValue(CVXR_beta)
        Rss <- mean((y - Kmatrix %*% alpha)^2)
        if (Rss_min > Rss) {
            alpha_min <- alpha
        }
        Rss_min <- min(Rss, Rss_min)
    }
    if (any(is.na( alpha)) | (result$status != "optimal")) {
        print("isnaalpha2")
        result <- solve(problem)
        alpha <- result$getValue(CVXR_beta)
        Rss <- mean((y - Kmatrix %*% alpha)^2)
        if (Rss_min > Rss) {
            alpha_min <- alpha
        }
        Rss_min <- min(Rss, Rss_min)
    }
    if(print_loss){
        print("FittingRSS")
        print(mean((y - Kmatrix %*% alpha_min)^2))
    }
    return(alpha_min)
}

