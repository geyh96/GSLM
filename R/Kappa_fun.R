#' The function select the RKHS ridge penalty parameter and threshold via stability selection proposed by SUN(2013)
#' 
#' @param X is the covariates matrix.
#' @param Y is the response vector.
#' @param Loss loss function defined: L2, Check,Huber,SVM_Classification, Logistic
#' @param kappa_Tau is the quantile level we choose. tau = 0.5 by default. for other loss, tau is useless.
#' @param lambdas the vector  contains the candiate of the RKHS ridge parameter lambda.
#' @param Kernel only Gaussian Kernel supported now.
#' @param kappa_time the times of the data split and kappa calculated
#' @param kappa_thres the candidates of the thresholding value
#' @param d_kapparatio the alpha_n in the paper SUN (2013) algorithm Step4
#' @param isKernelScale parameter for generating the kernel matrix deployed by the kernlab package. default FALSE and the performance is satisfactory.
#' @return return a list contains the lambda and threshold selected.
#' @export
Kappa_DSL <- function(X,
                      Y,
                      Loss = "L2",
                      kappa_Tau = 0.5,
                      lambdas = seq(0, 2, 0.5),
                      Kernel = "Gaussian",
                      kappa_time = 5,
                      kappa_thres = 10^(seq(-1.7, -1.2, 0.1)),
                      d_kapparatio = 1,
                      isKernelScale = FALSE) {
  N <- length(Y)
  P <- ncol(X)

  MatrixKappa <- matrix(0, length(lambdas), length(kappa_thres))
  rownames(MatrixKappa) <- paste("lambda", round(lambdas, 2), sep = "")
  colnames(MatrixKappa) <- paste("thres", round(kappa_thres, 2), sep = "")

  for (loop_lambda in 1:length(lambdas)) {
    ## For loop for the lambdas
    ilambdas <- lambdas[loop_lambda]
    set.seed(loop_lambda * 100)
    print("KappaLambdaTuning")
    print(loop_lambda)

    repeattime <- 1
    result1 <- matrix(rep(0), kappa_time, P)
    result2 <- matrix(rep(0), kappa_time, P)

    while (repeattime <= kappa_time) {
      index1 <- sample(1:N, N / 2, replace = F)

      x.1 <- X[index1, ]
      y.1 <- Y[index1]
      x.2 <- X[-index1, ]
      y.2 <- Y[-index1]
      ################# Group 1 ###############################

      result1[repeattime, ] <- DSL_base_screening(
                                  x.1,
                                  y.1,
                                  Loss = Loss,
                                  lambda = ilambdas,
                                  Tau = kappa_Tau,
                                  Kernel = Kernel,
                                  isKernelScale = isKernelScale,
                                  print_loss=FALSE)$gradient

      result2[repeattime, ] <- DSL_base_screening(
                                x.2,
                                y.2,
                                Loss = Loss,
                                lambda = ilambdas,
                                Tau = kappa_Tau,
                                Kernel = Kernel,
                                isKernelScale = isKernelScale,
                                print_loss=FALSE)$gradient
      repeattime <- repeattime + 1
    }
    ################### Comparing ############################
    MatrixKappa[loop_lambda, ] <- sapply(
                                      kappa_thres,
                                      Kapp.cal,
                                      kap.re1 = result1,
                                      kap.re2 = result2
                                    )
  }

  MatrixKappa2 <- MatrixKappa
  maxind <- max(which(MatrixKappa2 == max(MatrixKappa2), arr.ind = T)[, 2])
  candi_thres <- max(MatrixKappa2) * d_kapparatio

  thresind <- which(MatrixKappa2 >= candi_thres)[length(which(MatrixKappa2 >= candi_thres))]
  # thresind <- which(MatrixKappa2 >= candi_thres)[1]

  lambdaIndex <- 1
  thresIndex <- thresind

  print("Kappa_DSL is printing")
  print(lambdas[lambdaIndex])
  print(kappa_thres[thresIndex])
  print(MatrixKappa[lambdaIndex, thresIndex])
  print(MatrixKappa)

  resultList <- list(
    MatrixKappa = MatrixKappa,
    Lambda_Selected = lambdas[lambdaIndex],
    Thres_Selected = kappa_thres[thresIndex]
  )
  return(resultList)
}



#' The function select the RKHS ridge penalty parameter and threshold via stability selection proposed by SUN(2013)
#' 
#' @param X is the covariates matrix.
#' @param Y is the response vector.
#' @param Ind vector contains the main effect
#' @param isStrongHeredity boolean indicator
#' @param isInteraction boolean indicator
#' @param Loss loss function defined: L2, Check,Huber,SVM_Classification, Logistic
#' @param kappa_Tau is the quantile level we choose. tau = 0.5 by default. for other loss, tau is useless.
#' @param lambdas the vector  contains the candiate of the RKHS ridge parameter lambda.
#' @param Kernel only Gaussian Kernel supported now.
#' @param kappa_time the times of the data split and kappa calculated
#' @param kappa_thres the candidates of the thresholding value
#' @param d_kapparatio the alpha_n in the paper SUN (2013)
#' @param isKernelScale parameter for generating the kernel matrix deployed by the kernlab package. default FALSE and the performance is satisfactory.
#' @param print_loss a logic value to decide whether print loss or not
#' @return return a list contains the lambda and threshold selected.
#' @export
Kappa_DSL2 <- function(X,
                      Y,
                      Ind,
                      isStrongHeredity = TRUE,                      
                      isInteraction = TRUE,
                      Loss = "L2",
                      kappa_Tau = 0.5,
                      lambdas = seq(0, 2, 0.5),
                      Kernel = "Gaussian",
                      kappa_time = 5,
                      kappa_thres = 10^(seq(-1.7, -1.2, 0.1)),
                      d_kapparatio = 0.95,
                      isKernelScale = FALSE,
                      print_loss=FALSE) {
  # if(0){
  #   Y=y
  #   X=x
  #   kappa_thres=kappa_thres2
  #   lambdas=1e-6
  #   Ind=Ind1
  #   Kernel = "Gaussian"
  #   print_loss= TRUE
  #   Loss="L2"
  # }
  N <- length(Y)
  P <- ncol(X)

  MatrixKappa <- matrix(0, length(lambdas), length(kappa_thres))
  rownames(MatrixKappa) <- paste("lambda", round(lambdas, 2), sep = "")
  colnames(MatrixKappa) <- paste("thres", round(kappa_thres, 2), sep = "")

  ####################Find Dimension
      if(isInteraction){
            temp = DSL_base_order2(
                        X,
                        Y,
                      isStrongHeredity = isStrongHeredity,
                      Loss = "L2" ,
                      lambda = 1,
                      Tau = kappa_Tau,
                      Ind = Ind,
                      Kernel = Kernel,
                      isKernelScale=isKernelScale,
                      print_loss=print_loss
                          )$result_interaction
          d = length(temp)
      }
      if(!isInteraction){
            temp = DSL_base_order2(
                        X,
                        Y,
                      isStrongHeredity = isStrongHeredity,
                      Loss = "L2" ,
                      lambda = 1,
                      Tau = kappa_Tau,
                      Ind = Ind,
                      Kernel = Kernel,
                      isKernelScale=isKernelScale,
                      print_loss=print_loss
                          )$result_Order2
          d = length(temp)
      }
  ##################################


  for (loop_lambda in 1:length(lambdas)) {
    ## For loop for the lambdas
    ilambdas <- lambdas[loop_lambda]
    set.seed(loop_lambda * 100)
    cat("\n")
    print("KappaLambdaTuning")
    print(loop_lambda)

    repeattime <- 1
    result1 <- matrix(rep(0), kappa_time, d)
    result2 <- matrix(rep(0), kappa_time, d)

    while (repeattime <= kappa_time) {
      index1 <- sample(1:N, N / 2, replace = F)

      x.1 <- X[index1, ]
      y.1 <- Y[index1]
      x.2 <- X[-index1, ]
      y.2 <- Y[-index1]
      ################# interaction  ###############################
      if(isInteraction){
            result1[repeattime, ] <- DSL_base_order2(
                                            x.1,
                                            y.1,
                                          isStrongHeredity = isStrongHeredity,
                                          Loss = Loss ,
                                          lambda = ilambdas,
                                          Tau = kappa_Tau,
                                          Ind = Ind,
                                          Kernel = Kernel,
                                          isKernelScale=isKernelScale,
                                          print_loss=print_loss
                                              )$result_interaction

            result2[repeattime, ] <- DSL_base_order2(
                                            x.2,
                                            y.2,
                                          isStrongHeredity = isStrongHeredity,
                                          Loss = Loss ,
                                          lambda = ilambdas,
                                          Tau = kappa_Tau,
                                          Ind = Ind,
                                          Kernel = Kernel,
                                          isKernelScale=isKernelScale,
                                          print_loss=print_loss
                                              )$result_interaction
      
      }

      if(!isInteraction){
      result1[repeattime, ] <- DSL_base_order2(
                                      x.1,
                                      y.1,
                                    isStrongHeredity = isStrongHeredity,
                                    Loss = Loss ,
                                    lambda = ilambdas,
                                    Tau = kappa_Tau,
                                    Ind = Ind,
                                    Kernel = Kernel,
                                    isKernelScale=isKernelScale,
                                    print_loss=print_loss
                                        )$result_order2

      result2[repeattime, ] <- DSL_base_order2(
                                      x.2,
                                      y.2,
                                    isStrongHeredity = isStrongHeredity,
                                    Loss = Loss ,
                                    lambda = ilambdas,
                                    Tau = kappa_Tau,
                                    Ind = Ind,
                                    Kernel = Kernel,
                                    isKernelScale=isKernelScale,
                                    print_loss=print_loss
                                        )$result_order2
      
      }

      repeattime <- repeattime + 1
    }
    ################### Comparing ############################
    MatrixKappa[loop_lambda, ] <- sapply(
                                      kappa_thres,
                                      Kapp.cal,
                                      kap.re1 = result1,
                                      kap.re2 = result2
                                    )
  }

  MatrixKappa2 <- MatrixKappa
  maxind <- max(which(MatrixKappa2 == max(MatrixKappa2), arr.ind = T)[, 2])
  candi_thres <- max(MatrixKappa2) * d_kapparatio
  if(0){
    d_kapparatio = 1
  }
  # thresind <- which(MatrixKappa2 >= candi_thres)[length(which(MatrixKappa2 >= candi_thres))]

  thresind <- which(MatrixKappa2 >= candi_thres, arr.ind = T)[1,2]
  lambdaIndex <- which(MatrixKappa2 >= candi_thres, arr.ind = T)[1,1]
  # lambdaIndex <- 1
  thresIndex <- thresind

  print("Kappa_DSL is printing")
  print(lambdas[lambdaIndex])
  print(kappa_thres[thresIndex])
  print(MatrixKappa[lambdaIndex, thresIndex])
  print(MatrixKappa)

  resultList <- list(
    MatrixKappa = MatrixKappa,
    Lambda_Selected = lambdas[lambdaIndex],
    Thres_Selected = kappa_thres[thresIndex]
  )
  return(resultList)
}
