#' The function implement the kappa selection procedure and the variable selection procedure following the paper.
#' 
#' @param X is the covariates matrix.
#' @param Y is the response vector.
#' @param Loss the loss choosed among L2 Check Huber Logistic SVM-Classification
#' @param lambdas is the tunable parameter of RKHS ridge penalty.
#' @param Tau is the quantile level we choose. tau = 0.5 by default. for other loss, tau is useless.
#' @param Kernel only Gaussian Kernel supported now.
#' @param Kappa_time the number of stability selection loop
#' @param Kappa_thres the candidates of threshold selected for variable selection
#' @param Kappa_multiplier the actually used is the selected threshold times Kappa_multiplier. Suggest by the Sun 2013
#' @param Kappa_ratio the alpha suggested in  the step 4 of algorithm from Sun 2013.
#' @param isKernelScale parameter for generating the kernel matrix deployed by the kernlab package. default FALSE and the performance is satisfactory.
#' @param only_gradient a logic value to decide return a list of gradient or more information
#' @examples
#' run.n=500
#' run.p=1000
#' run.rho=0
#'if(1){
#'      UW <- matrix(runif(run.n * run.p, -0.5, 0.5), run.n, run.p)
#'      VW <- matrix(runif(run.n,  -0.5, 0.5), run.n, 1)
#'      x <- (UW + run.rho * VW %x% t(rep(1, run.p))) / (1 + run.rho)
#'      y0=4*(x[,1] - x[,2] + x[,3] - x[,4] ) - 6*(x[,1]*x[,2] + x[,2]*x[,3] - x[,3]*x[,4] )
#'      y=y0 + rnorm(run.n, 0, 1)
#'    }
#'
#' print("Begin DSL1")
#' DSL_result1_screening =  GSLM::DSL_screening(
#'             x, 
#'             y,
#'             Loss = "L2", 
#'             lambdas = 1e-5,  
#'             Tau = 0.5, 
#'             Kernel = "Gaussian",
#'             Kappa_time=10,
#'             Kappa_thres = 10^seq(-4,0,0.1))
#' rank(DSL_result1_screening$measure)[1:5]
#' Ind1 = which(unname(DSL_result1_screening$screening_result)==1)
#' iselected_DSL1=Ind1
#' print('which(Ind1==1)')
#' print(Ind1)   
#'     # Ind1=c(1,2,3,4,5)
#' if(length(Ind1)>1){
#' ###########################################begin KappaInteraction
#' Kappa_result1 = GSLM::Kappa_DSL2(
#'                 x,
#'                 y,
#'                 Ind1,
#'                 isStrongHeredity = TRUE,
#'                 isInteraction = TRUE,
#'                 Loss = "L2",
#'                 kappa_Tau = 0.5,
#'                 lambdas = 10^seq(-5,0,1),
#'                 Kernel = "Gaussian",
#'                 kappa_time = 10 ,
#'                 kappa_thres =10^seq(-4,1.5,0.1),
#'                 d_kapparatio = 1,
#'                 isKernelScale = FALSE)
#' names(Kappa_result1)
#' Kappa_result1$Lambda_Selected
#' Kappa_result1$Thres_Selected
#' DSL_result1_interaction = GSLM::DSL_base_order2(
#'         x, 
#'         y, 
#'         isStrongHeredity = TRUE,
#'         Loss = "L2",
#'         lambda = Kappa_result1$Lambda_Selected,
#'         Tau = 0.5,
#'         Ind = Ind1,
#'         Kernel ="Gaussian",
#'         isKernelScale=FALSE,
#'         print_loss=TRUE
#'             )$gradient_mat_interaction

#' iselected2_DSL1  = GSLM::get_order2_vec(DSL_result1_interaction>Kappa_result1$Thres_Selected)
#' }
#' print(iselected2_DSL1)
#' @export
DSL_screening=function(X, 
                        Y,
                       Loss = "L2", 
                       lambdas = 10^(seq(-3,1,0.3)),  
                       Tau = 0.5, 
                       Kernel = "Gaussian",
                       Kappa_time=10,
                       Kappa_thres=10^(seq(-6,0,0.1)),
                       Kappa_ratio=1,
                        Kappa_multiplier = sqrt(2), 
                       isKernelScale=FALSE,
                       only_gradient=TRUE){
    P=dim(X)[2]
    N=dim(X)[1]
    Kappa_Store <- Kappa_DSL(X = X,
                             Y = Y,
                             Loss = Loss,
                             kappa_Tau =Tau,
                             lambdas =lambdas,
                             Kernel = Kernel,
                             kappa_time = Kappa_time, 
                             kappa_thres = Kappa_thres,
                             d_kapparatio = Kappa_ratio,
                             isKernelScale = isKernelScale)
    print('Kappa_Store$Lambda_Selected')
    print(Kappa_Store$Lambda_Selected)
    print('Kappa_Store$Thres_Selected')
    print(Kappa_Store$Thres_Selected)

    result_screening <- DSL_base_screening(
                    X, 
                    Y,
                    Loss = Loss, 
                    lambda = Kappa_Store$Lambda_Selected,  
                    Tau = Tau, 
                    Kernel = Kernel, 
                    isKernelScale=isKernelScale ,
                    print_loss=TRUE,
                    only_gradient=only_gradient)
    imeasure_DSL = result_screening$gradient
    iselected_DSL = which(imeasure_DSL>Kappa_multiplier * Kappa_Store$Thres_Selected)
    deselected_DSL = which(imeasure_DSL< Kappa_multiplier * Kappa_Store$Thres_Selected)
                    #     gradient = gradient_screening,
                    # Kmatrix = result_alpha$Kmatrix,
                    # sigma = result_alpha$sigma,
    print(iselected_DSL)
    screening_result=rep(0,P)
    names(screening_result)=1:P
    screening_result[iselected_DSL]=1
    
    return(
        list(screening_result = screening_result,
             measure = imeasure_DSL,
             Kmatrix = result_screening$Kmatrix,
             sigma = result_screening$sigma
            )
        )
}

#' The function calculate the parameter of RKHS regression for L2, Quantile, Huber, Hinge and Logistic loss
#' 
#' @param X is the covariates matrix.
#' @param Y is the response vector.
#' @param Loss the loss choosed among L2 Check Huber Logistic SVM-Classification
#' @param lambda is the tunable parameter of RKHS ridge penalty.
#' @param Tau is the quantile level we choose. tau = 0.5 by default. for other loss, tau is useless.
#' @param Kernel only Gaussian Kernel supported now.
#' @param isKernelScale parameter for generating the kernel matrix deployed by the kernlab package. default FALSE and the performance is satisfactory.
#' @param print_loss a logic value to decide whether print loss or not
#' @export
DSL_calculate_alpha <- function(X, 
                    Y, 
                    Loss = "L2",
                    lambda = 1,
                    Tau = 0.5,
                    Kernel = "Gaussian",
                    isKernelScale=FALSE,
                    print_loss=TRUE
                    ){
    # X=x
    N <- dim(X)[1]
    P <- dim(X)[2]
    Xind=1:P
    if (Kernel == "Gaussian") {
        sigma_Kernel <- 1 / sigest(X, scaled = isKernelScale, frac = 1)[2]
        K <- kernelMatrix(rbfdot(sigma = 0.5 / sigma_Kernel), X)
    } else {
        assert(
            "10.15.2020:Now only the Gaussian Kernel is supported", 
            {Kernel == "Gaussian"}
                )
    }
    if (Loss == "L2") {
        DSL_alpha = CVXR_L2(y = Y,
                            K = K,
                            lambda = lambda,
                            print_loss=print_loss)
    }
    if (Loss == "Check") {
        # DSL_alpha <- CVXR_L2(y = Y,
        #                             K = K, 

        #                             quant_Lambda = lambda ,
        #                             print_loss)
        DSL_alpha <- CVXR_Quantile(y = Y,
                            K = K, 
                            tau = Tau, 
                            lambda = lambda ,
                            print_loss)

                }
    if(Loss=="SVM_Classification"){
                DSL_alpha=CVXR_SVC( Y, 
                                    K, 
                                    lambda,
                                    print_loss=print_loss)
                }
    if(Loss == "Logistic"){
                DSL_alpha = CVXR_Logistic(Y, 
                                          K, 
                                          lambda,
                                          print_loss)
                }
    if (Loss == "Huber") {
        DSL_alpha <- CVXR_Huber(y = Y, 
                                K = K, 
                                lambda = lambda ,
                                print_loss)
                        }
    return(list(
            DSL_alpha = DSL_alpha,
            sigma = sigma_Kernel,
            Kmatrix = K@.Data
        ))
}


#' The function calculate the gradient of RKHS regression for L2, Quantile, Huber, Hinge and Logistic loss. The combination of the DSL_calculate_alpha and DSL_calculate_gradient taht is a Cpp function in the Calculate_GradienKeh.cpp
#' 
#' @param X is the covariates matrix.
#' @param Y is the response vector.
#' @param Loss the loss choosed among L2 Check Huber Logistic SVM-Classification
#' @param lambda is the tunable parameter of RKHS ridge penalty.
#' @param Tau is the quantile level we choose. tau = 0.5 by default. for other loss, tau is useless.
#' @param Kernel only Gaussian Kernel supported now.
#' @param isKernelScale parameter for generating the kernel matrix deployed by the kernlab package. default FALSE and the performance is satisfactory.
#' @param print_loss a logic value to decide whether print loss or not
#' @param only_gradient only a list contains gradient returned or a list contains more information returned.
#' @export
DSL_base_screening <- function(
                    X, 
                    Y, 
                    Loss = "L2",
                    lambda = 1,
                    Tau = 0.5,
                    Kernel = "Gaussian",
                    isKernelScale=FALSE,
                    print_loss=TRUE,
                    only_gradient=TRUE){
        cat("Now calculate_alpha")
        result_alpha = DSL_calculate_alpha(
                        X,
                        Y,
                        Loss = Loss,
                        lambda = lambda,
                        Tau = Tau,
                        Kernel = Kernel,
                        isKernelScale = isKernelScale,
                        print_loss = FALSE
                        )
        cat("Now calculate_gradient")
        gradient_screening <- DSL_calculate_gradient(
                                    result_alpha$DSL_alpha,
                                    X,
                                    result_alpha$Kmatrix,
                                    result_alpha$sigma
                                    )
        if(only_gradient){
        
            return(
                list(
                    gradient = gradient_screening
                    ))
        
        }else{
        
            return(
                list(
                    gradient = gradient_screening,
                    Kmatrix = result_alpha$Kmatrix,
                    sigma = result_alpha$sigma,
                    alpha = result_alpha$DSL_alpha
                    ))
        
            }
}



#' Under strong heredity setting, the function calculate the second order gradient of RKHS regression for L2, Quantile, Huber, Hinge and Logistic loss
#' 
#' @param XX is the covariates matrix.
#' @param Y is the response vector.
#' @param Loss loss function defined: L2, Check,Huber,SVM_Classification, Logistic
#' @param lambda is the tunable parameter of RKHS ridge penalty.
#' @param Tau is the quantile level we choose. tau = 0.5 by default. for other loss, tau is useless.
#' @param Ind the vector of the index contains the main effect selected.
#' @param Kernel only Gaussian Kernel supported now.
#' @param isKernelScale parameter for generating the kernel matrix deployed by the kernlab package. default FALSE and the performance is satisfactory.
#' @param print_loss a logic value to decide whether print loss or not
#' @export
DSL_base_strong_heredity <- function(
                    XX, 
                    Y, 
                    Loss = "L2",
                    lambda = 1e-4,
                    Tau = 0.5,
                    Ind = NA,
                    Kernel = "Gaussian",
                    isKernelScale=FALSE,
                    print_loss=TRUE){

        # recalculate the alpha needed
        N = dim(XX)[1]
        P = dim(XX)[2]
        if(any(is.na(Ind))){
            Ind = c(1:P)
            Ind = as.integer(Ind)
            # Ind = 1:5
        }
        
        cat("Now calculate_alpha")
        result_alpha = DSL_calculate_alpha(
                        XX,
                        Y,
                        Loss = Loss,
                        lambda = lambda,
                        Tau = Tau,
                        Kernel = Kernel,
                        isKernelScale = isKernelScale,
                        print_loss = print_loss
                    )
        cat("Now calculate_gradient")
        gradient2 = Calculate_Gradient_order2_Cpp(
                                result_alpha$DSL_alpha,
                                XX, 
                                result_alpha$Kmatrix,
                                result_alpha$sigma,
                                Ind,
                                Ind)
    return(gradient2)
}

#' Under weak heredity setting, the function calculate the second order gradient of RKHS regression for L2, Quantile, Huber, Hinge and Logistic loss
#' 
#' @param X is the covariates matrix.
#' @param Y is the response vector.
#' @param Loss loss function defined: L2, Check,Huber,SVM_Classification, Logistic
#' @param lambda is the tunable parameter of RKHS ridge penalty.
#' @param Tau is the quantile level we choose. tau = 0.5 by default. for other loss, tau is useless.
#' @param Ind the vector of the index contains the main effect selected.
#' @param Kernel only Gaussian Kernel supported now.
#' @param isKernelScale parameter for generating the kernel matrix deployed by the kernlab package. default FALSE and the performance is satisfactory.
#' @param print_loss a logic value to decide whether print loss or not
#' @export
DSL_base_weak_heredity <- function(
                            X, 
                            Y, 
                            Loss = "L2",
                            lambda = 1,
                            Tau = 0.5,
                            Ind = NA,
                            Kernel = "Gaussian",
                            isKernelScale=FALSE,
                            print_loss=TRUE){
        N = dim(X)[1]
        P = dim(X)[2]
        cat("Now calculate_alpha")
        result_alpha = DSL_calculate_alpha(
                        X,
                        Y,
                        Loss = Loss,
                        lambda = lambda,
                        Tau = Tau,
                        Kernel = Kernel,
                        isKernelScale = isKernelScale,
                        print_loss = print_loss
                    )
        cat("Now calculate_gradient")
        if(any(is.na(Ind))){
            Ind = c(1:P)
        }
        Ind1 = as.integer(Ind)
        Ind2 = as.integer((1:P))
        d1 = length(Ind1)
        d2 = length(Ind2)
        cat("Now calculate_gradient")
        gradient2 = Calculate_Gradient_order2_Cpp(
                                result_alpha$DSL_alpha,
                                X,
                                result_alpha$Kmatrix,
                                result_alpha$sigma,
                                Ind1,
                                Ind2)
        # dim(gradient2)
        return(gradient2)
}


#' The function calculate the second order gradient of RKHS regression for L2, Quantile, Huber, Hinge and Logistic loss
#' 
#' @param X is the covariates matrix.
#' @param Y is the response vector.
#' @param isStrongHeredity boolean indicator
#' @param Loss loss function defined: L2, Check,Huber,SVM_Classification, Logistic
#' @param lambda is the tunable parameter of RKHS ridge penalty.
#' @param Tau is the quantile level we choose. tau = 0.5 by default. for other loss, tau is useless.
#' @param Ind the vector of the index contains the main effect selected.
#' @param Kernel only Gaussian Kernel supported now.
#' @param isKernelScale parameter for generating the kernel matrix deployed by the kernlab package. default FALSE and the performance is satisfactory.
#' @param print_loss a logic value to decide whether print loss or not
#' @export
DSL_base_order2=function(
                    X, 
                    Y, 
                    isStrongHeredity = TRUE,
                    Loss = "L2",
                    lambda = 1e-4,
                    Tau = 0.5,
                    Ind = NA,
                    Kernel = "Gaussian",
                    isKernelScale=FALSE,
                    print_loss=TRUE
                        ){

    if(isStrongHeredity){

        assert(
            "Strong Heredity does not permit Ind = NA",
                {   !any(is.na(Ind))  }
                )

        # Ind = 1:10
        XX = X[,Ind]
        # dim(XX)
        d1 = length(Ind)
        gradient_mat = DSL_base_strong_heredity(
                                    XX,
                                    Y,
                                    Loss =  Loss,
                                    lambda = lambda,
                                    Tau = Tau,
                                    Kernel = Kernel,
                                    isKernelScale=isKernelScale,
                                    print_loss=print_loss)
        Ind1 = as.integer(Ind)
        Ind2 = as.integer(Ind)
        rownames(gradient_mat)=paste("X", Ind1, sep="")
        colnames(gradient_mat)=paste("X", Ind2, sep="")
        
    }

    if(!isStrongHeredity){
        gradient_mat = DSL_base_weak_heredity(
                            X, 
                            Y, 
                            Loss = Loss,
                            lambda = lambda,
                            Tau = Tau,
                            Ind = as.integer(Ind),
                            Kernel = Kernel,
                            isKernelScale=isKernelScale,
                            print_loss=print_loss)
        P = dim(X)[2]
        dim(gradient_mat)
        Ind1 = as.integer(Ind)
        Ind2 = as.integer((1:P))
        rownames(gradient_mat)=paste("X", Ind1, sep="")
        colnames(gradient_mat)=paste("X", Ind2, sep="")
    }
        # P=20
    #split the interaction and quadratic term

        d1 = length(Ind1)
        d2 = length(Ind2)
        Boolean_Ind1 = kronecker(matrix(Ind1,d1,1),matrix(rep(1,length(d2)),1,d2))

        Boolean_Ind2 = kronecker(matrix(Ind2,1,d2),matrix(rep(1,length(d2)),d1,1))
        gradient_mat_order2 = gradient_mat * (Boolean_Ind1==Boolean_Ind2)
        gradient_mat_interaction = gradient_mat * (Boolean_Ind1!=Boolean_Ind2)
        result_order2=gradient_mat[Boolean_Ind1==Boolean_Ind2]
        result_interaction=gradient_mat[Boolean_Ind1!=Boolean_Ind2]
        # gradient_mat[1:10,1:10]>1
        return(
            list(
                gradient_mat_order2 = gradient_mat_order2,
                gradient_mat_interaction = gradient_mat_interaction,
                result_order2 = result_order2,
                result_interaction = result_interaction,
                Ind1 = Ind1,
                Ind2=Ind2
            )
        )
}