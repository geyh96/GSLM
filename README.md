# GSLM
GSLM is a R packages designed for structure learning based on the unstructured kernel based M-regression.
It provide functions to implement the variable selection and interaction selection for various Loss functions,like  L2 loss, Check loss, Huber loss for continuous response and Hinge loss(SVM) and Logistic loss for binary response.


# GSLM
Install the package from github via the $devtools$ package
```
devtools::install_github("geyh96/GSLM")
```
```
##example
run.n=500
run.p=1000
run.rho=0
if(1){
     UW <- matrix(runif(run.n * run.p, -0.5, 0.5), run.n, run.p)
     VW <- matrix(runif(run.n,  -0.5, 0.5), run.n, 1)
     x <- (UW + run.rho * VW \%x\% t(rep(1, run.p))) / (1 + run.rho)
     y0=4*(x[,1] - x[,2] + x[,3] - x[,4] ) - 6*(x[,1]*x[,2] + x[,2]*x[,3] - x[,3]*x[,4] )
     y=y0 + rnorm(run.n, 0, 1)
   }

print("Begin DSL1")
DSL_result1_screening =  GSLM::DSL_screening(
            x, 
            y,
            Loss = "L2", 
            lambdas = 1e-5,  
            Tau = 0.5, 
            Kernel = "Gaussian",
            Kappa_time=10,
            Kappa_thres = 10^seq(-4,0,0.1))
rank(DSL_result1_screening$measure)[1:5]
Ind1 = which(unname(DSL_result1_screening$screening_result)==1)
iselected_DSL1=Ind1
print('which(Ind1==1)')
print(Ind1)   
    # Ind1=c(1,2,3,4,5)
if(length(Ind1)>1){
###########################################begin KappaInteraction
Kappa_result1 = GSLM::Kappa_DSL2(
                x,
                y,
                Ind1,
                isStrongHeredity = TRUE,
                isInteraction = TRUE,
                Loss = "L2",
                kappa_Tau = 0.5,
                lambdas = 10^seq(-5,0,1),
                Kernel = "Gaussian",
                kappa_time = 10 ,
                kappa_thres =10^seq(-4,1.5,0.1),
                d_kapparatio = 1,
                isKernelScale = FALSE)
names(Kappa_result1)
Kappa_result1$Lambda_Selected
Kappa_result1$Thres_Selected
DSL_result1_interaction = GSLM::DSL_base_order2(
        x, 
        y, 
        isStrongHeredity = TRUE,
        Loss = "L2",
        lambda = Kappa_result1$Lambda_Selected,
        Tau = 0.5,
        Ind = Ind1,
        Kernel ="Gaussian",
        isKernelScale=FALSE,
        print_loss=TRUE
            )$gradient_mat_interaction
iselected2_DSL1  = GSLM::get_order2_vec(DSL_result1_interaction>Kappa_result1$Thres_Selected)
}
print(iselected2_DSL1)

```