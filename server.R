#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(plotly)
library(DT)
library(patchwork)
library(latex2exp)
library(RColorBrewer)
library(htmltools)


# Define server logic required to draw power plot
shinyServer(function(input, output, session) {
  options(shiny.sanitize.errors = FALSE)
  #### FUNCTIONS ####
  var_hte <- function(m, oicc, cicc, var_y, var_x, var_w){
    top <- var_y * (1 - oicc) * ( 1 + (m - 1) * oicc )
    bottom <- m * var_w * var_x * ( 1 + (m - 2) * oicc - (m - 1) * cicc * oicc )
    
    var_result <- top/bottom
    return(var_result)
  }
  
  power_hte <- function(oicc, cicc, m, n, var_y, var_x, var_w=0.25, d, a=0.05){
    z <- qnorm( 1 - a/2 )
    num_sqrt <- d^2 * var_w * var_x * ( 1 + (m - 2) * oicc - (m - 1) * cicc * oicc )
    denom_sqrt <- var_y * (1 - oicc) * ( 1 + (m - 1) * oicc )
    
    # dividing effect size d by var_hte and subtracting critical value #
    inside <- sqrt( n * m * num_sqrt/denom_sqrt ) - z
    
    power_result <- pnorm(inside)
    return(power_result)
  }
  
  n_hte <- function(oicc, cicc, m, power=0.9, var_y, var_x, var_w=0.25, d, a=0.05){
    
    s2_hte <- var_hte(m, oicc, cicc, var_y, var_x, var_w)
    
    z_typeI <- qnorm(1-a/2)
    z_power <- qnorm(power)
    
    n_result <- ( s2_hte * (z_typeI + z_power)^2 )/d^2
    
    z <- qnorm( 1 - a/2 )
    num_sqrt <- d^2 * var_w * var_x * ( 1 + (m - 2) * oicc - (m - 1) * cicc * oicc )
    denom_sqrt <- var_y * (1 - oicc) * ( 1 + (m - 1) * oicc )
    
    # dividing effect size d by var_hte and subtracting critical value #
    inside <- sqrt( ceiling(n_result) * m * num_sqrt/denom_sqrt ) - z
    
    power_result <- pnorm(inside)
    
    return( list(n=ceiling(n_result), power_emp = power_result ) )
  }
  
  power_ate <- function(oicc, m, n, var_y, var_w=0.25, d, a=0.05){
    z <- qnorm( 1 - a/2 )
    num_sqrt <- d^2 * var_w
    denom_sqrt <- var_y * (1 - oicc) * ( 1 + (m - 1) * oicc )
    
    inside <- sqrt( n * m * num_sqrt/denom_sqrt ) - z
    
    power_result <- pnorm(inside)
    return(power_result)
  }
  
  var_hte_three <- function(ns, m, pw, var_y, var_x, a0, a1, r0, r1, rand){
    
    # eigen values for outcome #
    l1 <- 1 - a0
    l2 <- 1 + (m-1)*a0 - m*a1
    l3 <- 1 + (m-1)*a0 + (ns-1)*m*a1
    
    # eigen values for effect modifier #
    e1 <- 1 - r0
    e2 <- 1 + (m-1)*r0 - m*r1
    e3 <- 1 + (m-1)*r0 + (ns-1)*m*r1
    
    # randomization variance #
    var_w <- pw*(1-pw)
    
    if(rand=="cluster"){
      var_4 <- var_y/(var_w * var_x) * (ns*m)/( (1/l3)*e3 + (ns-1)*(1/l2)*e2 + ns*(m-1)*(1/l1)*e1 )
    }else if(rand=="subcluster"){
      var_4 <- var_y/(var_w * var_x) * (m)/( m*(1/l1) - (1 + (m-1)*r0)*((1/l1) - (1/l2)) )
    }else if(rand=="individ"){
      var_4 <- var_y/(var_w*var_x) * l1
    }
    
    if(any(c(l1,l2,l3) < 0)){
      stop("Variance is not positive definite - please enter new outcome correlation parameters.")
    }else if(any(c(e1,e2,e3)<0)) {
      stop("Variance is not positive definite - please enter new covariate correlation parameters.")
    }else{
      return(var_4)#list(var=var_hte, tau=tau_w))  
    }
    
  }
  
  power_hte_three <- function(nc, ns, m, pw, var_y, var_x, a0, a1_a0, r0, r1_r0, rand, d, a=0.05){
    z <- qnorm(1 - a/2)
    a1 <- a1_a0*a0
    r1 <- r1_r0*r0
    
    var_4 <- var_hte_three(ns, m, pw, var_y, var_x, a0, a1, r0, r1, rand)
    
    inside <- abs(d)/sqrt(var_4/(nc*ns*m)) - z
    
    power_result <- pnorm(inside)
    return(power_result)
  }
  
  nc_hte_three <- function(m, ns, power=0.9, pw, var_y, var_x, a0, a1_a0, r0, r1_r0, rand, d, a=0.05){
    
    a1 <- a1_a0*a0
    r1 <- r1_r0*r0
    
    s2_hte <- var_hte_three(ns, m, pw, var_y, var_x, a0, a1, r0, r1, rand)
    
    z_typeI <- qnorm(1-a/2)
    z_power <- qnorm(power)
    
    nc_result <- (( s2_hte * (z_typeI + z_power)^2 )/d^2)/(m*ns)
    
    inside <- abs(d)/sqrt(s2_hte/(ceiling(nc_result)*ns*m)) - z_typeI
    
    power_result <- pnorm(inside)
    
    return( list(nc=ceiling(nc_result), power_emp=power_result) )
  }
  
  var_hte_swd_cs <- function(n,m,J,var_y,var_x,a0,a1,r0,r1){
    ## n clusters
    ## m people per cluster-period
    ## J periods
    ## var_y outcome variance
    ## var_x covariate variance
    ## a0 within-period outcome ICC
    ## a1 between-period outcome ICC
    ## r0 within-period covariate ICC
    ## r1 between-period covariate ICC
    
    
    # eigenvalues for outcome #
    l1 <- 1 - a0
    l2 <- 1 + (m-1)*a0 - m*a1
    l3 <- 1 + (m-1)*a0 + (J-1)*m*a1
    
    # eigenvalues for covariate #
    e1 <- 1 - r0
    e2 <- 1 + (m-1)*r0 - m*r1
    e3 <- 1 + (m-1)*r0 + (J-1)*m*r1
    
    # intervention sequence #
    W0 <- matrix(0, ncol=J, nrow=J-1)
    W0[upper.tri(W0)] <- 1
    c_per_seq <- n/(J-1)
    rows_rep <- rep(1:nrow(W0),each=c_per_seq)
    W_full <- W0[rows_rep,]
    
    # covariance of intervention vector #
    # WW_sum <- NULL
    # for(i in seq(n)){
    #   WW <- W_full[i,] %*% t(W_full[i,])
    # 
    #   WW_sum <- WW_sum + WW
    # 
    # }
    U <- sum(W_full)
    V <- sum(rowSums(W_full)^2)
    W <- sum(colSums(W_full)^2)
    
    # W_sum <- colSums(W_full)
    # Wt_sum <- rowSums(W_full)
    
    # omega <- (1/n)*WW_sum - ( (1/n)*W_sum )%*%( (1/n)*Wt_sum )
    omega_11 <- 1/(n^2)*(n*V-U^2)
    omega_tr <- 1/(n^2)*(n*U - W)
    # omega_tr <- sum(diag(omega))
    # omega_11 <- t(rep(1,J)) %*% omega %*% rep(1,J)
    
    tau_w <- solve( (J-1)*omega_tr ) * ( omega_11 - omega_tr)
    theta <- J*(m-1)*(e1/l1) + (J-1)*(e2/l2) + (e3/l3)
    
    var_hte <- ((var_y/var_x)/(n*omega_tr)) * J^2/( (J-1)*(1-tau_w)*(e3-e2)*(solve(l2)-solve(l3)) + J*theta )
    
    if(any(c(l1,l2,l3) < 0)){
      stop("Variance is not positive definite - please enter new outcome correlation parameters.")
    }else if(any(c(e1,e2,e3)<0)) {
      stop("Variance is not positive definite - please enter new covariate correlation parameters.")
    }else{
      return(var_hte)#list(var=var_hte, tau=tau_w))  
    }
    
  }
  
  var_hte_swd_cc <- function(n,m,J,var_y,var_x,a0,a1,a2,r0){
    
    ## n clusters
    ## m people per cluster-period
    ## J periods
    ## var_y outcome variance
    ## var_x covariate variance
    ## a0 within-period outcome ICC
    ## a1 between-period outcome ICC
    ## a2 within-individual outcome ICC
    ## r0 within-period covariate ICC
    
    # eigenvalues for outcome #
    t1 <- 1-a0+a1-a2
    t2 <- 1-a0-(J-1)*(a1-a2)
    t3 <- 1+(m-1)*(a0-a1)-a2
    t4 <- 1+(m-1)*a0 + (J-1)*(m-1)*a1 + (J-1)*a2
    
    #eigenvalues for covariate #
    eta1 <- 1-r0
    eta2 <- 1+(m-1)*r0
    
    # intervention sequence #
    W0 <- matrix(0, ncol=J, nrow=J-1)
    W0[upper.tri(W0)] <- 1
    c_per_seq <- n/(J-1)
    rows_rep <- rep(1:nrow(W0),each=c_per_seq)
    W_full <- W0[rows_rep,]
    
    # covariance of intervention vector #
    U <- sum(W_full)
    V <- sum(rowSums(W_full)^2)
    W <- sum(colSums(W_full)^2)
    
    omega_11 <- 1/(n^2)*(n*V-U^2)
    omega_tr <- 1/(n^2)*(n*U - W)
    
    tau_w <- solve( (J-1)*omega_tr ) * ( omega_11 - omega_tr)
    theta <- (m-1)*eta1/t2 + eta2/t4
    
    var_hte <- ((var_y/var_x)/(n*omega_tr)) * J/((J-1)*(1-tau_w)*(((1/t3) - (1/t4))*eta2 + (m-1)*((1/t1) - (1/t2))*eta1) + J*theta)
    
    if(any(c(t1,t2,t3,t4) < 0)){
      stop("Variance is not positive definite - please enter new outcome correlation parameters.")
    }else if(any(c(eta1,eta2)<0)) {
      stop("Variance is not positive definite - please enter new covariate correlation parameters.")
    }else{
      return(var_hte)#list(var=var_hte, tau=tau_w))  
    }
    
  }
  
  power_hte_swd <- function(n, m, J, var_y, var_x, a0, a1_a0, a2=NULL, r0, r1_r0=NULL, d, cohort, a=0.05){
    z <- qnorm(1 - a/2)
    a1 <- a1_a0*a0
    
    if(cohort=="cross"){
      r1 <- r1_r0*r0
      
      var_4 <- var_hte_swd_cs(n,m,J,var_y,var_x,a0,a1,r0,r1)#var_hte_swd_cs_billy(n=n, m=m, J=J, var_y=var_y, var_x=var_x, d=d, a0=a0, a1=a1, r0=r0, r1=r1, a=a)#
    }else if (cohort=="closed"){
      var_4 <- var_hte_swd_cc(n,m,J,var_y,var_x,a0,a1,a2,r0)
    }
    
    
    
    inside <- abs(d)/sqrt(var_4) - z
    
    power_result <- pnorm(inside)
    return(power_result)
  }
  
  ns_hte_swd <- function(m,J,var_y,var_x,a0, a1_a0, a2=NULL, r0, r1_r0=NULL, d, cohort, a=0.05, power=0.8){
    a1 <- a1_a0*a0
    
    z <- qnorm(1 - a/2)
    
    power_emp <- 0
    nca <- 0
    
    while(power_emp < power){
      nca <- nca + 1
      n_result <- (J-1)*nca
      
      if(cohort=="cross"){
        
        r1 <- r1_r0*r0
        
        s2_hte <- var_hte_swd_cs(n_result,m,J,var_y,var_x,a0,a1,r0,r1)#var_hte_swd_cs_billy(n=n, m=m, J=J, var_y=var_y, var_x=var_x, d=d, a0=a0, a1=a1, r0=r0, r1=r1, a=a)#
      }else if (cohort=="closed"){
        s2_hte <- var_hte_swd_cc(n_result,m,J,var_y,var_x,a0,a1,a2,r0)
      }
      
      power_emp <- pnorm( abs(d)/sqrt(s2_hte) - z )
      
      #n_result <- (( s2_hte * (z_typeI + z_power)^2 )/d^2)
      
    }
    
    
    
    return( list(ns=n_result/(J-1), power_emp=power_emp) )
    
  }
  
  var_hte_ownDes_cs <- function(desmat,n_seq,m,var_y,var_x,a0,a1,r0,r1){
    
    ## desmat design matrix
    ## n_seq clusters per sequence
    ## m people per cluster-period
    ## var_y outcome variance
    ## var_x covariate variance
    ## a0 within-period outcome ICC
    ## a1 between-period outcome ICC
    ## r0 within-period covariate ICC
    ## r1 between-period covariate ICC
    
    # determine number of periods #
    J <- ncol(desmat)
    n <- n_seq*nrow(desmat)
    
    # eigenvalues for outcome #
    l1 <- 1 - a0
    l2 <- 1 + (m-1)*a0 - m*a1
    l3 <- 1 + (m-1)*a0 + (J-1)*m*a1
    
    # eigenvalues for covariate #
    e1 <- 1 - r0
    e2 <- 1 + (m-1)*r0 - m*r1
    e3 <- 1 + (m-1)*r0 + (J-1)*m*r1
    
    # intervention sequence #
    # W0 <- matrix(0, ncol=J, nrow=J-1)
    # W0[upper.tri(W0)] <- 1
    c_per_seq <- n_seq
    rows_rep <- rep(1:nrow(desmat),each=c_per_seq)
    W_full <- desmat[rows_rep,]
    
    if( any(is.na(W_full)) ) stop("Design matrix cannot have missing cluster-periods. Upload new design matrix.")
    
    ## CAN'T JUST PUT NA.RM=T TO HAVE NON-OBSERVED CLUSTER-PERIODS ##
    U <- sum(W_full)
    V <- sum(rowSums(W_full)^2)
    W <- sum(colSums(W_full)^2)
    
    # W_sum <- colSums(W_full)
    # Wt_sum <- rowSums(W_full)
    
    # omega <- (1/n)*WW_sum - ( (1/n)*W_sum )%*%( (1/n)*Wt_sum )
    omega_11 <- 1/(n^2)*(n*V-U^2)
    omega_tr <- 1/(n^2)*(n*U - W)
    # omega_tr <- sum(diag(omega))
    # omega_11 <- t(rep(1,J)) %*% omega %*% rep(1,J)
    
    tau_w <- solve( (J-1)*omega_tr ) * ( omega_11 - omega_tr)
    theta <- J*(m-1)*(e1/l1) + (J-1)*(e2/l2) + (e3/l3)
    
    var_hte <- ((var_y/var_x)/(n*omega_tr)) * J^2/( (J-1)*(1-tau_w)*(e3-e2)*(solve(l2)-solve(l3)) + J*theta )
    
    if(any(c(l1,l2,l3) < 0)){
      stop("Variance is not positive definite - please enter new outcome correlation parameters.")
    }else if(any(c(e1,e2,e3)<0)) {
      stop("Variance is not positive definite - please enter new covariate correlation parameters.")
    }else{
      return(var_hte)#list(var=var_hte, tau=tau_w))  
    }
    
  }
  
  var_hte_ownDes_cc <- function(desmat, n_seq,m,var_y,var_x,a0,a1,a2,r0){
    
    ## desmat design matrix
    ## n_seq clusters per sequence
    ## m people per cluster-period
    ## J periods
    ## var_y outcome variance
    ## var_x covariate variance
    ## a0 within-period outcome ICC
    ## a1 between-period outcome ICC
    ## a2 within-individual outcome ICC
    ## r0 within-period covariate ICC
    
    # determine number of periods #
    J <- ncol(desmat)
    n <- n_seq*nrow(desmat)
    
    # eigenvalues for outcome #
    t1 <- 1-a0+a1-a2
    t2 <- 1-a0-(J-1)*(a1-a2)
    t3 <- 1+(m-1)*(a0-a1)-a2
    t4 <- 1+(m-1)*a0 + (J-1)*(m-1)*a1 + (J-1)*a2
    
    #eigenvalues for covariate #
    eta1 <- 1-r0
    eta2 <- 1+(m-1)*r0
    
    # intervention sequence #
    # W0 <- matrix(0, ncol=J, nrow=J-1)
    # W0[upper.tri(W0)] <- 1
    c_per_seq <- n_seq
    rows_rep <- rep(1:nrow(desmat),each=c_per_seq)
    W_full <- desmat[rows_rep,]
    
    if( any(is.na(W_full)) ) stop("Design matrix cannot have missing cluster-periods. Upload new design matrix.")
    
    # covariance of intervention vector #
    U <- sum(W_full)
    V <- sum(rowSums(W_full)^2)
    W <- sum(colSums(W_full)^2)
    
    omega_11 <- 1/(n^2)*(n*V-U^2)
    omega_tr <- 1/(n^2)*(n*U - W)
    
    tau_w <- solve( (J-1)*omega_tr ) * ( omega_11 - omega_tr)
    theta <- (m-1)*eta1/t2 + eta2/t4
    
    var_hte <- ((var_y/var_x)/(n*omega_tr)) * J/((J-1)*(1-tau_w)*(((1/t3) - (1/t4))*eta2 + (m-1)*((1/t1) - (1/t2))*eta1) + J*theta)
    
    if(any(c(t1,t2,t3,t4) < 0)){
      stop("Variance is not positive definite - please enter new outcome correlation parameters.")
    }else if(any(c(eta1,eta2)<0)) {
      stop("Variance is not positive definite - please enter new covariate correlation parameters.")
    }else{
      return(var_hte)#list(var=var_hte, tau=tau_w))  
    }
    
  }
  
  power_hte_ownDes <- function(desmat, n_seq, m, var_y, var_x, a0, a1_a0, a2=NULL, r0, r1_r0=NULL, d, cohort, a=0.05){
    z <- qnorm(1 - a/2)
    a1 <- a1_a0*a0
    
    if(cohort=="cross"){
      r1 <- r1_r0*r0
      
      var_4 <- var_hte_ownDes_cs(desmat,n_seq,m,var_y,var_x,a0,a1,r0,r1)#var_hte_swd_cs_billy(n=n, m=m, J=J, var_y=var_y, var_x=var_x, d=d, a0=a0, a1=a1, r0=r0, r1=r1, a=a)#
    }else if (cohort=="closed"){
      var_4 <- var_hte_ownDes_cc(desmat,n_seq,m,var_y,var_x,a0,a1,a2,r0)
    }
    
    
    
    inside <- abs(d)/sqrt(var_4) - z
    
    power_result <- pnorm(inside)
    return(power_result)
  }
  
  ns_hte_ownDes <- function(desmat,m,var_y,var_x,a0, a1_a0, a2=NULL, r0, r1_r0=NULL, d, cohort, a=0.05, power=0.8){
    a1 <- a1_a0*a0
    
    z <- qnorm(1 - a/2)
    
    power_emp <- 0
    nca <- 0
    
    J <- ncol(desmat)
    
    while(power_emp < power){
      nca <- nca + 1
      #n_result <- nrow(desmat)*nca
      
      if(cohort=="cross"){
        
        r1 <- r1_r0*r0
        
        s2_hte <- var_hte_ownDes_cs(desmat, nca,m,var_y,var_x,a0,a1,r0,r1)#var_hte_swd_cs_billy(n=n, m=m, J=J, var_y=var_y, var_x=var_x, d=d, a0=a0, a1=a1, r0=r0, r1=r1, a=a)#
      }else if (cohort=="closed"){
        s2_hte <- var_hte_ownDes_cc(desmat,nca,m,var_y,var_x,a0,a1,a2,r0)
      }
      
      power_emp <- pnorm( abs(d)/sqrt(s2_hte) - z )
      
      #n_result <- (( s2_hte * (z_typeI + z_power)^2 )/d^2)
      
    }
    
    
    
    return( list(ns=nca, power_emp=power_emp) )
    
  }
  
  var_hte_crossover_cs <- function(n_seq,m,J,var_y,var_x,a0,a1,r0,r1){
    
    ## n_seq clusters per sequence
    ## m people per cluster-period
    ## J periods
    ## var_y outcome variance
    ## var_x covariate variance
    ## a0 within-period outcome ICC
    ## a1 between-period outcome ICC
    ## r0 within-period covariate ICC
    ## r1 between-period covariate ICC
    
    # determine total clusters - crossovers always have 2 sequences#
    n <- n_seq*2
    
    # eigenvalues for outcome #
    l1 <- 1 - a0
    l2 <- 1 + (m-1)*a0 - m*a1
    l3 <- 1 + (m-1)*a0 + (J-1)*m*a1
    
    # eigenvalues for covariate #
    e1 <- 1 - r0
    e2 <- 1 + (m-1)*r0 - m*r1
    e3 <- 1 + (m-1)*r0 + (J-1)*m*r1
    
    # intervention sequence #
    if( J %% 2 == 1 ){
      W0 <- matrix(data=c(rep(c(0,1),J)), nrow=2, ncol=J, byrow=TRUE)
    }else if( J %% 2 == 0 ){
      W0 <- matrix(data=c(rep(c(0,1),J/2), rep(c(1,0), J/2)), nrow=2, ncol=J, byrow=TRUE)
    } 
    c_per_seq <- n_seq
    rows_rep <- rep(1:nrow(W0),each=c_per_seq)
    W_full <- W0[rows_rep,]
    
    U <- sum(W_full)
    V <- sum(rowSums(W_full)^2)
    W <- sum(colSums(W_full)^2)
    
    omega_11 <- 1/(n^2)*(n*V-U^2)
    omega_tr <- 1/(n^2)*(n*U - W)
    
    tau_w <- solve( (J-1)*omega_tr ) * ( omega_11 - omega_tr)
    theta <- J*(m-1)*(e1/l1) + (J-1)*(e2/l2) + (e3/l3)
    
    var_hte <- ((var_y/var_x)/(n*omega_tr)) * J^2/( (J-1)*(1-tau_w)*(e3-e2)*(solve(l2)-solve(l3)) + J*theta )
    
    if(any(c(l1,l2,l3) < 0)){
      stop("Variance is not positive definite - please enter new outcome correlation parameters.")
    }else if(any(c(e1,e2,e3)<0)) {
      stop("Variance is not positive definite - please enter new covariate correlation parameters.")
    }else{
      return(var_hte)#list(var=var_hte, tau=tau_w))  
    }
    
  }
  
  var_hte_crossover_cc <- function(n_seq,m,J,var_y,var_x,a0,a1,a2,r0){
    
    ## n_seq clusters per sequence
    ## m people per cluster-period
    ## J periods
    ## var_y outcome variance
    ## var_x covariate variance
    ## a0 within-period outcome ICC
    ## a1 between-period outcome ICC
    ## a2 within-individual outcome ICC
    ## r0 within-period covariate ICC
    
    # determine total number of clusters - crossovers always have 2 #
    n <- n_seq*2
    
    # eigenvalues for outcome #
    t1 <- 1-a0+a1-a2
    t2 <- 1-a0-(J-1)*(a1-a2)
    t3 <- 1+(m-1)*(a0-a1)-a2
    t4 <- 1+(m-1)*a0 + (J-1)*(m-1)*a1 + (J-1)*a2
    
    #eigenvalues for covariate #
    eta1 <- 1-r0
    eta2 <- 1+(m-1)*r0
    
    # intervention sequence #
    if( J %% 2 == 1 ){
      W0 <- matrix(data=c(rep(c(0,1),J)), nrow=2, ncol=J, byrow=TRUE)
    }else if( J %% 2 == 0 ){
      W0 <- matrix(data=c(rep(c(0,1),J/2), rep(c(1,0), J/2)), nrow=2, ncol=J, byrow=TRUE)
    } 
    c_per_seq <- n_seq
    rows_rep <- rep(1:nrow(W0),each=c_per_seq)
    W_full <- W0[rows_rep,]
    
    # covariance of intervention vector #
    U <- sum(W_full)
    V <- sum(rowSums(W_full)^2)
    W <- sum(colSums(W_full)^2)
    
    omega_11 <- 1/(n^2)*(n*V-U^2)
    omega_tr <- 1/(n^2)*(n*U - W)
    
    tau_w <- solve( (J-1)*omega_tr ) * ( omega_11 - omega_tr)
    theta <- (m-1)*eta1/t2 + eta2/t4
    
    var_hte <- ((var_y/var_x)/(n*omega_tr)) * J/((J-1)*(1-tau_w)*(((1/t3) - (1/t4))*eta2 + (m-1)*((1/t1) - (1/t2))*eta1) + J*theta)
    
    if(any(c(t1,t2,t3,t4) < 0)){
      stop("Variance is not positive definite - please enter new outcome correlation parameters.")
    }else if(any(c(eta1,eta2)<0)) {
      stop("Variance is not positive definite - please enter new covariate correlation parameters.")
    }else{
      return(var_hte)#list(var=var_hte, tau=tau_w))  
    }
    
  }
  
  power_hte_crossover <- function(n_seq, m, J, var_y, var_x, a0, a1_a0, a2=NULL, r0, r1_r0=NULL, d, cohort, a=0.05){
    z <- qnorm(1 - a/2)
    a1 <- a1_a0*a0
    
    if(cohort=="cross"){
      r1 <- r1_r0*r0
      
      var_4 <- var_hte_crossover_cs(n_seq,m,J,var_y,var_x,a0,a1,r0,r1)#var_hte_swd_cs_billy(n=n, m=m, J=J, var_y=var_y, var_x=var_x, d=d, a0=a0, a1=a1, r0=r0, r1=r1, a=a)#
    }else if (cohort=="closed"){
      var_4 <- var_hte_crossover_cc(n_seq,m,J,var_y,var_x,a0,a1,a2,r0)
    }
    
    
    
    inside <- abs(d)/sqrt(var_4) - z
    
    power_result <- pnorm(inside)
    return(power_result)
  }
  
  ns_hte_crossover <- function(m,J, var_y,var_x,a0, a1_a0, a2=NULL, r0, r1_r0=NULL, d, cohort, a=0.05, power=0.8){
    a1 <- a1_a0*a0
    
    z <- qnorm(1 - a/2)
    
    power_emp <- 0
    nca <- 0
    
    
    while(power_emp < power){
      nca <- nca + 1
      #n_result <- nrow(desmat)*nca
      
      if(cohort=="cross"){
        
        r1 <- r1_r0*r0
        
        s2_hte <- var_hte_crossover_cs(nca,m,J,var_y,var_x,a0,a1,r0,r1)#var_hte_swd_cs_billy(n=n, m=m, J=J, var_y=var_y, var_x=var_x, d=d, a0=a0, a1=a1, r0=r0, r1=r1, a=a)#
      }else if (cohort=="closed"){
        s2_hte <- var_hte_crossover_cc(nca,m,J,var_y,var_x,a0,a1,a2,r0)
      }
      
      power_emp <- pnorm( abs(d)/sqrt(s2_hte) - z )
      
      #n_result <- (( s2_hte * (z_typeI + z_power)^2 )/d^2)
      
    }
    
    
    
    return( list(ns=nca, power_emp=power_emp) )
    
  }
  
  var_irgt <- function(m1,m0,n1,n0, oicc1,oicc0, cicc=0, var_y1,var_y0, var_x){
    n <- n1+n0
    pi <- n1/n
    
    top_trt <- var_y1 * (1-oicc1) * (1 + (m1-1) * oicc1)
    bottom_trt <- var_x * pi * m1 * (1 + (m1-2) * oicc1 - (m1-1) * cicc * oicc1 )
    
    top_ctrl <- var_y0 * (1-oicc0) * (1 + (m0-1) * oicc0)
    bottom_ctrl <- var_x * (1-pi) * m0 * (1 + (m0-2) * oicc0 - (m0-1) * cicc * oicc0 )
    
    var_result <- top_trt/bottom_trt + top_ctrl/bottom_ctrl
    return(var_result)
  }
  
  power_irgt <- function(m1,m0,n1,n0, oicc1,oicc0, cicc=0, var_y1,var_y0, var_x, d, a=0.05){
    n <- n1+n0
    
    z <- qnorm(1 - a/2)
    
    var_4 <- var_irgt(m1,m0,n1,n0, oicc1,oicc0, cicc, var_y1,var_y0, var_x)
    
    inside <- abs(d)/sqrt(var_4/n) - z
    
    power_result <- pnorm(inside)
    return(power_result)
  }
  
  n_irgt <- function(m1,m0, oicc1,oicc0, cicc=0, var_y1,var_y0, var_x, d, pi_samp=NULL, pi=NULL, a=0.05, power=0.8){
    
    if(is.null(pi) & is.null(pi_samp)){
      
      pi_result <- 1/( 1 + sqrt( (var_y0 * m1 * (1-oicc0) * (1 + (m0-1) * oicc0) * (1 + (m1-2) * oicc1 - (m1-1) * cicc * oicc1))/(var_y1 * m0 * (1-oicc1) * (1 + (m1-1)* oicc1) * (1 + (m0-2) * oicc0 - (m0-1) * cicc*oicc0)) ) )
      
    }else if(!is.null(pi_samp)){
      
      pi_result <- (1-pi_samp)*m0/(pi_samp*m1+(1-pi_samp)*m0)
      
    }else{
      pi_result <- pi
    }
    
    
    top_trt <- var_y1 * (1-oicc1) * (1 + (m1-1) * oicc1)
    bottom_trt <- var_x * pi_result * m1 * (1 + (m1-2) * oicc1 - (m1-1) * cicc * oicc1 )
    
    top_ctrl <- var_y0 * (1-oicc0) * (1 + (m0-1) * oicc0)
    bottom_ctrl <- var_x * (1-pi_result) * m0 * (1 + (m0-2) * oicc0 - (m0-1) * cicc * oicc0 )
    
    s2_hte <- top_trt/bottom_trt + top_ctrl/bottom_ctrl
    
    z_typeI <- qnorm(1-a/2)
    z_power <- qnorm(power)
    
    n_result <- ceiling( ( s2_hte * (z_typeI + z_power)^2 )/d^2 )
    n1_result <- ceiling( n_result*pi_result )
    n0_result <- n_result - n1_result
    
    inside <- abs(d)/sqrt(s2_hte/n_result) - z_typeI
    
    power_result <- pnorm(inside)
    
    return( list(n=n_result, n1=n1_result, n0=n0_result, pi=pi_result, power_emp = power_result ) )
  }
  
  designMatrix <- function(design, periods=NULL, steps=NULL, file=NULL){
    if(design %in% c("parallel","three_level", "irgt", "het_two")){
      
      matrix(c(0,1),nrow=2,ncol=1)
      
    }else if(design == "parallel_m"){
    
        matrix(data=c(rep(c(0,1),periods)), nrow=2, ncol=periods, byrow=FALSE)
    
    }else if(design == "crossover_2" | design == "crossover_m"){
      
      if( periods %% 2 == 1 ){
        matrix(data=c(rep(c(0,1),periods)), nrow=2, ncol=periods, byrow=TRUE)
      }else if( periods %% 2 == 0 ){
        matrix(data=c(rep(c(0,1),periods/2), rep(c(1,0), periods/2)), nrow=2, ncol=periods, byrow=TRUE)
      } 
      
    }else if(design == "SWD"){
      
      W0 <- matrix(0, ncol=periods, nrow=steps)
      W0[upper.tri(W0)] <- 1
      W0
      
    }else if(design == "upload"){
      read.csv(file$datapath, header=FALSE)
    }
  }
  
  #### OUTPUT ####
  ## upload own observes ##
  file1 <- reactive({input$file1})
  
  observeEvent(input$reset, {
    file1 <- NULL
    reset('file1')
  })
  
  
  ## update UI sidebar options ##
  trial_react <- reactive({input$trial})
  cohort_react <- reactive({input$cohort})
  observeEvent({input$trial#trial_react()
    input$cohort}, {
      #cohort_react()}, {
      
      if(input$trial == "parallel" || input$trial == "irgt" || input$trial== "het_two"){
        # update plot display options #
        updateRadioButtons(session, "plot_display",
                           label="Plot display",
                           choiceNames = c("Cluster size vs Power",
                                           "Number of clusters vs Power",
                                           "Cluster size vs number of clusters"),#,
                           #"Number of clusters vs Cluster size (fixed power)"),
                           choiceValues=c("m_v_power", "n_v_power", "fixed_power")
        )
        
      }else if(input$trial == "three_level"){
        # update plot display options #
        updateRadioButtons(session, "plot_display",
                           label="Plot display",
                           choiceNames = c("Individuals vs Power (fixed clusters and subclusters)",
                                           #"Subclusters vs Power (fixed clusters and individuals)",
                                           "Clusters vs Power (fixed subclusters and individuals)",
                                           #"Subclusters vs Clusters (fixed power and individuals)",
                                           "Individuals vs Clusters (fixed power and subclusters)"#,
                                           #"Individuals vs Subclusters (fixed power and clusters)"
                           ),#,
                           #"Number of clusters vs Cluster size (fixed power)"),
                           choiceValues=c("m_v_power", #"ns_v_power",
                                          "nc_v_power", #"ns_v_nc",
                                          "m_v_nc"#,"m_v_ns"
                           )
        )
        
      }else if((input$trial == 'parallel_m' || input$trial == 'crossover_2' || input$trial == 'crossover_m' || input$trial == "SWD" || input$trial == "upload" ) & input$cohort == "cross"){
        updateRadioButtons(session, "plot_display",
                           label="Plot display",
                           choiceNames = c("Cluster size (per period) vs Power",
                                           "Number of clusters (per sequence) vs Power",
                                           "Cluster size (per period) vs number of clusters (per sequence)"),#,
                           #"Number of clusters vs Cluster size (fixed power)"),
                           choiceValues=c("m_v_power", "n_v_power", "fixed_power")
        )
        
      }else if((input$trial == 'parallel_m' || input$trial == 'crossover_2' || input$trial == 'crossover_m' || input$trial == "SWD" || input$trial == "upload" ) & input$cohort == "closed"){
        updateRadioButtons(session, "plot_display",
                           label="Plot display",
                           choiceNames = c("Total cluster size vs Power",
                                           "Number of clusters (per sequence) vs Power",
                                           "Total cluster size vs number of clusters (per sequence)"),#,
                           #"Number of clusters vs Cluster size (fixed power)"),
                           choiceValues=c("m_v_power", "n_v_power", "fixed_power")
        )
        
      }
    })
  
  # restrict irgt and heterogeneous CRTs to only continuous outcomes #
  observeEvent(input$trial,{#trial_react(),{
    if(input$trial == 'irgt' || input$trial == 'het_two'){
      updateRadioButtons(session, "outcome",
                         label="Outcome type",
                         choiceNames = c("Continuous"),
                         choiceValues=c("continuous")
      )
    }
  })
  # update choices for within-plot consistency in het CRT #
  observeEvent({input$trial
    input$icc_constant},{#trial_react(),{
      if(input$trial == 'het_two' & input$icc_constant == 'oicc1'){
        updateRadioButtons(session, "icc_constant_within",
                           label="ICC to stay constant within plot",
                           choiceNames = c(#"Treatment-arm outcome ICC",
                             "Control-arm outcome ICC",
                             "Covariate ICC"),
                           choiceValues = c(#"oicc1",
                             "oicc0",
                             "cicc"))
      }else if(input$trial == 'het_two' & input$icc_constant == 'oicc0'){
        updateRadioButtons(session, "icc_constant_within",
                           label="ICC to stay constant within plot",
                           choiceNames = c("Treatment-arm outcome ICC",
                                           #"Control-arm outcome ICC",
                                           "Covariate ICC"),
                           choiceValues = c("oicc1",
                                            #"oicc0",
                                            "cicc")) 
      }else if(input$trial == 'het_two' & input$icc_constant == 'cicc'){
        updateRadioButtons(session, "icc_constant_within",
                           label="ICC to stay constant within plot",
                           choiceNames = c("Treatment-arm outcome ICC",
                                           "Control-arm outcome ICC"),#,
                           #"Covariate ICC"),
                           choiceValues = c("oicc1",
                                            "oicc0"))
        #"cicc")) 
      }
    })
  
  ## power over number of clusters ##
  outcome_type_react <- reactive({input$outcome})
  covar_type_react <- reactive({input$covar})
  prop_control_react <- reactive({input$prop_control})
  prop_trt_react <- reactive({input$prop_trt}) 
  prop_covar_react <- reactive({input$prop_covar})
  sd_outcome_react <- reactive({input$sd_outcome})
  sd_outcome1_react <- reactive({input$sd_outcome1})
  sd_outcome0_react <- reactive({input$sd_outcome0})
  sd_covar_react <- reactive({input$sd_covar})
  plot_display_react <- reactive({input$plot_display})
  
  ## parallel reacts ##
  #icc_display_react <- reactive({input$icc_display})
  oicc_min <- reactive({input$oicc_range[1]})
  oicc_max <- reactive({input$oicc_range[2]})
  oicc_est <- reactive({input$oicc_est})
  cicc_min <- reactive({input$cicc_range[1]})
  cicc_max <- reactive({input$cicc_range[2]})
  cicc_est <- reactive({input$cicc_est})
  
  sensitivity_react <- reactive({input$icc_sensitivity})
  
  ## three-level reacts ##
  oicc_wsub_min_three <- reactive({input$oicc_wsub_range_three[1]})
  oicc_wsub_max_three <- reactive({input$oicc_wsub_range_three[2]})
  oicc_wsub_est_three <- reactive({input$oicc_wsub_est_three})
  
  oicc_ratio_min_three <- reactive({input$oicc_ratio_range_three[1]})
  oicc_ratio_max_three <- reactive({input$oicc_ratio_range_three[2]})
  oicc_ratio_est_three <- reactive({input$oicc_ratio_est_three})
  
  cicc_wsub_min_three <- reactive({input$cicc_wsub_range_three[1]})
  cicc_wsub_max_three <- reactive({input$cicc_wsub_range_three[2]})
  cicc_wsub_est_three <- reactive({input$cicc_wsub_est_three})
  
  cicc_ratio_min_three <- reactive({input$cicc_ratio_range_three[1]})
  cicc_ratio_max_three <- reactive({input$cicc_ratio_range_three[2]})
  cicc_ratio_est_three <- reactive({input$cicc_ratio_est_three})
  
  sensitivity_three_react <- reactive({input$icc_sensitivity_three})
  
  ## crossovers & SW-CRT & upload reacts ##
  cohort <- reactive({input$cohort})
  
  oicc_wperiod_min_swd <- reactive({input$oicc_wperiod_range_swd[1]})
  oicc_wperiod_max_swd <- reactive({input$oicc_wperiod_range_swd[2]})
  oicc_wperiod_est_swd <- reactive({input$oicc_wperiod_est_swd})
  
  oicc_wperiod_min_swd_cc <- reactive({input$oicc_wperiod_range_swd_cc[1]})
  oicc_wperiod_max_swd_cc <- reactive({input$oicc_wperiod_range_swd_cc[2]})
  oicc_wperiod_est_swd_cc <- reactive({input$oicc_wperiod_est_swd_cc})
  
  oicc_ratio_min_swd <- reactive({input$oicc_ratio_range_swd[1]})
  oicc_ratio_max_swd <- reactive({input$oicc_ratio_range_swd[2]})
  oicc_ratio_est_swd <- reactive({input$oicc_ratio_est_swd})
  
  oicc_ratio_min_swd_cc <- reactive({input$oicc_ratio_range_swd_cc[1]})
  oicc_ratio_max_swd_cc <- reactive({input$oicc_ratio_range_swd_cc[2]})
  oicc_ratio_est_swd_cc <- reactive({input$oicc_ratio_est_swd_cc})
  
  oicc_windiv_min_swd_cc <- reactive({input$oicc_windiv_range_swd_cc[1]})
  oicc_windiv_max_swd_cc <- reactive({input$oicc_windiv_range_swd_cc[2]})
  oicc_windiv_est_swd_cc <- reactive({input$oicc_windiv_est_swd_cc})
  
  cicc_wperiod_min_swd <- reactive({input$cicc_wperiod_range_swd[1]})
  cicc_wperiod_max_swd <- reactive({input$cicc_wperiod_range_swd[2]})
  cicc_wperiod_est_swd <- reactive({input$cicc_wperiod_est_swd})
  
  cicc_wperiod_min_swd_cc <- reactive({input$cicc_wperiod_range_swd_cc[1]})
  cicc_wperiod_max_swd_cc <- reactive({input$cicc_wperiod_range_swd_cc[2]})
  cicc_wperiod_est_swd_cc <- reactive({input$cicc_wperiod_est_swd_cc})
  
  cicc_ratio_min_swd <- reactive({input$cicc_ratio_range_swd[1]})
  cicc_ratio_max_swd <- reactive({input$cicc_ratio_range_swd[2]})
  cicc_ratio_est_swd <- reactive({input$cicc_ratio_est_swd})
  
  sensitivity_swd_react <- reactive({input$icc_sensitivity_swd})
  
  ## irgt reactives ##
  clustering_irgt <- reactive({input$control_cluster})
  
  oicc_trt_min_irgt <- reactive({input$oicc_trt_range_irgt[1]})
  oicc_trt_max_irgt <- reactive({input$oicc_trt_range_irgt[2]})
  oicc_trt_est_irgt <- reactive({input$oicc_trt_est_irgt})
  
  oicc_ctrl_min_irgt <- reactive({input$oicc_ctrl_range_irgt[1]})
  oicc_ctrl_max_irgt <- reactive({input$oicc_ctrl_range_irgt[2]})
  oicc_ctrl_est_irgt <- reactive({input$oicc_ctrl_est_irgt})
  
  sensitivity_irgt_react <- reactive({input$icc_sensitivity_irgt})
  
  ## het crt reactives ##
  oicc_trt_min_het <- reactive({input$oicc_trt_range_het[1]})
  oicc_trt_max_het <- reactive({input$oicc_trt_range_het[2]})
  oicc_trt_est_het <- reactive({input$oicc_trt_est_het})
  
  oicc_ctrl_min_het <- reactive({input$oicc_ctrl_range_het[1]})
  oicc_ctrl_max_het <- reactive({input$oicc_ctrl_range_het[2]})
  oicc_ctrl_est_het <- reactive({input$oicc_ctrl_est_het})
  
  cicc_min_het <- reactive({input$cicc_range_het[1]})
  cicc_max_het <- reactive({input$cicc_range_het[2]})
  cicc_est_het <- reactive({input$cicc_est_het})
  
  sensitivity_het_react <- reactive({input$icc_sensitivity_het})
  
  # window reactive #
  window_width <- reactive({input$dimension[2]})
  
  
  ## power over cluster size ##
  output$powerPlot_hte <- renderPlotly({
    
    
    
    colors_plot <- c("#66c2a5","#fc8d62","#8da0cb")
    
    # determine trial type #
    #### PARALLEL ####
    if(trial_react()=="parallel"){
      
      # determine effect sizes and variances depending on outcome/covariate type #
      if(outcome_type_react() =="continuous"){
        var_y <- (sd_outcome_react())^2
        d <- input$mean_diff_HTE
        
        if(covar_type_react()=="continuous"){
          # continuous outcome and covar power
          var_x <- (sd_covar_react())^2
          
        }else if(covar_type_react() == "binary"){
          # continuous outcome binary covar power
          var_x <- (prop_covar_react())*(1-prop_covar_react())
        }
        
      }else if(outcome_type_react() == "binary"){
        var_y <- ((prop_control_react()*(1-prop_control_react())) + (prop_trt_react()*(1-prop_trt_react())))/2
        #d <- (input$prop_trt - input$prop_control)
        
        if(covar_type_react() == "continuous"){
          # binary outcome continuous covar power
          var_x <-(sd_covar_react())^2
          
        }else if(covar_type_react() == "binary"){
          # binary outcome and covar power #
          var_x <- (prop_covar_react())*(1-prop_covar_react())
          
        }
      }
      
      if(plot_display_react() == "m_v_power"){
        m_range <- seq(input$m_range[1],input$m_range[2])
        if(sensitivity_react() == "est_only"){
          df_power <- expand.grid(n=input$n, m=m_range,
                                  oicc=#c(oicc_min(),
                                    oicc_est(),
                                  # oicc_max()),
                                  cicc=#c(cicc_min(),
                                    cicc_est(),
                                  #   cicc_max()),
                                  var_y=var_y, #(input$sd_outcome)^2,
                                  var_x=var_x,#(input$sd_covar)^2,
                                  var_w=(input$w)*(1-input$w),
                                  d=input$mean_diff_HTE, a=input$sig)
        }else{
          df_power <- expand.grid(n=input$n, m=m_range,
                                  oicc=c(oicc_min(), oicc_est(), oicc_max()),
                                  cicc=c(cicc_min(), cicc_est(), cicc_max()),
                                  var_y=var_y, #(input$sd_outcome)^2,
                                  var_x=var_x,#(input$sd_covar)^2,
                                  var_w=(input$w)*(1-input$w),
                                  d=input$mean_diff_HTE, a=input$sig)
        }
        
        
        power_hte_col <- rep(NA, nrow(df_power))
        for(i in seq(nrow(df_power))){
          power_hte_col[i] <- power_hte(oicc=df_power[i,"oicc"], cicc=df_power[i,"cicc"],
                                        m=df_power[i,"m"],n=df_power[i,"n"],
                                        var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                        var_w=df_power[i,"var_w"],
                                        d=df_power[i,"d"], a=df_power[i,"a"])
        }
        
        df_power <- cbind(df_power, power_hte_col)
        
        if(sensitivity_react() == "est_only"){
          p_est <- df_power %>%
            as.data.frame() %>%
            dplyr::filter(oicc == oicc_est(),#oicc_unique[2],
                          cicc == cicc_est()#cicc_unique[2]
            ) %>%
            mutate(cicc=factor(cicc)) %>%
            plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                    linetype=~cicc, color=~cicc,legendgroup=~cicc,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("outcome ICC: ", oicc, "<br>covariate ICC: ", cicc,
                                 "<br>Clusters (n):", n, "<br>Cluster size (m): ", m,
                                 "<br>HTE power: ", round(power_hte_col,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("o-ICC = ", oicc_est()),
                   xaxis=list(title="Cluster size (m)"),
                   yaxis=list(title="HTE Power"),
                   legend=list(title=list(text="covariate ICC")),
                   margin=0.001)
          
          subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                  margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
            layout(title = list(
              text='Cluster size vs HTE power',
              font=list(size=17)
            ),
            #margin=list(pad=50),
            annotations = list(
              list(
                x = 0.485,
                y = 1.03,
                text = paste0("<i>Assumed outcome ICC (", oicc_est(),") and covariate ICC (", cicc_est(),")</i>"),
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = FALSE, font=list(size=13)
              )
            ),
            legend=list(orientation="h",
                        yanchor="center",
                        y=0.25,
                        x=0.5)
            )
          
        }else if(sensitivity_react() == "sensitivity"){
          if(input$icc_display == "oICC_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_x$")
            
            oicc_unique <- sort(unique(df_power[,"oicc"]))
            cicc_unique <- sort(unique(df_power[,"cicc"]))
            p <- vector(mode="list", length=length(oicc_unique))
            
            for(i in seq(length(oicc_unique))){
              
              if(i != 3){
                p[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(cicc=factor(cicc)) %>%
                  dplyr::filter(oicc == oicc_unique[i]) %>%
                  plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                          linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("outcome ICC: ", oicc, "<br>covariate ICC: ", cicc,
                                       "<br>Clusters (n):", n, "<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC = ", oicc_unique[i]),
                         xaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="covariate ICC")),
                         margin=0.01)
              }else{
                p[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(cicc=factor(cicc)) %>%
                  dplyr::filter(oicc == oicc_unique[i]) %>%
                  plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                          linetype=~cicc, color=~cicc,legendgroup=~cicc,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("outcome ICC: ", oicc, "<br>covariate ICC: ", cicc,
                                       "<br>Clusters (n):", n, "<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC = ", oicc_unique[i]),
                         xaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="covariate ICC")),
                         margin=0.01)
              }
              
            }
            
            
            subplot(p[[1]],p[[2]],p[[3]], nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size vs HTE power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.03,
                  text = paste0("<i>Minimum outcome ICC (", oicc_unique[1],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE,
                  font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.43,
                  text = paste0("<i>Maximum outcome ICC (", oicc_unique[3],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.03,
                  text = paste0("<i>Assumed outcome ICC (", oicc_unique[2],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
            
            
            #p[[1]]
            
          }else if(input$icc_display == "cICC_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
            cicc_unique <- sort(unique(df_power[,"cicc"]))
            p <- vector(mode="list", length=length(cicc_unique))
            
            for(i in seq(length(cicc_unique))){
              
              if(i != 3){
                p[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(oicc=factor(oicc)) %>%
                  dplyr::filter(cicc == cicc_unique[i]) %>%
                  plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc,
                          linetype=~oicc, color=~oicc,legendgroup=~oicc, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("outcome ICC: ", oicc, "<br>covariate ICC: ", cicc,
                                       "<br>Clusters (n):", n, "<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)), height=input$dimension[2]*0.8) %>%
                  layout(title=paste("c-ICC = ", cicc_unique[i]),
                         xaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="outcome ICC")),
                         margin=0.01)
              }else{
                p[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(oicc=factor(oicc)) %>%
                  dplyr::filter(cicc == cicc_unique[i]) %>%
                  plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc,
                          linetype=~oicc, color=~oicc,legendgroup=~oicc,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("outcome ICC: ", oicc, "<br>covariate ICC: ", cicc,
                                       "<br>Clusters (n):", n, "<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)), height=input$dimension[2]*0.8) %>%
                  layout(title=paste("c-ICC = ", cicc_unique[i]),
                         xaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="outcome ICC")),
                         margin=0.01)
              }
              
            }
            
            subplot(p[[1]],p[[2]],p[[3]], nrows=2,
                    margin = 0.07, titleX=T, titleY=T) %>%
              layout(title = list(
                text='Cluster size vs HTE power',
                font=list(size=17)
              ),
              #margin = list(t=50,b=50,pad=20),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.03,
                  text = paste0("Mimimum covariate ICC (", cicc_unique[1],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.43,
                  text = paste0("<i>Maximum covariate ICC (", cicc_unique[2],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.03,
                  text = paste0("<i>Assumed covariate ICC (", cicc_unique[3],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
            
          }# end constant icc if/else
          
          
        }# end sensitivity if/else
        
        
        
      }else if(plot_display_react() == "n_v_power"){
        
        n_range <- seq(input$n_range[1],input$n_range[2])
        if(sensitivity_react() == "est_only"){
          df_power <- expand.grid(n=n_range, m=input$m,
                                  oicc=#c(oicc_min(), 
                                    oicc_est(),
                                  # oicc_max()),
                                  cicc=#c(cicc_min(),
                                    cicc_est(),
                                  # cicc_max()),
                                  var_y=var_y,
                                  var_x=var_x,
                                  var_w=(input$w)*(1-input$w),
                                  d=input$mean_diff_HTE, a=input$sig)
        }else{
          df_power <- expand.grid(n=n_range, m=input$m,
                                  oicc=c(oicc_min(), oicc_est(), oicc_max()),
                                  cicc=c(cicc_min(), cicc_est(), cicc_max()),
                                  var_y=var_y,
                                  var_x=var_x,
                                  var_w=(input$w)*(1-input$w),
                                  d=input$mean_diff_HTE, a=input$sig)
        }
        
        
        power_hte_col <- rep(NA, nrow(df_power))
        for(i in seq(nrow(df_power))){
          power_hte_col[i] <- power_hte(oicc=df_power[i,"oicc"], cicc=df_power[i,"cicc"],
                                        m=df_power[i,"m"],n=df_power[i,"n"],
                                        var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                        var_w=df_power[i,"var_w"],
                                        d=df_power[i,"d"], a=df_power[i,"a"])
        }
        
        df_power <- cbind(df_power, power_hte_col)
        
        if(sensitivity_react() == "est_only"){
          
          p_est <- df_power %>%
            as.data.frame() %>%
            dplyr::filter(oicc == oicc_est(),#oicc_unique[2],
                          cicc == cicc_est()#cicc_unique[2]
            ) %>%
            mutate(cicc=factor(cicc)) %>%
            plot_ly(x=~n,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                    linetype=~cicc, color=~cicc,legendgroup=~cicc,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("outcome ICC: ", oicc, "<br>covariate ICC: ", cicc,
                                 "<br>Clusters (n):", n, "<br>Cluster size (m): ", m,
                                 "<br>HTE power: ", round(power_hte_col,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("o-ICC = ", oicc_est()),
                   xaxis=list(title="Number of Clusters (n)"),
                   yaxis=list(title="HTE Power"),
                   legend=list(title=list(text="covariate ICC")),
                   margin=0.001)
          
          subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                  margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
            layout(title = list(
              text='Number of clusters vs HTE power',
              font=list(size=17)
            ),
            #margin=list(pad=50),
            annotations = list(
              list(
                x = 0.485,
                y = 1.03,
                text = paste0("<i>Assumed outcome ICC (", oicc_est(),") and covariate ICC (", cicc_est(),")</i>"),
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = FALSE, font=list(size=13)
              )
            ),
            legend=list(orientation="h",
                        yanchor="center",
                        y=0.25,
                        x=0.5)
            )
          
        }else if(sensitivity_react() == "sensitivity"){
          
          if(input$icc_display == "oICC_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_x$")
            
            oicc_unique <- sort(unique(df_power[,"oicc"]))
            cicc_unique <- sort(unique(df_power[,"cicc"]))
            p <- vector(mode="list", length=length(oicc_unique))
            
            for(i in seq(length(oicc_unique))){
              
              if(i != 3){
                p[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(cicc=factor(cicc)) %>%
                  dplyr::filter(oicc == oicc_unique[i]) %>%
                  plot_ly(x=~n,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                          linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("outcome ICC: ", oicc, "<br>covariate ICC: ", cicc,
                                       "<br>Clusters (n):", n, "<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)), height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC = ", oicc_unique[i]),
                         xaxis=list(title="Number of clusters (n)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="covariate ICC")),
                         margin=0.01)
              }else{
                p[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(cicc=factor(cicc)) %>%
                  dplyr::filter(oicc == oicc_unique[i]) %>%
                  plot_ly(x=~n,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                          linetype=~cicc, color=~cicc,legendgroup=~cicc,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("outcome ICC: ", oicc, "<br>covariate ICC: ", cicc,
                                       "<br>Clusters (n):", n, "<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)), height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC = ", oicc_unique[i]),
                         xaxis=list(title="Number of clusters (n)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="covariate ICC")),
                         margin=0.01)
              }
              
            }
            
            subplot(p[[1]],p[[2]],p[[3]], nrows=2,
                    margin = 0.07, titleX=T, titleY=T) %>%
              layout(title = list(
                text='Number of clusters vs HTE power',
                font=list(size=17)
              ),
              #margin = list(t=50,b=50,pad=20),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.03,
                  text = paste0("<i>Minimum outcome ICC (", oicc_unique[1],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.43,
                  text = paste0("<i>Maximum outcome ICC (", oicc_unique[2],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.03,
                  text = paste0("<i>Assumed outcome ICC (", oicc_unique[3],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
            
            
          }else if(input$icc_display == "cICC_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
            cicc_unique <- sort(unique(df_power[,"cicc"]))
            p <- vector(mode="list", length=length(cicc_unique))
            
            for(i in seq(length(cicc_unique))){
              
              if(i != 3){
                p[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(oicc=factor(oicc)) %>%
                  dplyr::filter(cicc == cicc_unique[i]) %>%
                  plot_ly(x=~n,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc,
                          linetype=~oicc, color=~oicc,legendgroup=~oicc, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("outcome ICC: ", oicc, "<br>covariate ICC: ", cicc,
                                       "<br>Clusters (n):", n, "<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)), height=input$dimension[2]*0.8) %>%
                  layout(title=paste("c-ICC = ", cicc_unique[i]),
                         xaxis=list(title="Number of clusters (n)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="outcome ICC")),
                         margin=0.01)
              }else{
                p[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(oicc=factor(oicc)) %>%
                  dplyr::filter(cicc == cicc_unique[i]) %>%
                  plot_ly(x=~n,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc,
                          linetype=~oicc, color=~oicc,legendgroup=~oicc,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("outcome ICC: ", oicc, "<br>covariate ICC: ", cicc,
                                       "<br>Clusters (n):", n, "<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)), height=input$dimension[2]*0.8) %>%
                  layout(title=paste("c-ICC = ", cicc_unique[i]),
                         xaxis=list(title="Number of clusters (n)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="outcome ICC")),
                         margin=0.01)
              }
              
            }
            
            subplot(p[[1]],p[[2]],p[[3]], nrows=2,
                    margin = 0.07, titleX=T, titleY=T) %>%
              layout(title = list(
                text='Number of clusters vs HTE power',
                font=list(size=17)
              ),
              #margin = list(t=50,b=50,pad=20),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.03,
                  text = paste0("<i>Minimum covariate ICC (", cicc_unique[1],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.43,
                  text = paste0("<i>Maximum covariate ICC (", cicc_unique[2],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.03,
                  text = paste0("<i>Assumed covariate ICC (", cicc_unique[3],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
            
          }
          
        }# end sensitivity if/else
        
      }else if(plot_display_react() == "fixed_power"){
        
        m_range <- seq(input$m_range[1],input$m_range[2])
        if(sensitivity_react() == "est_only"){
          df_n <- expand.grid(power=input$power, m=m_range,
                              oicc=#c(oicc_min(),
                                oicc_est(),
                              # oicc_max()),
                              cicc=#c(cicc_min(),
                                cicc_est(),
                              # cicc_max()),
                              var_y=var_y,#(input$sd_outcome)^2,
                              var_x=var_x,#(input$sd_covar)^2,
                              var_w=(input$w)*(1-input$w),
                              d=input$mean_diff_HTE, a=input$sig)
        }else{
          df_n <- expand.grid(power=input$power, m=m_range,
                              oicc=c(oicc_min(), oicc_est(), oicc_max()),
                              cicc=c(cicc_min(), cicc_est(), cicc_max()),
                              var_y=var_y,#(input$sd_outcome)^2,
                              var_x=var_x,#(input$sd_covar)^2,
                              var_w=(input$w)*(1-input$w),
                              d=input$mean_diff_HTE, a=input$sig)
        }
        
        
        
        n_hte_col <- matrix(NA, nrow(df_n), ncol=2)
        for(i in seq(nrow(df_n))){
          n_hte_col[i,] <- unlist( n_hte(oicc=df_n[i,"oicc"], cicc=df_n[i,"cicc"],
                                         m=df_n[i,"m"],power=df_n[i,"power"],
                                         var_y=df_n[i,"var_y"], var_x=df_n[i,"var_x"],
                                         var_w=df_n[i,"var_w"],d=df_n[i,"d"],
                                         a=df_n[i,"a"]) )
          # power_hte(oicc=df_power[i,"oicc"], cicc=df_power[i,"cicc"],
          #                             m=df_power[i,"m"],n=df_power[i,"n"],
          #                             var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
          #                             var_w=df_power[i,"var_w"],
          #                             d=df_power[i,"d"], a=df_power[i,"a"])
        }
        
        df_n <- cbind(df_n, n_hte_col)
        colnames(df_n)[(ncol(df_n)-1):ncol(df_n)] <- c("n","power_emp")
        max_n <- max(df_n[,"n"])
        
        if(sensitivity_react() == "est_only"){
          p_est <- df_n %>%
            as.data.frame() %>%
            dplyr::filter(oicc == oicc_est(),#oicc_unique[2],
                          cicc == cicc_est()#cicc_unique[2]
            ) %>%
            mutate(cicc=factor(cicc)) %>%
            plot_ly(x=~m,y=~n, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                    linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("outcome ICC: ", oicc, "<br>covariate ICC: ", cicc,
                                 "<br>Clusters (n):", n, "<br>Cluster size (m): ", m,
                                 "<br>HTE power: ", round(power_emp,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("o-ICC = ", oicc_est()),
                   xaxis=list(title="Cluster size (m)"),
                   yaxis=list(title="Number of clusters (n)",
                              range=list(0,max_n)),
                   legend=list(title=list(text="covariate ICC")),
                   margin=0.001)
          
          subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                  margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
            layout(title = list(
              text='Cluster size vs number of clusters',
              font=list(size=17)
            ),
            #margin=list(pad=50),
            annotations = list(
              list(
                x = 0.485,
                y = 1.03,
                text = paste0("<i>Assumed outcome ICC (", oicc_est(),") and covariate ICC (", cicc_est(),")</i>"),
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = FALSE, font=list(size=13)
              )
            ),
            legend=list(orientation="h",
                        yanchor="center",
                        y=0.25,
                        x=0.5)
            )
        }else if(sensitivity_react() == "sensitivity"){
          if(input$icc_display == "oICC_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_x$")
            
            oicc_unique <- sort(unique(df_n[,"oicc"]))
            cicc_unique <- sort(unique(df_n[,"cicc"]))
            p <- vector(mode="list", length=length(oicc_unique))
            
            for(i in seq(length(oicc_unique))){
              
              if(i != 3){
                p[[i]] <- df_n %>%
                  as.data.frame() %>%
                  mutate(cicc=factor(cicc)) %>%
                  dplyr::filter(oicc == oicc_unique[i]) %>%
                  plot_ly(x=~m,y=~n, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                          linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("outcome ICC: ", oicc, "<br>covariate ICC: ", cicc,
                                       "<br>Clusters (n):", n, "<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_emp,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC = ", oicc_unique[i]),
                         xaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="Number of clusters (n)",
                                    range=list(0,max_n)),
                         legend=list(title=list(text="covariate ICC")),
                         margin=0.01)
              }else{
                p[[i]] <- df_n %>%
                  as.data.frame() %>%
                  mutate(cicc=factor(cicc)) %>%
                  dplyr::filter(oicc == oicc_unique[i]) %>%
                  plot_ly(x=~m,y=~n, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                          linetype=~cicc, color=~cicc,legendgroup=~cicc,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("outcome ICC: ", oicc, "<br>covariate ICC: ", cicc,
                                       "<br>Clusters (n):", n, "<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_emp,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC = ", oicc_unique[i]),
                         xaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="Number of clusters (n)",
                                    range=list(0,max_n)),
                         legend=list(title=list(text="covariate ICC")),
                         margin=0.01)
              }
              
            }
            
            subplot(p[[1]],p[[2]],p[[3]], nrows=2,
                    margin = 0.07, titleX=T, titleY=T) %>%
              layout(title = list(
                text='Cluster size vs number of clusters',
                font=list(size=17)
              ),
              #margin = list(t=50,b=50,pad=20),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.03,
                  text = paste0("<i>Minimum outcome ICC (", oicc_unique[1], ")"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.43,
                  text = paste0("<i>Maximum outcome ICC (", oicc_unique[2], ")"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.03,
                  text = paste0("<i>Assumed outcome ICC (", oicc_unique[3],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
            
          }else if(input$icc_display == "cICC_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
            cicc_unique <- sort(unique(df_n[,"cicc"]))
            p <- vector(mode="list", length=length(cicc_unique))
            
            for(i in seq(length(cicc_unique))){
              
              if(i != 3){
                p[[i]] <- df_n %>%
                  as.data.frame() %>%
                  mutate(oicc=factor(oicc)) %>%
                  dplyr::filter(cicc == cicc_unique[i]) %>%
                  plot_ly(x=~m,y=~n, type='scatter', mode='lines', line=list(width=3), name=~oicc,
                          linetype=~oicc, color=~oicc,legendgroup=~oicc, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("outcome ICC: ", oicc, "<br>covariate ICC: ", cicc,
                                       "<br>Clusters (n):", n, "<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_emp,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("c-ICC = ", cicc_unique[i]),
                         xaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="Number of clusters (n)",
                                    range=list(0,max_n)),
                         legend=list(title=list(text="outcome ICC")),
                         margin=0.01)
              }else{
                p[[i]] <- df_n %>%
                  as.data.frame() %>%
                  mutate(oicc=factor(oicc)) %>%
                  dplyr::filter(cicc == cicc_unique[i]) %>%
                  plot_ly(x=~m,y=~n, type='scatter', mode='lines', line=list(width=3), name=~oicc,
                          linetype=~oicc, color=~oicc,legendgroup=~oicc,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("outcome ICC: ", oicc, "<br>covariate ICC: ", cicc,
                                       "<br>Clusters (n):", n, "<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_emp,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("c-ICC = ", cicc_unique[i]),
                         xaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="Number of clusters (n)",
                                    range=list(0,max_n)),
                         legend=list(title=list(text="outcome ICC")),
                         margin=0.01)
              }
              
            }
            
            subplot(p[[1]],p[[2]],p[[3]], nrows=2,
                    margin = 0.07, titleX=T, titleY=T) %>%
              layout(title = list(
                text='Cluster size vs number of clusters',
                font=list(size=17)
              ),
              #margin = list(t=50,b=50,pad=20),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.03,
                  text = paste0("<i>Minimum covariate ICC (", cicc_unique[1],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.43,
                  text = paste0("<i>Maximum covariate ICC (", cicc_unique[2],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.03,
                  text = paste0("<i>Assumed covariate ICC (", cicc_unique[3],")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
            
          }# end constant ICC if/else
        }# end sensitivity if/else
        
        
      }# end plot display if/else
      
      #### THREE LEVEL ####
    }else if(trial_react() == "three_level"){
      
      # determine effect sizes and variances depending on outcome/covariate type #
      if(outcome_type_react() =="continuous"){
        var_y <- (sd_outcome_react())^2
        d <- input$mean_diff_HTE
        
        if(covar_type_react()=="continuous"){
          # continuous outcome and covar power
          var_x <- (sd_covar_react())^2
          
        }else if(covar_type_react() == "binary"){
          # continuous outcome binary covar power
          var_x <- (prop_covar_react())*(1-prop_covar_react())
        }
        
      }else if(outcome_type_react() == "binary"){
        var_y <- ((prop_control_react()*(1-prop_control_react())) + (prop_trt_react()*(1-prop_trt_react())))/2
        #d <- (input$prop_trt - input$prop_control)
        
        if(covar_type_react() == "continuous"){
          # binary outcome continuous covar power
          var_x <-(sd_covar_react())^2
          
        }else if(covar_type_react() == "binary"){
          # binary outcome and covar power #
          var_x <- (prop_covar_react())*(1-prop_covar_react())
          
        }
      }
      
      if(plot_display_react() == "m_v_power"){
        m_range <- seq(input$m_three_range[1],input$m_three_range[2])
        if(sensitivity_three_react() == "est_only"){
          df_power <- expand.grid(nc=input$nc, ns=input$ns, m=m_range,
                                  a0=#c(oicc_wsub_min_three(),
                                    oicc_wsub_est_three(),
                                  # oicc_wsub_max_three()),
                                  a1_a0=#c(oicc_ratio_min_three(),
                                    oicc_ratio_est_three(),
                                  # oicc_ratio_max_three()),
                                  r0=#c(cicc_wsub_min_three(),
                                    cicc_wsub_est_three(),
                                  # cicc_wsub_max_three()),
                                  r1_r0=#c(cicc_ratio_min_three(),
                                    cicc_ratio_est_three(),
                                  # cicc_ratio_max_three()),
                                  var_y=var_y, #(input$sd_outcome)^2,
                                  var_x=var_x,#(input$sd_covar)^2,
                                  pw=input$w,
                                  d=input$mean_diff_HTE, a=input$sig,
                                  rand=input$randomization_three)
        }else{
          df_power <- expand.grid(nc=input$nc, ns=input$ns, m=m_range,
                                  a0=c(oicc_wsub_min_three(),
                                       oicc_wsub_est_three(),
                                       oicc_wsub_max_three()),
                                  a1_a0=c(oicc_ratio_min_three(),
                                          oicc_ratio_est_three(),
                                          oicc_ratio_max_three()),
                                  r0=c(cicc_wsub_min_three(),
                                       cicc_wsub_est_three(),
                                       cicc_wsub_max_three()),
                                  r1_r0=c(cicc_ratio_min_three(),
                                          cicc_ratio_est_three(),
                                          cicc_ratio_max_three()),
                                  var_y=var_y, #(input$sd_outcome)^2,
                                  var_x=var_x,#(input$sd_covar)^2,
                                  pw=input$w,
                                  d=input$mean_diff_HTE, a=input$sig,
                                  rand=input$randomization_three)
        }
        
        
        power_hte_col <- rep(NA, nrow(df_power))
        for(i in seq(nrow(df_power))){
          power_hte_col[i] <- power_hte_three(nc=df_power[i,"nc"], ns=df_power[i,"ns"],
                                              m=df_power[i,"m"],
                                              pw=df_power[i,"pw"],
                                              var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                              a0=df_power[i,"a0"], a1_a0=df_power[i,"a1_a0"],
                                              r0=df_power[i,"r0"], r1_r0=df_power[i,"r1_r0"],
                                              rand=df_power[i,"rand"],
                                              d=df_power[i,"d"], a=df_power[i,"a"])
        }
        
        df_power <- cbind(df_power, power_hte_col)
        
        if(sensitivity_three_react() == "est_only"){
          p_est <- df_power %>%
            as.data.frame() %>%
            dplyr::filter(r0 == cicc_wsub_est_three(),
                          r1_r0 == cicc_ratio_est_three(),
                          a0 == oicc_wsub_est_three(),
                          a1_a0 == oicc_ratio_est_three()) %>%
            mutate(r1_r0=factor(r1_r0)) %>%
            plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                    linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("within-subcluster outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                 "; <br>within-subcluster covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                 "<br>Clusters (nc):", nc, "; Subclusters (ns):", ns,"<br>Cluster size (m): ", m,
                                 "<br>HTE power: ", round(power_hte_col,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("within-subcluster covariate ICC = ", cicc_wsub_est_three()),
                   #subtitle=paste("within-subcluster outcome ICC = ", oicc_wsub_est_three(), ", outcome ICC ratio = ", oicc_ratio_est_three()),
                   xaxis=list(title="Cluster size (m)"),
                   yaxis=list(title="HTE Power"),
                   legend=list(title=list(text="covariate ICC ratio")),
                   margin=0.001)
          
          subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                  margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
            layout(title = list(
              text='Cluster size vs HTE power',
              font=list(size=17)
            ),
            #margin=list(pad=50),
            annotations = list(
              list(
                x = 0.485,
                y = 1.03,#0.97,
                text = paste0("<i>Assumed a0 (", oicc_wsub_est_three(), "), a1/a0 (", oicc_ratio_est_three(),"), r0 (", cicc_wsub_est_three(),") and r1/r0 (", cicc_ratio_est_three(),")</i>"),
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = FALSE, font=list(size=13)
              )
            ),
            legend=list(orientation="h",
                        yanchor="center",
                        y=0.25,
                        x=0.5)
            )
        }else if(sensitivity_three_react() == "sensitivity"){
          if(input$icc_display_three == "oICC_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_x$")
            
            a0_unique <- sort(unique(df_power[,"a0"]))
            a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
            r0_unique <- sort(unique(df_power[,"r0"]))
            r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
            p_o <- p_c <- vector(mode="list", length=length(a0_unique))
            
            for(i in seq(length(r0_unique))){
              
              if(i != 3){
                # outcome ICCs #
                p_o[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(r1_r0=factor(r1_r0)) %>%
                  dplyr::filter(a0 == oicc_wsub_est_three(),
                                a1_a0 == oicc_ratio_est_three(),
                                r0 == r0_unique[i]) %>%
                  plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                          linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("within-subcluster outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                       "; <br>within-subcluster covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                       "<br>Clusters (nc):", nc, "; Subclusters (ns):", ns,"<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("within-subcluster covariate ICC = ", r0_unique[i]),
                         #subtitle=paste("within-subcluster outcome ICC = ", oicc_wsub_est_three(), ", outcome ICC ratio = ", oicc_ratio_est_three()),
                         xaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="covariate CAC")),
                         margin=0.01)
                
              }else{
                p_o[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(r1_r0=factor(r1_r0)) %>%
                  dplyr::filter(a0 == oicc_wsub_est_three(),
                                a1_a0 == oicc_ratio_est_three(),
                                r0 == r0_unique[i]) %>%
                  plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                          linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("within-subcluster outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                       "; <br>within-subcluster covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                       "<br>Clusters (nc):", nc, "; Subclusters (ns):", ns,"<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("within-subcluster covariate ICC = ", r0_unique[i]),
                         #subtitle=paste("within-subcluster outcome ICC = ", oicc_wsub_est_three(), ", outcome ICC ratio = ", oicc_ratio_est_three()),
                         xaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="covariate CAC")),
                         margin=0.01)
              }
              
            }
            
            subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                    margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size vs HTE power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.03,
                  #yshift=-30,
                  text = paste0("<i>Assumed a0 (", oicc_wsub_est_three(), ") and a1/a0 (", oicc_ratio_est_three(),"), minimum r0 (", cicc_wsub_min_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE#,
                  #font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.43,
                  text = paste0("<i>Assumed a0 (", oicc_wsub_est_three(), ") and a1/a0 (", oicc_ratio_est_three(),"), and maximum r0 (", cicc_wsub_max_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.03,
                  text = paste0("<i>Assumed a0 (", oicc_wsub_est_three(), "), a1/a0 (", oicc_ratio_est_three(),"), and r0 (", cicc_wsub_est_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
            
          }else if(input$icc_display_three == "cICC_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
            a0_unique <- sort(unique(df_power[,"a0"]))
            a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
            r0_unique <- sort(unique(df_power[,"r0"]))
            r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
            p_o <- p_c <- vector(mode="list", length=length(a0_unique))
            
            for(i in seq(length(a0_unique))){
              
              if(i != 3){
                # outcome ICCs #
                p_c[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(a1_a0=factor(a1_a0)) %>%
                  dplyr::filter(r0 == cicc_wsub_est_three(),
                                r1_r0 == cicc_ratio_est_three(),
                                a0 == a0_unique[i]) %>%
                  plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                          linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("within-subcluster outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                       "; <br>within-subcluster covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                       "<br>Clusters (nc):", nc, "; Subclusters (ns):", ns,"<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("within-subcluster outcome ICC = ", a0_unique[i]),
                         #subtitle=paste("within-subcluster covariate ICC = ", cicc_wsub_est_three(), ", covariate ICC ratio = ", cicc_ratio_est_three()),
                         xaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="outcome CAC")),
                         margin=0.01)
                
              }else{
                p_c[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(a1_a0=factor(a1_a0)) %>%
                  dplyr::filter(r0 == cicc_wsub_est_three(),
                                r1_r0 == cicc_ratio_est_three(),
                                a0 == a0_unique[i]) %>%
                  plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                          linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("within-subcluster outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                       "; <br>within-subcluster covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                       "<br>Clusters (nc):", nc, "; Subclusters (ns):", ns,"<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("within-subcluster outcome ICC = ", a0_unique[i]),
                         #subtitle=paste("within-subcluster covariate ICC = ", cicc_wsub_est_three(), ", covariate ICC ratio = ", cicc_ratio_est_three()),
                         xaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="outcome CAC")),
                         margin=0.01)
              }
              
            }
            
            subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                    margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size vs HTE power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.03,
                  #yshift=-30,
                  text = paste0("<i>Assumed r0 (", cicc_wsub_est_three(), ") and r1/r0 (", cicc_ratio_est_three(),"), minimum a0 (", oicc_wsub_min_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE#,
                  #font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.43,
                  text = paste0("<i>Assumed r0 (", cicc_wsub_est_three(), ") and r1/r0 (", cicc_ratio_est_three(),"), and maximum a0 (", oicc_wsub_max_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.03,
                  text = paste0("<i>Assumed r0 (", cicc_wsub_est_three(), "), r1/r0 (", cicc_ratio_est_three(),"), and a0 (", oicc_wsub_est_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
            
          }
          
        }# end sensitivity if/else
        
        
      }else if(plot_display_react() == "nc_v_power"){
        
        nc_range <- seq(input$nc_range[1],input$nc_range[2])
        if(sensitivity_three_react() == "est_only"){
          df_power <- expand.grid(nc=nc_range, ns=input$ns, m=input$m_three,
                                  a0=#c(oicc_wsub_min_three(),
                                    oicc_wsub_est_three(),
                                  # oicc_wsub_max_three()),
                                  a1_a0=#c(oicc_ratio_min_three(),
                                    oicc_ratio_est_three(),
                                  # oicc_ratio_max_three()),
                                  r0=#c(cicc_wsub_min_three(),
                                    cicc_wsub_est_three(),
                                  # cicc_wsub_max_three()),
                                  r1_r0=#c(cicc_ratio_min_three(),
                                    cicc_ratio_est_three(),
                                  # cicc_ratio_max_three()),
                                  var_y=var_y, #(input$sd_outcome)^2,
                                  var_x=var_x,#(input$sd_covar)^2,
                                  pw=input$w,
                                  d=input$mean_diff_HTE, a=input$sig,
                                  rand=input$randomization_three)
        }else{
          df_power <- expand.grid(nc=nc_range, ns=input$ns, m=input$m_three,
                                  a0=c(oicc_wsub_min_three(),
                                       oicc_wsub_est_three(),
                                       oicc_wsub_max_three()),
                                  a1_a0=c(oicc_ratio_min_three(),
                                          oicc_ratio_est_three(),
                                          oicc_ratio_max_three()),
                                  r0=c(cicc_wsub_min_three(),
                                       cicc_wsub_est_three(),
                                       cicc_wsub_max_three()),
                                  r1_r0=c(cicc_ratio_min_three(),
                                          cicc_ratio_est_three(),
                                          cicc_ratio_max_three()),
                                  var_y=var_y, #(input$sd_outcome)^2,
                                  var_x=var_x,#(input$sd_covar)^2,
                                  pw=input$w,
                                  d=input$mean_diff_HTE, a=input$sig,
                                  rand=input$randomization_three)
        }
        
        
        power_hte_col <- rep(NA, nrow(df_power))
        for(i in seq(nrow(df_power))){
          power_hte_col[i] <- power_hte_three(nc=df_power[i,"nc"], ns=df_power[i,"ns"],
                                              m=df_power[i,"m"],
                                              pw=df_power[i,"pw"],
                                              var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                              a0=df_power[i,"a0"], a1_a0=df_power[i,"a1_a0"],
                                              r0=df_power[i,"r0"], r1_r0=df_power[i,"r1_r0"],
                                              rand=df_power[i,"rand"],
                                              d=df_power[i,"d"], a=df_power[i,"a"])
        }
        
        df_power <- cbind(df_power, power_hte_col)
        
        if(sensitivity_three_react() == "est_only"){
          p_est <- df_power %>%
            as.data.frame() %>%
            dplyr::filter(r0 == cicc_wsub_est_three(),
                          r1_r0 == cicc_ratio_est_three(),
                          a0 == oicc_wsub_est_three(),
                          a1_a0 == oicc_ratio_est_three()) %>%
            mutate(r1_r0=factor(r1_r0)) %>%
            plot_ly(x=~nc,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                    linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("within-subcluster outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                 "; <br>within-subcluster covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                 "<br>Clusters (nc):", nc, "; Subclusters (ns):", ns,"<br>Cluster size (m): ", m,
                                 "<br>HTE power: ", round(power_hte_col,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("within-subcluster covariate ICC = ", cicc_wsub_est_three()),
                   #subtitle=paste("within-subcluster outcome ICC = ", oicc_wsub_est_three(), ", outcome ICC ratio = ", oicc_ratio_est_three()),
                   xaxis=list(title="Number of clusters (nc)"),
                   yaxis=list(title="HTE Power"),
                   legend=list(title=list(text="covariate ICC ratio")),
                   margin=0.001)
          
          subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                  margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
            layout(title = list(
              text='Number of clusters vs HTE Power',
              font=list(size=17)
            ),
            #margin=list(pad=50),
            annotations = list(
              list(
                x = 0.485,
                y = 1.03,#0.97,
                text = paste0("<i>Assumed a0 (", oicc_wsub_est_three(), "), a1/a0 (", oicc_ratio_est_three(),"), r0 (", cicc_wsub_est_three(),") and r1/r0 (", cicc_ratio_est_three(),")</i>"),
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = FALSE, font=list(size=13)
              )
            ),
            legend=list(orientation="h",
                        yanchor="center",
                        y=0.25,
                        x=0.5)
            )
        }else if(sensitivity_three_react() == "sensitivity"){
          if(input$icc_display_three == "oICC_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_x$")
            
            a0_unique <- sort(unique(df_power[,"a0"]))
            a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
            r0_unique <- sort(unique(df_power[,"r0"]))
            r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
            p_o <- p_c <- vector(mode="list", length=length(a0_unique))
            
            for(i in seq(length(r0_unique))){
              
              if(i != 3){
                # outcome ICCs #
                p_o[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(r1_r0=factor(r1_r0)) %>%
                  dplyr::filter(a0 == oicc_wsub_est_three(),
                                a1_a0 == oicc_ratio_est_three(),
                                r0 == r0_unique[i]) %>%
                  plot_ly(x=~nc,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                          linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("within-subcluster outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                       "; <br>within-subcluster covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                       "<br>Clusters (nc):", nc, "; Subclusters (ns):", ns,"<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("within-subcluster covariate ICC = ", r0_unique[i]),
                         #subtitle=paste("within-subcluster outcome ICC = ", oicc_wsub_est_three(), ", outcome ICC ratio = ", oicc_ratio_est_three()),
                         xaxis=list(title="Number of clusters (nc)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="covariate CAC")),
                         margin=0.01)
                
              }else{
                p_o[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(r1_r0=factor(r1_r0)) %>%
                  dplyr::filter(a0 == oicc_wsub_est_three(),
                                a1_a0 == oicc_ratio_est_three(),
                                r0 == r0_unique[i]) %>%
                  plot_ly(x=~nc,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                          linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("within-subcluster outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                       "; <br>within-subcluster covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                       "<br>Clusters (nc):", nc, "; Subclusters (ns):", ns,"<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("within-subcluster covariate ICC = ", r0_unique[i]),
                         #subtitle=paste("within-subcluster outcome ICC = ", oicc_wsub_est_three(), ", outcome ICC ratio = ", oicc_ratio_est_three()),
                         xaxis=list(title="Number of clusters (nc)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="covariate CAC")),
                         margin=0.01)
              }
              
            }
            
            subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                    margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
              layout(title = list(
                text='Number of clusters vs HTE Power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.03,
                  #yshift=-30,
                  text = paste0("<i>Assumed a0 (", oicc_wsub_est_three(), ") and a1/a0 (", oicc_ratio_est_three(),"), minimum r0 (", cicc_wsub_min_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE#,
                  #font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.43,
                  text = paste0("<i>Assumed a0 (", oicc_wsub_est_three(), ") and a1/a0 (", oicc_ratio_est_three(),"), and maximum r0 (", cicc_wsub_max_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.03,
                  text = paste0("<i>Assumed a0 (", oicc_wsub_est_three(), "), a1/a0 (", oicc_ratio_est_three(),"), and r0 (", cicc_wsub_est_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
            
          }else if(input$icc_display_three == "cICC_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
            a0_unique <- sort(unique(df_power[,"a0"]))
            a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
            r0_unique <- sort(unique(df_power[,"r0"]))
            r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
            p_o <- p_c <- vector(mode="list", length=length(a0_unique))
            
            for(i in seq(length(a0_unique))){
              
              if(i != 3){
                # outcome ICCs #
                p_c[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(a1_a0=factor(a1_a0)) %>%
                  dplyr::filter(r0 == cicc_wsub_est_three(),
                                r1_r0 == cicc_ratio_est_three(),
                                a0 == a0_unique[i]) %>%
                  plot_ly(x=~nc,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                          linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("within-subcluster outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                       "; <br>within-subcluster covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                       "<br>Clusters (nc):", nc, "; Subclusters (ns):", ns,"<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("within-subcluster outcome ICC = ", a0_unique[i]),
                         #subtitle=paste("within-subcluster covariate ICC = ", cicc_wsub_est_three(), ", covariate ICC ratio = ", cicc_ratio_est_three()),
                         xaxis=list(title="Number of clusters (nc)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="outcome CAC")),
                         margin=0.01)
                
              }else{
                p_c[[i]] <- df_power %>%
                  as.data.frame() %>%
                  mutate(a1_a0=factor(a1_a0)) %>%
                  dplyr::filter(r0 == cicc_wsub_est_three(),
                                r1_r0 == cicc_ratio_est_three(),
                                a0 == a0_unique[i]) %>%
                  plot_ly(x=~nc,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                          linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("within-subcluster outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                       "; <br>within-subcluster covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                       "<br>Clusters (nc):", nc, "; Subclusters (ns):", ns,"<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("within-subcluster outcome ICC = ", a0_unique[i]),
                         #subtitle=paste("within-subcluster covariate ICC = ", cicc_wsub_est_three(), ", covariate ICC ratio = ", cicc_ratio_est_three()),
                         xaxis=list(title="Number of clusters (nc)"),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="outcome CAC")),
                         margin=0.01)
              }
              
            }
            
            subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                    margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
              layout(title = list(
                text='Number of clusters vs HTE Power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.03,
                  #yshift=-30,
                  text = paste0("<i>Assumed r0 (", cicc_wsub_est_three(), ") and r1/r0 (", cicc_ratio_est_three(),"), minimum a0 (", oicc_wsub_min_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE#,
                  #font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.43,
                  text = paste0("<i>Assumed r0 (", cicc_wsub_est_three(), ") and r1/r0 (", cicc_ratio_est_three(),"), and maximum a0 (", oicc_wsub_max_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.03,
                  text = paste0("<i>Assumed r0 (", cicc_wsub_est_three(), "), r1/r0 (", cicc_ratio_est_three(),"), and a0 (", oicc_wsub_est_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
            
          }
          
        }# end sensitivity if/else
        
      }else if(plot_display_react() == "m_v_nc"){
        
        
        m_range <- seq(input$m_three_range[1],input$m_three_range[2])
        
        if(sensitivity_three_react() == "est_only"){
          df_nc <- expand.grid(power=input$power_three, m=m_range,
                               ns=input$ns,
                               a0=#c(oicc_wsub_min_three(),
                                 oicc_wsub_est_three(),
                               # oicc_wsub_max_three()),
                               a1_a0=#c(oicc_ratio_min_three(),
                                 oicc_ratio_est_three(),
                               # oicc_ratio_max_three()),
                               r0=#c(cicc_wsub_min_three(),
                                 cicc_wsub_est_three(),
                               # cicc_wsub_max_three()),
                               r1_r0=#c(cicc_ratio_min_three(),
                                 cicc_ratio_est_three(),
                               # cicc_ratio_max_three()),
                               var_y=var_y, #(input$sd_outcome)^2,
                               var_x=var_x,#(input$sd_covar)^2,
                               pw=input$w,
                               d=input$mean_diff_HTE, a=input$sig,
                               rand=input$randomization_three)
        }else{
          df_nc <- expand.grid(power=input$power_three, m=m_range,
                               ns=input$ns,
                               a0=c(oicc_wsub_min_three(),
                                    oicc_wsub_est_three(),
                                    oicc_wsub_max_three()),
                               a1_a0=c(oicc_ratio_min_three(),
                                       oicc_ratio_est_three(),
                                       oicc_ratio_max_three()),
                               r0=c(cicc_wsub_min_three(),
                                    cicc_wsub_est_three(),
                                    cicc_wsub_max_three()),
                               r1_r0=c(cicc_ratio_min_three(),
                                       cicc_ratio_est_three(),
                                       cicc_ratio_max_three()),
                               var_y=var_y, #(input$sd_outcome)^2,
                               var_x=var_x,#(input$sd_covar)^2,
                               pw=input$w,
                               d=input$mean_diff_HTE, a=input$sig,
                               rand=input$randomization_three)
        }
        
        
        nc_hte_col <- matrix(NA, nrow(df_nc), ncol=2)
        for(i in seq(nrow(df_nc))){
          nc_hte_col[i,] <- unlist(nc_hte_three(m=df_nc[i,"m"],
                                                ns=df_nc[i,"ns"],
                                                power=df_nc[i,"power"],
                                                pw=df_nc[i,"pw"],
                                                var_y=df_nc[i,"var_y"], var_x=df_nc[i,"var_x"],
                                                a0=df_nc[i,"a0"], a1_a0=df_nc[i,"a1_a0"],
                                                r0=df_nc[i,"r0"], r1_r0=df_nc[i,"r1_r0"],
                                                rand=df_nc[i,"rand"],
                                                d=df_nc[i,"d"], a=df_nc[i,"a"]))
        }
        
        df_nc <- cbind(df_nc, nc_hte_col)
        colnames(df_nc)[(ncol(df_nc)-1):ncol(df_nc)] <- c("nc","power_emp")
        
        if(sensitivity_three_react() == "est_only"){
          p_est <- df_nc %>%
            as.data.frame() %>%
            dplyr::filter(r0 == cicc_wsub_est_three(),
                          r1_r0 == cicc_ratio_est_three(),
                          a0 == oicc_wsub_est_three(),
                          a1_a0 == oicc_ratio_est_three()) %>%
            mutate(r1_r0=factor(r1_r0)) %>%
            plot_ly(x=~m,y=~nc, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                    linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("within-subcluster outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                 "; <br>within-subcluster covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                 "<br>Clusters (nc):", nc, "; Subclusters (ns):", ns,"<br>Cluster size (m): ", m,
                                 "<br>HTE power: ", round(power_emp,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("within-subcluster covariate ICC = ", cicc_wsub_est_three()),
                   #subtitle=paste("within-subcluster outcome ICC = ", oicc_wsub_est_three(), ", outcome ICC ratio = ", oicc_ratio_est_three()),
                   xaxis=list(title="Cluster size (m)"),
                   yaxis=list(title="Number of clusters (nc)"),
                   legend=list(title=list(text="covariate ICC ratio")),
                   margin=0.001)
          
          subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                  margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
            layout(title = list(
              text='Cluster size vs number of clusters',
              font=list(size=17)
            ),
            #margin=list(pad=50),
            annotations = list(
              list(
                x = 0.485,
                y = 1.03,#0.97,
                text = paste0("<i>Assumed a0 (", oicc_wsub_est_three(), "), a1/a0 (", oicc_ratio_est_three(),"), r0 (", cicc_wsub_est_three(),") and r1/r0 (", cicc_ratio_est_three(),")</i>"),
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = FALSE, font=list(size=13)
              )
            ),
            legend=list(orientation="h",
                        yanchor="center",
                        y=0.25,
                        x=0.5)
            )
        }else if(sensitivity_three_react() == "sensitivity"){
          if(input$icc_display_three == "oICC_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_x$")
            
            a0_unique <- sort(unique(df_nc[,"a0"]))
            a1_a0_unique <- sort(unique(df_nc[,"a1_a0"]))
            r0_unique <- sort(unique(df_nc[,"r0"]))
            r1_r0_unique <- sort(unique(df_nc[,"r1_r0"]))
            p_o <- p_c <- vector(mode="list", length=length(a0_unique))
            
            for(i in seq(length(r0_unique))){
              
              if(i != 3){
                # outcome ICCs #
                p_o[[i]] <- df_nc %>%
                  as.data.frame() %>%
                  mutate(r1_r0=factor(r1_r0)) %>%
                  dplyr::filter(a0 == oicc_wsub_est_three(),
                                a1_a0 == oicc_ratio_est_three(),
                                r0 == r0_unique[i]) %>%
                  plot_ly(x=~m,y=~nc, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                          linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("within-subcluster outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                       "; <br>within-subcluster covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                       "<br>Clusters (nc):", nc, "; Subclusters (ns):", ns,"<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_emp,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("within-subcluster covariate ICC = ", r0_unique[i]),
                         #subtitle=paste("within-subcluster outcome ICC = ", oicc_wsub_est_three(), ", outcome ICC ratio = ", oicc_ratio_est_three()),
                         xaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="Number of clusters (nc)"),
                         legend=list(title=list(text="covariate CAC")),
                         margin=0.01)
                
              }else{
                p_o[[i]] <- df_nc %>%
                  as.data.frame() %>%
                  mutate(r1_r0=factor(r1_r0)) %>%
                  dplyr::filter(a0 == oicc_wsub_est_three(),
                                a1_a0 == oicc_ratio_est_three(),
                                r0 == r0_unique[i]) %>%
                  plot_ly(x=~m,y=~nc, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                          linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("within-subcluster outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                       "; <br>within-subcluster covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                       "<br>Clusters (nc):", nc, "; Subclusters (ns):", ns,"<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_emp,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("within-subcluster covariate ICC = ", r0_unique[i]),
                         #subtitle=paste("within-subcluster outcome ICC = ", oicc_wsub_est_three(), ", outcome ICC ratio = ", oicc_ratio_est_three()),
                         xaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="Number of clusters (nc)"),
                         legend=list(title=list(text="covariate CAC")),
                         margin=0.01)
              }
              
            }
            
            subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                    margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size vs number of clusters',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.03,
                  #yshift=-30,
                  text = paste0("<i>Assumed a0 (", oicc_wsub_est_three(), ") and a1/a0 (", oicc_ratio_est_three(),"), minimum r0 (", cicc_wsub_min_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE#,
                  #font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.43,
                  text = paste0("<i>Assumed a0 (", oicc_wsub_est_three(), ") and a1/a0 (", oicc_ratio_est_three(),"), and maximum r0 (", cicc_wsub_max_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.03,
                  text = paste0("<i>Assumed a0 (", oicc_wsub_est_three(), "), a1/a0 (", oicc_ratio_est_three(),"), and r0 (", cicc_wsub_est_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
            
          }else if(input$icc_display_three == "cICC_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
            a0_unique <- sort(unique(df_nc[,"a0"]))
            a1_a0_unique <- sort(unique(df_nc[,"a1_a0"]))
            r0_unique <- sort(unique(df_nc[,"r0"]))
            r1_r0_unique <- sort(unique(df_nc[,"r1_r0"]))
            p_o <- p_c <- vector(mode="list", length=length(a0_unique))
            
            for(i in seq(length(a0_unique))){
              
              if(i != 3){
                # outcome ICCs #
                p_c[[i]] <- df_nc %>%
                  as.data.frame() %>%
                  mutate(a1_a0=factor(a1_a0)) %>%
                  dplyr::filter(r0 == cicc_wsub_est_three(),
                                r1_r0 == cicc_ratio_est_three(),
                                a0 == a0_unique[i]) %>%
                  plot_ly(x=~m,y=~nc, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                          linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("within-subcluster outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                       "; <br>within-subcluster covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                       "<br>Clusters (nc):", nc, "; Subclusters (ns):", ns,"<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_emp,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("within-subcluster outcome ICC = ", a0_unique[i]),
                         #subtitle=paste("within-subcluster covariate ICC = ", cicc_wsub_est_three(), ", covariate ICC ratio = ", cicc_ratio_est_three()),
                         yaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="Number of clusters (nc)"),
                         legend=list(title=list(text="outcome CAC")),
                         margin=0.01)
                
              }else{
                p_c[[i]] <- df_nc %>%
                  as.data.frame() %>%
                  mutate(a1_a0=factor(a1_a0)) %>%
                  dplyr::filter(r0 == cicc_wsub_est_three(),
                                r1_r0 == cicc_ratio_est_three(),
                                a0 == a0_unique[i]) %>%
                  plot_ly(x=~m,y=~nc, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                          linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("within-subcluster outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                       "; <br>within-subcluster covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                       "<br>Clusters (nc):", nc, "; Subclusters (ns):", ns,"<br>Cluster size (m): ", m,
                                       "<br>HTE power: ", round(power_emp,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("within-subcluster outcome ICC = ", a0_unique[i]),
                         #subtitle=paste("within-subcluster covariate ICC = ", cicc_wsub_est_three(), ", covariate ICC ratio = ", cicc_ratio_est_three()),
                         yaxis=list(title="Cluster size (m)"),
                         yaxis=list(title="Number of clusters (nc)"),
                         legend=list(title=list(text="outcome CAC")),
                         margin=0.01)
              }
              
            }
            
            subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                    margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size vs number of clusters',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.03,
                  #yshift=-30,
                  text = paste0("<i>Assumed r0 (", cicc_wsub_est_three(), ") and r1/r0 (", cicc_ratio_est_three(),"), minimum a0 (", oicc_wsub_min_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE#,
                  #font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.43,
                  text = paste0("<i>Assumed r0 (", cicc_wsub_est_three(), ") and r1/r0 (", cicc_ratio_est_three(),"), and maximum a0 (", oicc_wsub_max_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.03,
                  text = paste0("<i>Assumed r0 (", cicc_wsub_est_three(), "), r1/r0 (", cicc_ratio_est_three(),"), and a0 (", oicc_wsub_est_three(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
            
          }
        }# end sensitivity if/else
        
        
      }# end plot display if/else
      
      #### SWD ####
    }else if(trial_react() == "SWD"){
      ## NEED TO DO  closed- cohort option ##
      # determine effect sizes and variances depending on outcome/covariate type #
      if(outcome_type_react() =="continuous"){
        var_y <- (sd_outcome_react())^2
        d <- input$mean_diff_HTE
        
        if(covar_type_react()=="continuous"){
          # continuous outcome and covar power
          var_x <- (sd_covar_react())^2
          
        }else if(covar_type_react() == "binary"){
          # continuous outcome binary covar power
          var_x <- (prop_covar_react())*(1-prop_covar_react())
        }
        
      }else if(outcome_type_react() == "binary"){
        var_y <- ((prop_control_react()*(1-prop_control_react())) + (prop_trt_react()*(1-prop_trt_react())))/2
        #d <- (input$prop_trt - input$prop_control)
        
        if(covar_type_react() == "continuous"){
          # binary outcome continuous covar power
          var_x <-(sd_covar_react())^2
          
        }else if(covar_type_react() == "binary"){
          # binary outcome and covar power #
          var_x <- (prop_covar_react())*(1-prop_covar_react())
          
        }
      }
      
      if(plot_display_react() == "m_v_power"){
        m_range <- seq(input$m_swd_range[1],input$m_swd_range[2])
        J <- input$J_1 + 1
        n <- input$ns_swd*input$J_1
        
        if(cohort()=="cross"){
          
          if(sensitivity_swd_react() == "est_only"){
            df_power <- expand.grid(n=n, m=m_range, J=J,
                                    a0=#c(oicc_wperiod_min_swd(),
                                      oicc_wperiod_est_swd(),
                                    #oicc_wperiod_max_swd()),
                                    a1_a0=#c(oicc_ratio_min_swd(),
                                      oicc_ratio_est_swd(),
                                    # oicc_ratio_max_swd()),
                                    r0=#c(cicc_wperiod_min_swd(),
                                      cicc_wperiod_est_swd(),
                                    # cicc_wperiod_max_swd()),
                                    r1_r0=#c(cicc_ratio_min_swd(),
                                      cicc_ratio_est_swd(),
                                    # cicc_ratio_max_swd()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
          }else{
            df_power <- expand.grid(n=n, m=m_range, J=J,
                                    a0=c(oicc_wperiod_min_swd(),
                                         oicc_wperiod_est_swd(),
                                         oicc_wperiod_max_swd()),
                                    a1_a0=c(oicc_ratio_min_swd(),
                                            oicc_ratio_est_swd(),
                                            oicc_ratio_max_swd()),
                                    r0=c(cicc_wperiod_min_swd(),
                                         cicc_wperiod_est_swd(),
                                         cicc_wperiod_max_swd()),
                                    r1_r0=c(cicc_ratio_min_swd(),
                                            cicc_ratio_est_swd(),
                                            cicc_ratio_max_swd()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
          }
          
          
          power_hte_col <- rep(NA, nrow(df_power))
          for(i in seq(nrow(df_power))){
            power_hte_col[i] <- power_hte_swd(n=df_power[i,"n"],
                                              m=df_power[i,"m"], J=df_power[i,"J"],
                                              var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                              a0=df_power[i,"a0"], a1_a0=df_power[i,"a1_a0"],
                                              r0=df_power[i,"r0"], r1_r0=df_power[i,"r1_r0"],
                                              cohort=df_power[i,"cohort"],
                                              d=df_power[i,"d"], a=df_power[i,"a"])
          }
          
          df_power <- cbind(df_power, power_hte_col)
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_power %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd(),
                            r1_r0 == cicc_ratio_est_swd(),
                            a0 == oicc_wperiod_est_swd(),
                            a1_a0 == oicc_ratio_est_swd()) %>%
              mutate(r1_r0=factor(r1_r0)) %>%
              plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                      linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                   "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                   "<br>Clusters (total):", n, "; Clusters (per sequence):", n/(J-1),"<br>Cluster size (per period): ", m,
                                   "<br>Sequences (steps): ", (J-1),
                                   "<br>HTE power: ", round(power_hte_col,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-period covariate ICC = ", cicc_wperiod_est_swd()),
                     #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Cluster size (per period)"),
                     yaxis=list(title="HTE Power"),
                     legend=list(title=list(text="covariate ICC ratio")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size (per period) vs HTE power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), r0 (", cicc_wperiod_est_swd(),") and r1/r0 (", cicc_ratio_est_swd(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_swd == "oICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
              
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", n/(J-1),"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", n/(J-1),"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), minimum r0 (", cicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), and maximum r0 (", cicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), and r0 (", cicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", n/(J-1),"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-subcluster covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", n/(J-1),"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-subcluster covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), minimum a0 (", oicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), and maximum a0 (", oicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), "), r1/r0 (", cicc_ratio_est_swd(),"), and a0 (", oicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }else if(cohort() == "closed"){
          if(sensitivity_swd_react() == "est_only"){
            df_power <- expand.grid(n=n, m=m_range, J=J,
                                    a0=#c(oicc_wperiod_min_swd_cc(),
                                      oicc_wperiod_est_swd_cc(),
                                    # oicc_wperiod_max_swd_cc()),
                                    a1_a0=#c(oicc_ratio_min_swd_cc(),
                                      oicc_ratio_est_swd_cc(),
                                    # oicc_ratio_max_swd_cc()),
                                    a2=#c(oicc_windiv_min_swd_cc(),
                                      oicc_windiv_est_swd_cc(),
                                    # oicc_windiv_max_swd_cc()),
                                    r0=#c(cicc_wperiod_min_swd_cc(),
                                      cicc_wperiod_est_swd_cc(),
                                    # cicc_wperiod_max_swd_cc()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
          }else{
            df_power <- expand.grid(n=n, m=m_range, J=J,
                                    a0=c(oicc_wperiod_min_swd_cc(),
                                         oicc_wperiod_est_swd_cc(),
                                         oicc_wperiod_max_swd_cc()),
                                    a1_a0=c(oicc_ratio_min_swd_cc(),
                                            oicc_ratio_est_swd_cc(),
                                            oicc_ratio_max_swd_cc()),
                                    a2=c(oicc_windiv_min_swd_cc(),
                                         oicc_windiv_est_swd_cc(),
                                         oicc_windiv_max_swd_cc()),
                                    r0=c(cicc_wperiod_min_swd_cc(),
                                         cicc_wperiod_est_swd_cc(),
                                         cicc_wperiod_max_swd_cc()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
          }
          
          
          power_hte_col <- rep(NA, nrow(df_power))
          for(i in seq(nrow(df_power))){
            power_hte_col[i] <- power_hte_swd(n=df_power[i,"n"],
                                              m=df_power[i,"m"], J=df_power[i,"J"],
                                              var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                              a0=df_power[i,"a0"], a1_a0=df_power[i,"a1_a0"],
                                              a2=df_power[i,"a2"],
                                              r0=df_power[i,"r0"],
                                              cohort=df_power[i,"cohort"],
                                              d=df_power[i,"d"], a=df_power[i,"a"])
            #if(is.nan(power_hte_col[i])) print(i)
          }
          
          df_power <- cbind(df_power, power_hte_col)
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_power %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                            a0 == oicc_wperiod_est_swd_cc(),
                            a1_a0 == oicc_ratio_est_swd_cc(),
                            a2 == oicc_windiv_est_swd_cc()) %>%
              mutate(r0=factor(r0)) %>%
              plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                      linetype=~r0, color=~r0,legendgroup=~r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0, "; <br>within-individual outcome ICC (a2): ", a2,
                                   "; <br>within-period covariate ICC (r0): ", r0,
                                   "<br>Clusters (total):", n, "; Clusters (per sequence):", n/(J-1),"<br>Cluster size (per period): ", m,
                                   "<br>Sequences (steps): ", (J-1),
                                   "<br>HTE power: ", round(power_hte_col,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-individual outcome ICC = ", oicc_windiv_est_swd_cc()),
                     #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Cluster size (per period)"),
                     yaxis=list(title="HTE Power"),
                     legend=list(title=list(text="covariate ICC")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size (per period) vs HTE power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), a2 (", oicc_windiv_est_swd_cc(), "), and r0 (", cicc_wperiod_est_swd(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_swd == "oICC_constant"){
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              a2_unique <- sort(unique(df_power[,"a2"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", n/(J-1),"<br>Cluster size (total): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate ICC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", n/(J-1),"<br>Cluster size (total): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate ICC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), minimum r0 (", cicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), maximum r0 (", cicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), a2 (", oicc_windiv_est_swd_cc(),"), and r0 (", cicc_wperiod_est_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              a2_unique <- sort(unique(df_power[,"a2"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", n/(J-1),"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-subcluster covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", n/(J-1),"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-subcluster covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), minimum a0 (", oicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), and maximum a0 (", oicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), "), a2 (", oicc_windiv_est_swd_cc(),"), and a0 (", oicc_wperiod_est_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }#end cross/closed if/else
        
      }else if(plot_display_react() == "n_v_power"){
        
        ns_range <- seq(input$ns_swd_range[1],input$ns_swd_range[2])
        J <- input$J_1 + 1
        n_range <- ns_range*input$J_1
        
        if(cohort() == "cross"){
          
          if(sensitivity_swd_react() == "est_only"){
            df_power <- expand.grid(n=n_range, m=input$m_swd,J=J,
                                    a0=#c(oicc_wperiod_min_swd(),
                                      oicc_wperiod_est_swd(),
                                    # oicc_wperiod_max_swd()),
                                    a1_a0=#c(oicc_ratio_min_swd(),
                                      oicc_ratio_est_swd(),
                                    #  oicc_ratio_max_swd()),
                                    r0=#c(cicc_wperiod_min_swd(),
                                      cicc_wperiod_est_swd(),
                                    # cicc_wperiod_max_swd()),
                                    r1_r0=#c(cicc_ratio_min_swd(),
                                      cicc_ratio_est_swd(),
                                    # cicc_ratio_max_swd()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
            
          }else{
            df_power <- expand.grid(n=n_range, m=input$m_swd,J=J,
                                    a0=c(oicc_wperiod_min_swd(),
                                         oicc_wperiod_est_swd(),
                                         oicc_wperiod_max_swd()),
                                    a1_a0=c(oicc_ratio_min_swd(),
                                            oicc_ratio_est_swd(),
                                            oicc_ratio_max_swd()),
                                    r0=c(cicc_wperiod_min_swd(),
                                         cicc_wperiod_est_swd(),
                                         cicc_wperiod_max_swd()),
                                    r1_r0=c(cicc_ratio_min_swd(),
                                            cicc_ratio_est_swd(),
                                            cicc_ratio_max_swd()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
            
          }
          
          power_hte_col <- rep(NA, nrow(df_power))
          for(i in seq(nrow(df_power))){
            power_hte_col[i] <- power_hte_swd(n=df_power[i,"n"],
                                              m=df_power[i,"m"], J=df_power[i,"J"],
                                              var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                              a0=df_power[i,"a0"], a1_a0=df_power[i,"a1_a0"],
                                              r0=df_power[i,"r0"], r1_r0=df_power[i,"r1_r0"],
                                              cohort=df_power[i,"cohort"],
                                              d=df_power[i,"d"], a=df_power[i,"a"])
          }
          
          df_power <- cbind(df_power, power_hte_col)
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_power %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd(),
                            r1_r0 == cicc_ratio_est_swd(),
                            a0 == oicc_wperiod_est_swd(),
                            a1_a0 == oicc_ratio_est_swd()) %>%
              mutate(r1_r0=factor(r1_r0),
                     ns=n/(J-1)) %>%
              plot_ly(x=~ns,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                      linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                   "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                   "<br>Clusters (total):", n, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                   "<br>Sequences (steps): ", (J-1),
                                   "<br>HTE power: ", round(power_hte_col,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-period covariate ICC = ", cicc_wperiod_est_swd()),
                     #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Number of clusters (per sequence)"),
                     yaxis=list(title="HTE Power"),
                     legend=list(title=list(text="covariate ICC ratio")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Number of clusters (per sequence) vs HTE Power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), r0 (", cicc_wperiod_est_swd(),") and r1/r0 (", cicc_ratio_est_swd(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_three == "oICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0),
                           ns=n/(J-1)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~ns,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0),
                           ns=n/(J-1)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~ns,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters (per sequence) vs HTE Power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), minimum r0 (", cicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), and maximum r0 (", cicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), and r0 (", cicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0),
                           ns=n/(J-1)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~ns,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0),
                           ns=n/(J-1)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~ns,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters (per sequence) vs HTE Power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), minimum a0 (", oicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), and maximum a0 (", oicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), "), r1/r0 (", cicc_ratio_est_swd(),"), and a0 (", oicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }else if (cohort()== "closed"){
          
          if(sensitivity_swd_react() == "est_only"){
            df_power <- expand.grid(n=n_range, m=input$m_swd,J=J,
                                    a0=#c(oicc_wperiod_min_swd(),
                                      oicc_wperiod_est_swd_cc(),
                                    # oicc_wperiod_max_swd()),
                                    a1_a0=#c(oicc_ratio_min_swd(),
                                      oicc_ratio_est_swd_cc(),
                                    #  oicc_ratio_max_swd()),
                                    a2=#c(oicc_windiv_min_swd_cc(),
                                      oicc_windiv_est_swd_cc(),
                                    # oicc_windiv_max_swd_cc()),
                                    r0=#c(cicc_wperiod_min_swd_cc(),
                                      cicc_wperiod_est_swd_cc(),
                                    # cicc_wperiod_max_swd_cc()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
            
          }else{
            df_power <- expand.grid(n=n_range, m=input$m_swd,J=J,
                                    a0=c(oicc_wperiod_min_swd_cc(),
                                         oicc_wperiod_est_swd_cc(),
                                         oicc_wperiod_max_swd_cc()),
                                    a1_a0=c(oicc_ratio_min_swd_cc(),
                                            oicc_ratio_est_swd_cc(),
                                            oicc_ratio_max_swd_cc()),
                                    a2=c(oicc_windiv_min_swd_cc(),
                                         oicc_windiv_est_swd_cc(),
                                         oicc_windiv_max_swd_cc()),
                                    r0=c(cicc_wperiod_min_swd_cc(),
                                         cicc_wperiod_est_swd_cc(),
                                         cicc_wperiod_max_swd_cc()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
            
          }
          
          # df_power <- df_power <- cbind(df_power, a1=df_power$a0*df_power$a1_a0)
          # df_power <- df_power[which(df_power$a1 <=  df_power$a0),]
          
          power_hte_col <- rep(NA, nrow(df_power))
          for(i in seq(nrow(df_power))){
            power_hte_col[i] <- power_hte_swd(n=df_power[i,"n"],
                                              m=df_power[i,"m"], J=df_power[i,"J"],
                                              var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                              a0=df_power[i,"a0"], a1_a0=df_power[i,"a1_a0"],
                                              a2=df_power[i,"a2"], r0=df_power[i,"r0"],
                                              cohort=df_power[i,"cohort"],
                                              d=df_power[i,"d"], a=df_power[i,"a"])
          }
          
          df_power <- cbind(df_power, power_hte_col)
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_power %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                            a0 == oicc_wperiod_est_swd_cc(),
                            a1_a0 == oicc_ratio_est_swd_cc(),
                            a2 == oicc_windiv_est_swd_cc()) %>%
              mutate(r0=factor(r0),
                     ns=n/(J-1)) %>%
              plot_ly(x=~ns,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                      linetype=~r0, color=~r0,legendgroup=~r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0, "; <br>within-individual outcome ICC (a2): ", a2,
                                   "; <br>within-period covariate ICC (r0): ", r0,
                                   "<br>Clusters (total):", n, "; Clusters (per sequence):", n/(J-1),"<br>Cluster size (per period): ", m,
                                   "<br>Sequences (steps): ", (J-1),
                                   "<br>HTE power: ", round(power_hte_col,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-individual outcome ICC = ", oicc_windiv_est_swd_cc()),
                     #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Number of clusters (per sequence)"),
                     yaxis=list(title="HTE Power"),
                     legend=list(title=list(text="covariate ICC ratio")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Number of clusters (per sequence) vs HTE Power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), a2 (", oicc_windiv_est_swd_cc(), "), and r0 (", cicc_wperiod_est_swd_cc(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_three == "oICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              a2_unique <- sort(unique(df_power[,"a2"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0),
                           ns=n/(J-1)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~ns,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", n/(J-1),"<br>Cluster size (total): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0),
                           ns=n/(J-1)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~ns,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", n/(J-1),"<br>Cluster size (total): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters (per sequence) vs HTE Power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), minimum r0 (", cicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), maximum r0 (", cicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), a2 (", oicc_windiv_est_swd_cc(),"), and r0 (", cicc_wperiod_est_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              a2_unique <- sort(unique(df_power[,"a2"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0),
                           ns=n/(J-1)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~ns,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", n/(J-1),"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0),
                           ns=n/(J-1)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~ns,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n, "; Clusters (per sequence):", n/(J-1),"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", (J-1),
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters (per sequence) vs HTE Power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), minimum a0 (", oicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), and maximum a0 (", oicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), "), a2 (", oicc_windiv_est_swd_cc(),"), and a0 (", oicc_wperiod_est_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
        }#end cross/closed
        
      }else if(plot_display_react() == "fixed_power"){
        
        m_range <- seq(input$m_swd_range[1],input$m_swd_range[2])
        J <- input$J_1 + 1
        
        if(cohort() == "cross"){
          if(sensitivity_swd_react() == "est_only"){
            df_ns <- expand.grid(power=input$power_swd, m=m_range, J=J,
                                 a0=#c(oicc_wperiod_min_swd(),
                                   oicc_wperiod_est_swd(),
                                 # oicc_wperiod_max_swd()),
                                 a1_a0=#c(oicc_ratio_min_swd(),
                                   oicc_ratio_est_swd(),
                                 # oicc_ratio_max_swd()),
                                 r0=#c(cicc_wperiod_min_swd(),
                                   cicc_wperiod_est_swd(),
                                 # cicc_wperiod_max_swd()),
                                 r1_r0=#c(cicc_ratio_min_swd(),
                                   cicc_ratio_est_swd(),
                                 # cicc_ratio_max_swd()),
                                 var_y=var_y, #(input$sd_outcome)^2,
                                 var_x=var_x,#(input$sd_covar)^2,
                                 cohort=cohort(),
                                 d=input$mean_diff_HTE, a=input$sig)
            
          }else{
            df_ns <- expand.grid(power=input$power_swd, m=m_range, J=J,
                                 a0=c(oicc_wperiod_min_swd(),
                                      oicc_wperiod_est_swd(),
                                      oicc_wperiod_max_swd()),
                                 a1_a0=c(oicc_ratio_min_swd(),
                                         oicc_ratio_est_swd(),
                                         oicc_ratio_max_swd()),
                                 r0=c(cicc_wperiod_min_swd(),
                                      cicc_wperiod_est_swd(),
                                      cicc_wperiod_max_swd()),
                                 r1_r0=c(cicc_ratio_min_swd(),
                                         cicc_ratio_est_swd(),
                                         cicc_ratio_max_swd()),
                                 var_y=var_y, #(input$sd_outcome)^2,
                                 var_x=var_x,#(input$sd_covar)^2,
                                 cohort=cohort(),
                                 d=input$mean_diff_HTE, a=input$sig)
            
          }
          
          
          ns_hte_col <- matrix(NA, nrow(df_ns), ncol=2)
          for(i in seq(nrow(df_ns))){
            ns_hte_col[i,] <- unlist(ns_hte_swd(m=df_ns[i,"m"],
                                                J=df_ns[i,"J"],
                                                var_y=df_ns[i,"var_y"], var_x=df_ns[i,"var_x"],
                                                a0=df_ns[i,"a0"], a1_a0=df_ns[i,"a1_a0"],
                                                r0=df_ns[i,"r0"], r1_r0=df_ns[i,"r1_r0"],
                                                d=df_ns[i,"d"],
                                                cohort=df_ns[i,"cohort"],
                                                a=df_ns[i,"a"],
                                                power=df_ns[i,"power"]))
          }
          
          df_ns <- cbind(df_ns, ns_hte_col)
          colnames(df_ns)[(ncol(df_ns)-1):ncol(df_ns)] <- c("ns", "power_emp")
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_ns %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd(),
                            r1_r0 == cicc_ratio_est_swd(),
                            a0 == oicc_wperiod_est_swd(),
                            a1_a0 == oicc_ratio_est_swd()) %>%
              mutate(r1_r0=factor(r1_r0)) %>%
              plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                      linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                   "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                   "<br>Clusters (total):", ns*(J-1), "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                   "<br>Sequences (steps): ", input$J_1,
                                   "<br>HTE power: ", round(power_emp,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-period covariate ICC = ", cicc_wperiod_est_swd()),
                     #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Cluster size (per period)"),
                     yaxis=list(title="Number of clusters (per sequence)"),
                     legend=list(title=list(text="covariate ICC ratio")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size (per period) vs number of clusters (per sequence)',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), r0 (", cicc_wperiod_est_swd(),") and r1/r0 (", cicc_ratio_est_swd(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_swd == "oICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              
              a0_unique <- sort(unique(df_ns[,"a0"]))
              a1_a0_unique <- sort(unique(df_ns[,"a1_a0"]))
              r0_unique <- sort(unique(df_ns[,"r0"]))
              r1_r0_unique <- sort(unique(df_ns[,"r1_r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", ns*(J-1), "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", input$J_1,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", ns*(J-1), "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", input$J_1,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs number of clusters (per sequence)',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), minimum r0 (", cicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), and maximum r0 (", cicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), and r0 (", cicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_ns[,"a0"]))
              a1_a0_unique <- sort(unique(df_ns[,"a1_a0"]))
              r0_unique <- sort(unique(df_ns[,"r0"]))
              r1_r0_unique <- sort(unique(df_ns[,"r1_r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", ns*(J-1), "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", input$J_1,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           yaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", ns*(J-1), "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", input$J_1,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           yaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs number of clusters (per sequence)',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), minimum a0 (", oicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), and maximum a0 (", oicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), "), r1/r0 (", cicc_ratio_est_swd(),"), and a0 (", oicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }else if(cohort() == "closed"){
          if(sensitivity_swd_react() == "est_only"){
            df_ns <- expand.grid(power=input$power_swd, m=m_range, J=J,
                                 a0=#c(oicc_wperiod_min_swd(),
                                   oicc_wperiod_est_swd_cc(),
                                 # oicc_wperiod_max_swd()),
                                 a1_a0=#c(oicc_ratio_min_swd(),
                                   oicc_ratio_est_swd_cc(),
                                 # oicc_ratio_max_swd()),
                                 a2=#c(oicc_windiv_min_swd_cc(),
                                   oicc_windiv_est_swd_cc(),
                                 # oicc_windiv_max_swd_cc()),
                                 r0=#c(cicc_wperiod_min_swd_cc(),
                                   cicc_wperiod_est_swd_cc(),
                                 # cicc_wperiod_max_swd_cc()),
                                 var_y=var_y, #(input$sd_outcome)^2,
                                 var_x=var_x,#(input$sd_covar)^2,
                                 cohort=cohort(),
                                 d=input$mean_diff_HTE, a=input$sig)
            
          }else{
            df_ns <- expand.grid(power=input$power_swd, m=m_range, J=J,
                                 a0=c(oicc_wperiod_min_swd_cc(),
                                      oicc_wperiod_est_swd_cc(),
                                      oicc_wperiod_max_swd_cc()),
                                 a1_a0=c(oicc_ratio_min_swd_cc(),
                                         oicc_ratio_est_swd_cc(),
                                         oicc_ratio_max_swd_cc()),
                                 a2=c(oicc_windiv_min_swd_cc(),
                                      oicc_windiv_est_swd_cc(),
                                      oicc_windiv_max_swd_cc()),
                                 r0=c(cicc_wperiod_min_swd_cc(),
                                      cicc_wperiod_est_swd_cc(),
                                      cicc_wperiod_max_swd_cc()),
                                 var_y=var_y, #(input$sd_outcome)^2,
                                 var_x=var_x,#(input$sd_covar)^2,
                                 cohort=cohort(),
                                 d=input$mean_diff_HTE, a=input$sig)
            
          }
          
          
          ns_hte_col <- matrix(NA, nrow(df_ns), ncol=2)
          for(i in seq(nrow(df_ns))){
            ns_hte_col[i,] <- unlist(ns_hte_swd(m=df_ns[i,"m"],
                                                J=df_ns[i,"J"],
                                                var_y=df_ns[i,"var_y"], var_x=df_ns[i,"var_x"],
                                                a0=df_ns[i,"a0"], a1_a0=df_ns[i,"a1_a0"],
                                                a2=df_ns[i,"a2"], r0=df_ns[i,"r0"],
                                                d=df_ns[i,"d"],
                                                cohort=df_ns[i,"cohort"],
                                                a=df_ns[i,"a"],
                                                power=df_ns[i,"power"]))
          }
          
          df_ns <- cbind(df_ns, ns_hte_col)
          colnames(df_ns)[(ncol(df_ns)-1):ncol(df_ns)] <- c("ns", "power_emp")
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_ns %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                            a0 == oicc_wperiod_est_swd_cc(),
                            a1_a0 == oicc_ratio_est_swd_cc(),
                            a2 == oicc_windiv_est_swd_cc()) %>%
              mutate(r0=factor(r0)) %>%
              plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r0,
                      linetype=~r0, color=~r0,legendgroup=~r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0, "; <br>within-individual outcome ICC (a2): ", a2,
                                   "; <br>within-period covariate ICC (r0): ", r0,
                                   "<br>Clusters (total):", ns*(J-1), "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                   "<br>Sequences (steps): ", input$J_1,
                                   "<br>HTE power: ", round(power_emp,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-individual outcome ICC = ", oicc_windiv_est_swd_cc()),
                     #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Cluster size (per period)"),
                     yaxis=list(title="Number of clusters (per sequence)"),
                     legend=list(title=list(text="covariate ICC ratio")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size (per period) vs number of clusters (per sequence)',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), a2 (", oicc_windiv_est_swd_cc(), "), and r0 (", cicc_wperiod_est_swd_cc(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_swd == "oICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              
              a0_unique <- sort(unique(df_ns[,"a0"]))
              a1_a0_unique <- sort(unique(df_ns[,"a1_a0"]))
              a2_unique <- sort(unique(df_ns[,"a2"]))
              r0_unique <- sort(unique(df_ns[,"r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", ns*(J-1), "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", input$J_1,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", ns*(J-1), "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", input$J_1,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs number of clusters (per sequence)',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), minimum r0 (", cicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), maximum r0 (", cicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), a2 (", oicc_windiv_est_swd_cc(),"), and r0 (", cicc_wperiod_est_swd_cc(),")</i>"),                      xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_ns[,"a0"]))
              a1_a0_unique <- sort(unique(df_ns[,"a1_a0"]))
              a2_unique <- sort(unique(df_ns[,"a2"]))
              r0_unique <- sort(unique(df_ns[,"r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", ns*(J-1), "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", input$J_1,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           yaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", ns*(J-1), "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", input$J_1,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           yaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs number of clusters (per sequence)',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), minimum a0 (", oicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), and maximum a0 (", oicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), "), a2 (", oicc_windiv_est_swd_cc(),"), and a0 (", oicc_wperiod_est_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }#end cross/closed
        
        
      }#end fixed power
      
      #### UPLOAD OWN DESIGN ####
    }else if(trial_react() == 'parallel_m' | trial_react() == "upload"){
      
      if(is.null(file1()) & trial_react()=="upload") stop("User needs to upload design matrix before the function can continue")
      
      if(trial_react() == "upload"){
        desmat <- read.csv(file1()$datapath, header=FALSE) 
      }
      
      if(trial_react() == "parallel_m"){
        desmat <- designMatrix(design = trial_react(), periods = input$J)
      }
      
      # determine effect sizes and variances depending on outcome/covariate type #
      if(outcome_type_react() =="continuous"){
        var_y <- (sd_outcome_react())^2
        d <- input$mean_diff_HTE
        
        if(covar_type_react()=="continuous"){
          # continuous outcome and covar power
          var_x <- (sd_covar_react())^2
          
        }else if(covar_type_react() == "binary"){
          # continuous outcome binary covar power
          var_x <- (prop_covar_react())*(1-prop_covar_react())
        }
        
      }else if(outcome_type_react() == "binary"){
        var_y <- ((prop_control_react()*(1-prop_control_react())) + (prop_trt_react()*(1-prop_trt_react())))/2
        #d <- (input$prop_trt - input$prop_control)
        
        if(covar_type_react() == "continuous"){
          # binary outcome continuous covar power
          var_x <-(sd_covar_react())^2
          
        }else if(covar_type_react() == "binary"){
          # binary outcome and covar power #
          var_x <- (prop_covar_react())*(1-prop_covar_react())
          
        }
      }
      
      J <- ncol(desmat)
      seqs <- nrow(desmat)
      
      if(plot_display_react() == "m_v_power"){
        m_range <- seq(input$m_swd_range[1],input$m_swd_range[2])
        n_seq <- input$ns_swd
        #J <- input$J_1 + 1
        #n <- input$ns_swd*input$J_1
        
        if(cohort()=="cross"){
          
          if(sensitivity_swd_react() == "est_only"){
            df_power <- expand.grid(n_seq=n_seq, m=m_range, J=J,seqs=seqs,
                                    a0=#c(oicc_wperiod_min_swd(),
                                      oicc_wperiod_est_swd(),
                                    #oicc_wperiod_max_swd()),
                                    a1_a0=#c(oicc_ratio_min_swd(),
                                      oicc_ratio_est_swd(),
                                    # oicc_ratio_max_swd()),
                                    r0=#c(cicc_wperiod_min_swd(),
                                      cicc_wperiod_est_swd(),
                                    # cicc_wperiod_max_swd()),
                                    r1_r0=#c(cicc_ratio_min_swd(),
                                      cicc_ratio_est_swd(),
                                    # cicc_ratio_max_swd()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
          }else{
            df_power <- expand.grid(n_seq=n_seq, m=m_range, J=J,seqs=seqs,
                                    a0=c(oicc_wperiod_min_swd(),
                                         oicc_wperiod_est_swd(),
                                         oicc_wperiod_max_swd()),
                                    a1_a0=c(oicc_ratio_min_swd(),
                                            oicc_ratio_est_swd(),
                                            oicc_ratio_max_swd()),
                                    r0=c(cicc_wperiod_min_swd(),
                                         cicc_wperiod_est_swd(),
                                         cicc_wperiod_max_swd()),
                                    r1_r0=c(cicc_ratio_min_swd(),
                                            cicc_ratio_est_swd(),
                                            cicc_ratio_max_swd()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
          }
          
          
          power_hte_col <- rep(NA, nrow(df_power))
          for(i in seq(nrow(df_power))){
            power_hte_col[i] <- power_hte_ownDes(desmat, n_seq=df_power[i,"n_seq"],
                                                 m=df_power[i,"m"],# J=df_power[i,"J"],
                                                 var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                                 a0=df_power[i,"a0"], a1_a0=df_power[i,"a1_a0"],
                                                 r0=df_power[i,"r0"], r1_r0=df_power[i,"r1_r0"],
                                                 cohort=df_power[i,"cohort"],
                                                 d=df_power[i,"d"], a=df_power[i,"a"])
          }
          
          df_power <- cbind(df_power, power_hte_col)
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_power %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd(),
                            r1_r0 == cicc_ratio_est_swd(),
                            a0 == oicc_wperiod_est_swd(),
                            a1_a0 == oicc_ratio_est_swd()) %>%
              mutate(r1_r0=factor(r1_r0)) %>%
              plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                      linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                   "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                   "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                   "<br>Sequences (steps): ", seqs,
                                   "<br>HTE power: ", round(power_hte_col,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-period covariate ICC = ", cicc_wperiod_est_swd()),
                     #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Cluster size (per period)"),
                     yaxis=list(title="HTE Power"),
                     legend=list(title=list(text="covariate ICC ratio")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size (per period) vs HTE power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), r0 (", cicc_wperiod_est_swd(),") and r1/r0 (", cicc_ratio_est_swd(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_swd == "oICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
              
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), minimum r0 (", cicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), and maximum r0 (", cicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), and r0 (", cicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-subcluster covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-subcluster covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), minimum a0 (", oicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), and maximum a0 (", oicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), "), r1/r0 (", cicc_ratio_est_swd(),"), and a0 (", oicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }else if(cohort() == "closed"){
          if(sensitivity_swd_react() == "est_only"){
            df_power <- expand.grid(n_seq=n_seq, m=m_range, J=J,seqs=seqs,
                                    a0=#c(oicc_wperiod_min_swd_cc(),
                                      oicc_wperiod_est_swd_cc(),
                                    # oicc_wperiod_max_swd_cc()),
                                    a1_a0=#c(oicc_ratio_min_swd_cc(),
                                      oicc_ratio_est_swd_cc(),
                                    # oicc_ratio_max_swd_cc()),
                                    a2=#c(oicc_windiv_min_swd_cc(),
                                      oicc_windiv_est_swd_cc(),
                                    # oicc_windiv_max_swd_cc()),
                                    r0=#c(cicc_wperiod_min_swd_cc(),
                                      cicc_wperiod_est_swd_cc(),
                                    # cicc_wperiod_max_swd_cc()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
          }else{
            df_power <- expand.grid(n_seq=n_seq, m=m_range, J=J,seqs=seqs,
                                    a0=c(oicc_wperiod_min_swd_cc(),
                                         oicc_wperiod_est_swd_cc(),
                                         oicc_wperiod_max_swd_cc()),
                                    a1_a0=c(oicc_ratio_min_swd_cc(),
                                            oicc_ratio_est_swd_cc(),
                                            oicc_ratio_max_swd_cc()),
                                    a2=c(oicc_windiv_min_swd_cc(),
                                         oicc_windiv_est_swd_cc(),
                                         oicc_windiv_max_swd_cc()),
                                    r0=c(cicc_wperiod_min_swd_cc(),
                                         cicc_wperiod_est_swd_cc(),
                                         cicc_wperiod_max_swd_cc()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
          }
          
          
          power_hte_col <- rep(NA, nrow(df_power))
          for(i in seq(nrow(df_power))){
            power_hte_col[i] <- power_hte_ownDes(desmat,n_seq=df_power[i,"n_seq"],
                                                 m=df_power[i,"m"],# J=df_power[i,"J"],
                                                 var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                                 a0=df_power[i,"a0"], a1_a0=df_power[i,"a1_a0"],
                                                 a2=df_power[i,"a2"],
                                                 r0=df_power[i,"r0"],
                                                 cohort=df_power[i,"cohort"],
                                                 d=df_power[i,"d"], a=df_power[i,"a"])
            # if(is.nan(power_hte_col[i])) print(i)
          }
          
          df_power <- cbind(df_power, power_hte_col)
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_power %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                            a0 == oicc_wperiod_est_swd_cc(),
                            a1_a0 == oicc_ratio_est_swd_cc(),
                            a2 == oicc_windiv_est_swd_cc()) %>%
              mutate(r0=factor(r0)) %>%
              plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                      linetype=~r0, color=~r0,legendgroup=~r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0, "; <br>within-individual outcome ICC (a2): ", a2,
                                   "; <br>within-period covariate ICC (r0): ", r0,
                                   "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                   "<br>Sequences (steps): ", seqs,
                                   "<br>HTE power: ", round(power_hte_col,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-individual outcome ICC = ", oicc_windiv_est_swd_cc()),
                     #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Cluster size (per period)"),
                     yaxis=list(title="HTE Power"),
                     legend=list(title=list(text="covariate ICC")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size (per period) vs HTE power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), a2 (", oicc_windiv_est_swd_cc(), "), and r0 (", cicc_wperiod_est_swd(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_swd == "oICC_constant"){
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              a2_unique <- sort(unique(df_power[,"a2"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (total): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate ICC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (total): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate ICC")),
                           margin=0.01)
                }
                
              } 
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), minimum r0 (", cicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), maximum r0 (", cicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), a2 (", oicc_windiv_est_swd_cc(),"), and r0 (", cicc_wperiod_est_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              a2_unique <- sort(unique(df_power[,"a2"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-subcluster covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-subcluster covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), minimum a0 (", oicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), and maximum a0 (", oicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), "), a2 (", oicc_windiv_est_swd_cc(),"), and a0 (", oicc_wperiod_est_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }#end cross/closed if/else
        
      }else if(plot_display_react() == "n_v_power"){
        
        ns_range <- seq(input$ns_swd_range[1],input$ns_swd_range[2])
        #J <- input$J_1 + 1
        #n_range <- ns_range*input$J_1
        
        if(cohort() == "cross"){
          
          if(sensitivity_swd_react() == "est_only"){
            df_power <- expand.grid(n_seq=ns_range, m=input$m_swd,J=J,seqs=seqs,
                                    a0=#c(oicc_wperiod_min_swd(),
                                      oicc_wperiod_est_swd(),
                                    # oicc_wperiod_max_swd()),
                                    a1_a0=#c(oicc_ratio_min_swd(),
                                      oicc_ratio_est_swd(),
                                    #  oicc_ratio_max_swd()),
                                    r0=#c(cicc_wperiod_min_swd(),
                                      cicc_wperiod_est_swd(),
                                    # cicc_wperiod_max_swd()),
                                    r1_r0=#c(cicc_ratio_min_swd(),
                                      cicc_ratio_est_swd(),
                                    # cicc_ratio_max_swd()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
            
          }else{
            df_power <- expand.grid(n_seq=ns_range, m=input$m_swd,J=J,seqs=seqs,
                                    a0=c(oicc_wperiod_min_swd(),
                                         oicc_wperiod_est_swd(),
                                         oicc_wperiod_max_swd()),
                                    a1_a0=c(oicc_ratio_min_swd(),
                                            oicc_ratio_est_swd(),
                                            oicc_ratio_max_swd()),
                                    r0=c(cicc_wperiod_min_swd(),
                                         cicc_wperiod_est_swd(),
                                         cicc_wperiod_max_swd()),
                                    r1_r0=c(cicc_ratio_min_swd(),
                                            cicc_ratio_est_swd(),
                                            cicc_ratio_max_swd()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
            
          }
          
          power_hte_col <- rep(NA, nrow(df_power))
          for(i in seq(nrow(df_power))){
            power_hte_col[i] <- power_hte_ownDes(desmat,n_seq=df_power[i,"n_seq"],
                                                 m=df_power[i,"m"],# J=df_power[i,"J"],
                                                 var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                                 a0=df_power[i,"a0"], a1_a0=df_power[i,"a1_a0"],
                                                 r0=df_power[i,"r0"], r1_r0=df_power[i,"r1_r0"],
                                                 cohort=df_power[i,"cohort"],
                                                 d=df_power[i,"d"], a=df_power[i,"a"])
          }
          
          df_power <- cbind(df_power, power_hte_col)
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_power %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd(),
                            r1_r0 == cicc_ratio_est_swd(),
                            a0 == oicc_wperiod_est_swd(),
                            a1_a0 == oicc_ratio_est_swd()) %>%
              mutate(r1_r0=factor(r1_r0)) %>%
              plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                      linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                   "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                   "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                   "<br>Sequences (steps): ", seqs,
                                   "<br>HTE power: ", round(power_hte_col,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-period covariate ICC = ", cicc_wperiod_est_swd()),
                     #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Number of clusters (per sequence)"),
                     yaxis=list(title="HTE Power"),
                     legend=list(title=list(text="covariate ICC ratio")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Number of clusters (per sequence) vs HTE Power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), r0 (", cicc_wperiod_est_swd(),") and r1/r0 (", cicc_ratio_est_swd(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_three == "oICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters (per sequence) vs HTE Power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), minimum r0 (", cicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), and maximum r0 (", cicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), and r0 (", cicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters (per sequence) vs HTE Power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), minimum a0 (", oicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), and maximum a0 (", oicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), "), r1/r0 (", cicc_ratio_est_swd(),"), and a0 (", oicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }else if (cohort()== "closed"){
          
          if(sensitivity_swd_react() == "est_only"){
            df_power <- expand.grid(n_seq=ns_range, m=input$m_swd,J=J,seqs=seqs,
                                    a0=#c(oicc_wperiod_min_swd(),
                                      oicc_wperiod_est_swd_cc(),
                                    # oicc_wperiod_max_swd()),
                                    a1_a0=#c(oicc_ratio_min_swd(),
                                      oicc_ratio_est_swd_cc(),
                                    #  oicc_ratio_max_swd()),
                                    a2=#c(oicc_windiv_min_swd_cc(),
                                      oicc_windiv_est_swd_cc(),
                                    # oicc_windiv_max_swd_cc()),
                                    r0=#c(cicc_wperiod_min_swd_cc(),
                                      cicc_wperiod_est_swd_cc(),
                                    # cicc_wperiod_max_swd_cc()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
            
          }else{
            df_power <- expand.grid(n_seq=ns_range, m=input$m_swd,J=J,seqs=seqs,
                                    a0=c(oicc_wperiod_min_swd_cc(),
                                         oicc_wperiod_est_swd_cc(),
                                         oicc_wperiod_max_swd_cc()),
                                    a1_a0=c(oicc_ratio_min_swd_cc(),
                                            oicc_ratio_est_swd_cc(),
                                            oicc_ratio_max_swd_cc()),
                                    a2=c(oicc_windiv_min_swd_cc(),
                                         oicc_windiv_est_swd_cc(),
                                         oicc_windiv_max_swd_cc()),
                                    r0=c(cicc_wperiod_min_swd_cc(),
                                         cicc_wperiod_est_swd_cc(),
                                         cicc_wperiod_max_swd_cc()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
            
          }
          
          # df_power <- df_power <- cbind(df_power, a1=df_power$a0*df_power$a1_a0)
          # df_power <- df_power[which(df_power$a1 <=  df_power$a0),]
          
          power_hte_col <- rep(NA, nrow(df_power))
          for(i in seq(nrow(df_power))){
            power_hte_col[i] <- power_hte_ownDes(desmat,n_seq=df_power[i,"n_seq"],
                                                 m=df_power[i,"m"], #J=df_power[i,"J"],
                                                 var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                                 a0=df_power[i,"a0"], a1_a0=df_power[i,"a1_a0"],
                                                 a2=df_power[i,"a2"], r0=df_power[i,"r0"],
                                                 cohort=df_power[i,"cohort"],
                                                 d=df_power[i,"d"], a=df_power[i,"a"])
          }
          
          df_power <- cbind(df_power, power_hte_col)
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_power %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                            a0 == oicc_wperiod_est_swd_cc(),
                            a1_a0 == oicc_ratio_est_swd_cc(),
                            a2 == oicc_windiv_est_swd_cc()) %>%
              mutate(r0=factor(r0)) %>%
              plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                      linetype=~r0, color=~r0,legendgroup=~r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0, "; <br>within-individual outcome ICC (a2): ", a2,
                                   "; <br>within-period covariate ICC (r0): ", r0,
                                   "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                   "<br>Sequences (steps): ", seqs,
                                   "<br>HTE power: ", round(power_hte_col,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-individual outcome ICC = ", oicc_windiv_est_swd_cc()),
                     #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Number of clusters (per sequence)"),
                     yaxis=list(title="HTE Power"),
                     legend=list(title=list(text="covariate ICC ratio")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Number of clusters (per sequence) vs HTE Power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), a2 (", oicc_windiv_est_swd_cc(), "), and r0 (", cicc_wperiod_est_swd_cc(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_three == "oICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              a2_unique <- sort(unique(df_power[,"a2"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (total): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (total): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters (per sequence) vs HTE Power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), minimum r0 (", cicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), maximum r0 (", cicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), a2 (", oicc_windiv_est_swd_cc(),"), and r0 (", cicc_wperiod_est_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              a2_unique <- sort(unique(df_power[,"a2"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*seqs, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters (per sequence) vs HTE Power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), minimum a0 (", oicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), and maximum a0 (", oicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), "), a2 (", oicc_windiv_est_swd_cc(),"), and a0 (", oicc_wperiod_est_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }#end cross/closed
        
      }else if(plot_display_react() == "fixed_power"){
        
        m_range <- seq(input$m_swd_range[1],input$m_swd_range[2])
        #J <- input$J_1 + 1
        
        if(cohort() == "cross"){
          if(sensitivity_swd_react() == "est_only"){
            df_ns <- expand.grid(power=input$power_swd, m=m_range, J=J,seqs=seqs,
                                 a0=#c(oicc_wperiod_min_swd(),
                                   oicc_wperiod_est_swd(),
                                 # oicc_wperiod_max_swd()),
                                 a1_a0=#c(oicc_ratio_min_swd(),
                                   oicc_ratio_est_swd(),
                                 # oicc_ratio_max_swd()),
                                 r0=#c(cicc_wperiod_min_swd(),
                                   cicc_wperiod_est_swd(),
                                 # cicc_wperiod_max_swd()),
                                 r1_r0=#c(cicc_ratio_min_swd(),
                                   cicc_ratio_est_swd(),
                                 # cicc_ratio_max_swd()),
                                 var_y=var_y, #(input$sd_outcome)^2,
                                 var_x=var_x,#(input$sd_covar)^2,
                                 cohort=cohort(),
                                 d=input$mean_diff_HTE, a=input$sig)
            
          }else{
            df_ns <- expand.grid(power=input$power_swd, m=m_range, J=J,seqs=seqs,
                                 a0=c(oicc_wperiod_min_swd(),
                                      oicc_wperiod_est_swd(),
                                      oicc_wperiod_max_swd()),
                                 a1_a0=c(oicc_ratio_min_swd(),
                                         oicc_ratio_est_swd(),
                                         oicc_ratio_max_swd()),
                                 r0=c(cicc_wperiod_min_swd(),
                                      cicc_wperiod_est_swd(),
                                      cicc_wperiod_max_swd()),
                                 r1_r0=c(cicc_ratio_min_swd(),
                                         cicc_ratio_est_swd(),
                                         cicc_ratio_max_swd()),
                                 var_y=var_y, #(input$sd_outcome)^2,
                                 var_x=var_x,#(input$sd_covar)^2,
                                 cohort=cohort(),
                                 d=input$mean_diff_HTE, a=input$sig)
            
          }
          
          
          ns_hte_col <- matrix(NA, nrow(df_ns), ncol=2)
          for(i in seq(nrow(df_ns))){
            ns_hte_col[i,] <- unlist(ns_hte_ownDes(desmat, m=df_ns[i,"m"],
                                                   #J=df_ns[i,"J"],
                                                   var_y=df_ns[i,"var_y"], var_x=df_ns[i,"var_x"],
                                                   a0=df_ns[i,"a0"], a1_a0=df_ns[i,"a1_a0"],
                                                   r0=df_ns[i,"r0"], r1_r0=df_ns[i,"r1_r0"],
                                                   d=df_ns[i,"d"],
                                                   cohort=df_ns[i,"cohort"],
                                                   a=df_ns[i,"a"],
                                                   power=df_ns[i,"power"]))
          }
          
          df_ns <- cbind(df_ns, ns_hte_col)
          colnames(df_ns)[(ncol(df_ns)-1):ncol(df_ns)] <- c("ns", "power_emp")
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_ns %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd(),
                            r1_r0 == cicc_ratio_est_swd(),
                            a0 == oicc_wperiod_est_swd(),
                            a1_a0 == oicc_ratio_est_swd()) %>%
              mutate(r1_r0=factor(r1_r0)) %>%
              plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                      linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                   "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                   "<br>Clusters (total):", ns*seqs, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                   "<br>Sequences (steps): ", seqs,
                                   "<br>HTE power: ", round(power_emp,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-period covariate ICC = ", cicc_wperiod_est_swd()),
                     #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Cluster size (per period)"),
                     yaxis=list(title="Number of clusters (per sequence)"),
                     legend=list(title=list(text="covariate ICC ratio")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size (per period) vs number of clusters (per sequence)',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), r0 (", cicc_wperiod_est_swd(),") and r1/r0 (", cicc_ratio_est_swd(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_swd == "oICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              
              a0_unique <- sort(unique(df_ns[,"a0"]))
              a1_a0_unique <- sort(unique(df_ns[,"a1_a0"]))
              r0_unique <- sort(unique(df_ns[,"r0"]))
              r1_r0_unique <- sort(unique(df_ns[,"r1_r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", ns*seqs, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", ns*seqs, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs number of clusters (per sequence)',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), minimum r0 (", cicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), and maximum r0 (", cicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), and r0 (", cicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_ns[,"a0"]))
              a1_a0_unique <- sort(unique(df_ns[,"a1_a0"]))
              r0_unique <- sort(unique(df_ns[,"r0"]))
              r1_r0_unique <- sort(unique(df_ns[,"r1_r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", ns*seqs, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           yaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", ns*seqs, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           yaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs number of clusters (per sequence)',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), minimum a0 (", oicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), and maximum a0 (", oicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), "), r1/r0 (", cicc_ratio_est_swd(),"), and a0 (", oicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }else if(cohort() == "closed"){
          if(sensitivity_swd_react() == "est_only"){
            df_ns <- expand.grid(power=input$power_swd, m=m_range, J=J,seqs=seqs,
                                 a0=#c(oicc_wperiod_min_swd(),
                                   oicc_wperiod_est_swd_cc(),
                                 # oicc_wperiod_max_swd()),
                                 a1_a0=#c(oicc_ratio_min_swd(),
                                   oicc_ratio_est_swd_cc(),
                                 # oicc_ratio_max_swd()),
                                 a2=#c(oicc_windiv_min_swd_cc(),
                                   oicc_windiv_est_swd_cc(),
                                 # oicc_windiv_max_swd_cc()),
                                 r0=#c(cicc_wperiod_min_swd_cc(),
                                   cicc_wperiod_est_swd_cc(),
                                 # cicc_wperiod_max_swd_cc()),
                                 var_y=var_y, #(input$sd_outcome)^2,
                                 var_x=var_x,#(input$sd_covar)^2,
                                 cohort=cohort(),
                                 d=input$mean_diff_HTE, a=input$sig)
            
          }else{
            df_ns <- expand.grid(power=input$power_swd, m=m_range, J=J,seqs=seqs,
                                 a0=c(oicc_wperiod_min_swd_cc(),
                                      oicc_wperiod_est_swd_cc(),
                                      oicc_wperiod_max_swd_cc()),
                                 a1_a0=c(oicc_ratio_min_swd_cc(),
                                         oicc_ratio_est_swd_cc(),
                                         oicc_ratio_max_swd_cc()),
                                 a2=c(oicc_windiv_min_swd_cc(),
                                      oicc_windiv_est_swd_cc(),
                                      oicc_windiv_max_swd_cc()),
                                 r0=c(cicc_wperiod_min_swd_cc(),
                                      cicc_wperiod_est_swd_cc(),
                                      cicc_wperiod_max_swd_cc()),
                                 var_y=var_y, #(input$sd_outcome)^2,
                                 var_x=var_x,#(input$sd_covar)^2,
                                 cohort=cohort(),
                                 d=input$mean_diff_HTE, a=input$sig)
            
          }
          
          
          ns_hte_col <- matrix(NA, nrow(df_ns), ncol=2)
          for(i in seq(nrow(df_ns))){
            ns_hte_col[i,] <- unlist(ns_hte_ownDes(desmat,m=df_ns[i,"m"],
                                                   #J=df_ns[i,"J"],
                                                   var_y=df_ns[i,"var_y"], var_x=df_ns[i,"var_x"],
                                                   a0=df_ns[i,"a0"], a1_a0=df_ns[i,"a1_a0"],
                                                   a2=df_ns[i,"a2"], r0=df_ns[i,"r0"],
                                                   d=df_ns[i,"d"],
                                                   cohort=df_ns[i,"cohort"],
                                                   a=df_ns[i,"a"],
                                                   power=df_ns[i,"power"]))
          }
          
          df_ns <- cbind(df_ns, ns_hte_col)
          colnames(df_ns)[(ncol(df_ns)-1):ncol(df_ns)] <- c("ns", "power_emp")
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_ns %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                            a0 == oicc_wperiod_est_swd_cc(),
                            a1_a0 == oicc_ratio_est_swd_cc(),
                            a2 == oicc_windiv_est_swd_cc()) %>%
              mutate(r0=factor(r0)) %>%
              plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r0,
                      linetype=~r0, color=~r0,legendgroup=~r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0, "; <br>within-individual outcome ICC (a2): ", a2,
                                   "; <br>within-period covariate ICC (r0): ", r0,
                                   "<br>Clusters (total):", ns*seqs, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                   "<br>Sequences (steps): ", seqs,
                                   "<br>HTE power: ", round(power_emp,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-individual outcome ICC = ", oicc_windiv_est_swd_cc()),
                     #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Cluster size (per period)"),
                     yaxis=list(title="Number of clusters (per sequence)"),
                     legend=list(title=list(text="covariate ICC ratio")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size (per period) vs number of clusters (per sequence)',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), a2 (", oicc_windiv_est_swd_cc(), "), and r0 (", cicc_wperiod_est_swd_cc(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_swd == "oICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              
              a0_unique <- sort(unique(df_ns[,"a0"]))
              a1_a0_unique <- sort(unique(df_ns[,"a1_a0"]))
              a2_unique <- sort(unique(df_ns[,"a2"]))
              r0_unique <- sort(unique(df_ns[,"r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", ns*seqs, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", ns*seqs, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs number of clusters (per sequence)',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), minimum r0 (", cicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), maximum r0 (", cicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), a2 (", oicc_windiv_est_swd_cc(),"), and r0 (", cicc_wperiod_est_swd_cc(),")</i>"),                      xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_ns[,"a0"]))
              a1_a0_unique <- sort(unique(df_ns[,"a1_a0"]))
              a2_unique <- sort(unique(df_ns[,"a2"]))
              r0_unique <- sort(unique(df_ns[,"r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", ns*seqs, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           yaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", ns*seqs, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Sequences (steps): ", seqs,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           yaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs number of clusters (per sequence)',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), minimum a0 (", oicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), and maximum a0 (", oicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), "), a2 (", oicc_windiv_est_swd_cc(),"), and a0 (", oicc_wperiod_est_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }#end cross/closed
        
        
      }#end fixed power
      
      #### TWO-PERIOD AND MULTI-PERIOD CROSSOVERS ####
    }else if(trial_react() == "crossover_2" | trial_react() == "crossover_m"){
      # determine effect sizes and variances depending on outcome/covariate type #
      if(outcome_type_react() =="continuous"){
        var_y <- (sd_outcome_react())^2
        d <- input$mean_diff_HTE
        
        if(covar_type_react()=="continuous"){
          # continuous outcome and covar power
          var_x <- (sd_covar_react())^2
          
        }else if(covar_type_react() == "binary"){
          # continuous outcome binary covar power
          var_x <- (prop_covar_react())*(1-prop_covar_react())
        }
        
      }else if(outcome_type_react() == "binary"){
        var_y <- ((prop_control_react()*(1-prop_control_react())) + (prop_trt_react()*(1-prop_trt_react())))/2
        #d <- (input$prop_trt - input$prop_control)
        
        if(covar_type_react() == "continuous"){
          # binary outcome continuous covar power
          var_x <-(sd_covar_react())^2
          
        }else if(covar_type_react() == "binary"){
          # binary outcome and covar power #
          var_x <- (prop_covar_react())*(1-prop_covar_react())
          
        }
      }
      
      if(trial_react() == "crossover_2"){
        J <- 2
      }else if(trial_react() == "crossover_m"){
        J <- input$J
      }
      
      if(plot_display_react() == "m_v_power"){
        m_range <- seq(input$m_swd_range[1],input$m_swd_range[2])
        n_seq <- input$ns_swd
        #J <- input$J_1 + 1
        #n <- input$ns_swd*input$J_1
        
        if(cohort()=="cross"){
          
          if(sensitivity_swd_react() == "est_only"){
            df_power <- expand.grid(n_seq=n_seq, m=m_range, J=J,
                                    a0=#c(oicc_wperiod_min_swd(),
                                      oicc_wperiod_est_swd(),
                                    #oicc_wperiod_max_swd()),
                                    a1_a0=#c(oicc_ratio_min_swd(),
                                      oicc_ratio_est_swd(),
                                    # oicc_ratio_max_swd()),
                                    r0=#c(cicc_wperiod_min_swd(),
                                      cicc_wperiod_est_swd(),
                                    # cicc_wperiod_max_swd()),
                                    r1_r0=#c(cicc_ratio_min_swd(),
                                      cicc_ratio_est_swd(),
                                    # cicc_ratio_max_swd()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
          }else{
            df_power <- expand.grid(n_seq=n_seq, m=m_range, J=J,
                                    a0=c(oicc_wperiod_min_swd(),
                                         oicc_wperiod_est_swd(),
                                         oicc_wperiod_max_swd()),
                                    a1_a0=c(oicc_ratio_min_swd(),
                                            oicc_ratio_est_swd(),
                                            oicc_ratio_max_swd()),
                                    r0=c(cicc_wperiod_min_swd(),
                                         cicc_wperiod_est_swd(),
                                         cicc_wperiod_max_swd()),
                                    r1_r0=c(cicc_ratio_min_swd(),
                                            cicc_ratio_est_swd(),
                                            cicc_ratio_max_swd()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
          }
          
          
          power_hte_col <- rep(NA, nrow(df_power))
          for(i in seq(nrow(df_power))){
            power_hte_col[i] <- power_hte_crossover(n_seq=df_power[i,"n_seq"],
                                                    m=df_power[i,"m"], J=df_power[i,"J"],
                                                    var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                                    a0=df_power[i,"a0"], a1_a0=df_power[i,"a1_a0"],
                                                    r0=df_power[i,"r0"], r1_r0=df_power[i,"r1_r0"],
                                                    cohort=df_power[i,"cohort"],
                                                    d=df_power[i,"d"], a=df_power[i,"a"])
          }
          
          df_power <- cbind(df_power, power_hte_col)
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_power %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd(),
                            r1_r0 == cicc_ratio_est_swd(),
                            a0 == oicc_wperiod_est_swd(),
                            a1_a0 == oicc_ratio_est_swd()) %>%
              mutate(r1_r0=factor(r1_r0)) %>%
              plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                      linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                   "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                   "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                   "<br>Periods (J): ", J,
                                   "<br>HTE power: ", round(power_hte_col,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-period covariate ICC = ", cicc_wperiod_est_swd()),
                     #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Cluster size (per period)"),
                     yaxis=list(title="HTE Power"),
                     legend=list(title=list(text="covariate ICC ratio")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size (per period) vs HTE power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), r0 (", cicc_wperiod_est_swd(),") and r1/r0 (", cicc_ratio_est_swd(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_swd == "oICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
              
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), minimum r0 (", cicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), and maximum r0 (", cicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), and r0 (", cicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-subcluster covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-subcluster covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), minimum a0 (", oicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), and maximum a0 (", oicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), "), r1/r0 (", cicc_ratio_est_swd(),"), and a0 (", oicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }else if(cohort() == "closed"){
          if(sensitivity_swd_react() == "est_only"){
            df_power <- expand.grid(n_seq=n_seq, m=m_range, J=J,
                                    a0=#c(oicc_wperiod_min_swd_cc(),
                                      oicc_wperiod_est_swd_cc(),
                                    # oicc_wperiod_max_swd_cc()),
                                    a1_a0=#c(oicc_ratio_min_swd_cc(),
                                      oicc_ratio_est_swd_cc(),
                                    # oicc_ratio_max_swd_cc()),
                                    a2=#c(oicc_windiv_min_swd_cc(),
                                      oicc_windiv_est_swd_cc(),
                                    # oicc_windiv_max_swd_cc()),
                                    r0=#c(cicc_wperiod_min_swd_cc(),
                                      cicc_wperiod_est_swd_cc(),
                                    # cicc_wperiod_max_swd_cc()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
          }else{
            df_power <- expand.grid(n_seq=n_seq, m=m_range, J=J,
                                    a0=c(oicc_wperiod_min_swd_cc(),
                                         oicc_wperiod_est_swd_cc(),
                                         oicc_wperiod_max_swd_cc()),
                                    a1_a0=c(oicc_ratio_min_swd_cc(),
                                            oicc_ratio_est_swd_cc(),
                                            oicc_ratio_max_swd_cc()),
                                    a2=c(oicc_windiv_min_swd_cc(),
                                         oicc_windiv_est_swd_cc(),
                                         oicc_windiv_max_swd_cc()),
                                    r0=c(cicc_wperiod_min_swd_cc(),
                                         cicc_wperiod_est_swd_cc(),
                                         cicc_wperiod_max_swd_cc()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
          }
          
          
          power_hte_col <- rep(NA, nrow(df_power))
          for(i in seq(nrow(df_power))){
            power_hte_col[i] <- power_hte_crossover(n_seq=df_power[i,"n_seq"],
                                                    m=df_power[i,"m"], J=df_power[i,"J"],
                                                    var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                                    a0=df_power[i,"a0"], a1_a0=df_power[i,"a1_a0"],
                                                    a2=df_power[i,"a2"],
                                                    r0=df_power[i,"r0"],
                                                    cohort=df_power[i,"cohort"],
                                                    d=df_power[i,"d"], a=df_power[i,"a"])
            # if(is.nan(power_hte_col[i])) print(i)
          }
          
          df_power <- cbind(df_power, power_hte_col)
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_power %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                            a0 == oicc_wperiod_est_swd_cc(),
                            a1_a0 == oicc_ratio_est_swd_cc(),
                            a2 == oicc_windiv_est_swd_cc()) %>%
              mutate(r0=factor(r0)) %>%
              plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                      linetype=~r0, color=~r0,legendgroup=~r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0, "; <br>within-individual outcome ICC (a2): ", a2,
                                   "; <br>within-period covariate ICC (r0): ", r0,
                                   "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                   "<br>Periods (J): ", J,
                                   "<br>HTE power: ", round(power_hte_col,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-individual outcome ICC = ", oicc_windiv_est_swd_cc()),
                     #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Cluster size (per period)"),
                     yaxis=list(title="HTE Power"),
                     legend=list(title=list(text="covariate ICC")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size (per period) vs HTE power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), a2 (", oicc_windiv_est_swd_cc(), "), and r0 (", cicc_wperiod_est_swd(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_swd == "oICC_constant"){
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              a2_unique <- sort(unique(df_power[,"a2"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (total): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate ICC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (total): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-subcluster outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate ICC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), minimum r0 (", cicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), maximum r0 (", cicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), a2 (", oicc_windiv_est_swd_cc(),"), and r0 (", cicc_wperiod_est_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              a2_unique <- sort(unique(df_power[,"a2"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br> Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-subcluster covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-subcluster covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), minimum a0 (", oicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), and maximum a0 (", oicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), "), a2 (", oicc_windiv_est_swd_cc(),"), and a0 (", oicc_wperiod_est_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
        }#end cross/closed if/else
        
      }else if(plot_display_react() == "n_v_power"){
        
        ns_range <- seq(input$ns_swd_range[1],input$ns_swd_range[2])
        #J <- input$J_1 + 1
        #n_range <- ns_range*input$J_1
        
        if(cohort() == "cross"){
          
          if(sensitivity_swd_react() == "est_only"){
            df_power <- expand.grid(n_seq=ns_range, m=input$m_swd,J=J,
                                    a0=#c(oicc_wperiod_min_swd(),
                                      oicc_wperiod_est_swd(),
                                    # oicc_wperiod_max_swd()),
                                    a1_a0=#c(oicc_ratio_min_swd(),
                                      oicc_ratio_est_swd(),
                                    #  oicc_ratio_max_swd()),
                                    r0=#c(cicc_wperiod_min_swd(),
                                      cicc_wperiod_est_swd(),
                                    # cicc_wperiod_max_swd()),
                                    r1_r0=#c(cicc_ratio_min_swd(),
                                      cicc_ratio_est_swd(),
                                    # cicc_ratio_max_swd()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
            
          }else{
            df_power <- expand.grid(n_seq=ns_range, m=input$m_swd,J=J,
                                    a0=c(oicc_wperiod_min_swd(),
                                         oicc_wperiod_est_swd(),
                                         oicc_wperiod_max_swd()),
                                    a1_a0=c(oicc_ratio_min_swd(),
                                            oicc_ratio_est_swd(),
                                            oicc_ratio_max_swd()),
                                    r0=c(cicc_wperiod_min_swd(),
                                         cicc_wperiod_est_swd(),
                                         cicc_wperiod_max_swd()),
                                    r1_r0=c(cicc_ratio_min_swd(),
                                            cicc_ratio_est_swd(),
                                            cicc_ratio_max_swd()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
            
          }
          
          power_hte_col <- rep(NA, nrow(df_power))
          for(i in seq(nrow(df_power))){
            power_hte_col[i] <- power_hte_crossover(n_seq=df_power[i,"n_seq"],
                                                    m=df_power[i,"m"], J=df_power[i,"J"],
                                                    var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                                    a0=df_power[i,"a0"], a1_a0=df_power[i,"a1_a0"],
                                                    r0=df_power[i,"r0"], r1_r0=df_power[i,"r1_r0"],
                                                    cohort=df_power[i,"cohort"],
                                                    d=df_power[i,"d"], a=df_power[i,"a"])
          }
          
          df_power <- cbind(df_power, power_hte_col)
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_power %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd(),
                            r1_r0 == cicc_ratio_est_swd(),
                            a0 == oicc_wperiod_est_swd(),
                            a1_a0 == oicc_ratio_est_swd()) %>%
              mutate(r1_r0=factor(r1_r0)) %>%
              plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                      linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                   "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                   "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                   "<br>Periods (J): ", J,
                                   "<br>HTE power: ", round(power_hte_col,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-period covariate ICC = ", cicc_wperiod_est_swd()),
                     #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Number of clusters (per sequence)"),
                     yaxis=list(title="HTE Power"),
                     legend=list(title=list(text="covariate ICC ratio")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Number of clusters (per sequence) vs HTE Power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), r0 (", cicc_wperiod_est_swd(),") and r1/r0 (", cicc_ratio_est_swd(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_three == "oICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters (per sequence) vs HTE Power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), minimum r0 (", cicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), and maximum r0 (", cicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), and r0 (", cicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              r1_r0_unique <- sort(unique(df_power[,"r1_r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters (per sequence) vs HTE Power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), minimum a0 (", oicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), and maximum a0 (", oicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), "), r1/r0 (", cicc_ratio_est_swd(),"), and a0 (", oicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }else if (cohort()== "closed"){
          
          if(sensitivity_swd_react() == "est_only"){
            df_power <- expand.grid(n_seq=ns_range, m=input$m_swd,J=J,
                                    a0=#c(oicc_wperiod_min_swd(),
                                      oicc_wperiod_est_swd_cc(),
                                    # oicc_wperiod_max_swd()),
                                    a1_a0=#c(oicc_ratio_min_swd(),
                                      oicc_ratio_est_swd_cc(),
                                    #  oicc_ratio_max_swd()),
                                    a2=#c(oicc_windiv_min_swd_cc(),
                                      oicc_windiv_est_swd_cc(),
                                    # oicc_windiv_max_swd_cc()),
                                    r0=#c(cicc_wperiod_min_swd_cc(),
                                      cicc_wperiod_est_swd_cc(),
                                    # cicc_wperiod_max_swd_cc()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
            
          }else{
            df_power <- expand.grid(n_seq=ns_range, m=input$m_swd,J=J,
                                    a0=c(oicc_wperiod_min_swd_cc(),
                                         oicc_wperiod_est_swd_cc(),
                                         oicc_wperiod_max_swd_cc()),
                                    a1_a0=c(oicc_ratio_min_swd_cc(),
                                            oicc_ratio_est_swd_cc(),
                                            oicc_ratio_max_swd_cc()),
                                    a2=c(oicc_windiv_min_swd_cc(),
                                         oicc_windiv_est_swd_cc(),
                                         oicc_windiv_max_swd_cc()),
                                    r0=c(cicc_wperiod_min_swd_cc(),
                                         cicc_wperiod_est_swd_cc(),
                                         cicc_wperiod_max_swd_cc()),
                                    var_y=var_y, #(input$sd_outcome)^2,
                                    var_x=var_x,#(input$sd_covar)^2,
                                    cohort=cohort(),
                                    d=input$mean_diff_HTE, a=input$sig)
            
          }
          
          # df_power <- df_power <- cbind(df_power, a1=df_power$a0*df_power$a1_a0)
          # df_power <- df_power[which(df_power$a1 <=  df_power$a0),]
          
          power_hte_col <- rep(NA, nrow(df_power))
          for(i in seq(nrow(df_power))){
            power_hte_col[i] <- power_hte_crossover(n_seq=df_power[i,"n_seq"],
                                                    m=df_power[i,"m"], J=df_power[i,"J"],
                                                    var_y=df_power[i,"var_y"], var_x=df_power[i,"var_x"],
                                                    a0=df_power[i,"a0"], a1_a0=df_power[i,"a1_a0"],
                                                    a2=df_power[i,"a2"], r0=df_power[i,"r0"],
                                                    cohort=df_power[i,"cohort"],
                                                    d=df_power[i,"d"], a=df_power[i,"a"])
          }
          
          df_power <- cbind(df_power, power_hte_col)
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_power %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                            a0 == oicc_wperiod_est_swd_cc(),
                            a1_a0 == oicc_ratio_est_swd_cc(),
                            a2 == oicc_windiv_est_swd_cc()) %>%
              mutate(r0=factor(r0)) %>%
              plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                      linetype=~r0, color=~r0,legendgroup=~r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0, "; <br>within-individual outcome ICC (a2): ", a2,
                                   "; <br>within-period covariate ICC (r0): ", r0,
                                   "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                   "<br>Periods (J): ", J,
                                   "<br>HTE power: ", round(power_hte_col,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-individual outcome ICC = ", oicc_windiv_est_swd_cc()),
                     #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Number of clusters (per sequence)"),
                     yaxis=list(title="HTE Power"),
                     legend=list(title=list(text="covariate ICC ratio")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Number of clusters (per sequence) vs HTE Power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), a2 (", oicc_windiv_est_swd_cc(), "), and r0 (", cicc_wperiod_est_swd_cc(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_three == "oICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              a2_unique <- sort(unique(df_power[,"a2"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (total): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (total): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters (per sequence) vs HTE Power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), minimum r0 (", cicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), maximum r0 (", cicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), a2 (", oicc_windiv_est_swd_cc(),"), and r0 (", cicc_wperiod_est_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_power[,"a0"]))
              a1_a0_unique <- sort(unique(df_power[,"a1_a0"]))
              a2_unique <- sort(unique(df_power[,"a2"]))
              r0_unique <- sort(unique(df_power[,"r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_power %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~n_seq,y=~power_hte_col, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", n_seq*2, "; Clusters (per sequence):", n_seq,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           xaxis=list(title="Number of clusters (per sequence)"),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters (per sequence) vs HTE Power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), minimum a0 (", oicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), and maximum a0 (", oicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), "), a2 (", oicc_windiv_est_swd_cc(),"), and a0 (", oicc_wperiod_est_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }#end cross/closed
        
      }else if(plot_display_react() == "fixed_power"){
        
        m_range <- seq(input$m_swd_range[1],input$m_swd_range[2])
        #J <- input$J_1 + 1
        
        if(cohort() == "cross"){
          if(sensitivity_swd_react() == "est_only"){
            df_ns <- expand.grid(power=input$power_swd, m=m_range, J=J,
                                 a0=#c(oicc_wperiod_min_swd(),
                                   oicc_wperiod_est_swd(),
                                 # oicc_wperiod_max_swd()),
                                 a1_a0=#c(oicc_ratio_min_swd(),
                                   oicc_ratio_est_swd(),
                                 # oicc_ratio_max_swd()),
                                 r0=#c(cicc_wperiod_min_swd(),
                                   cicc_wperiod_est_swd(),
                                 # cicc_wperiod_max_swd()),
                                 r1_r0=#c(cicc_ratio_min_swd(),
                                   cicc_ratio_est_swd(),
                                 # cicc_ratio_max_swd()),
                                 var_y=var_y, #(input$sd_outcome)^2,
                                 var_x=var_x,#(input$sd_covar)^2,
                                 cohort=cohort(),
                                 d=input$mean_diff_HTE, a=input$sig)
            
          }else{
            df_ns <- expand.grid(power=input$power_swd, m=m_range, J=J,
                                 a0=c(oicc_wperiod_min_swd(),
                                      oicc_wperiod_est_swd(),
                                      oicc_wperiod_max_swd()),
                                 a1_a0=c(oicc_ratio_min_swd(),
                                         oicc_ratio_est_swd(),
                                         oicc_ratio_max_swd()),
                                 r0=c(cicc_wperiod_min_swd(),
                                      cicc_wperiod_est_swd(),
                                      cicc_wperiod_max_swd()),
                                 r1_r0=c(cicc_ratio_min_swd(),
                                         cicc_ratio_est_swd(),
                                         cicc_ratio_max_swd()),
                                 var_y=var_y, #(input$sd_outcome)^2,
                                 var_x=var_x,#(input$sd_covar)^2,
                                 cohort=cohort(),
                                 d=input$mean_diff_HTE, a=input$sig)
            
          }
          
          
          ns_hte_col <- matrix(NA, nrow(df_ns), ncol=2)
          for(i in seq(nrow(df_ns))){
            ns_hte_col[i,] <- unlist(ns_hte_crossover(m=df_ns[i,"m"],
                                                      J=df_ns[i,"J"],
                                                      var_y=df_ns[i,"var_y"], var_x=df_ns[i,"var_x"],
                                                      a0=df_ns[i,"a0"], a1_a0=df_ns[i,"a1_a0"],
                                                      r0=df_ns[i,"r0"], r1_r0=df_ns[i,"r1_r0"],
                                                      d=df_ns[i,"d"],
                                                      cohort=df_ns[i,"cohort"],
                                                      a=df_ns[i,"a"],
                                                      power=df_ns[i,"power"]))
          }
          
          df_ns <- cbind(df_ns, ns_hte_col)
          colnames(df_ns)[(ncol(df_ns)-1):ncol(df_ns)] <- c("ns", "power_emp")
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_ns %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd(),
                            r1_r0 == cicc_ratio_est_swd(),
                            a0 == oicc_wperiod_est_swd(),
                            a1_a0 == oicc_ratio_est_swd()) %>%
              mutate(r1_r0=factor(r1_r0)) %>%
              plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                      linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                   "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                   "<br>Clusters (total):", ns*2, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                   "<br>Periods (J): ", J,
                                   "<br>HTE power: ", round(power_emp,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-period covariate ICC = ", cicc_wperiod_est_swd()),
                     #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Cluster size (per period)"),
                     yaxis=list(title="Number of clusters (per sequence)"),
                     legend=list(title=list(text="covariate ICC ratio")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size (per period) vs number of clusters (per sequence)',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), r0 (", cicc_wperiod_est_swd(),") and r1/r0 (", cicc_ratio_est_swd(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_swd == "oICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              
              a0_unique <- sort(unique(df_ns[,"a0"]))
              a1_a0_unique <- sort(unique(df_ns[,"a1_a0"]))
              r0_unique <- sort(unique(df_ns[,"r0"]))
              r1_r0_unique <- sort(unique(df_ns[,"r1_r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", ns*2, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(r1_r0=factor(r1_r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd(),
                                  a1_a0 == oicc_ratio_est_swd(),
                                  r0 == r0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r1_r0,
                            linetype=~r1_r0, color=~r1_r0,legendgroup=~r1_r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", ns*2, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period covariate ICC = ", r0_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs number of clusters (per sequence)',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), minimum r0 (", cicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), ") and a1/a0 (", oicc_ratio_est_swd(),"), and maximum r0 (", cicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd(), "), a1/a0 (", oicc_ratio_est_swd(),"), and r0 (", cicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_ns[,"a0"]))
              a1_a0_unique <- sort(unique(df_ns[,"a1_a0"]))
              r0_unique <- sort(unique(df_ns[,"r0"]))
              r1_r0_unique <- sort(unique(df_ns[,"r1_r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", ns*2, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           yaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd(),
                                  r1_r0 == cicc_ratio_est_swd(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,
                                         "; <br>within-period covariate ICC (r0): ", r0, "; <br>covariate CAC (r1/r0): ", r1_r0,
                                         "<br>Clusters (total):", ns*2, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           yaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs number of clusters (per sequence)',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), minimum a0 (", oicc_wperiod_min_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), ") and r1/r0 (", cicc_ratio_est_swd(),"), and maximum a0 (", oicc_wperiod_max_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd(), "), r1/r0 (", cicc_ratio_est_swd(),"), and a0 (", oicc_wperiod_est_swd(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }else if(cohort() == "closed"){
          if(sensitivity_swd_react() == "est_only"){
            df_ns <- expand.grid(power=input$power_swd, m=m_range, J=J,
                                 a0=#c(oicc_wperiod_min_swd(),
                                   oicc_wperiod_est_swd_cc(),
                                 # oicc_wperiod_max_swd()),
                                 a1_a0=#c(oicc_ratio_min_swd(),
                                   oicc_ratio_est_swd_cc(),
                                 # oicc_ratio_max_swd()),
                                 a2=#c(oicc_windiv_min_swd_cc(),
                                   oicc_windiv_est_swd_cc(),
                                 # oicc_windiv_max_swd_cc()),
                                 r0=#c(cicc_wperiod_min_swd_cc(),
                                   cicc_wperiod_est_swd_cc(),
                                 # cicc_wperiod_max_swd_cc()),
                                 var_y=var_y, #(input$sd_outcome)^2,
                                 var_x=var_x,#(input$sd_covar)^2,
                                 cohort=cohort(),
                                 d=input$mean_diff_HTE, a=input$sig)
            
          }else{
            df_ns <- expand.grid(power=input$power_swd, m=m_range, J=J,
                                 a0=c(oicc_wperiod_min_swd_cc(),
                                      oicc_wperiod_est_swd_cc(),
                                      oicc_wperiod_max_swd_cc()),
                                 a1_a0=c(oicc_ratio_min_swd_cc(),
                                         oicc_ratio_est_swd_cc(),
                                         oicc_ratio_max_swd_cc()),
                                 a2=c(oicc_windiv_min_swd_cc(),
                                      oicc_windiv_est_swd_cc(),
                                      oicc_windiv_max_swd_cc()),
                                 r0=c(cicc_wperiod_min_swd_cc(),
                                      cicc_wperiod_est_swd_cc(),
                                      cicc_wperiod_max_swd_cc()),
                                 var_y=var_y, #(input$sd_outcome)^2,
                                 var_x=var_x,#(input$sd_covar)^2,
                                 cohort=cohort(),
                                 d=input$mean_diff_HTE, a=input$sig)
            
          }
          
          
          ns_hte_col <- matrix(NA, nrow(df_ns), ncol=2)
          for(i in seq(nrow(df_ns))){
            ns_hte_col[i,] <- unlist(ns_hte_crossover(m=df_ns[i,"m"],
                                                      J=df_ns[i,"J"],
                                                      var_y=df_ns[i,"var_y"], var_x=df_ns[i,"var_x"],
                                                      a0=df_ns[i,"a0"], a1_a0=df_ns[i,"a1_a0"],
                                                      a2=df_ns[i,"a2"], r0=df_ns[i,"r0"],
                                                      d=df_ns[i,"d"],
                                                      cohort=df_ns[i,"cohort"],
                                                      a=df_ns[i,"a"],
                                                      power=df_ns[i,"power"]))
          }
          
          df_ns <- cbind(df_ns, ns_hte_col)
          colnames(df_ns)[(ncol(df_ns)-1):ncol(df_ns)] <- c("ns", "power_emp")
          
          if(sensitivity_swd_react() == "est_only"){
            p_est <- df_ns %>%
              as.data.frame() %>%
              dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                            a0 == oicc_wperiod_est_swd_cc(),
                            a1_a0 == oicc_ratio_est_swd_cc(),
                            a2 == oicc_windiv_est_swd_cc()) %>%
              mutate(r0=factor(r0)) %>%
              plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r0,
                      linetype=~r0, color=~r0,legendgroup=~r0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0, "; <br>within-individual outcome ICC (a2): ", a2,
                                   "; <br>within-period covariate ICC (r0): ", r0,
                                   "<br>Clusters (total):", ns*2, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                   "<br>Periods (J): ", J,
                                   "<br>HTE power: ", round(power_emp,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("within-individual outcome ICC = ", oicc_windiv_est_swd_cc()),
                     #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                     xaxis=list(title="Cluster size (per period)"),
                     yaxis=list(title="Number of clusters (per sequence)"),
                     legend=list(title=list(text="covariate ICC ratio")),
                     margin=0.001)
            
            subplot(p_est,# nrows=2, widths = c(0.5,0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size (per period) vs number of clusters (per sequence)',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.485,
                  y = 1.03,#0.97,
                  text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), a2 (", oicc_windiv_est_swd_cc(), "), and r0 (", cicc_wperiod_est_swd_cc(),")</i>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else if(sensitivity_swd_react() == "sensitivity"){
            if(input$icc_display_swd == "oICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              
              a0_unique <- sort(unique(df_ns[,"a0"]))
              a1_a0_unique <- sort(unique(df_ns[,"a1_a0"]))
              a2_unique <- sort(unique(df_ns[,"a2"]))
              r0_unique <- sort(unique(df_ns[,"r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(r0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_o[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", ns*2, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                  
                }else{
                  p_o[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(r0=factor(r0)) %>%
                    dplyr::filter(a0 == oicc_wperiod_est_swd_cc(),
                                  a1_a0 == oicc_ratio_est_swd_cc(),
                                  a2 == a2_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~r0,
                            linetype=~r0, color=~r0,legendgroup=~r0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", ns*2, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-individual outcome ICC = ", a2_unique[i]),
                           #subtitle=paste("within-period outcome ICC = ", oicc_wperiod_est_swd(), ", outcome ICC ratio = ", oicc_ratio_est_swd()),
                           xaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="covariate CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_o[[1]], p_o[[2]],p_o[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs number of clusters (per sequence)',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), minimum r0 (", cicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), and a2 (", oicc_windiv_est_swd_cc(),"), maximum r0 (", cicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed a0 (", oicc_wperiod_est_swd_cc(), "), a1/a0 (", oicc_ratio_est_swd_cc(),"), a2 (", oicc_windiv_est_swd_cc(),"), and r0 (", cicc_wperiod_est_swd_cc(),")</i>"),                      xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }else if(input$icc_display_swd == "cICC_constant"){
              #legend_title <- latex2exp::TeX("$\\rho_{y|x}$")
              a0_unique <- sort(unique(df_ns[,"a0"]))
              a1_a0_unique <- sort(unique(df_ns[,"a1_a0"]))
              a2_unique <- sort(unique(df_ns[,"a2"]))
              r0_unique <- sort(unique(df_ns[,"r0"]))
              p_o <- p_c <- vector(mode="list", length=length(a0_unique))
              
              for(i in seq(length(a0_unique))){
                
                if(i != 3){
                  # outcome ICCs #
                  p_c[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", ns*2, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           yaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                  
                }else{
                  p_c[[i]] <- df_ns %>%
                    as.data.frame() %>%
                    mutate(a1_a0=factor(a1_a0)) %>%
                    dplyr::filter(r0 == cicc_wperiod_est_swd_cc(),
                                  a2 == oicc_windiv_est_swd_cc(),
                                  a0 == a0_unique[i]) %>%
                    plot_ly(x=~m,y=~ns, type='scatter', mode='lines', line=list(width=3), name=~a1_a0,
                            linetype=~a1_a0, color=~a1_a0,legendgroup=~a1_a0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("within-period outcome ICC (a0): ", a0, "; <br>outcome CAC (a1/a0): ", a1_a0,"; <br>within-individual outcome ICC (a2): ", a2,
                                         "; <br>within-period covariate ICC (r0): ", r0,
                                         "<br>Clusters (total):", ns*2, "; Clusters (per sequence):", ns,"<br>Cluster size (per period): ", m,
                                         "<br>Periods (J): ", J,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("within-period outcome ICC = ", a0_unique[i]),
                           #subtitle=paste("within-period covariate ICC = ", cicc_wperiod_est_swd(), ", covariate ICC ratio = ", cicc_ratio_est_swd()),
                           yaxis=list(title="Cluster size (per period)"),
                           yaxis=list(title="Number of clusters (per sequence)"),
                           legend=list(title=list(text="outcome CAC")),
                           margin=0.01)
                }
                
              }
              
              subplot(p_c[[1]], p_c[[2]],p_c[[3]], nrows=2, widths = c(0.5,0.5),
                      margin = 0.09, titleX=T, titleY=T) %>% #, list(b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size (per period) vs number of clusters (per sequence)',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    #yshift=-30,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), minimum a0 (", oicc_wperiod_min_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE#,
                    #font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), ") and a2 (", oicc_windiv_est_swd_cc(),"), and maximum a0 (", oicc_wperiod_max_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste0("<i>Assumed r0 (", cicc_wperiod_est_swd_cc(), "), a2 (", oicc_windiv_est_swd_cc(),"), and a0 (", oicc_wperiod_est_swd_cc(),")</i>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }
          }# end sensitivity if/else
          
          
        }#end cross/closed
        
        
      }#end fixed power
      
      #### IRGT ####
    }else if( trial_react() == "irgt"){
      if(clustering_irgt() == "indiv"){
        m0_fix <- 1
      }
      
      # determine effect sizes and variances depending on outcome/covariate type #
      if(outcome_type_react() =="continuous"){
        var_y1 <- (sd_outcome1_react())^2
        var_y0 <- (sd_outcome0_react())^2
        d <- input$mean_diff_HTE
        
        if(covar_type_react()=="continuous"){
          # continuous outcome and covar power
          var_x <- (sd_covar_react())^2
          
        }else if(covar_type_react() == "binary"){
          # continuous outcome binary covar power
          var_x <- (prop_covar_react())*(1-prop_covar_react())
        }
        
      }
      
      if(plot_display_react() == "m_v_power"){
        m1_range <- seq(input$m1_slide[1],input$m1_slide[2])
        m1_fix <- input$m1_fix
        
        if(clustering_irgt() == "cluster"){
          m0_range <- seq(input$m0_slide[1],input$m0_slide[2])
          m0_fix <- input$m0_fix
        }
        
        if(sensitivity_irgt_react() == "est_only"){
          # m1 on x-axis #
          df1_power <- expand.grid(n1=input$n1_fix, m1=m1_range,
                                   n0=input$n0_fix, m0=m0_fix,
                                   oicc1=oicc_trt_est_irgt(),
                                   oicc0=oicc_ctrl_est_irgt(),
                                   var_y1=var_y1,
                                   var_y0=var_y0,
                                   var_x=var_x,
                                   d=input$mean_diff_HTE, a=input$sig)
          # print(m0_fix)
          
          if(clustering_irgt() == "cluster"){
            # m0 on x-axis #
            df0_power <- expand.grid(n1=input$n1_fix, m1=m1_fix,
                                     n0=input$n0_fix, m0=m0_range,
                                     oicc1=oicc_trt_est_irgt(),
                                     oicc0=oicc_ctrl_est_irgt(),
                                     var_y1=var_y1,
                                     var_y0=var_y0,
                                     var_x=var_x,
                                     d=input$mean_diff_HTE, a=input$sig)
          }
          
        }else{
          # m1 on x-axis #
          df1_power <- expand.grid(n1=input$n1_fix, m1=m1_range,
                                   n0=input$n0_fix, m0=m0_fix,
                                   oicc1=c(oicc_trt_min_irgt(), oicc_trt_est_irgt(),
                                           oicc_trt_max_irgt()),
                                   oicc0=c(oicc_ctrl_min_irgt(), oicc_ctrl_est_irgt(),
                                           oicc_ctrl_max_irgt()),
                                   var_y1=var_y1,
                                   var_y0=var_y0,
                                   var_x=var_x,
                                   d=input$mean_diff_HTE, a=input$sig)
          if(clustering_irgt() == "cluster"){
            # m0 on x-axis #
            df0_power <- expand.grid(n1=input$n1_fix, m1=m1_fix,
                                     n0=input$n0_fix, m0=m0_range,
                                     oicc1=c(oicc_trt_min_irgt(), oicc_trt_est_irgt(),
                                             oicc_trt_max_irgt()),
                                     oicc0=c(oicc_ctrl_min_irgt(), oicc_ctrl_est_irgt(),
                                             oicc_ctrl_max_irgt()),
                                     var_y1=var_y1,
                                     var_y0=var_y0,
                                     var_x=var_x,
                                     d=input$mean_diff_HTE, a=input$sig)
          }
        }
        
        power1_hte_col <- rep(NA, nrow(df1_power))
        
        for(i in seq(nrow(df1_power))){
          power1_hte_col[i] <- power_irgt(m1=df1_power[i,"m1"],
                                          m0=df1_power[i,"m0"],
                                          n1=df1_power[i,"n1"],
                                          n0=df1_power[i,"n0"],
                                          oicc1=df1_power[i,"oicc1"],
                                          oicc0=df1_power[i,"oicc0"],
                                          var_y1=df1_power[i,"var_y1"],
                                          var_y0=df1_power[i,"var_y0"],
                                          var_x=df1_power[i,"var_x"],
                                          d=df1_power[i,"d"], a=df1_power[i,"a"])
          
        }
        
        if(clustering_irgt() == "cluster"){
          power0_hte_col <- rep(NA, nrow(df0_power))
          
          for(i in seq(nrow(df0_power))){
            power0_hte_col[i] <- power_irgt(m1=df0_power[i,"m1"],
                                            m0=df0_power[i,"m0"],
                                            n1=df0_power[i,"n1"],
                                            n0=df0_power[i,"n0"],
                                            oicc1=df0_power[i,"oicc1"],
                                            oicc0=df0_power[i,"oicc0"],
                                            var_y1=df0_power[i,"var_y1"],
                                            var_y0=df0_power[i,"var_y0"],
                                            var_x=df0_power[i,"var_x"],
                                            d=df0_power[i,"d"], a=df0_power[i,"a"])
          }
          df0_power <- cbind(df0_power, power0_hte_col)
          
        }
        
        df1_power <- cbind(df1_power, power1_hte_col)
        
        if(sensitivity_irgt_react() == "est_only"){
          # m1 on x-axis #
          p1_est <- df1_power %>%
            as.data.frame() %>%
            dplyr::filter(oicc1 == oicc_trt_est_irgt(),#oicc_unique[2],
                          oicc0 == oicc_ctrl_est_irgt()#cicc_unique[2]
            ) %>%
            mutate(oicc0=factor(oicc0)) %>%
            plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                    linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("Treatment outcome ICC: ", oicc1,
                                 "<br>Control outcome ICC: ", oicc0,
                                 "<br>Treatment clusters (n1):", n1,
                                 "<br>Treatment cluster size (m1): ", m1,
                                 "<br>Control clusters (n0):", n0,
                                 "<br>Control cluster size (m0): ", m0,
                                 "<br>HTE power: ", round(power1_hte_col,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("o-ICC = ", oicc_trt_est_irgt()),
                   xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                         standoff=10)),
                   yaxis=list(title="HTE Power"),
                   legend=list(title=list(text="Control-arm outcome ICC")),
                   margin=0.001)
          
          if(clustering_irgt() == "cluster"){
            
            # m0 on x-axis #
            p0_est <- df0_power %>%
              as.data.frame() %>%
              dplyr::filter(oicc1 == oicc_trt_est_irgt(),#oicc_unique[2],
                            oicc0 == oicc_ctrl_est_irgt()#cicc_unique[2]
              ) %>%
              mutate(oicc0=factor(oicc0)) %>%
              plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                      linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("Treatment outcome ICC: ", oicc1,
                                   "<br>Control outcome ICC: ", oicc0,
                                   "<br>Treatment clusters (n1):", n1,
                                   "<br>Treatment cluster size (m1): ", m1,
                                   "<br>Control clusters (n0):", n0,
                                   "<br>Control cluster size (m0): ", m0,
                                   "<br>HTE power: ", round(power0_hte_col,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("o-ICC = ", oicc_est()),
                     xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                           standoff=10)),
                     yaxis=list(title="HTE Power"),
                     legend=list(title=list(text="Control-arm outcome ICC")),
                     margin=0.001)
            
            subplot(p1_est, p0_est,#nrows=1, widths = c(0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size vs HTE power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.0,
                  text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Assumed treatment-arm outcome ICC (", oicc_trt_est_irgt(),") and control-arm outcome ICC (", oicc_ctrl_est_irgt(),")</i>"),
                                       width=0.8*getOption("width")/2), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.0,
                  text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Assumed treatment-arm outcome ICC (", oicc_trt_est_irgt(),") and control-arm outcome ICC (", oicc_ctrl_est_irgt(),")</i>"),
                                       width=0.8*getOption("width")/2), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else{
            subplot(p1_est, #p0_est,#nrows=1, widths = c(0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Treatment-arm cluster size vs HTE power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.5,
                  y = 1.0,
                  text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc_trt_est_irgt(),") and control-arm outcome ICC (", oicc_ctrl_est_irgt(),")</i>"),
                                       width=0.75*getOption("width")), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }# end clustering if/else
          
        }else if(sensitivity_irgt_react() == "sensitivity"){
          if(input$icc_display_irgt == "oICC1_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_x$")
            oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
            oicc0_unique <- sort(unique(df1_power[,"oicc0"]))
            p1 <- p0 <- vector(mode="list", length=length(oicc1_unique))
            
            for(i in seq(length(oicc1_unique))){
              
              if(i != 3){
                # m1 on x-axis #
                p1[[i]] <- df1_power %>%
                  as.data.frame() %>%
                  mutate(oicc0=factor(oicc0)) %>%
                  dplyr::filter(oicc1 == oicc1_unique[i]) %>%
                  plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                          linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("Treatment outcome ICC: ", oicc1,
                                       "<br>Control outcome ICC: ", oicc0,
                                       "<br>Treatment clusters (n1):", n1,
                                       "<br>Treatment cluster size (m1): ", m1,
                                       "<br>Control clusters (n0):", n0,
                                       "<br>Control cluster size (m0): ", m0,
                                       "<br>HTE power: ", round(power1_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                         xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                               standoff=10)),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="Control-arm outcome ICC")),
                         margin=0.01)
                
                if(clustering_irgt() == "cluster"){
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i]) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.01)
                }
              }else{
                
                if(clustering_irgt() == "cluster"){
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i]) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i]) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                }else{
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i]) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                }
              }
              
            }# end row loop
            
            
            
            if(clustering_irgt() == "cluster"){
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.6,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.27,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.6,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.27,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }else{
              # m1 on x-axis #
              subplot(p1[[1]],
                      p1[[2]],
                      p1[[3]],
                      nrows=2,# heights = c(0.33,0.33,0.33),
                      margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Treatment-arm cluster size vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    text = paste(strwrap(paste0("<i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }# end clustering if/else
          }else if(input$icc_display_irgt == "oICC0_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_x$")
            oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
            oicc0_unique <- sort(unique(df1_power[,"oicc0"]))
            p1 <- p0 <- vector(mode="list", length=length(oicc0_unique))
            
            for(i in seq(length(oicc1_unique))){
              
              if(i != 3){
                # m1 on x-axis #
                p1[[i]] <- df1_power %>%
                  as.data.frame() %>%
                  mutate(oicc1=factor(oicc1)) %>%
                  dplyr::filter(oicc0 == oicc0_unique[i]) %>%
                  plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                          linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("Treatment outcome ICC: ", oicc1,
                                       "<br>Control outcome ICC: ", oicc0,
                                       "<br>Treatment clusters (n1):", n1,
                                       "<br>Treatment cluster size (m1): ", m1,
                                       "<br>Control clusters (n0):", n0,
                                       "<br>Control cluster size (m0): ", m0,
                                       "<br>HTE power: ", round(power1_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                         xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                               standoff=10)),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="Treatment-arm outcome ICC")),
                         margin=0.01)
                if(clustering_irgt() == "cluster"){
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i]) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.01)
                }
              }else{
                
                if(clustering_irgt() == "cluster"){
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i]) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i]) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                }else{
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i]) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                }
              }
              
            }# end row loop
            
            if(clustering_irgt() == "cluster"){
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.6,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.27,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.6,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.27,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }else{
              # m1 on x-axis #
              subplot(p1[[1]],#p0[[1]],
                      p1[[2]],#p0[[2]],
                      p1[[3]],#p0[[3]],
                      nrows=2,# heights = c(0.33,0.33,0.33),
                      margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Treatment-arm cluster size vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    text = paste(strwrap(paste0("<i>Minimum control-arm outcome ICC (", oicc0_unique[1],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
            }# end clustering if/else
            
          }# oicc1 constant end here
        }# end sensitivity
        
      }else if(plot_display_react() == "n_v_power"){
        
        n1_range <- seq(input$n1_slide[1],input$n1_slide[2])
        n1_fix <- input$n1_fix
        
        n0_range <- seq(input$m0_slide[1],input$m0_slide[2])
        n0_fix <- input$m0_fix
        
        if(sensitivity_irgt_react() == "est_only"){
          # n1 on x-axis #
          df1_power <- expand.grid(n1=n1_range, m1=input$m1_fix,
                                   n0=input$n0_fix, m0=input$m0_fix,
                                   oicc1=oicc_trt_est_irgt(),
                                   oicc0=oicc_ctrl_est_irgt(),
                                   var_y1=var_y1,
                                   var_y0=var_y0,
                                   var_x=var_x,
                                   d=input$mean_diff_HTE, a=input$sig)
          
          # n0 on x-axis #
          df0_power <- expand.grid(n1=n1_fix, m1=input$m1_fix,
                                   n0=n0_range, m0=input$m0_fix,
                                   oicc1=oicc_trt_est_irgt(),
                                   oicc0=oicc_ctrl_est_irgt(),
                                   var_y1=var_y1,
                                   var_y0=var_y0,
                                   var_x=var_x,
                                   d=input$mean_diff_HTE, a=input$sig)
          
        }else{
          # n1 on x-axis #
          df1_power <- expand.grid(n1=n1_range, m1=input$m1_fix,
                                   n0=n0_fix, m0=input$m0_fix,
                                   oicc1=c(oicc_trt_min_irgt(), oicc_trt_est_irgt(),
                                           oicc_trt_max_irgt()),
                                   oicc0=c(oicc_ctrl_min_irgt(), oicc_ctrl_est_irgt(),
                                           oicc_ctrl_max_irgt()),
                                   var_y1=var_y1,
                                   var_y0=var_y0,
                                   var_x=var_x,
                                   d=input$mean_diff_HTE, a=input$sig)
          
          # n0 on x-axis #
          df0_power <- expand.grid(n1=n1_fix, m1=input$m1_fix,
                                   n0=n0_range, m0=input$m0_fix,
                                   oicc1=c(oicc_trt_min_irgt(), oicc_trt_est_irgt(),
                                           oicc_trt_max_irgt()),
                                   oicc0=c(oicc_ctrl_min_irgt(), oicc_ctrl_est_irgt(),
                                           oicc_ctrl_max_irgt()),
                                   var_y1=var_y1,
                                   var_y0=var_y0,
                                   var_x=var_x,
                                   d=input$mean_diff_HTE, a=input$sig)
        }
        
        power1_hte_col <- rep(NA, nrow(df1_power))
        
        for(i in seq(nrow(df1_power))){
          power1_hte_col[i] <- power_irgt(m1=df1_power[i,"m1"],
                                          m0=df1_power[i,"m0"],
                                          n1=df1_power[i,"n1"],
                                          n0=df1_power[i,"n0"],
                                          oicc1=df1_power[i,"oicc1"],
                                          oicc0=df1_power[i,"oicc0"],
                                          var_y1=df1_power[i,"var_y1"],
                                          var_y0=df1_power[i,"var_y0"],
                                          var_x=df1_power[i,"var_x"],
                                          d=df1_power[i,"d"], a=df1_power[i,"a"])
          
        }
        
        power0_hte_col <- rep(NA, nrow(df0_power))
        
        for(i in seq(nrow(df0_power))){
          power0_hte_col[i] <- power_irgt(m1=df0_power[i,"m1"],
                                          m0=df0_power[i,"m0"],
                                          n1=df0_power[i,"n1"],
                                          n0=df0_power[i,"n0"],
                                          oicc1=df0_power[i,"oicc1"],
                                          oicc0=df0_power[i,"oicc0"],
                                          var_y1=df0_power[i,"var_y1"],
                                          var_y0=df0_power[i,"var_y0"],
                                          var_x=df0_power[i,"var_x"],
                                          d=df0_power[i,"d"], a=df0_power[i,"a"])
        }
        df0_power <- cbind(df0_power, power0_hte_col)
        
        df1_power <- cbind(df1_power, power1_hte_col)
        
        if(sensitivity_irgt_react() == "est_only"){
          # n1 on x-axis #
          p1_est <- df1_power %>%
            as.data.frame() %>%
            dplyr::filter(oicc1 == oicc_trt_est_irgt(),#oicc_unique[2],
                          oicc0 == oicc_ctrl_est_irgt()#cicc_unique[2]
            ) %>%
            mutate(oicc0=factor(oicc0)) %>%
            plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                    linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("Treatment outcome ICC: ", oicc1,
                                 "<br>Control outcome ICC: ", oicc0,
                                 "<br>Treatment clusters (n1):", n1,
                                 "<br>Treatment cluster size (m1): ", m1,
                                 "<br>Control clusters (n0):", n0,
                                 "<br>Control cluster size (m0): ", m0,
                                 "<br>HTE power: ", round(power1_hte_col,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("o-ICC = ", oicc_trt_est_irgt()),
                   xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                         standoff=10)),
                   yaxis=list(title="HTE Power"),
                   legend=list(title=list(text="Control-arm outcome ICC")),
                   margin=0.001)
          
          # n0 on x-axis #
          p0_est <- df0_power %>%
            as.data.frame() %>%
            dplyr::filter(oicc1 == oicc_trt_est_irgt(),#oicc_unique[2],
                          oicc0 == oicc_ctrl_est_irgt()#cicc_unique[2]
            ) %>%
            mutate(oicc0=factor(oicc0)) %>%
            plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                    linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("Treatment outcome ICC: ", oicc1,
                                 "<br>Control outcome ICC: ", oicc0,
                                 "<br>Treatment clusters (n1):", n1,
                                 "<br>Treatment cluster size (m1): ", m1,
                                 "<br>Control clusters (n0):", n0,
                                 "<br>Control cluster size (m0): ", m0,
                                 "<br>HTE power: ", round(power0_hte_col,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("o-ICC = ", oicc_est()),
                   xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                         standoff=10)),
                   yaxis=list(title="HTE Power"),
                   legend=list(title=list(text="Control-arm outcome ICC")),
                   margin=0.001)
          
          subplot(p1_est, p0_est,#nrows=1, widths = c(0.5),
                  margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
            layout(title = list(
              text='Cluster size vs HTE power',
              font=list(size=17)
            ),
            #margin=list(pad=50),
            annotations = list(
              list(
                x = 0.23,
                y = 1.0,
                text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Assumed treatment-arm outcome ICC (", oicc_trt_est_irgt(),") and control-arm outcome ICC (", oicc_ctrl_est_irgt(),")</i>"),
                                     width=0.8*getOption("width")/2), collapse="<br>"),
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = FALSE, font=list(size=13)
              ),
              list(
                x = 0.77,
                y = 1.0,
                text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Assumed treatment-arm outcome ICC (", oicc_trt_est_irgt(),") and control-arm outcome ICC (", oicc_ctrl_est_irgt(),")</i>"),
                                     width=0.8*getOption("width")/2), collapse="<br>"),
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = FALSE, font=list(size=13)
              )
            ),
            legend=list(orientation="h",
                        yanchor="center",
                        y=0.25,
                        x=0.5)
            )
          
        }else if(sensitivity_irgt_react() == "sensitivity"){
          if(input$icc_display_irgt == "oICC1_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_x$")
            oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
            oicc0_unique <- sort(unique(df1_power[,"oicc0"]))
            p1 <- p0 <- vector(mode="list", length=length(oicc1_unique))
            
            for(i in seq(length(oicc1_unique))){
              
              if(i != 3){
                # n1 on x-axis #
                p1[[i]] <- df1_power %>%
                  as.data.frame() %>%
                  mutate(oicc0=factor(oicc0)) %>%
                  dplyr::filter(oicc1 == oicc1_unique[i]) %>%
                  plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                          linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("Treatment outcome ICC: ", oicc1,
                                       "<br>Control outcome ICC: ", oicc0,
                                       "<br>Treatment clusters (n1):", n1,
                                       "<br>Treatment cluster size (m1): ", m1,
                                       "<br>Control clusters (n0):", n0,
                                       "<br>Control cluster size (m0): ", m0,
                                       "<br>HTE power: ", round(power1_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                         xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                               standoff=10)),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="Control-arm outcome ICC")),
                         margin=0.01)
                
                # n0 on x-axis #
                p0[[i]] <- df0_power %>%
                  as.data.frame() %>%
                  mutate(oicc0=factor(oicc0)) %>%
                  dplyr::filter(oicc1 == oicc1_unique[i]) %>%
                  plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                          linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("Treatment outcome ICC: ", oicc1,
                                       "<br>Control outcome ICC: ", oicc0,
                                       "<br>Treatment clusters (n1):", n1,
                                       "<br>Treatment cluster size (m1): ", m1,
                                       "<br>Control clusters (n0):", n0,
                                       "<br>Control cluster size (m0): ", m0,
                                       "<br>HTE power: ", round(power0_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                         xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                               standoff=10)),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="Control-arm outcome ICC")),
                         margin=0.01)
              }else{
                
                # n1 on x-axis #
                p1[[i]] <- df1_power %>%
                  as.data.frame() %>%
                  mutate(oicc0=factor(oicc0)) %>%
                  dplyr::filter(oicc1 == oicc1_unique[i]) %>%
                  plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                          linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("Treatment outcome ICC: ", oicc1,
                                       "<br>Control outcome ICC: ", oicc0,
                                       "<br>Treatment clusters (n1):", n1,
                                       "<br>Treatment cluster size (m1): ", m1,
                                       "<br>Control clusters (n0):", n0,
                                       "<br>Control cluster size (m0): ", m0,
                                       "<br>HTE power: ", round(power1_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                         xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                               standoff=10)),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="Control-arm outcome ICC")),
                         margin=0.001)
                
                # n0 on x-axis #
                p0[[i]] <- df0_power %>%
                  as.data.frame() %>%
                  mutate(oicc0=factor(oicc0)) %>%
                  dplyr::filter(oicc1 == oicc1_unique[i]) %>%
                  plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                          linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("Treatment outcome ICC: ", oicc1,
                                       "<br>Control outcome ICC: ", oicc0,
                                       "<br>Treatment clusters (n1):", n1,
                                       "<br>Treatment cluster size (m1): ", m1,
                                       "<br>Control clusters (n0):", n0,
                                       "<br>Control cluster size (m0): ", m0,
                                       "<br>HTE power: ", round(power0_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                         xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                               standoff=10)),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="Control-arm outcome ICC")),
                         margin=0.001)
                
              }
              
            }# end row loop
            
            # n1 on x-axis #
            subplot(p1[[1]],p0[[1]],
                    p1[[2]],p0[[2]],
                    p1[[3]],p0[[3]],
                    nrows=3,# heights = c(0.33,0.33,0.33),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size vs HTE power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.0,
                  text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],")</i>"),
                                       width=0.75*getOption("width")), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE,
                  font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.6,
                  text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],")</i>"),
                                       width=0.75*getOption("width")), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.27,
                  text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],")</i>"),
                                       width=0.75*getOption("width")), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.0,
                  text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],")</i>"),
                                       width=0.75*getOption("width")), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE,
                  font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 0.6,
                  text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],")</i>"),
                                       width=0.75*getOption("width")), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 0.27,
                  text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],")</i>"),
                                       width=0.75*getOption("width")), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=-0.1,
                          x=0.6)
              )
            
          }else if(input$icc_display_irgt == "oICC0_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_x$")
            oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
            oicc0_unique <- sort(unique(df1_power[,"oicc0"]))
            p1 <- p0 <- vector(mode="list", length=length(oicc0_unique))
            
            for(i in seq(length(oicc1_unique))){
              
              if(i != 3){
                # n1 on x-axis #
                p1[[i]] <- df1_power %>%
                  as.data.frame() %>%
                  mutate(oicc1=factor(oicc1)) %>%
                  dplyr::filter(oicc0 == oicc0_unique[i]) %>%
                  plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                          linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("Treatment outcome ICC: ", oicc1,
                                       "<br>Control outcome ICC: ", oicc0,
                                       "<br>Treatment clusters (n1):", n1,
                                       "<br>Treatment cluster size (m1): ", m1,
                                       "<br>Control clusters (n0):", n0,
                                       "<br>Control cluster size (m0): ", m0,
                                       "<br>HTE power: ", round(power1_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                         xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                               standoff=10)),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="Treatment-arm outcome ICC")),
                         margin=0.01)
                
                # n0 on x-axis #
                p0[[i]] <- df0_power %>%
                  as.data.frame() %>%
                  mutate(oicc1=factor(oicc1)) %>%
                  dplyr::filter(oicc0 == oicc0_unique[i]) %>%
                  plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                          linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("Treatment outcome ICC: ", oicc1,
                                       "<br>Control outcome ICC: ", oicc0,
                                       "<br>Treatment clusters (n1):", n1,
                                       "<br>Treatment cluster size (m1): ", m1,
                                       "<br>Control clusters (n0):", n0,
                                       "<br>Control cluster size (m0): ", m0,
                                       "<br>HTE power: ", round(power0_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                         xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                               standoff=10)),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="Treatment-arm outcome ICC")),
                         margin=0.01)
                
              }else{
                
                # n1 on x-axis #
                p1[[i]] <- df1_power %>%
                  as.data.frame() %>%
                  mutate(oicc1=factor(oicc1)) %>%
                  dplyr::filter(oicc0 == oicc0_unique[i]) %>%
                  plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                          linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("Treatment outcome ICC: ", oicc1,
                                       "<br>Control outcome ICC: ", oicc0,
                                       "<br>Treatment clusters (n1):", n1,
                                       "<br>Treatment cluster size (m1): ", m1,
                                       "<br>Control clusters (n0):", n0,
                                       "<br>Control cluster size (m0): ", m0,
                                       "<br>HTE power: ", round(power1_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                         xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                               standoff=10)),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="Treatment-arm outcome ICC")),
                         margin=0.001)
                
                # n0 on x-axis #
                p0[[i]] <- df0_power %>%
                  as.data.frame() %>%
                  mutate(oicc1=factor(oicc1)) %>%
                  dplyr::filter(oicc0 == oicc0_unique[i]) %>%
                  plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                          linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("Treatment outcome ICC: ", oicc1,
                                       "<br>Control outcome ICC: ", oicc0,
                                       "<br>Treatment clusters (n1):", n1,
                                       "<br>Treatment cluster size (m1): ", m1,
                                       "<br>Control clusters (n0):", n0,
                                       "<br>Control cluster size (m0): ", m0,
                                       "<br>HTE power: ", round(power0_hte_col,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                         xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                               standoff=10)),
                         yaxis=list(title="HTE Power"),
                         legend=list(title=list(text="Treatment-arm outcome ICC")),
                         margin=0.001)
              }
              
            }# end row loop
            
            # n1 on x-axis #
            subplot(p1[[1]],p0[[1]],
                    p1[[2]],p0[[2]],
                    p1[[3]],p0[[3]],
                    nrows=3,# heights = c(0.33,0.33,0.33),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size vs HTE power',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.0,
                  text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],")</i>"),
                                       width=0.75*getOption("width")), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE,
                  font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.6,
                  text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],")</i>"),
                                       width=0.75*getOption("width")), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.23,
                  y = 0.27,
                  text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],")</i>"),
                                       width=0.75*getOption("width")), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.0,
                  text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],")</i>"),
                                       width=0.75*getOption("width")), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE,
                  font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 0.6,
                  text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],")</i>"),
                                       width=0.75*getOption("width")), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 0.27,
                  text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],")</i>"),
                                       width=0.75*getOption("width")), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=-0.1,
                          x=0.6)
              )
            
          }# oicc1 constant end here
        }# end sensitivity
        
      }else if(plot_display_react() == "fixed_power"){
        ## NEED TO DO THIS WHOLE THING FOR POWER ##
        m1_range <- seq(input$m1_slide[1],input$m1_slide[2])
        m1_fix <- input$m1_fix
        
        if(clustering_irgt() == "cluster"){
          m0_range <- seq(input$m0_slide[1],input$m0_slide[2])
          m0_fix <- input$m0_fix
        }
        
        if(sensitivity_irgt_react() == "est_only"){
          # m1 on x-axis #
          df1_n <- expand.grid(#n1=input$n1_fix, 
            m1=m1_range,
            #n0=input$n0_fix,
            m0=m0_fix,
            power=input$power_irgt,
            oicc1=oicc_trt_est_irgt(),
            oicc0=oicc_ctrl_est_irgt(),
            var_y1=var_y1,
            var_y0=var_y0,
            var_x=var_x,
            d=input$mean_diff_HTE, a=input$sig)
          # print(m0_fix)
          
          if(clustering_irgt() == "cluster"){
            # m0 on x-axis #
            df0_n <- expand.grid(#n1=input$n1_fix,
              m1=m1_fix,
              #n0=input$n0_fix, 
              m0=m0_range,
              power=input$power_irgt,
              oicc1=oicc_trt_est_irgt(),
              oicc0=oicc_ctrl_est_irgt(),
              var_y1=var_y1,
              var_y0=var_y0,
              var_x=var_x,
              d=input$mean_diff_HTE, a=input$sig)
          }
          
        }else{
          # m1 on x-axis #
          df1_n <- expand.grid(#n1=input$n1_fix, 
            m1=m1_range,
            #n0=input$n0_fix, 
            m0=m0_fix,
            power=input$power_irgt,
            oicc1=c(oicc_trt_min_irgt(), oicc_trt_est_irgt(),
                    oicc_trt_max_irgt()),
            oicc0=c(oicc_ctrl_min_irgt(), oicc_ctrl_est_irgt(),
                    oicc_ctrl_max_irgt()),
            var_y1=var_y1,
            var_y0=var_y0,
            var_x=var_x,
            d=input$mean_diff_HTE, a=input$sig)
          if(clustering_irgt() == "cluster"){
            # m0 on x-axis #
            df0_n <- expand.grid(#n1=input$n1_fix, 
              m1=m1_fix,
              #n0=input$n0_fix, 
              m0=m0_range,
              power=input$power_irgt,
              oicc1=c(oicc_trt_min_irgt(), oicc_trt_est_irgt(),
                      oicc_trt_max_irgt()),
              oicc0=c(oicc_ctrl_min_irgt(), oicc_ctrl_est_irgt(),
                      oicc_ctrl_max_irgt()),
              var_y1=var_y1,
              var_y0=var_y0,
              var_x=var_x,
              d=input$mean_diff_HTE, a=input$sig)
          }
        }
        
        n1_hte_col <- matrix(NA, ncol=5, nrow=nrow(df1_n))
        
        for(i in seq(nrow(df1_n))){
          n1_hte_col[i,] <- unlist( n_irgt(m1=df1_n[i,"m1"],
                                           m0=df1_n[i,"m0"],
                                           #n1=df1_n[i,"n1"],
                                           #n0=df1_n[i,"n0"],
                                           oicc1=df1_n[i,"oicc1"],
                                           oicc0=df1_n[i,"oicc0"],
                                           var_y1=df1_n[i,"var_y1"],
                                           var_y0=df1_n[i,"var_y0"],
                                           var_x=df1_n[i,"var_x"],
                                           d=df1_n[i,"d"], a=df1_n[i,"a"],
                                           power=df1_n[i,"power"]) )
          
        }
        
        if(clustering_irgt() == "cluster"){
          n0_hte_col <- matrix(NA, ncol=5, nrow=nrow(df0_n))
          
          for(i in seq(nrow(df0_n))){
            n0_hte_col[i,] <- unlist( n_irgt(m1=df0_n[i,"m1"],
                                             m0=df0_n[i,"m0"],
                                             #n1=df0_n[i,"n1"],
                                             #n0=df0_n[i,"n0"],
                                             oicc1=df0_n[i,"oicc1"],
                                             oicc0=df0_n[i,"oicc0"],
                                             var_y1=df0_n[i,"var_y1"],
                                             var_y0=df0_n[i,"var_y0"],
                                             var_x=df0_n[i,"var_x"],
                                             d=df0_n[i,"d"], a=df0_n[i,"a"],
                                             power=df0_n[i,"power"]) )
          }
          colnames(n1_hte_col) <- c("n", "n1", "n0","pi","power_emp")
          
          df0_n <- cbind(df0_n, n0_hte_col)
          colnames(df0_n)[11:15] <- c("n", "n1", "n0","pi","power_emp")
          
        }
        
        df1_n <- cbind(df1_n, n1_hte_col)
        colnames(df1_n)[11:15] <- c("n", "n1", "n0","pi","power_emp")
        
        ## do a line break in the y-axis title ##
        if(sensitivity_irgt_react() == "est_only"){
          # m1 on x-axis, n1 on y-axis #
          p1_est <- df1_n %>%
            as.data.frame() %>%
            dplyr::filter(oicc1 == oicc_trt_est_irgt(),#oicc_unique[2],
                          oicc0 == oicc_ctrl_est_irgt()#cicc_unique[2]
            ) %>%
            mutate(oicc0=factor(oicc0)) %>%
            plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                    linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("Treatment outcome ICC: ", oicc1,
                                 "<br>Control outcome ICC: ", oicc0,
                                 "<br>Treatment clusters (n1):", n1,
                                 "<br>Treatment cluster size (m1): ", m1,
                                 "<br>Control clusters (n0):", n0,
                                 "<br>Control cluster size (m0): ", m0,
                                 "<br>HTE power: ", round(power_emp,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("o-ICC = ", oicc_trt_est_irgt()),
                   xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                         standoff=10)),
                   yaxis=list(title="Treatment-arm\n number of clusters (n1)"),
                   legend=list(title=list(text="Control-arm outcome ICC")),
                   margin=0.001)
          
          if(clustering_irgt() == "cluster"){
            
            # m0 on x-axis #
            p0_est <- df0_n %>%
              as.data.frame() %>%
              dplyr::filter(oicc1 == oicc_trt_est_irgt(),#oicc_unique[2],
                            oicc0 == oicc_ctrl_est_irgt()#cicc_unique[2]
              ) %>%
              mutate(oicc0=factor(oicc0)) %>%
              plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                      linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                      colors=colors_plot,
                      hoverinfo="text",
                      text=~paste0("Treatment outcome ICC: ", oicc1,
                                   "<br>Control outcome ICC: ", oicc0,
                                   "<br>Treatment clusters (n1):", n1,
                                   "<br>Treatment cluster size (m1): ", m1,
                                   "<br>Control clusters (n0):", n0,
                                   "<br>Control cluster size (m0): ", m0,
                                   "<br>HTE power: ", round(power_emp,4)),
                      height=input$dimension[2]*0.8) %>%
              layout(title=paste("o-ICC = ", oicc_est()),
                     xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                           standoff=10)),
                     yaxis=list(title="Control-arm\n number of clusters (n0)"),
                     legend=list(title=list(text="Control-arm outcome ICC")),
                     margin=0.001)
            
            subplot(p1_est, p0_est,#nrows=1, widths = c(0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Cluster size vs number of clusters',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.23,
                  y = 1.0,
                  text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Assumed treatment-arm outcome ICC (", oicc_trt_est_irgt(),") and control-arm outcome ICC (", oicc_ctrl_est_irgt(),")</i>"),
                                       width=0.8*getOption("width")/2), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                ),
                list(
                  x = 0.77,
                  y = 1.0,
                  text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Assumed treatment-arm outcome ICC (", oicc_trt_est_irgt(),") and control-arm outcome ICC (", oicc_ctrl_est_irgt(),")</i>"),
                                       width=0.8*getOption("width")/2), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }else{
            subplot(p1_est, #p0_est,#nrows=1, widths = c(0.5),
                    margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
              layout(title = list(
                text='Treatment-arm cluster size vs number of clusters',
                font=list(size=17)
              ),
              #margin=list(pad=50),
              annotations = list(
                list(
                  x = 0.5,
                  y = 1.0,
                  text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc_trt_est_irgt(),") and control-arm outcome ICC (", oicc_ctrl_est_irgt(),")</i>"),
                                       width=0.75*getOption("width")), collapse="<br>"),
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "bottom",
                  showarrow = FALSE, font=list(size=13)
                )
              ),
              legend=list(orientation="h",
                          yanchor="center",
                          y=0.25,
                          x=0.5)
              )
          }# end clustering if/else
          
        }else if(sensitivity_irgt_react() == "sensitivity"){
          if(input$icc_display_irgt == "oICC1_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_x$")
            oicc1_unique <- sort(unique(df1_n[,"oicc1"]))
            oicc0_unique <- sort(unique(df1_n[,"oicc0"]))
            p1 <- p0 <- vector(mode="list", length=length(oicc1_unique))
            
            for(i in seq(length(oicc1_unique))){
              
              if(i != 3){
                # m1 on x-axis #
                p1[[i]] <- df1_n %>%
                  as.data.frame() %>%
                  mutate(oicc0=factor(oicc0)) %>%
                  dplyr::filter(oicc1 == oicc1_unique[i]) %>%
                  plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                          linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("Treatment outcome ICC: ", oicc1,
                                       "<br>Control outcome ICC: ", oicc0,
                                       "<br>Treatment clusters (n1):", n1,
                                       "<br>Treatment cluster size (m1): ", m1,
                                       "<br>Control clusters (n0):", n0,
                                       "<br>Control cluster size (m0): ", m0,
                                       "<br>HTE power: ", round(power_emp,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                         xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                               standoff=10)),
                         yaxis=list(title="Treatment-arm\n number of clusters (n1)"),
                         legend=list(title=list(text="Control-arm outcome ICC")),
                         margin=0.01)
                
                if(clustering_irgt() == "cluster"){
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i]) %>%
                    plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm\n number of clusters (n0)"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.01)
                }
              }else{
                
                if(clustering_irgt() == "cluster"){
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i]) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm\n number of clusters (n1)"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i]) %>%
                    plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm\n number of clusters (n0)"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                }else{
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i]) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm\n number of clusters (n1)"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                }
              }
              
            }# end row loop
            
            
            
            if(clustering_irgt() == "cluster"){
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs number of clusters',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.6,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.27,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.6,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.27,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }else{
              # m1 on x-axis #
              subplot(p1[[1]],
                      p1[[2]],
                      p1[[3]],
                      nrows=2,# heights = c(0.33,0.33,0.33),
                      margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Treatment-arm cluster size vs number of clusters',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    text = paste(strwrap(paste0("<i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
              
            }# end clustering if/else
          }else if(input$icc_display_irgt == "oICC0_constant"){
            #legend_title <- latex2exp::TeX("$\\rho_x$")
            oicc1_unique <- sort(unique(df1_n[,"oicc1"]))
            oicc0_unique <- sort(unique(df1_n[,"oicc0"]))
            p1 <- p0 <- vector(mode="list", length=length(oicc0_unique))
            
            for(i in seq(length(oicc1_unique))){
              
              if(i != 3){
                # m1 on x-axis #
                p1[[i]] <- df1_n %>%
                  as.data.frame() %>%
                  mutate(oicc1=factor(oicc1)) %>%
                  dplyr::filter(oicc0 == oicc0_unique[i]) %>%
                  plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                          linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                          colors=colors_plot,
                          hoverinfo="text",
                          text=~paste0("Treatment outcome ICC: ", oicc1,
                                       "<br>Control outcome ICC: ", oicc0,
                                       "<br>Treatment clusters (n1):", n1,
                                       "<br>Treatment cluster size (m1): ", m1,
                                       "<br>Control clusters (n0):", n0,
                                       "<br>Control cluster size (m0): ", m0,
                                       "<br>HTE power: ", round(power_emp,4)),
                          height=input$dimension[2]*0.8) %>%
                  layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                         xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                               standoff=10)),
                         yaxis=list(title="Treatment-arm\n number of clusters (n1)"),
                         legend=list(title=list(text="Treatment-arm outcome ICC")),
                         margin=0.01)
                if(clustering_irgt() == "cluster"){
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i]) %>%
                    plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm number\n of clusters (n0)"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.01)
                }
              }else{
                
                if(clustering_irgt() == "cluster"){
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i]) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm number\n of clusters (n0)"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i]) %>%
                    plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm\n number of clusters (n0)"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                }else{
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i]) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm number of clusters (n1)"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                }
              }
              
            }# end row loop
            
            if(clustering_irgt() == "cluster"){
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs number of clusters',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.6,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.27,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.6,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.27,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }else{
              # m1 on x-axis #
              subplot(p1[[1]],#p0[[1]],
                      p1[[2]],#p0[[2]],
                      p1[[3]],#p0[[3]],
                      nrows=2,# heights = c(0.33,0.33,0.33),
                      margin = 0.07, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Treatment-arm cluster size vs number of clusters',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.03,
                    text = paste(strwrap(paste0("<i>Minimum control-arm outcome ICC (", oicc0_unique[1],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.43,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.03,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],")</i>"),
                                         width=0.75*getOption("width")), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=0.25,
                            x=0.5)
                )
            }# end clustering if/else
            
          }# oicc1 constant end here
        }# end sensitivity
      }#end plot_display if/else
      
      
      #### HET CLUSTER ####
    }else if( trial_react() == "het_two"){
      
      # determine effect sizes and variances depending on outcome/covariate type #
      if(outcome_type_react() =="continuous"){
        var_y1 <- (sd_outcome1_react())^2
        var_y0 <- (sd_outcome0_react())^2
        d <- input$mean_diff_HTE
        
        if(covar_type_react()=="continuous"){
          # continuous outcome and covar power
          var_x <- (sd_covar_react())^2
          
        }else if(covar_type_react() == "binary"){
          # continuous outcome binary covar power
          var_x <- (prop_covar_react())*(1-prop_covar_react())
        }
        
      }
      
      if(plot_display_react() == "m_v_power"){
        m1_range <- seq(input$m1_slide_het[1],input$m1_slide_het[2])
        m1_fix <- input$m1_fix_het
        
        m0_range <- seq(input$m0_slide_het[1],input$m0_slide_het[2])
        m0_fix <- input$m0_fix_het
        
        if(sensitivity_het_react() == "est_only"){
          # m1 on x-axis #
          df1_power <- expand.grid(n1=input$n1_fix_het, m1=m1_range,
                                   n0=input$n0_fix_het, m0=m0_fix,
                                   oicc1=oicc_trt_est_het(),
                                   oicc0=oicc_ctrl_est_het(),
                                   cicc=cicc_est_het(),
                                   var_y1=var_y1,
                                   var_y0=var_y0,
                                   var_x=var_x,
                                   d=input$mean_diff_HTE, a=input$sig)
          # print(m0_fix)
          
          # m0 on x-axis #
          df0_power <- expand.grid(n1=input$n1_fix_het, m1=m1_fix,
                                   n0=input$n0_fix_het, m0=m0_range,
                                   oicc1=oicc_trt_est_het(),
                                   oicc0=oicc_ctrl_est_het(),
                                   cicc=cicc_est_het(),
                                   var_y1=var_y1,
                                   var_y0=var_y0,
                                   var_x=var_x,
                                   d=input$mean_diff_HTE, a=input$sig)
          
        }else{
          # m1 on x-axis #
          df1_power <- expand.grid(n1=input$n1_fix_het, m1=m1_range,
                                   n0=input$n0_fix_het, m0=m0_fix,
                                   oicc1=c(oicc_trt_min_het(), oicc_trt_est_het(),
                                           oicc_trt_max_het()),
                                   oicc0=c(oicc_ctrl_min_het(), oicc_ctrl_est_het(),
                                           oicc_ctrl_max_het()),
                                   cicc=c(cicc_min_het(), cicc_est_het(),
                                          cicc_max_het()),
                                   var_y1=var_y1,
                                   var_y0=var_y0,
                                   var_x=var_x,
                                   d=input$mean_diff_HTE, a=input$sig)
          
          # m0 on x-axis #
          df0_power <- expand.grid(n1=input$n1_fix_het, m1=m1_fix,
                                   n0=input$n0_fix_het, m0=m0_range,
                                   oicc1=c(oicc_trt_min_het(), oicc_trt_est_het(),
                                           oicc_trt_max_het()),
                                   oicc0=c(oicc_ctrl_min_het(), oicc_ctrl_est_het(),
                                           oicc_ctrl_max_het()),
                                   cicc=c(cicc_min_het(), cicc_est_het(),
                                          cicc_max_het()),
                                   var_y1=var_y1,
                                   var_y0=var_y0,
                                   var_x=var_x,
                                   d=input$mean_diff_HTE, a=input$sig)
        }
        
        power1_hte_col <- rep(NA, nrow(df1_power))
        power0_hte_col <- rep(NA, nrow(df0_power))
        
        for(i in seq(nrow(df1_power))){
          power1_hte_col[i] <- power_irgt(m1=df1_power[i,"m1"],
                                          m0=df1_power[i,"m0"],
                                          n1=df1_power[i,"n1"],
                                          n0=df1_power[i,"n0"],
                                          oicc1=df1_power[i,"oicc1"],
                                          oicc0=df1_power[i,"oicc0"],
                                          cicc=df1_power[i,"cicc"],
                                          var_y1=df1_power[i,"var_y1"],
                                          var_y0=df1_power[i,"var_y0"],
                                          var_x=df1_power[i,"var_x"],
                                          d=df1_power[i,"d"], a=df1_power[i,"a"])
          power0_hte_col[i] <- power_irgt(m1=df0_power[i,"m1"],
                                          m0=df0_power[i,"m0"],
                                          n1=df0_power[i,"n1"],
                                          n0=df0_power[i,"n0"],
                                          oicc1=df0_power[i,"oicc1"],
                                          oicc0=df0_power[i,"oicc0"],
                                          cicc=df0_power[i,"cicc"],
                                          var_y1=df0_power[i,"var_y1"],
                                          var_y0=df0_power[i,"var_y0"],
                                          var_x=df0_power[i,"var_x"],
                                          d=df0_power[i,"d"], a=df0_power[i,"a"])
          
        }
        
        df0_power <- cbind(df0_power, power0_hte_col)
        df1_power <- cbind(df1_power, power1_hte_col)
        
        if(sensitivity_het_react() == "est_only"){
          # m1 on x-axis #
          p1_est <- df1_power %>%
            as.data.frame() %>%
            dplyr::filter(oicc1 == oicc_trt_est_het(),#oicc_unique[2],
                          oicc0 == oicc_ctrl_est_het(),#cicc_unique[2]
                          cicc == cicc_est_het()#cicc_unique[2]
            ) %>%
            mutate(oicc0=factor(oicc0)) %>%
            plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                    linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("Treatment outcome ICC: ", oicc1,
                                 "<br>Control outcome ICC: ", oicc0,
                                 "<br>Covariate ICC: ", cicc,
                                 "<br>Treatment clusters (n1):", n1,
                                 "<br>Treatment cluster size (m1): ", m1,
                                 "<br>Control clusters (n0):", n0,
                                 "<br>Control cluster size (m0): ", m0,
                                 "<br>HTE power: ", round(power1_hte_col,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("o-ICC = ", oicc_trt_est_het()),
                   xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                         standoff=10)),
                   yaxis=list(title="HTE Power"),
                   legend=list(title=list(text="Control-arm outcome ICC")),
                   margin=0.001)
          
          # m0 on x-axis #
          p0_est <- df0_power %>%
            as.data.frame() %>%
            dplyr::filter(oicc1 == oicc_trt_est_het(),#oicc_unique[2],
                          oicc0 == oicc_ctrl_est_het(),#cicc_unique[2]
                          cicc == cicc_est_het()#cicc_unique[2]
            ) %>%
            mutate(oicc0=factor(oicc0)) %>%
            plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                    linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("Treatment outcome ICC: ", oicc1,
                                 "<br>Control outcome ICC: ", oicc0,
                                 "<br>Covariate ICC: ", cicc,
                                 "<br>Treatment clusters (n1):", n1,
                                 "<br>Treatment cluster size (m1): ", m1,
                                 "<br>Control clusters (n0):", n0,
                                 "<br>Control cluster size (m0): ", m0,
                                 "<br>HTE power: ", round(power0_hte_col,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("o-ICC = ", oicc_est()),
                   xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                         standoff=10)),
                   yaxis=list(title="HTE Power"),
                   legend=list(title=list(text="Control-arm outcome ICC")),
                   margin=0.001)
          
          subplot(p1_est, p0_est,#nrows=1, widths = c(0.5),
                  margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
            layout(title = list(
              text='Cluster size vs HTE power',
              font=list(size=17)
            ),
            #margin=list(pad=50),
            annotations = list(
              list(
                x = 0.23,
                y = 1.0,
                text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Assumed treatment-arm outcome ICC (", oicc_trt_est_het(),"), control-arm outcome ICC (", oicc_ctrl_est_het(),") and covariate ICC (", cicc_est_het(), ")</i>"),
                                     width=0.8*getOption("width")/2), collapse="<br>"),
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = FALSE, font=list(size=13)
              ),
              list(
                x = 0.77,
                y = 1.0,
                text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Assumed treatment-arm outcome ICC (", oicc_trt_est_het(),"), control-arm outcome ICC (", oicc_ctrl_est_het(),") and covariate ICC (", cicc_est_het(), ")</i>"),
                                     width=0.8*getOption("width")/2), collapse="<br>"),
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = FALSE, font=list(size=13)
              )
            ),
            legend=list(orientation="h",
                        yanchor="center",
                        y=0.25,
                        x=0.5)
            )
          
        }else if(sensitivity_het_react() == "sensitivity"){
          if(input$icc_constant_within == "oicc1"){
            if(input$icc_constant == "cicc"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
              oicc0_unique <- sort(unique(df1_power[,"oicc0"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc1_unique))
              
              for(i in seq(length(oicc1_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],"), assumed covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],") and covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],"), assumed covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],"), assumed covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],") and covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],"), assumed covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }else if(input$icc_constant == "oicc0"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
              cicc_unique <- sort(unique(df1_power[,"cicc"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc1_unique))
              
              for(i in seq(length(cicc_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  oicc0 == oicc_ctrl_est_het()) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  oicc0 == oicc_ctrl_est_het()) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  oicc0 == oicc_ctrl_est_het()) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  oicc0 == oicc_ctrl_est_het()) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.001)
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }# end across plot constant
            
            
          }else if(input$icc_constant_within == "oicc0"){
            if(input$icc_constant == "cicc"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
              oicc0_unique <- sort(unique(df1_power[,"oicc0"]))
              cicc_unique <- sort(unique(df1_power[,"cicc"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc0_unique))
              
              for(i in seq(length(oicc1_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                  
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }else if(input$icc_constant == "oicc1"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
              oicc0_unique <- sort(unique(df1_power[,"oicc0"]))
              cicc_unique <- sort(unique(df1_power[,"cicc"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc0_unique))
              
              for(i in seq(length(cicc_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.001)
                  
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }# end across-plot constant
            
            
          }else if(input$icc_constant_within == "cicc"){
            if(input$icc_constant == "oicc0"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
              oicc0_unique <- sort(unique(df1_power[,"oicc0"]))
              cicc_unique <- sort(unique(df1_power[,"cicc"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc0_unique))
              
              for(i in seq(length(oicc1_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc_ctrl_est_het(),
                                  cicc == cicc_unique[i]) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc_ctrl_est_het(),
                                  cicc == cicc_unique[i]) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc_ctrl_est_het(),
                                  cicc == cicc_unique[i]) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc_ctrl_est_het(),
                                  cicc == cicc_unique[i]) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                  
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum covariate ICC (", cicc_unique[1],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed covariate ICC (", cicc_unique[2],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum covariate ICC (", cicc_unique[3],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum covariate ICC (", cicc_unique[1],"), assumed control-arm outcomeICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed covariate ICC (", cicc_unique[2],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum covariate ICC (", cicc_unique[3],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }else if(input$icc_constant == "oicc1"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
              oicc0_unique <- sort(unique(df1_power[,"oicc0"]))
              cicc_unique <- sort(unique(df1_power[,"cicc"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc0_unique))
              
              for(i in seq(length(cicc_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(cicc == cicc_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(cicc == cicc_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(cicc == cicc_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(cicc == cicc_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                  
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum covariate ICC (", cicc_unique[1],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed covariate ICC (", cicc_unique[2],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum covariate ICC (", cicc_unique[3],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum covariate ICC (", cicc_unique[1],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed covariate ICC (", cicc_unique[2],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum covariate ICC (", cicc_unique[3],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }# end across-plot constant
            
            
          }# end with-plot constant
        }# end sensitivity
        
      }else if(plot_display_react() == "n_v_power"){
        
        n1_range <- seq(input$n1_slide_het[1],input$n1_slide_het[2])
        n1_fix <- input$n1_fix_het
        
        n0_range <- seq(input$n0_slide_het[1],input$n0_slide_het[2])
        n0_fix <- input$n0_fix_het
        
        if(sensitivity_het_react() == "est_only"){
          # m1 on x-axis #
          df1_power <- expand.grid(n1=n1_range, m1=input$m1_fix_het,
                                   n0=n0_fix, m0=input$m0_fix_het,
                                   oicc1=oicc_trt_est_het(),
                                   oicc0=oicc_ctrl_est_het(),
                                   cicc=cicc_est_het(),
                                   var_y1=var_y1,
                                   var_y0=var_y0,
                                   var_x=var_x,
                                   d=input$mean_diff_HTE, a=input$sig)
          # print(m0_fix)
          
          # m0 on x-axis #
          df0_power <- expand.grid(n1=n1_fix, m1=input$m1_fix_het,
                                   n0=n0_range, m0=input$m0_fix_het,
                                   oicc1=oicc_trt_est_het(),
                                   oicc0=oicc_ctrl_est_het(),
                                   cicc=cicc_est_het(),
                                   var_y1=var_y1,
                                   var_y0=var_y0,
                                   var_x=var_x,
                                   d=input$mean_diff_HTE, a=input$sig)
          
        }else{
          # m1 on x-axis #
          df1_power <- expand.grid(n1=n1_range, m1=input$m1_fix_het,
                                   n0=n0_fix, m0=input$m0_fix_het,
                                   oicc1=c(oicc_trt_min_het(), oicc_trt_est_het(),
                                           oicc_trt_max_het()),
                                   oicc0=c(oicc_ctrl_min_het(), oicc_ctrl_est_het(),
                                           oicc_ctrl_max_het()),
                                   cicc=c(cicc_min_het(), cicc_est_het(),
                                          cicc_max_het()),
                                   var_y1=var_y1,
                                   var_y0=var_y0,
                                   var_x=var_x,
                                   d=input$mean_diff_HTE, a=input$sig)
          
          # m0 on x-axis #
          df0_power <- expand.grid(n1=n1_fix, m1=input$m1_fix_het,
                                   n0=n0_range, m0=input$m0_fix_het,
                                   oicc1=c(oicc_trt_min_het(), oicc_trt_est_het(),
                                           oicc_trt_max_het()),
                                   oicc0=c(oicc_ctrl_min_het(), oicc_ctrl_est_het(),
                                           oicc_ctrl_max_het()),
                                   cicc=c(cicc_min_het(), cicc_est_het(),
                                          cicc_max_het()),
                                   var_y1=var_y1,
                                   var_y0=var_y0,
                                   var_x=var_x,
                                   d=input$mean_diff_HTE, a=input$sig)
        }
        
        power1_hte_col <- rep(NA, nrow(df1_power))
        power0_hte_col <- rep(NA, nrow(df0_power))
        
        for(i in seq(nrow(df1_power))){
          power1_hte_col[i] <- power_irgt(m1=df1_power[i,"m1"],
                                          m0=df1_power[i,"m0"],
                                          n1=df1_power[i,"n1"],
                                          n0=df1_power[i,"n0"],
                                          oicc1=df1_power[i,"oicc1"],
                                          oicc0=df1_power[i,"oicc0"],
                                          cicc=df1_power[i,"cicc"],
                                          var_y1=df1_power[i,"var_y1"],
                                          var_y0=df1_power[i,"var_y0"],
                                          var_x=df1_power[i,"var_x"],
                                          d=df1_power[i,"d"], a=df1_power[i,"a"])
          power0_hte_col[i] <- power_irgt(m1=df0_power[i,"m1"],
                                          m0=df0_power[i,"m0"],
                                          n1=df0_power[i,"n1"],
                                          n0=df0_power[i,"n0"],
                                          oicc1=df0_power[i,"oicc1"],
                                          oicc0=df0_power[i,"oicc0"],
                                          cicc=df0_power[i,"cicc"],
                                          var_y1=df0_power[i,"var_y1"],
                                          var_y0=df0_power[i,"var_y0"],
                                          var_x=df0_power[i,"var_x"],
                                          d=df0_power[i,"d"], a=df0_power[i,"a"])
          
        }
        
        df0_power <- cbind(df0_power, power0_hte_col)
        df1_power <- cbind(df1_power, power1_hte_col)
        
        if(sensitivity_het_react() == "est_only"){
          # m1 on x-axis #
          p1_est <- df1_power %>%
            as.data.frame() %>%
            dplyr::filter(oicc1 == oicc_trt_est_het(),#oicc_unique[2],
                          oicc0 == oicc_ctrl_est_het(),#cicc_unique[2]
                          cicc == cicc_est_het()#cicc_unique[2]
            ) %>%
            mutate(oicc0=factor(oicc0)) %>%
            plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                    linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("Treatment outcome ICC: ", oicc1,
                                 "<br>Control outcome ICC: ", oicc0,
                                 "<br>Covariate ICC: ", cicc,
                                 "<br>Treatment clusters (n1):", n1,
                                 "<br>Treatment cluster size (m1): ", m1,
                                 "<br>Control clusters (n0):", n0,
                                 "<br>Control cluster size (m0): ", m0,
                                 "<br>HTE power: ", round(power1_hte_col,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("o-ICC = ", oicc_trt_est_het()),
                   xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                         standoff=10)),
                   yaxis=list(title="HTE Power"),
                   legend=list(title=list(text="Control-arm outcome ICC")),
                   margin=0.001)
          
          # m0 on x-axis #
          p0_est <- df0_power %>%
            as.data.frame() %>%
            dplyr::filter(oicc1 == oicc_trt_est_het(),#oicc_unique[2],
                          oicc0 == oicc_ctrl_est_het(),#cicc_unique[2]
                          cicc == cicc_est_het()#cicc_unique[2]
            ) %>%
            mutate(oicc0=factor(oicc0)) %>%
            plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                    linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("Treatment outcome ICC: ", oicc1,
                                 "<br>Control outcome ICC: ", oicc0,
                                 "<br>Covariate ICC: ", cicc,
                                 "<br>Treatment clusters (n1):", n1,
                                 "<br>Treatment cluster size (m1): ", m1,
                                 "<br>Control clusters (n0):", n0,
                                 "<br>Control cluster size (m0): ", m0,
                                 "<br>HTE power: ", round(power0_hte_col,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("o-ICC = ", oicc_est()),
                   xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                         standoff=10)),
                   yaxis=list(title="HTE Power"),
                   legend=list(title=list(text="Control-arm outcome ICC")),
                   margin=0.001)
          
          subplot(p1_est, p0_est,#nrows=1, widths = c(0.5),
                  margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
            layout(title = list(
              text='Number of clusters vs HTE power',
              font=list(size=17)
            ),
            #margin=list(pad=50),
            annotations = list(
              list(
                x = 0.23,
                y = 1.0,
                text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Assumed treatment-arm outcome ICC (", oicc_trt_est_het(),"), control-arm outcome ICC (", oicc_ctrl_est_het(),") and covariate ICC (", cicc_est_het(), ")</i>"),
                                     width=0.8*getOption("width")/2), collapse="<br>"),
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = FALSE, font=list(size=13)
              ),
              list(
                x = 0.77,
                y = 1.0,
                text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Assumed treatment-arm outcome ICC (", oicc_trt_est_het(),"), control-arm outcome ICC (", oicc_ctrl_est_het(),") and covariate ICC (", cicc_est_het(),")</i>"),
                                     width=0.8*getOption("width")/2), collapse="<br>"),
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = FALSE, font=list(size=13)
              )
            ),
            legend=list(orientation="h",
                        yanchor="center",
                        y=0.25,
                        x=0.5)
            )
          
        }else if(sensitivity_het_react() == "sensitivity"){
          if(input$icc_constant_within == "oicc1"){
            if(input$icc_constant == "cicc"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
              oicc0_unique <- sort(unique(df1_power[,"oicc0"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc1_unique))
              
              for(i in seq(length(oicc1_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],"), assumed covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],") and covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],"), assumed covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],"), assumed covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],") and covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],"), assumed covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }else if(input$icc_constant == "oicc0"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
              cicc_unique <- sort(unique(df1_power[,"cicc"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc1_unique))
              
              for(i in seq(length(cicc_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  oicc0 == oicc_ctrl_est_het()) %>%
                    plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  oicc0 == oicc_ctrl_est_het()) %>%
                    plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  oicc0 == oicc_ctrl_est_het()) %>%
                    plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  oicc0 == oicc_ctrl_est_het()) %>%
                    plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.001)
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }# end across plot constant
            
            
          }else if(input$icc_constant_within == "oicc0"){
            if(input$icc_constant == "cicc"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
              oicc0_unique <- sort(unique(df1_power[,"oicc0"]))
              cicc_unique <- sort(unique(df1_power[,"cicc"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc0_unique))
              
              for(i in seq(length(oicc1_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                  
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }else if(input$icc_constant == "oicc1"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
              oicc0_unique <- sort(unique(df1_power[,"oicc0"]))
              cicc_unique <- sort(unique(df1_power[,"cicc"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc0_unique))
              
              for(i in seq(length(cicc_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.001)
                  
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }# end across-plot constant
            
            
          }else if(input$icc_constant_within == "cicc"){
            if(input$icc_constant == "oicc0"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
              oicc0_unique <- sort(unique(df1_power[,"oicc0"]))
              cicc_unique <- sort(unique(df1_power[,"cicc"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc0_unique))
              
              for(i in seq(length(oicc1_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc_ctrl_est_het(),
                                  cicc == cicc_unique[i]) %>%
                    plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc_ctrl_est_het(),
                                  cicc == cicc_unique[i]) %>%
                    plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc_ctrl_est_het(),
                                  cicc == cicc_unique[i]) %>%
                    plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc_ctrl_est_het(),
                                  cicc == cicc_unique[i]) %>%
                    plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                  
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum covariate ICC (", cicc_unique[1],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed covariate ICC (", cicc_unique[2],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum covariate ICC (", cicc_unique[3],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum covariate ICC (", cicc_unique[1],"), assumed control-arm outcomeICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed covariate ICC (", cicc_unique[2],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum covariate ICC (", cicc_unique[3],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }else if(input$icc_constant == "oicc1"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_power[,"oicc1"]))
              oicc0_unique <- sort(unique(df1_power[,"oicc0"]))
              cicc_unique <- sort(unique(df1_power[,"cicc"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc0_unique))
              
              for(i in seq(length(cicc_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(cicc == cicc_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(cicc == cicc_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(cicc == cicc_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~n1,y=~power1_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power1_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm number of clusters (n1)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_power %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(cicc == cicc_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~n0,y=~power0_hte_col, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power0_hte_col,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm number of clusters (n0)",
                                                 standoff=10)),
                           yaxis=list(title="HTE Power"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                  
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Number of clusters vs HTE power',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum covariate ICC (", cicc_unique[1],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed covariate ICC (", cicc_unique[2],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum covariate ICC (", cicc_unique[3],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum covariate ICC (", cicc_unique[1],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed covariate ICC (", cicc_unique[2],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum covariate ICC (", cicc_unique[3],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }# end across-plot constant
            
            
          }# end with-plot constant
        }# end sensitivity
        
      }else if(plot_display_react() == "fixed_power"){
        
        m1_range <- seq(input$m1_slide_het[1],input$m1_slide_het[2])
        m1_fix <- input$m1_fix_het
        
        m0_range <- seq(input$m0_slide_het[1],input$m0_slide_het[2])
        m0_fix <- input$m0_fix_het
        
        if(sensitivity_het_react() == "est_only"){
          # m1 on x-axis #
          df1_n <- expand.grid(#n1=input$n1_fix_het, 
            m1=m1_range,
            #n0=input$n0_fix_het,
            m0=m0_fix,
            power=input$power_het,
            oicc1=oicc_trt_est_het(),
            oicc0=oicc_ctrl_est_het(),
            cicc=cicc_est_het(),
            var_y1=var_y1,
            var_y0=var_y0,
            var_x=var_x,
            d=input$mean_diff_HTE, a=input$sig)
          # print(m0_fix)
          
          # m0 on x-axis #
          df0_n <- expand.grid(#n1=n1_fix,
            m1=m1_fix,
            #n0=n0_range,
            m0=m0_range,
            power=input$power_het,
            oicc1=oicc_trt_est_het(),
            oicc0=oicc_ctrl_est_het(),
            cicc=cicc_est_het(),
            var_y1=var_y1,
            var_y0=var_y0,
            var_x=var_x,
            d=input$mean_diff_HTE, a=input$sig)
          
        }else{
          # m1 on x-axis #
          df1_n <- expand.grid(#n1=n1_range,
            m1=m1_range,
            #n0=n0_fix,
            m0=m0_fix,
            power=input$power_het,
            oicc1=c(oicc_trt_min_het(), oicc_trt_est_het(),
                    oicc_trt_max_het()),
            oicc0=c(oicc_ctrl_min_het(), oicc_ctrl_est_het(),
                    oicc_ctrl_max_het()),
            cicc=c(cicc_min_het(), cicc_est_het(),
                   cicc_max_het()),
            var_y1=var_y1,
            var_y0=var_y0,
            var_x=var_x,
            d=input$mean_diff_HTE, a=input$sig)
          
          # m0 on x-axis #
          df0_n <- expand.grid(#n1=n1_fix,
            m1=m1_fix,
            #n0=n0_range,
            m0=m0_range,
            power=input$power_het,
            oicc1=c(oicc_trt_min_het(), oicc_trt_est_het(),
                    oicc_trt_max_het()),
            oicc0=c(oicc_ctrl_min_het(), oicc_ctrl_est_het(),
                    oicc_ctrl_max_het()),
            cicc=c(cicc_min_het(), cicc_est_het(),
                   cicc_max_het()),
            var_y1=var_y1,
            var_y0=var_y0,
            var_x=var_x,
            d=input$mean_diff_HTE, a=input$sig)
        }
        
        n1_hte_col <- matrix(NA, nrow=nrow(df1_n), ncol=5)
        n0_hte_col <- matrix(NA, nrow=nrow(df0_n), ncol=5)
        
        for(i in seq(nrow(df1_n))){
          n1_hte_col[i,] <- unlist(n_irgt(m1=df1_n[i,"m1"],
                                          m0=df1_n[i,"m0"],
                                          #n1=df1_n[i,"n1"],
                                          #n0=df1_n[i,"n0"],
                                          oicc1=df1_n[i,"oicc1"],
                                          oicc0=df1_n[i,"oicc0"],
                                          cicc=df1_n[i,"cicc"],
                                          var_y1=df1_n[i,"var_y1"],
                                          var_y0=df1_n[i,"var_y0"],
                                          var_x=df1_n[i,"var_x"],
                                          d=df1_n[i,"d"], a=df1_n[i,"a"],
                                          power=df1_n[i,"power"]))
          n0_hte_col[i,] <- unlist(n_irgt(m1=df0_n[i,"m1"],
                                          m0=df0_n[i,"m0"],
                                          #n1=df0_n[i,"n1"],
                                          #n0=df0_n[i,"n0"],
                                          oicc1=df0_n[i,"oicc1"],
                                          oicc0=df0_n[i,"oicc0"],
                                          cicc=df0_n[i,"cicc"],
                                          var_y1=df0_n[i,"var_y1"],
                                          var_y0=df0_n[i,"var_y0"],
                                          var_x=df0_n[i,"var_x"],
                                          d=df0_n[i,"d"], a=df0_n[i,"a"],
                                          power=df0_n[i,"power"]))
          
        }
        
        df0_n <- cbind(df0_n, n0_hte_col)
        df1_n <- cbind(df1_n, n1_hte_col)
        
        colnames(df0_n)[12:16] <- c("n", "n1", "n0", "pi", "power_emp")
        colnames(df1_n)[12:16] <- c("n", "n1", "n0", "pi", "power_emp")
        
        if(sensitivity_het_react() == "est_only"){
          # m1 on x-axis #
          p1_est <- df1_n %>%
            as.data.frame() %>%
            dplyr::filter(oicc1 == oicc_trt_est_het(),#oicc_unique[2],
                          oicc0 == oicc_ctrl_est_het(),#cicc_unique[2]
                          cicc == cicc_est_het()#cicc_unique[2]
            ) %>%
            mutate(oicc0=factor(oicc0)) %>%
            plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                    linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("Treatment outcome ICC: ", oicc1,
                                 "<br>Control outcome ICC: ", oicc0,
                                 "<br>Covariate ICC: ", cicc,
                                 "<br>Treatment clusters (n1):", n1,
                                 "<br>Treatment cluster size (m1): ", m1,
                                 "<br>Control clusters (n0):", n0,
                                 "<br>Control cluster size (m0): ", m0,
                                 "<br>HTE power: ", round(power_emp,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("o-ICC = ", oicc_trt_est_het()),
                   xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                         standoff=10)),
                   yaxis=list(title="Treatment-arm number of clusters (n1)"),
                   legend=list(title=list(text="Control-arm outcome ICC")),
                   margin=0.001)
          
          # m0 on x-axis #
          p0_est <- df0_n %>%
            as.data.frame() %>%
            dplyr::filter(oicc1 == oicc_trt_est_het(),#oicc_unique[2],
                          oicc0 == oicc_ctrl_est_het(),#cicc_unique[2]
                          cicc == cicc_est_het()#cicc_unique[2]
            ) %>%
            mutate(oicc0=factor(oicc0)) %>%
            plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                    linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                    colors=colors_plot,
                    hoverinfo="text",
                    text=~paste0("Treatment outcome ICC: ", oicc1,
                                 "<br>Control outcome ICC: ", oicc0,
                                 "<br>Covariate ICC: ", cicc,
                                 "<br>Treatment clusters (n1):", n1,
                                 "<br>Treatment cluster size (m1): ", m1,
                                 "<br>Control clusters (n0):", n0,
                                 "<br>Control cluster size (m0): ", m0,
                                 "<br>HTE power: ", round(power_emp,4)),
                    height=input$dimension[2]*0.8) %>%
            layout(title=paste("o-ICC = ", oicc_est()),
                   xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                         standoff=10)),
                   yaxis=list(title="Control-arm number of clusters (n0)"),
                   legend=list(title=list(text="Control-arm outcome ICC")),
                   margin=0.001)
          
          subplot(p1_est, p0_est,#nrows=1, widths = c(0.5),
                  margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
            layout(title = list(
              text='Cluster size vs number of clusters',
              font=list(size=17)
            ),
            #margin=list(pad=50),
            annotations = list(
              list(
                x = 0.23,
                y = 1.0,
                text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Assumed treatment-arm outcome ICC (", oicc_trt_est_het(),"), control-arm outcome ICC (", oicc_ctrl_est_het(),") and covariate ICC (", cicc_est_het(), ")</i>"),
                                     width=0.8*getOption("width")/2), collapse="<br>"),
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = FALSE, font=list(size=13)
              ),
              list(
                x = 0.77,
                y = 1.0,
                text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Assumed treatment-arm outcome ICC (", oicc_trt_est_het(),"), control-arm outcome ICC (", oicc_ctrl_est_het(),") and covariate ICC (", cicc_est_het(),")</i>"),
                                     width=0.8*getOption("width")/2), collapse="<br>"),
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = FALSE, font=list(size=13)
              )
            ),
            legend=list(orientation="h",
                        yanchor="center",
                        y=0.25,
                        x=0.5)
            )
          
        }else if(sensitivity_het_react() == "sensitivity"){
          if(input$icc_constant_within == "oicc1"){
            if(input$icc_constant == "cicc"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_n[,"oicc1"]))
              oicc0_unique <- sort(unique(df1_n[,"oicc0"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc1_unique))
              
              for(i in seq(length(oicc1_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size(m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm <br>number of <br>clusters (n1)"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm <br>number of <br>clusters (n0)"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm <br>number of <br>clusters (n1)"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm <br>number of <br>clusters (n0)"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs number of clusters',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],"), assumed covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],") and covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],"), assumed covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],"), assumed covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],") and covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],"), assumed covariate ICC (", cicc_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }else if(input$icc_constant == "oicc0"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_n[,"oicc1"]))
              cicc_unique <- sort(unique(df1_n[,"cicc"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc1_unique))
              
              for(i in seq(length(cicc_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  oicc0 == oicc_ctrl_est_het()) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm <br>number of <br>clusters (n1)"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  oicc0 == oicc_ctrl_est_het()) %>%
                    plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm <br>number of <br>clusters (n0)"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  oicc0 == oicc_ctrl_est_het()) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm <br>number of <br>clusters (n1)"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc1 == oicc1_unique[i],
                                  oicc0 == oicc_ctrl_est_het()) %>%
                    plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC1 = ", oicc1_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm <br>number of <br>clusters (n0)"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.001)
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs number of clusters',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum treatment-arm outcome ICC (", oicc1_unique[1],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed treatment-arm outcome ICC (", oicc1_unique[2],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum treatment-arm outcome ICC (", oicc1_unique[3],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(),")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }# end across plot constant
            
            
          }else if(input$icc_constant_within == "oicc0"){
            if(input$icc_constant == "cicc"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_n[,"oicc1"]))
              oicc0_unique <- sort(unique(df1_n[,"oicc0"]))
              cicc_unique <- sort(unique(df1_n[,"cicc"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc0_unique))
              
              for(i in seq(length(oicc1_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm <br>number of <br>clusters (n1)"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm <br>number of <br>clusters (n0)"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm <br>number of <br>clusters (n1)"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  cicc == cicc_est_het()) %>%
                    plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm <br>number of <br>clusters (n0)"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                  
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs number of clusters',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],"), assumed covariate ICC (", cicc_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }else if(input$icc_constant == "oicc1"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_n[,"oicc1"]))
              oicc0_unique <- sort(unique(df1_n[,"oicc0"]))
              cicc_unique <- sort(unique(df1_n[,"cicc"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc0_unique))
              
              for(i in seq(length(cicc_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm <br>number of <br>clusters (n1)"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm <br>number of <br>clusters (n0)"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm <br>number of <br>clusters (n1)"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(cicc=factor(cicc)) %>%
                    dplyr::filter(oicc0 == oicc0_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m0,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~cicc,
                            linetype=~cicc, color=~cicc,legendgroup=~cicc,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm <br>number of <br>clusters (n0)"),
                           legend=list(title=list(text="Covariate ICC")),
                           margin=0.001)
                  
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs number of clusters',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum control-arm outcome ICC (", oicc0_unique[1],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), "</i>)"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed control-arm outcome ICC (", oicc0_unique[2],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum control-arm outcome ICC (", oicc0_unique[3],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }# end across-plot constant
            
            
          }else if(input$icc_constant_within == "cicc"){
            if(input$icc_constant == "oicc0"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_n[,"oicc1"]))
              oicc0_unique <- sort(unique(df1_n[,"oicc0"]))
              cicc_unique <- sort(unique(df1_n[,"cicc"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc0_unique))
              
              for(i in seq(length(oicc1_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc_ctrl_est_het(),
                                  cicc == cicc_unique[i]) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm <br>number of <br>clusters (n1)"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc_ctrl_est_het(),
                                  cicc == cicc_unique[i]) %>%
                    plot_ly(x=~m0,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm <br>number of <br>clusters (n0)"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc_ctrl_est_het(),
                                  cicc == cicc_unique[i]) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm <br>number of <br>clusters (n1)"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(oicc1=factor(oicc1)) %>%
                    dplyr::filter(oicc0 == oicc_ctrl_est_het(),
                                  cicc == cicc_unique[i]) %>%
                    plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~oicc1,
                            linetype=~oicc1, color=~oicc1,legendgroup=~oicc1,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm <br>number of <br>clusters (n0)"),
                           legend=list(title=list(text="Treatment-arm outcome ICC")),
                           margin=0.001)
                  
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs number of clusters',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum covariate ICC (", cicc_unique[1],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed covariate ICC (", cicc_unique[2],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum covariate ICC (", cicc_unique[3],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum covariate ICC (", cicc_unique[1],"), assumed control-arm outcomeICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed covariate ICC (", cicc_unique[2],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum covariate ICC (", cicc_unique[3],"), assumed control-arm outcome ICC (", oicc_ctrl_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }else if(input$icc_constant == "oicc1"){
              #legend_title <- latex2exp::TeX("$\\rho_x$")
              oicc1_unique <- sort(unique(df1_n[,"oicc1"]))
              oicc0_unique <- sort(unique(df1_n[,"oicc0"]))
              cicc_unique <- sort(unique(df1_n[,"cicc"]))
              p1 <- p0 <- vector(mode="list", length=length(oicc0_unique))
              
              for(i in seq(length(cicc_unique))){
                
                if(i != 3){
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(cicc == cicc_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm <br>number of <br>clusters (n1)"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.01)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(cicc == cicc_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0, showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm <br>number of <br>clusters (n0)"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.01)
                }else{
                  
                  # m1 on x-axis #
                  p1[[i]] <- df1_n %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(cicc == cicc_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m1,y=~n1, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,showlegend=F,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Treatment-arm cluster size (m1)",
                                                 standoff=10)),
                           yaxis=list(title="Treatment-arm <br>number of <br>clusters (n1)"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                  
                  # m0 on x-axis #
                  p0[[i]] <- df0_n %>%
                    as.data.frame() %>%
                    mutate(oicc0=factor(oicc0)) %>%
                    dplyr::filter(cicc == cicc_unique[i],
                                  oicc1 == oicc_trt_est_het()) %>%
                    plot_ly(x=~m0,y=~n0, type='scatter', mode='lines', line=list(width=3), name=~oicc0,
                            linetype=~oicc0, color=~oicc0,legendgroup=~oicc0,
                            colors=colors_plot,
                            hoverinfo="text",
                            text=~paste0("Treatment outcome ICC: ", oicc1,
                                         "<br>Control outcome ICC: ", oicc0,
                                         "<br>Covariate ICC: ", cicc,
                                         "<br>Treatment clusters (n1):", n1,
                                         "<br>Treatment cluster size (m1): ", m1,
                                         "<br>Control clusters (n0):", n0,
                                         "<br>Control cluster size (m0): ", m0,
                                         "<br>HTE power: ", round(power_emp,4)),
                            height=input$dimension[2]*0.8) %>%
                    layout(title=paste("o-ICC0 = ", oicc0_unique[i]),
                           xaxis=list(title=list(text="Control-arm cluster size (m0)",
                                                 standoff=10)),
                           yaxis=list(title="Control-arm <br>number of <br>clusters (n0)"),
                           legend=list(title=list(text="Control-arm outcome ICC")),
                           margin=0.001)
                  
                }
                
              }# end row loop
              
              # m1 on x-axis #
              subplot(p1[[1]],p0[[1]],
                      p1[[2]],p0[[2]],
                      p1[[3]],p0[[3]],
                      nrows=3,# heights = c(0.33,0.33,0.33),
                      margin = 0.09, titleX=T, titleY=T) %>% #list(t=50,b=50,pad=50)) %>%
                layout(title = list(
                  text='Cluster size vs number of clusters',
                  font=list(size=17)
                ),
                #margin=list(pad=50),
                annotations = list(
                  list(
                    x = 0.23,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Treatment arm</b> <br><i>Minimum covariate ICC (", cicc_unique[1],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed covariate ICC (", cicc_unique[2],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.23,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum covariate ICC (", cicc_unique[3],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 1.02,
                    text = paste(strwrap(paste0("<b>Control arm</b> <br><i>Minimum covariate ICC (", cicc_unique[1],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE,
                    font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.59,
                    text = paste(strwrap(paste0("<i>Assumed covariate ICC (", cicc_unique[2],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  ),
                  list(
                    x = 0.77,
                    y = 0.26,
                    text = paste(strwrap(paste0("<i>Maximum covariate ICC (", cicc_unique[3],"), assumed treatment-arm outcome ICC (", oicc_trt_est_het(), ")</i>"),
                                         width=0.8*getOption("width")/2), collapse="<br>"),
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "bottom",
                    showarrow = FALSE, font=list(size=13)
                  )
                ),
                legend=list(orientation="h",
                            yanchor="center",
                            y=-0.1,
                            x=0.6)
                )
            }# end across-plot constant
            
            
          }# end with-plot constant
        }# end sensitivity
      }# end of plot type if/else
      
      
    }# end of trial type if/else
    
    
    
    
  })#end of output plot
  
  #### DESIGN MATRIX OUTPUT ####
  output$design_matrix <- renderTable({
    
    if(trial_react() =="crossover_2"){
      desmat_display <- designMatrix(design = trial_react(), periods = 2)
    }else if(trial_react() =="crossover_m"){
      desmat_display <- designMatrix(design = trial_react(), periods = input$J)
    }else if(trial_react() == "SWD"){
      desmat_display <- designMatrix(design = trial_react(), periods = input$J_1+1, steps = input$J_1)
    }else if(trial_react() == "upload"){
      if(is.null(file1())) stop("User needs to upload design matrix before the function can continue")
      desmat_display <- designMatrix(design = trial_react(), file=file1())
    }else if(trial_react() =="parallel_m"){
      desmat_display <- designMatrix(design = trial_react(), periods = input$J)
    }else{
      desmat_display <- designMatrix(design = trial_react())
    }
    
    head(desmat_display, n=nrow(desmat_display))
    
  },digits=0)
  
  
})


