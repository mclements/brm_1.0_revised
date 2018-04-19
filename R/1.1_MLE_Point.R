
max.likelihood = function(param, y, x, va, vb, alpha.start, beta.start, weights, 
    max.step, thres, pa, pb) {
    
    startpars = c(alpha.start, beta.start)
    
    getProb = if (param == "RR") getProbRR else getProbRD
    
    ## negative log likelihood function
    neg.log.likelihood = function(pars) {
        alpha = pars[1:pa]
        beta = pars[(pa + 1):(pa + pb)]
        p0p1 = getProb(va %*% alpha, vb %*% beta)
        p0 = p0p1[, 1];   p1 = p0p1[, 2]
        
        return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) * weights[x == 0] + 
            (y[x == 0]) * log(p0[x == 0]) * weights[x == 0]) - sum((1 - y[x == 
            1]) * log(1 - p1[x == 1]) * weights[x == 1] + (y[x == 1]) * log(p1[x == 
            1]) * weights[x == 1]))
    }
    
    neg.log.likelihood.alpha = function(alpha){
        p0p1  = getProb(va %*% alpha, vb %*% beta)
        p0    = p0p1[,1];  p1 = p0p1[,2]
        
        return(-sum((1-y[x==0])*log(1-p0[x==0])*weights[x==0] +
                        (y[x==0])*log(p0[x==0])*weights[x==0]) -
                   sum((1-y[x==1])*log(1-p1[x==1])*weights[x==1] +
                           (y[x==1])*log(p1[x==1])*weights[x==1]))  
    }
    
    neg.log.likelihood.beta = function(beta){
        p0p1 = getProb(va %*% alpha, vb %*% beta)
        p0    = p0p1[,1];  p1 = p0p1[,2]
        
        return(-sum((1-y[x==0])*log(1-p0[x==0])*weights[x==0] +
                        (y[x==0])*log(p0[x==0])*weights[x==0]) -
                   sum((1-y[x==1])*log(1-p1[x==1])*weights[x==1] +
                           (y[x==1])*log(p1[x==1])*weights[x==1]))  
    }
    
    
    ## Optimization 

    opt = stats::optim(c(alpha.start,beta.start), neg.log.likelihood, control=list(reltol=thres)) # add hessian=TRUE?
    opt$convergence = (opt$convergence == 0) # change cf. optim()
    
    return(opt)
}

