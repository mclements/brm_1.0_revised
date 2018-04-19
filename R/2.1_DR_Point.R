

dr.estimate.onestep = function(param, y, x, va, vb, alpha.start, beta, pscore, 
    wt, weights, max.step, thres, message) {
    
    startpars = c(alpha.start)  # pars only contain alpha
    ## DR estimation equation^2
    if (param == "RR") {
        dr.objective = function(pars) {
            p0 = t(mapply(getProbScalarRR, va %*% pars, vb %*% beta))[, 1]
            H.alpha = as.vector(y * exp(-x * (va %*% pars)))
            tmp = t(va) %*% (wt * (x - pscore) * (H.alpha - p0) * weights)
            return(sum(tmp^2))
        }
    }
    if (param == "RD") {
        dr.objective = function(pars) {
            p0 = t(mapply(getProbScalarRD, va %*% pars, vb %*% beta))[, 1]
            H.alpha = y - x * tanh(va %*% pars)
            tmp = t((H.alpha - p0) * (x - pscore)) %*% (va * wt * weights)
            return(sum(tmp^2))
        }
    }
    
    Diff = function(x,y) abs(x-y)/(x+thres)
    step = 0;    value.old = diff = thres + 1
    while(diff > thres & value.old > thres & step < max.step){
        step = step + 1
        opt = stats::optim(startpars,dr.objective,control=list(maxit=max.step))
        diff = Diff(opt$value,value.old)
        value.old = opt$value
        startpars = opt$par
        if(message & step %% 10 == 0){
            cat("This is the ", step, "th step. The optimum value is ",
                opt$value," after ",opt$counts[1]," iterations \n",sep="")   
        }
    }
    
    opt = list(par = opt$par, convergence = (step < max.step), 
               value = opt$value)
    
    return(opt)
}



dr.estimate.noiterate = function(param, y, x, va, vb, vc, alpha.ml, beta.ml, 
    gamma, optimal, weights, max.step, thres, alpha.start, message) {
    
    pscore = as.vector(expit(vc %*% gamma))
    
    if (optimal == TRUE) {
        if (param == "RR") {
            p0 = t(mapply(getProbScalarRR, va %*% alpha.ml, vb %*% beta.ml))[, 
                1]
            wt = as.vector(1/(1 - p0 + (1 - pscore) * (exp(-va %*% alpha.ml) - 
                1)))
        }
        if (param == "RD") {
            p0 = t(mapply(getProbScalarRD, va %*% alpha.ml, vb %*% beta.ml))[, 
                1]
            rho = as.vector(tanh(va %*% alpha.ml))
            wt = (1 - rho) * (1 + rho)/(p0 * (1 - p0) + rho * (1 - pscore) * 
                (1 - 2 * p0 - rho))
        }
    } else {
        wt = rep(1, length(pscore))
    }
    
    if (is.null(alpha.start)) 
        alpha.start = alpha.ml
    
    alpha.dr.opt = dr.estimate.onestep(param, y, x, va, vb, alpha.start, beta.ml, 
        pscore, wt, weights, max.step, thres, message)
    
    # if(MESSAGE){ print(paste('DR One Step: ',' Alpha:
    # ',paste(round(alpha.dr,5),collapse=', '),' Beta:
    # ',paste(round(beta.ml,5),collapse=', '))) }
    
    return(alpha.dr.opt)
} 
