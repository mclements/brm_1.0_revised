
#' Calculate risks from arctanh RD and log OP
#' 
#' @param atanhrd arctanh of risk difference
#' 
#' @param logop log of odds product
#' 
#' @details The \eqn{log OP} is defined as \eqn{log OP = log[(P(y=1|x=0)/P(y=0|x=0))*(P(y=1|x=1)/P(y=0|x=1))]}. 
#' The inverse hyperbolic tangent function \code{arctanh} is defined as \eqn{arctanh(z) = [log(1+z) - log(1-z)] / 2}. 
#' 
#' @return a vector \eqn{(P(y=1|x=0),P(y=1|x=1))}
#' 
#' @examples getProbScalarRD(0,0)
#' 
#' set.seed(0)
#' logrr = rnorm(10,0,1)
#' logop = rnorm(10,0,1)
#' probs = mapply(getProbScalarRD, logrr, logop)
#' rownames(probs) = c("P(y=1|x=0)","P(y=1|x=1)")
#' probs
#' 
#' @export


getProbScalarRD = function(atanhrd, logop) {
    
    if(length(atanhrd) == 2){
        logop   = atanhrd[2]
        atanhrd = atanhrd[1]
    }
    
    if (logop > 350) {
        if (atanhrd < 0) {
            p0 = 1
            p1 = p0 + tanh(atanhrd)
        } else {
            p1 = 1
            p0 = p1 - tanh(atanhrd)
        }
    } else {
        ## not on boundary logop = 0; solving linear equations
        if (same(logop, 0)) {
            p0 = 0.5 * (1 - tanh(atanhrd))
        } else {
            p0 = (-(exp(logop) * (tanh(atanhrd) - 2) - tanh(atanhrd)) - sqrt((exp(logop) * 
                (tanh(atanhrd) - 2) - tanh(atanhrd))^2 + 4 * exp(logop) * 
                (1 - tanh(atanhrd)) * (1 - exp(logop))))/(2 * (exp(logop) - 
                1))
        }
        p1 = p0 + tanh(atanhrd)
    }
    return(c(p0, p1))
} 
