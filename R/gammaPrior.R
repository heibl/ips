
realParameter <- function(x){
  
  id <- get("counter", envir = x$environment)
  out <- paste("RealParameter", id$realParameter, sep = ".")
  id$realParameter <- id$realParameter + 1
  assign("counter", id, envir = x$environment)
  return(out)
}

uniformPrior <- function(x, name = "distr", upper = "Infinity"){
  
  id <- get("counter", envir = x$environment)
  uniform <- list(xmlNode("Uniform", 
                          attrs = c(id = paste0("Uniform.", id$uniform),
                                    name = name,
                                    upper = upper)))
  id$uniform <- id$uniform + 1
  assign("counter", id, envir = x$environment)
  return(uniform)
}

gammaPrior <- function(x, alpha, beta){
  
  id <- get("counter", envir = x$environment)
  alpha <- xmlNode("parameter",
                   attrs = c(estimate = "false",
                             id = paste0("RealParameter.", id$realParameter),
                             name = "alpha",
                             spec = "parameter.RealParameter"),
                   value = alpha)
  id$realParameter <- id$realParameter + 1
  
  beta <- xmlNode("parameter",
                  attrs = c(estimate = "false",
                            id = paste0("RealParameter.", id$realParameter),
                            name = "beta",
                            spec = "parameter.RealParameter"),
                  value = beta)
  id$realParameter <- id$realParameter + 1
  
  gamma <- xmlNode("Gamma", attrs = c(id = paste0("Gamma.", id$gamma),
                             name = "distr"),
          .children = list(alpha, beta))
  id$gamma <- id$gamma + 1
  
  assign("counter", id, envir = x$environment)
  return(gamma)
}

LogNormalPrior <- function(x, M, S){
  
  id <- get("counter", envir = x$environment)
  M <- xmlNode("parameter",
                   attrs = c(estimate = "false",
                             id = paste0("RealParameter.", id$realParameter),
                             name = "M",
                             spec = "parameter.RealParameter"),
                   value = M)
  id$realParameter <- id$realParameter + 1
  
  S <- xmlNode("parameter",
                  attrs = c(estimate = "false",
                            id = paste0("RealParameter.", id$realParameter),
                            name = "S",
                            spec = "parameter.RealParameter"),
                  value = S)
  id$realParameter <- id$realParameter + 1
  
  LogNormal <- xmlNode("LogNormal", attrs = c(id = paste0("LogNormalDistributionModel.", id$logNormal),
                                      name = "distr"),
                   .children = list(M, S))
  id$logNormal <- id$logNormal + 1
  
  assign("counter", id, envir = x$environment)
  return(LogNormal)
}




