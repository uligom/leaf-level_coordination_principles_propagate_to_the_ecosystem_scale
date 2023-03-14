## Implementation of bootstrapped eigenvector
## Follows Peres-Neto et al (2003)
netoboot <- function (x, permutations = 1000, ...)
{
    ## x is a data.frame
    obs <- x$c1
    pcnull <- ade4::dudi.pca(x, ...)
    res <- pcnull$c1
    out <- matrix(0, nrow = nrow(res), ncol = ncol(res))
    colnames(out) <- names(res)
    rownames(out) <- rownames(res)
    N <- nrow(x)
    out <- vector(mode = "list", length = permutations)
    for (i in 1:permutations) {
        pc <- dudi.pca(x[sample(N, replace = TRUE), ], ...)
        
        r <-  cor(pcnull$c1, pc$c1)
        k <- apply(abs(r), 2, which.max)
        reve <- sign(diag(r[k,]))
        sol <- pc$c1[ ,k]
        sol <- sweep(sol, 2, reve, "*")
        out[[i]] <- sol
        ## out <- out + ifelse(res > 0, sol <=  0, sol >= 0)
    }
    #"out/permutations
    out
}
# debugonce(netoboot)

plotboot <- function(obs, sim, mfrow =  n2mfrow(ncol(obs))){
    ndim <- ncol(obs)
    par(mfrow = mfrow)
    for(i in 1:ndim){
        df.sim <- t(do.call("cbind", lapply(sim, function(x) x[,i]))) 
        colnames(df.sim) <- rownames(obs)
        boxplot(df.sim, main = paste("Axis", i))
        abline(h = 0, col = 'red', lty = 3)
    }
}

computePval_peres <- function(obs, sim, mfrow =  n2mfrow(ncol(obs))){
  ndim <- ncol(obs)
  #par(mfrow = mfrow)
  pvalue <- matrix(NA, nrow = ndim, ncol = nrow(obs))
  for(i in 1:ndim){
    df.sim <- t(do.call("cbind", lapply(sim, function(x) x[,i])))
    colnames(df.sim) <- rownames(obs)
    
    posval <- apply(df.sim, 2, function(i) sum(i >= 0))
    posval[which(obs[,i] >= 0)] <- NA
    negval <- apply(df.sim, 2, function(i) sum(i <= 0))
    posval[which(obs[,i] >= 0)] <- negval[which(obs[,i] >= 0)]
    #-- Calculate pval vector
    
    pvalue[i,] <- (posval + 1)/(dim(df.sim)[1]+1)

  }
  colnames(pvalue) <- rownames(obs)
  return(pvalue)
}

# p95 <- aaa %>% summarise_all(quantile(., 0.95))
# p05 <- aaa %>% summarise_all(quantile(., 0.05))
# mean <- a1 <- aaa %>% summarise_all(mean = mean(., na.rm = FALSE))
# apply(aaa, 2, function(i) sum(i > 0))
# 
# p05<-apply(aaa, 2, function(i) quantile(i, 0.05))
# p95<-apply(aaa, 2, function(i) quantile(i, 0.95))
# meanp<-apply(aaa, 2, function(i) mean(i))
# sdp<-apply(aaa, 2, function(i) sd(i))
# sdp/abs(meanp)

## selection of a proportion of data
subselect.data <- function(df, percent = 0.95){
    df0 <- scale(df, scale = FALSE)
    di <- rowSums(df0^2)
    thresh <- quantile(di, percent)
    idx <- which(di>thresh)
    if(length(idx) > 0)
        df <- df[-idx,]
    return(df)
}


## permutation models
nullmodel <- function(tab, model){
    if(model == 1){
        res <- apply(tab, 2, function(x) runif(nrow(tab), min = min(x), max = max(x)))
    }
    
    if(model == 2){
        res <- scale(matrix(rnorm(ncol(tab) * nrow(tab)), nrow(tab), ncol(tab)))
    }
    
    if(model == 3){
        res <- apply(tab,2,sample)
    }
    
    if(model == 4){
        corM <- cor(tab)
        res <- scale(scale(matrix(rnorm(ncol(tab) * nrow(tab)), nrow(tab), ncol(tab)))%*%chol(corM))
    }
    
    return(res)
}

## concentration analysis
myfunc <- function(pc){
    nbreaks <- 11
    ngroups <- nbreaks - 1
    idx <- 1:ngroups
    thegrid <- expand.grid(LeafArea = idx, Nmass = idx, LMA = idx, Height = idx, SeedMass = idx, StemDensity = idx)
    pc.cut <- sapply(as.data.frame(pc), function(x) cut(x, breaks = seq(-5, 5, length = nbreaks), labels = FALSE))
    whichcell <- function(x) 1 + sum((x - 1) * 10^(0:(length(x)-1)))
    idx <- apply(pc.cut, 1, whichcell)
    n <- rep(0, nrow(thegrid))
    n[as.numeric(names(table(idx)))] <- table(idx)
    res <- cbind(thegrid, n = n)[as.numeric(names(table(idx))),]
    res
}

plot.curves <- function(obs, null, main = ""){
    null.cum <- lapply(null, function(x) cumsum(sort(x$n, dec = T)) / sum(x$n))
    obs.cum <- cumsum(sort(obs$n, dec = T)) / sum(obs$n)
    nmax <- max(length(obs.cum), sapply(null.cum, length))
    res.null <- matrix(1, nmax, length(null.cum))
    for(j in 1:length(null.cum))
        res.null[1:length(null.cum[[j]]), j] <- null.cum[[j]]
    obs.tmp <- rep(1, nmax)
    obs.tmp[1:length(obs.cum)] <- obs.cum
    res <- apply(res.null, 1, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
    res <- rbind(obs.tmp, res)
    plot(1:ncol(res), res[1, ], ylim = range(res), 
         ty = "n", xlab = "number of cells", ylab = "Proportion of species", main = main)
    polygon(c(1:ncol(res), rev(1:ncol(res))), 
            c(res[2, ], rev(res[4, ])), col = "lightgrey", 
            border = "grey")
    lines(1:ncol(res), res[3, ], col = "darkgrey", 
          ty = "l", lty = 2)
    points(1:ncol(res), res[1, ], ty = "l", lwd = 3)
    print(paste("N_10 (obs):", which(res[1,] >0.1)[1]))
    print(paste("N_10 (null):", which(res[1,] >0.5)[1]))
    print(paste("N_50 (obs):", which(res[3,] >0.1)[1]))
    print(paste("N_50 (null):", which(res[3,] >0.5)[1]))
}

