svyquantile<-function(x,design,quantiles,...) UseMethod("svyquantile", design)

qq<-function (x, w, probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE, 
    type = 7, ...) 
{
    ## if (is.factor(x)) {
    ##     if (is.ordered(x)) {
    ##         if (!any(type == c(1L, 3L))) 
    ##             stop("'type' must be 1 or 3 for ordered factors")
    ##     }
    ##     else stop("factors are not allowed")
    ##     lx <- levels(x)
    ## }
    ## else lx <- NULL
    if (na.rm) 
        x <- x[!is.na(x)]
    else if (anyNA(x)) 
        stop("missing values and NaN's not allowed if 'na.rm' is FALSE")
    eps <- 100 * .Machine$double.eps
    if (any((p.ok <- !is.na(probs)) & (probs < -eps | probs > 
        1 + eps))) 
        stop("'probs' outside [0,1]")
    n <- length(x)
    if (na.p <- any(!p.ok)) {
        o.pr <- probs
        probs <- probs[p.ok]
        probs <- pmax(0, pmin(1, probs))
    }
    np <- length(probs)
    if (n > 0 && np > 0) {
        if (type == 7) {
            index <- 1 + (n - 1) * probs
            lo <- floor(index)
            hi <- ceiling(index)
            x <- sort(x, partial = unique(c(lo, hi)))
            qs <- x[lo]
            i <- which(index > lo & x[hi] != qs)
            h <- (index - lo)[i]
            qs[i] <- (1 - h) * qs[i] + h * x[hi[i]]
        }
        else {
            if (type <= 3) {
                nppm <- if (type == 3) 
                  n * probs - 0.5
                else n * probs
                j <- floor(nppm)
                h <- switch(type, (nppm > j), ((nppm > j) + 1)/2, 
                  (nppm != j) | ((j%%2L) == 1L))
            }
            else {
                switch(type - 3, {
                  a <- 0
                  b <- 1
                }, a <- b <- 0.5, a <- b <- 0, a <- b <- 1, a <- b <- 1/3, 
                  a <- b <- 3/8)
                fuzz <- 4 * .Machine$double.eps
                nppm <- a + probs * (n + 1 - a - b)
                j <- floor(nppm + fuzz)
                h <- nppm - j
                if (any(sml <- abs(h) < fuzz)) 
                  h[sml] <- 0
            }
            x <- sort(x, partial = unique(c(1, j[j > 0L & j <= 
                n], (j + 1)[j > 0L & j < n], n)))
            x <- c(x[1L], x[1L], x, x[n], x[n])
            qs <- x[j + 2L]
            qs[h == 1] <- x[j + 3L][h == 1]
            other <- (0 < h) & (h < 1) & (x[j + 2L] != x[j + 
                3L])
            if (any(other)) 
                qs[other] <- ((1 - h) * x[j + 2L] + h * x[j + 
                  3L])[other]
        }
    }
    else {
        qs <- rep(NA_real_, np)
    }
    if (is.character(lx)) 
        qs <- factor(qs, levels = seq_along(lx), labels = lx, 
            ordered = TRUE)
    if (names && np > 0L) {
        names(qs) <- format_perc(probs)
    }
    if (na.p) {
        o.pr[p.ok] <- qs
        names(o.pr) <- rep("", length(o.pr))
        names(o.pr)[p.ok] <- names(qs)
        o.pr
    }
    else qs
}
