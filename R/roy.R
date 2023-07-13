#' roy
#'
#' teste roy e bose para comparacoes multiplas
#'
#' @param mod modelo
#' @param var variavel
#' @param trt tratamento
#' @param alpha nivel alfa
#' @export
roy <- function(mod, var, trt, alpha = 0.05) {
  data <- data.frame(trt = mod[["model"]][[trt]], y = mod[["model"]][[1]][,var])
  means <- data %>% group_by(trt) %>% summarise(y = mean(y)) %>% as.data.frame

  H <- summary(mod)$SS[[trt]]
  E <- summary(mod)$SS$Residuals

  treatments <- means$trt
  ntr <- nlevels(treatments)
  meanstrt <- means$y

  q <- summary(mod)$stats['trat','Df']
  p <- nrow(H)
  s <- min(c(p,q))
  nE <- mod$df.residual
  ml <- 0.5*(abs(p-q)-1)
  nl <- 0.5*(nE-p-1)
  r <- unique(xtabs( ~ trt, data = data))
  sqr <- E[var,var]
  theta <- uniroot(function(x) {rootWishart::doubleWishart(x, p = p, n = q, m = nE, type = "multiple") - (1-alpha)}, c(0,1))$root
  dms <- sqrt((theta/(1-theta)*sqr/r*2))

  comb <- utils::combn(ntr, 2)
  nn <- ncol(comb)
  dif <- rep(0, nn)
  sig <- NULL
  pvalue <- rep(0, nn)
  for (k in 1:nn) {
    i <- comb[1, k]
    j <- comb[2, k]
    dif[k] <- unlist(means[i, 2][[1]] - means[j, 2][[1]])
    pvalue[k] <- round(1 - rootWishart::doubleWishart((dif[k]^2/(dif[k]^2+(sqr*2/r))), p = p, n = q, m = nE, type = "multiple"),4)
    sig[k] <- " "
    if (pvalue[k] <= 0.001)
      sig[k] <- "***"
    else if (pvalue[k] <= 0.01)
      sig[k] <- "**"
    else if (pvalue[k] <= 0.05)
      sig[k] <- "*"
    else if (pvalue[k] <= 0.1)
      sig[k] <- "."
  }
  tr.i <- levels(means$trt)[comb[1, ]]
  tr.j <- levels(means$trt)[comb[2, ]]
  comparison <- data.frame(difference = as.vector(dif), pvalue = pvalue, signif. = sig)
  rownames(comparison) <- paste(tr.i, tr.j, sep = " - ")

  diffmat <- function(x){
    D = matrix(as.numeric(NA), NROW(x), NROW(x))
    for (i in 1:NROW(x)){
      d = x[[i]] - x[-i]
      D[i,-i] = d
    }
    if (!all(is.na(diag(D)))){
      stop("Not all diagonal elements zero")
    }
    diag(D) = 0
    if (!is.null(names(x))) colnames(D) = rownames(D) = names(x)
    return(D)
  }
  contr <- abs(diffmat(means$y))
  n <- length(meanstrt)
  z <- data.frame(treatments, meanstrt)
  letras <- c(letters[1:26])
  w <- z[order(z[, 2], decreasing = TRUE), ]
  M <- rep("", n)
  k <- 1
  k1 <- 0
  j <- 1
  i <- 1
  cambio <- n
  cambio1 <- 0
  chequeo = 0
  M[1] <- letras[k]
  q <- as.numeric(rownames(w))
  while (j < n) {
    chequeo <- chequeo + 1
    if (chequeo > n)
      break
    for (i in j:n) {
      s <- contr[q[i], q[j]] < dms
      if (s) {
        if (
          sapply(M[i], function(x)
          { y <- sub(" +$", "", x)
          p1 <- nchar(y)
          cc <- substr(y, p1, p1)
          return(cc)}) != letras[k])
          M[i] <- paste(M[i], letras[k], sep = "")
      }
      else {
        k <- k + 1
        cambio <- i
        cambio1 <- 0
        ja <- j
        for (jj in cambio:n) M[jj] <- paste(M[jj], "",
                                            sep = "")
        M[cambio] <- paste(M[cambio], letras[k], sep = "")
        for (v in ja:cambio) {
          if (contr[q[v], q[cambio]] <= dms) {
            j <- j + 1
            cambio1 <- 1
          }
          else break
        }
        break
      }
    }
    if (cambio1 == 0)
      j <- j + 1
  }
  w <- data.frame(w, stat = M)
  trt <- as.character(w$treatments)
  meanstrt <- as.numeric(w$meanstrt)
  output <- list(contrasts = comparison, group = data.frame(trt, meanstrt, groups = M), dms = dms, theta = theta)

  return(output)
}
