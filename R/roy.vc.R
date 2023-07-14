#' roy.vc
#'
#' teste roy e bose para comparacoes multiplas na primeira variavel canonica
#'
#' @param mod modelo
#' @param escores vetor com os escores
#' @param trat vetor com os tratamentos
#' @param nrep numero de repeticoes
#' @param alpha nivel alfa
#' @export
#'
roy.vc <- function(mod, escores, trat, nrep, alpha = 0.05) {
meanstrt <- as.vector(xtabs(escores ~ trat)/nrep)
treatment <- unique(trat)
trat <- as.factor(trat)
ecores <- escores
r <- nrep
nometrat <- names(mod[["model"]])[2]
# ANOVA dos escores
m1 <- aov(escores ~ trat)
an_es <- anova(m1)

# roy e bose
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
contr <- abs(diffmat(meanstrt))

H <- summary(mod)$SS[[nometrat]]
E <- summary(mod)$SS$Residuals
nE <- df.residual(mod)
q <- summary(mod)$stats[nometrat,'Df']
p <- nrow(H)
s <- min(c(p,q))
nE <- mod$df.residual
alpha <- 0.05
theta <- uniroot(function(x) {
  doubleWishart(x, p = p, n = q, m = nE, type = "multiple") - (1-alpha)},
  c(0,1))$root
sqr <- anova(m1)['Residuals','Sum Sq']
dms <- sqrt(theta/(1-theta)*sqr/r*2)

n <- length(meanstrt)
z <- data.frame(treatment, meanstrt)
letras <- c(letters[1:26])
w <- z[order(z[, 2], decreasing = FALSE), ]
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
trt <- as.character(w$treatment)
meanstrt <- as.numeric(w$meanstrt)
output <- list(anova.escores = an_es,
               group = data.frame(trat = trt, means = meanstrt, groups = M),
               theta = theta,
               dms = dms)
return(output)
}
