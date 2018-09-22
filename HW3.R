# random sample
set.seed(20180909)
rdm_sample <- rcauchy(10, 5, 1)
rdm_sample

log_llh <- function(theta, sample_X){
  lsum <- 0
  for (i in 1:length(sample_X)){
    lsum <- lsum + -log(pi) - log(1 + (theta - sample_X[i])^2)
  }
  lsum
}

#plot of loglikelihood value
x <- seq(0, 10, 0.01)
y <- sapply(x, log_llh, sample_X = rdm_sample)
data_l <- as.data.frame(cbind(y, x))
library(ggplot2)
ggplot(data_l, aes(x, y)) + geom_smooth() +
  labs(title = expression(paste("Loglikelihood of sample against ", theta)),
       x = expression(theta), y = "loglikelihood") +
  theme(plot.title = element_text(hjust = 0.5))

#first derivative of loglikelihood

first_derv <- function(theta, sample_X) {
  first_derv <- 0
  for (i in 1:length(sample_X)){
    first_derv <- first_derv -2 *
      ((theta - sample_X[i])/(1 + (theta - sample_X[i])^2))
  }
  first_derv
}

#second derivative of loglikelihood

second_derv <- function(theta, sample_X) {
  second_derv <- 0
  for (i in 1:length(sample_X)){
    second_derv <- second_derv -2 * 
      ((1-(theta - sample_X[i])^2)/(1 + (theta - sample_X[i])^2)^2)
  }
  second_derv
}

#Newton–Raphson method
newton <- function(init, pre=1e-50, maxrun=200) {
  n <- 1
  xt <- init
  while (n<maxrun){
    fx <- first_derv(xt, rdm_sample)
    fx_d <- second_derv(xt, rdm_sample)
    if (fx == 0) {break}
    ht <- -fx/fx_d
    xt1 <- xt + ht
    if (abs(xt1-xt) < pre) {break}
    xt <- xt1
    n <- n+1
  }
return(c(root = xt, iter = n))
}

init <- seq(-10, 30, 0.5)
result <- as.data.frame(matrix(0, nrow = length(init), ncol = 3))
for (i in 1:length(init) ) {
  result[i,1] <- paste("Initial = ", init[i])
  result[i,2:3] <- newton(init[i])
}
colnames(result) <- (c("Initial", "Root", "# of iterations"))
library(pander)
pander(result)

# plot
result <- cbind(init, result)
ggplot(result, aes(init, Root)) +
  geom_line() + geom_point() + 
  labs(title = expression(paste("Root vs. ", theta)),
       x = expression(paste("Initial value of ", theta)), y = "Root") +
  theme(plot.title = element_text(hjust = 0.5))

# Improve it by halving the steps
# Improved Newton–Raphson method
newton2 <- function(init, pre=1e-50, maxrun=200) {
  n <- 1
  xt <- init
  while (n<maxrun){
    fx <- first_derv(xt, rdm_sample)
    fx_d <- second_derv(xt, rdm_sample)
    if (fx == 0) {break}
    ht <- -fx/fx_d
    xt1 <- xt + ht/2
    if (abs(xt1-xt) < pre) {break}
    xt <- xt1
    n <- n+1
  }
  return(c(root = xt, iter = n))
}

init <- seq(-10, 30, 0.5)
result2 <- as.data.frame(matrix(0, nrow = length(init), ncol = 3))
for (i in 1:length(init) ) {
  result2[i,1] <- paste("Initial = ", init[i])
  result2[i,2:3] <- newton2(init[i])
}
colnames(result2) <- (c("Initial", "Root", "# of iterations"))
library(pander)
pander(result2)

# plot
result2 <- cbind(init, result2)
ggplot(result2, aes(init, Root)) +
  geom_line() + geom_point() + 
  labs(title = expression(paste("Root vs. ", theta, " (Improved)")),
       x = expression(paste("Initial value of ", theta)), y = "Root") +
  theme(plot.title = element_text(hjust = 0.5))

#fixed-point iteration
fix_pnt <- function(init, alpha, pre=1e-50, maxrun=200) {
  n <- 1
  x <- init
  while (n<maxrun){
    fx <- first_derv(x, rdm_sample)
    if (fx == 0) {break}
    Gx <- x + alpha*fx
    if (abs(Gx-x) < pre) {break}
    x <- Gx
    n <- n+1
  }
  return(c(root = x, iter = n))
}

init <- seq(-10, 30, 0.5)
alpha <- c(1, 0.64, 0.25)
result3 <- as.data.frame(matrix(0, nrow = length(init), ncol = 7))
for (i in 1:length(init) ) {
  result3[i,1] <- paste("Init.=", init[i])
  result3[i,2:3] <- fix_pnt(init[i], alpha[1])
  result3[i,4:5] <- fix_pnt(init[i], alpha[2])
  result3[i,6:7] <- fix_pnt(init[i], alpha[3])
}
colnames(result3) <- c("Initial", paste("Root (alpha=",alpha[1],")"),
                       paste0("# of iterations (alpha=",alpha[1],")"),
                       paste0("Root (alpha=",alpha[2],")"),
                       paste0("# of iterations (alpha=",alpha[2],")"),
                       paste0("Root (alpha=",alpha[3],")"),
                       paste0("# of iterations (alpha=",alpha[3],")") )
library(pander)
pander(result3, style="rmarkdown", split.table=Inf, split.cells=Inf)

# plot
result3_plot <- cbind(init, result3)
colnames(result3_plot)[c(3,5,7)] <- c("y1","y2","y3")
ggplot(result3_plot, aes(init)) +
  geom_point(aes(y = y1, colour = "var0")) + 
  geom_point(aes(y = y2, colour = "var1")) +
  geom_point(aes(y = y3, colour = "var2")) +
  labs(title = expression(paste("Root vs. ", theta)),
       x = expression(paste("Initial value of ", theta)), y = "Root") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_colour_discrete(breaks=c("var0", "var1", "var2"),
                      labels=c("alpha = 1", "alpha = 0.64", "alpha = 0.25"))

# Fisher Scoring and Newton-Raphson method
# Fisher Scoring
fisher <- function(init, pre=1e-10, maxrun=200) {
  n <- 1
  Ix <- 10/2
  xt <- init
  while (n<maxrun){
    fx <- first_derv(xt, rdm_sample)
    if (fx == 0) {break}
    xt1 <- xt + fx/Ix
    if (abs(xt1-xt) < pre) {break}
    xt <- xt1
    n <- n+1
  }
  return(c(root = xt, iter = n))
}

init <- seq(-10, 30, 0.5)
result4 <- as.data.frame(matrix(0, nrow = length(init), ncol = 5))
options(digits = 8)
for (i in 1:length(init) ) {
  result4[i,1] <- paste("Initial = ", init[i])
  result4[i,2:3] <- fisher(init[i])
  result4[i,4:5] <- newton(result4[i,2])
}
colnames(result4) <- c("Initial", "Root(Fisher Scoring)",
                       "# of iterations(Fisher Scoring)", 
                       "Root(Newton-Raphson)",
                       "# of iterations (Newton-Raphson)")
library(pander)
pander(result4, style="rmarkdown",split.table=Inf, split.cells=Inf)

# plot
result4_plot <- cbind(init, result4)
colnames(result4_plot)[c(3,5)] <- c("y1","y2")
ggplot(result4_plot, aes(init)) +
  geom_point(aes(y = y1, colour = "var0")) + 
  geom_point(aes(y = y2, colour = "var1")) +
  labs(title = expression(paste("Root vs. ", theta)),
       x = expression(paste("Initial value of ", theta)), y = "Root") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_colour_discrete(breaks=c("var0", "var1"),
                        labels=c("Fisher Scoring", "Fisher & Newton"))









