
#Exercise J-4.1
x = c(1997,907,904,32)

EM <- function(theta, x, maxit, tolerr){
  theta_star = (-1657 + sqrt(3728689))/7680
  for (it in 1:maxit) {
    E = x[1]*(theta/(2+theta))
    theta_tilde <- (E + x[4])/(E + x[2] + x[3] + x[4])
    mod_rel_err <- max(abs((theta_tilde - theta) / max(1,theta_tilde)))
    convergence_ratio = abs(theta_tilde-theta_star)/abs(theta-theta_star)
    print(sprintf('it = %3.0f   theta = %12.12f     MRE=%2.1e     Convergence Ratio = %3.3e',
                  it, theta_tilde , mod_rel_err, convergence_ratio), quote = FALSE)
    if(mod_rel_err < tolerr) {
      break
    }
    theta = theta_tilde
  }
}

a <- EM(0.02, x, 200, 1e-6)



# Exercise J-4.2


# Part a
y <- read.table('C:/Users/mikej/Desktop/CSUF/Math 534/Homework_5_Part_I/ExJ42.txt')
y <- y[2:nrow(y),]
y <- as.numeric(y)
y <- as.matrix(y)

hist(y,breaks = 20,freq = FALSE)
lines(density(y),col='red')


#Part c
alpha = 0.3; beta=0.5
mu1=1; mu2=2; mu3=3
sigma2 = 1
theta0 <- c(0.3,0.5,1,2,3,1)

EM_mixture <- function(y, theta, maxit, tolerr){
  n <- length(y)
  for (it in 1:maxit) {
    #E-Step
    alpha <- theta[1]
    beta <- theta[2]
    pi <- c(alpha, beta, 1-alpha-beta)
    mu1 <- theta[3]
    mu2 <- theta[4]
    mu3 <- theta[5]
    sigma2 <- theta[6]
    f1 <- dnorm(y, mean = mu1, sd=sqrt(sigma2))
    f2 <- dnorm(y, mean = mu2, sd=sqrt(sigma2))
    f3 <- dnorm(y, mean = mu3, sd=sqrt(sigma2))
    N1 <- f1*pi[1]
    N2 <- f2*pi[2]
    N3 <- f3*pi[3]
    D <- N1 + N2 + N3
    Post1 <- N1/D #Posterior probability of belonging to group 1
    Post2 <- N2/D #Posterior probability of belonging to group 2
    Post3 <- N3/D #Posterior probability of belonging to group 3
    Z <- cbind(Post1,Post2,Post3)
    
    # M-Step
    alpha_tilde <- ((1 - beta)*sum(Z[,1]))/(sum(Z[,1]) + sum(Z[,3]))
    beta_tilde <-  ((1 - alpha)*sum(Z[,2]))/(sum(Z[,2]) + sum(Z[,3]))
    mu1_tilde <- sum(Z[,1]*y)/sum(Z[,1])
    mu2_tilde <- sum(Z[,2]*y)/sum(Z[,2])
    mu3_tilde <- sum(Z[,3]*y)/sum(Z[,3])
    sigma2_tilde <- (sum(Z[,1]*((y - matrix(rep(mu1, times=n), nrow=n))^2)) + 
                       sum(Z[,2]*((y - matrix(rep(mu2, times=n), nrow=n)))^2) + 
                       sum(Z[,3]*((y - matrix(rep(mu3, times=n), nrow=n))^2)))/n
    
    #Compute the log-likelihood
    log_likelihood <- sum(Z[,1]*(log(dnorm(y,mu1_tilde,sqrt(sigma2_tilde))) +log(alpha_tilde))) +
      sum(Z[,2]*(log(dnorm(y,mu2_tilde,sqrt(sigma2_tilde))) +log(beta_tilde))) +
      sum(Z[,3]*(log(dnorm(y,mu3_tilde,sqrt(sigma2_tilde))) +log(1-alpha_tilde-beta_tilde)))
    
    #Convergence criteria
    theta_tilde <- c(alpha_tilde, beta_tilde, mu1_tilde, mu2_tilde, mu3_tilde, sigma2_tilde)
    mod_rel_err <- max(abs((theta_tilde - theta) / max(1,theta_tilde)))
    print(sprintf('it = %3.0f   log_likelihood = %12.12f     MRE=%2.1e',
                  it, log_likelihood , mod_rel_err), quote = FALSE)
    theta <- theta_tilde
    if(mod_rel_err < tolerr) {
      break
    }
  }
  return(list(alpha = theta[1], beta = theta[2], mu1 = theta[3], 
       mu2 = theta[4], mu3 = theta[5], sigma2 = theta[6], Z = Z))
  print(sprintf('alpha = %1.4f   beta = %1.4f     mu1=%1.4f     mu2=%1.4f     mu3=%1.4f     sigma2=%1.4f',
                alpha, beta, mu1, mu2, mu3, sigma2), quote = FALSE)
}



theta <- EM_mixture(y, theta0, maxit=200, tolerr = 1e-6)

# Part d
hist(y,breaks = 20,freq = FALSE)
lines(density(y),col='red')
x1 <- seq(min(y), max(y), length.out = 100)
y1 <- theta$alpha*dnorm(x1, theta$mu1, sqrt(theta$sigma2)) +
  theta$beta*dnorm(x1, theta$mu2, sqrt(theta$sigma2)) +
  (1-theta$alpha-theta$beta)*dnorm(x1, theta$mu3, sqrt(theta$sigma2))
lines(x1,y1, col='blue')

# Part e
n=length(y)
class = numeric(n)
for (i in 1:n) {
  class[i] = which(theta$Z[i,] == max(theta$Z[i,]))
}
cases <- data.frame(case = 1:n, class = class)
library(ggplot2)
ggplot(data = cases) +
  geom_point(mapping = aes(x = case, y = class, color = factor(class))) +
  scale_color_manual("class", values = c("red", "green", "blue"))

