#This is a problem from chapter 2 of this book https://www.amazon.com.br/An%C3%A1lise-Sobreviv%C3%AAncia-Aplicada-Ant%C3%B4nio-Colosimo/dp/8521203845
#It is about a cheese company that would like to test if its product would last longer in one of two different packages, A and B.
#The last 4 observations in each group are censored. I will apply a Bayesian parametric model with a data augmentation step,
#to try to estimate the means of the two populations based on their samples, at the end a Bayesian equivalent of a hypothesis
#test is performed to compare the two means.
#Group A
y1 <- c(31, 40, 43, 44, 46, 46, 47, 48, 48, 49, 50, 50, 60, 60, 60, 60, 60, 60, 60, 60)
censored1 <- c(rep(0, 16), rep(1, 4))
#Group B
y2 <- c(48, 48, 49, 49, 49, 49, 50, 50, 50, 50, 53, 53, 54, 54, 54, 55, 55, 55, 55, 55)
censored2 <- c(rep(0, 16), rep(1, 4))

A <- cbind(censored1,y1)
B <- cbind(censored2,y2)
######################################################################################

updateTheta <- function(a,b,Y){
  n <- length(Y)
  y_bar <- mean(Y)
  a.post <- a + n
  b.post <- b + n * y_bar
  return(rgamma(1,a.post,b.post))
}

######################################################################################
#MCMC using data augmentation 
Niter <- 10000
n_miss1 <- sum(censored1)
n_miss2 <- sum(censored2)
Theta.out <- array(NA, dim = c(Niter,2))
Theta.out[1,] <- c(0.3,0.5)

for(i in 1:Niter){
  Theta.out[i,1] <- updateTheta(0.001,0.001,A[,2])
  Theta.out[i,2] <- updateTheta(0.001,0.001,B[,2])
  #Because of the memoryless property of the exponential distribution
  A[which(A[,1]==1),2] <- 60 + rexp(n_miss1,Theta.out[i,1])
  B[which(B[,1]==1),2] <- 55 + rexp(n_miss2,Theta.out[i,2])
  print(i)
}

Nburn <- 1000

plot(1/Theta.out[seq(Nburn+1,Niter),1],type='l')
plot(1/Theta.out[seq(Nburn+1,Niter),2],type='l')
hist(1/Theta.out[seq(Nburn+1,Niter),1],main="Mean time of group A",probability=T)
hist(1/Theta.out[seq(Nburn+1,Niter),2],main="Mean time of group B",probability=T)
mean(1/Theta.out[seq(Nburn+1,Niter),1])
mean(1/Theta.out[seq(Nburn+1,Niter),2])

#Is Theta1 different from Theta2 ?
Mean_times <- 1/Theta.out[seq(Nburn+1,Niter),]
hist(Mean_times[,1] - Mean_times[,2])
#Using a bayesian hypothesis test as described in 
#https://stats.stackexchange.com/questions/130389/bayesian-equivalent-of-two-sample-t-test
#we can conclude that there is no substantial difference in the mean of the two populations,
#once this is not much better than throwing a fair coin.
mean(Mean_times[,1] - Mean_times[,2] < 0) 

