#
#     R-code to generate data for Assignment 1.  This code MUST be
#     placed at the start of your own R-script.  You must edit
#     the argument to the set.seed( ) function to fit your own
#     registration number
#
set.seed(3031)
#
#     Generate sample A times with censorng
#
LenA <- 100
TimeA <- sort(round(rexp(LenA, 0.5), digits = 10))
CensorA <- rbinom(LenA,1,0.5)
cbind(TimeA, CensorA)
#
#     Generate sample B times with censorng
#
LenB <- 100
TimeB <- sort(round(rexp(LenB, 0.5), digits = 10))
CensorB <- rbinom(LenB,1,0.6)
cbind(TimeB, CensorB)
#
#     Combine data
#
Time <- c(TimeA, TimeB)
Censor <- c(CensorA, CensorB)
Group <- factor(c(rep(1,LenA), rep(2,LenB)))
cbind(Time, Censor, Group)

####################################################
# Please insert your R code after this line
####################################################
library(survival)

par(mfrow=c(1,2))   #fit two graph side by side
KM.survA <- survfit(Surv(TimeA, CensorA)~1,conf.type="plain") #fitted Cox model for drug A
names(KM.survA)
summary(KM.survA)   #Kaplan-Meier estimate of the survival function for drug A
#Graph of Kaplan-Meier etimate of survivor function for drug A
plot(KM.survA,main="Time to discontinuation of drug A",xlab="Time (in weeks)",ylab="Survival function with 95% CI")

KM.survB <- survfit(Surv(TimeB, CensorB)~1,conf.type="plain") #fitted Cox model for drug B
names(KM.survB)
summary(KM.survB)   #Kaplan-Meier estimate of the survival function for drug B
#Graph of Kaplan-Meier etimate of survivor function for drug B
plot(KM.survB,main="Time to discontinuation of drug B",xlab="Time (in weeks)",ylab="Survival function with 95% CI")

par(mfrow=c(1,1))   #fit original size
#fit both Kaplan-Meier plots into a graph
plot(KM.survA,main="Time to discontinuation of drugs",xlim=c(0,10),xlab="Time (in weeks)",ylab="Survival function with 95% CI",col="red")
lines(KM.survB,col="blue")
legend(7, 0.9, legend=c("Drug A", "Drug B"), col=c("red", "blue"), lty=1:1, cex=0.7)

Cox.fit <- coxph(Surv(Time,Censor) ~ Group)   #Cox fit proportional hazards regression model
summary(Cox.fit)

A <- summary(KM.survA)  #summary of survivor function of drug A
B <- summary(KM.survB)  #summary of survivor function of drug B

#partial log-likelihood function
partloglik = function(beta){
likeA <- numeric()    #empty vector for individual likelihood of patient taking drug A
likeB <- numeric()    #empty vector for individual likelihood of patient taking drug B
multA = 1   #initial value of multA
multB = 1   #initial value of multB

#partial likelihood for survivor function of A
for(i in 1:length(A$time)){
  likeA[i]=1/(length(which(TimeA>=A$time[i]))+(length(which(TimeB>=A$time[i])))*exp(beta))
  multA=multA*likeA[i]
}

#partial likelihood for survivor function of B
for(j in 1:length(B$time)){
  likeB[j]=exp(beta)/(length(which(TimeA>=B$time[j]))+(length(which(TimeB>=B$time[j])))*exp(beta))
  multB=multB*likeB[j]
}

#partial log-likelihood
logfunction = log(multA*multB)
return(logfunction)
}

partial = rep(0,1000)   #empty vector for beta values
y <- seq(0.28, 0.31, length = 1000)   #sequence
  
#obtain value from function
for(g in 1:length(y)){
  partial[g]=partloglik(y[g])
}

#plot graph of partial log-likelihood function
plot(y, partial, main = "Partial log-likelihood", type = "l", xlab = "beta", ylab = "logL", col = "red")

which.max(partial)    #maximum point of partial log-likelihood graph
Beta = y[which.max(partial)] ; Beta #maximum likelihood estimator of Beta

#score function
score <- function(x,h){
  h=10^(-4)   #small h
  (partloglik(x+h/2) - partloglik(x-h/2))/h     #first derivative principle
}

#function to plot score function graph
plotscore = function(Left, Right){
score1=rep(0,1000)  #empty vector for score values
x <- seq(Left, Right, length = 1000)  #sequence

#obtain values from function
for(k in 1:length(x)){
  score1[k]=score(x[k])
}

#plot graph of score function
plot(x, score1, main = "Score function", type = "l", xlab = "beta", ylab = "Score", col = "red", lwd = 2)
grid()
abline(h = 0, lty = 2)  #horizontal line at y=0
}
par(mfrow=c(1,2))
plotscore(-2,2) #plot graph
abline(v=Cox.fit$coefficients)  #vertical line at maximum likelihood estimator of Beta
text(0.29,40,labels=round(Beta,3))  #add label to line
plotscore(0.29,0.3) #plot graph
abline(v=Cox.fit$coefficients)  #vertical line at maximum likelihood estimator of Beta

############WALD TEST########################
#function of second derivative
Deriv2 <- function(x,h){
    (partloglik(x+h) - 2*partloglik(x) + partloglik(x-h))/h^2
}
Info <- -Deriv2(Beta, 10^(-4))  #information function
Var <- 1/Info   #variance
Var
zsq <- ((Beta - 0)^2)/Var ; zsq

############LRT##############################
Lambda <- -2 *(partloglik(0) - partloglik(Beta))

############Score Test######################
U0 <- score(0, 10^(-4))
I0 <- -Deriv2(0, 10^(-4))
Score <- U0^2/I0

#display values
Tests <- c(zsq, Lambda, Score)
names(Tests) <- c("Wald", "LRT", "Score")
Tests
