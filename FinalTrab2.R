library(survival)
library(dplyr)
data(colon)
head(colon)
attach(colon)

#O banco de dados contém duas linhas por pessoa, a variável etype refere-se a recor-
#rência ou à morte , vamos primeiro análisar para o caso de morte, isto é, etype==2

death <- filter(colon, etype == 2)
head(death)
detach(colon); attach(death)
agec <- ifelse(age > 65, 1, 0)
nodesc <- ifelse(nodes > 4, 1, 0)
table(rx)
table(sex)
table(obstruct)
table(perfor)
table(adhere)
table(differ)
table(extent)
table(surg)
table(node4)
table(etype)
table(agec)
table(nodesc)

###############Ajustando sobrevivência por Kapla-Meier
ekm <- survfit(Surv(time, status)~1)
summary(ekm)
plot(ekm$time, ekm$surv, type = 'l', main = 'Kaplan-Meier S(t)', xlab = 'tempos',
     ylab = 'S(t)')

#############Vamos separar a variável sex para verificar se há alguma diferença na
#sobrevivência de homens e mulheres
library(survminer)
ggsurvplot(survfit(Surv(time, status)~sex, data = death), pval = T)
#Não parece haver diferença significativa na sobrevivência de homens e mulheres

####Vamos verificar se os dados se ajustam bem a algum dos modelos paramétricos
#usados anterioremente
summary(ekm)
st <- ekm$surv
temp <- ekm$time
invst <- qnorm(st)
par(mfrow = c(1, 3))
plot(temp, -log(st), pch = 16, xlab = "Tempos", ylab = "-log(S(t))",
     main = "Exponencial")
plot(log(temp), log(-log(st)), pch = 16, xlab = "log(tempos)",
     ylab = "log(-log(S(t)))", main = "Weibull")
plot(log(temp), invst, pch = 16, xlab = "log(tempos)",
     ylab = expression(Phi^-1*(S(t))), main = "Log-normal")

#Log normal está mais próximo de uma reta, vamos olhar a fç de sob dos modelos
#versus a fç de KM

###Exponencial
adj11 <- survreg(Surv(time, status) ~ 1, dist = "exponential")
alpha <- exp(adj11$coefficients[1])
ste <- exp(-tm/alpha)

###Weibull
adj22 <- survreg(Surv(time, status) ~ 1, dist = "weibull")
alpha2 <- exp(adj22$coefficients[1])
gama <- 1/adj22$scale
stw <- exp(-(tm/alpha2)^gama)

###Log-Normal
adj33 <- survreg(Surv(time, status) ~ 1, dist = "lognorm")
stln <-  pnorm((-log(tm) + adj33$coefficients[1])/adj33$scale)

par(mfrow = c(1, 3))

plot(st, ste, pch = 16, ylim = c(0, 1), xlab = "S(t): Kaplan-Meier",
     ylab = "S(t): Exponential", main = "Exponencial")
lines(c(0, 1), c(0, 1), type = "l", lty = 1)

plot(st, stw, pch=16, ylim=range(c(0.0,1)), xlim=range(c(0,1)),
     xlab = "S(t): Kaplan-Meier", ylab="S(t): Weibull", main = "Weibull")
lines(c(0,1),c(0,1),type="l",lty=1)

plot(st,stln,pch=16,ylim=range(c(0.0,1)), xlim=range(c(0,1)),
     xlab="S(t): Kaplan-Meier", ylab="S(t): log-normal", main = "Log-Norm")
lines(c(0,1),c(0,1),type="l",lty=1)

#De fato o Log normal aparente ser o melhor dos modelos para esses dados vamos
#obter o teste da razão de verossimilhança para obter uma medida númerica de qual
#seria o melhor modelo compando todos eles com a gama generalizada

require(flexsurv)
adj4 <- flexsurvreg(Surv(time, status) ~ 1, dist = "gengamma")


###Teste Razão de Verossimilhança - forçando a barra é possível não rejeitar o 
#modelo
###exponencial ao nível de 5%

###EXponencial
TRV1 <- 2 * (adj4$loglik - adj11$loglik[2])
pexp <- 1 - pchisq(TRV1, df = 2)

###Weibull
TRV2 <- 2 * (adj4$loglik - adj22$loglik[2])
pwei <- 1 - pchisq(TRV2, df = 1)

###Log-Normal
TRV3 <- 2 * (adj4$loglik - adj33$loglik[2])
plog <- 1 - pchisq(TRV3, df = 1)

logl <- c(adj4$loglik, adj11$loglik[2], adj22$loglik[2], adj33$loglik[2])
trv <- c("", round(TRV1, digits = 4), round(TRV2, digits = 4), 
         round(TRV3, digits = 4))
pv <- c("", round(pexp, digits = 4), round(pwei, digits = 4),
        round(plog, digits = 4))


t <- data.frame(logl, trv, pv, row.names = c("Gama Generalizada", "Exponencial",
                "Weibull", "Log-Normal"))
colnames(t) <- c("Log-Verossimilhança", "TRV", "p-valor")
t

###Todos os modelos foram rejeitados

###Modelo de Cox
#################################################################################
####SELEÇÃO DE VARIÁVEIS MODELO DE COX###
#################################################################################

#################################################################################
####PASSO 1###
#################################################################################

#MODELO NULO
fit0 <- coxph(Surv(time, status) ~ 1, data = death, x = T, 
              method = "breslow")
summary(fit0)

#Categoria Lev+5FU significante ao nível de .001
fit1 <- coxph(Surv(time, status) ~ as.factor(rx), data = death, x = T, 
              method = "breslow")
summary(fit1)

#Não significante
fit2 <- coxph(Surv(time, status) ~ sex, data = death, x = T, 
              method = "breslow")
summary(fit2)

#É significante ao nível de .01
fit3 <- coxph(Surv(time, status) ~ obstruct, data = death, x = T, 
              method = "breslow")
summary(fit3)

#Não significante
fit4 <- coxph(Surv(time, status) ~ perfor, data = death, x = T, 
              method = "breslow")
summary(fit4)

#Significante ao nível de .001
fit5 <- coxph(Surv(time, status) ~ adhere, data = death, x = T, 
              method = "breslow")
summary(fit5)

#Significativo ao nível de .001
fit6 <- coxph(Surv(time, status) ~ nodesc, data = death, x = T, 
              method = "breslow")
summary(fit6)

#Significativo ao nível de .001
fit7 <- coxph(Surv(time, status) ~ as.factor(differ), data = death, x = T, 
              method = "breslow")
summary(fit7)

#Significativo ao nível de .001
fit8 <- coxph(Surv(time, status) ~ as.factor(extent), data = death, x = T, 
              method = "breslow")
summary(fit8)

#Significativo ao nível de .001
fit9 <- coxph(Surv(time, status) ~ surg, data = death, x = T, 
              method = "breslow")
summary(fit9)

logv <- c(fit0$loglik, fit1$loglik[2], fit2$loglik[2], fit3$loglik[2], 
          fit4$loglik[2],
          fit5$loglik[2], fit6$loglik[2], fit7$loglik[2], fit8$loglik[2], 
          fit9$loglik[2])
vloglik <- c("", 12.1, 0.02, 5.05, 0.33, 6.04, 86.3, 15.2, 29.2, 5.01)
p.valor <- c("", 0.00231, 0.887, 0.0246, 0.565, 0.014, 0, 0.000487, 2.03e-06,
             0.0253)

frame <- data.frame(logv, vloglik, p.valor, row.names = c('nulo', 'rx', 'sex',
                                                          'obstruct', 'perfor', 'adhere', 'nodesc', 'differ', 'extent', 
                                                          'surg'))
colnames(frame) <- c('log-veros', 'TRV', 'p-valor')

#################################################################################
####PASSO 2###
#################################################################################

#########################Variáveis .1 de significância##########################
fit <- coxph(Surv(time, status) ~ as.factor(rx) + obstruct + adhere + 
             as.factor(differ) + as.factor(extent) + surg + nodesc, x = T,
             data = death, method = 'breslow')
summary(fit)
#######################Retirando adhere#########################
fit2 <- coxph(Surv(time, status) ~ as.factor(rx) + obstruct + as.factor(differ) + 
                as.factor(extent) + surg + nodesc, x = T, data = death, 
              method = 'breslow')
summary(fit2)
######################Retirando obstruct#######################
fit3 <- coxph(Surv(time, status) ~ as.factor(rx) + adhere + as.factor(differ) + 
                as.factor(extent) + surg + nodesc , x = T, data = death,
              method = 'breslow')
summary(fit3)
#####################Retirando differ#########################
fit4 <- coxph(Surv(time, status) ~ as.factor(rx) + obstruct + adhere + 
              as.factor(extent) + surg + nodesc, x = T, data = death, 
              method = 'breslow')
summary(fit4)
###################Retirando extent###########################
fit5 <- coxph(Surv(time, status) ~ as.factor(rx) + obstruct + adhere + 
              as.factor(differ) + surg + nodesc, x = T, data = death, 
              method = 'breslow')
summary(fit5)
##################Retirando surg#############################
fit6 <- coxph(Surv(time, status) ~ as.factor(rx) + obstruct + adhere + 
              as.factor(differ) + as.factor(extent) + nodesc, x = T, data = death, 
              method = 'breslow')
summary(fit6)
###############Retirando rx#################################
fit9 <- coxph(Surv(time, status) ~ obstruct + adhere + as.factor(differ) + 
              as.factor(extent) + surg + nodesc, x = T, data = death, 
              method = 'breslow')
summary(fit9)

###############Retirando nodesc#################################
fit10 <- coxph(Surv(time, status) ~ as.factor(rx) + obstruct + adhere + 
             as.factor(differ) + as.factor(extent) + surg, x = T, data = death, 
              method = 'breslow')
summary(fit10)

#################################################################################
########Excluir apenas variável adhere############
#################################################################################
mod <- coxph(Surv(time, status) ~ as.factor(rx) + obstruct + as.factor(extent) +
             surg + nodesc, x = T, data = death, method = 'breslow')
summary(mod)

################################################################################
####Verificando adequação do modelo de Cox####
################################################################################
####################RX
adj <-coxph(Surv(time[rx=='Lev'], status[rx=='Lev'])~1,data=death,x=T,
            method="breslow")
ss<- survfit(adj)
s0<-round(ss$surv,digits=5)
H0<- -log(s0)
plot(ss$time,log(H0),xlab="Tempos",
     ylab=expression(log(Lambda[0]*(t))), bty="n",type="s")
adj2 <-coxph(Surv(time[rx== 'Lev+5FU'], status[rx== 'Lev+5FU'])~1,data=death,x=T,
             method="breslow")
ss<- survfit(adj2)
s0<-round(ss$surv,digits=5)
H0<- -log(s0)
lines(ss$time,log(H0),type="s",lty=2)
legend(1500,-2,lty=c(2,1),c("rx = Lev ","rx = Lev+5FU"),lwd=1,bty="n",cex=0.7)
title("Tratamento (rx)")

################OBSTRUCT
adj <-coxph(Surv(time[obstruct==0], status[obstruct==0])~1,data=death,x=T,
            method="breslow")
ss<- survfit(adj)
s0<-round(ss$surv,digits=5)
H0<- -log(s0)
plot(ss$time,log(H0),xlab="Tempos",
     ylab=expression(log(Lambda[0]*(t))), bty="n",type="s")
adj2 <-coxph(Surv(time[obstruct==1], status[obstruct==1])~1,data=death,x=T,
             method="breslow")
ss<- survfit(adj2)
s0<-round(ss$surv,digits=5)
H0<- -log(s0)
lines(ss$time,log(H0),type="s",lty=2)
legend(1500,-3,lty=c(2,1),c("obstruct = 0","obstruct = 1"),lwd=1,bty="n",cex=0.7)
title("Obstrução do Colon")

#############################Para variável extent (Viola Hipóstese)
adj <-coxph(Surv(time[extent==2], status[extent==2])~1,data=death,x=T,
            method="breslow")
ss<- survfit(adj)
s0<-round(ss$surv,digits=5)
H0<- -log(s0)
plot(ss$time,log(H0),xlab="Tempos",
     ylab=expression(log(Lambda[0]*(t))), bty="n",type='s')
adj2 <-coxph(Surv(time[extent==3], status[extent==3])~1,data=death,x=T,
             method="breslow")
ss<- survfit(adj2)
s0<-round(ss$surv,digits=5)
H0<- -log(s0)
lines(ss$time,log(H0),type='S',lty=2)
adj3 <-coxph(Surv(time[extent==4], status[extent==4])~1,data=death,x=T,
             method="breslow")
ss<- survfit(adj3)
s0<-round(ss$surv,digits=5)
H0<- -log(s0)
lines(ss$time,log(H0),type='S',lty=3)
legend(1500,-3,lty=c(1,2,3),c("extent = 2","extent = 3", "extent = 4"),
       lwd=1,bty="n",cex=0.7)
title("Extensão do Tumor (extent)")

###########################Para variável surg
adj <-coxph(Surv(time[surg==0], status[surg==0])~1,data=death,x=T,
            method="breslow")
ss<- survfit(adj)
s0<-round(ss$surv,digits=5)
H0<- -log(s0)
plot(ss$time,log(H0),xlab="Tempos",
     ylab=expression(log(Lambda[0]*(t))), bty="n",type="s")
adj2 <-coxph(Surv(time[surg==1], status[surg==1])~1,data=death,x=T,
             method="breslow")
ss<- survfit(adj2)
s0<-round(ss$surv,digits=5)
H0<- -log(s0)
lines(ss$time,log(H0),type="s",lty=2)
legend(1500,-3,lty=c(2,1),c("surg = 0","surg = 1"),lwd=1,bty="n",cex=0.7)
title("Tempo Cirurgia (surg)")

##############################Para variável nodesc
adj <-coxph(Surv(time[nodesc==0], status[nodesc==0])~1,data=death,x=T,
            method="breslow")
ss<- survfit(adj)
s0<-round(ss$surv,digits=5)
H0<- -log(s0)
plot(ss$time,log(H0),xlab="Tempos",
     ylab=expression(log(Lambda[0]*(t))), bty="n",type="s")
adj2 <-coxph(Surv(time[nodesc==1], status[nodesc==1])~1,data=death,x=T,
             method="breslow")
ss<- survfit(adj2)
s0<-round(ss$surv,digits=5)
H0<- -log(s0)
lines(ss$time,log(H0),type="s",lty=2)
legend(1500,-3,lty=c(2,1),c("nodesc = 0","nodesc = 1"),lwd=1,bty="n",cex=0.7)
title("+4 nódulos linfáticos")

######Resíduos de Cox-Snell#######
resid(mod, type = "scaledsch")
cox.zph(mod, transform = "identity")
par(mfrow = c(3, 2))
plot(cox.zph(mod))

####Resíduos Martingal e Deviance
par(mfrow=c(1,2))
rd<-resid(mod,type="deviance")       # resíduos deviance
rm<-resid(mod,type="martingale")     # resíduos martingal
pl<-mod$linear.predictors
plot(pl,rm, xlab="Preditor linear", ylab="Resíduo martingal", pch=16)
plot(pl,rd,  xlab="Preditor linear", ylab="Resíduo deviance" , pch=16)
