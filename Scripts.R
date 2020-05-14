rm(list = ls())
library(maxLik)
x=c(125,18,20,34)

logvero=function(teta){
  theta=teta[1]
  x[1]*log(2+theta)+(x[2]+x[3])*log(1-theta)+x[4]*log(theta)
}
theta.hat=maxLik(logLik = logvero,start = c(theta=0.5))
coef(theta.hat);theta.hat

##-----------------------------------------------------------------
#Aplicação 1
erro=0.00001;dif=0.1;teta.m=0.5; it=0
while (dif>erro) {
  teta.m1=(68+159*teta.m)/(144+197*teta.m)
  dif=(teta.m1-teta.m);
  it=it+1;it;print(dif); print(teta.m1)
  teta.m=teta.m1
}
##-----------------------------------------------------------------
#Aplicação 2 Mistura de normais
rm(list = ls())
n=100
m1=100
m2=120
s1=4
s2=6
pr=0.7
w=rbinom(n,1,pr)
x=rnorm(n,m1,s1)*(w==0)+rnorm(n,m2,s2)*(w==1) 
hist(x)
#funcao para calcular densidade
fc=function(x,m1,m2,s1,s2,prob){
  fx=(1-prob)*dnorm(x,m1,s1)+prob*dnorm(x,m2,s2)
  return(fx)
}
grad=seq(60,200,by=0.1)
fx=fc(grad,m1,m2,s1,s2,pr)
#plot(grad,fx,type="l")
hist(x,ylim = range(fx),prob=T)
lines(grad,fx)

#param[1] Media 1
#param[2] Media 2
#param[3] variancia 1
#param[4] variancia 2
#param[5] epsilon (peso da componente 2)

#y=c(195,166,188,195,179,198,161,179,200,191)
param=c(180,100,15,20,0.5) #chute inicial
e=0.0001
erro=0.5; it=0
while(erro>e){
  param0=param
  #Passo E
  part1 = (1-param[5])*dnorm(x,param[1],param[3])
  part2 = param[5]*dnorm(x,param[2],param[4])
  gam = part2/(part1+part2)
  #passo M
  aux= c(sum( (1-gam)*x)/sum(1-gam),sum(gam*x)/sum(gam))
  param[1] = min(aux)
  param[2] = max(aux)
  param[3] = sqrt(sum( ((1-gam)*(x-param[1])^2) )/sum(1-gam))
  param[4] = sqrt(sum( (gam*(x-param[2])^2) )/sum(gam))
  param[5] = mean(gam)
  #erro = sum((param0-param)^2)/(sum(param))^2
  erro = max(abs(param0-param)/(abs(param0)+0.001))  
  it=it+1
}
round(param,3); it

##-----------------------------------------------------------------
#Exemplo probit
rm(list=ls())
resultado = matrix(0,nrow=1000,ncol=3)
for(jj in 1:1000){
  n=500
  beta0=0.2
  beta1=1.2
  beta2=-0.6
  x1=rbinom(n,1,0.5)
  x2=rnorm(n)
  pred=beta0+beta1*x1+beta2*x2
  prob=pnorm(pred)
  y=rbinom(n,1,prob)

  x=cbind(1,x1,x2)
  beta=array(0,c(3,1))    
  e = 0.0001 ;it=0 ; erro = 0.5
  while(erro > e){
    beta.ini = beta
    xb = x%*%beta
    EZ_aux1= xb - ( dnorm(-xb)/pnorm(-xb) )      # se Zi <= 0
    EZ_aux2= xb + ( dnorm(-xb)/(1-pnorm(-xb)) )  # se Zi > 0
    EZ = EZ_aux1*(y==0) + EZ_aux2*(y==1)        # vetor de esperancas 
    beta = solve( t(x)%*%x ) %*% t(x)%*%EZ      # estimativas na iteracao m
    erro = max(abs(beta.ini - beta)/(abs(beta.ini)+0.001))
    it = it+1 
  }
#it;round(t(beta),3); cbind(beta0,beta1,beta2)
resultado[jj,]=t(beta)
}
mean(resultado[,1]);mean(resultado[,2]);mean(resultado[,2])

x11()
par(mfrow=c(1,3))
hist(resultado[,1])
abline(v=beta0,col=2)
hist(resultado[,2])
abline(v=beta1,col=2)
hist(resultado[,3])
abline(v=beta2,col=2)
##-----------------------------------------------------------------
#MCEM da Aplicação 1
n=1000
erro=0.00001;dif=0.1;
teta.m=0.5; it=0
while (dif>erro) {
  teta.aux= (teta.m/(2+teta.m) )
  z2=rbinom(n,125,teta.aux)
  EZ=mean(z2)
  teta.m1 = (34+EZ)/(72+EZ)
  dif=abs( (teta.m1-teta.m)/teta.m1 );
  teta.m=teta.m1
  print(c(it,teta.m))
  it=it+1
}
##------------------------------------------------------------------------
#Exemplo probit
rm(list=ls())
library(msm)
  n=500
  m=c(rep(10,30),rep(50,20),rep(100,20),rep(500,15),rep(1000,15))
  r=100
  beta0=0.2
  beta1=1.2
  beta2=-0.6
  x1=rbinom(n,1,0.5)
  x2=rnorm(n)
  pred=beta0+beta1*x1+beta2*x2
  prob=pnorm(pred)
  y=rbinom(n,1,prob)
  x=cbind(1,x1,x2)
  beta=array(0,c(3,1))    
  e = 0.0001 ;it=0 ; erro = 0.5 ; result=array(0,c(r,3))
#  while(erro > e){
for(k in 1:r){
    beta.ini = beta
    xb = x%*%beta
    Z_aux1= t(sapply(1:n, function(i){rtnorm(m[k], mean=xb[i], sd=1, lower=-Inf, upper=0)}) ) # se Zi <= 0
    Z_aux2= t(sapply(1:n, function(i){rtnorm(m[k], mean=xb[i], sd=1, lower=0, upper= Inf)}) ) # se Zi > 0
    EZ_aux1=rowMeans(Z_aux1)
    EZ_aux2=rowMeans(Z_aux2)
    EZ = EZ_aux1*(y==0) + EZ_aux2*(y==1)        # vetor de esperancas 
    beta = solve( t(x)%*%x ) %*% t(x)%*%EZ      # estimativas na iteracao m
    #erro = max(abs(beta.ini - beta)/(abs(beta.ini)+0.001))
    result[k,]=t(beta) ; #print(k)
  }
#  it;round(t(beta),3); cbind(beta0,beta1,beta2)

