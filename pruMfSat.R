satw <- scan(file="uf50-01.cnf", what = character())
nvar=as.numeric(satw[27])
nk=as.numeric(satw[28])
satw=as.numeric(satw[29:(nk*4+28)])
print("num restricciones")
print(nk)
print("num de variables")
print(nvar)

w <- matrix(nrow = nk, ncol = nvar)
w[,] = 0
r=1
ls = length(satw)
for(s in 1:ls){
  if(satw[s] == 0) {
    r = r+1
  } else{
    n = abs(satw[s])
    w[r,n] = sign(satw[s])
  }
}
x = vector(mode = "numeric", length=nvar)
m = vector(mode = "numeric", length=nvar)
mu = vector(mode = "numeric", length=nk)
g = vector(mode = "numeric", length=nk)
negations = vector(mode = "numeric", length=nk)
### Initial mu vector: ###
mu[] = 0
mnviol=10000000000
print(mu)
for(i in 1:100000){
  m = 1 / (1 + exp(as.vector(1-mu%*%w)))+.00000001
  xb = round(m)
  x=xb
  falsos=which(x==0)
  for(i in 1:length(falsos)){
    x[falsos[i]]=-1
  }
  rclau=as.vector(w%*%x)
  #print(rclau)
  nviol=0
  for(i in 1:nk){
    if(rclau[i]==-3){
      #print("violada")
      nviol=nviol+1
      mu[i]=(mu[i]-(1/10^4)*(as.vector(w[i,]%*%x)))
    #  print(mu[i])
    }
    else{
      #print("satisfecha")
    #m[i]=mu[i]+1/10^3
     #mu[i]=abs(mu[i]-.5)
     #mu[i]=(mu[i]+(1/10^10)*abs(as.vector(w[i,]%*%x)))
    }
  #  print(mu)
  }
  #print("Initial mean field solution:")
  #print(m)
  #print(x)
  #print("---")
  ### Violated constraints: ###
  if(nviol<mnviol){
    mx=x
    mnviol=nviol
    print("Number of violated constraints:")
    print(mnviol)
  #  print(mu)
    bmu=mu
  }

}
#print(bmu)
print(mu)
