repeticiones=1
numinst=10
nPB=3  ##almacenar los n mejores individuos en PB
sink("pru-feb17-Res.txt")
#libRsrary(lpSolve)
###
EvaluaBeneficio <- function(x){
		return(t(x) %*% b)
}

###
EvaluaFactibilidad <- function(x){
		sumF =  w %*% x
		factibilidad=(pm-sumF)

		if(sum(factibilidad>=0)==m){
			return(1)#factible
		}
		else{
			return(0)#infactible
		}
}

generadorSol <- function(imu){
	summuw= as.vector(imu%*%w)
	sol = 1 / (1 + exp((-b+summuw)/Temp^2))
	return(sol)
}

generaFq<-function(sol,mmu){
	Fq=0
	ss=0
	ent=0#entropía
	Fq=mmu%*%(pm-w %*% sol)
	ben=EvaluaBeneficio(sol)
	Fq=ben+Fq
	for(i in 1:n){
		ss=(1-sol[i])*log(1-sol[i])+sol[i]*log(sol[i])
		if(ss!="NaN"){ent=ent+ss}
	}
	Fq=Fq+ent
	return(Fq)
}


	optimiza=function(){
			for(i in 1:m){
				dir.binario[i]="<="
			}
			ld= re[(n*m+n+4):(n*m+n+m+3)]
			MKP<- lp ("max",b, w, dir.binario, ld, compute.sens=TRUE, binary.vec=c(1:n))
			resoptim[1]=MKP$objval
			resoptim[2:(n+1)]=MKP$solution
			return(resoptim)
	}



###############################################################################################################
casos = c("OR5x100", "OR5x250", "OR5x500", "OR10x100", "OR10x250", "OR10x500",
"OR30x100", "OR30x250", "OR30x500")
densidad=c("-0.25_","-0.50_","-0.75_")
optimos=scan("optimosMKP.txt")
nn=0
#ExpResprom nombre del archivo que contiene los promedios
archivo = paste("abril16-prom.txt", sep="")
etiq=c("Instancia", "PROMgapcm","PROMgapmu","PROMtcm","PROMtm","PROMtt","PROMgapccm","PROMgapcmu","SDgapcm","SDgapmu","SDtcm","SDtm","SDtt","SDgapccm","SDgapcmu")
col = length(etiq)
resultados=matrix(data = 0, nrow = repeticiones, ncol = (col-1)/2)
estadisticas = vector(mode = "numeric", length=col)
PROMdensidad=matrix(data = 0, nrow = numinst, ncol = col-1)
ESTdensidad=matrix(data = 0, nrow = length(densidad), ncol = col)
write("RESULTADOS PROMEDIO MKP",archivo, sep=" ", append=FALSE, ncolumns=col)
for (cc in 1:length(casos)){#inicio de casos
		caso = casos[cc]
		nombre=c("\n",caso)

		write(nombre,archivo, sep=" ", append=TRUE, ncolumns=col)
		write(etiq,archivo, sep=" ", append=TRUE, ncolumns=col)
		for (ii in 1:length(densidad)){#densidad
			den=densidad[ii]
				for(jj in 1:numinst){# número de instancia de 1 a 10
					nn=nn+1
					instancia=paste(caso,"/",caso,den,jj,".dat",sep="")
					instn=paste(caso,den,jj,".dat",sep="")
					re <- scan(instancia)
					#re <- scan("ins5-2.txt")
					n = re[1] ## núm de variables ##
					m = re[2] ## núm de restricciones
					opt = re[3] ## valor óptimo
					b = vector(mode = "numeric", length=n)
					w = matrix(data = NA, nrow = m, ncol = n)
					pm = vector(mode = "numeric", length=m)
					mui = vector(mode = "numeric", length=m)
					muf = vector(mode = "numeric", length=m)
					muc = vector(mode = "numeric", length=m)
					mu = vector(mode = "numeric", length=m)
					SumF = matrix(data = NA, nrow = m, ncol = 1)
					smuw = vector(mode = "numeric", length=m)
					cero = vector(mode = "numeric", length=m)
					sum = vector(mode = "numeric", length=m)
					sol = vector(mode = "numeric", length=n)
					A = vector(mode = "numeric", length=n)
					PB = matrix(data = 0, nrow = nPB, ncol = n)
					pos_unos  = vector(mode = "numeric", length=n)
					pos_ceros = vector(mode = "numeric", length=n)
					resbl = vector(mode = "numeric", length=n+1)
					resoptim = vector(mode = "numeric", length=n+1)
					solmuestra = vector(mode = "numeric", length=n)
					variabilidad = vector(mode = "numeric", length=n)
					varia = vector(mode = "numeric", length=n)

					dir.binario = vector(mode = "numeric", length=m)
					ld = vector(mode = "numeric", length=m)
					pm = re[(n*m+n+4):(n*m+n+m+3)]
					b = re[4:(n+3)]
					aux = n+4
					for(i in 1:m){
						for(j in 1:n){
							w[i,j] = re[aux]
							aux= aux + 1
						}
					}

					for(kk in 1:repeticiones){##repeticiones


##########################
#Inicio del algoritmo de CM
###########################

#############################
Temp=1
alpha=1/10^(7)
mu = vector(mode = "numeric", length=m)
priomu = vector(mode = "numeric", length=m)
con=0
mejorval=0
#print(mu)
set.seed(Sys.time())
ptm <- proc.time()
ptmo = ptm
while(TRUE){
	scmc=generadorSol(mu)
	x = round(scmc)
	factibilidad=(pm-w %*% x)
	if(sum(factibilidad>=0)==m){
		#break
		con=con+1
		val=EvaluaBeneficio(x)
		if(val>mejorval){
			mejorval=val
			gen=con
		}
		if(con==50){
			break
		}
	}
	for(i in 1:m){#factible
		if(factibilidad[i]>=0){
			mu[i]=mu[i]-mu[i]*.01
		}
		else{#infactible
			mu[i]=mu[i]-alpha*(pm[i]-w[i,] %*% x)
		}
	}

}
ptm = proc.time() - ptmo
tcm=ptm[1]

#####################################################
cat(instn, mejorval,gen,"\n")
#####################################################

}#repeticiones

}#numinst
}#densidades


}#casos
