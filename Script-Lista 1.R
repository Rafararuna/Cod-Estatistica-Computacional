# Questão 1:

n <- 10000
u <- runif(n)
x <- tan(pi*u-pi/2)+1

ggplot()+
  geom_histogram(aes(x, y=stat(density)), fill = "black", color = "white")+
  labs(x="x", y='Frequência')+
  ggtitle("                                            Histograma da f(x)")+
  stat_function(fun = function(x) (1/pi*(1/(1+(x-1)^2))), color = "red")+
  xlim(-10,10)+
  theme_minimal()
##################################################################

# Questão 2:

n <- 100
a <- 2
b <- 2
u <- runif(n)
x <- b/((1-u)^(1/a))
print(x)

seq <- seq(0, 100, 0.1)
hist(x, col = "black", border = "white", probability = TRUE,
     main = "Histograma da f(x)", ylab = "Densidade")
lines(seq, 8/seq^3, col = "red")
###################################################################

# Questão 3:

## Temos que Y ~ U[0,1] é a nossa densidade candidata e X ~ Beta é a nossa densidade alvo
## Primeiro, temos que achar o M, e para isso fazemos: f(x)/g(x) = dbeta/dunif <= M ; 
#  porém, sabamos que g(x) = dunif = 1, logo dbeta <= M, ou seja, vamos realizar apenas o optmize do dbeta para achar  M

f = function(x) {
  dbeta(x, 3, 5)
}
optimize(f, interval = c(0,1), maximum = TRUE)$objective #M = 2.3045

## Agora, geramos nossas NPA's. Porém, antes, é importante mostrar que:

#U <= (1/M) * f(y)/g(y) ; porém, g(y) é igual a 1, logo U <= (1/M) * f(y), ou seja U * M <= f(y), ou seja U[0,M] <= f(y)

## Sabendo disso, podemos gerar as NPA's:

N <- 1000
alpha <- 3
beta <- 5
M <- 2.3

u <- runif(N, min = 0, max=M) #geramos a uniforme U[0,M]
y <- runif(N) #candidata
x <- y[u<dbeta(y, alpha, beta)]

w <- seq(0,1,0.01)
hist(x,probability = TRUE, ylim = c(0,3), col = "black", border = "white", 
     main = "Histograma da f(x)", ylab = "Densidade")
lines(w,dbeta(w,alpha,beta),col = c("red"))
######################################################################

#Questão 4:

n <- 1000
v <- 6 #numero de normais padrao

x_v <- matrix(rnorm(n*v),n,v)^2

### Método 1: Soma por linhas da matriz anterior
y <- rowSums(x_v) #esse somatório das linhas gera os valores da qui-quadrado com 6 g.l
mean(y) #tem que ser próximo de 6, pois a media da qui-quadrado é igual a "v"
var(y) #tem que ser próximo de 12, pois a variancia da qui-quadrado é igual a "2*v"

##ou

### Método 2
y <- apply(x_v, MARGIN=1, FUN = sum) #vector de dimensão n ; o argumento MARGIN define se vc quer somar linha ( = 1) ou coluna ( = 2)
mean(y) #deverá ser próximo de v=6
var(y) #deverá ser próximo de 2*v = 12

hist(y,probability = TRUE, ylim = c(0,0.15), col = "black", border = "white", 
     main = "Histograma da f(x)", ylab = "Densidade")



n <- 1000
v <- 6 #numero de normais padrao

x_v <- matrix(runif(n*v),n,v)

x_v2 <- vector()
for (i in 1:nrow(x_v)) {
  y <- x_v[i,1]*x_v[i,2]*x_v[i,3]
  x_v2 <- append(x_v2, y)
}

x_v3 <- -2*log(x_v2)
mean(x_v3)
var(x_v3)

hist(x_v3,probability = TRUE, ylim = c(0,0.15),col = "black", border = "white", 
     main = "Histograma da f(x)", ylab = "Densidade")