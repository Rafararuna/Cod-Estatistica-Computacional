dados <- c(2.8660833, 2.9974002, 3.7438357, 1.8547118, 1.8265233, 3.4525441,
       1.3621581, 3.8131878, 2.0969880, 2.7967733, 1.7850382, 1.5230772,
       1.3722042, 1.1855542, 3.1024554, 1.0848315, 2.2334891, 1.3668661,
       1.9122130, 1.7089584, 0.7841663, 1.9528258, 4.0061368, 2.0407105,
       0.8405855, 5.1794931, 1.6179411, 4.3789515, 1.4596143, 2.1202899)

# 1) Estima��o dos parametros usando o m�todo BFGS atrav�s da fun��o optim:

## Definindo vari�veis globais:

n <- length(dados) # tamanho da amostra
alpha <- 0.5 # valor verdadeiro de alpha
beta <- 2 # valor verdadeiro de beta

df <- data.frame( "n" = c("30"), "alpha" = c(0.5), "beta" = c(2.0))

library(knitr)
knitr::kable(df,main = "Vari�veis Globais",align = "c")


## Implementando a fun��o de log-verossimilhan�a:

logL <- function(theta,x) {
  
  n     <- length(dados)
  alpha <- theta[1]
  beta  <- theta[2]
  
  loglike <- n/(alpha^2) - 1/(2*alpha^2) * sum(x/beta+beta/x) -
    n*log(alpha) - (n/2)*log(beta) + sum(log(x+beta))
  
  return(-loglike) # a fun��o optim, por padr�o, faz a minimiza��o, por isso colocamos o sinal de -, para que ela fa�a a maximiza��o
}

## Temos que os valores iniciais s�o:

alpha_0 <-  0.1
beta_0 <- 1.0

## Implementando a fun��o optim usando o m�todo BFGS:

resulta <- optim(c(alpha_0,beta_0), logL, x = dados, method = "BFGS")
resulta$par # estimativas de alpha e beta respectimente


## Temos que os valores iniciais s�o:

alpha_0 <-  1.0
beta_0 <- 1.0

## Implementando a fun��o optim usando o m�todo BFGS:

resulta <- optim(c(alpha_0,beta_0), logL, x = dados, method = "BFGS")
resulta$par # estimativas de alpha e beta respectimente


# 2) Estima��o dos parametros usando o m�todo de Newton:

NR <- function(x, start, tol = 1e-10){
  
  theta <- start
  ch <- 1
  
  while(ch > tol){ 
    
    n <- length(dados)
    alpha <- theta[1]
    beta <- theta[2]
    
    l_alpha <- (-2*n/alpha^3) + ((1/alpha^3)*(sum((x/beta) + (beta/x)))) - n/alpha
    
    l_beta <- (1/(2*(alpha^2))*(sum((x/(beta^2))-(1/x)))) - (n/(2*beta)) + (sum(1/(x+beta)))
    
    gradiente <- c(l_alpha,l_beta)  
    
    l_alpha_alpha <- ((6*n)/(alpha^4)) - ((3/(alpha^4))*(sum((x/beta) + (beta/x)))) + (n/(alpha^2))
    
    l_alpha_beta <- (1/(alpha^3))*(sum((1/x) - (x/(beta^2))))
    
    l_beta_alpha <- (1/(alpha^3))*(sum((1/x) - (x/(beta^2))))
    
    l_beta_beta <- ((-1/((alpha^2)*(beta^3)))*(sum(x))) + (n/(2*(beta^2))) - (sum(1/((x+beta)^2))) 
    
    hessiana <- matrix(c( l_alpha_alpha,l_beta_alpha,l_alpha_beta,l_beta_beta), nrow = 2)
    
    theta = theta - gradiente%*%solve(hessiana)
    
    ch <- sqrt(t(gradiente)%*%gradiente)
    
  }
  
  return(theta)  
  
}  

v_iniciais <- c(0.1,1.0)

NR(dados, v_iniciais)
