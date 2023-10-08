library(astsa)


## Gerando um série
set.seed(2023)
N = 2^10
t= 1:N
tendencia  =  15 - 50*(t/N-0.5) + 50*(t/N-0.5)^2 + 300*(t/N-0.5)^3 + 5*(t/N-0.4)^4
serie = tendencia + 2.0*sin(2*pi*t*1/40) + 4.0*cos(2*pi*t*1/40) +
  -2.5*sin(2*pi*t*1/24) + 2.5*cos(2*pi*t*1/24) +
  3.0*sin(2*pi*t*1/150) + 1.0*cos(2*pi*t*1/150) + 
  arima.sim(n = N, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)),
            sd = 2.5)
par(mfrow = c(2,1))
plot(serie, type = "l")

idx = sample(1:N, 250)
serie[idx] = serie[idx] + 1.0 * serie[idx] * sample(c(-1,1), 250, T)
serie[idx] = serie[50:70] * seq(2., 2.5, length.out=11)

plot(serie, type = "l")


# removendo a tendência

est_tendencia = lowess(serie, f = 0.25)
serie_Detrend = serie - est_tendencia$y

par(mfrow = c(2,1))
plot(serie, type = "l")
lines(est_tendencia$y, col = 2, lwd = 2)
plot(serie_Detrend, type = "l")



## Estimando os periodogramas cru:

par(mfrow = c(1,1))
period_cru = astsa::mvspec(serie_Detrend, detrend = F)
1/period_cru$freq[order(period_cru$spec, decreasing = T)[1:5]]


## Estimando periodograma suavisado:

kernel_01 = kernel("daniell", 7)
kernel_02 = kernel("modified.daniell", c(7,7))
kernel_03 = kernel("daniell", c(5,5))
kernel_04 = kernel("daniell", c(5,6,6,5))
kernel_05 = kernel("modified.daniell", c(5,6,5))
kernel_06 = kernel("modified.daniell", c(3,5,7,5,3))

par(mfrow = c(3,2))
plot(kernel_01)
plot(kernel_02)
plot(kernel_03)
plot(kernel_04)
plot(kernel_05)
plot(kernel_06)


par(mfrow = c(3,2))
period_suav01 = astsa::mvspec(serie_Detrend, detrend = F, kernel = kernel_01, taper = 0.0)
period_suav02 = astsa::mvspec(serie_Detrend, detrend = F, kernel = kernel_02, taper = 0.0)
period_suav03 = astsa::mvspec(serie_Detrend, detrend = F, kernel = kernel_03, taper = 0.0)
period_suav04 = astsa::mvspec(serie_Detrend, detrend = F, kernel = kernel_04, taper = 0.0)
period_suav05 = astsa::mvspec(serie_Detrend, detrend = F, kernel = kernel_05, taper = 0.0)
period_suav06 = astsa::mvspec(serie_Detrend, detrend = F, kernel = kernel_06, taper = 0.0)

## o parâmetro taper aplica o filtro "Cosine Taper" à série no dominio do tempo e 
## depois o periodograma é calculado. É um tipo suavisador de covariancia. ver stoffer ed4 pag 201 e morettin 2ed pag 246

x = serie_Detrend 
p = 0.40
nr = length(x)

m   = floor(nr * p)
w   = 0.5 * (1 - cos(pi * seq.int(1, 2 * m - 1, by = 2)/(2 * m)))
nuc = c(w, rep_len(1, nr - 2 * m), rev(w))
x = nuc * x

x == c(spec.taper(serie_Detrend, p = 0.4))
u2 = (1 - (5/8) * p * 2)

par(mfrow = c(3,1))
plot(serie_Detrend, type = "l")
plot(nuc, type = "l", main = paste0("Cossine Bell, p = ", p))
plot(x, type = "l")

par(mfrow = c(2,1))
period_suav01 = astsa::mvspec(x, demean = F, detrend = F, taper = 0.0, plot = F)
period_suav02 = astsa::mvspec(serie_Detrend, demean = F, detrend = F, taper = 0.4, plot = F)

plot(period_suav01$spec/u2 ~ period_suav01$freq, type = "l")
plot(period_suav02$spec ~ period_suav02$freq, type = "l")
period_suav02$spec ==  period_suav01$spec/u2 


par(mfrow = c(3,2))
period_suav01 = astsa::mvspec(serie_Detrend, detrend = F, taper = 0.5)
period_suav02 = astsa::mvspec(serie_Detrend, detrend = F, taper = 0.4)
period_suav03 = astsa::mvspec(serie_Detrend, detrend = F, taper = 0.3)
period_suav04 = astsa::mvspec(serie_Detrend, detrend = F, taper = 0.2)
period_suav05 = astsa::mvspec(serie_Detrend, detrend = F, taper = 0.1)
period_suav06 = astsa::mvspec(serie_Detrend, detrend = F, taper = 0.0)


par(mfrow = c(3,1))
fac = acf(serie_Detrend, lag.max = 512, plot = F)
fac$acf
period_Cov = astsa::mvspec(fac$acf, demean = F, detrend = F, taper = 0.0)
idx = order(period_Cov$details[,3], decreasing = T)
period_Cov$details[idx[1:5], 2]

k = kernel("daniell", m=c(7))
ser2 = kernapply(c(fac$acf), k)
period_SuvCov = astsa::mvspec(ser2, demean = F, detrend = F, taper = 0.0)
idx = order(period_SuvCov$details[,3], decreasing = T)
period_SuvCov$details[idx[1:5], 2]

k = kernel("modified.daniell", m=c(6, 6))
ser2 = kernapply(c(fac$acf), k)
period_SuvCov = astsa::mvspec(ser2, demean = F, detrend = F, taper = 0.0)
idx = order(period_SuvCov$details[,3], decreasing = T)
period_SuvCov$details[idx[1:5], 2]

## COMBINANDO AMBOS:

par(mfrow = c(3,2))
period_suav01 = astsa::mvspec(serie_Detrend, detrend = F, kernel = kernel_01, taper = 0.5)
period_suav02 = astsa::mvspec(serie_Detrend, detrend = F, kernel = kernel_02, taper = 0.5)
period_suav03 = astsa::mvspec(serie_Detrend, detrend = F, kernel = kernel_03, taper = 0.5)
period_suav04 = astsa::mvspec(serie_Detrend, detrend = F, kernel = kernel_04, taper = 0.5)
period_suav05 = astsa::mvspec(serie_Detrend, detrend = F, kernel = kernel_05, taper = 0.5)
period_suav06 = astsa::mvspec(serie_Detrend, detrend = F, kernel = kernel_06, taper = 0.5)


